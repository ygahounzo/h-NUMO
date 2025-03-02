!---------------------------------------------------------------------!
!>@brief This routine writes out the all 8 Global Variables: Density, velocity, potential temperature,
!> water vapor, cloud water, and rain water in VTK format
!> This data is then read into Paraview or VisIT
!>@author James F. Kelly 5 Feb. 2011
!>
!>@date S. Gopalakrishnan
!> 26 Mar. 2011
!> 
!>@date Simone Marras
!> 21 May 2013
!>
!>Adding binary output version
!>@date F.X. Giraldo
!>@date February 27 2014 to include Reference Pressure
!>@date April 23 2014 to include Temperature
!>@date May 1 1014 to include Radius (for Spherical Cases)
!>
!>@date A. Muller
!>June 2014 to include outasciimaya for CSV output to be fed to Maya Autodesk
!>@date F.X. Giraldo on October 16, 2016: rewritten to remove all LES, VMS, SMAG, NAZ stuff and to only handle LAV
!---------------------------------------------------------------------!


subroutine outvtk_g_binary_mlswe(q,qb,fname,time)

    use mod_basis, only: ngl, nglx, ngly, nglz, is_2d

    use mod_constants, only: pi, earth_radius, gravity
  
    use mod_global_grid, only: coord_g, intma_g, npoin_g, nelem_g, ncol_g

    use mod_grid, only: npoin, intma, coord, nelem
  
    use mod_initial, only: nvar, nvar_diag, kvector, rho_layers, bathymetry

    use mod_input, only: nelx, nely, nelz, nopx, nopy, nopz, out_type, &
        eqn_set, format_vtk, space_method, is_mlswe

    use mod_mpi_utilities, only: irank, irank0
  
    use mod_parallel, only: nproc, num_send_recv_total

    use mod_vtk_binary
  
    implicit none

    !global variables
    real, intent(in) :: q(nvar,npoin), qb(4,npoin)
    real, intent(in) :: time
    character, intent(in) :: fname*100

    !local variables
    integer ie, i, j, k, nglm1, nglm13, ii, jj, kk
    integer ncells, nsize
    integer ip, ip_g, iproc
    integer Nn, Ne, Nnodi, ncon
    integer E_IO
    real T_k, theta_k, P_k, pi_k, rho_k, E_k
    real T_ref, theta_ref, P_ref, pi_ref, rho_ref, E_ref
    real u_k, v_k, w_k
    real u2_k, v2_k, w2_k
    real u2_ref, v2_ref, w2_ref
    real z
    real time_value
    real sound, Mach
    real c, rho, theta
    integer displs1(nproc)
    integer ierr

    !allocatable local arrays
    real,   dimension(:,:,:), allocatable                  :: q_l
    real,   dimension(:,:),   allocatable                  :: q_g, coord_dg_gathered,qb_g
    real,   dimension(:),     allocatable                  :: km
    real,   dimension(:),     allocatable                  :: x_uns
    real,   dimension(:),     allocatable                  :: y_uns
    real,   dimension(:),     allocatable                  :: z_uns
    integer,dimension(:),     allocatable                  :: eltype
    integer,dimension(:),     allocatable                  :: conn
    real,   dimension(:),     allocatable                  :: var_uns_grid, var_uns_grid_ref

    character*72 :: cbuf
    character*12 :: output_format
    character*24 :: fnp
    integer elemType, l
    integer AllocateStatus
    real :: xfactor, yfactor, zfactor
    logical:: is_cgc
    
    ! is cgc
    is_cgc = .false.
    if(space_method == 'cgc') is_cgc = .true.

    if (irank == irank0) then
        allocate(q_g(nvar,npoin_g), coord_dg_gathered(3,npoin_g), &
            x_uns(npoin_g), y_uns(npoin_g), z_uns(npoin_g), eltype(nelem_g*max(nglx-1,1)*max(ngly-1,1)*max(nglz-1,1)), &
            conn(9*nelem_g*max(nglx-1,1)*max(ngly-1,1)*max(nglz-1,1)), var_uns_grid(npoin_g), var_uns_grid_ref(npoin_g),&
            qb_g(4,npoin_g),stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - OUTVTK_G_BINARY **"
    end if ! irank==irank0
  
    ! Gather Data onto Head node
    call gather_data(q_g,q,nvar)
    call gather_data(qb_g,qb,4)
    
    call gather_data(coord_dg_gathered,coord,3)
    
    ! copy sections of the diagnostic array
    if (irank == irank0) then

        nglm13 = max(nglx-1,1)*max(ngly-1,1)*max(nglz-1,1)   !Number of cells per element
        ncells = nelem_g*nglm13      !Total number of cells

        !
        ! Open VTK file
        !
        call vtk_ini(output_format = format_vtk, &
            filename      = fname,               &
            title         = 'NUMO3d data',     &
            mesh_topology = 'UNSTRUCTURED_GRID', &
            time_value    = time)

        x_uns(:) = coord_dg_gathered(1,:)
        y_uns(:) = coord_dg_gathered(2,:)
        z_uns(:) = coord_dg_gathered(3,:)

        !
        ! Write coordinates to file:
        !
        call vtk_geo_unst_R8(npoin_g, x_uns, y_uns, z_uns)
     
        ! Cell type => Hex or Quads in case of 2D
        if(is_2d) then
            do i = 1,ncells
                eltype(i) = 9
            end do
            nsize = 5*ncells
        else
            do i = 1,ncells
                eltype(i) = 12
            end do
            nsize = 9*ncells
        endif
     
        !
        ! Write connectivity to file:
        !
        l=1

        do ie =1,nelem_g
            do i = 1,max(nglx-1,1)
                do j = 1,max(ngly-1,1)
                    do k = 1,max(nglz-1,1)
                        ii=min(i+1,nglx)
                        jj=min(j+1,ngly)
                        kk=min(k+1,nglz)
                
                        if(nglx == 1) then
                            if(is_cgc) then
                                conn(l)=4
                                conn(l+1)=(intma_g( i, j, k,ie)  -1)
                                conn(l+2)=(intma_g( i,jj, k,ie)  -1)
                                conn(l+3)=(intma_g( i,jj,kk,ie)  -1)
                                conn(l+4)=(intma_g( i, j,kk,ie)  -1)
                                l = l + 5
                            else
                                conn(l)=4
                                conn(l+1)=(intma( i, j, k,ie)   -1)
                                conn(l+2)=(intma( i,jj, k,ie)   -1)
                                conn(l+3)=(intma( i,jj,kk,ie)   -1)
                                conn(l+4)=(intma( i, j,kk,ie)   -1)
                                l = l + 5
                            endif
                        else if(ngly == 1) then
                            if(is_cgc) then
                                conn(l)=4
                                conn(l+1)=(intma_g( i, j, k,ie)  -1)
                                conn(l+2)=(intma_g(ii, j, k,ie)  -1)
                                conn(l+3)=(intma_g(ii, j,kk,ie)  -1)
                                conn(l+4)=(intma_g( i, j,kk,ie)  -1)
                                l = l + 5
                            else
                                conn(l)=4
                                conn(l+1)=(intma( i, j, k,ie)   -1)
                                conn(l+2)=(intma(ii, j, k,ie)   -1)
                                conn(l+3)=(intma(ii, j,kk,ie)   -1)
                                conn(l+4)=(intma( i, j,kk,ie)   -1)
                                l = l + 5
                            endif
                        else if(nglz == 1) then
                            if(is_cgc) then
                                conn(l)=4
                                conn(l+1)=(intma_g( i, j, k,ie)  -1)
                                conn(l+2)=(intma_g(ii, j, k,ie)  -1)
                                conn(l+3)=(intma_g(ii,jj, k,ie)  -1)
                                conn(l+4)=(intma_g( i,jj, k,ie)  -1)
                                l = l + 5
                            else
                                conn(l)=4
                                conn(l+1)=(intma( i, j, k,ie)   -1)
                                conn(l+2)=(intma(ii, j, k,ie)   -1)
                                conn(l+3)=(intma(ii,jj, k,ie)   -1)
                                conn(l+4)=(intma( i,jj, k,ie)   -1)
                                l = l + 5
                            endif
                        else
                            if(is_cgc) then
                                conn(l)=8
                                conn(l+1)=(intma_g( i, j, k,ie)  -1)
                                conn(l+2)=(intma_g(ii, j, k,ie)  -1)
                                conn(l+3)=(intma_g(ii,jj, k,ie)  -1)
                                conn(l+4)=(intma_g( i,jj, k,ie)  -1)
                                conn(l+5)=(intma_g( i, j,kk,ie)  -1)
                                conn(l+6)=(intma_g(ii, j,kk,ie)  -1)
                                conn(l+7)=(intma_g(ii,jj,kk,ie)  -1)
                                conn(l+8)=(intma_g( i,jj,kk,ie)  -1)
                                l = l + 9
                            else
                                conn(l)=8
                                conn(l+1)=(intma( i, j, k,ie)   -1)
                                conn(l+2)=(intma(ii, j, k,ie)   -1)
                                conn(l+3)=(intma(ii,jj, k,ie)   -1)
                                conn(l+4)=(intma( i,jj, k,ie)   -1)
                                conn(l+5)=(intma( i, j,kk,ie)   -1)
                                conn(l+6)=(intma(ii, j,kk,ie)   -1)
                                conn(l+7)=(intma(ii,jj,kk,ie)   -1)
                                conn(l+8)=(intma( i,jj,kk,ie)   -1)
                                l = l + 9
                            endif
                        endif
                
                    end do
                end do
            end do
        end do
        ncon = l-1

        !
        ! Write connectivity to file
        !
        call vtk_con(ncells, nsize, ncon, conn, eltype)

        !
        ! Write data to file:
        !
        call vtk_dat(npoin_g, 'NODE')

        do i=1,npoin_g
            var_uns_grid(i) = q_g(1,i)
        end do
        call vtk_var_scal_R8(npoin_g,'h',var_uns_grid)
        
        !Write bathymetry below reference level
        ! do i=1,npoin
        ! var_uns_grid(i) = q_g(4,i)
        ! end do
        ! call vtk_var_scal_R8(npoin,'elevation',var_uns_grid)


        !Write u-velo :
        do i=1,npoin_g
            var_uns_grid(i) = q_g(2,i)
        end do
        
        call vtk_var_scal_R8(npoin_g,'u',var_uns_grid)

        !Write v-velo :
        do i=1,npoin_g
            var_uns_grid(i) = q_g(3,i)
        end do

        call vtk_var_scal_R8(npoin_g,'v',var_uns_grid)

        !Write velocity vectors:
        do i=1,npoin_g
            x_uns(i) = q_g(2,i)
            y_uns(i) = q_g(3,i)
            z_uns(i) = 0.0
        end do
        call vtk_var_vect_R8('VECT',npoin_g,'MOMENTUM',x_uns, y_uns, z_uns)

        do i=1,npoin_g
            var_uns_grid(i) = qb_g(1,i)
        end do
        call vtk_var_scal_R8(npoin_g,'pb',var_uns_grid)
            
        !Write bathymetry below reference level
        do i=1,npoin_g
            var_uns_grid(i) = q_g(5,i)
        end do
        call vtk_var_scal_R8(npoin_g,'SSH',var_uns_grid)
    
        !Write u-velo :
        do i=1,npoin_g
            var_uns_grid(i) = qb_g(3,i) / qb_g(1,i)
        end do
        
        call vtk_var_scal_R8(npoin_g,'ub',var_uns_grid)

        !Write v-velo :
        do i=1,npoin_g
            var_uns_grid(i) = qb_g(4,i) / qb_g(1,i)
        end do

        call vtk_var_scal_R8(npoin_g,'vb',var_uns_grid)
    
        !Close the file
        call vtk_end()
     
        !deallocate global arrays
        deallocate(q_g,qb_g)
        deallocate(var_uns_grid, var_uns_grid_ref)
     
    end if !irank0
  
end subroutine outvtk_g_binary_mlswe


subroutine outvtk_g_binary_mlswe_global(q,qb,qprime,fname,time)

    use mod_basis, only: ngl, nglx, ngly, nglz, is_2d

    use mod_constants, only: pi, earth_radius, gravity
  
    use mod_global_grid, only: coord_g, intma_g, npoin_g, nelem_g, ncol_g

    use mod_grid, only: npoin, intma, coord, nelem
  
    use mod_initial, only: nvar, nvar_diag, kvector, rho_layers, bathymetry

    use mod_input, only: nelx, nely, nelz, nopx, nopy, nopz, out_type, &
        eqn_set, format_vtk, space_method, is_mlswe

    use mod_mpi_utilities, only: irank, irank0
  
    use mod_parallel, only: nproc, num_send_recv_total

    use mod_vtk_binary
  
    implicit none

    !global variables
    real, intent(in) :: q(3,npoin), qb(3,npoin), qprime(3,npoin)
    real, intent(in) :: time
    character, intent(in) :: fname*100

    !local variables
    integer ie, i, j, k, nglm1, nglm13, ii, jj, kk
    integer ncells, nsize
    integer ip, ip_g, iproc
    integer Nn, Ne, Nnodi, ncon
    integer E_IO
    real T_k, theta_k, P_k, pi_k, rho_k, E_k
    real T_ref, theta_ref, P_ref, pi_ref, rho_ref, E_ref
    real u_k, v_k, w_k
    real u2_k, v2_k, w2_k
    real u2_ref, v2_ref, w2_ref
    real z
    real time_value
    real sound, Mach
    real c, rho, theta
    integer displs1(nproc)
    integer ierr

    !allocatable local arrays
    real,   dimension(:,:,:), allocatable                  :: q_l
    real,   dimension(:,:),   allocatable                  :: q_g, coord_dg_gathered,qb_g, qprime_g
    real,   dimension(:),     allocatable                  :: km
    real,   dimension(:),     allocatable                  :: x_uns
    real,   dimension(:),     allocatable                  :: y_uns
    real,   dimension(:),     allocatable                  :: z_uns
    integer,dimension(:),     allocatable                  :: eltype
    integer,dimension(:),     allocatable                  :: conn
    real,   dimension(:),     allocatable                  :: var_uns_grid, var_uns_grid_ref

    character*72 :: cbuf
    character*12 :: output_format
    character*24 :: fnp
    integer elemType, l
    integer AllocateStatus
    real :: xfactor, yfactor, zfactor
    logical:: is_cgc
    

    if (irank == irank0) then
        allocate(q_g(3,npoin_g), coord_dg_gathered(3,npoin_g), &
            x_uns(npoin_g), y_uns(npoin_g), z_uns(npoin_g), eltype(nelem_g*max(nglx-1,1)*max(ngly-1,1)*max(nglz-1,1)), &
            conn(9*nelem_g*max(nglx-1,1)*max(ngly-1,1)*max(nglz-1,1)), var_uns_grid(npoin_g),&
            qb_g(4,npoin_g),qprime_g(3,npoin_g), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - OUTVTK_G_BINARY **"
    end if ! irank==irank0
  
    ! Gather Data onto Head node
    call gather_data(q_g,q,3)
    call gather_data(qb_g,qb,4)
    call gather_data(qprime_g,qprime,3)
    
    call gather_data(coord_dg_gathered,coord,3)
    
    ! copy sections of the diagnostic array
    if (irank == irank0) then

        nglm13 = max(nglx-1,1)*max(ngly-1,1)*max(nglz-1,1)   !Number of cells per element
        ncells = nelem_g*nglm13      !Total number of cells

        !
        ! Open VTK file
        !
        call vtk_ini(output_format = format_vtk, &
            filename      = fname,               &
            title         = 'NUMO3d data',     &
            mesh_topology = 'UNSTRUCTURED_GRID', &
            time_value    = time)

        x_uns(:) = coord_dg_gathered(1,:)
        y_uns(:) = coord_dg_gathered(2,:)
        z_uns(:) = coord_dg_gathered(3,:)

        !
        ! Write coordinates to file:
        !
        call vtk_geo_unst_R8(npoin_g, x_uns, y_uns, z_uns)
     
        ! Cell type => Hex or Quads in case of 2D
        if(is_2d) then
            do i = 1,ncells
                eltype(i) = 9
            end do
            nsize = 5*ncells
        else
            do i = 1,ncells
                eltype(i) = 12
            end do
            nsize = 9*ncells
        endif
     
        !
        ! Write connectivity to file:
        !
        l=1

        do ie =1,nelem_g
            do i = 1,max(nglx-1,1)
                do j = 1,max(ngly-1,1)
                    do k = 1,max(nglz-1,1)
                        ii=min(i+1,nglx)
                        jj=min(j+1,ngly)
                        kk=min(k+1,nglz)
                
                        if(nglx == 1) then
                            conn(l)=4
                            conn(l+1)=(intma( i, j, k,ie)   -1)
                            conn(l+2)=(intma( i,jj, k,ie)   -1)
                            conn(l+3)=(intma( i,jj,kk,ie)   -1)
                            conn(l+4)=(intma( i, j,kk,ie)   -1)
                            l = l + 5
                        else if(ngly == 1) then
                            conn(l)=4
                            conn(l+1)=(intma( i, j, k,ie)   -1)
                            conn(l+2)=(intma(ii, j, k,ie)   -1)
                            conn(l+3)=(intma(ii, j,kk,ie)   -1)
                            conn(l+4)=(intma( i, j,kk,ie)   -1)
                            l = l + 5
                        else if(nglz == 1) then
                            conn(l)=4
                            conn(l+1)=(intma( i, j, k,ie)   -1)
                            conn(l+2)=(intma(ii, j, k,ie)   -1)
                            conn(l+3)=(intma(ii,jj, k,ie)   -1)
                            conn(l+4)=(intma( i,jj, k,ie)   -1)
                            l = l + 5
                        else
                            conn(l)=8
                            conn(l+1)=(intma( i, j, k,ie)   -1)
                            conn(l+2)=(intma(ii, j, k,ie)   -1)
                            conn(l+3)=(intma(ii,jj, k,ie)   -1)
                            conn(l+4)=(intma( i,jj, k,ie)   -1)
                            conn(l+5)=(intma( i, j,kk,ie)   -1)
                            conn(l+6)=(intma(ii, j,kk,ie)   -1)
                            conn(l+7)=(intma(ii,jj,kk,ie)   -1)
                            conn(l+8)=(intma( i,jj,kk,ie)   -1)
                            l = l + 9
                        endif
                
                    end do
                end do
            end do
        end do
        ncon = l-1

        !
        ! Write connectivity to file
        !
        call vtk_con(ncells, nsize, ncon, conn, eltype)

        !
        ! Write data to file:
        !
        call vtk_dat(npoin_g, 'NODE')

        ! do i=1,npoin_g
        !     var_uns_grid(i) = q_g(1,i)
        ! end do
        call vtk_var_scal_R8(npoin_g,'dp_df',q_g(1,:))
        call vtk_var_scal_R8(npoin_g,'dpp_df',qprime_g(1,:))
        call vtk_var_scal_R8(npoin_g,'pbpert_df',qb_g(1,:))
        !Write velocity vectors:
        x_uns(:) = q_g(2,:)
        y_uns(:) = q_g(3,:)
        call vtk_var_vect2D_R8('VECT',npoin_g,'VELOCITY',x_uns, y_uns)

        !Write velocity vectors uprime, vprime
        x_uns(:) = qprime_g(2,:)
        y_uns(:) = qprime_g(3,:)
        call vtk_var_vect2D_R8('VECT',npoin_g,'VELOPRIME',x_uns, y_uns)

        !Write velocity vectors uprime, vprime
        x_uns(:) = qb_g(2,:)
        y_uns(:) = qb_g(3,:)
        call vtk_var_vect2D_R8('VECT',npoin_g,'VELOBARO',x_uns, y_uns)
    
        !Close the file
        call vtk_end()
     
        !deallocate global arrays
        deallocate(q_g,qb_g,qprime_g)
        deallocate(var_uns_grid)
     
    end if !irank0
  
end subroutine outvtk_g_binary_mlswe_global