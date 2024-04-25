!----------------------------------------------------------------------!
!>@brief This module builds the Grid Geometry: INTMA, COORD, etc.
!>@author  Francis X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>
!>@date FX Giraldo to add sigma_g          on Jun 2014
!>
!>@date S. Marras  to add bdy_flg_g        on Jun 2014
!>
!----------------------------------------------------------------------!
module mod_global_grid

    use mod_basis, only: nopx, nopy, nopz, nglx, ngly, nglz, &
        xgl, xglx, xgly, xglz, FACE_LEN, &
        nqx, nqy, nqz, xnq, xnqx, xnqy, xnqz

    use mod_constants, only: earth_radius, tol, gravity

    use mod_input, only : equations, icase, nelx, nely, nelz, &
        xdims, ydims, &
        x_boundary, y_boundary, z_boundary, &
        ztop_in => ztop, zbot_in => zbottom, geometry_type, space_method

    use mod_types, only : r8

    public :: &
        mod_global_grid_create, &
        mod_global_grid_bcs, &
        mod_global_grid_dimensions, &
        mod_global_grid_destroy, &
        mod_global_grid_out, &
        mod_global_grid_init_coord, &
        coord_g_cg, &
        coord_g, &
        ! coord_g_cg_q, & ! added by Yao
        ! coord_g_q, &   ! added by Yao
        sigma_g, &
        index_g, &
        intma_g, &
        ! intma_g_q, &   ! added by Yao
        intma_lev_g, &
        intma_s, &
        ele_col_g, &
        bsido_g, &
        iperiodic_g, &
        flag_periodic_g, &
        nface_int, &
        nface_g,face_g, &
        npoin_g_cg, npoin_g, nelem_g, &
        npoin_g_cg_q, npoin_g_q, &   ! This is the number of quadrature points in the global grid
        nboun_g, &          !This is the global number of boundary elements
        nboun_poin_g, &     !This is the global number of boundary points
        ncol_g, &
        nelem_s, npoin_r, npoin_s, &
        nelem_r, &
        xperiodic, yperiodic, zperiodic, &
        nx, ny, nz, &
        xmin, xmax, ymin, ymax, zmin, zmax, &
        iboundary, &
        saved_coord, saved_sigma, saved_index, &
        saved_intma, &
        saved_bsido, &
        saved_face

    private
    !-----------------------------------------------------------------------
    !module variables and parameters
    real,    dimension(:,:),      allocatable :: coord_g_cg, coord_g
    real,    dimension(:),        allocatable :: sigma_g
    integer, dimension(:,:),      allocatable :: index_g, face_g
    integer, dimension(:,:,:,:),  allocatable :: intma_g, intma_lev_g
    integer, dimension(:,:,:),    allocatable :: intma_s
    integer, dimension(:,:),      allocatable :: bsido_g
    integer, dimension(:),        allocatable :: ele_col_g
    real :: xmin, xmax, ymin, ymax, zmin, zmax
    logical xperiodic, yperiodic, zperiodic

    integer :: iboundary(6)
    integer :: nface_int, nface_g
    integer :: npoin_g_cg, npoin_g, nelem_g, nelem_r, npoin_g_cg_q, npoin_g_q
    integer :: nboun_g
    integer :: nboun_poin_g
    integer :: ncol_g
    integer :: nelem_s, npoin_r, npoin_s
    integer :: nx, ny, nz

    real,dimension(:,:,:),       allocatable :: saved_coord
    real,dimension(:,:),         allocatable :: saved_sigma
    integer,dimension(:,:,:),    allocatable :: saved_index, saved_bsido, saved_face
    integer,dimension(:,:,:,:,:),allocatable :: saved_intma
  !-----------------------------------------------------------------------

contains

    !-----------------------------------------------------------------------
    subroutine mod_global_grid_create()

        implicit none

        if (geometry_type == 'cube') then
            call mod_global_grid_cube_create()
        else
            print *, "Invalid Geometry Type"
            stop
        end if

        !Initialize coord.
        call mod_global_grid_init_coord()

        !create faces
        call create_face(face_g,nface_g,FACE_LEN)

    end subroutine mod_global_grid_create

    !-----------------------------------------------------------------------
    subroutine mod_global_grid_init_coord()

        implicit none

        integer i,j,k,e,ip,ip1
        
        logical:: is_cgc
        
        is_cgc = .false.
        if(space_method == 'cgc') is_cgc = .true.

        !Initialize coord.
        do i = 1,nglx
            do j = 1,ngly
                do k = 1,nglz
                    do e = 1,nelem_g
                        ip = intma_g(i,j,k,e)
                        if(is_cgc) then
                            ip1 = ip
                        else
                            ip1 = (e-1) * nglz * ngly * nglx + &
                                (k-1) * ngly * nglx + &
                                (j-1) * nglx + &
                                (i-1) &
                                + 1
                        endif
                        coord_g(:,ip1) = coord_g_cg(:,ip)
                    enddo
                enddo
            enddo
        enddo

    end subroutine mod_global_grid_init_coord

    !-----------------------------------------------------------------------
    subroutine mod_global_grid_cube_create()
        
        implicit none

        !local arrays
        integer :: AllocateStatus
        integer :: i, j
        real :: xc, yc, zc, ac, hc
        real :: x, y, z, ztop
        real, dimension(:), allocatable :: zsurf
        
        real f0, beta, gg, rho, tau

        !choose either CG or DG sotrage
        nx=nelx*nopx + 1
        ny=nely*nopy + 1
        nz=nelz*nopz + 1

        !Local Sizes
        nelem_s = nelx*nely
        npoin_s = nx*ny

        !Global sizes
        nelem_g = nelx*nely*nelz
        npoin_g_cg = nx*ny*nz
        npoin_g_cg_q = ((nqx-1)*nelx + 1)*((nqy-1)*nely + 1)*((nqz-1)*nelz + 1)

        if(space_method == 'cgc') then
            npoin_g = npoin_g_cg
        else
            npoin_g = nelem_g * (nopx + 1) * (nopy + 1) * (nopz + 1)
            npoin_g_q = nelem_g * nqx * nqy * nqz
        endif
        ncol_g  = npoin_g_cg/nz
    
        nboun_g = 0
        if(nopx > 0) nboun_g = nboun_g + 2*nely*nelz
        if(nopy > 0) nboun_g = nboun_g + 2*nelx*nelz
        if(nopz > 0) nboun_g = nboun_g + 2*nelx*nely
    
        nboun_poin_g = 8
        if(nopx > 0) nboun_poin_g = nboun_poin_g + 2*(ny-2)*(nz-2) + 4*(nx-2)
        if(nopy > 0) nboun_poin_g = nboun_poin_g + 2*(nx-2)*(nz-2) + 4*(ny-2)
        if(nopz > 0) nboun_poin_g = nboun_poin_g + 2*(nx-2)*(ny-2) + 4*(nz-2)

        nface_g = 0
        if(nopx > 0) nface_g = nface_g + (nelx + 1)*nely*nelz
        if(nopy > 0) nface_g = nface_g + (nely + 1)*nelx*nelz
        if(nopz > 0) nface_g = nface_g + (nelz + 1)*nelx*nely
 
        !initial grid & bcs
        iboundary(1:2)=z_boundary(1:2)
        iboundary(3:4)=y_boundary(1:2)
        iboundary(5:6)=x_boundary(1:2)
        
        xmin=xdims(1) ;  xmax=xdims(2)
        ymin=ydims(1) ;  ymax=ydims(2)
        zmin=zbot_in  ;  zmax=ztop_in
        
        !Determine Periodicity
        xperiodic=.false.
        yperiodic=.false.
        zperiodic=.false.

        !THIS Works: Note that iboundary == 30 in order to let COUNT_FACES count the right number 
        !that CREATE_FACE counts since Periodic Faces are counted twice since they are seen by CREATE_FACE 
        !as boundary faces but in reality are interior faces.
        if (iboundary(1) == 3 .and. iboundary(2) == 3) zperiodic=.true. 
        if (iboundary(3) == 3 .and. iboundary(4) == 3) yperiodic=.true. 
        if (iboundary(5) == 3 .and. iboundary(6) == 3) xperiodic=.true. 

        print*,' Boundary Conditions are: '
        do i=1,6
            print*,' i, iboundary = ',i,iboundary(i)
        end do
        print*,' xperiodic = ',xperiodic
        print*,' yperiodic = ',yperiodic
        print*,' zperiodic = ',zperiodic

        allocate( coord_g_cg(3,npoin_g_cg),coord_g(3,npoin_g), face_g(FACE_LEN,nface_g), &
            sigma_g(npoin_g_cg), index_g(2,npoin_g_cg),  &
            intma_g(nglx,ngly,nglz,nelem_g), &
            bsido_g(6,nboun_g), &
            zsurf(npoin_g_cg), &
            ele_col_g(nelem_g), &
            ! coord_g_cg_q(3,npoin_g_cg_q),coord_g_q(3,npoin_g_q), &   ! added by Yao
            ! intma_g_q(nqx,nqy,nqz,nelem_g), &  ! added by Yao
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Grid **"

        !Initialize allocated arrays:
        coord_g_cg      = 0.0
        coord_g         = 0.0
        sigma_g         = 0.0
        index_g         = 0
        intma_g         = 0
        bsido_g         = 0
        zsurf           = 0.0
        ele_col_g       = 0
        ! coord_g_cg_q      = 0.0
        ! coord_g_q         = 0.0
        ! intma_g_q         = 0
        
        !Construct Grid
        call create_grid_cube(coord_g_cg,index_g,intma_g,ele_col_g,bsido_g,npoin_g_cg,npoin_g,nelem_g,nboun_g,&
            xglx,xgly,xglz,nglx,ngly,nglz,nelx,nely,nelz,xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,&
            iboundary,xperiodic,yperiodic,zperiodic)

end subroutine mod_global_grid_cube_create

!-----------------------------------------------------------------------!
subroutine mod_global_grid_bcs()

    use mpi

    implicit none

    integer :: ierr

    call mpi_bcast(xperiodic,1,mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(yperiodic,1,mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(zperiodic,1,mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(iboundary,6,mpi_integer,0,mpi_comm_world,ierr)

end subroutine mod_global_grid_bcs

!-----------------------------------------------------------------------!
subroutine mod_global_grid_dimensions()

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    integer :: ierr

    call mpi_bcast(xmin,1,MPI_PRECISION,0,mpi_comm_world,ierr)
    call mpi_bcast(xmax,1,MPI_PRECISION,0,mpi_comm_world,ierr)
    call mpi_bcast(ymin,1,MPI_PRECISION,0,mpi_comm_world,ierr)
    call mpi_bcast(ymax,1,MPI_PRECISION,0,mpi_comm_world,ierr)
    call mpi_bcast(zmin,1,MPI_PRECISION,0,mpi_comm_world,ierr)
    call mpi_bcast(zmax,1,MPI_PRECISION,0,mpi_comm_world,ierr)

end subroutine mod_global_grid_dimensions

!-----------------------------------------------------------------------!
subroutine mod_global_grid_destroy()
    deallocate(coord_g_cg,coord_g,intma_g,bsido_g)
end subroutine mod_global_grid_destroy

!-----------------------------------------------------------------------!
subroutine mod_global_grid_out

    !local 
    integer :: ivar, ip, ie, k, mk, j, mj, i, mi
    character :: text*72
    character :: fname*72

    !Constants
    !    nelemt=nelem*(ngl-1)*(ngl-1)

    !Open GKS file
    fname = 'globalgrid.gri'
    open(1,file=fname)


    text=" NPOIN  NELEM NBOUN NOPX NOPY NOPZ "
    write(1,'(a)')text
    write(1,'(6(i7,1x))')npoin_g_cg,nelem_g,nbsido_g,nopx,nopy,nopz

    text=" COORDINATES "
    write(1,'(a)')text
    do ip=1,npoin_g_cg
        write(1, *) ip, coord_g_cg(1,ip), coord_g_cg(2,ip), coord_g_cg(3,ip)
    end do !ip 

    text=" ELEMENT CONNECTIVITY "
    write(1,'(a)')text
    do ie=1,nelem_g
        do k = 1,nglz
            mk = (k - 1)*nglx*ngly
            do j=1,ngly
                mj=(j-1)*nglx
                do i=1,nglx
                    mi=mk + mj + i
                    write(1,'(3(i7,1x))')ie,mi,intma_g(i,j,k,ie)
                end do !i
            end do !j
        end do ! k
    end do !ip   

end subroutine mod_global_grid_out

end module mod_global_grid
