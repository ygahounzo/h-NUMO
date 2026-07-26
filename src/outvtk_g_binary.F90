!---------------------------------------------------------------------!
!>@brief This routine writes out the all Global Variables in VTK format
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
!>@ modified by Yao Gahounzo
!> 20 May 2024
!---------------------------------------------------------------------!

subroutine outvtk_g_binary_mlswe(G, inp, b, init, gg, par, q, qb, abs_vort, pot_vort, fname, time)

    use mod_mpi_utilities, only: irank, irank0, MPI_PRECISION
    use mod_parallel,      only: parallel_CS
    use mod_global_grid,   only: grid_global
    use mod_vtk_binary
    use mod_grid,          only: grid
    use mod_basis,         only: basis
    use mod_input,         only: input
    use mod_initial,       only: initial

    implicit none

    type(grid),        intent(in) :: G
    type(input),       intent(in) :: inp
    type(basis),       intent(in) :: b
    type(initial),     intent(in) :: init
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par

    real, intent(in) :: q(init%nvar, G%npoin), qb(inp%nvar_btp, G%npoin)
    real, intent(in) :: abs_vort(G%npoin), pot_vort(G%npoin)
    real, intent(in) :: time
    character, intent(in) :: fname*100
    logical :: has_w

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
    integer displs1(par%nproc)
    integer ierr

    real,    dimension(:,:,:), allocatable :: q_l
    real,    dimension(:,:),   allocatable :: q_g, coord_dg_gathered, qb_g
    real,    dimension(:),     allocatable :: abs_vort_g, pot_vort_g
    real,    dimension(:),     allocatable :: km
    real,    dimension(:),     allocatable :: x_uns, y_uns, z_uns
    integer, dimension(:),     allocatable :: eltype, conn
    integer, dimension(:,:),   allocatable :: conn_local, conn_g_pts
    real,    dimension(:),     allocatable :: var_uns_grid, var_uns_grid_ref

    character*72  :: cbuf
    character*12  :: output_format
    character*24  :: fnp
    integer elemType, l
    integer AllocateStatus
    integer ncells_local
    real :: xfactor, yfactor, zfactor

    has_w = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex

    ! Connectivity is built from each rank's own LOCAL grid and gathered with
    ! a point-index offset via gather_connectivity, rather than read from
    ! gg%intma_g: gg%intma_g is only ever populated by mod_global_grid_create,
    ! which nothing in this p4est-based build calls, so it is always
    ! unallocated here. Indexing it (or indexing local G%intma out to the
    ! global element count gg%nelem_g) read out-of-bounds/uninitialized
    ! memory, corrupting the CELLS block (or segfaulting).
    nglm13       = max(b%nglx-1,1)*max(b%ngly-1,1)*max(b%nglz-1,1)
    ncells       = gg%nelem_g*nglm13
    ncells_local = G%nelem*nglm13

    allocate(conn_local(b%CELL_CHILDREN, max(ncells_local,1)))

    l = 0
    do ie = 1, G%nelem
        do i = 1, max(b%nglx-1,1)
            do j = 1, max(b%ngly-1,1)
                do k = 1, max(b%nglz-1,1)
                    ii = min(i+1, b%nglx)
                    jj = min(j+1, b%ngly)
                    kk = min(k+1, b%nglz)
                    l = l + 1

                    if(b%nglx == 1) then
                        conn_local(1,l) = (G%intma( i, j, k,ie) - 1)
                        conn_local(2,l) = (G%intma( i,jj, k,ie) - 1)
                        conn_local(3,l) = (G%intma( i,jj,kk,ie) - 1)
                        conn_local(4,l) = (G%intma( i, j,kk,ie) - 1)
                    else if(b%ngly == 1) then
                        conn_local(1,l) = (G%intma( i, j, k,ie) - 1)
                        conn_local(2,l) = (G%intma(ii, j, k,ie) - 1)
                        conn_local(3,l) = (G%intma(ii, j,kk,ie) - 1)
                        conn_local(4,l) = (G%intma( i, j,kk,ie) - 1)
                    else if(b%nglz == 1) then
                        conn_local(1,l) = (G%intma( i, j, k,ie) - 1)
                        conn_local(2,l) = (G%intma(ii, j, k,ie) - 1)
                        conn_local(3,l) = (G%intma(ii,jj, k,ie) - 1)
                        conn_local(4,l) = (G%intma( i,jj, k,ie) - 1)
                    else
                        conn_local(1,l) = (G%intma( i, j, k,ie) - 1)
                        conn_local(2,l) = (G%intma(ii, j, k,ie) - 1)
                        conn_local(3,l) = (G%intma(ii,jj, k,ie) - 1)
                        conn_local(4,l) = (G%intma( i,jj, k,ie) - 1)
                        conn_local(5,l) = (G%intma( i, j,kk,ie) - 1)
                        conn_local(6,l) = (G%intma(ii, j,kk,ie) - 1)
                        conn_local(7,l) = (G%intma(ii,jj,kk,ie) - 1)
                        conn_local(8,l) = (G%intma( i,jj,kk,ie) - 1)
                    endif
                end do
            end do
        end do
    end do

    if (irank == irank0) allocate(conn_g_pts(b%CELL_CHILDREN, ncells))
    call gather_connectivity(b, par, conn_g_pts, conn_local, ncells_local, ncells)
    deallocate(conn_local)

    if (irank == irank0) then
        allocate(q_g(init%nvar, gg%npoin_g), coord_dg_gathered(3, gg%npoin_g), &
            x_uns(gg%npoin_g), y_uns(gg%npoin_g), z_uns(gg%npoin_g),           &
            eltype(ncells), conn(ncells*(b%CELL_CHILDREN+1)),                  &
            var_uns_grid(gg%npoin_g), var_uns_grid_ref(gg%npoin_g),             &
            qb_g(inp%nvar_btp, gg%npoin_g),                                    &
            abs_vort_g(gg%npoin_g), pot_vort_g(gg%npoin_g), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - OUTVTK_G_BINARY **"
    end if

    ! Gather Data onto Head node
    call gather_data(G, inp, gg, par, q_g, q, init%nvar)
    call gather_data(G, inp, gg, par, qb_g, qb, inp%nvar_btp)
    call gather_data(G, inp, gg, par, coord_dg_gathered, G%coord, 3)
    call gather_data(G, inp, gg, par, abs_vort_g, abs_vort, 1)
    call gather_data(G, inp, gg, par, pot_vort_g, pot_vort, 1)

    if (irank == irank0) then

        call vtk_ini(output_format = inp%format_vtk,   &
            filename      = fname,                     &
            title         = 'NUMO3d data',             &
            mesh_topology = 'UNSTRUCTURED_GRID',       &
            time_value    = time)

        x_uns(:) = coord_dg_gathered(1,:)
        y_uns(:) = coord_dg_gathered(2,:)
        z_uns(:) = coord_dg_gathered(3,:)

        call vtk_geo_unst_R8(gg%npoin_g, x_uns, y_uns, z_uns)

        if(b%is_2d) then
            do i = 1, ncells
                eltype(i) = 9
            end do
            nsize = 5*ncells
        else
            do i = 1, ncells
                eltype(i) = 12
            end do
            nsize = 9*ncells
        endif

        ! Build the flat VTK connectivity list (cell-size prefix + point
        ! indices per cell) from the already-gathered, correctly globally-
        ! numbered conn_g_pts.
        l = 1
        do i = 1, ncells
            conn(l) = b%CELL_CHILDREN
            conn(l+1:l+b%CELL_CHILDREN) = conn_g_pts(1:b%CELL_CHILDREN, i)
            l = l + b%CELL_CHILDREN + 1
        end do
        ncon = l-1

        call vtk_con(ncells, nsize, ncon, conn, eltype)
        call vtk_dat(gg%npoin_g, 'NODE')

        do i = 1, gg%npoin_g
            var_uns_grid(i) = q_g(1,i)
        end do
        call vtk_var_scal_R8(gg%npoin_g, 'h', var_uns_grid)

        do i = 1, gg%npoin_g
            var_uns_grid(i) = q_g(2,i)
        end do
        call vtk_var_scal_R8(gg%npoin_g, 'u', var_uns_grid)

        do i = 1, gg%npoin_g
            var_uns_grid(i) = q_g(3,i)
        end do
        call vtk_var_scal_R8(gg%npoin_g, 'v', var_uns_grid)

        if (has_w) then
            do i = 1, gg%npoin_g
                var_uns_grid(i) = q_g(4,i)
            end do
            call vtk_var_scal_R8(gg%npoin_g, 'w', var_uns_grid)
        end if

        do i = 1, gg%npoin_g
            x_uns(i) = q_g(2,i)
            y_uns(i) = q_g(3,i)
            z_uns(i) = 0.0
            if (has_w) z_uns(i) = q_g(4,i)
        end do
        call vtk_var_vect_R8('VECT', gg%npoin_g, 'MOMENTUM', x_uns, y_uns, z_uns)

        do i = 1, gg%npoin_g
            var_uns_grid(i) = qb_g(1,i)
        end do
        call vtk_var_scal_R8(gg%npoin_g, 'pb', var_uns_grid)

        do i = 1, gg%npoin_g
            var_uns_grid(i) = q_g(5,i)
        end do
        call vtk_var_scal_R8(gg%npoin_g, 'SSH', var_uns_grid)

        do i = 1, gg%npoin_g
            var_uns_grid(i) = qb_g(3,i) / qb_g(1,i)
        end do
        call vtk_var_scal_R8(gg%npoin_g, 'ub', var_uns_grid)

        do i = 1, gg%npoin_g
            var_uns_grid(i) = qb_g(4,i) / qb_g(1,i)
        end do
        call vtk_var_scal_R8(gg%npoin_g, 'vb', var_uns_grid)

        if (has_w) then
            do i = 1, gg%npoin_g
                var_uns_grid(i) = qb_g(5,i) / qb_g(1,i)
            end do
            call vtk_var_scal_R8(gg%npoin_g, 'wb', var_uns_grid)
        end if

        if (inp%lwrite_vorticity) then
            call vtk_var_scal_R8(gg%npoin_g, 'AbsVorticity', abs_vort_g)
            call vtk_var_scal_R8(gg%npoin_g, 'PotVorticity', pot_vort_g)
        end if

        call vtk_end()

        deallocate(q_g, qb_g)
        deallocate(abs_vort_g, pot_vort_g)
        deallocate(var_uns_grid, var_uns_grid_ref)
        deallocate(conn_g_pts)

    end if

end subroutine outvtk_g_binary_mlswe


subroutine outvtk_g_binary_mlswe_global(G, inp, b, gg, par, q, qb, qprime, fname, time)

    use mod_mpi_utilities, only: irank, irank0, MPI_PRECISION
    use mod_parallel,      only: parallel_CS
    use mod_global_grid,   only: grid_global
    use mod_vtk_binary
    use mod_grid,          only: grid
    use mod_basis,         only: basis
    use mod_input,         only: input

    implicit none

    type(grid),        intent(in) :: G
    type(input),       intent(in) :: inp
    type(basis),       intent(in) :: b
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par

    real, intent(in) :: q(3, G%npoin), qb(3, G%npoin), qprime(3, G%npoin)
    real, intent(in) :: time
    character, intent(in) :: fname*100

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
    integer displs1(par%nproc)
    integer ierr

    real,    dimension(:,:,:), allocatable :: q_l
    real,    dimension(:,:),   allocatable :: q_g, coord_dg_gathered, qb_g, qprime_g
    real,    dimension(:),     allocatable :: km
    real,    dimension(:),     allocatable :: x_uns, y_uns, z_uns
    integer, dimension(:),     allocatable :: eltype, conn
    real,    dimension(:),     allocatable :: var_uns_grid

    character*72  :: cbuf
    character*12  :: output_format
    character*24  :: fnp
    integer elemType, l
    integer AllocateStatus
    real :: xfactor, yfactor, zfactor
    logical :: is_cgc

    if (irank == irank0) then
        allocate(q_g(3, gg%npoin_g), coord_dg_gathered(3, gg%npoin_g),         &
            x_uns(gg%npoin_g), y_uns(gg%npoin_g), z_uns(gg%npoin_g),           &
            eltype(gg%nelem_g*max(b%nglx-1,1)*max(b%ngly-1,1)*max(b%nglz-1,1)), &
            conn(9*gg%nelem_g*max(b%nglx-1,1)*max(b%ngly-1,1)*max(b%nglz-1,1)), &
            var_uns_grid(gg%npoin_g),                                           &
            qb_g(4, gg%npoin_g), qprime_g(3, gg%npoin_g), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - OUTVTK_G_BINARY **"
    end if

    ! Gather Data onto Head node
    call gather_data(G, inp, gg, par, q_g, q, 3)
    call gather_data(G, inp, gg, par, qb_g, qb, 4)
    call gather_data(G, inp, gg, par, qprime_g, qprime, 3)
    call gather_data(G, inp, gg, par, coord_dg_gathered, G%coord, 3)

    if (irank == irank0) then

        nglm13 = max(b%nglx-1,1)*max(b%ngly-1,1)*max(b%nglz-1,1)
        ncells = gg%nelem_g*nglm13

        call vtk_ini(output_format = inp%format_vtk,   &
            filename      = fname,                     &
            title         = 'NUMO3d data',             &
            mesh_topology = 'UNSTRUCTURED_GRID',       &
            time_value    = time)

        x_uns(:) = coord_dg_gathered(1,:)
        y_uns(:) = coord_dg_gathered(2,:)
        z_uns(:) = coord_dg_gathered(3,:)

        call vtk_geo_unst_R8(gg%npoin_g, x_uns, y_uns, z_uns)

        if(b%is_2d) then
            do i = 1, ncells
                eltype(i) = 9
            end do
            nsize = 5*ncells
        else
            do i = 1, ncells
                eltype(i) = 12
            end do
            nsize = 9*ncells
        endif

        l = 1

        do ie = 1, gg%nelem_g
            do i = 1, max(b%nglx-1,1)
                do j = 1, max(b%ngly-1,1)
                    do k = 1, max(b%nglz-1,1)
                        ii = min(i+1, b%nglx)
                        jj = min(j+1, b%ngly)
                        kk = min(k+1, b%nglz)

                        if(b%nglx == 1) then
                            conn(l)   = 4
                            conn(l+1) = (G%intma( i, j, k,ie) - 1)
                            conn(l+2) = (G%intma( i,jj, k,ie) - 1)
                            conn(l+3) = (G%intma( i,jj,kk,ie) - 1)
                            conn(l+4) = (G%intma( i, j,kk,ie) - 1)
                            l = l + 5
                        else if(b%ngly == 1) then
                            conn(l)   = 4
                            conn(l+1) = (G%intma( i, j, k,ie) - 1)
                            conn(l+2) = (G%intma(ii, j, k,ie) - 1)
                            conn(l+3) = (G%intma(ii, j,kk,ie) - 1)
                            conn(l+4) = (G%intma( i, j,kk,ie) - 1)
                            l = l + 5
                        else if(b%nglz == 1) then
                            conn(l)   = 4
                            conn(l+1) = (G%intma( i, j, k,ie) - 1)
                            conn(l+2) = (G%intma(ii, j, k,ie) - 1)
                            conn(l+3) = (G%intma(ii,jj, k,ie) - 1)
                            conn(l+4) = (G%intma( i,jj, k,ie) - 1)
                            l = l + 5
                        else
                            conn(l)   = 8
                            conn(l+1) = (G%intma( i, j, k,ie) - 1)
                            conn(l+2) = (G%intma(ii, j, k,ie) - 1)
                            conn(l+3) = (G%intma(ii,jj, k,ie) - 1)
                            conn(l+4) = (G%intma( i,jj, k,ie) - 1)
                            conn(l+5) = (G%intma( i, j,kk,ie) - 1)
                            conn(l+6) = (G%intma(ii, j,kk,ie) - 1)
                            conn(l+7) = (G%intma(ii,jj,kk,ie) - 1)
                            conn(l+8) = (G%intma( i,jj,kk,ie) - 1)
                            l = l + 9
                        endif

                    end do
                end do
            end do
        end do
        ncon = l-1

        call vtk_con(ncells, nsize, ncon, conn, eltype)
        call vtk_dat(gg%npoin_g, 'NODE')

        call vtk_var_scal_R8(gg%npoin_g, 'dp_df', q_g(1,:))
        call vtk_var_scal_R8(gg%npoin_g, 'dpp_df', qprime_g(1,:))
        call vtk_var_scal_R8(gg%npoin_g, 'pbpert_df', qb_g(1,:))

        x_uns(:) = q_g(2,:)
        y_uns(:) = q_g(3,:)
        call vtk_var_vect2D_R8('VECT', gg%npoin_g, 'VELOCITY', x_uns, y_uns)

        x_uns(:) = qprime_g(2,:)
        y_uns(:) = qprime_g(3,:)
        call vtk_var_vect2D_R8('VECT', gg%npoin_g, 'VELOPRIME', x_uns, y_uns)

        x_uns(:) = qb_g(2,:)
        y_uns(:) = qb_g(3,:)
        call vtk_var_vect2D_R8('VECT', gg%npoin_g, 'VELOBARO', x_uns, y_uns)

        call vtk_end()

        deallocate(q_g, qb_g, qprime_g)
        deallocate(var_uns_grid)

    end if

end subroutine outvtk_g_binary_mlswe_global
