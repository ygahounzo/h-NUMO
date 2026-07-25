!---------------------------------------------------------------------!
!> @brief This module contains subroutines that read output files to restart simulation
!>
!> @author Yao Gahounzo on 01/2024
!---------------------------------------------------------------------!

module mod_restart

    use mod_grid,    only: grid
    use mod_basis,   only: basis
    use mod_input,   only: input
    use mod_initial, only: initial

    implicit none

    public :: read_mlswe, restart_mlswe, write_restart_nc

    contains

!---------------------------------------------------------------------!
!> @brief Writes the full prognostic state needed to restart the run
!> to a numbered NetCDF file under RESTART/ (created on rank 0 if
!> missing). Called periodically (every inp%write_restart_time) and
!> once, unconditionally, at the end of the run.
!---------------------------------------------------------------------!
subroutine write_restart_nc(G, inp, gg, par, init, q_df, qb, irestart_nc, time, itime)

    use mod_grid,          only: grid
    use mod_input,         only: input
    use mod_global_grid,   only: grid_global
    use mod_parallel,      only: parallel_CS
    use mod_initial,       only: initial
    use mod_constants,     only: gravity
    use mod_mpi_utilities, only: irank, irank0
    use netcdf

    implicit none

    type(grid),        intent(in) :: G
    type(input),       intent(in) :: inp
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par
    type(initial),     intent(in) :: init

    real, intent(in)    :: q_df(inp%nvar_bcl, G%npoin, inp%nlayers)
    real, intent(in)    :: qb(inp%nvar_btp, G%npoin)
    integer, intent(in) :: irestart_nc
    real, intent(in)    :: time
    integer, intent(in) :: itime

    character(len=*), parameter :: restart_dir = 'RESTART'
    character(len=*), parameter :: restart_log = 'RESTART/hNUMO.res'
    character(len=25) :: restart_fname
    character(len=25) :: hdr_filename = 'FILENAME'
    logical :: log_exists
    integer :: log_unit

    character(len=200) :: fn
    character(len=4)   :: num
    logical :: has_w
    integer :: j, k, iloop, ncid, dimids_2d(2)
    real :: q(5, G%npoin, inp%nlayers)
    real :: q_gg(5, gg%npoin_g, inp%nlayers), ql(5, G%npoin), zbot_g(gg%npoin_g)
    real :: coord_dg_gathered(3, gg%npoin_g), qb_g(inp%nvar_btp, gg%npoin_g), q_g(5, gg%npoin_g)
    real, dimension(G%npoin,    inp%nlayers+1) :: mslwe_elevation
    real, dimension(gg%npoin_g, inp%nlayers+1) :: eta

    character (len = *), parameter :: TIME_NAME      = "time"
    character (len = *), parameter :: NK_NAME        = "nlayers"
    character (len = *), parameter :: NPOIN_NAME     = "npoin"
    character (len = *), parameter :: DT_NAME        = "dt"
    character (len = *), parameter :: DT_BTP_NAME    = "dt_btp"
    character (len = *), parameter :: X_NAME         = "x"
    character (len = *), parameter :: Y_NAME         = "y"
    character (len = *), parameter :: PB_NAME        = "pb"
    character (len = *), parameter :: PBUB_NAME      = "pbub"
    character (len = *), parameter :: PBVB_NAME      = "pbvb"
    character (len = *), parameter :: PBWB_NAME      = "pbwb"
    character (len = *), parameter :: H_NAME         = "h"
    character (len = *), parameter :: U_NAME         = "u"
    character (len = *), parameter :: V_NAME         = "v"
    character (len = *), parameter :: W_NAME         = "w"
    character (len = *), parameter :: ETA_NAME       = "eta"
    character (len = *), parameter :: INTERFACE_NAME = "zi"

    integer :: time_dimid, npoin_dimid, nlayers_dimid, zi_dimid
    integer :: dt_varid, dt_btp_varid, pb_varid
    integer :: pbub_varid, pbvb_varid, pbwb_varid, h_varid, u_varid, v_varid, w_varid, e_varid
    integer :: x_varid, y_varid

    has_w = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

    do k = 1, inp%nlayers
        q(1,:,k) = (init%alpha_mlswe(k)/gravity)*q_df(1,:,k)
        q(2,:,k) = q_df(2,:,k) / q_df(1,:,k)
        q(3,:,k) = q_df(3,:,k) / q_df(1,:,k)
        if (has_w) then
            q(4,:,k) = q_df(4,:,k) / q_df(1,:,k)
        else
            q(4,:,k) = q_df(1,:,k)
        end if
    end do

    mslwe_elevation = 0.0
    mslwe_elevation(:, inp%nlayers+1) = init%zbot_df

    do k = inp%nlayers, 1, -1
        mslwe_elevation(:,k) = mslwe_elevation(:,k+1) + q(1,:,k)
    end do

    q(5,:,1) = mslwe_elevation(:,1)
    q(5,:,2:inp%nlayers) = mslwe_elevation(:,2:inp%nlayers)

    ! Gather Data onto Head node
    do k = 1, inp%nlayers
        ql = q(:,:,k)
        call gather_data(G, inp, gg, par, q_g, ql, 5)
        q_gg(:,:,k) = q_g(:,:)
    enddo

    call gather_data(G, inp, gg, par, qb_g,               qb,             inp%nvar_btp)
    call gather_data(G, inp, gg, par, coord_dg_gathered,   G%coord,        3)
    call gather_data(G, inp, gg, par, zbot_g,              init%zbot_df,   1)

    eta(:, 1:inp%nlayers) = q_gg(5,:,1:inp%nlayers)
    eta(:, inp%nlayers+1) = zbot_g

    if (irank == irank0) then

        call execute_command_line('mkdir -p ' // restart_dir)

        write(num,'(i4)') irestart_nc
        iloop = 3 - int(log10(real(irestart_nc)))
        do j = 1, iloop
            num(j:j) = '0'
        end do

        restart_fname = 'restart_mlswe_' // num // '.nc'
        fn = trim(restart_dir) // '/' // trim(restart_fname)

        ! Create NetCDF file
        call check(nf90_create(trim(fn), nf90_clobber+nf90_64bit_offset, ncid))

        ! Define dimensions
        call check(nf90_def_dim(ncid, TIME_NAME,      NF90_UNLIMITED,   time_dimid))
        call check(nf90_def_dim(ncid, NPOIN_NAME,     gg%npoin_g,       npoin_dimid))
        call check(nf90_def_dim(ncid, NK_NAME,        inp%nlayers,      nlayers_dimid))
        call check(nf90_def_dim(ncid, INTERFACE_NAME, inp%nlayers+1,    zi_dimid))

        call check(nf90_put_att(ncid, NF90_GLOBAL, "filename", trim(fn)))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "npoin", "Number of points in the mesh"))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "zi", "Number of interfaces"))

        ! Name each variable
        call check(nf90_def_var(ncid, DT_NAME,     NF90_DOUBLE, 1,          dt_varid))
        call check(nf90_put_att(ncid, dt_varid,     "name", "Baroclinic time step"))
        call check(nf90_put_att(ncid, dt_varid,     "units", "seconds"))
        call check(nf90_def_var(ncid, DT_BTP_NAME, NF90_DOUBLE, 1,          dt_btp_varid))
        call check(nf90_put_att(ncid, dt_btp_varid, "name", "Barotropic time step"))
        call check(nf90_put_att(ncid, dt_btp_varid, "units", "seconds"))
        call check(nf90_def_var(ncid, X_NAME, NF90_DOUBLE, npoin_dimid, x_varid))
        call check(nf90_put_att(ncid, x_varid, "name", "cartesian coordinates"))
        call check(nf90_put_att(ncid, x_varid, "axis", "X"))
        call check(nf90_def_var(ncid, Y_NAME, NF90_DOUBLE, npoin_dimid, y_varid))
        call check(nf90_put_att(ncid, y_varid, "name", "cartesian coordinates"))
        call check(nf90_put_att(ncid, y_varid, "axis", "Y"))
        call check(nf90_def_var(ncid, PB_NAME,   NF90_DOUBLE, npoin_dimid, pb_varid))
        call check(nf90_put_att(ncid, pb_varid,   "name", "Barotropic pressure pb"))
        call check(nf90_put_att(ncid, pb_varid,   "units", "N/m^2"))
        call check(nf90_def_var(ncid, PBUB_NAME, NF90_DOUBLE, npoin_dimid, pbub_varid))
        call check(nf90_put_att(ncid, pbub_varid, "name", "Barotropic u-momentum"))
        call check(nf90_put_att(ncid, pbub_varid, "units", "kg.m/s"))
        call check(nf90_def_var(ncid, PBVB_NAME, NF90_DOUBLE, npoin_dimid, pbvb_varid))
        call check(nf90_put_att(ncid, pbvb_varid, "name", "Barotropic v-momentum"))
        call check(nf90_put_att(ncid, pbvb_varid, "units", "kg.m/s"))
        if (has_w) then
            call check(nf90_def_var(ncid, PBWB_NAME, NF90_DOUBLE, npoin_dimid, pbwb_varid))
            call check(nf90_put_att(ncid, pbwb_varid, "name", "Barotropic w-momentum"))
            call check(nf90_put_att(ncid, pbwb_varid, "units", "kg.m/s"))
        end if
        dimids_2d = (/npoin_dimid, nlayers_dimid/)
        call check(nf90_def_var(ncid, H_NAME, NF90_DOUBLE, dimids_2d, h_varid))
        call check(nf90_put_att(ncid, h_varid, "name", "Layer thickness"))
        call check(nf90_put_att(ncid, h_varid, "units", "m"))
        call check(nf90_def_var(ncid, U_NAME, NF90_DOUBLE, dimids_2d, u_varid))
        call check(nf90_put_att(ncid, u_varid, "name", "Baroclinic u-velocity"))
        call check(nf90_put_att(ncid, u_varid, "units", "m/s"))
        call check(nf90_def_var(ncid, V_NAME, NF90_DOUBLE, dimids_2d, v_varid))
        call check(nf90_put_att(ncid, v_varid, "name", "Baroclinic v-velocity"))
        call check(nf90_put_att(ncid, v_varid, "units", "m/s"))
        if (has_w) then
            call check(nf90_def_var(ncid, W_NAME, NF90_DOUBLE, dimids_2d, w_varid))
            call check(nf90_put_att(ncid, w_varid, "name", "Baroclinic w-velocity"))
            call check(nf90_put_att(ncid, w_varid, "units", "m/s"))
        end if
        dimids_2d = (/npoin_dimid, zi_dimid/)
        call check(nf90_def_var(ncid, ETA_NAME, NF90_DOUBLE, dimids_2d, e_varid))
        call check(nf90_put_att(ncid, e_varid, "name", "Interface Height Relative to Mean Sea Level"))
        call check(nf90_put_att(ncid, e_varid, "units", "m"))

        ! End define mode
        call check(nf90_enddef(ncid))

        ! Write data
        call check(nf90_put_var(ncid, dt_varid,     inp%dt))
        call check(nf90_put_var(ncid, dt_btp_varid, inp%dt_btp))
        call check(nf90_put_var(ncid, x_varid,    coord_dg_gathered(1,:)))
        call check(nf90_put_var(ncid, y_varid,    coord_dg_gathered(2,:)))
        call check(nf90_put_var(ncid, pb_varid,   qb_g(1,:)))
        call check(nf90_put_var(ncid, pbub_varid, qb_g(3,:)))
        call check(nf90_put_var(ncid, pbvb_varid, qb_g(4,:)))
        if (has_w) call check(nf90_put_var(ncid, pbwb_varid, qb_g(5,:)))
        call check(nf90_put_var(ncid, h_varid, q_gg(1,:,:)))
        call check(nf90_put_var(ncid, u_varid, q_gg(2,:,:)))
        call check(nf90_put_var(ncid, v_varid, q_gg(3,:,:)))
        if (has_w) call check(nf90_put_var(ncid, w_varid, q_gg(4,:,:)))
        call check(nf90_put_var(ncid, e_varid, eta(:,:)))

        ! Close NetCDF file
        call check(nf90_close(ncid))

        print *, "Restart written: ", trim(fn)

        ! Log the restart index against the simulation time, so the user
        ! knows what TIME_INITIAL to set in numo3d.in to resume from it.
        inquire(file=restart_log, exist=log_exists)
        if (.not. log_exists) then
            open(newunit=log_unit, file=restart_log, status='new', action='write')
            write(log_unit,'(a)') '# h-NUMO restart index'
            write(log_unit,'(a)') '# To resume: set TIME_INITIAL to the value below and IRESTART_FILE_NUMBER to INDEX in numo3d.in'
            write(log_unit,'(a)') '# (RESTART/restart_mlswe_####.nc is picked up automatically regardless of OUT_TYPE).'
            write(log_unit,'(a)')
            write(log_unit,'(a6,2x,a25,1x,a24,2x,a10)') 'INDEX', hdr_filename, 'TIME_INITIAL[time_scale]', 'ITIME'
        else
            open(newunit=log_unit, file=restart_log, status='old', action='write', position='append')
        end if

        write(log_unit,'(i6,2x,a25,1x,f24.6,2x,i10)') irestart_nc, restart_fname, time/inp%time_scale, itime

        close(log_unit)

    end if

    contains

    !-----------------------!
    !>@brief Check status
    !-----------------------!
    subroutine check(status)
        integer, intent(in) :: status

        if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check

end subroutine write_restart_nc

subroutine restart_mlswe(G, inp, b, init, q_df, qb_df, qp_df_out, fname)

    use mod_constants, only: gravity

    implicit none

    type(grid),    intent(in) :: G
    type(input),   intent(in) :: inp
    type(basis),   intent(in) :: b
    type(initial), intent(in) :: init

    character(len=*), intent(in) :: fname

    real, dimension(4, G%npoin),              intent(out) :: qb_df
    real, dimension(3, G%npoin, inp%nlayers), intent(out) :: q_df
    real, dimension(5, G%npoin, inp%nlayers), intent(out) :: qp_df_out

    real, dimension(G%npoin, inp%nlayers+1) :: mslwe_elevation
    real, dimension(3, G%npoin)             :: qb_df_read
    real, dimension(3, G%npoin, inp%nlayers) :: q_df_read

    integer :: k

    call read_mlswe(G, inp, b, q_df_read, qb_df_read, fname)

    ! Barotropic variables at the dofs (nodal points)
    qb_df(1,:)   = qb_df_read(1,:)
    qb_df(3:4,:) = qb_df_read(2:3,:)
    qb_df(2,:)   = qb_df(1,:) - init%pbprime_df(:)

    ! Baroclinic variables
    do k = 1, inp%nlayers
        q_df(1,:,k) = (gravity / init%alpha_mlswe(k)) * q_df_read(1,:,k)
        q_df(2,:,k) = q_df_read(2,:,k) * q_df(1,:,k)
        q_df(3,:,k) = q_df_read(3,:,k) * q_df(1,:,k)
    end do

    ! Prepare output variables
    do k = 1, inp%nlayers
        qp_df_out(1,:,k) = (init%alpha_mlswe(k) / gravity) * q_df(1,:,k)
        qp_df_out(2,:,k) = q_df(2,:,k) / q_df(1,:,k)
        qp_df_out(3,:,k) = q_df(3,:,k) / q_df(1,:,k)
    end do

    mslwe_elevation(:, inp%nlayers+1) = init%zbot_df

    do k = inp%nlayers, 1, -1
        mslwe_elevation(:,k) = mslwe_elevation(:,k+1) + qp_df_out(1,:,k)
    end do

    qp_df_out(4,:,1) = qb_df(3,:)
    qp_df_out(4,:,2) = qb_df(3,:)

    qp_df_out(5,:,1) = mslwe_elevation(:,1)
    qp_df_out(5,:,2) = mslwe_elevation(:,2)

end subroutine restart_mlswe

subroutine read_mlswe(G, inp, b, q_df, qb_df, fname)

    use mpi

    implicit none

    type(grid),  intent(in) :: G
    type(input), intent(in) :: inp
    type(basis), intent(in) :: b

    real :: q_df(3, G%npoin, inp%nlayers), qb_df(3, G%npoin)
    character fname*100, fname_proc*200, proc_num*10
    character grid_file*100

    real, dimension(:,:),   allocatable :: qb_grp
    real, dimension(:,:,:), allocatable :: q_grp

    integer :: i, j, k
    integer :: ncells, nproc_grp
    integer :: ivar, ifactor, ip, ic, ie, e, ndof, im, mp, me
    integer :: irank, ierr, nproc

    integer :: out_group, ngroups
    integer :: igroup, ig_rank, buffer
    integer :: npoin_grp, nelem_grp, ncells_grp

    integer, dimension(:), allocatable :: npoin_l_grp, nelem_l_grp
    integer :: npoin_l_max_grp, nelem_l_max_grp, pos

    ngroups = 1

    call mpi_comm_rank(mpi_comm_world, irank, ierr)
    call mpi_comm_size(mpi_comm_world, nproc, ierr)

    call set_groups(out_group, ngroups, igroup, ig_rank)
    call mpi_comm_size(out_group, nproc_grp, ierr)

    ncells = G%nelem * b%nsubcells

    if(ig_rank == 0) allocate(npoin_l_grp(nproc_grp), nelem_l_grp(nproc_grp))

    call mpi_gather(G%npoin, 1, mpi_integer, npoin_l_grp, 1, mpi_integer, 0, out_group, ierr)
    call mpi_gather(G%nelem, 1, mpi_integer, nelem_l_grp, 1, mpi_integer, 0, out_group, ierr)

    if(ig_rank == 0) then

        npoin_grp = sum(npoin_l_grp)
        nelem_grp = sum(nelem_l_grp)

        npoin_l_max_grp = maxval(npoin_l_grp)
        nelem_l_max_grp = maxval(nelem_l_grp)

        ncells_grp = nelem_grp * b%nsubcells
        allocate(q_grp(3, npoin_grp, inp%nlayers))
        allocate(qb_grp(3, npoin_grp))

        pos = index(fname, '.nc')

        if (pos > 0 .and. pos + 2 == len_trim(fname)) then
            call load_data_mlswe_nc(fname, qb_grp, q_grp, npoin_grp, inp%nlayers)
        else
            call load_data_mlswe(fname, qb_grp, q_grp, npoin_grp, inp%nlayers)
        end if

    end if

    call scatter_data_to(qb_grp, qb_df, 3, G%npoin, npoin_grp, npoin_l_grp, &
        npoin_l_max_grp, nproc_grp, out_group)

    do k = 1, inp%nlayers
        call scatter_data_to(q_grp(:,:,k), q_df(:,:,k), 3, G%npoin, npoin_grp, &
            npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    end do

    if(ig_rank == 0) then
        deallocate(q_grp, qb_grp)
        deallocate(npoin_l_grp, nelem_l_grp)
    end if

end subroutine read_mlswe

subroutine load_data_mlswe(name_data_file, qb_df_g, q_df_g, npoin_g, nlayers)

    character(len=*), intent(in) :: name_data_file
    integer, intent(in) :: npoin_g, nlayers

    real, dimension(3, npoin_g),          intent(out) :: qb_df_g
    real, dimension(3, npoin_g, nlayers), intent(out) :: q_df_g

    real, dimension(2, npoin_g)       :: coord_g
    real, dimension(npoin_g, nlayers+1) :: z_g

    integer :: count, nk, npoin_gg, i, j
    real :: dt, dt_btp

    open(unit=10, file=name_data_file, status='old', action='read')

    read(10, *) nk
    read(10, *) npoin_gg
    read(10, *) dt
    read(10, *) dt_btp

    read(10, *) (coord_g(:, i), i = 1, npoin_gg)

    do j = 1, 3
        read(10, *) (qb_df_g(j, i), i = 1, npoin_gg)
    end do

    read(10, *) (q_df_g(1, :, j), j = 1, nk)
    read(10, *) (q_df_g(2, :, j), j = 1, nk)
    read(10, *) (q_df_g(3, :, j), j = 1, nk)

    read(10, *) (z_g(:, j), j = 1, nk + 1)

    close(10)
end subroutine load_data_mlswe

!---------------------------------------------------------------------!
!> @brief NetCDF counterpart of load_data_mlswe: reads the pb/pbub/pbvb
!> and per-layer h/u/v fields written by write_restart_nc back into the
!> same gathered layout load_data_mlswe produces (w is not needed to
!> restart the run and is not read back, matching the txt path).
!---------------------------------------------------------------------!
subroutine load_data_mlswe_nc(name_data_file, qb_df_g, q_df_g, npoin_g, nlayers)

    use netcdf

    implicit none

    character(len=*), intent(in) :: name_data_file
    integer, intent(in) :: npoin_g, nlayers

    real, dimension(3, npoin_g),          intent(out) :: qb_df_g
    real, dimension(3, npoin_g, nlayers), intent(out) :: q_df_g

    integer :: ncid, varid

    call check(nf90_open(trim(name_data_file), nf90_nowrite, ncid))

    call check(nf90_inq_varid(ncid, 'pb',   varid))
    call check(nf90_get_var(ncid, varid, qb_df_g(1,:)))
    call check(nf90_inq_varid(ncid, 'pbub', varid))
    call check(nf90_get_var(ncid, varid, qb_df_g(2,:)))
    call check(nf90_inq_varid(ncid, 'pbvb', varid))
    call check(nf90_get_var(ncid, varid, qb_df_g(3,:)))

    call check(nf90_inq_varid(ncid, 'h', varid))
    call check(nf90_get_var(ncid, varid, q_df_g(1,:,:)))
    call check(nf90_inq_varid(ncid, 'u', varid))
    call check(nf90_get_var(ncid, varid, q_df_g(2,:,:)))
    call check(nf90_inq_varid(ncid, 'v', varid))
    call check(nf90_get_var(ncid, varid, q_df_g(3,:,:)))

    call check(nf90_close(ncid))

    contains

    subroutine check(status)
        integer, intent(in) :: status

        if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check

end subroutine load_data_mlswe_nc

!---------------------------------------------------------------------!
!> @brief This subroutine creates communicators for ngroups of
!> processors (group_comm) and returns group number (igroup) and
!> rank within the group (ig_rank)
!>
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine set_groups(group_comm, ngroups, igroup, ig_rank)

    use mpi

    implicit none

    integer :: group_comm, ngroups

    integer :: irank, ierr, i
    integer :: igroup, ig_rank, group_size, nproc

    call mpi_comm_rank(mpi_comm_world, irank, ierr)
    call mpi_comm_size(mpi_comm_world, nproc, ierr)

    group_size = ceiling(nproc / real(ngroups))

    if(group_size*(ngroups-1) == nproc) ngroups = ngroups - 1

    igroup  = int(irank / group_size)
    ig_rank = mod(irank, group_size)

    call mpi_comm_split(mpi_comm_world, igroup, irank, group_comm, ierr)

end subroutine set_groups

end module mod_restart
