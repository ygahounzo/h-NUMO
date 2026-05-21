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

    public :: read_mlswe, restart_mlswe

    contains

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
            ! call load_data_mlswe_nc(fname, qb_grp, q_grp, npoin_grp, inp%nlayers)
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
