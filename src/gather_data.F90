subroutine gather_data(G, inp, gg, par, q_g, q, nvar)

    use mod_grid,          only: grid
    use mod_input,         only: input
    use mod_global_grid,   only: grid_global
    use mod_parallel,      only: parallel_CS
    use mpi
    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    type(grid),        intent(in) :: G
    type(input),       intent(in) :: inp
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par

    integer :: nvar
    real :: q(nvar, G%npoin)
    real :: q_g(nvar, gg%npoin_g)
    integer :: displs1(par%nproc)
    integer :: ierr, irank, iproc, ip, ip_g
    real, dimension(:,:,:), allocatable :: q_l
    integer, dimension(:), allocatable  :: recvcount

    call mpi_comm_rank(mpi_comm_world, irank, ierr)
    allocate(recvcount(1:par%nproc))

    if (irank == 0) then
        allocate(q_l(nvar, par%npoin_l_max, par%nproc))
        do iproc = 1, par%nproc
            displs1(iproc)  = nvar*(iproc - 1)*par%npoin_l_max
            recvcount(iproc) = nvar*par%npoin_l(iproc)
        end do
    else
        do iproc = 1, par%nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(q, nvar*G%npoin, MPI_PRECISION, &
        q_l, recvcount, displs1, MPI_PRECISION, 0, mpi_comm_world, ierr)

    if (irank == 0) then
        if(inp%space_method == 'cgc') then
            do iproc = 1, par%nproc
                do ip = 1, par%npoin_l(iproc)
                    ip_g = par%local_global_poin_l(ip, iproc)
                    q_g(:, ip_g) = q_l(:, ip, iproc)
                end do
            end do
        else
            ip_g = 0
            do iproc = 1, par%nproc
                do ip = 1, par%npoin_l(iproc)
                    ip_g = ip_g + 1
                    q_g(:, ip_g) = q_l(:, ip, iproc)
                end do
            end do
        endif
        deallocate(q_l)
    end if

    deallocate(recvcount)

end subroutine gather_data

subroutine gather_data_column(G, gg, par, q_g, q)

    use mod_grid,          only: grid, node_column
    use mod_global_grid,   only: grid_global
    use mod_parallel,      only: parallel_CS
    use mpi
    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    type(grid),        intent(in) :: G
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par

    real, intent(in)  :: q(G%ncol)
    real, intent(out) :: q_g(gg%ncol_g)
    integer           :: displs1(par%nproc)
    integer           :: ncol_l
    integer           :: ierr, irank, iproc, ip, ip_g

    real, dimension(:,:), allocatable :: q_l
    integer, dimension(:), allocatable :: recvcount

    call mpi_comm_rank(mpi_comm_world, irank, ierr)
    allocate(recvcount(1:par%nproc))

    if (irank == 0) then
        allocate(q_l(par%ncol_l_max, par%nproc))
        do iproc = 1, par%nproc
            displs1(iproc)   = (iproc - 1)*par%ncol_l_max
            recvcount(iproc) = par%npoin_l(iproc)/gg%nz !ncol_l(iproc) <--> npoin_l(iproc)/nz
        end do
    else
        do iproc = 1, par%nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(q, G%ncol, MPI_PRECISION, &
        q_l, recvcount, displs1, MPI_PRECISION, 0, mpi_comm_world, ierr)

    if (irank == 0) then
        do iproc = 1, par%nproc
            ncol_l = par%npoin_l(iproc)/gg%nz
            do ip = 1, ncol_l
                ip_g = node_column(G, ip, 1)
                q_g(ip_g) = q_l(ip, iproc)
            end do
        end do
        deallocate(q_l)
    end if
    deallocate(recvcount)

end subroutine gather_data_column


subroutine gather_data_column_p4est(G, gg, par, q_g, q, nvar)
    !
    ! Added by Simone marras.
    ! It's a modification of
    ! gather_data_column_v0. This version uses also nvar.
    !
    ! October 2014
    !
    ! Modified by MA Kopera to not use local_global_poin_l

    use mod_grid,          only: grid
    use mod_global_grid,   only: grid_global
    use mod_parallel,      only: parallel_CS
    use mpi
    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    type(grid),        intent(in) :: G
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par

    integer, intent(in)  :: nvar
    real,    intent(in)  :: q(nvar, G%ncol)
    real,    intent(out) :: q_g(nvar, gg%ncol_g)

    integer              :: displs1(par%nproc)
    integer              :: ierr, irank, iproc, ip, ip_g, ip_aux

    real, dimension(:,:,:), allocatable :: q_l
    integer, dimension(:), allocatable  :: recvcount

    call mpi_comm_rank(mpi_comm_world, irank, ierr)
    allocate(recvcount(1:par%nproc))

    if (irank == 0) then
        allocate(q_l(nvar, par%ncol_l_max, par%nproc))
        do iproc = 1, par%nproc
            displs1(iproc)   = nvar*(iproc - 1)*par%ncol_l_max
            recvcount(iproc) = nvar*par%ncol_l(iproc)
        end do
    else
        do iproc = 1, par%nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(q, nvar*G%ncol, MPI_PRECISION, &
        q_l, recvcount, displs1, MPI_PRECISION, 0, mpi_comm_world, ierr)

    if (irank == 0) then
        ip_g = 0
        do iproc = 1, par%nproc
            do ip = 1, par%ncol_l(iproc)
                ip_g = ip_g + 1
                q_g(:, ip_g) = q_l(:, ip, iproc)
            end do
        end do
        deallocate(q_l)
    end if

    deallocate(recvcount)

end subroutine gather_data_column_p4est


!>@author M.A.Kopera
!>Added on 01/06/2015 to gather data without local_global
subroutine gather_data_p4est(G, gg, par, q_g, q, nvar)

    use mod_grid,          only: grid
    use mod_global_grid,   only: grid_global
    use mod_parallel,      only: parallel_CS
    use mpi
    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    type(grid),        intent(in) :: G
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par

    integer :: nvar
    real :: q(nvar, G%npoin)
    real :: q_g(nvar, gg%npoin_g)
    integer :: displs1(par%nproc)
    integer :: ierr, irank, iproc, ip, ip_g
    real, dimension(:,:,:), allocatable :: q_l
    integer, dimension(:), allocatable  :: recvcount

    call mpi_comm_rank(mpi_comm_world, irank, ierr)
    allocate(recvcount(1:par%nproc))

    if (irank == 0) then
        allocate(q_l(nvar, par%npoin_l_max, par%nproc))
        do iproc = 1, par%nproc
            displs1(iproc)   = nvar*(iproc - 1)*par%npoin_l_max
            recvcount(iproc) = nvar*par%npoin_l(iproc)
        end do
    else
        do iproc = 1, par%nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(q, nvar*G%npoin, MPI_PRECISION, &
        q_l, recvcount, displs1, MPI_PRECISION, 0, mpi_comm_world, ierr)

    if (irank == 0) then
        ip_g = 0
        do iproc = 1, par%nproc
            do ip = 1, par%npoin_l(iproc)
                ip_g = ip_g + 1
                q_g(:, ip_g) = q_l(:, ip, iproc)
            end do
        end do
        deallocate(q_l)
    end if
    deallocate(recvcount)

end subroutine gather_data_p4est

!>@author M.A.Kopera
!>Added on 11/2016 to scatter data without local_global
subroutine scatter_data_p4est(G, gg, par, q_g, q, nvar)

    use mod_grid,          only: grid
    use mod_global_grid,   only: grid_global
    use mod_parallel,      only: parallel_CS
    use mpi
    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    type(grid),        intent(in) :: G
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par

    integer :: nvar
    real :: q(nvar, G%npoin)
    real :: q_g(nvar, gg%npoin_g)
    integer :: displs1(par%nproc)
    integer :: ierr, irank, iproc, ip, ip_g
    real, dimension(:,:,:), allocatable :: q_l
    integer, dimension(:), allocatable  :: sndcount

    call mpi_comm_rank(mpi_comm_world, irank, ierr)
    allocate(sndcount(1:par%nproc))

    if (irank == 0) then
        allocate(q_l(nvar, par%npoin_l_max, par%nproc))
        do iproc = 1, par%nproc
            displs1(iproc)  = nvar*(iproc - 1)*par%npoin_l_max
            sndcount(iproc) = nvar*par%npoin_l(iproc)
        end do
    else
        do iproc = 1, par%nproc
            sndcount(iproc) = 0
        end do
    end if

    if (irank == 0) then
        ip_g = 0
        do iproc = 1, par%nproc
            do ip = 1, par%npoin_l(iproc)
                ip_g = ip_g + 1
                q_l(:, ip, iproc) = q_g(:, ip_g)
            end do
        end do
    end if

    ! Scatter all the data from the head node
    call mpi_scatterv(q_l, sndcount, displs1, MPI_PRECISION, &
        q, nvar*G%npoin, MPI_PRECISION, 0, mpi_comm_world, ierr)

    if(irank == 0) deallocate(q_l)

    deallocate(sndcount)

end subroutine scatter_data_p4est


!---------------------------------------------------------------------!
!> @brief Added on 01/06/2015 to gather data without local_global
!> from a specified communicator
!> It is a bit rusty for now since it needs a lot of input
!>
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine gather_data_from(q_g, q, nvar, npoin, npoin_g, npoin_l, npoin_l_max, nproc, mpi_comm_custom)

    use mpi
    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    integer :: mpi_comm_custom
    integer :: npoin_l_max, npoin_l(nproc)
    integer :: nvar, npoin, npoin_g, nproc
    real :: q(nvar, npoin)
    real :: q_g(nvar, npoin_g)
    integer :: ierr, irank, iproc, ip, ip_g
    real, dimension(:,:,:), allocatable :: q_l
    integer, dimension(:), allocatable  :: recvcount, displs1

    call mpi_comm_rank(mpi_comm_custom, irank, ierr)

    allocate(recvcount(1:nproc), displs1(1:nproc))
    if (irank == 0) then
        allocate(q_l(nvar, npoin_l_max, nproc))
        do iproc = 1, nproc
            displs1(iproc)   = nvar*(iproc - 1)*npoin_l_max
            recvcount(iproc) = nvar*npoin_l(iproc)
        end do
    else
        do iproc = 1, nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(q, nvar*npoin, MPI_PRECISION, &
        q_l, recvcount, displs1, MPI_PRECISION, 0, mpi_comm_custom, ierr)

    if (irank == 0) then
        ip_g = 0
        do iproc = 1, nproc
            do ip = 1, npoin_l(iproc)
                ip_g = ip_g + 1
                q_g(:, ip_g) = q_l(:, ip, iproc)
            end do
        end do
        deallocate(q_l)
    end if

    deallocate(recvcount, displs1)

end subroutine gather_data_from

!---------------------------------------------------------------------!
!> @brief Added on 03/05/2015 to scatter data without local_global
!> in a specified communicator
!> It is a bit rusty for now since it needs a lot of input
!>
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine scatter_data_to(q_g, q, nvar, npoin, npoin_g, npoin_l, npoin_l_max, nproc, mpi_comm_custom)

    use mpi
    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    integer :: mpi_comm_custom
    integer :: npoin_l_max, npoin_l(nproc)
    integer :: nvar, npoin, npoin_g, nproc
    real :: q(nvar, npoin)
    real :: q_g(nvar, npoin_g)
    integer :: ierr, irank, iproc, ip, ip_g
    real, dimension(:,:,:), allocatable :: q_l
    integer, dimension(:), allocatable  :: sndcount, displs1

    call mpi_comm_rank(mpi_comm_custom, irank, ierr)

    allocate(sndcount(1:nproc), displs1(1:nproc))
    if (irank == 0) then
        allocate(q_l(nvar, npoin_l_max, nproc))
        do iproc = 1, nproc
            displs1(iproc)  = nvar*(iproc - 1)*npoin_l_max
            sndcount(iproc) = nvar*npoin_l(iproc)
        end do
    else
        do iproc = 1, nproc
            sndcount(iproc) = 0
        end do
    end if

    if (irank == 0) then
        ip_g = 0
        do iproc = 1, nproc
            do ip = 1, npoin_l(iproc)
                ip_g = ip_g + 1
                q_l(:, ip, iproc) = q_g(:, ip_g)
            end do
        end do
    end if

    ! Scatter all the data from the head node
    call mpi_scatterv(q_l, sndcount, displs1, MPI_PRECISION, &
        q, nvar*npoin, MPI_PRECISION, 0, mpi_comm_custom, ierr)

    if (irank == 0) then
        deallocate(q_l)
    end if

    deallocate(sndcount, displs1)

end subroutine scatter_data_to


!>@author M.A.Kopera
!>Added on 01/28/2015 to gather integer data without local_global using offset
!>This is used to offset connectivity while assembling data
subroutine gather_connectivity(b, par, conn_g, conn, ncells, ncells_g)

    use mod_basis,    only: basis
    use mod_parallel, only: parallel_CS
    use mpi

    implicit none

    type(basis),       intent(in) :: b
    type(parallel_CS), intent(in) :: par

    integer :: ncells, ncells_g
    integer :: conn(b%CELL_CHILDREN, ncells)
    integer :: conn_g(b%CELL_CHILDREN, ncells_g)
    integer :: displs1(par%nproc)
    integer :: ierr, irank, iproc, ip, ip_g
    integer, dimension(:,:,:), allocatable :: conn_l
    integer, dimension(:),     allocatable :: recvcount
    integer :: offset

    integer :: ncells_l_max
    integer, dimension(:), allocatable :: ncells_l

    call mpi_comm_rank(mpi_comm_world, irank, ierr)
    allocate(recvcount(1:par%nproc), ncells_l(par%nproc))

    if (irank == 0) then
        ncells_l     = par%nelem_l * b%nsubcells !number of cells per element = (ngl-1)*(ngl-1)*(ngl-1)
        ncells_l_max = par%nelem_l_max * b%nsubcells

        allocate(conn_l(b%CELL_CHILDREN, ncells_l_max, par%nproc))

        do iproc = 1, par%nproc
            displs1(iproc)   = b%CELL_CHILDREN*(iproc - 1)*ncells_l_max
            recvcount(iproc) = b%CELL_CHILDREN*ncells_l(iproc)
        end do
    else
        do iproc = 1, par%nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(conn, b%CELL_CHILDREN*ncells, mpi_integer, &
        conn_l, recvcount, displs1, mpi_integer, 0, mpi_comm_world, ierr)

    if (irank == 0) then
        ip_g   = 0
        offset = 0
        do iproc = 1, par%nproc
            do ip = 1, ncells_l(iproc)
                ip_g = ip_g + 1
                conn_g(:, ip_g) = conn_l(:, ip, iproc) + offset
            end do
            offset = offset + par%npoin_l(iproc)
        end do
        deallocate(conn_l)
    end if
    deallocate(recvcount, ncells_l)

end subroutine gather_connectivity

!---------------------------------------------------------------------!
!> @brief Added on 01/30/2015 to gather integer data without
!> local_global using offset from a specified communicator
!> This is used to offset connectivity while assembling data
!> It is a bit rusty for now since it needs a lot of input
!>
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine gather_connectivity_from(b, conn_g, conn, ncells, ncells_g, npoin_l, npoin_l_max, &
            nproc, nelem_l, nelem_l_max, mpi_comm_custom)

    use mod_basis, only: basis
    use mpi

    implicit none

    type(basis), intent(in) :: b

    integer :: ncells, ncells_g
    integer :: conn(b%CELL_CHILDREN, ncells)
    integer :: conn_g(b%CELL_CHILDREN, ncells_g)
    integer :: nproc, nelem_l_max, npoin_l_max
    integer :: nelem_l(nproc), npoin_l(nproc)
    integer :: mpi_comm_custom

    integer :: displs1(nproc)
    integer :: ierr, irank, iproc, ip, ip_g
    integer, dimension(:,:,:), allocatable :: conn_l
    integer, dimension(:),     allocatable :: recvcount
    integer :: offset

    integer :: ncells_l_max
    integer, dimension(:), allocatable :: ncells_l

    call mpi_comm_rank(mpi_comm_custom, irank, ierr)
    allocate(recvcount(1:nproc), ncells_l(nproc))

    if (irank == 0) then
        ncells_l     = nelem_l * b%nsubcells !number of cells per element = (ngl-1)*(ngl-1)*(ngl-1)
        ncells_l_max = nelem_l_max * b%nsubcells

        allocate(conn_l(b%CELL_CHILDREN, ncells_l_max, nproc))

        do iproc = 1, nproc
            displs1(iproc)   = b%CELL_CHILDREN*(iproc - 1)*ncells_l_max
            recvcount(iproc) = b%CELL_CHILDREN*ncells_l(iproc)
        end do
    else
        do iproc = 1, nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(conn, b%CELL_CHILDREN*ncells, mpi_integer, &
        conn_l, recvcount, displs1, mpi_integer, 0, mpi_comm_custom, ierr)

    if (irank == 0) then
        ip_g   = 0
        offset = 0
        do iproc = 1, nproc
            do ip = 1, ncells_l(iproc)
                ip_g = ip_g + 1
                conn_g(:, ip_g) = conn_l(:, ip, iproc) + offset
            end do
            offset = offset + npoin_l(iproc)
        end do
        deallocate(conn_l)
    end if
    deallocate(recvcount, ncells_l)

end subroutine gather_connectivity_from
