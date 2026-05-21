!----------------------------------------------------------------------!
!>@brief This module stores the MPI COMMUNICATOR arrays
!>@author Francis X. Giraldo
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@date Modified May 2, 2018
!----------------------------------------------------------------------!
module mod_mpi_communicator
  
    use mpi
    use mod_parallel

    type :: mpi_communicator
        integer, dimension(:),   allocatable :: ireq
        integer, dimension(:,:), allocatable :: status
        integer :: nreq, ierr
    end type mpi_communicator

    public :: &
        mod_mpi_communicator_create, mpi_communicator
    private
  
    !module variables and parameters
    integer, dimension(:),   allocatable :: ireq
    integer, dimension(:,:), allocatable :: status
  
contains
  
    subroutine mod_mpi_communicator_create(par, mpic)
    
        implicit none

        type(parallel_CS), intent(in) :: par
        type(mpi_communicator), intent(inout) :: mpic
    
        integer :: AllocateStatus

        if(allocated(mpic%ireq)) deallocate(mpic%ireq, mpic%status)
        allocate( mpic%ireq(2*par%num_nbh), mpic%status(mpi_status_size,2*par%num_nbh), &
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - MOD_MPI_COMMUNICATOR_CREATE **"

    end subroutine mod_mpi_communicator_create

end module mod_mpi_communicator
