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

    use mod_parallel, only: num_nbh

    public :: &
        mod_mpi_communicator_create, &
        ireq, status, nreq, ierr
    private
  
    !module variables and parameters
    integer, dimension(:),   allocatable :: ireq
    integer, dimension(:,:), allocatable :: status
  
contains
  
    subroutine mod_mpi_communicator_create()
    
        implicit none
    
        integer :: AllocateStatus

        if(allocated(ireq)) deallocate(ireq, status)
        allocate( ireq(2*num_nbh), status(mpi_status_size,2*num_nbh), &
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - MOD_MPI_COMMUNICATOR_CREATE **"

    end subroutine mod_mpi_communicator_create

end module mod_mpi_communicator
