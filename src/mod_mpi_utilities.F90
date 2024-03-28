!---------------------------------------------------------------------!
!>@brief MPI utility functions 
!---------------------------------------------------------------------!
module mod_mpi_utilities
  
    use mpi
    use mod_types, only : r8
  
    implicit none
    private
  
    !--- PUBLIC MODULE-LEVLE VARIABLES --!
    integer, parameter :: irank0=0
    integer :: irank, numproc
    !--- END PUBLIC MODULE-LEVLE VARIABLES --!
  
    !--- PRIVATE MODULE-LEVLE VARIABLES --!
    integer :: ierr_mpi
    logical :: is_mpi_initialized=.false.
    !--- END PRIVATE MODULE-LEVLE VARIABLES --!
#ifdef SINGLE  
    integer, parameter, public :: MPI_PRECISION =  MPI_REAL
#else
    integer, parameter, public :: MPI_PRECISION =  MPI_DOUBLE_PRECISION
#endif
    public :: irank0, irank, numproc
    public :: initialize_mpi_util, finalize_mpi_util, task_sync, exit_all, &
        wtime
  
contains
  
    subroutine initialize_mpi_util(ierr)
        integer, intent(out), optional :: ierr
    
        call MPI_INIT(ierr_mpi)
        if(ierr_mpi /= MPI_SUCCESS) then
            print *,'MPI_INIT returned error: ',ierr_mpi
            call exit_all(-99)
        end if
    
        is_mpi_initialized = .true.
    
        call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr_mpi)
        if(ierr_mpi /= MPI_SUCCESS) then
            print *,'MPI_COMM_RANK returned error: ',ierr_mpi
            call exit_all(-99)
        end if
    
        call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierr_mpi)
        if(ierr_mpi /= MPI_SUCCESS) then
            print *,'MPI_COMM_SIZE returned error: ',ierr_mpi
            call exit_all(-99)
        end if
    
        if(present(ierr)) ierr = ierr_mpi
    end subroutine initialize_mpi_util
  
    subroutine finalize_mpi_util(ierr)
        integer, intent(out), optional :: ierr
    
        if(.not.is_mpi_initialized) return
    
        call task_sync()
    
        call MPI_FINALIZE(ierr_mpi)
        if(ierr_mpi /= MPI_SUCCESS) then
            print *,'MPI_COMM_RANK returned error: ',ierr_mpi
            call exit_all(-99)
        end if
    
        ierr=0
    
    end subroutine finalize_mpi_util
  
    subroutine task_sync()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)
        if(ierr_mpi /= MPI_SUCCESS) then
            print *,'MPI_BARRIER returned error: ',ierr_mpi
            call exit_all(-99)
        end if
    end subroutine task_sync
  
    function wtime()
        real(kind=r8) wtime
        call MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)
        wtime=MPI_WTIME()
    end function wtime
  
    subroutine exit_all(exit_code)
        integer, intent(in) :: exit_code
    
        if(.not.is_mpi_initialized) then
            call exit(-99)
        end if
    
        call MPI_ABORT(MPI_COMM_WORLD, exit_code, ierr_mpi)
    end subroutine exit_all
  
end module mod_mpi_utilities
