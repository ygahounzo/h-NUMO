!----------------------------------------------------------------------!
!>@brief This module defines functions that can be used for measuring the
!> efficiency of the code by using hardware counters with the help
!> of the function library PAPI
!>@author A. Mueller
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943-5216
!----------------------------------------------------------------------!
module mod_papi
  
    use mod_mpi_utilities, only : irank, irank0

    use mpi

    implicit none
  
    public :: papi_start, papi_update, papi_print, papi_write, papi_stop
  
    private
  
  !-----------------------------------------------------------------------
  
  !module variables and parameters
#ifdef usepapi
 ! include '/home/amueller/opt/papi/5.1.0/include/f90papi.h'
  real, dimension(3) :: papi_results, papi_results_g
  integer events(3), numevents, ierr, cache_line_size
  integer(KIND=8), dimension(3) :: counts
  character*(PAPI_MAX_STR_LEN) errorstring
#endif
  
contains
  
    !-----------------------------------------------------------------------
    subroutine papi_start()
  
        implicit none
     
#ifdef usepapi
     numevents = 3
     events(1) = PAPI_TOT_CYC !total number of clock cycles
     events(2) = PAPI_FP_OPS !total number of floating point operations
     events(3) = PAPI_L2_DCM !L2 cache misses
       !get cache line size: (works on Hamming)
     open(1,file='/sys/devices/system/cpu/cpu0/cache/index1/coherency_line_size')
     read(1,*)cache_line_size
     close(1)
       !start PAPI counters
     call PAPIF_start_counters(events, numevents, ierr)
     if ( ierr .ne. PAPI_OK ) then
        call PAPIF_perror(ierr, errorstring, PAPI_MAX_STR_LEN)
        print *, 'PAPI-ERROR:'
        print *, errorstring
     endif
     papi_results = 0.0
#endif  
    end subroutine papi_start
  
    !-----------------------------------------------------------------------
    subroutine papi_update()
  
        implicit none
     
#ifdef usepapi
     counts = 0
     call PAPIF_read_counters(counts, numevents, ierr)
     if ( ierr .ne. PAPI_OK ) then
        call PAPIF_perror(ierr, errorstring, PAPI_MAX_STR_LEN)
        print *, 'PAPI-ERROR:'
        print *, errorstring
     endif
     counts(3) = counts(3)*cache_line_size
     papi_results = papi_results + real(counts)
#endif  
    end subroutine papi_update
  
    !-----------------------------------------------------------------------
    subroutine papi_print()
  
        implicit none
     
#ifdef usepapi
       !collect PAPI counters:
     call mpi_reduce(papi_results,papi_results_g,3,MPI_PRECISION,mpi_sum,0,mpi_comm_world,ierr)
     if (irank == irank0) then
        write(*,'(a)')"PAPI: cycles rank0, total cycles, flops per cycle, memory bytes per cycle"
        write(*,'(4(e26.16,1x))')papi_results(1),papi_results_g(1),papi_results_g(2)/papi_results_g(1),papi_results_g(3)/papi_results_g(1)
     end if
#endif  
    end subroutine papi_print
  
    !-----------------------------------------------------------------------
    subroutine papi_write()
  
        implicit none
     
#ifdef usepapi
     !to be written
#endif  
    end subroutine papi_write
  
    !-----------------------------------------------------------------------
    subroutine papi_stop()
  
        implicit none
     
#ifdef usepapi
     call PAPIF_stop_counters(counts, numevents, ierr)
     if ( ierr .ne. PAPI_OK ) then
        call PAPIF_perror(ierr, errorstring, PAPI_MAX_STR_LEN)
        print *, 'PAPI-ERROR:'
        print *, errorstring
     endif
#endif  
    end subroutine papi_stop
  
end module mod_papi
