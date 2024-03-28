!-------------------------------
!> @brief Applies filter
!-------------------------------
module mod_filter
  
    use mod_grid, only: npoin

    use mod_initial, only: nvar

    use mod_parallel, only: num_send_recv_total
  
    public :: &
        mod_filter_create_rhs, &
        b, &
        b_data
  
    private

    !-----------------------------------------------------------------------
    real, dimension(:,:),   allocatable :: b
    real, dimension(:,:),   allocatable :: b_data
  !-----------------------------------------------------------------------
  
contains
  
    !-----------------------------------------------------------------------
    subroutine mod_filter_create_rhs()
    
        implicit none
    
        integer AllocateStatus
   
        if(allocated(b)) deallocate(b, b_data)
        allocate( b(nvar,npoin), b_data(nvar,num_send_recv_total), stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - MOD_FILTER_CREATE_RHS **"

    end subroutine mod_filter_create_rhs
  
end module mod_filter
