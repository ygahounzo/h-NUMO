!==========================================================================================
!>@brief Apply DSS (Direct Stiffness Summation) and optionally multiply by inverse mass.
!> Note: parallel communication part not yet implemented in OO style.
!>@author Daniel S. Abdi and F.X. Giraldo
!===========================================================================================

subroutine create_global_rhs(rhs, recv_data, mvar, imass)

  implicit none

  integer, intent(in) :: mvar, imass
  real, intent(inout) :: rhs(mvar, *)
  real, intent(out)   :: recv_data(mvar, *)

end subroutine create_global_rhs
