!==========================================================================================
!>@brief A better name for this routine would be APPLY_DSS. 
!>This routine communicates RHS Matrix and Possibly Multiply by Inverse mass matrix = DSS. 
!>Note: this routine does both the DSS on-processor and across processors.
!>@author Daniel S. Abdi and F.X. Giraldo
!===========================================================================================

subroutine create_global_rhs(rhs,recv_data,mvar,imass)

  use mod_basis, only: nglx, ngly, nglz

  use mod_grid, only: npoin, npoin_cg, intma, intma_cg, nelem, is_cg_coupled, EToNC

  use mod_input, only: space_method, is_non_conforming_flg

  use mod_metrics, only: massinv

  use mod_parallel, only: num_send_recv_total, nbh_send_recv

  implicit none

  !Global Variables
  integer, intent(in) :: mvar, imass
  real, intent(inout) :: rhs(mvar,npoin)
  real, intent(out)   :: recv_data(mvar, num_send_recv_total)

  !Local Variables
  integer :: ip, ip_g, ivar, ibb, i, j, k, e
  integer :: Np

  Np = nglx * ngly * nglz

  !Multiply by Inverse Mass Matrix
  if (imass == 1 .and. (is_non_conforming_flg == 0 .or. space_method .ne. 'cgd')) then
    do ivar = 1,mvar
      rhs(ivar,:) = rhs(ivar,:)*massinv(:)
    end do
  end if

end subroutine create_global_rhs

subroutine make_field_continuous(q, nflds)
  use mod_metrics, only: mass
  use mod_ref, only: recv_data
  use mod_input, only: is_non_conforming_flg, space_method
  use mod_grid, only: npoin

  implicit none

  integer, intent(in) :: nflds
  real, intent(inout) :: q(nflds,npoin)
  integer :: ivar

  if(space_method == 'cgd') then
    do ivar = 1,nflds
      q(ivar, :) = mass(:) * q(ivar, :)
    enddo
    call create_global_rhs(q,recv_data,nflds,1)
  else
    stop "Make Field Continuous Currently only works for cgd"
  endif

end subroutine make_field_continuous
