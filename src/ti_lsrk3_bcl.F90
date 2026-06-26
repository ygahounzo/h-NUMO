! ===========================================================================================================================
! LSRK3 (Wicker & Skamarock, 2002) time integration for the baroclinic (BCL) mode:
!   phi^(n+1/3) = phi^n + (dt/3) * F(phi^n)
!   phi^(n+1/2) = phi^n + (dt/2) * F(phi^(n+1/3))
!   phi^(n+1)   = phi^n + dt     * F(phi^(n+1/2))
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
! ==========================================================================================================================

subroutine ti_lsrk3_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, q_df, qb_df)

  use mod_grid,              only: grid
  use mod_input,             only: input
  use mod_basis,             only: basis
  use mod_face,              only: face_CS
  use mod_parallel,          only: parallel_CS
  use mod_initial,           only: initial
  use mod_ref,               only: mref
  use mod_mpi_communicator,  only: mpi_communicator
  use mod_metrics,           only: metrics
  use mod_tensor,            only: tensor_CS
  use mod_variables,         only: bcl_CS, btp_CS
  use mod_splitting,         only: create_rhs_bcl
  use mod_rk_mlswe,          only: ti_barotropic_ssprk_mlswe
  use mod_barotropic_terms,  only: btp_bcl_coeffs_qdf
  use mod_layer_terms,       only: extract_qprime_df_face, layer_mom_boundary_df, extract_velocity
  use mod_initial_mlswe,    only: poslimiter

  implicit none

  type(grid),             intent(in)    :: G
  type(input),            intent(in)    :: inp
  type(basis),            intent(in)    :: b
  type(face_CS),          intent(in)    :: mf
  type(parallel_CS),      intent(in)    :: par
  type(initial),          intent(inout) :: init
  type(mref),             intent(inout) :: ref
  type(mpi_communicator), intent(inout) :: mpic
  type(metrics),          intent(in)    :: mt
  type(tensor_CS),        intent(in)    :: tsp
  type(bcl_CS),           intent(inout) :: bcl
  type(btp_CS),           intent(inout) :: btp

  real, dimension(4,G%npoin),             intent(inout) :: qb_df
  real, dimension(3,G%npoin,inp%nlayers), intent(inout) :: q_df

  integer :: k, ik
  real :: dtt
  real, parameter :: lsrk3_beta(3) = (/ 1.0/3.0, 1.0/2.0, 1.0 /)

  ! Save stage-0 qb_df on GPU for the ES branch reset each stage.
  !$acc kernels present(bcl%qbp_df, qb_df)
  bcl%qbp_df = qb_df
  !$acc end kernels

  ! Initialise stage array on device.
  !$acc kernels present(bcl%q0_df, q_df)
  bcl%q0_df = q_df
  !$acc end kernels

  do ik = 1, 3

    dtt = inp%dt * lsrk3_beta(ik)

    ! GPU: compute baroclinic primed variables from the current stage state.
    call extract_qprime_df_face(G, inp, init, bcl%qprime_df, q_df, qb_df)
    ! GPU: compute bcl/btp coefficient arrays (fused single kernel).
    call btp_bcl_coeffs_qdf(G, inp, b, tsp, bcl, btp, bcl%qprime_df)

    if (ik == 1 .and. inp%rk_bcl_FS) then
      ! FS: BTP at stage 1 only, full N_btp.
      call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, &
                                      qb_df, bcl%qprime_df)

    elseif (.not. inp%rk_bcl_FS) then
      ! ES: restore qb_df to the stage-0 value on GPU.
      !$acc kernels present(qb_df, bcl%qbp_df)
      qb_df      = bcl%qbp_df
      !$acc end kernels
      init%N_btp = max(1, ceiling(dtt / inp%dt_btp))
      call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, &
                                      qb_df, bcl%qprime_df)
    endif

    ! GPU: build baroclinic RHS; rhs_bcl downloaded to host inside but device copy valid.
    call create_rhs_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, bcl%rhs_bcl, bcl%qprime_df, q_df)

    ! GPU LSRK3 update: q = q0 + dt * rhs.
    !$acc kernels present(q_df, bcl%q0_df, bcl%rhs_bcl)
    do k = 1, inp%nlayers
      q_df(1,:,k) = bcl%q0_df(1,:,k) + dtt * bcl%rhs_bcl(1,:,k)
      q_df(2,:,k) = bcl%q0_df(2,:,k) + dtt * bcl%rhs_bcl(2,:,k)
      q_df(3,:,k) = bcl%q0_df(3,:,k) + dtt * bcl%rhs_bcl(3,:,k)
    end do
    !$acc end kernels

    ! GPU: wall BC — gang over faces, atomic updates for corner nodes.
    call layer_mom_boundary_df(G, inp, b, mf, q_df)

    ! GPU: extract baroclinic velocity (removes barotropic component).
    call extract_velocity(G, inp, bcl%uv_df, q_df, qb_df)

    ! GPU: reconstruct momentum from corrected velocity.
    !$acc kernels present(q_df, bcl%uv_df)
    do k = 1, inp%nlayers
      q_df(2,:,k) = bcl%uv_df(1,:,k) * q_df(1,:,k)
      q_df(3,:,k) = bcl%uv_df(2,:,k) * q_df(1,:,k)
    end do
    !$acc end kernels

    ! GPU: Zhang-Shu positivity limiter (element-local, no cross-element races).
    call poslimiter(b, G, inp, mt, q_df, init%alpha_mlswe)

  end do

end subroutine ti_lsrk3_bcl
