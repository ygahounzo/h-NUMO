! ===========================================================================================================================
! This module contains the routines for the predictor-corrector (for baroclinic) and the RK35 time integration (for barotropic) methods
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: October 27, 2023
! ==========================================================================================================================

subroutine ti_rk3_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, q_df, qb_df)

  ! q: layer variable dp, u*dp, v*dp at quad points and their face values: q_face
  ! q_df : layer variable dp, u*dp, v*dp at nodal (dof) point
  ! qprime: value dp', u' and v' at quad points and their face values: qprime_face
  ! qb : barotopic variable pb, pb_pert = pb'*eta, ub*pb, vb*pb at quad points and their face values: qb_face
  ! qb_df: : barotopic variable pb, pb_pert = pb'*eta, ub*pb, vb*pb at nodal points
  ! qprime_df: value dp', u' and v' at nodal points
  ! qp_df_out: output variable, thickness h_k, velocity u_k,v_k, free surface ssh

  use mod_grid,              only: grid
  use mod_input,             only: input
  use mod_basis,             only: basis
  use mod_face,              only: face_CS
  use mod_parallel,          only: parallel_CS
  use mod_initial,           only: initial
  use mod_ref,               only: mref
  use mod_mpi_communicator,  only: mpi_communicator
  use mod_metrics,           only: metrics
  use mod_tensor,         only: tensor_CS
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
  type(mpi_communicator), intent(inout)    :: mpic
  type(metrics),          intent(in)    :: mt
  type(tensor_CS),        intent(in)    :: tsp
  type(bcl_CS),           intent(inout) :: bcl
  type(btp_CS),           intent(inout) :: btp

  real, dimension(4,G%npoin),             intent(inout) :: qb_df
  real, dimension(3,G%npoin,inp%nlayers), intent(inout) :: q_df

  real, dimension(4,G%npoin)             :: qbp_df
  integer :: k, ik
  real, dimension(3,G%npoin,inp%nlayers) :: q0_df, q1_df
  real, dimension(2,G%npoin,inp%nlayers) :: uv_df
  real :: dtt, a1, a2
  !$acc declare create(q0_df, q1_df, uv_df)

  ! qbp_df saves the initial qb_df on host for the ES branch reset each stage.
  qbp_df = qb_df

  ! Upload q_df and qb_df; keep them resident for the full integrator.
  !$acc data copyin(q_df, qb_df)

  ! Initialise stage arrays on device.
  !$acc kernels present(q0_df, q1_df, q_df)
  q0_df = q_df
  q1_df = q_df
  !$acc end kernels

  do ik = 1, inp%kstages_bcl

    ! Read SSPRK coefficients on CPU; they become firstprivate scalars in kernels.
    dtt = inp%dt * init%ssprk_beta_bcl(ik)
    a1  = init%ssprk_a_bcl(ik,1)
    a2  = init%ssprk_a_bcl(ik,2)

    ! GPU: compute baroclinic primed variables from the current stage state.
    ! Writes bcl%qprime_df on device — no host upload needed before create_rhs_bcl.
    call extract_qprime_df_face(G, inp, init, bcl%qprime_df, q1_df, qb_df)

    ! GPU: compute bcl/btp coefficient arrays (dpp_graduv, dpprime_visc, etc.).
    call btp_bcl_coeffs_qdf(G, inp, b, tsp, bcl, btp, bcl%qprime_df)

    if (ik == 1 .and. inp%rk_bcl_FS) then
      ! FS: BTP at stage 1 only, full N_btp.
      ! ti_barotropic_ssprk_mlswe uploads qb_df internally before use.
      call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, &
                                      qb_df, bcl%qprime_df)

    elseif (.not. inp%rk_bcl_FS) then
      ! ES: restore qb_df to the stage-0 value on host; BTP re-uploads it.
      qb_df      = qbp_df
      init%N_btp = max(1, ceiling(dtt / inp%dt_btp))
      call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, &
                                      qb_df, bcl%qprime_df)
    endif
    ! After ti_barotropic_ssprk_mlswe: qb_df host and device are both current.

    ! GPU: build baroclinic RHS; rhs_bcl downloaded to host inside but device copy valid.
    call create_rhs_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, bcl%rhs_bcl, bcl%qprime_df, q1_df)

    ! GPU SSPRK combination: q_df = a1*q0 + a2*q1 + dt*rhs.
    !$acc kernels present(q_df, q0_df, q1_df, bcl%rhs_bcl)
    do k = 1, inp%nlayers
      q_df(1,:,k) = a1*q0_df(1,:,k) + a2*q1_df(1,:,k) + dtt*bcl%rhs_bcl(1,:,k)
      q_df(2,:,k) = a1*q0_df(2,:,k) + a2*q1_df(2,:,k) + dtt*bcl%rhs_bcl(2,:,k)
      q_df(3,:,k) = a1*q0_df(3,:,k) + a2*q1_df(3,:,k) + dtt*bcl%rhs_bcl(3,:,k)
    end do
    !$acc end kernels

    ! Wall BC: potential race at shared corner nodes prevents GPU port.
    ! Download q_df, apply on CPU, re-upload.
    !$acc update host(q_df)
    call layer_mom_boundary_df(G, inp, b, mf, q_df)
    !$acc update device(q_df)

    ! GPU: extract baroclinic velocity (removes barotropic component).
    call extract_velocity(G, inp, uv_df, q_df, qb_df)

    ! GPU: reconstruct momentum from corrected velocity.
    !$acc kernels present(q_df, uv_df)
    do k = 1, inp%nlayers
      q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
      q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
    end do
    !$acc end kernels

    ! GPU: Zhang-Shu positivity limiter (element-local, no cross-element races).
    call poslimiter(b, G, inp, mt, q_df, init%alpha_mlswe)

    ! GPU: save stage result for next stage.
    !$acc kernels present(q1_df, q_df)
    q1_df = q_df
    !$acc end kernels

  end do

  ! Download final state for the caller.
  !$acc update host(q_df, qb_df)
  !$acc end data

end subroutine ti_rk3_bcl
