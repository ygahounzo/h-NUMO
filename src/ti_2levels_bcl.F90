! ===========================================================================================================================
! This module contains the routines for the predictor-corrector (for baroclinic) and the RK35 time integration (for barotropic) methods
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: October 27, 2023
! ==========================================================================================================================

subroutine ti_2levels_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, q_df, qb_df)

   use mod_grid,             only: grid
   use mod_basis,            only: basis
   use mod_input,            only: input
   use mod_initial,          only: initial
   use mod_metrics,          only: metrics
   use mod_face,             only: face_CS
   use mod_parallel,         only: parallel_CS
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator
   use mod_tensor,           only: tensor_CS
   use mod_variables,        only: btp_CS, bcl_CS
   use mod_splitting,        only: create_rhs_bcl, rhs_momentum
   use mod_rk_mlswe,         only: ti_barotropic_ssprk_mlswe
   use mod_barotropic_terms, only: btp_bcl_coeffs_qdf
   use mod_layer_terms,      only: extract_qprime_df_face, layer_mom_boundary_df, extract_velocity
   use mod_create_rhs_mlswe, only: layer_mass_rhs

   implicit none

   type(grid),             intent(in)    :: G
   type(input),            intent(in)    :: inp
   type(basis),            intent(in)    :: b
   type(face_CS),          intent(in)    :: mf
   type(parallel_CS),      intent(in)    :: par
   type(btp_CS),           intent(inout) :: btp
   type(bcl_CS),           intent(inout) :: bcl
   type(initial),          intent(in)    :: init
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic
   type(tensor_CS),        intent(in)    :: tsp
   type(metrics),          intent(in)    :: mt

   ! NOTE: this two-level (predictor-corrector) scheme is not yet sphere-aware —
   ! rhs_momentum/layer_momentum_rhs only ever compute u,v-momentum (rhs_mom
   ! stays 2-wide); q_df(4,...) (w) is never advanced by this integrator.
   real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(inout) :: q_df
   real, dimension(inp%nvar_btp, G%npoin),              intent(inout) :: qb_df

   real, dimension(inp%nvar_btp, G%npoin)              :: qb_df_n
   real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers) :: qprime_df_n, qprime_df_pred, q_df_pred
   real, dimension(inp%nvar_bcl-1, G%npoin, inp%nlayers) :: uv_df
   real, dimension(2, G%npoin, inp%nlayers) :: rhs_mom
   real, dimension(G%npoin, inp%nlayers)    :: rhs_dp
   real, dimension(G%npoin)                 :: dp_norm_inv
   integer :: k

   ! ==================== Prediction step =================================
   ! Forward Euler predictor: advance all baroclinic variables by dt using
   ! the RHS evaluated at t=n.

   qb_df_n = qb_df

   call extract_qprime_df_face(G, inp, init, bcl%qprime_df, q_df, qb_df_n)

   qprime_df_n = bcl%qprime_df

   call btp_bcl_coeffs_qdf(G, inp, b, tsp, bcl, btp, bcl%qprime_df, init%alpha_mlswe, init%pbprime_df)

   !$acc update device(bcl%qprime_df)
   call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, qb_df, bcl%qprime_df, inp%dt_btp)
   call create_rhs_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, bcl%rhs_bcl, bcl%qprime_df, q_df)

   q_df_pred = q_df + inp%dt*bcl%rhs_bcl

   call layer_mom_boundary_df(G, inp, b, mf, init, q_df_pred)

   ! Enforce barotropic-baroclinic velocity consistency on the predicted state.
   call extract_velocity(G, inp, uv_df, q_df_pred, qb_df)
   q_df_pred(2,:,:) = uv_df(1,:,:) * q_df_pred(1,:,:)
   q_df_pred(3,:,:) = uv_df(2,:,:) * q_df_pred(1,:,:)

   ! ==================== Correction step =================================
   ! Re-integrate the barotropic from t=n using the
   ! average of the t=n and predicted qprime, then update
   ! continuity and momentum separately.

   qb_df = qb_df_n

   call extract_qprime_df_face(G, inp, init, qprime_df_pred, q_df_pred, qb_df)

   ! Centre-in-time average of the baroclinic layer-thickness forcing.
   bcl%qprime_df = 0.5*(qprime_df_pred + bcl%qprime_df)

   call btp_bcl_coeffs_qdf(G, inp, b, tsp, bcl, btp, bcl%qprime_df, init%alpha_mlswe, init%pbprime_df)
   !$acc update device(bcl%qprime_df)
   call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, qb_df, bcl%qprime_df, inp%dt_btp)

   ! Layer continuity equation
   call layer_mass_rhs(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs_dp, bcl%qprime_df)

   ! Update layer thickness
   q_df(1,:,:) = q_df(1,:,:) + inp%dt*rhs_dp

   dp_norm_inv = init%pbprime_df / sum(q_df(1,:,:), dim=2)

   do k = 1, inp%nlayers
      bcl%qprime_df(1,:,k) = q_df(1,:,k) * dp_norm_inv
   end do

   ! Average of dp' for the pressure gradient in the momentum RHS.
   bcl%qprime_df(1,:,:) = 0.5*(bcl%qprime_df(1,:,:) + qprime_df_n(1,:,:))

   ! Layer momentum equation
   call rhs_momentum(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs_mom, bcl%qprime_df, q_df)

   ! Update momentum
   q_df(2,:,:) = q_df(2,:,:) + inp%dt*rhs_mom(1,:,:)
   q_df(3,:,:) = q_df(3,:,:) + inp%dt*rhs_mom(2,:,:)

   call layer_mom_boundary_df(G, inp, b, mf, init, q_df)

   ! Enforce barotropic-baroclinic velocity consistency on the corrected state.
   call extract_velocity(G, inp, uv_df, q_df, qb_df)
   q_df(2,:,:) = uv_df(1,:,:) * q_df(1,:,:)
   q_df(3,:,:) = uv_df(2,:,:) * q_df(1,:,:)

end subroutine ti_2levels_bcl
