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
   use mod_splitting,        only: create_rhs_bcl
   use mod_rk_mlswe,         only: ti_barotropic_ssprk_mlswe
   use mod_barotropic_terms, only: btp_bcl_coeffs_qdf
   use mod_layer_terms,      only: extract_qprime_df_face, layer_mom_boundary_df, extract_velocity
   use mod_create_rhs_mlswe, only: layer_mass_rhs
   use mod_initial_mlswe,    only: check_layer_thickness, poslimiter

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

   ! Sphere-aware: (nvar_bcl==4), q_df(2:4,...) holds
   ! 3D Cartesian tangent-plane momentum (Px,Py,Pz); component 4 ("w" in name
   ! only).
   real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(inout) :: q_df
   real, dimension(inp%nvar_btp, G%npoin),              intent(inout) :: qb_df

   real, dimension(inp%nvar_btp, G%npoin)              :: qb_df_n
   real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers) :: qprime_df_n, qprime_df_pred, q_df_pred
   real, dimension(inp%nvar_bcl-1, G%npoin, inp%nlayers) :: uv_df
   real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers) :: q_df_temp
   real, dimension(G%npoin, inp%nlayers)    :: rhs_dp
   real, dimension(G%npoin)                 :: dp_norm_inv
   real, dimension(G%npoin)                 :: tempu, tempv
   integer :: k, I
   logical :: has_w

   ! Semi-implicit Coriolis rotation coefficients/locals (both stages, always
   ! applied — Coriolis is never added explicitly in the RHS).
   ! Cartesian (2D beta-plane): init%fdt2_bcl, init%a_bcl, init%b_bcl (mod_initial_mlswe.F90).
   real :: rl, rm, rn, omg, Px, Py, Pz, d_inv

   has_w = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

   ! ==================== Prediction step =================================
   ! Forward Euler predictor: advance all baroclinic variables by dt using
   ! the RHS evaluated at t=n.

   qb_df_n = qb_df

   call extract_qprime_df_face(G, inp, b, mt, tsp, init, bcl, bcl%qprime_df, q_df, qb_df_n)

   qprime_df_n = bcl%qprime_df

   call btp_bcl_coeffs_qdf(G, inp, b, tsp, bcl, btp, bcl%qprime_df, init%alpha_mlswe, init%pbprime_df)

   !$acc update device(bcl%qprime_df)
   call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, qb_df, bcl%qprime_df, inp%dt_btp)
   call create_rhs_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, bcl%rhs_bcl, bcl%qprime_df, q_df)

   q_df_pred = q_df + inp%dt*bcl%rhs_bcl

   ! Semi-implicit Coriolis rotation (predictor stage), applied unconditionally
   ! — Coriolis is never added explicitly in create_rhs_bcl's RHS
   if (.not. has_w) then
      ! Cartesian 2D beta-plane rotation.
      do k = 1, inp%nlayers
         tempu(:) = q_df_pred(2,:,k) + init%fdt2_bcl(:)*q_df(3,:,k)
         tempv(:) = q_df_pred(3,:,k) - init%fdt2_bcl(:)*q_df(2,:,k)
         q_df_pred(2,:,k) = init%a_bcl(:)*tempu(:) + init%b_bcl(:)*tempv(:)
         q_df_pred(3,:,k) = -init%b_bcl(:)*tempu(:) + init%a_bcl(:)*tempv(:)
      end do
   else
      ! Spherical 3D rotation, (I + omega*K_r)^-1 where K_r is the
      ! cross-product matrix for the local radial unit vector init%kvector.
      !$acc update host(q_df_pred)
      do k = 1, inp%nlayers
         do I = 1, G%npoin
            rl    = init%kvector(1,I)
            rm    = init%kvector(2,I)
            rn    = init%kvector(3,I)
            omg   = init%fdt2_bcl(I)
            Px    = q_df_pred(2,I,k) + omg*(rn*q_df(3,I,k) - rm*q_df(4,I,k))
            Py    = q_df_pred(3,I,k) + omg*(rl*q_df(4,I,k) - rn*q_df(2,I,k))
            Pz    = q_df_pred(4,I,k) + omg*(rm*q_df(2,I,k) - rl*q_df(3,I,k))
            d_inv = 1.0 / (1.0 + omg**2)
            q_df_pred(2,I,k) = (Px - omg*(rm*Pz - rn*Py)) * d_inv
            q_df_pred(3,I,k) = (Py - omg*(rn*Px - rl*Pz)) * d_inv
            q_df_pred(4,I,k) = (Pz - omg*(rl*Py - rm*Px)) * d_inv
         end do
      end do
      !$acc update device(q_df_pred)
   end if

   call layer_mom_boundary_df(G, inp, b, mf, init, q_df_pred)

   ! Floor thickness before extract_velocity to prevent division by ~1e-20.
   call poslimiter(b, G, inp, mt, q_df_pred, init%alpha_mlswe)

   ! Enforce barotropic-baroclinic velocity consistency on the predicted state.
   call extract_velocity(G, inp, b, mt, tsp, init, bcl, uv_df, q_df_pred, qb_df)
   q_df_pred(2,:,:) = uv_df(1,:,:) * q_df_pred(1,:,:)
   q_df_pred(3,:,:) = uv_df(2,:,:) * q_df_pred(1,:,:)
   if (has_w) q_df_pred(4,:,:) = uv_df(3,:,:) * q_df_pred(1,:,:)

   ! ==================== Correction step =================================
   ! Re-integrate the barotropic from t=n using the
   ! average of the t=n and predicted qprime, then update
   ! continuity and momentum separately.

   qb_df = qb_df_n

   call extract_qprime_df_face(G, inp, b, mt, tsp, init, bcl, qprime_df_pred, q_df_pred, qb_df)

   ! Centre-in-time average of the baroclinic layer-thickness forcing.
   bcl%qprime_df = 0.5*(qprime_df_pred + bcl%qprime_df)

   call btp_bcl_coeffs_qdf(G, inp, b, tsp, bcl, btp, bcl%qprime_df, init%alpha_mlswe, init%pbprime_df)
   !$acc update device(bcl%qprime_df)
   call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, qb_df, bcl%qprime_df, inp%dt_btp)

   ! Layer continuity equation
   call layer_mass_rhs(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs_dp, bcl%qprime_df)

   ! Update layer thickness
   q_df(1,:,:) = q_df(1,:,:) + inp%dt*rhs_dp

   call check_layer_thickness(b, G, inp, q_df, 'ti_2levels_bcl', 0)

   dp_norm_inv = init%pbprime_df / sum(q_df(1,:,:), dim=2)

   do k = 1, inp%nlayers
      bcl%qprime_df(1,:,k) = q_df(1,:,k) * dp_norm_inv
   end do

   ! Average of dp' for the pressure gradient in the momentum RHS.
   bcl%qprime_df(1,:,:) = 0.5*(bcl%qprime_df(1,:,:) + qprime_df_n(1,:,:))

   ! Layer momentum equation.
   ! Reuses create_rhs_bcl. Row 1
   ! (mass) of the result is discarded here — the dp update already happened via
   ! layer_mass_rhs above.
   call create_rhs_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, bcl%rhs_bcl, bcl%qprime_df, q_df)

   ! Explicit momentum update, kept in q_df_temp (not written into q_df yet):
   ! q_df(2,3,(4),:,:) still holds its time-n value and is needed below as the
   ! Coriolis rotation's cross-term source.
   q_df_temp(2,:,:) = q_df(2,:,:) + inp%dt*bcl%rhs_bcl(2,:,:)
   q_df_temp(3,:,:) = q_df(3,:,:) + inp%dt*bcl%rhs_bcl(3,:,:)
   if (has_w) q_df_temp(4,:,:) = q_df(4,:,:) + inp%dt*bcl%rhs_bcl(4,:,:)

   ! Semi-implicit Coriolis rotation (corrector stage) — same formulas as
   ! the predictor stage above, applied unconditionally.
   if (.not. has_w) then
      do k = 1, inp%nlayers
         tempu(:) = q_df_temp(2,:,k) + init%fdt2_bcl(:)*q_df(3,:,k)
         tempv(:) = q_df_temp(3,:,k) - init%fdt2_bcl(:)*q_df(2,:,k)
         q_df(2,:,k) = init%a_bcl(:)*tempu(:) + init%b_bcl(:)*tempv(:)
         q_df(3,:,k) = -init%b_bcl(:)*tempu(:) + init%a_bcl(:)*tempv(:)
      end do
   else
      !$acc update host(q_df_temp)
      do k = 1, inp%nlayers
         do I = 1, G%npoin
            rl    = init%kvector(1,I)
            rm    = init%kvector(2,I)
            rn    = init%kvector(3,I)
            omg   = init%fdt2_bcl(I)
            Px    = q_df_temp(2,I,k) + omg*(rn*q_df(3,I,k) - rm*q_df(4,I,k))
            Py    = q_df_temp(3,I,k) + omg*(rl*q_df(4,I,k) - rn*q_df(2,I,k))
            Pz    = q_df_temp(4,I,k) + omg*(rm*q_df(2,I,k) - rl*q_df(3,I,k))
            d_inv = 1.0 / (1.0 + omg**2)
            q_df(2,I,k) = (Px - omg*(rm*Pz - rn*Py)) * d_inv
            q_df(3,I,k) = (Py - omg*(rn*Px - rl*Pz)) * d_inv
            q_df(4,I,k) = (Pz - omg*(rl*Py - rm*Px)) * d_inv
         end do
      end do
      !$acc update device(q_df)
   end if

   call layer_mom_boundary_df(G, inp, b, mf, init, q_df)

   ! Floor thickness before extract_velocity to prevent division by ~1e-20.
   call poslimiter(b, G, inp, mt, q_df, init%alpha_mlswe)

   ! Enforce barotropic-baroclinic velocity consistency on the corrected state.
   call extract_velocity(G, inp, b, mt, tsp, init, bcl, uv_df, q_df, qb_df)
   q_df(2,:,:) = uv_df(1,:,:) * q_df(1,:,:)
   q_df(3,:,:) = uv_df(2,:,:) * q_df(1,:,:)
   if (has_w) q_df(4,:,:) = uv_df(3,:,:) * q_df(1,:,:)

end subroutine ti_2levels_bcl
