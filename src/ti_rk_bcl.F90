! ===========================================================================================================================
! This module contains the routines for the predictor-corrector (for baroclinic) and the RK35 time integration (for barotropic) methods
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: October 27, 2023
! ==========================================================================================================================

subroutine ti_rk_bcl(q_df, qb_df0)

    use mod_splitting,       only: create_rhs_bcl, rhs_momentum
    use mod_input,           only: nlayers, method_visc, dt
    use mod_grid,            only: npoin, npoin_q, nface
    use mod_constants,       only: gravity
    use mod_initial,         only: alpha_mlswe, zbot_df, pbprime_df
    use mod_basis,           only: nq, ngl
    use mod_rk_mlswe,        only: ti_barotropic_ssprk_mlswe
    use mod_variables,       only: one_plus_eta_df, dpprime_visc, dpprime_visc_q, qprime_df, qb_df
    use mod_barotropic_terms, only: btp_bcl_coeffs_qdf
    use mod_layer_terms,     only: extract_qprime_df_face, layer_mom_boundary_df, extract_velocity
    use mod_initial_mlswe,   only: poslimiter
    use mod_create_rhs_mlswe, only: layer_mass_rhs

    implicit none

    real, dimension(4,npoin),          intent(inout) :: qb_df0
    real, dimension(3,npoin,nlayers),  intent(inout) :: q_df

    real, dimension(4,npoin)           :: qbp_df
    real, dimension(3,npoin,nlayers)   :: qprime_df0, qprime_df2, q_df2, rhs
    real, dimension(2,npoin,nlayers)   :: uv_df, rhs_mom
    real, dimension(npoin,nlayers)     :: rhs_dp
    real, dimension(npoin)             :: inv_ope   ! per-node reciprocal
    integer :: k

    qb_df0 = qb_df

    ! ==================== Prediction step =================================

    qbp_df = qb_df

    call extract_qprime_df_face(qprime_df, q_df, qbp_df)

    qprime_df0 = qprime_df

    call btp_bcl_coeffs_qdf(qprime_df)
    !!$acc update device(qprime_df)
    call ti_barotropic_ssprk_mlswe(qb_df, qprime_df)

    call create_rhs_bcl(rhs, qprime_df, q_df)

    ! Forward Euler predictor step: q_df2 = q_df + dt*rhs
    q_df2 = q_df + dt*rhs

    call layer_mom_boundary_df(q_df2)

    ! Extract velocity for the predictor step and applied consistency bewteen btp and bcl
    call extract_velocity(uv_df, q_df2, qb_df)

    ! Recompute the momentum
    q_df2(2,:,:) = uv_df(1,:,:) * q_df2(1,:,:)
    q_df2(3,:,:) = uv_df(2,:,:) * q_df2(1,:,:)

    ! ==================== Correction step =================================

    qb_df = qbp_df

    call extract_qprime_df_face(qprime_df2, q_df2, qb_df)

    qprime_df = 0.5*(qprime_df2 + qprime_df)

    call btp_bcl_coeffs_qdf(qprime_df)
    !!$acc update device(qprime_df)
    call ti_barotropic_ssprk_mlswe(qb_df, qprime_df)

    ! Layer continuity equation
    call layer_mass_rhs(rhs_dp, qprime_df)

    ! Update layer thickness: q_df(1,:,:) = q_df(1,:,:) + dt*rhs_dp
    q_df(1,:,:) = q_df(1,:,:) + dt*rhs_dp

    inv_ope = pbprime_df / sum(q_df(1,:,:), dim=2)   ! store reciprocal

    do k = 1, nlayers
        qprime_df(1,:,k) = q_df(1,:,k) * inv_ope  
    end do

    ! Average the dpprime for the momentum equation
    qprime_df(1,:,:) = 0.5*(qprime_df(1,:,:) + qprime_df0(1,:,:))

    ! Layer momentum equation
    call rhs_momentum(rhs_mom, qprime_df, q_df)

    ! Update momentum
    q_df(2,:,:) = q_df(2,:,:) + dt*rhs_mom(1,:,:)
    q_df(3,:,:) = q_df(3,:,:) + dt*rhs_mom(2,:,:)

    call layer_mom_boundary_df(q_df)

    ! Extract velocity for the predictor step and applied consistency bewteen btp and bcl
    call extract_velocity(uv_df, q_df, qb_df)

    ! Recompute the momentum
    q_df(2,:,:) = uv_df(1,:,:) * q_df(1,:,:)
    q_df(3,:,:) = uv_df(2,:,:) * q_df(1,:,:)

    qb_df0 = qb_df

end subroutine ti_rk_bcl
