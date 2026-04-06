! ===========================================================================================================================
! This module contains the routines for the predictor-corrector (for baroclinic) and the RK35 time integration (for barotropic) methods
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: October 27, 2023
! ==========================================================================================================================

subroutine ti_rk3_bcl(q_df, qb_df)

    ! q: layer variable dp, u*dp, v*dp at quad points and their face values: q_face
    ! q_df : layer variable dp, u*dp, v*dp at nodal (dof) point
    ! qprime: value dp', u' and v' at quad points and their face values: qprime_face
    ! qb : barotopic variable pb, pb_pert = pb'*eta, ub*pb, vb*pb at quad points and their face values: qb_face
    ! qb_df: : barotopic variable pb, pb_pert = pb'*eta, ub*pb, vb*pb at nodal points 
    ! qprime_df: value dp', u' and v' at nodal points
    ! qp_df_out: output variable, thickness h_k, velocity u_k,v_k, free surface ssh

    use mod_splitting, only: create_rhs_bcl
    use mod_input, only: nlayers, method_visc, kstages, dt
    use mod_grid, only: npoin, npoin_q, nface
    use mod_constants, only: gravity
    use mod_initial, only: alpha_mlswe, zbot_df, ssprk_a, ssprk_beta
    use mod_basis, only: nq, ngl
    use mod_rk_mlswe, only: ti_barotropic_ssprk_mlswe
    use mod_variables, only: one_plus_eta_df, dpprime_visc, dpprime_visc_q
    use mod_barotropic_terms, only: btp_bcl_coeffs_qdf
    use mod_layer_terms, only: extract_qprime_df_face, layer_mom_boundary_df, extract_velocity
    use mod_layer_terms, only: limit_layer_mass, limit_layer_momentum
    use mod_initial_mlswe, only: poslimiter

    implicit none

    real, dimension(4,npoin), intent(inout) :: qb_df
    real, dimension(3,npoin,nlayers), intent(inout) :: q_df
    
    real, dimension(4,npoin) :: qbp_df
    real, dimension(3,npoin,nlayers) :: qprime_df, qprime_df2, q_df2, q0_df, q1_df, rhs, q2_df
    real, dimension(2,npoin,nlayers) :: uv_df
    integer :: k, ik
    real :: dtt
    
    q0_df = q_df
    q1_df = q_df
    q2_df = 0.0
    ! ==================== Prediction step =================================

    do ik = 1, 3

        dtt = dt*ssprk_beta(ik)

        call extract_qprime_df_face(qprime_df, q1_df, qb_df)
        call btp_bcl_coeffs_qdf(qprime_df)

        call ti_barotropic_ssprk_mlswe(qb_df, qprime_df)

        call create_rhs_bcl(rhs, qprime_df, q1_df)
        
        do k = 1, nlayers
            q_df(1,:,k) = ssprk_a(ik,1)*q0_df(1,:,k) + ssprk_a(ik,2)*q1_df(1,:,k) + ssprk_a(ik,3)*q2_df(1,:,k) + dtt*rhs(1,:,k)
            q_df(2,:,k) = ssprk_a(ik,1)*q0_df(2,:,k) + ssprk_a(ik,2)*q1_df(2,:,k) + ssprk_a(ik,3)*q2_df(2,:,k) + dtt*rhs(2,:,k)
            q_df(3,:,k) = ssprk_a(ik,1)*q0_df(3,:,k) + ssprk_a(ik,2)*q1_df(3,:,k) + ssprk_a(ik,3)*q2_df(3,:,k) + dtt*rhs(3,:,k)
        end do

        call layer_mom_boundary_df(q_df)
        ! call poslimiter(q_df, alpha_mlswe)

        ! call limit_layer_mass(q_df)
        ! call limit_layer_momentum(q_df)

        ! Compute dpprime, uprime and vprime at the quad and nodal points
        call extract_velocity(uv_df, q_df, qb_df)

        do k = 1,nlayers
            q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
            q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
        end do

        call layer_mom_boundary_df(q_df)

        ! call limit_layer_mass(q_df)

        ! call poslimiter(q_df, alpha_mlswe)

        q1_df = q_df
        if(ik == 5 .and. k == 2) q2_df = q_df

    end do

end subroutine ti_rk3_bcl
