module mod_rk_mlswe

    ! ===========================================================================================================================
    ! This module contains the routines for the barotropic-baroclinic splitting
    !   Author: Yao Gahounzo 
    !   Date: March 27, 2023
    ! It contains the following routines:
    ! - ti_barotropic: barotropic substem for splitting system using ssprk(5,3) time integration (Ruuth, Steven. 2006)
    !
    ! ===========================================================================================================================

    use mod_initial, only: N_btp, pbprime_df, tau_wind, ssprk_a, ssprk_beta
    use mod_grid, only: npoin, npoin_q, nface
    use mod_basis, only: nqx, nqy, nqz, nq, ngl
    use mod_input, only: nlayers, dt_btp, ifilter, nlayers, kstages, method_visc
    use mod_create_rhs_mlswe, only: create_rhs_btp_momentum_new, create_rhs_btp_momentum_new1, create_rhs_btp_momentum_new2
    use mod_barotropic_terms, only: btp_bcl_coeffs, btp_evaluate_mom_dp, btp_evaluate_mom_dp_face, btp_evaluate_mom_dp_graduvdp_face

    use mod_barotropic_terms, only: btp_mom_boundary_df, compute_btp_terms, btp_laplacian_terms, btp_laplacian_terms_v1
    use mod_layer_terms, only: filter_mlswe
    use mod_laplacian_quad, only: create_laplacian_mlswe_v5, create_laplacian_mlswe_v6

    implicit none

    public :: ti_barotropic_rk_mlswe2, ti_barotropic_ssprk_mlswe, ti_barotropic_rk_mlswe3, ti_rk35_btp
    
    contains

    subroutine ti_barotropic_rk_mlswe2(one_plus_eta_out,one_plus_eta_edge_2_ave,uvb_ave, ope_ave, ope2_ave,&
        H_ave, Qu_ave, Qv_ave, Quv_ave, btp_mass_flux_ave, ope_ave_df, uvb_face_ave, &
        ope_face_ave, btp_mass_flux_face_ave, H_face_ave, &
        Qu_face_ave, Qv_face_ave, Quv_face_ave, tau_wind_ave, tau_bot_ave, qb,qb_face,&
        qb_df,qprime,qprime_face, qprime_df, flag_pred, uvb_df_ave)

        ! ===========================================================================================================================
        ! This subroutine predicts and corrects the barotropic quantities for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the barotropic variables (pb_df(or pbpert), pbub_df, pbvb_df) are stored in qb_df
        ! The quadrature points of the barotropic variables (pb, pbpert, ub, vb, pbub, pbvb) are stored in qb. pb = pbprime + pb_df
        ! The face values of the barotropic variables are stored in qb_face
        ! The quadrature points variables of(dpprime, uprime, vprime) are stored in qprime
        !
        ! During the prediction step, use the baroclinic quantities from baroclinic time level  n. 
        ! During the correction step, use averages of the predicted values and values from baroclinic time level n.
        ! 
        ! Return the next baroclinic time step values of the barotropic variables and the barotropic substep averages.
        ! ===========================================================================================================================


        implicit none

        real, dimension(4,npoin_q), intent(inout) :: qb
        real, dimension(4,2,nq,nface), intent(inout) :: qb_face
        real, dimension(4,npoin), intent(inout) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface, nlayers), intent(in) :: qprime_face
        real, dimension(nq,nface), intent(out) :: one_plus_eta_edge_2_ave
        real, dimension(npoin_q), intent(out) :: ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave
        real, dimension(2,npoin_q), intent(out) :: btp_mass_flux_ave, uvb_ave
        real, dimension(npoin), intent(out) :: ope_ave_df
        real, dimension(2,2,nq,nface), intent(out) :: uvb_face_ave
        real, dimension(2,nq,nface), intent(out) :: btp_mass_flux_face_ave, ope_face_ave
        real, dimension(nq,nface), intent(out) :: H_face_ave
        real, dimension(2,nq,nface), intent(out) :: Qu_face_ave, Qv_face_ave, Quv_face_ave
        real, dimension(npoin), intent(out) :: one_plus_eta_out
        real, dimension(2,npoin_q), intent(out) :: tau_wind_ave, tau_bot_ave
        integer, intent(in) :: flag_pred
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(2,npoin) :: uvb_df_ave

        real, dimension(npoin_q) :: one_plus_eta
        real, dimension(4,npoin) :: qb_df_pred
        real, dimension(4,npoin_q) :: qb_init
        real, dimension(4,2,nq,nface) :: qb_face_init
        real, dimension(2,npoin) :: rhs_visc_btp
        real, dimension(2,nq,nface) :: Qu_face, Qv_face, one_plus_eta_face, flux_edge, H_bcl_edge, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge
        real, dimension(nq,nface) :: H_face, one_plus_eta_edge_2, one_plus_eta_edge
        real, dimension(npoin_q) :: Quu, Qvv, Quv, H, Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, pbprime_visc
        real, dimension(npoin) :: pb_advec, tempu, tempv, one_plus_eta_df
        real, dimension(2,npoin_q) :: btp_mass_flux, tau_bot

        real, dimension(2,npoin_q,nlayers) :: grad_uprime, grad_vprime
        real, dimension(2,npoin_q) :: btp_dpuvprime, coriolis
        real, dimension(2,nq,nface) :: qb_com1,btp_dp_face
        real, dimension(2,2,nq,nface) :: qb_com2,btp_dpuv_face
        real, dimension(3,npoin) :: rhs_mom
        real, dimension(2,npoin) :: q_vic, uvb_df
        real, dimension(4,npoin) :: qb0_df, qb1_df, qb2_df
        real, dimension(4,2,nq,nface) :: grad_uvdp_face
        real, dimension(4,npoin_q) :: grad_uvdp

        integer :: mstep, I, Iq, iquad, iface, ilr, k, ik
        real, dimension(npoin) :: fdt_btp1, fdt2_btp1, a_btp1, b_btp1
        real :: N_inv, a0,a1,a2, beta, dtt

        one_plus_eta_edge_2_ave = 0.0
        uvb_ave  = 0.0
        ope_ave = 0.0
        btp_mass_flux_ave = 0.0
        H_ave = 0.0
        Qu_ave = 0.0
        Qv_ave = 0.0
        Quv_ave = 0.0
        ope_ave_df = 0.0
        uvb_face_ave  = 0.0
        ope_face_ave = 0.0
        btp_mass_flux_face_ave = 0.0
        H_face_ave = 0.0
        Qu_face_ave = 0.0
        Qv_face_ave = 0.0
        Quv_face_ave = 0.0
        tau_wind_ave = 0.0
        tau_bot_ave = 0.0
        rhs_visc_btp = 0.0
        ope2_ave = 0.0
        uvb_df_ave = 0.0

        ! Compute baroclinic coefficients in the barotropic momentum fluxes, barotropic pressure forcing, and barotropic
        ! horizontal viscosity terms.  These are needed for the barotropic momentum equation.

        call btp_bcl_coeffs(Q_uu_dp,Q_uv_dp,H_bcl,Q_vv_dp,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge,H_bcl_edge, qprime,qprime_face)

        qb_init = qb
        qb_face_init = qb_face
        qb0_df = qb_df

        call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb0_df, qprime(1,:,:))

        !call create_communicator_quad(qb_face_init,4)

        call create_communicator_quad_all(qb_face_init,grad_uvdp_face,4)

        ! ========== The time loop begins here ==================

        do mstep = 1, N_btp

            ! Communicate the grid cell face values of the barotropic variables to the neighboring processors.

            !call create_communicator_quad(qb_face_init,4)

            ik = 1
            a0 = 1.0
            a1 = 0.0
            a2 = 0.0
            beta = 0.377268915331368
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, flux_edge, btp_mass_flux, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face_init,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb0_df)

            ! Compute RHS viscosity terms

            if(method_visc > 0) then
                call create_laplacian_mlswe_v5(rhs_visc_btp,grad_uvdp,grad_uvdp_face)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb_init,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qb_face_init)

            ! Update barotropic variables

            qb_df(3,:) = qb0_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
            qb_df(4,:) = qb0_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))
            qb_df(2,:) = qb0_df(2,:) + dtt * rhs_mom(3,:)

            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            call btp_mom_boundary_df(qb_df(3:4,:))

            ! Use basis functions and the predicted degrees of freedom to compute predicted values of
            ! pb*ubar, pb*vbar, pbpert, and pb at endpoints and quadrature points in the grid cells.
            ! Then use division to compute ubar and vbar at those points.

            call btp_evaluate_mom_dp(qb,qb_df)
            call btp_evaluate_mom_dp_face(qb_face, qb)

            call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))

            ! Communication of the faces values within the processors
            !call create_communicator_quad(qb_face,4)
            call create_communicator_quad_all(qb_face,grad_uvdp_face,4)

            uvb_ave(1,:) = uvb_ave(1,:) + qb_init(3,:)/qb_init(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb_init(4,:)/qb_init(1,:)

            uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face_init(3,:,:,:)/qb_face_init(1,:,:,:)
            uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face_init(4,:,:,:)/qb_face_init(1,:,:,:)

            ope_ave = ope_ave + one_plus_eta;
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
            ope2_ave = ope2_ave + one_plus_eta**2

            uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

            ! ============ Step 2 of the RK3 scheme ============

            qb1_df = qb_df

            ik = 2
            a0 = 0.0
            a1 = 1.0
            a2 = 0.0
            beta = 0.377268915331368
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, flux_edge, btp_mass_flux, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc > 0) then
                call create_laplacian_mlswe_v5(rhs_visc_btp,grad_uvdp,grad_uvdp_face)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qb_face)

            ! Update barotropic variables

            qb_df(3,:) = qb1_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
            qb_df(4,:) = qb1_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))
            qb_df(2,:) = qb1_df(2,:) + dtt * rhs_mom(3,:)

            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            call btp_mom_boundary_df(qb_df(3:4,:))

            ! Use basis functions and the predicted degrees of freedom to compute predicted values of
            ! pb*ubar, pb*vbar, pbpert, and pb at endpoints and quadrature points in the grid cells.
            ! Then use division to compute ubar and vbar at those points.

            call btp_evaluate_mom_dp(qb,qb_df)
            call btp_evaluate_mom_dp_face(qb_face, qb)

            call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))

            ! Communication of the faces values within the processors
            !call create_communicator_quad(qb_face,4)
            call create_communicator_quad_all(qb_face,grad_uvdp_face,4)

            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face(3,:,:,:)/qb_face(1,:,:,:)
            uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face(4,:,:,:)/qb_face(1,:,:,:)

            ope_ave = ope_ave + one_plus_eta;
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
            ope2_ave = ope2_ave + one_plus_eta**2

            uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

            ! ============ Step 3 of the RK3 scheme ============

            qb1_df = qb_df
            qb2_df = qb_df
            ik = 3
            a0 = 0.355909775063327
            a1 = 1.0-a0
            a2 = 0.0
            beta = 0.242995220537396
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, flux_edge, btp_mass_flux, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc > 0) then
                call create_laplacian_mlswe_v5(rhs_visc_btp,grad_uvdp,grad_uvdp_face)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qb_face)

            ! Update barotropic variables

            qb_df(3,:) = a0*qb0_df(3,:) + a1*qb1_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
            qb_df(4,:) = a0*qb0_df(4,:) + a1*qb1_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))
            qb_df(2,:) = a0*qb0_df(2,:) + a1*qb1_df(2,:) + dtt * rhs_mom(3,:)

            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            call btp_mom_boundary_df(qb_df(3:4,:))

            ! Use basis functions and the predicted degrees of freedom to compute predicted values of
            ! pb*ubar, pb*vbar, pbpert, and pb at endpoints and quadrature points in the grid cells.
            ! Then use division to compute ubar and vbar at those points.

            call btp_evaluate_mom_dp(qb,qb_df)
            call btp_evaluate_mom_dp_face(qb_face, qb)

            call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))

            ! Communication of the faces values within the processors
            !call create_communicator_quad(qb_face,4)
            call create_communicator_quad_all(qb_face,grad_uvdp_face,4)

            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face(3,:,:,:)/qb_face(1,:,:,:)
            uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face(4,:,:,:)/qb_face(1,:,:,:)

            ope_ave = ope_ave + one_plus_eta;
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
            ope2_ave = ope2_ave + one_plus_eta**2

            uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

            ! ============ Step 4 of the RK3 scheme ============

            qb1_df = qb_df

            ik = 4
            a0 = 0.367933791638137
            a1 = 1.0-a0
            a2 = 0.0
            beta = 0.238458932846290
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, flux_edge, btp_mass_flux, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc > 0) then
                call create_laplacian_mlswe_v5(rhs_visc_btp,grad_uvdp,grad_uvdp_face)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qb_face)

            ! Update barotropic variables

            qb_df(3,:) = a0*qb0_df(3,:) + a1*qb1_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
            qb_df(4,:) = a0*qb0_df(4,:) + a1*qb1_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))
            qb_df(2,:) = a0*qb0_df(2,:) + a1*qb1_df(2,:) + dtt * rhs_mom(3,:)

            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            call btp_mom_boundary_df(qb_df(3:4,:))

            ! Use basis functions and the predicted degrees of freedom to compute predicted values of
            ! pb*ubar, pb*vbar, pbpert, and pb at endpoints and quadrature points in the grid cells.
            ! Then use division to compute ubar and vbar at those points.

            call btp_evaluate_mom_dp(qb,qb_df)
            call btp_evaluate_mom_dp_face(qb_face, qb)

            call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))

            ! Communication of the faces values within the processors
            !call create_communicator_quad(qb_face,4)
            call create_communicator_quad_all(qb_face,grad_uvdp_face,4)

            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face(3,:,:,:)/qb_face(1,:,:,:)
            uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face(4,:,:,:)/qb_face(1,:,:,:)

            ope_ave = ope_ave + one_plus_eta;
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
            ope2_ave = ope2_ave + one_plus_eta**2

            uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

            ! ============ Step 5 of the RK3 scheme ============

            qb1_df = qb_df

            ik = 5
            a0 = 0.0
            a1 = 0.762406163401431
            a2 = 1.0-a1
            beta = 0.287632146308408
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, flux_edge, btp_mass_flux, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc > 0) then
                call create_laplacian_mlswe_v5(rhs_visc_btp,grad_uvdp,grad_uvdp_face)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qb_face)

            ! Update barotropic variables

            qb_df(3,:) = a1*qb1_df(3,:) + a2*qb2_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
            qb_df(4,:) = a1*qb1_df(4,:) + a2*qb2_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))
            qb_df(2,:) = a1*qb1_df(2,:) + a2*qb2_df(2,:) + dtt * rhs_mom(3,:)

            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            call btp_mom_boundary_df(qb_df(3:4,:))

            if(mstep == N_btp .and. ifilter > 0 .and. flag_pred == 0) then 
                call filter_mlswe(qb_df(3:4,:),2)
                call btp_mom_boundary_df(qb_df(3:4,:))
                call filter_mlswe(qb_df(1,:),1)
                qb_df(2,:) = qb_df(1,:) - pbprime_df(:)
            end if

            call btp_evaluate_mom_dp(qb,qb_df)
            call btp_evaluate_mom_dp_face(qb_face, qb)

            !qb_com2(1,:,:,:) = qb_face(3,:,:,:)/qb_face(1,:,:,:)
            !qb_com2(2,:,:,:) = qb_face(4,:,:,:)/qb_face(1,:,:,:)
            !call create_communicator_quad(qb_com2,2)

            if(mstep == N_btp) then 
                qb_com2(1,:,:,:) = qb_face(3,:,:,:)/qb_face(1,:,:,:)
                qb_com2(2,:,:,:) = qb_face(4,:,:,:)/qb_face(1,:,:,:)
                call create_communicator_quad(qb_com2,2)
                uvb_face_ave = uvb_face_ave + qb_com2(1:2,:,:,:)
            else 
                !call create_communicator_quad(qb_face,4)
                call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))
                call create_communicator_quad_all(qb_face,grad_uvdp_face,4)
                uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face(3,:,:,:)/qb_face(1,:,:,:)
                uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face(4,:,:,:)/qb_face(1,:,:,:)
            end if

            ! Accumulate sums for time averaging
            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            !uvb_face_ave = uvb_face_ave + qb_com2(1:2,:,:,:)

            ope_ave = ope_ave + one_plus_eta; tau_bot_ave = tau_bot_ave + tau_bot
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
            ope2_ave = ope2_ave + one_plus_eta**2

            tau_wind_ave = tau_wind_ave + tau_wind

            uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

            qb_init = qb
            qb_face_init = qb_face
            qb0_df = qb_df

        end do

        ! Update the barotropic variables at the baroclinic time level n+1

        one_plus_eta_out = one_plus_eta_df

        N_inv = 1.0 / real(5*N_btp)

        uvb_ave = N_inv*uvb_ave

        ope_ave = N_inv*ope_ave
        H_ave = N_inv*H_ave 
        Qu_ave = N_inv*Qu_ave
        Qv_ave = N_inv*Qv_ave 
        Quv_ave = N_inv*Quv_ave 
        btp_mass_flux_ave = N_inv*btp_mass_flux_ave

        uvb_face_ave = N_inv*uvb_face_ave 

        ope_face_ave = N_inv*ope_face_ave 
        H_face_ave = N_inv*H_face_ave 
        Qu_face_ave = N_inv*Qu_face_ave 
        Qv_face_ave = N_inv*Qv_face_ave 
        Quv_face_ave = N_inv*Quv_face_ave 
        btp_mass_flux_face_ave = N_inv*btp_mass_flux_face_ave 

        one_plus_eta_edge_2_ave = N_inv*one_plus_eta_edge_2_ave 
        ope2_ave = N_inv*ope2_ave

        ope_ave_df = N_inv*ope_ave_df
        uvb_df_ave = N_inv*uvb_df_ave

        tau_wind_ave = tau_wind_ave / real(N_btp)
        tau_bot_ave = tau_bot_ave / real(N_btp)

    end subroutine ti_barotropic_rk_mlswe2

    subroutine ti_barotropic_rk_mlswe3(one_plus_eta_out,one_plus_eta_edge_2_ave,uvb_ave, ope_ave, ope2_ave,&
        H_ave, Qu_ave, Qv_ave, Quv_ave, btp_mass_flux_ave, ope_ave_df, uvb_face_ave, &
        ope_face_ave, btp_mass_flux_face_ave, H_face_ave, &
        Qu_face_ave, Qv_face_ave, Quv_face_ave, tau_wind_ave, tau_bot_ave, qb,qb_face,&
        qb_df,qprime,qprime_face, qprime_df, flag_pred, uvb_df_ave)

        ! ===========================================================================================================================
        ! This subroutine predicts and corrects the barotropic quantities for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the barotropic variables (pb_df(or pbpert), pbub_df, pbvb_df) are stored in qb_df
        ! The quadrature points of the barotropic variables (pb, pbpert, ub, vb, pbub, pbvb) are stored in qb. pb = pbprime + pb_df
        ! The face values of the barotropic variables are stored in qb_face
        ! The quadrature points variables of(dpprime, uprime, vprime) are stored in qprime
        !
        ! During the prediction step, use the baroclinic quantities from baroclinic time level  n. 
        ! During the correction step, use averages of the predicted values and values from baroclinic time level n.
        ! 
        ! Return the next baroclinic time step values of the barotropic variables and the barotropic substep averages.
        ! ===========================================================================================================================

        use mod_barotropic_terms, only: btp_extract_df_face, btp_interpolate_face
        use mod_create_rhs_mlswe, only:  create_rhs_btp_dynamics_volume_new1, Apply_btp_fluxes_new
        use mod_laplacian_quad, only: btp_create_laplacian


        implicit none

        real, dimension(4,npoin_q), intent(inout) :: qb
        real, dimension(4,2,nq,nface), intent(inout) :: qb_face
        real, dimension(4,npoin), intent(inout) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface, nlayers), intent(in) :: qprime_face
        real, dimension(nq,nface), intent(out) :: one_plus_eta_edge_2_ave
        real, dimension(npoin_q), intent(out) :: ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave
        real, dimension(2,npoin_q), intent(out) :: btp_mass_flux_ave, uvb_ave
        real, dimension(npoin), intent(out) :: ope_ave_df
        real, dimension(2,2,nq,nface), intent(out) :: uvb_face_ave
        real, dimension(2,nq,nface), intent(out) :: btp_mass_flux_face_ave, ope_face_ave
        real, dimension(nq,nface), intent(out) :: H_face_ave
        real, dimension(2,nq,nface), intent(out) :: Qu_face_ave, Qv_face_ave, Quv_face_ave
        real, dimension(npoin), intent(out) :: one_plus_eta_out
        real, dimension(2,npoin_q), intent(out) :: tau_wind_ave, tau_bot_ave
        integer, intent(in) :: flag_pred
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(2,npoin) :: uvb_df_ave

        real, dimension(npoin_q) :: one_plus_eta
        real, dimension(4,npoin) :: qb_df_pred
        real, dimension(4,npoin_q) :: qb_init
        real, dimension(4,2,nq,nface) :: qb_face_init, qb_face2
        real, dimension(4,2,ngl,nface) :: qb_df_face, qb_df_face1
        real, dimension(2,npoin) :: rhs_visc_btp
        real, dimension(2,nq,nface) :: Qu_face, Qv_face, one_plus_eta_face, flux_edge, H_bcl_edge, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge
        real, dimension(nq,nface) :: H_face, one_plus_eta_edge_2, one_plus_eta_edge
        real, dimension(npoin_q) :: Quu, Qvv, Quv, H, Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, pbprime_visc
        real, dimension(npoin) :: pb_advec, tempu, tempv, one_plus_eta_df
        real, dimension(2,npoin_q) :: btp_mass_flux, tau_bot

        real, dimension(2,npoin_q,nlayers) :: grad_uprime, grad_vprime
        real, dimension(2,npoin_q) :: btp_dpuvprime, coriolis
        real, dimension(2,nq,nface) :: qb_com1,btp_dp_face
        real, dimension(2,2,nq,nface) :: qb_com2,btp_dpuv_face
        real, dimension(3,npoin) :: rhs_mom
        real, dimension(2,npoin) :: q_vic, uvb_df
        real, dimension(4,npoin) :: qb0_df, qb1_df, qb2_df
        real, dimension(4,2,nq,nface) :: grad_uvdp_face
        real, dimension(4,npoin_q) :: grad_uvdp

        integer :: mstep, I, Iq, iquad, iface, ilr, k, ik
        real, dimension(npoin) :: fdt_btp1, fdt2_btp1, a_btp1, b_btp1
        real :: N_inv, a0,a1,a2, beta, dtt

        one_plus_eta_edge_2_ave = 0.0
        uvb_ave  = 0.0
        ope_ave = 0.0
        btp_mass_flux_ave = 0.0
        H_ave = 0.0
        Qu_ave = 0.0
        Qv_ave = 0.0
        Quv_ave = 0.0
        ope_ave_df = 0.0
        uvb_face_ave  = 0.0
        ope_face_ave = 0.0
        btp_mass_flux_face_ave = 0.0
        H_face_ave = 0.0
        Qu_face_ave = 0.0
        Qv_face_ave = 0.0
        Quv_face_ave = 0.0
        tau_wind_ave = 0.0
        tau_bot_ave = 0.0
        rhs_visc_btp = 0.0
        ope2_ave = 0.0
        uvb_df_ave = 0.0

        ! Compute baroclinic coefficients in the barotropic momentum fluxes, barotropic pressure forcing, and barotropic
        ! horizontal viscosity terms.  These are needed for the barotropic momentum equation.

        call btp_bcl_coeffs(Q_uu_dp,Q_uv_dp,H_bcl,Q_vv_dp,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge,H_bcl_edge, qprime,qprime_face)

        qb_init = qb
        qb_face_init = qb_face
        qb0_df = qb_df

        ! Extract the face values for the nodal points and communicate them to the neighboring processors.
        !call btp_extract_df_face(qb_df_face, qb_df)
        !call create_communicator_df(qb_df_face,4)

        ! Interpolate the nodal face values to the quadrature points 
        !call btp_interpolate_face(qb_face, qb_df_face)

        ! ========== The time loop begins here ==================

        do mstep = 1, N_btp

            ik = 1
            a0 = 1.0
            a1 = 0.0
            a2 = 0.0
            beta = 0.377268915331368
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, flux_edge, btp_mass_flux, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face_init,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb0_df)

            ! Compute RHS viscosity terms

            if(method_visc > 0) then
                ! Compute the gradient of uprime and vrime and multiple by dpprime. 
                ! These values will be used in the LDG method for the Laplacian term
                
                !call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))
                !call create_laplacian_mlswe_v6(rhs_visc_btp,grad_uvdp,grad_uvdp_face)

                call btp_create_laplacian(rhs_visc_btp,qprime,qprime_face,qb,qb_face, qprime_df(1,:,:))
            end if

            ! Compute RHS for the barotropic

            !call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb_init,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qb_face_init)
            
            call create_rhs_btp_dynamics_volume_new1(rhs_mom, H, qb_init, Quu, Quv, Qvv,tau_bot, btp_mass_flux)
            call Apply_btp_fluxes_new(rhs_mom,H_face,Qu_face,Qv_face, qb_face_init, flux_edge)

            ! Update barotropic variables

            qb_df(3,:) = qb0_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
            qb_df(4,:) = qb0_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))
            qb_df(2,:) = qb0_df(2,:) + dtt * rhs_mom(3,:)

            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            call btp_mom_boundary_df(qb_df(3:4,:))

            ! Use basis functions and the predicted degrees of freedom to compute predicted values of
            ! pb*ubar, pb*vbar, pbpert, and pb at endpoints and quadrature points in the grid cells.
            ! Then use division to compute ubar and vbar at those points.

            call btp_evaluate_mom_dp(qb,qb_df)

            ! Extract the face values for the nodal points and communicate them to the neighboring processors.
            call btp_extract_df_face(qb_df_face, qb_df)
            call create_communicator_df(qb_df_face,4)

            ! Interpolate the nodal face values to the quadrature points 
            call btp_interpolate_face(qb_face, qb_df_face)

            uvb_ave(1,:) = uvb_ave(1,:) + qb_init(3,:)/qb_init(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb_init(4,:)/qb_init(1,:)

            uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face_init(3,:,:,:)/qb_face_init(1,:,:,:)
            uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face_init(4,:,:,:)/qb_face_init(1,:,:,:)

            ope_ave = ope_ave + one_plus_eta;
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
            ope2_ave = ope2_ave + one_plus_eta**2

            uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

            ! ============ Step 2 of the RK3 scheme ============

            qb1_df = qb_df

            ik = 2
            a0 = 0.0
            a1 = 1.0
            a2 = 0.0
            beta = 0.377268915331368
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, flux_edge, btp_mass_flux, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc > 0) then
                ! Compute the gradient of uprime and vrime and multiple by dpprime. 
                ! These values will be used in the LDG method for the Laplacian term
                
                !call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))
                !call create_laplacian_mlswe_v6(rhs_visc_btp,grad_uvdp,grad_uvdp_face)

                call btp_create_laplacian(rhs_visc_btp,qprime,qprime_face,qb,qb_face, qprime_df(1,:,:))
            end if

            ! Compute RHS for the barotropic

            !call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qb_face)

            call create_rhs_btp_dynamics_volume_new1(rhs_mom, H, qb, Quu, Quv, Qvv,tau_bot, btp_mass_flux)
            call Apply_btp_fluxes_new(rhs_mom,H_face,Qu_face,Qv_face, qb_face, flux_edge)

            ! Update barotropic variables

            qb_df(3,:) = qb1_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
            qb_df(4,:) = qb1_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))
            qb_df(2,:) = qb1_df(2,:) + dtt * rhs_mom(3,:)

            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            call btp_mom_boundary_df(qb_df(3:4,:))

            ! Use basis functions and the predicted degrees of freedom to compute predicted values of
            ! pb*ubar, pb*vbar, pbpert, and pb at endpoints and quadrature points in the grid cells.
            ! Then use division to compute ubar and vbar at those points.

            call btp_evaluate_mom_dp(qb,qb_df)

            ! Extract the face values for the nodal points and communicate them to the neighboring processors.
            call btp_extract_df_face(qb_df_face, qb_df)
            call create_communicator_df(qb_df_face,4)

            ! Interpolate the nodal face values to the quadrature points 
            call btp_interpolate_face(qb_face, qb_df_face)

            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face(3,:,:,:)/qb_face(1,:,:,:)
            uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face(4,:,:,:)/qb_face(1,:,:,:)

            ope_ave = ope_ave + one_plus_eta;
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
            ope2_ave = ope2_ave + one_plus_eta**2

            uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

            ! ============ Step 3 of the RK3 scheme ============

            qb1_df = qb_df
            qb2_df = qb_df
            ik = 3
            a0 = 0.355909775063327
            a1 = 1.0-a0
            a2 = 0.0
            beta = 0.242995220537396
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, flux_edge, btp_mass_flux, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc > 0) then
                ! Compute the gradient of uprime and vrime and multiple by dpprime. 
                ! These values will be used in the LDG method for the Laplacian term
                
                !call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))
                !call create_laplacian_mlswe_v6(rhs_visc_btp,grad_uvdp,grad_uvdp_face)

                call btp_create_laplacian(rhs_visc_btp,qprime,qprime_face,qb,qb_face, qprime_df(1,:,:))
            end if

            ! Compute RHS for the barotropic

            !call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qb_face)
            call create_rhs_btp_dynamics_volume_new1(rhs_mom, H, qb, Quu, Quv, Qvv,tau_bot, btp_mass_flux)
            call Apply_btp_fluxes_new(rhs_mom,H_face,Qu_face,Qv_face, qb_face, flux_edge)

            ! Update barotropic variables

            qb_df(3,:) = a0*qb0_df(3,:) + a1*qb1_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
            qb_df(4,:) = a0*qb0_df(4,:) + a1*qb1_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))
            qb_df(2,:) = a0*qb0_df(2,:) + a1*qb1_df(2,:) + dtt * rhs_mom(3,:)

            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            call btp_mom_boundary_df(qb_df(3:4,:))

            ! Use basis functions and the predicted degrees of freedom to compute predicted values of
            ! pb*ubar, pb*vbar, pbpert, and pb at endpoints and quadrature points in the grid cells.
            ! Then use division to compute ubar and vbar at those points.

            call btp_evaluate_mom_dp(qb,qb_df)

            ! Extract the face values for the nodal points and communicate them to the neighboring processors.
            call btp_extract_df_face(qb_df_face, qb_df)
            call create_communicator_df(qb_df_face,4)

            ! Interpolate the nodal face values to the quadrature points 
            call btp_interpolate_face(qb_face, qb_df_face)

            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face(3,:,:,:)/qb_face(1,:,:,:)
            uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face(4,:,:,:)/qb_face(1,:,:,:)

            ope_ave = ope_ave + one_plus_eta;
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
            ope2_ave = ope2_ave + one_plus_eta**2

            uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

            ! ============ Step 4 of the RK3 scheme ============

            qb1_df = qb_df

            ik = 4
            a0 = 0.367933791638137
            a1 = 1.0-a0
            a2 = 0.0
            beta = 0.238458932846290
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, flux_edge, btp_mass_flux, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc > 0) then
                ! Compute the gradient of uprime and vrime and multiple by dpprime. 
                ! These values will be used in the LDG method for the Laplacian term
                
                !call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))
                !call create_laplacian_mlswe_v6(rhs_visc_btp,grad_uvdp,grad_uvdp_face)

                call btp_create_laplacian(rhs_visc_btp,qprime,qprime_face,qb,qb_face, qprime_df(1,:,:))
            end if

            ! Compute RHS for the barotropic

            !call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qb_face)
            call create_rhs_btp_dynamics_volume_new1(rhs_mom, H, qb, Quu, Quv, Qvv,tau_bot, btp_mass_flux)
            call Apply_btp_fluxes_new(rhs_mom,H_face,Qu_face,Qv_face, qb_face, flux_edge)

            ! Update barotropic variables

            qb_df(3,:) = a0*qb0_df(3,:) + a1*qb1_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
            qb_df(4,:) = a0*qb0_df(4,:) + a1*qb1_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))
            qb_df(2,:) = a0*qb0_df(2,:) + a1*qb1_df(2,:) + dtt * rhs_mom(3,:)

            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            call btp_mom_boundary_df(qb_df(3:4,:))

            ! Use basis functions and the predicted degrees of freedom to compute predicted values of
            ! pb*ubar, pb*vbar, pbpert, and pb at endpoints and quadrature points in the grid cells.
            ! Then use division to compute ubar and vbar at those points.

            call btp_evaluate_mom_dp(qb,qb_df)

            ! Extract the face values for the nodal points and communicate them to the neighboring processors.
            call btp_extract_df_face(qb_df_face, qb_df)
            call create_communicator_df(qb_df_face,4)

            ! Interpolate the nodal face values to the quadrature points 
            call btp_interpolate_face(qb_face, qb_df_face)

            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face(3,:,:,:)/qb_face(1,:,:,:)
            uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face(4,:,:,:)/qb_face(1,:,:,:)

            ope_ave = ope_ave + one_plus_eta;
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
            ope2_ave = ope2_ave + one_plus_eta**2

            uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

            ! ============ Step 5 of the RK3 scheme ============

            qb1_df = qb_df

            ik = 5
            a0 = 0.0
            a1 = 0.762406163401431
            a2 = 1.0-a1
            beta = 0.287632146308408
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, flux_edge, btp_mass_flux, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc > 0) then
                ! Compute the gradient of uprime and vrime and multiple by dpprime. 
                ! These values will be used in the LDG method for the Laplacian term
                
                !call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))
                !call create_laplacian_mlswe_v6(rhs_visc_btp,grad_uvdp,grad_uvdp_face)

                call btp_create_laplacian(rhs_visc_btp,qprime,qprime_face,qb,qb_face, qprime_df(1,:,:))
            end if

            ! Compute RHS for the barotropic

            !call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qb_face)
            call create_rhs_btp_dynamics_volume_new1(rhs_mom, H, qb, Quu, Quv, Qvv,tau_bot, btp_mass_flux)
            call Apply_btp_fluxes_new(rhs_mom,H_face,Qu_face,Qv_face, qb_face, flux_edge)

            ! Update barotropic variables

            qb_df(3,:) = a1*qb1_df(3,:) + a2*qb2_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
            qb_df(4,:) = a1*qb1_df(4,:) + a2*qb2_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))
            qb_df(2,:) = a1*qb1_df(2,:) + a2*qb2_df(2,:) + dtt * rhs_mom(3,:)

            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            call btp_mom_boundary_df(qb_df(3:4,:))

            if(mstep == N_btp .and. ifilter > 0 .and. flag_pred == 0) then 
                call filter_mlswe(qb_df(3:4,:),2)
                call btp_mom_boundary_df(qb_df(3:4,:))
                call filter_mlswe(qb_df(1,:),1)
                qb_df(2,:) = qb_df(1,:) - pbprime_df(:)
            end if

            call btp_evaluate_mom_dp(qb,qb_df)

            ! Extract the face values for the nodal points and communicate them to the neighboring processors.
            call btp_extract_df_face(qb_df_face, qb_df)

            call create_communicator_df(qb_df_face,4)

            !if(mstep == N_btp) then 
            !    qb_df_face1 = qb_df_face
            !end if

            ! Interpolate the nodal face values to the quadrature points 
            call btp_interpolate_face(qb_face, qb_df_face)

            uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face(3,:,:,:)/qb_face(1,:,:,:)
            uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face(4,:,:,:)/qb_face(1,:,:,:)

            ! Accumulate sums for time averaging
            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            !uvb_face_ave = uvb_face_ave + qb_com2(1:2,:,:,:)

            ope_ave = ope_ave + one_plus_eta; tau_bot_ave = tau_bot_ave + tau_bot
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
            ope2_ave = ope2_ave + one_plus_eta**2

            tau_wind_ave = tau_wind_ave + tau_wind

            uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

            qb_init = qb
            qb_face_init = qb_face
            qb0_df = qb_df

        end do

        ! Interpolate the nodal face values to the quadrature points 
        !call btp_interpolate_face(qb_face, qb_df_face1)

        ! Update the barotropic variables at the baroclinic time level n+1

        one_plus_eta_out = one_plus_eta_df

        N_inv = 1.0 / real(5*N_btp)

        uvb_ave = N_inv*uvb_ave

        ope_ave = N_inv*ope_ave
        H_ave = N_inv*H_ave 
        Qu_ave = N_inv*Qu_ave
        Qv_ave = N_inv*Qv_ave 
        Quv_ave = N_inv*Quv_ave 
        btp_mass_flux_ave = N_inv*btp_mass_flux_ave

        uvb_face_ave = N_inv*uvb_face_ave 

        ope_face_ave = N_inv*ope_face_ave 
        H_face_ave = N_inv*H_face_ave 
        Qu_face_ave = N_inv*Qu_face_ave 
        Qv_face_ave = N_inv*Qv_face_ave 
        Quv_face_ave = N_inv*Quv_face_ave 
        btp_mass_flux_face_ave = N_inv*btp_mass_flux_face_ave 

        one_plus_eta_edge_2_ave = N_inv*one_plus_eta_edge_2_ave 
        ope2_ave = N_inv*ope2_ave

        ope_ave_df = N_inv*ope_ave_df
        uvb_df_ave = N_inv*uvb_df_ave

        tau_wind_ave = tau_wind_ave / real(N_btp)
        tau_bot_ave = tau_bot_ave / real(N_btp)

    end subroutine ti_barotropic_rk_mlswe3

    subroutine ti_rk35_btp(one_plus_eta_out,one_plus_eta_edge_2_ave,uvb_ave, ope_ave, ope2_ave,&
        H_ave, Qu_ave, Qv_ave, Quv_ave, btp_mass_flux_ave, ope_ave_df, uvb_face_ave, &
        ope_face_ave, btp_mass_flux_face_ave, H_face_ave, &
        Qu_face_ave, Qv_face_ave, Quv_face_ave, tau_wind_ave, tau_bot_ave, qb,qb_face,&
        qb_df,qprime,qprime_face, qprime_df, flag_pred, uvb_df_ave)

        ! ===========================================================================================================================
        ! This subroutine predicts and corrects the barotropic quantities for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the barotropic variables (pb_df(or pbpert), pbub_df, pbvb_df) are stored in qb_df
        ! The quadrature points of the barotropic variables (pb, pbpert, ub, vb, pbub, pbvb) are stored in qb. pb = pbprime + pb_df
        ! The face values of the barotropic variables are stored in qb_face
        ! The quadrature points variables of(dpprime, uprime, vprime) are stored in qprime
        !
        ! During the prediction step, use the baroclinic quantities from baroclinic time level  n. 
        ! During the correction step, use averages of the predicted values and values from baroclinic time level n.
        ! 
        ! Return the next baroclinic time step values of the barotropic variables and the barotropic substep averages.
        ! ===========================================================================================================================

        use mod_barotropic_terms, only: btp_extract_df_face, btp_interpolate_face
        use mod_create_rhs_mlswe, only:  create_rhs_btp_dynamics_volume_new1, Apply_btp_fluxes_new
        use mod_laplacian_quad, only: btp_create_laplacian

        implicit none

        real, dimension(4,npoin_q), intent(inout) :: qb
        real, dimension(4,2,nq,nface), intent(inout) :: qb_face
        real, dimension(4,npoin), intent(inout) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface, nlayers), intent(in) :: qprime_face
        real, dimension(nq,nface), intent(out) :: one_plus_eta_edge_2_ave
        real, dimension(npoin_q), intent(out) :: ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave
        real, dimension(2,npoin_q), intent(out) :: btp_mass_flux_ave, uvb_ave
        real, dimension(npoin), intent(out) :: ope_ave_df
        real, dimension(2,2,nq,nface), intent(out) :: uvb_face_ave
        real, dimension(2,nq,nface), intent(out) :: btp_mass_flux_face_ave, ope_face_ave
        real, dimension(nq,nface), intent(out) :: H_face_ave
        real, dimension(2,nq,nface), intent(out) :: Qu_face_ave, Qv_face_ave, Quv_face_ave
        real, dimension(npoin), intent(out) :: one_plus_eta_out
        real, dimension(2,npoin_q), intent(out) :: tau_wind_ave, tau_bot_ave
        integer, intent(in) :: flag_pred
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(2,npoin) :: uvb_df_ave

        real, dimension(npoin_q) :: one_plus_eta
        real, dimension(2,npoin) :: rhs_visc_btp
        real, dimension(2,nq,nface) :: Qu_face, Qv_face, one_plus_eta_face, flux_edge, H_bcl_edge, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge
        real, dimension(nq,nface) :: H_face, one_plus_eta_edge_2, one_plus_eta_edge
        real, dimension(npoin_q) :: Quu, Qvv, Quv, H, Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl
        real, dimension(npoin) :: one_plus_eta_df
        real, dimension(2,npoin_q) :: btp_mass_flux, tau_bot

        real, dimension(4,2,nq,nface) :: qb_face2
        real, dimension(3,npoin) :: rhs_mom
        real, dimension(4,npoin) :: qb0_df, qb1_df
        real, dimension(4,2,nq,nface) :: grad_uvdp_face
        real, dimension(4,npoin_q) :: grad_uvdp
        real, dimension(4,npoin_q) :: qb_init, qbp
        real, dimension(4,2,nq,nface) :: qb_face_init, qbp_face
        real, dimension(4,2,ngl,nface) :: qb_df_face

        integer :: mstep, I, Iq, iquad, iface, ilr, k, ik
        real :: N_inv, dtt

        one_plus_eta_edge_2_ave = 0.0
        uvb_ave  = 0.0
        ope_ave = 0.0
        btp_mass_flux_ave = 0.0
        H_ave = 0.0
        Qu_ave = 0.0
        Qv_ave = 0.0
        Quv_ave = 0.0
        ope_ave_df = 0.0
        uvb_face_ave  = 0.0
        ope_face_ave = 0.0
        btp_mass_flux_face_ave = 0.0
        H_face_ave = 0.0
        Qu_face_ave = 0.0
        Qv_face_ave = 0.0
        Quv_face_ave = 0.0
        tau_wind_ave = 0.0
        tau_bot_ave = 0.0
        rhs_visc_btp = 0.0
        ope2_ave = 0.0
        uvb_df_ave = 0.0

        ! Compute baroclinic coefficients in the barotropic momentum fluxes, barotropic pressure forcing, and barotropic
        ! horizontal viscosity terms.  These are needed for the barotropic momentum equation.

        call btp_bcl_coeffs(Q_uu_dp,Q_uv_dp,H_bcl,Q_vv_dp,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge,H_bcl_edge, qprime,qprime_face)

        qb_init = qb
        qb_face_init = qb_face
        qb0_df = qb_df

        qb1_df = 0.0

        ! ========== The time loop begins here ==================

        do mstep = 1, N_btp

            do ik=1,kstages

                dtt = dt_btp*ssprk_beta(ik)

                qbp = qb
                qbp_face = qb_face

                if(ik == 1) then 
                    qbp = qb_init
                    qbp_face = qb_face_init
                end if 
        
                ! ******** Barotropic mass & momentum equation ***********

                call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                    one_plus_eta_face, flux_edge, btp_mass_flux, qbp,Q_uu_dp,Q_uv_dp,Q_vv_dp,qbp_face,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                    H_bcl, H_bcl_edge, qb_df)

                ! Compute RHS viscosity terms

                if(method_visc > 0) then
                    ! Compute the gradient of uprime and vrime and multiple by dpprime. 
                    ! These values will be used in the LDG method for the Laplacian term
                    
                    call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))
                    call create_laplacian_mlswe_v6(rhs_visc_btp,grad_uvdp,grad_uvdp_face)

                    !call btp_create_laplacian(rhs_visc_btp,qprime,qprime_face,qbp,qbp_face, qprime_df(1,:,:))
                end if

                ! Compute RHS for the barotropic

                call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qbp,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qbp_face)

                !Update Solution

                qb_df(2,:) = ssprk_a(ik,1)*qb0_df(2,:) + ssprk_a(ik,2)*qb_df(2,:) + ssprk_a(ik,3)*qb1_df(2,:) + dtt*rhs_mom(3,:)
                qb_df(3,:) = ssprk_a(ik,1)*qb0_df(3,:) + ssprk_a(ik,2)*qb_df(3,:) + ssprk_a(ik,3)*qb1_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
                qb_df(4,:) = ssprk_a(ik,1)*qb0_df(4,:) + ssprk_a(ik,2)*qb_df(4,:) + ssprk_a(ik,3)*qb1_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))

                qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

                !Apply Boundary Conditions
                call btp_mom_boundary_df(qb_df(3:4,:))
        
                !Update

                !Store the 2nd stage for SSP(5,3)
                if (ik == 2) then
                    qb1_df = qb_df
                end if

                if(mstep == N_btp .and. ifilter > 0 .and. flag_pred == 0) then 
                    call filter_mlswe(qb_df(3:4,:),2)
                    call btp_mom_boundary_df(qb_df(3:4,:))
                    call filter_mlswe(qb_df(1,:),1)
                    qb_df(2,:) = qb_df(1,:) - pbprime_df(:)
                end if

                ! Extract the face values for the nodal points and communicate them to the neighboring processors.
                call btp_extract_df_face(qb_df_face, qb_df)
                call create_communicator_df(qb_df_face,4)

                ! Interpolate the nodal face values to the quadrature points 
                call btp_evaluate_mom_dp(qb,qb_df)
                call btp_interpolate_face(qb_face, qb_df_face)

                if(ik == 1) then 

                    uvb_ave(1,:) = uvb_ave(1,:) + qbp(3,:)/qbp(1,:)
                    uvb_ave(2,:) = uvb_ave(2,:) + qbp(4,:)/qbp(1,:)

                    uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qbp_face(3,:,:,:)/qbp_face(1,:,:,:)
                    uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qbp_face(4,:,:,:)/qbp_face(1,:,:,:)
                else 
                    uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
                    uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

                    uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face(3,:,:,:)/qb_face(1,:,:,:)
                    uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face(4,:,:,:)/qb_face(1,:,:,:)
                end if

                ope_ave = ope_ave + one_plus_eta;
                btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
                H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
                Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

                ope_face_ave = ope_face_ave + one_plus_eta_face
                btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
                H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

                one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
                ope2_ave = ope2_ave + one_plus_eta**2

                uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
                uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

                tau_bot_ave = tau_bot_ave + tau_bot

            end do !ik

            qb0_df = qb_df
            qb_init = qb
            qb_face_init = qb_face

            tau_wind_ave = tau_wind_ave + tau_wind

        end do

        ! Update the barotropic variables at the baroclinic time level n+1

        one_plus_eta_out = one_plus_eta_df

        N_inv = 1.0 / real(kstages*N_btp)

        uvb_ave = N_inv*uvb_ave
        ope_ave = N_inv*ope_ave
        H_ave = N_inv*H_ave 
        Qu_ave = N_inv*Qu_ave
        Qv_ave = N_inv*Qv_ave 
        Quv_ave = N_inv*Quv_ave 
        btp_mass_flux_ave = N_inv*btp_mass_flux_ave
        uvb_face_ave = N_inv*uvb_face_ave 
        ope_face_ave = N_inv*ope_face_ave 
        H_face_ave = N_inv*H_face_ave 
        Qu_face_ave = N_inv*Qu_face_ave 
        Qv_face_ave = N_inv*Qv_face_ave 
        Quv_face_ave = N_inv*Quv_face_ave 
        btp_mass_flux_face_ave = N_inv*btp_mass_flux_face_ave 
        one_plus_eta_edge_2_ave = N_inv*one_plus_eta_edge_2_ave 
        ope2_ave = N_inv*ope2_ave
        ope_ave_df = N_inv*ope_ave_df
        uvb_df_ave = N_inv*uvb_df_ave
        tau_wind_ave = tau_wind_ave / real(N_btp)
        tau_bot_ave = tau_bot_ave*N_inv

    end subroutine ti_rk35_btp

    subroutine ti_barotropic_ssprk_mlswe(one_plus_eta_out,one_plus_eta_edge_2_ave,uvb_ave, ope_ave, ope2_ave,&
        H_ave, Qu_ave, Qv_ave, Quv_ave, btp_mass_flux_ave, ope_ave_df, uvb_face_ave, &
        ope_face_ave, btp_mass_flux_face_ave, H_face_ave, &
        Qu_face_ave, Qv_face_ave, Quv_face_ave, tau_wind_ave, tau_bot_ave, qb,qb_face,&
        qb_df,qprime,qprime_face, qprime_df, flag_pred, uvb_df_ave)

        ! ===========================================================================================================================
        ! This subroutine predicts and corrects the barotropic quantities for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the barotropic variables (pb_df(or pbpert), pbub_df, pbvb_df) are stored in qb_df
        ! The quadrature points of the barotropic variables (pb, pbpert, ub, vb, pbub, pbvb) are stored in qb. pb = pbprime + pb_df
        ! The face values of the barotropic variables are stored in qb_face
        ! The quadrature points variables of(dpprime, uprime, vprime) are stored in qprime
        !
        ! During the prediction step, use the baroclinic quantities from baroclinic time level  n. 
        ! During the correction step, use averages of the predicted values and values from baroclinic time level n.
        ! 
        ! Return the next baroclinic time step values of the barotropic variables and the barotropic substep averages.
        ! ===========================================================================================================================

        implicit none

        real, dimension(4,npoin_q), intent(inout) :: qb
        real, dimension(4,2,nq,nface), intent(inout) :: qb_face
        real, dimension(4,npoin), intent(inout) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface, nlayers), intent(in) :: qprime_face
        real, dimension(nq,nface), intent(out) :: one_plus_eta_edge_2_ave
        real, dimension(npoin_q), intent(out) :: ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave
        real, dimension(2,npoin_q), intent(out) :: btp_mass_flux_ave, uvb_ave
        real, dimension(npoin), intent(out) :: ope_ave_df
        real, dimension(2,2,nq,nface), intent(out) :: uvb_face_ave
        real, dimension(2,nq,nface), intent(out) :: btp_mass_flux_face_ave, ope_face_ave
        real, dimension(nq,nface), intent(out) :: H_face_ave
        real, dimension(2,nq,nface), intent(out) :: Qu_face_ave, Qv_face_ave, Quv_face_ave
        real, dimension(npoin), intent(out) :: one_plus_eta_out
        real, dimension(2,npoin_q), intent(out) :: tau_wind_ave, tau_bot_ave
        integer, intent(in) :: flag_pred
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(2,npoin) :: uvb_df_ave

        real, dimension(npoin_q) :: one_plus_eta
        real, dimension(2,npoin) :: rhs_visc_btp
        real, dimension(2,nq,nface) :: Qu_face, Qv_face, one_plus_eta_face, flux_edge, H_bcl_edge, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge
        real, dimension(nq,nface) :: H_face, one_plus_eta_edge_2, one_plus_eta_edge
        real, dimension(npoin_q) :: Quu, Qvv, Quv, H, Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl
        real, dimension(npoin) :: one_plus_eta_df
        real, dimension(2,npoin_q) :: btp_mass_flux, tau_bot

        real, dimension(4,2,nq,nface) :: qb_face2
        real, dimension(3,npoin) :: rhs_mom
        real, dimension(4,npoin) :: qb0_df, qb1_df
        real, dimension(4,2,nq,nface) :: grad_uvdp_face
        real, dimension(4,npoin_q) :: grad_uvdp

        integer :: mstep, I, Iq, iquad, iface, ilr, k, ik
        real :: N_inv, dtt

        one_plus_eta_edge_2_ave = 0.0
        uvb_ave  = 0.0
        ope_ave = 0.0
        btp_mass_flux_ave = 0.0
        H_ave = 0.0
        Qu_ave = 0.0
        Qv_ave = 0.0
        Quv_ave = 0.0
        ope_ave_df = 0.0
        uvb_face_ave  = 0.0
        ope_face_ave = 0.0
        btp_mass_flux_face_ave = 0.0
        H_face_ave = 0.0
        Qu_face_ave = 0.0
        Qv_face_ave = 0.0
        Quv_face_ave = 0.0
        tau_wind_ave = 0.0
        tau_bot_ave = 0.0
        rhs_visc_btp = 0.0
        ope2_ave = 0.0
        uvb_df_ave = 0.0

        ! Compute baroclinic coefficients in the barotropic momentum fluxes, barotropic pressure forcing, and barotropic
        ! horizontal viscosity terms.  These are needed for the barotropic momentum equation.

        call btp_bcl_coeffs(Q_uu_dp,Q_uv_dp,H_bcl,Q_vv_dp,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge,H_bcl_edge, qprime,qprime_face)

        qb0_df = qb_df

        ! Communicate interface values to the neighboring processors

        call create_communicator_quad(qb_face,4)

        qb1_df = 0.0

        ! ========== The time loop begins here ==================

        do mstep = 1, N_btp

            do ik=1,kstages

                dtt = dt_btp*ssprk_beta(ik)
        
                ! ******** Barotropic mass & momentum equation ***********

                call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                    one_plus_eta_face, flux_edge, btp_mass_flux, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                    H_bcl, H_bcl_edge, qb_df)

                ! Compute RHS viscosity terms

                if(method_visc > 0) then
                    ! Compute the gradient of uprime and vrime and multiple by dpprime. 
                    ! These values will be used in the LDG method for the Laplacian term
                    
                    call btp_laplacian_terms(grad_uvdp,grad_uvdp_face,qprime_df,qb_df, qprime(1,:,:))
                    call create_laplacian_mlswe_v6(rhs_visc_btp,grad_uvdp,grad_uvdp_face)
                end if

                ! Compute RHS for the barotropic

                call create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_btp,btp_mass_flux, flux_edge,qb_face)
        
                !Update Solution

                qb_df(2,:) = ssprk_a(ik,1)*qb0_df(2,:) + ssprk_a(ik,2)*qb_df(2,:) + ssprk_a(ik,3)*qb1_df(2,:) + dtt*rhs_mom(3,:)
                qb_df(3,:) = ssprk_a(ik,1)*qb0_df(3,:) + ssprk_a(ik,2)*qb_df(3,:) + ssprk_a(ik,3)*qb1_df(3,:) + dtt*(rhs_mom(1,:) + rhs_visc_btp(1,:))
                qb_df(4,:) = ssprk_a(ik,1)*qb0_df(4,:) + ssprk_a(ik,2)*qb_df(4,:) + ssprk_a(ik,3)*qb1_df(4,:) + dtt*(rhs_mom(2,:) + rhs_visc_btp(2,:))

                qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

                !Apply Boundary Conditions
                call btp_mom_boundary_df(qb_df(3:4,:))
        
                !Update

                !Store the 2nd stage for SSP(5,3)
                if (kstages == 5 .and. ik == 2) then
                    qb1_df = qb_df
                end if

                if(mstep == N_btp .and. ifilter > 0 .and. flag_pred == 0) then 
                    call filter_mlswe(qb_df(3:4,:),2)
                    call btp_mom_boundary_df(qb_df(3:4,:))
                    call filter_mlswe(qb_df(1,:),1)
                    qb_df(2,:) = qb_df(1,:) - pbprime_df(:)
                end if
    
                call btp_evaluate_mom_dp(qb,qb_df)
                call btp_evaluate_mom_dp_face(qb_face, qb)

                if(mstep == N_btp) then 
                    qb_face2 = qb_face
                end if

                ! Communication of the faces values within the processors
                call create_communicator_quad(qb_face,4)

                uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
                uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

                uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face(3,:,:,:)/qb_face(1,:,:,:)
                uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face(4,:,:,:)/qb_face(1,:,:,:)

                ope_ave = ope_ave + one_plus_eta;
                btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
                H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
                Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

                ope_face_ave = ope_face_ave + one_plus_eta_face
                btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
                H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

                one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
                ope2_ave = ope2_ave + one_plus_eta**2

                uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
                uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

                tau_bot_ave = tau_bot_ave + tau_bot

            end do !ik

            qb0_df = qb_df

            tau_wind_ave = tau_wind_ave + tau_wind

        end do

        qb_face = qb_face2

        ! Update the barotropic variables at the baroclinic time level n+1

        one_plus_eta_out = one_plus_eta_df

        N_inv = 1.0 / real(kstages*N_btp)

        uvb_ave = N_inv*uvb_ave

        ope_ave = N_inv*ope_ave
        H_ave = N_inv*H_ave 
        Qu_ave = N_inv*Qu_ave
        Qv_ave = N_inv*Qv_ave 
        Quv_ave = N_inv*Quv_ave 
        btp_mass_flux_ave = N_inv*btp_mass_flux_ave

        uvb_face_ave = N_inv*uvb_face_ave 

        ope_face_ave = N_inv*ope_face_ave 
        H_face_ave = N_inv*H_face_ave 
        Qu_face_ave = N_inv*Qu_face_ave 
        Qv_face_ave = N_inv*Qv_face_ave 
        Quv_face_ave = N_inv*Quv_face_ave 
        btp_mass_flux_face_ave = N_inv*btp_mass_flux_face_ave 

        one_plus_eta_edge_2_ave = N_inv*one_plus_eta_edge_2_ave 
        ope2_ave = N_inv*ope2_ave

        ope_ave_df = N_inv*ope_ave_df
        uvb_df_ave = N_inv*uvb_df_ave

        tau_wind_ave = tau_wind_ave / real(N_btp)
        tau_bot_ave = tau_bot_ave*N_inv

    end subroutine ti_barotropic_ssprk_mlswe

end module mod_rk_mlswe