module mod_rk_mlswe

    ! ===========================================================================================================================
    ! This module contains the routines for the barotropic-baroclinic splitting
    !   Author: Yao Gahounzo 
    !   Date: March 27, 2023
    ! It contains the following routines:
    ! - ti_barotropic: barotropic substem for splitting system using SSPRK35 time integration (see Higdon et al. 2005)
    !
    ! ===========================================================================================================================
        
    implicit none

    public :: ti_barotropic_rk_mlswe
    
    contains

    subroutine ti_barotropic_rk_mlswe(qb,qb_face,qb_df,qprime,qprime_face, qprime_df, flag_pred)

        ! ===========================================================================================================================
        ! This subroutine predicts and corrects the barotropic quantities for the splitting system using SSPRK35 time integration
        ! The nodal points or degree of freedom of the barotropic variables (pb_df(or pbpert), pbub_df, pbvb_df) are stored in qb_df
        ! The quadrature points of the barotropic variables (pb, pbpert, ub, vb, pbub, pbvb) are stored in qb. pb = pbprime + pb_df
        ! The face values of the barotropic variables are stored in qb_face
        ! The quadrature points variables of(dpprime, uprime, vprime) are stored in qprime
        ! 
        ! Return the next baroclinic time step values of the barotropic variables and the barotropic substep averages.
        ! ===========================================================================================================================

        use mod_initial, only: N_btp, pbprime_df, fdt_btp, fdt2_btp, a_btp, b_btp, tau_wind
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nqx, nqy, nqz, nq
        use mod_input, only: nlayers, dt_btp, ifilter, nlayers
        use mod_create_rhs_mlswe, only: create_rhs_btp_momentum_new1
        use mod_barotropic_terms, only: btp_mass_advection_terms, btp_bcl_coeffs, &
                                        btp_evaluate_mom_dp, btp_evaluate_mom_dp_face

        use mod_input, only: method_visc
        use mod_barotropic_terms, only: btp_mom_boundary_df, evaluate_quprime, compute_btp_terms, &
                                        evaluate_quprime2, evaluate_visc_terms
        use mod_layer_terms, only: filter_mlswe
        use mod_laplacian_quad, only: btp_create_laplacian_v1, btp_create_laplacian
        use mod_variables, only: one_plus_eta_edge_2_ave, ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave, &
                                ope_ave_df, uvb_face_ave, btp_mass_flux_face_ave, ope_face_ave, H_face_ave, &
                                Qu_face_ave, Qv_face_ave, Quv_face_ave, one_plus_eta_out, tau_wind_ave, tau_bot_ave, &
                                btp_mass_flux_ave, uvb_ave, one_plus_eta, &
                                Qu_face, Qv_face, one_plus_eta_face, flux_edge, H_bcl_edge, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, &
                                H_face, one_plus_eta_edge_2, one_plus_eta_edge, &
                                Quu, Qvv, Quv, H, Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, one_plus_eta_df, &
                                btp_mass_flux, tau_bot, uvb_ave_df


        implicit none

        real, dimension(4,npoin_q), intent(inout) :: qb
        real, dimension(4,2,nq,nface), intent(inout) :: qb_face
        real, dimension(4,npoin), intent(inout) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface, nlayers), intent(in) :: qprime_face
        integer, intent(in) :: flag_pred
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df

        !real, dimension(npoin_q) :: one_plus_eta
        real, dimension(4,npoin) :: qb_df_pred
        real, dimension(4,npoin_q) :: qb_init
        real, dimension(4,2,nq,nface) :: qb_face_init
        real, dimension(2,npoin) :: rhs_visc_btp

        real, dimension(2,2,nq,nface) :: qb_com2
        real, dimension(3,npoin) :: rhs_mom
        real, dimension(4,npoin) :: qb0_df, qb1_df, qb2_df

        integer :: mstep, I, Iq, iquad, iface, ilr, k, ik
        real :: N_inv, a0,a1,a2, beta, dtt

        one_plus_eta_edge_2_ave = 0.0
        uvb_ave  = 0.0
        uvb_ave_df  = 0.0
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

        ! Compute baroclinic coefficients in the barotropic momentum fluxes, barotropic pressure forcing, and barotropic
        ! horizontal viscosity terms.  These are needed for the barotropic momentum equation.

        call btp_bcl_coeffs(qprime,qprime_face)

        qb_init = qb
        qb_face_init = qb_face
        qb0_df = qb_df

        ! ========== The time loop begins here ==================

        do mstep = 1, N_btp

            ! Communicate the grid cell face values of the barotropic variables to the neighboring processors.

            call create_communicator_quad(qb_face_init,4)

            ik = 1
            a0 = 1.0
            a1 = 0.0
            a2 = 0.0
            beta = 0.377268915331368
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(qb_init,qb_face_init, qprime, qb0_df)

            ! Compute RHS viscosity terms

            if(method_visc == 1) then
                call btp_create_laplacian_v1(rhs_visc_btp,qprime,qprime_face,qb_init,qb_face_init)
            elseif(method_visc == 2) then
                !call btp_create_laplacian(rhs_visc_btp,qprime_df,qb0_df,qprime(1,:,:), qprime_face, qb_face)
                call btp_create_laplacian(rhs_visc_btp,qprime_df,qb0_df)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_momentum_new1(rhs_mom,qb_init,qb_face_init)

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

            ! Communication of the faces values within the processors
            call create_communicator_quad(qb_face,4)

            uvb_ave(1,:) = uvb_ave(1,:) + qb_init(3,:)/qb_init(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb_init(4,:)/qb_init(1,:)

            uvb_ave_df(1,:) = uvb_ave_df(1,:) + qb0_df(3,:)/qb0_df(1,:)
            uvb_ave_df(2,:) = uvb_ave_df(2,:) + qb0_df(4,:)/qb0_df(1,:)

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

            ! ============ Step 2 of the RK3 scheme ============

            qb1_df = qb_df

            ik = 2
            a0 = 0.0
            a1 = 1.0
            a2 = 0.0
            beta = 0.377268915331368
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(qb,qb_face, qprime, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc == 1) then
                call btp_create_laplacian_v1(rhs_visc_btp,qprime,qprime_face,qb,qb_face)
            elseif(method_visc == 2) then
                !call btp_create_laplacian(rhs_visc_btp,qprime_df,qb1_df,qprime(1,:,:), qprime_face, qb_face)
                call btp_create_laplacian(rhs_visc_btp,qprime_df,qb1_df)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_momentum_new1(rhs_mom,qb,qb_face)

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

            ! Communication of the faces values within the processors
            call create_communicator_quad(qb_face,4)

            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            uvb_ave_df(1,:) = uvb_ave_df(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_ave_df(2,:) = uvb_ave_df(2,:) + qb_df(4,:)/qb_df(1,:)

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

            call compute_btp_terms(qb,qb_face, qprime, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc == 1) then
                call btp_create_laplacian_v1(rhs_visc_btp,qprime,qprime_face,qb,qb_face)
            elseif(method_visc == 2) then
                !call btp_create_laplacian(rhs_visc_btp,qprime_df,qb1_df,qprime(1,:,:), qprime_face, qb_face)
                call btp_create_laplacian(rhs_visc_btp,qprime_df,qb1_df)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_momentum_new1(rhs_mom,qb,qb_face)

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

            ! Communication of the faces values within the processors
            call create_communicator_quad(qb_face,4)

            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            uvb_ave_df(1,:) = uvb_ave_df(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_ave_df(2,:) = uvb_ave_df(2,:) + qb_df(4,:)/qb_df(1,:)

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

            ! ============ Step 4 of the RK3 scheme ============

            qb1_df = qb_df

            ik = 4
            a0 = 0.367933791638137
            a1 = 1.0-a0
            a2 = 0.0
            beta = 0.238458932846290
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(qb,qb_face, qprime, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc == 1) then
                call btp_create_laplacian_v1(rhs_visc_btp,qprime,qprime_face,qb,qb_face)
            elseif(method_visc == 2) then
                !call btp_create_laplacian(rhs_visc_btp,qprime_df,qb1_df,qprime(1,:,:), qprime_face, qb_face)
                call btp_create_laplacian(rhs_visc_btp,qprime_df,qb1_df)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_momentum_new1(rhs_mom,qb,qb_face)

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

            ! Communication of the faces values within the processors
            call create_communicator_quad(qb_face,4)

            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            uvb_ave_df(1,:) = uvb_ave_df(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_ave_df(2,:) = uvb_ave_df(2,:) + qb_df(4,:)/qb_df(1,:)

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

            ! ============ Step 5 of the RK3 scheme ============

            qb1_df = qb_df

            ik = 5
            a0 = 0.0
            a1 = 0.762406163401431
            a2 = 1.0-a1
            beta = 0.287632146308408
            dtt = dt_btp*beta

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(qb,qb_face, qprime, qb1_df)

            ! Compute RHS viscosity terms

            if(method_visc == 1) then
                call btp_create_laplacian_v1(rhs_visc_btp,qprime,qprime_face,qb,qb_face)
            elseif(method_visc == 2) then
                !call btp_create_laplacian(rhs_visc_btp,qprime_df,qb1_df,qprime(1,:,:), qprime_face, qb_face)
                call btp_create_laplacian(rhs_visc_btp,qprime_df,qb1_df)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_momentum_new1(rhs_mom,qb,qb_face)

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

            qb_com2(1,:,:,:) = qb_face(3,:,:,:)/qb_face(1,:,:,:)
            qb_com2(2,:,:,:) = qb_face(4,:,:,:)/qb_face(1,:,:,:)
            call create_communicator_quad(qb_com2,2)

            ! Accumulate sums for time averaging
            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            uvb_ave_df(1,:) = uvb_ave_df(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_ave_df(2,:) = uvb_ave_df(2,:) + qb_df(4,:)/qb_df(1,:)

            uvb_face_ave = uvb_face_ave + qb_com2(1:2,:,:,:)

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

            qb_init = qb
            qb_face_init = qb_face
            qb0_df = qb_df

        end do

        ! Update the barotropic variables at the baroclinic time level n+1

        one_plus_eta_out = one_plus_eta_df

        N_inv = 1.0 / real(5*N_btp)

        uvb_ave = N_inv*uvb_ave
        uvb_ave_df = N_inv*uvb_ave_df

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

        tau_wind_ave = tau_wind_ave / real(N_btp)
        tau_bot_ave = tau_bot_ave / real(N_btp)

    end subroutine ti_barotropic_rk_mlswe

end module mod_rk_mlswe