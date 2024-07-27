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

    public :: ti_barotropic_rk_mlswe, ti_barotropic_ssprk_mlswe
    
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

        use mod_initial, only: N_btp, pbprime_df, fdt_btp, fdt2_btp, a_btp, b_btp, tau_wind, one_over_pbprime_df, ssprk_a, ssprk_beta
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nqx, nqy, nqz, nq
        use mod_input, only: nlayers, dt_btp, ifilter, nlayers, kstages, method_visc
        use mod_rhs_btp, only: create_rhs_btp
        use mod_barotropic_terms, only: btp_mass_advection_terms, btp_bcl_coeffs, &
                                        btp_evaluate_mom_dp, btp_evaluate_mom_dp_face, btp_bcl_grad_coeffs

        use mod_barotropic_terms, only: btp_mom_boundary_df, compute_btp_terms, btp_interpolate_avg_v1
        use mod_laplacian_quad, only: btp_create_laplacian_v1, btp_create_laplacian, btp_create_laplacian_v3
        use mod_variables, only: one_plus_eta_edge_2_ave, ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave, &
                                ope_ave_df, uvb_face_ave, btp_mass_flux_face_ave, ope_face_ave, H_face_ave, &
                                Qu_face_ave, Qv_face_ave, Quv_face_ave, one_plus_eta_out, tau_wind_ave, tau_bot_ave, &
                                btp_mass_flux_ave, uvb_ave, one_plus_eta, &
                                Qu_face, Qv_face, one_plus_eta_face, flux_edge, H_bcl_edge, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, &
                                H_face, one_plus_eta_edge_2, one_plus_eta_edge, &
                                Quu, Qvv, Quv, H, Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, one_plus_eta_df, &
                                btp_mass_flux, tau_bot, uvb_ave_df

        use mod_variables, only: graduvb_face_ave, graduvb_ave


        implicit none

        real, dimension(4,npoin_q), intent(inout) :: qb
        real, dimension(4,2,nq,nface), intent(inout) :: qb_face
        real, dimension(4,npoin), intent(inout) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface, nlayers), intent(in) :: qprime_face
        integer, intent(in) :: flag_pred
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        !real, dimension(npoin_q) :: one_plus_eta
        real, dimension(4,npoin) :: qb_df_pred
        real, dimension(4,npoin_q) :: qb_init
        real, dimension(4,2,nq,nface) :: qb_face_init
        real, dimension(2,npoin) :: rhs_visc_btp, uvdp_ave

        real, dimension(2,2,nq,nface) :: qb_com2
        real, dimension(3,npoin) :: rhs
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

        graduvb_face_ave = 0.0
        graduvb_ave = 0.0

        ! Compute baroclinic coefficients in the barotropic momentum fluxes, barotropic pressure forcing, and barotropic
        ! horizontal viscosity terms.  These are needed for the barotropic momentum equation.

        call btp_bcl_coeffs(qprime,qprime_face)
        if (method_visc > 0) call btp_bcl_grad_coeffs(qprime_df)

        qb2_df = 0.0

        ! ========== The time loop begins here ==================

        do mstep = 1, N_btp

            qb0_df = qb_df
            qb1_df = qb_df

            do ik=1,kstages

                dtt = dt_btp*ssprk_beta(ik)

                ope_ave_df = ope_ave_df + (1.0 + qb1_df(2,:) * one_over_pbprime_df(:))

                ! Compute RHS for the barotropic

                call create_rhs_btp(rhs,qb,qb1_df,qprime,qb_face)

                ! Update barotropic variables

                qb_df(2,:) = ssprk_a(ik,1)*qb0_df(2,:) + ssprk_a(ik,2)*qb1_df(2,:) + ssprk_a(ik,3)*qb2_df(2,:) + dtt*rhs(1,:)
                qb_df(3,:) = ssprk_a(ik,1)*qb0_df(3,:) + ssprk_a(ik,2)*qb1_df(3,:) + ssprk_a(ik,3)*qb2_df(3,:) + dtt*rhs(2,:)
                qb_df(4,:) = ssprk_a(ik,1)*qb0_df(4,:) + ssprk_a(ik,2)*qb1_df(4,:) + ssprk_a(ik,3)*qb2_df(4,:) + dtt*rhs(3,:)

                qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

                call btp_mom_boundary_df(qb_df(3:4,:))

                !Update
                qb1_df = qb_df
        
                !Store the 2nd stage for SSP(5,3)
                if (kstages == 5 .and. ik == 2) then
                    qb2_df = qb_df
                end if

                call btp_evaluate_mom_dp(qb,qb_df)
                call btp_evaluate_mom_dp_face(qb_face, qb)

                uvb_ave_df(1,:) = uvb_ave_df(1,:) + qb_df(3,:)/qb_df(1,:)
                uvb_ave_df(2,:) = uvb_ave_df(2,:) + qb_df(4,:)/qb_df(1,:)

            end do 

            tau_wind_ave = tau_wind_ave + tau_wind

        end do

        ! Compute time averages of the various quantities listed earlier,
        ! over all barotropic substeps of the baroclinic time interval.   

        one_plus_eta_out = one_plus_eta_df

        N_inv = 1.0 / real(kstages*N_btp)

        !uvb_ave = N_inv*uvb_ave
        uvb_ave_df = N_inv*uvb_ave_df

        ope_ave = N_inv*ope_ave
        H_ave = N_inv*H_ave 
        Qu_ave = N_inv*Qu_ave
        Qv_ave = N_inv*Qv_ave 
        Quv_ave = N_inv*Quv_ave 
        btp_mass_flux_ave = N_inv*btp_mass_flux_ave

        !uvb_face_ave = N_inv*uvb_face_ave 

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
        tau_bot_ave = N_inv*tau_bot_ave

        graduvb_face_ave = N_inv*graduvb_face_ave 
        graduvb_ave = N_inv*graduvb_ave

        call btp_interpolate_avg_v1()

    end subroutine ti_barotropic_rk_mlswe


    subroutine ti_barotropic_ssprk_mlswe(qb_df,qprime,qprime_face, qprime_df)

        ! ===========================================================================================================================
        ! This subroutine predicts and corrects the barotropic quantities for the splitting system using SSPRK35 time integration
        ! The nodal points or degree of freedom of the barotropic variables (pb_df(or pbpert), pbub_df, pbvb_df) are stored in qb_df
        ! The quadrature points of the barotropic variables (pb, pbpert, ub, vb, pbub, pbvb) are stored in qb. pb_df = pbprime + pb_df
        ! The face values of the barotropic variables are stored in qb_face
        ! The quadrature points variables of(dpprime, uprime, vprime) are stored in qprime
        ! 
        ! Return the next baroclinic time step values of the barotropic variables and the barotropic substep averages.
        ! ===========================================================================================================================

        use mod_initial, only: N_btp, pbprime_df, fdt_btp, fdt2_btp, a_btp, b_btp, tau_wind, one_over_pbprime_df, ssprk_a, ssprk_beta
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nqx, nqy, nqz, nq
        use mod_input, only: nlayers, dt_btp, ifilter, nlayers, kstages, method_visc
        use mod_rhs_btp, only: create_rhs_btp_v1, create_rhs_btp_v2
        use mod_barotropic_terms, only: btp_bcl_coeffs, btp_bcl_grad_coeffs

        use mod_barotropic_terms, only: btp_mom_boundary_df, btp_interpolate_avg_v1
        use mod_variables, only: one_plus_eta_edge_2_ave, ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave, &
                                ope_ave_df, uvb_face_ave, btp_mass_flux_face_ave, ope_face_ave, H_face_ave, &
                                Qu_face_ave, Qv_face_ave, Quv_face_ave, one_plus_eta_out, tau_wind_ave, tau_bot_ave, &
                                btp_mass_flux_ave, uvb_ave, uvb_ave_df

        use mod_variables, only: graduvb_face_ave, graduvb_ave


        implicit none

        real, dimension(4,npoin), intent(inout) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface, nlayers), intent(in) :: qprime_face
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        real, dimension(3,npoin) :: rhs
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
        ope2_ave = 0.0

        graduvb_face_ave = 0.0
        graduvb_ave = 0.0

        ! Compute baroclinic coefficients in the barotropic momentum fluxes, barotropic pressure forcing, and barotropic
        ! horizontal viscosity terms.  These are needed for the barotropic momentum equation.

        call btp_bcl_coeffs(qprime,qprime_face)
        if (method_visc > 0) call btp_bcl_grad_coeffs(qprime_df)

        qb2_df = 0.0

        ! ========== The time loop begins here ==================

        do mstep = 1, N_btp

            qb0_df = qb_df
            qb1_df = qb_df

            do ik=1,kstages

                dtt = dt_btp*ssprk_beta(ik)

                ope_ave_df = ope_ave_df + (1.0 + qb1_df(2,:) * one_over_pbprime_df(:))

                ! Compute RHS for the barotropic

                !call create_rhs_btp_v1(rhs,qb1_df,qprime)
                call create_rhs_btp_v2(rhs,qb1_df,qprime)

                ! Update barotropic variables

                qb_df(2,:) = ssprk_a(ik,1)*qb0_df(2,:) + ssprk_a(ik,2)*qb1_df(2,:) + ssprk_a(ik,3)*qb2_df(2,:) + dtt*rhs(1,:)
                qb_df(3,:) = ssprk_a(ik,1)*qb0_df(3,:) + ssprk_a(ik,2)*qb1_df(3,:) + ssprk_a(ik,3)*qb2_df(3,:) + dtt*rhs(2,:)
                qb_df(4,:) = ssprk_a(ik,1)*qb0_df(4,:) + ssprk_a(ik,2)*qb1_df(4,:) + ssprk_a(ik,3)*qb2_df(4,:) + dtt*rhs(3,:)

                qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

                call btp_mom_boundary_df(qb_df(3:4,:))

                !Update
                qb1_df = qb_df
        
                !Store the 2nd stage for SSP(5,3)
                if (kstages == 5 .and. ik == 2) then
                    qb2_df = qb_df
                end if

                uvb_ave_df(1,:) = uvb_ave_df(1,:) + qb_df(3,:)/qb_df(1,:)
                uvb_ave_df(2,:) = uvb_ave_df(2,:) + qb_df(4,:)/qb_df(1,:)

            end do 

            tau_wind_ave = tau_wind_ave + tau_wind

        end do

        ! Compute time averages of the various quantities listed earlier,
        ! over all barotropic substeps of the baroclinic time interval.   

        !one_plus_eta_out = one_plus_eta_df

        N_inv = 1.0 / real(kstages*N_btp)

        !uvb_ave = N_inv*uvb_ave
        uvb_ave_df = N_inv*uvb_ave_df

        ope_ave = N_inv*ope_ave
        H_ave = N_inv*H_ave 
        Qu_ave = N_inv*Qu_ave
        Qv_ave = N_inv*Qv_ave 
        Quv_ave = N_inv*Quv_ave 
        btp_mass_flux_ave = N_inv*btp_mass_flux_ave

        !uvb_face_ave = N_inv*uvb_face_ave 

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
        tau_bot_ave = N_inv*tau_bot_ave

        graduvb_face_ave = N_inv*graduvb_face_ave 
        graduvb_ave = N_inv*graduvb_ave

        call btp_interpolate_avg_v1()

    end subroutine ti_barotropic_ssprk_mlswe


end module mod_rk_mlswe