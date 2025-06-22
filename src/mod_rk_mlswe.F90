module mod_rk_mlswe

    ! ============================================================================================
    ! This module contains the routine for the barotropic solver
    !   Author: Yao Gahounzo 
    !   Date: March 27, 2023
    ! It contains the following routine:
    ! - ti_barotropic_ssprk_mlswe: based on SSPRK time integration 
    !   (see Higdon et al. 2005)
    !
    ! =============================================================================================
        
    implicit none

    public :: ti_barotropic_ssprk_mlswe
    
    contains

    subroutine ti_barotropic_ssprk_mlswe(qb_df,qprime_df)

        use mod_initial, only: N_btp, pbprime_df, tau_wind, one_over_pbprime_df, ssprk_a, ssprk_beta
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nqx, nqy, nqz, nq
        use mod_input, only: nlayers, dt_btp,nlayers, kstages, method_visc
        use mod_rhs_btp, only: create_rhs_btp
        use mod_barotropic_terms, only: btp_mom_boundary_df
        use mod_variables, only: one_plus_eta_edge_2_ave, ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, &
                                 ope2_ave, ope_ave_df, uvb_face_ave, btp_mass_flux_face_ave, &
                                 ope_face_ave, H_face_ave, Qu_face_ave, Qv_face_ave, Quv_face_ave, &
                                 one_plus_eta_out, tau_wind_ave, tau_bot_ave, &
                                 btp_mass_flux_ave, uvb_ave, uvb_ave_df

        use mod_variables, only: graduvb_face_ave, graduvb_ave

        implicit none

        real, dimension(4,npoin), intent(inout) :: qb_df
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

        ! Compute baroclinic coefficients in the barotropic momentum fluxes, barotropic 
        ! pressure forcing, and barotropic
        ! horizontal viscosity terms.  These are needed for the barotropic momentum equation.

        qb2_df = 0.0

        ! ========== The time loop begins here ==================

        do mstep = 1, N_btp

            qb0_df = qb_df
            qb1_df = qb_df

            do ik=1,kstages

                dtt = dt_btp*ssprk_beta(ik)
                ope_ave_df = ope_ave_df + (1.0 + qb1_df(2,:) * one_over_pbprime_df(:))
                uvb_ave_df(1,:) = uvb_ave_df(1,:) + qb1_df(3,:)/qb1_df(1,:)
                uvb_ave_df(2,:) = uvb_ave_df(2,:) + qb1_df(4,:)/qb1_df(1,:)

                ! Compute RHS for the barotropic
                call create_rhs_btp(rhs,qb1_df,qprime_df)

                ! Update barotropic variables

                qb_df(2,:) = ssprk_a(ik,1)*qb0_df(2,:) + ssprk_a(ik,2)*qb1_df(2,:) + &
                             ssprk_a(ik,3)*qb2_df(2,:) + dtt*rhs(1,:)
                qb_df(3,:) = ssprk_a(ik,1)*qb0_df(3,:) + ssprk_a(ik,2)*qb1_df(3,:) + &
                             ssprk_a(ik,3)*qb2_df(3,:) + dtt*rhs(2,:)
                qb_df(4,:) = ssprk_a(ik,1)*qb0_df(4,:) + ssprk_a(ik,2)*qb1_df(4,:) + &
                             ssprk_a(ik,3)*qb2_df(4,:) + dtt*rhs(3,:)

                qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

                call btp_mom_boundary_df(qb_df(3:4,:))

                !Update
                qb1_df = qb_df
        
                !Store the 2nd stage for SSP(5,3)
                if (kstages == 5 .and. ik == 2) qb2_df = qb_df

            end do 
            tau_wind_ave = tau_wind_ave + tau_wind

        end do

        ! Compute time averages of the various quantities listed earlier,
        ! over all barotropic substeps of the baroclinic time interval.   

        N_inv = 1.0 / real(kstages*N_btp)

        uvb_ave_df = N_inv*uvb_ave_df
        graduvb_face_ave = N_inv*graduvb_face_ave 
        graduvb_ave = N_inv*graduvb_ave
        tau_wind_ave = tau_wind_ave / real(N_btp)
        ope_ave_df = N_inv*ope_ave_df
        ope2_ave = N_inv*ope2_ave

        ope_ave = N_inv*ope_ave
        H_ave = N_inv*H_ave
        Qu_ave = N_inv*Qu_ave
        Qv_ave = N_inv*Qv_ave
        Quv_ave = N_inv*Quv_ave
        btp_mass_flux_ave = N_inv*btp_mass_flux_ave
        tau_bot_ave = N_inv*tau_bot_ave

        ope_face_ave = N_inv*ope_face_ave
        H_face_ave = N_inv*H_face_ave
        Qu_face_ave = N_inv*Qu_face_ave
        Qv_face_ave = N_inv*Qv_face_ave
        btp_mass_flux_face_ave = N_inv*btp_mass_flux_face_ave
        one_plus_eta_edge_2_ave = N_inv*one_plus_eta_edge_2_ave
        uvb_ave = N_inv*uvb_ave
        uvb_face_ave = N_inv*uvb_face_ave

    end subroutine ti_barotropic_ssprk_mlswe


end module mod_rk_mlswe
