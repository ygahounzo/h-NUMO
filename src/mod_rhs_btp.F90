! =======================================================================
!>@brief This subroutine builds the RHS vector for the DG method for MLSWE
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: April 03, 2023
! =======================================================================
module mod_rhs_btp


    use mod_constants, only: gravity

    use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad

    use mod_basis, only: nglx, ngly, nglz, npts, dpsiqx, dpsiqy, dpsiqz, nqx, nqy, nqz, &
        psiqx, psiqy, psiqz, nq, ngl

    use mod_grid, only: intma, npoin_q, npoin, nface

    use mod_laplacian_quad, only: btp_create_laplacian
    use mod_barotropic_terms, only: btp_extract_df

    use mod_input, only: nlayers, method_visc

    use mod_metrics, only: &
        ksiq_x, ksiq_y, ksiq_z, &
        etaq_x, etaq_y, etaq_z, &
        zetaq_x, zetaq_y, zetaq_z, &
        jacq, massinv

    public :: create_rhs_btp, create_rhs_btp2, create_rhs_btp3, &
                create_rhs_btp_momentum, create_rhs_btp_mom_mass, btp_mass_advection_rhs
             

contains

    subroutine create_rhs_btp(rhs,qb_df,qprime)

        implicit none

        real, dimension(3, npoin), intent(out) :: rhs
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime

        real, dimension(4, 2, ngl, nface) :: qb_df_face
        real, dimension(2,npoin) :: rhs_visc_btp

        call btp_extract_df(qb_df_face, qb_df)

        call btp_create_precommunicator(qb_df_face,4)

        call create_rhs_btp_volume_qdf(rhs, qb_df, qprime)

        call btp_create_postcommunicator(qb_df_face,4)

        call creat_btp_fluxes_qdf(rhs,qb_df_face)

        ! Compute RHS viscosity terms

        if(method_visc > 0) call btp_create_laplacian(rhs_visc_btp,qb_df)

        rhs(2,:) = rhs(2,:) + rhs_visc_btp(1,:)
        rhs(3,:) = rhs(3,:) + rhs_visc_btp(2,:)

    end subroutine create_rhs_btp

    subroutine create_rhs_btp2(rhs,qb_df,qprime_df)

        implicit none

        real, dimension(3, npoin), intent(out) :: rhs
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        real, dimension(4, 2, ngl, nface) :: qb_df_face
        real, dimension(2,npoin) :: rhs_visc_btp

        call btp_extract_df(qb_df_face, qb_df)

        call btp_create_precommunicator(qb_df_face,4)

        call create_rhs_btp_volume_qdf2(rhs, qb_df, qprime_df)

        call btp_create_postcommunicator(qb_df_face,4)

        call creat_btp_fluxes_qdf(rhs,qb_df_face)

        ! Compute RHS viscosity terms

        if(method_visc > 0) call btp_create_laplacian(rhs_visc_btp,qb_df)

        rhs(2,:) = rhs(2,:) + rhs_visc_btp(1,:)
        rhs(3,:) = rhs(3,:) + rhs_visc_btp(2,:)

    end subroutine create_rhs_btp2

    subroutine create_rhs_btp3(rhs,qb_df,qprime_df)

        implicit none

        real, dimension(3, npoin), intent(out) :: rhs
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        real, dimension(4, 2, ngl, nface) :: qb_df_face
        real, dimension(2,npoin) :: rhs_visc_btp

        call btp_extract_df(qb_df_face, qb_df)

        call btp_create_precommunicator(qb_df_face,4)

        call create_rhs_btp_volume_qdf3(rhs, qb_df, qprime_df)

        call btp_create_postcommunicator(qb_df_face,4)

        call creat_btp_fluxes_qdf3(rhs,qb_df_face)
        !call creat_btp_fluxes_qdf5(rhs,qb_df_face)
        !call creat_btp_fluxes_qdf4(rhs,qb_df_face)

        ! Compute RHS viscosity terms

        if(method_visc > 0) call btp_create_laplacian(rhs_visc_btp,qb_df)

        rhs(2,:) = rhs(2,:) + rhs_visc_btp(1,:)
        rhs(3,:) = rhs(3,:) + rhs_visc_btp(2,:)

    end subroutine create_rhs_btp3

    subroutine create_rhs_btp_momentum(rhs_mom,qb,qb_face)

        implicit none

        real, dimension(2, npoin), intent(out) :: rhs_mom
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face
        
        integer :: I

        call create_rhs_btp_dynamics_volume(rhs_mom,qb)

        call Apply_btp_fluxes(rhs_mom,qb_face)

    end subroutine create_rhs_btp_momentum

    subroutine create_rhs_btp_mom_mass(rhs,qb,qb_face)

        implicit none

        real, dimension(3, npoin), intent(out) :: rhs
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face

        call create_rhs_btp_dynamics_volume_all(rhs, qb)

        call Apply_btp_fluxes_all(rhs,qb_face)

    end subroutine create_rhs_btp_mom_mass

    subroutine btp_mass_advection_rhs(pb_advec)

        implicit none

        real, dimension(npoin), intent(out) :: pb_advec

        call create_rhs_btp_dynamics_volume_mass(pb_advec)
        call  Apply_btp_flux_mass(pb_advec)

    end subroutine btp_mass_advection_rhs

    subroutine create_rhs_btp_volume(rhs, qb, qprime)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, coriolis_quad, one_over_pbprime
        use mod_variables, only: tau_bot_ave, H_ave, Qu_ave, Quv_ave, Qv_ave, ope_ave, uvb_ave, btp_mass_flux_ave
        use mod_variables, only: Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl
        use mod_input, only: nlayers, cd_mlswe

        implicit none 

        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,npoin), intent(out) :: rhs

        real :: source_x, source_y, Hq, quux, quvxy, qvvy
        real :: wq, hi, dhdx, dhdy, speed1, tau_bot_u, tau_bot_v, ope, &
                dp, dpp, udp, vdp, ub, vb

        integer :: I, Iq, ip

        speed1 = cd_mlswe/gravity

        rhs = 0.0

        do Iq = 1, npoin_q

            dp = qb(1,Iq)
            dpp = qb(2,Iq)
            udp = qb(3,Iq)
            vdp = qb(4,Iq) 

            wq = wjac(Iq)

            ub = udp/dp; vb = vdp/dp
        
            tau_bot_u = speed1 * qprime(1,Iq,nlayers)*(qprime(2,Iq,nlayers) + ub)
            tau_bot_v = speed1 * qprime(1,Iq,nlayers)*(qprime(3,Iq,nlayers) + vb)

            source_x = coriolis_quad(Iq)*vdp + gravity*(tau_wind(1,Iq) - tau_bot_u) - &
                        gravity*dp*grad_zbot_quad(1,Iq)

            source_y = -coriolis_quad(Iq)*udp + gravity*(tau_wind(2,Iq) - tau_bot_v) - &
                        gravity*dp*grad_zbot_quad(2,Iq)

            ope = 1.0 + dpp * one_over_pbprime(Iq)
            Hq = (ope**2) * H_bcl(Iq)

            quux = ub * udp + ope * Q_uu_dp(Iq)
            quvxy = ub * vdp + ope * Q_uv_dp(Iq)
            qvvy = vb * vdp + ope * Q_vv_dp(Iq)

            H_ave(Iq) = H_ave(Iq) + Hq;  Qu_ave(Iq) = Qu_ave(Iq) + quux
            Qv_ave(Iq) = Qv_ave(Iq) + qvvy;  Quv_ave(Iq) = Quv_ave(Iq) + quvxy

            tau_bot_ave(1,Iq) = tau_bot_ave(1,Iq) + tau_bot_u
            tau_bot_ave(2,Iq) = tau_bot_ave(2,Iq) + tau_bot_v

            ope_ave(Iq) = ope_ave(Iq) + ope

            !uvb_ave(1,Iq) = uvb_ave(1,Iq) + ub
            !uvb_ave(2,Iq) = uvb_ave(2,Iq) + vb

            btp_mass_flux_ave(1,Iq) = btp_mass_flux_ave(1,Iq) + udp
            btp_mass_flux_ave(2,Iq) = btp_mass_flux_ave(2,Iq) + vdp

            do ip = 1,npts

                I = indexq(ip,Iq)
                hi = psih(ip,Iq)

                !Xi derivatives                                                                                
                dhdx = dpsidx(ip,Iq)
                !Eta derivatives
                dhdy = dpsidy(ip,Iq)

                rhs(1,I) = rhs(1,I) + wq*(dhdx*udp + dhdy*vdp)
                rhs(2,I) = rhs(2,I) + wq*(hi*source_x + dhdx*(Hq + quux) + quvxy*dhdy)
                rhs(3,I) = rhs(3,I) + wq*(hi*source_y + dhdx*quvxy + dhdy*(Hq + qvvy))

            end do
        end do

    end subroutine create_rhs_btp_volume

    subroutine create_rhs_btp_volume_qdf(rhs, qb_df, qprime)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, &
                                coriolis_quad, one_over_pbprime, alpha_mlswe
        use mod_variables, only: tau_bot_ave, H_ave, Qu_ave, Quv_ave, Qv_ave, ope_ave, uvb_ave, btp_mass_flux_ave
        use mod_variables, only: Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl
        use mod_input, only: nlayers, cd_mlswe, botfr

        implicit none 

        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,npoin), intent(out) :: rhs

        real :: source_x, source_y, Hq, quux, quvxy, qvvy
        real :: wq, hi, dhdx, dhdy, coef_fric, tau_bot_u, tau_bot_v, ope, &
                dp, dpp, udp, vdp, ub, vb, Pstress, Pbstress, ubot, vbot, speed

        integer :: I, Iq, ip

        !speed1 = cd_mlswe/gravity

        rhs = 0.0
        tau_bot_u = 0.0 ; tau_bot_v = 0.0
        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth 
                                                  ! at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth
                                                         ! at which bottom stress is reduced to 0

        do Iq = 1, npoin_q

            dp = 0.0; dpp = 0.0; udp = 0.0; vdp = 0.0

            do ip = 1,npts 
                I = indexq(ip,Iq)
                hi = psih(ip,Iq)

                dp = dp + hi*qb_df(1,I)
                dpp = dpp + hi*qb_df(2,I)
                udp = udp + hi*qb_df(3,I)
                vdp = vdp + hi*qb_df(4,I)
            end do 

            wq = wjac(Iq)
        
            ub = udp/dp; vb = vdp/dp

            if (botfr == 1) then
                
                ubot = qprime(2,Iq,nlayers) + ub
                vbot = qprime(3,Iq,nlayers) + vb

                tau_bot_u = (cd_mlswe/alpha_mlswe(nlayers))*ubot
                tau_bot_v = (cd_mlswe/alpha_mlswe(nlayers))*vbot
                
            elseif (botfr == 2) then

                ubot = qprime(2,Iq,nlayers) + ub
                vbot = qprime(3,Iq,nlayers) + vb
                speed = (cd_mlswe/gravity)*sqrt(ubot**2 + vbot**2)

                tau_bot_u = speed*ubot
                tau_bot_v = speed*vbot
            end if

            source_x = coriolis_quad(Iq)*vdp + gravity*(tau_wind(1,Iq) - tau_bot_u) - &
                        gravity*dp*grad_zbot_quad(1,Iq)

            source_y = -coriolis_quad(Iq)*udp + gravity*(tau_wind(2,Iq) - tau_bot_v) - &
                        gravity*dp*grad_zbot_quad(2,Iq)

            ope = 1.0 + dpp * one_over_pbprime(Iq)
            Hq = (ope**2) * H_bcl(Iq)

            quux = ub * udp + ope * Q_uu_dp(Iq)
            quvxy = ub * vdp + ope * Q_uv_dp(Iq)
            qvvy = vb * vdp + ope * Q_vv_dp(Iq)

            H_ave(Iq) = H_ave(Iq) + Hq;  Qu_ave(Iq) = Qu_ave(Iq) + quux
            Qv_ave(Iq) = Qv_ave(Iq) + qvvy;  Quv_ave(Iq) = Quv_ave(Iq) + quvxy

            tau_bot_ave(1,Iq) = tau_bot_ave(1,Iq) + tau_bot_u
            tau_bot_ave(2,Iq) = tau_bot_ave(2,Iq) + tau_bot_v

            ope_ave(Iq) = ope_ave(Iq) + ope

            btp_mass_flux_ave(1,Iq) = btp_mass_flux_ave(1,Iq) + udp
            btp_mass_flux_ave(2,Iq) = btp_mass_flux_ave(2,Iq) + vdp

            do ip = 1,npts

                I = indexq(ip,Iq)

                hi = psih(ip,Iq)

                !Xi derivatives
                dhdx = dpsidx(ip,Iq)
                !Eta derivatives
                dhdy = dpsidy(ip,Iq)

                rhs(1,I) = rhs(1,I) + wq*(dhdx*udp + dhdy*vdp)
                rhs(2,I) = rhs(2,I) + wq*(hi*source_x + dhdx*(Hq + quux) + quvxy*dhdy)
                rhs(3,I) = rhs(3,I) + wq*(hi*source_y + dhdx*quvxy + dhdy*(Hq + qvvy))

            end do
        end do

    end subroutine create_rhs_btp_volume_qdf

    subroutine create_rhs_btp_volume_qdf2(rhs, qb_df, qprime_df)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_df, tau_wind_df, psih, dpsidx,dpsidy, indexq, wjac, &
                                coriolis_df, one_over_pbprime_df, alpha_mlswe
        use mod_input, only: nlayers, cd_mlswe, botfr
        use mod_variables, only: tau_bot_ave, H_ave, Qu_ave, Quv_ave, Qv_ave, ope_ave, uvb_ave, btp_mass_flux_ave
        use mod_initial, only: grad_zbot_quad, one_over_pbprime
        use mod_variables, only: Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl

        implicit none 

        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(3,npoin), intent(out) :: rhs

        real :: source_xq, source_yq, Hq, quux, quvxy, qvvy
        real :: wq, hi, dhdx, dhdy, coef_fric, tau_bot_uq, tau_bot_vq, opeq, &
                udpq, vdpq, ub, vb, Pstress, Pbstress, ubot, vbot, speed, dpq, dppq
        real, dimension(npoin) :: Htmp, ope, tau_bot_u, tau_bot_v, source_x, source_y, quu, quv, qvv

        integer :: I, Iq, ip

        !speed1 = cd_mlswe/gravity

        rhs = 0.0
        tau_bot_u = 0.0 ; tau_bot_v = 0.0
        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth
                                                  ! at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth
                                                         ! at which bottom stress is reduced to 0

        call nodal_computation(Htmp, ope, tau_bot_u, tau_bot_v, source_x, source_y, quu, quv, qvv, qb_df, qprime_df)

        do Iq = 1, npoin_q

            udpq = 0.0; vdpq = 0.0
            dpq = 0.0; dppq = 0.0
            Hq = 0.0; tau_bot_uq = 0.0; tau_bot_vq = 0.0; quux = 0.0
            quvxy = 0.0; qvvy = 0.0; opeq = 0.0; source_xq = 0.0; source_yq = 0.0

            do ip = 1,npts
                I = indexq(ip,Iq)
                hi = psih(ip,Iq)
                
                dpq = dpq + hi*qb_df(1,I)
                dppq = dppq + hi*qb_df(2,I)
                udpq = udpq + hi*qb_df(3,I)
                vdpq = vdpq + hi*qb_df(4,I)

                source_xq = source_xq + hi*source_x(I) 
                source_yq = source_yq + hi*source_y(I)
                opeq = opeq + hi*ope(I)
                Hq = Hq + hi*Htmp(I)
                tau_bot_uq = tau_bot_uq + hi*tau_bot_u(I)
                tau_bot_vq = tau_bot_vq + hi*tau_bot_v(I)

                quux = quux + hi*quu(I)
                quvxy = quvxy + hi*quv(I)
                qvvy = qvvy + hi*qvv(I)
            end do

            wq = wjac(Iq)
        
            opeq = 1.0 + dppq * one_over_pbprime(Iq)
            Hq = (opeq**2) * H_bcl(Iq)
            H_ave(Iq) = H_ave(Iq) + Hq;

            Qu_ave(Iq) = Qu_ave(Iq) + quux
            Qv_ave(Iq) = Qv_ave(Iq) + qvvy;  Quv_ave(Iq) = Quv_ave(Iq) + quvxy
            tau_bot_ave(1,Iq) = tau_bot_ave(1,Iq) + tau_bot_uq
            tau_bot_ave(2,Iq) = tau_bot_ave(2,Iq) + tau_bot_vq
            ope_ave(Iq) = ope_ave(Iq) + opeq
            btp_mass_flux_ave(1,Iq) = btp_mass_flux_ave(1,Iq) + udpq
            btp_mass_flux_ave(2,Iq) = btp_mass_flux_ave(2,Iq) + vdpq

            do ip = 1,npts

                I = indexq(ip,Iq)

                hi = psih(ip,Iq)

                !Xi derivatives
                dhdx = dpsidx(ip,Iq)
                !Eta derivatives
                dhdy = dpsidy(ip,Iq)

                rhs(1,I) = rhs(1,I) + wq*(dhdx*udpq + dhdy*vdpq)
                rhs(2,I) = rhs(2,I) + wq*(hi*source_xq + dhdx*(Hq + quux) + quvxy*dhdy)
                rhs(3,I) = rhs(3,I) + wq*(hi*source_yq + dhdx*quvxy + dhdy*(Hq + qvvy))

            end do
        end do

    end subroutine create_rhs_btp_volume_qdf2

    subroutine create_rhs_btp_volume_qdf3(rhs, qb_df, qprime_df)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_df, tau_wind_df, coriolis_df, one_over_pbprime_df, alpha_mlswe
        use mod_initial, only: psih_df, dpsidx_df,dpsidy_df, index_df, wjac_df
        use mod_variables, only: tau_bot_ave_df, H_ave_df, Qu_ave_df, Quv_ave_df, Qv_ave_df, &
                                 ope_ave_df, uvb_ave_df, btp_mass_flux_ave_df
        use mod_variables, only: Q_uu_dp_df, Q_uv_dp_df, Q_vv_dp_df, H_bcl_df
        use mod_input, only: nlayers, cd_mlswe, botfr

        implicit none

        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(3,npoin), intent(inout) :: rhs

        real :: source_x, source_y, Hq, quux, quvxy, qvvy
        real :: wq, hi, dhdx, dhdy, coef_fric, tau_bot_u, tau_bot_v, ope, &
                dp, dpp, udp, vdp, ub, vb, Pstress, Pbstress, ubot, vbot, speed

        integer :: I, Iq, ip, I1

        !speed1 = cd_mlswe/gravity

        rhs = 0.0
        tau_bot_u = 0.0 ; tau_bot_v = 0.0
        Pstress = (gravity/alpha_mlswe(1)) * 50.0        ! pressure corresponding to 50m depth
                                                         ! at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth
                                                         ! at which bottom stress is reduced to 0

        !call nodal_computation(Htmp, ope, tau_bot_u, tau_bot_v, source_x, source_y, quu, quv, qvv, qb_df, qprime_df)

        do I = 1, npoin

            dp = qb_df(1,I)
            dpp = qb_df(2,I)
            udp = qb_df(3,I)
            vdp = qb_df(4,I)

            ub = udp/dp; vb = vdp/dp

            wq = wjac_df(I)

            if (botfr == 1) then

                ubot = qprime_df(2,I,nlayers) + ub
                vbot = qprime_df(3,I,nlayers) + vb

                tau_bot_u = (cd_mlswe/alpha_mlswe(nlayers))*ubot
                tau_bot_v = (cd_mlswe/alpha_mlswe(nlayers))*vbot

            elseif (botfr == 2) then

                ubot = qprime_df(2,I,nlayers) + ub
                vbot = qprime_df(3,I,nlayers) + vb
                speed = (cd_mlswe/gravity)*sqrt(ubot**2 + vbot**2)

                tau_bot_u = speed*ubot
                tau_bot_v = speed*vbot
            end if

            source_x = coriolis_df(I)*vdp + gravity*(tau_wind_df(1,I) - tau_bot_u) - &
                        gravity*dp*grad_zbot_df(1,I)

            source_y = -coriolis_df(I)*udp + gravity*(tau_wind_df(2,I) - tau_bot_v) - &
                        gravity*dp*grad_zbot_df(2,I)

            ope = 1.0 + dpp * one_over_pbprime_df(I)
            Hq = (ope**2) * H_bcl_df(I)

            quux = ub * udp + ope * Q_uu_dp_df(I)
            quvxy = ub * vdp + ope * Q_uv_dp_df(I)
            qvvy = vb * vdp + ope * Q_vv_dp_df(I)

            H_ave_df(I) = H_ave_df(I) + Hq;  Qu_ave_df(I) = Qu_ave_df(I) + quux
            Qv_ave_df(I) = Qv_ave_df(I) + qvvy;  Quv_ave_df(I) = Quv_ave_df(I) + quvxy
            tau_bot_ave_df(1,I) = tau_bot_ave_df(1,I) + tau_bot_u
            tau_bot_ave_df(2,I) = tau_bot_ave_df(2,I) + tau_bot_v
            !ope_ave_df(I) = ope_ave_df(I) + ope
            btp_mass_flux_ave_df(1,I) = btp_mass_flux_ave_df(1,I) + udp
            btp_mass_flux_ave_df(2,I) = btp_mass_flux_ave_df(2,I) + vdp

            do ip = 1,npts

                I1 = index_df(ip,I)

                hi = psih_df(ip,I)

                !Xi derivatives
                dhdx = dpsidx_df(ip,I)
                !Eta derivatives
                dhdy = dpsidy_df(ip,I)

                rhs(1,I1) = rhs(1,I1) + wq*(dhdx*udp + dhdy*vdp)
                rhs(2,I1) = rhs(2,I1) + wq*(hi*source_x + dhdx*(Hq + quux) + quvxy*dhdy)
                rhs(3,I1) = rhs(3,I1) + wq*(hi*source_y + dhdx*quvxy + dhdy*(Hq + qvvy))

            end do
        end do

    end subroutine create_rhs_btp_volume_qdf3

    subroutine nodal_computation(Htmp, ope, tau_bot_u, tau_bot_v, source_x, source_y, quu, quv, qvv, qb_df, qprime_df)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_df, tau_wind_df, coriolis_df, one_over_pbprime_df, alpha_mlswe
        use mod_variables, only: Q_uu_dp_df, Q_uv_dp_df, Q_vv_dp_df, H_bcl_df
        use mod_input, only: nlayers, cd_mlswe, botfr

        implicit none

        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(npoin), intent(out) :: Htmp, ope, tau_bot_u, tau_bot_v, source_x, source_y, &
                                                quu, quv, qvv

        real :: dp, dpp, udp, vdp, ub, vb, Pstress, Pbstress, ubot, vbot, speed

        integer :: I, Iq, ip

        !speed1 = cd_mlswe/gravity

        tau_bot_u = 0.0 ; tau_bot_v = 0.0
        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth at which bottom stress is reduced to 0

        do I = 1, npoin

            dp = qb_df(1,I)
            dpp = qb_df(2,I)
            udp = qb_df(3,I)
            vdp = qb_df(4,I)

            ub = udp/dp; vb = vdp/dp

            if (botfr == 1) then

                ubot = qprime_df(2,I,nlayers) + ub
                vbot = qprime_df(3,I,nlayers) + vb

                tau_bot_u(I) = (cd_mlswe/alpha_mlswe(nlayers))*ubot
                tau_bot_v(I) = (cd_mlswe/alpha_mlswe(nlayers))*vbot

            elseif (botfr == 2) then

                ubot = qprime_df(2,I,nlayers) + ub
                vbot = qprime_df(3,I,nlayers) + vb
                speed = (cd_mlswe/gravity)*sqrt(ubot**2 + vbot**2)

                tau_bot_u(I) = speed*ubot
                tau_bot_v(I) = speed*vbot
            end if

            source_x(I) = coriolis_df(I)*vdp + gravity*(tau_wind_df(1,I) - tau_bot_u(I)) - &
                        gravity*dp*grad_zbot_df(1,I)

            source_y(I) = -coriolis_df(I)*udp + gravity*(tau_wind_df(2,I) - tau_bot_v(I)) - &
                        gravity*dp*grad_zbot_df(2,I)

            ope(I) = 1.0 + dpp * one_over_pbprime_df(I)
            Htmp(I) = (ope(I)**2) * H_bcl_df(I)

            quu(I) = ub * udp + ope(I) * Q_uu_dp_df(I)
            quv(I) = ub * vdp + ope(I) * Q_uv_dp_df(I)
            qvv(I) = vb * vdp + ope(I) * Q_vv_dp_df(I)

       end do

    end subroutine nodal_computation

    subroutine creat_btp_fluxes_qdf(rhs,qb_df_face)

        use mod_basis, only: psiq, ngl
        use mod_grid, only:  npoin, intma, nface,face
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: H_face_ave,ope_face_ave,btp_mass_flux_face_ave,Qu_face_ave, Qv_face_ave, &
                                one_plus_eta_edge_2_ave, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, H_bcl_edge, uvb_face_ave

        use mod_initial, only: coeff_pbpert_L, coeff_pbub_LR, coeff_pbpert_R, &
                                one_over_pbprime_face, one_over_pbprime_edge, &
                                coeff_mass_pbub_L, coeff_mass_pbub_R, coeff_mass_pbpert_LR

        implicit none

        real, dimension(3, npoin), intent(inout)  :: rhs
        real, dimension(4, 2, ngl, nface), intent(in) :: qb_df_face
    
        integer :: iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_x, flux_y, flux, nxr, nyr, un
        real, dimension(nq) :: quu, quv, qvu, qvv, H_face_temp, flux_edge_x, flux_edge_y, one_plus_eta_edge
        real, dimension(nq) :: ul, ur, vl, vr
        real :: pU_L, pU_R, lamb, dispu, dispv, pbpert_edge
        real :: qbl(4,nq), qbr(4,nq)

        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            qbl = 0.0; qbr = 0.0

            do iquad = 1,nq

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)
        
                nxr = -nxl
                nyr = -nyl

                do n = 1, ngl

                    hi = psiq(n,iquad)

                    qbl(1:4,iquad) = qbl(1:4,iquad) + hi*qb_df_face(1:4,1,n,iface)

                    qbr(1:4,iquad) = qbr(1:4,iquad) + hi*qb_df_face(1:4,2,n,iface)

                end do 
        
                pU_L = nxl * qbl(3,iquad) + nyl * qbl(4,iquad)
                pU_R = nxr * qbr(3,iquad) + nyr * qbr(4,iquad)
        
                pbpert_edge = coeff_pbpert_L(iquad, iface) * qbl(2,iquad) &
                                            + coeff_pbpert_R(iquad, iface) * qbr(2,iquad) &
                                            + coeff_pbub_LR(iquad, iface) * (pU_L + pU_R)

                !one_plus_eta_edge(iquad) = 1.0 + pbpert_edge * one_over_pbprime_edge(iquad, iface)
                one_plus_eta_edge(iquad) = 1.0 + pbpert_edge * one_over_pbprime_face(1,n,iface)

                ! Compute mass fluxes at each element face.

                flux_edge_x(iquad) = coeff_mass_pbub_L(iquad,iface) * qbl(3,iquad) &
                                            + coeff_mass_pbub_R(iquad,iface) * qbr(3,iquad) &
                                            + coeff_mass_pbpert_LR(iquad,iface) * (nxl * qbl(2,iquad) + nxr * qbr(2,iquad))

                flux_edge_y(iquad) = coeff_mass_pbub_L(iquad,iface) * qbl(4,iquad) &
                                            + coeff_mass_pbub_R(iquad,iface) * qbr(4,iquad) &
                                            + coeff_mass_pbpert_LR(iquad,iface) * (nyl * qbl(2,iquad) + nyr * qbr(2,iquad))

            end do

            ul = qbl(3,:)/qbl(1,:); ur = qbr(3,:)/qbr(1,:)
            vl = qbl(4,:)/qbl(1,:); vr = qbr(4,:)/qbr(1,:)

            quu(:) = 0.5*(ul*qbl(3,:) + ur*qbr(3,:)) + one_plus_eta_edge(:) * Q_uu_dp_edge(:, iface)
            quv(:) = 0.5*(vl*qbl(3,:) + vr*qbr(3,:)) + one_plus_eta_edge(:) * Q_uv_dp_edge(:, iface)

            qvu(:) = 0.5*(ul*qbl(4,:) + ur*qbr(4,:)) + one_plus_eta_edge(:) * Q_uv_dp_edge(:, iface)
            qvv(:) = 0.5*(vl*qbl(4,:) + vr*qbr(4,:)) + one_plus_eta_edge(:) * Q_vv_dp_edge(:, iface)

            ! Compute pressure forcing H_face at each element face.

            H_face_temp(:) = (one_plus_eta_edge(:)**2) * H_bcl_edge(:,iface)

            ! Accumulate sums for time averaging

            btp_mass_flux_face_ave(1,:,iface) = btp_mass_flux_face_ave(1,:,iface) + flux_edge_x(:)
            btp_mass_flux_face_ave(2,:,iface) = btp_mass_flux_face_ave(2,:,iface) + flux_edge_y(:)

            H_face_ave(:,iface) = H_face_ave(:,iface) + H_face_temp(:)
            Qu_face_ave(1,:,iface) = Qu_face_ave(1,:,iface) + quu(:); Qu_face_ave(2,:,iface) = Qu_face_ave(2,:,iface) + quv(:)
            Qv_face_ave(1,:,iface) = Qv_face_ave(1,:,iface) + qvu(:); Qv_face_ave(2,:,iface) = Qv_face_ave(2,:,iface) + qvv(:)

            ope_face_ave(1,:,iface) = ope_face_ave(1,:,iface) + 1.0 + qbl(2,:) * one_over_pbprime_face(1,:,iface)
            ope_face_ave(2,:,iface) = ope_face_ave(2,:,iface) + 1.0 + qbr(2,:) * one_over_pbprime_face(2,:,iface)

            one_plus_eta_edge_2_ave(:,iface) = one_plus_eta_edge_2_ave(:,iface) + one_plus_eta_edge(:)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                H_kx = nxl*H_face_temp(iquad)
                H_ky = nyl*H_face_temp(iquad) 

                lamb = coeff_mass_pbpert_LR(iquad,iface)

                dispu = 0.5*lamb*(qbr(3,iquad) - qbl(3,iquad))
                dispv = 0.5*lamb*(qbr(4,iquad) - qbl(4,iquad))

                flux_x = nxl*quu(iquad) + nyl*quv(iquad) - dispu
                flux_y = nxl*qvu(iquad) + nyl*qvv(iquad) - dispv

                flux = nxl*flux_edge_x(iquad) + nyl*flux_edge_y(iquad)

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    rhs(1,I) = rhs(1,I) - wq*hi*flux
                    rhs(2,I) = rhs(2,I) - wq*hi*(H_kx + flux_x)
                    rhs(3,I) = rhs(3,I) - wq*hi*(H_ky + flux_y)
                    
                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        rhs(1,I) = rhs(1,I) + wq*hi*flux
                        rhs(2,I) = rhs(2,I) + wq*hi*(H_kx + flux_x)
                        rhs(3,I) = rhs(3,I) + wq*hi*(H_ky + flux_y)

                    end if

                end do
            end do

        end do

        rhs(1,:) = massinv(:)*rhs(1,:)
        rhs(2,:) = massinv(:)*rhs(2,:)
        rhs(3,:) = massinv(:)*rhs(3,:)
    
    end subroutine creat_btp_fluxes_qdf

    
    subroutine creat_btp_fluxes_qdf3(rhs,qb_df_face)

        use mod_basis, only: psi, ngl
        use mod_grid, only:  npoin, intma, nface,face
        use mod_face, only: imapl, imapr, normal_vector, jac_face
        use mod_variables, only: H_face_ave_df,ope_face_ave_df,btp_mass_flux_face_ave_df,Qu_face_ave_df, Qv_face_ave_df, &
                                one_plus_eta_edge_2_ave_df, Q_uu_dp_edge_df, Q_uv_dp_edge_df, Q_vv_dp_edge_df, &
                                H_bcl_edge_df

        use mod_initial, only: coeff_pbpert_L_df, coeff_pbub_LR_df, coeff_pbpert_R_df, &
                                one_over_pbprime_df_face, one_over_pbprime_edge, &
                                coeff_mass_pbub_L_df, coeff_mass_pbub_R_df, coeff_mass_pbpert_LR_df, &
                                pbprime_df, zbot_df
        use mod_constants, only: gravity

        implicit none

        real, dimension(3, npoin), intent(inout)  :: rhs
        real, dimension(4, 2, ngl, nface), intent(in) :: qb_df_face

        integer :: iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n, m
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_x, flux_y, flux, nxr, nyr, un
        real, dimension(ngl) :: quu, quv, qvu, qvv, H_face_temp, flux_edge_x, flux_edge_y, one_plus_eta_edge
        real, dimension(ngl) :: ul, ur, vl, vr
        real :: pU_L, pU_R, lamb, dispu, dispv, pbpert_edge
        real :: qbl(4,ngl), qbr(4,ngl), H

        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            qbl = 0.0; qbr = 0.0

            do n = 1,ngl

                nxl = normal_vector(1,n,1,iface)
                nyl = normal_vector(2,n,1,iface)

                nxr = -nxl ; nyr = -nyl

                qbl(1:4,n) = qb_df_face(1:4,1,n,iface)
                qbr(1:4,n) = qb_df_face(1:4,2,n,iface)

                pU_L = nxl * qbl(3,n) + nyl * qbl(4,n)
                pU_R = nxr * qbr(3,n) + nyr * qbr(4,n)

                pbpert_edge = coeff_pbpert_L_df(n, iface) * qbl(2,n) &
                                            + coeff_pbpert_R_df(n, iface) * qbr(2,n) &
                                            + coeff_pbub_LR_df(n, iface) * (pU_L + pU_R)

                one_plus_eta_edge(n) = 1.0 + pbpert_edge * one_over_pbprime_df_face(1,n,iface)

                ! Compute mass fluxes at each element face.

                flux_edge_x(n) = coeff_mass_pbub_L_df(n,iface) * qbl(3,n) &
                                            + coeff_mass_pbub_R_df(n,iface) * qbr(3,n) &
                                            + coeff_mass_pbpert_LR_df(n,iface) * &
                                              (nxl * qbl(2,n) + nxr * qbr(2,n))

                flux_edge_y(n) = coeff_mass_pbub_L_df(n,iface) * qbl(4,n) &
                                            + coeff_mass_pbub_R_df(n,iface) * qbr(4,n) &
                                            + coeff_mass_pbpert_LR_df(n,iface) * &
                                              (nyl * qbl(2,n) + nyr * qbr(2,n))

            end do

            ul = qbl(3,:)/qbl(1,:); ur = qbr(3,:)/qbr(1,:)
            vl = qbl(4,:)/qbl(1,:); vr = qbr(4,:)/qbr(1,:)

            quu(:) = 0.5*(ul*qbl(3,:) + ur*qbr(3,:)) + one_plus_eta_edge(:) * Q_uu_dp_edge_df(:, iface)
            quv(:) = 0.5*(vl*qbl(3,:) + vr*qbr(3,:)) + one_plus_eta_edge(:) * Q_uv_dp_edge_df(:, iface)

            qvu(:) = 0.5*(ul*qbl(4,:) + ur*qbr(4,:)) + one_plus_eta_edge(:) * Q_uv_dp_edge_df(:, iface)
            qvv(:) = 0.5*(vl*qbl(4,:) + vr*qbr(4,:)) + one_plus_eta_edge(:) * Q_vv_dp_edge_df(:, iface)

            ! Compute pressure forcing H_face at each element face.

            H_face_temp(:) = (one_plus_eta_edge(:)**2) * H_bcl_edge_df(:,iface)

            ! Accumulate sums for time averaging

            btp_mass_flux_face_ave_df(1,:,iface) = btp_mass_flux_face_ave_df(1,:,iface) + flux_edge_x(:)
            btp_mass_flux_face_ave_df(2,:,iface) = btp_mass_flux_face_ave_df(2,:,iface) + flux_edge_y(:)

            H_face_ave_df(:,iface) = H_face_ave_df(:,iface) + H_face_temp(:)
            Qu_face_ave_df(1,:,iface) = Qu_face_ave_df(1,:,iface) + quu(:);
            Qu_face_ave_df(2,:,iface) = Qu_face_ave_df(2,:,iface) + quv(:)
            Qv_face_ave_df(1,:,iface) = Qv_face_ave_df(1,:,iface) + qvu(:);
            Qv_face_ave_df(2,:,iface) = Qv_face_ave_df(2,:,iface) + qvv(:)

            ope_face_ave_df(1,:,iface) = ope_face_ave_df(1,:,iface) + (1.0 + qbl(2,:) * one_over_pbprime_df_face(1,:,iface))
            ope_face_ave_df(2,:,iface) = ope_face_ave_df(2,:,iface) + (1.0 + qbr(2,:) * one_over_pbprime_df_face(2,:,iface))

            one_plus_eta_edge_2_ave_df(:,iface) = one_plus_eta_edge_2_ave_df(:,iface) + one_plus_eta_edge(:)

            do n = 1, ngl

                wq = jac_face(n,1,iface)

                nxl = normal_vector(1,n,1,iface)
                nyl = normal_vector(2,n,1,iface)

                H_kx = nxl*H_face_temp(n)
                H_ky = nyl*H_face_temp(n)

                flux = nxl*flux_edge_x(n) + nyl*flux_edge_y(n)
                
                il = imapl(1,n,1,iface)
                jl = imapl(2,n,1,iface)
                kl = imapl(3,n,1,iface)

                I = intma(il,jl,kl,el)
                H = abs(zbot_df(I))

                lamb = coeff_mass_pbpert_LR_df(n,iface)
                !lamb = sqrt(gravity*H)
                dispu = 0.5*lamb*(qbr(3,n) - qbl(3,n))
                dispv = 0.5*lamb*(qbr(4,n) - qbl(4,n))

                flux_x = nxl*quu(n) + nyl*quv(n) !- dispu
                flux_y = nxl*qvu(n) + nyl*qvv(n) !- dispv

                rhs(1,I) = rhs(1,I) - wq*flux
                rhs(2,I) = rhs(2,I) - wq*(H_kx + flux_x)
                rhs(3,I) = rhs(3,I) - wq*(H_ky + flux_y)

                if(er > 0) then

                    ir = imapr(1,n,1,iface)
                    jr = imapr(2,n,1,iface)
                    kr = imapr(3,n,1,iface)

                    I = intma(ir,jr,kr,er)

                    rhs(1,I) = rhs(1,I) + wq*flux
                    rhs(2,I) = rhs(2,I) + wq*(H_kx + flux_x)
                    rhs(3,I) = rhs(3,I) + wq*(H_ky + flux_y)

                end if

            end do
        end do

        rhs(1,:) = massinv(:)*rhs(1,:)
        rhs(2,:) = massinv(:)*rhs(2,:)
        rhs(3,:) = massinv(:)*rhs(3,:)

    end subroutine creat_btp_fluxes_qdf3

     subroutine creat_btp_fluxes_qdf4(rhs,qb_df_face)

        use mod_basis, only: psi, ngl
        use mod_grid, only:  npoin, intma, nface,face
        use mod_face, only: imapl, imapr, normal_vector, jac_face
        use mod_variables, only: H_face_ave_df,ope_face_ave_df,btp_mass_flux_face_ave_df,Qu_face_ave_df, Qv_face_ave_df, &
                                one_plus_eta_edge_2_ave_df, Q_uu_dp_edge_df, Q_uv_dp_edge_df, Q_vv_dp_edge_df, &
                                H_bcl_edge_df

        use mod_initial, only: coeff_pbpert_L_df, coeff_pbub_LR_df, coeff_pbpert_R_df, &
                                one_over_pbprime_df_face, one_over_pbprime_edge, &
                                coeff_mass_pbub_L_df, coeff_mass_pbub_R_df, coeff_mass_pbpert_LR_df

        implicit none

        real, dimension(3, npoin), intent(inout)  :: rhs
        real, dimension(4, 2, ngl, nface), intent(in) :: qb_df_face

        integer :: iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n, m
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_x, flux_y, flux, nxr, nyr, un
        real, dimension(ngl) :: quu, quv, qvu, qvv, H_face_temp, flux_edge_x, flux_edge_y, one_plus_eta_edge
        real, dimension(ngl) :: ul, ur, vl, vr
        real :: pU_L, pU_R, lamb, dispu, dispv, pbpert_edge
        real :: qbl(4,ngl), qbr(4,ngl)

        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            qbl = 0.0; qbr = 0.0

            do n = 1,ngl

                nxl = normal_vector(1,n,1,iface)
                nyl = normal_vector(2,n,1,iface)

                nxr = -nxl ; nyr = -nyl

                qbl(1:4,n) = qb_df_face(1:4,1,n,iface)
                qbr(1:4,n) = qb_df_face(1:4,2,n,iface)

                pU_L = nxl * qbl(3,n) + nyl * qbl(4,n)
                pU_R = nxr * qbr(3,n) + nyr * qbr(4,n)

                pbpert_edge = coeff_pbpert_L_df(n, iface) * qbl(2,n) &
                                            + coeff_pbpert_R_df(n, iface) * qbr(2,n) &
                                            + coeff_pbub_LR_df(n, iface) * (pU_L + pU_R)

                one_plus_eta_edge(n) = 1.0 + pbpert_edge * one_over_pbprime_df_face(1,n,iface)

                ! Compute mass fluxes at each element face.

                flux_edge_x(n) = coeff_mass_pbub_L_df(n,iface) * qbl(3,n) &
                                            + coeff_mass_pbub_R_df(n,iface) * qbr(3,n) &
                                            + coeff_mass_pbpert_LR_df(n,iface) * &
                                              (nxl * qbl(2,n) + nxr * qbr(2,n))

                flux_edge_y(n) = coeff_mass_pbub_L_df(n,iface) * qbl(4,n) &
                                            + coeff_mass_pbub_R_df(n,iface) * qbr(4,n) &
                                            + coeff_mass_pbpert_LR_df(n,iface) * &
                                              (nyl * qbl(2,n) + nyr * qbr(2,n))

            end do

            ul = qbl(3,:)/qbl(1,:); ur = qbr(3,:)/qbr(1,:)
            vl = qbl(4,:)/qbl(1,:); vr = qbr(4,:)/qbr(1,:)

            quu(:) = 0.5*(ul*qbl(3,:) + ur*qbr(3,:)) + one_plus_eta_edge(:) * Q_uu_dp_edge_df(:, iface)
            quv(:) = 0.5*(vl*qbl(3,:) + vr*qbr(3,:)) + one_plus_eta_edge(:) * Q_uv_dp_edge_df(:, iface)

            qvu(:) = 0.5*(ul*qbl(4,:) + ur*qbr(4,:)) + one_plus_eta_edge(:) * Q_uv_dp_edge_df(:, iface)
            qvv(:) = 0.5*(vl*qbl(4,:) + vr*qbr(4,:)) + one_plus_eta_edge(:) * Q_vv_dp_edge_df(:, iface)

            ! Compute pressure forcing H_face at each element face.

            H_face_temp(:) = (one_plus_eta_edge(:)**2) * H_bcl_edge_df(:,iface)

            ! Accumulate sums for time averaging

            btp_mass_flux_face_ave_df(1,:,iface) = btp_mass_flux_face_ave_df(1,:,iface) + flux_edge_x(:)
            btp_mass_flux_face_ave_df(2,:,iface) = btp_mass_flux_face_ave_df(2,:,iface) + flux_edge_y(:)

            H_face_ave_df(:,iface) = H_face_ave_df(:,iface) + H_face_temp(:)
            Qu_face_ave_df(1,:,iface) = Qu_face_ave_df(1,:,iface) + quu(:);
            Qu_face_ave_df(2,:,iface) = Qu_face_ave_df(2,:,iface) + quv(:)
            Qv_face_ave_df(1,:,iface) = Qv_face_ave_df(1,:,iface) + qvu(:);
            Qv_face_ave_df(2,:,iface) = Qv_face_ave_df(2,:,iface) + qvv(:)

            ope_face_ave_df(1,:,iface) = ope_face_ave_df(1,:,iface) + (1.0 + qbl(2,:) * one_over_pbprime_df_face(1,:,iface))
            ope_face_ave_df(2,:,iface) = ope_face_ave_df(2,:,iface) + (1.0 + qbr(2,:) * one_over_pbprime_df_face(2,:,iface))

            one_plus_eta_edge_2_ave_df(:,iface) = one_plus_eta_edge_2_ave_df(:,iface) + one_plus_eta_edge(:)

            do n = 1, ngl

                wq = jac_face(n,1,iface)

                nxl = normal_vector(1,n,1,iface)
                nyl = normal_vector(2,n,1,iface)

                H_kx = nxl*H_face_temp(n)
                H_ky = nyl*H_face_temp(n)

                lamb = coeff_mass_pbpert_LR_df(n,iface)

                dispu = 0.5*lamb*(qbr(3,n) - qbl(3,n))
                dispv = 0.5*lamb*(qbr(4,n) - qbl(4,n))

                flux_x = nxl*quu(n) + nyl*quv(n) !- dispu
                flux_y = nxl*qvu(n) + nyl*qvv(n) !- dispv

                flux = nxl*flux_edge_x(n) + nyl*flux_edge_y(n)

                !do m = 1, ngl

                hi = 1.0 !psi(m,n)
                m = n

                il = imapl(1,m,1,iface)
                jl = imapl(2,m,1,iface)
                kl = imapl(3,m,1,iface)

                I = intma(il,jl,kl,el)

                rhs(1,I) = rhs(1,I) - wq*hi*flux
                rhs(2,I) = rhs(2,I) - wq*hi*(H_kx + flux_x)
                rhs(3,I) = rhs(3,I) - wq*hi*(H_ky + flux_y)

                if(er > 0) then

                    ir = imapr(1,m,1,iface)
                    jr = imapr(2,m,1,iface)
                    kr = imapr(3,m,1,iface)

                    I = intma(ir,jr,kr,er)

                    rhs(1,I) = rhs(1,I) + wq*hi*flux
                    rhs(2,I) = rhs(2,I) + wq*hi*(H_kx + flux_x)
                    rhs(3,I) = rhs(3,I) + wq*hi*(H_ky + flux_y)

                end if

                !end do

            end do

        end do

        rhs(1,:) = massinv(:)*rhs(1,:)
        rhs(2,:) = massinv(:)*rhs(2,:)
        rhs(3,:) = massinv(:)*rhs(3,:)

    end subroutine creat_btp_fluxes_qdf4

     subroutine creat_btp_fluxes_qdf5(rhs,qb_df_face)

        use mod_basis, only: psiq, ngl
        use mod_grid, only:  npoin, intma, nface,face
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: H_face_ave,ope_face_ave,btp_mass_flux_face_ave,Qu_face_ave, Qv_face_ave, &
                                one_plus_eta_edge_2_ave, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, H_bcl_edge, uvb_face_ave

        use mod_initial, only: coeff_pbpert_L, coeff_pbub_LR, coeff_pbpert_R, &
                                one_over_pbprime_face, one_over_pbprime_edge, &
                                coeff_mass_pbub_L, coeff_mass_pbub_R, coeff_mass_pbpert_LR

        implicit none

        real, dimension(3, npoin), intent(inout)  :: rhs
        real, dimension(4, 2, ngl, nface), intent(in) :: qb_df_face

        integer :: iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_x, flux_y, flux, nxr, nyr, un
        real, dimension(nq) :: quu, quv, qvu, qvv, H_face_temp, flux_edge_x, flux_edge_y, one_plus_eta_edge
        real, dimension(nq) :: ul, ur, vl, vr
        real :: pU_L, pU_R, lamb, dispu, dispv, pbpert_edge
        real :: qbl(4,nq), qbr(4,nq)

        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            qbl = 0.0; qbr = 0.0

            do iquad = 1,nq

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                nxr = -nxl
                nyr = -nyl

                do n = 1, ngl

                    hi = psiq(n,iquad)

                    qbl(1:4,iquad) = qbl(1:4,iquad) + hi*qb_df_face(1:4,1,n,iface)

                    qbr(1:4,iquad) = qbr(1:4,iquad) + hi*qb_df_face(1:4,2,n,iface)

                end do

                pU_L = nxl * qbl(3,iquad) + nyl * qbl(4,iquad)
                pU_R = nxr * qbr(3,iquad) + nyr * qbr(4,iquad)

                pbpert_edge = coeff_pbpert_L(iquad, iface) * qbl(2,iquad) &
                                            + coeff_pbpert_R(iquad, iface) * qbr(2,iquad) &
                                            + coeff_pbub_LR(iquad, iface) * (pU_L + pU_R)

                !one_plus_eta_edge(iquad) = 1.0 + pbpert_edge * one_over_pbprime_edge(iquad, iface)
                one_plus_eta_edge(iquad) = 1.0 + pbpert_edge * one_over_pbprime_face(1,n,iface)

                ! Compute mass fluxes at each element face.

                flux_edge_x(iquad) = coeff_mass_pbub_L(iquad,iface) * qbl(3,iquad) &
                                            + coeff_mass_pbub_R(iquad,iface) * qbr(3,iquad) &
                                            + coeff_mass_pbpert_LR(iquad,iface) * (nxl * qbl(2,iquad) + nxr * qbr(2,iquad))

                flux_edge_y(iquad) = coeff_mass_pbub_L(iquad,iface) * qbl(4,iquad) &
                                            + coeff_mass_pbub_R(iquad,iface) * qbr(4,iquad) &
                                            + coeff_mass_pbpert_LR(iquad,iface) * (nyl * qbl(2,iquad) + nyr * qbr(2,iquad))

            end do

            ul = qbl(3,:)/qbl(1,:); ur = qbr(3,:)/qbr(1,:)
            vl = qbl(4,:)/qbl(1,:); vr = qbr(4,:)/qbr(1,:)

            quu(:) = 0.5*(ul*qbl(3,:) + ur*qbr(3,:)) + one_plus_eta_edge(:) * Q_uu_dp_edge(:, iface)
            quv(:) = 0.5*(vl*qbl(3,:) + vr*qbr(3,:)) + one_plus_eta_edge(:) * Q_uv_dp_edge(:, iface)

            qvu(:) = 0.5*(ul*qbl(4,:) + ur*qbr(4,:)) + one_plus_eta_edge(:) * Q_uv_dp_edge(:, iface)
            qvv(:) = 0.5*(vl*qbl(4,:) + vr*qbr(4,:)) + one_plus_eta_edge(:) * Q_vv_dp_edge(:, iface)

            ! Compute pressure forcing H_face at each element face.

            H_face_temp(:) = (one_plus_eta_edge(:)**2) * H_bcl_edge(:,iface)

            ! Accumulate sums for time averaging

            btp_mass_flux_face_ave(1,:,iface) = btp_mass_flux_face_ave(1,:,iface) + flux_edge_x(:)
            btp_mass_flux_face_ave(2,:,iface) = btp_mass_flux_face_ave(2,:,iface) + flux_edge_y(:)

            H_face_ave(:,iface) = H_face_ave(:,iface) + H_face_temp(:)
            Qu_face_ave(1,:,iface) = Qu_face_ave(1,:,iface) + quu(:); Qu_face_ave(2,:,iface) = Qu_face_ave(2,:,iface) + quv(:)
            Qv_face_ave(1,:,iface) = Qv_face_ave(1,:,iface) + qvu(:); Qv_face_ave(2,:,iface) = Qv_face_ave(2,:,iface) + qvv(:)

            ope_face_ave(1,:,iface) = ope_face_ave(1,:,iface) + 1.0 + qbl(2,:) * one_over_pbprime_face(1,:,iface)
            ope_face_ave(2,:,iface) = ope_face_ave(2,:,iface) + 1.0 + qbr(2,:) * one_over_pbprime_face(2,:,iface)

            one_plus_eta_edge_2_ave(:,iface) = one_plus_eta_edge_2_ave(:,iface) + one_plus_eta_edge(:)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                H_kx = nxl*H_face_temp(iquad)
                H_ky = nyl*H_face_temp(iquad)

                lamb = coeff_mass_pbpert_LR(iquad,iface)

                dispu = 0.5*lamb*(qbr(3,iquad) - qbl(3,iquad))
                dispv = 0.5*lamb*(qbr(4,iquad) - qbl(4,iquad))

                flux_x = nxl*quu(iquad) + nyl*quv(iquad) - dispu
                flux_y = nxl*qvu(iquad) + nyl*qvv(iquad) - dispv

                flux = nxl*flux_edge_x(iquad) + nyl*flux_edge_y(iquad)

                do n = 1, ngl

                    hi = psiq(n,iquad)

                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    !rhs(1,I) = rhs(1,I) - wq*hi*flux
                    !rhs(2,I) = rhs(2,I) - wq*hi*(H_kx + flux_x)
                    !rhs(3,I) = rhs(3,I) - wq*hi*(H_ky + flux_y)

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        !rhs(1,I) = rhs(1,I) + wq*hi*flux
                        !rhs(2,I) = rhs(2,I) + wq*hi*(H_kx + flux_x)
                        !rhs(3,I) = rhs(3,I) + wq*hi*(H_ky + flux_y)

                    end if

                end do
            end do

        end do

        !rhs(1,:) = massinv(:)*rhs(1,:)
        !rhs(2,:) = massinv(:)*rhs(2,:)
        !rhs(3,:) = massinv(:)*rhs(3,:)

    end subroutine creat_btp_fluxes_qdf5

    subroutine create_rhs_btp_dynamics_volume(rhs_mom, qb)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, coriolis_quad
        use mod_variables, only: Quu, Quv, Qvv, tau_bot, H

        implicit none 

        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(2,npoin), intent(out) :: rhs_mom

        real :: source_x, source_y, Hq, quux, quvxy, qvvy

        real :: wq, hi, dhdx, dhdy

        integer :: I, Iq, ip

        rhs_mom = 0.0

        do Iq = 1, npoin_q

            source_x = gravity*(tau_wind(1,Iq) - tau_bot(1,Iq)) - &
                        gravity*qb(1,Iq)*grad_zbot_quad(1,Iq)

            source_y = gravity*(tau_wind(2,Iq) - tau_bot(2,Iq)) - &
                        gravity*qb(1,Iq)*grad_zbot_quad(2,Iq)

            Hq = H(Iq)
            quux = Quu(Iq)
            quvxy = Quv(Iq)
            qvvy = Qvv(Iq)

            wq = wjac(Iq)

            do ip = 1,npts

                I = indexq(ip,Iq)

                hi = psih(ip,Iq)

                !Xi derivatives
                dhdx = dpsidx(ip,Iq)
                !Eta derivatives
                dhdy = dpsidy(ip,Iq)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*(hi*source_x + dhdx*Hq + quux*dhdx + quvxy*dhdy)
                rhs_mom(2,I) = rhs_mom(2,I) + wq*(hi*source_y + dhdy*Hq + quvxy*dhdx + qvvy*dhdy)

            end do
        end do

    end subroutine create_rhs_btp_dynamics_volume

    subroutine create_rhs_btp_dynamics_volume_all(rhs,qb)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, coriolis_quad
        use mod_variables, only: btp_mass_flux, tau_bot, H, Quu, Quv, Qvv

        implicit none 

        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(3,npoin), intent(out) :: rhs

        real :: source_x, source_y, Hq, quux, quvxy, qvvy
        real :: wq, hi, dhdx, dhdy, var_u, var_v

        integer :: I, Iq, ip

        rhs = 0.0

        do Iq = 1, npoin_q

            wq = wjac(Iq)
            var_u = btp_mass_flux(1,Iq)
            var_v = btp_mass_flux(2,Iq)

            source_x = gravity*(tau_wind(1,Iq) - tau_bot(1,Iq)) - &
                        gravity*qb(1,Iq)*grad_zbot_quad(1,Iq)

            source_y = gravity*(tau_wind(2,Iq) - tau_bot(2,Iq)) - &
                        gravity*qb(1,Iq)*grad_zbot_quad(2,Iq)

            Hq = H(Iq)
            quux = Quu(Iq)
            quvxy = Quv(Iq)
            qvvy = Qvv(Iq)

            do ip = 1,npts

                I = indexq(ip,Iq)

                hi = psih(ip,Iq)

                !Xi derivatives
                dhdx = dpsidx(ip,Iq)
                !Eta derivatives
                dhdy = dpsidy(ip,Iq)

                rhs(1,I) = rhs(1,I) + wq*(dhdx*var_u + dhdy*var_v)
                rhs(2,I) = rhs(2,I) + wq*(hi*source_x + dhdx*Hq + quux*dhdx + quvxy*dhdy)
                rhs(3,I) = rhs(3,I) + wq*(hi*source_y + dhdy*Hq + quvxy*dhdx + qvvy*dhdy)

            end do
        end do

    end subroutine create_rhs_btp_dynamics_volume_all

    subroutine Apply_btp_fluxes(rhs_mom,qb_face)

        use mod_basis, only: nq, psiq, ngl
        use mod_grid, only:  npoin, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_input, only: nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe, coeff_mass_pbpert_LR, pbprime_face
        use mod_variables, only: H_face,Qu_face,Qv_face

        implicit none

        real, dimension(2, npoin), intent(inout)  :: rhs_mom
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face
    
        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, nq_i, nq_j
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_x, flux_y
        real, dimension(nq) :: quu, quv, qvu, qvv, Hface_l
        integer :: n, m
        real :: unl, unr, hbl, hbr, lambl, lambr, lamb, dispu, dispv, ubl, vbl, ubr, vbr

        quu = 0.0
        quv = 0.0
        qvu = 0.0
        qvv = 0.0
        Hface_l = 0.0

        do iface = 1, nface

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            Hface_l(:) = H_face(:,iface)

            quu(:) = Qu_face(1,:,iface)
            quv(:) = Qu_face(2,:,iface)

            qvu(:) = Qv_face(1,:,iface)
            qvv(:) = Qv_face(2,:,iface)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                H_kx = nxl*Hface_l(iquad)
                H_ky = nyl*Hface_l(iquad) 

                lamb = coeff_mass_pbpert_LR(iquad,iface)

                dispu = 0.5*lamb*(qb_face(3,2,iquad,iface) - qb_face(3,1,iquad,iface))
                dispv = 0.5*lamb*(qb_face(4,2,iquad,iface) - qb_face(4,1,iquad,iface))

                flux_x = nxl*quu(iquad) + nyl*quv(iquad) - dispu
                flux_y = nxl*qvu(iquad) + nyl*qvv(iquad) - dispv

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    rhs_mom(1,I) = rhs_mom(1,I) - wq*hi*(H_kx + flux_x)
                    rhs_mom(2,I) = rhs_mom(2,I) - wq*hi*(H_ky + flux_y)

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        rhs_mom(1,I) = rhs_mom(1,I) + wq*hi*(H_kx + flux_x)
                        rhs_mom(2,I) = rhs_mom(2,I) + wq*hi*(H_ky + flux_y)

                    end if

                end do
            end do

        end do

        rhs_mom(1,:) = massinv(:)*rhs_mom(1,:)
        rhs_mom(2,:) = massinv(:)*rhs_mom(2,:)
    
    end subroutine Apply_btp_fluxes

    subroutine Apply_btp_fluxes_all(rhs,qb_face)

        use mod_basis, only: nq, psiq, ngl
        use mod_grid, only:  npoin, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_input, only: nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe, coeff_mass_pbpert_LR, pbprime_face
        use mod_variables, only: H_face,Qu_face,Qv_face,flux_edge

        implicit none

        real, dimension(3, npoin), intent(inout)  :: rhs
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face
    
        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, nq_i, nq_j
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_x, flux_y, flux
        real, dimension(nq) :: quu, quv, qvu, qvv, Hface_l, var_x, var_y
        integer :: n, m
        real :: unl, unr, hbl, hbr, lambl, lambr, lamb, dispu, dispv, ubl, vbl, ubr, vbr

        quu = 0.0
        quv = 0.0
        qvu = 0.0
        qvv = 0.0
        Hface_l = 0.0

        do iface = 1, nface

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            Hface_l(:) = H_face(:,iface)

            quu(:) = Qu_face(1,:,iface)
            quv(:) = Qu_face(2,:,iface)

            qvu(:) = Qv_face(1,:,iface)
            qvv(:) = Qv_face(2,:,iface)

            var_x = flux_edge(1,:,iface)
            var_y = flux_edge(2,:,iface)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                H_kx = nxl*Hface_l(iquad)
                H_ky = nyl*Hface_l(iquad) 

                lamb = coeff_mass_pbpert_LR(iquad,iface)

                dispu = 0.5*lamb*(qb_face(3,2,iquad,iface) - qb_face(3,1,iquad,iface))
                dispv = 0.5*lamb*(qb_face(4,2,iquad,iface) - qb_face(4,1,iquad,iface))

                flux_x = nxl*quu(iquad) + nyl*quv(iquad) - dispu
                flux_y = nxl*qvu(iquad) + nyl*qvv(iquad) - dispv

                flux = nxl*var_x(iquad) + nyl*var_y(iquad)

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    rhs(1,I) = rhs(1,I) - wq*hi*flux
                    rhs(2,I) = rhs(2,I) - wq*hi*(H_kx + flux_x)
                    rhs(3,I) = rhs(3,I) - wq*hi*(H_ky + flux_y)

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        rhs(1,I) = rhs(1,I) + wq*hi*flux
                        rhs(2,I) = rhs(2,I) + wq*hi*(H_kx + flux_x)
                        rhs(3,I) = rhs(3,I) + wq*hi*(H_ky + flux_y)

                    end if

                end do
            end do

        end do

        rhs(1,:) = massinv(:)*rhs(1,:)
        rhs(2,:) = massinv(:)*rhs(2,:)
        rhs(3,:) = massinv(:)*rhs(3,:)
    
    end subroutine Apply_btp_fluxes_all

    subroutine create_rhs_btp_dynamics_volume_mass(pb_advec)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac
        use mod_variables, only: btp_mass_flux

        implicit none 

        real, dimension(npoin), intent(out) :: pb_advec


        real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, hi, dhde, dhdn
        real :: h_e, h_c, h_n, var_u, var_v

        integer :: k, e, iquad, jquad, kquad, l, n, m, I, Iq, ip

        pb_advec = 0.0

        do Iq = 1,npoin_q
            
            wq = wjac(Iq)

            var_u = btp_mass_flux(1,Iq)
            var_v = btp_mass_flux(2,Iq)

            do ip = 1,npts

                I = indexq(ip,Iq)

                pb_advec(I) = pb_advec(I) + wq*(dpsidx(ip,Iq)*var_u + dpsidy(ip,Iq)*var_v)

            end do
        end do

    end subroutine create_rhs_btp_dynamics_volume_mass

    subroutine Apply_btp_flux_mass(pb_advec)

        use mod_basis, only: nq, psiq, ngl
        use mod_grid, only:  npoin, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: flux_edge

        implicit none

        real, dimension(npoin), intent(inout) :: pb_advec
    
        integer :: k, iface, iquad, el, er, il, jl, ir, jr, I
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j, n, m
        real :: wq, nxl, nyl, flux
        real, dimension(nq) :: var_x, var_y
        real :: hi
    
        do iface = 1, nface

            !Store Left Side Variables

            el = face(7,iface)
            er = face(8,iface)

            var_x = flux_edge(1,:,iface)
            var_y = flux_edge(2,:,iface)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                flux = nxl*var_x(iquad) + nyl*var_y(iquad)

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    pb_advec(I) = pb_advec(I) - wq*hi*flux

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        pb_advec(I) = pb_advec(I) + wq*hi*flux

                    end if !er > 0

                end do !n
            end do !iquad

        end do

        pb_advec(:) = massinv(:)*pb_advec(:)
    
    end subroutine Apply_btp_flux_mass
    

end module mod_rhs_btp
