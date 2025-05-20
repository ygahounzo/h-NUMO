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

    subroutine create_rhs_btp(rhs,qb_df,qprime_df)

        implicit none

        real, dimension(3, npoin), intent(out) :: rhs
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        real, dimension(4, 2, ngl, nface) :: qb_df_face
        real, dimension(2,npoin) :: rhs_visc_btp

        call btp_extract_df(qb_df_face, qb_df)

        call btp_create_precommunicator(qb_df_face,4)

        call create_rhs_btp_volume_qdf(rhs, qb_df, qprime_df)

        call btp_create_postcommunicator(qb_df_face,4)

        call creat_btp_fluxes_qdf(rhs,qb_df_face)

        ! Compute RHS viscosity terms

        if(method_visc > 0) call btp_create_laplacian(rhs_visc_btp,qb_df)

        rhs(2,:) = rhs(2,:) + rhs_visc_btp(1,:)
        rhs(3,:) = rhs(3,:) + rhs_visc_btp(2,:)

    end subroutine create_rhs_btp

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

    subroutine create_rhs_btp_volume_qdf(rhs, qb_df, qprime_df)

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
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(3,npoin), intent(out) :: rhs

        real :: source_x, source_y, Hq, quux, quvxy, qvvy
        real :: wq, hi, dhdx, dhdy, coef_fric, tau_bot_u, tau_bot_v, ope, &
                dp, dpp, udp, vdp, ub, vb, Pstress, Pbstress, ubot, vbot, speed
        real :: pp, up, vp

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
            pp = 0.0; up = 0.0; vp = 0.0

            do ip = 1,npts 
                I = indexq(ip,Iq)
                hi = psih(ip,Iq)

                dp = dp + hi*qb_df(1,I)
                dpp = dpp + hi*qb_df(2,I)
                udp = udp + hi*qb_df(3,I)
                vdp = vdp + hi*qb_df(4,I)
                pp = pp + hi*qprime_df(1,I,nlayers)
                up = up + hi*qprime_df(2,I,nlayers)
                vp = vp + hi*qprime_df(3,I,nlayers)
            end do 

            wq = wjac(Iq)
            ub = udp/dp; vb = vdp/dp

            if (botfr == 1) then
                ubot = up + ub
                vbot = vp + vb
                speed = (cd_mlswe/gravity)*pp
                tau_bot_u = speed*ubot
                tau_bot_v = speed*vbot
            elseif (botfr == 2) then
                ubot = up + ub
                vbot = vp + vb
                speed = (cd_mlswe/alpha_mlswe(nlayers))*sqrt(ubot**2 + vbot**2)
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

     subroutine creat_btp_fluxes_qdf(rhs,qb_df_face)

        use mod_basis, only: psiq, ngl
        use mod_grid, only:  npoin, intma, nface,face
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: H_face_ave,ope_face_ave,btp_mass_flux_face_ave,Qu_face_ave, Qv_face_ave, &
                                one_plus_eta_edge_2_ave, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, H_bcl_edge

        use mod_initial, only: coeff_pbpert_L, coeff_pbub_LR, coeff_pbpert_R, &
                                one_over_pbprime_face, one_over_pbprime_edge, &
                                coeff_mass_pbub_L, coeff_mass_pbub_R, coeff_mass_pbpert_LR, &
                                alpha_mlswe, pbprime_df_face

        implicit none

        real, dimension(3, npoin), intent(inout)  :: rhs
        real, dimension(4, 2, ngl, nface), intent(in) :: qb_df_face

        integer :: iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_x, flux_y, flux, nxr, nyr, un
        real, dimension(nq) :: quu, quv, qvu, qvv, H_face_temp, flux_edge_x, flux_edge_y, one_plus_eta_edge
        real, dimension(nq) :: ul, ur, vl, vr, pbl, pbr
        real :: pU_L, pU_R, lamb, dispu, dispv, pbpert_edge
        real :: qbl(4,nq), qbr(4,nq), lam1, lam2, dispp

        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            qbl = 0.0; qbr = 0.0
            pbl = 0.0; pbr = 0.0

            do iquad = 1,nq

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                nxr = -nxl; nyr = -nyl

                do n = 1, ngl
                    hi = psiq(n,iquad)
                    qbl(1:4,iquad) = qbl(1:4,iquad) + hi*qb_df_face(1:4,1,n,iface)
                    qbr(1:4,iquad) = qbr(1:4,iquad) + hi*qb_df_face(1:4,2,n,iface)
                    pbl(iquad) = pbl(iquad) + hi*pbprime_df_face(1,n,iface)
                    pbr(iquad) = pbr(iquad) + hi*pbprime_df_face(2,n,iface)
                end do

                pU_L = nxl * qbl(3,iquad) + nyl * qbl(4,iquad)
                pU_R = nxr * qbr(3,iquad) + nyr * qbr(4,iquad)

                pbpert_edge = 0.5*(qbl(2,iquad) + qbr(2,iquad))
                one_plus_eta_edge(iquad) = 1.0 + (pbpert_edge/pbl(iquad))

                ! Compute mass fluxes at each element face.
                flux_edge_x(iquad) = 0.5 * (qbl(3,iquad) + qbr(3,iquad))
                flux_edge_y(iquad) = 0.5 * (qbl(4,iquad) + qbr(4,iquad))

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
            ope_face_ave(1,:,iface) = ope_face_ave(1,:,iface) + 1.0 + (qbl(2,:)/pbl)
            ope_face_ave(2,:,iface) = ope_face_ave(2,:,iface) + 1.0 + (qbr(2,:)/pbr)
            one_plus_eta_edge_2_ave(:,iface) = one_plus_eta_edge_2_ave(:,iface) + one_plus_eta_edge(:)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                H_kx = nxl*H_face_temp(iquad)
                H_ky = nyl*H_face_temp(iquad)

                !lamb = coeff_mass_pbpert_LR(iquad,iface)
                lam1 = abs(nxl*ul(iquad) + nyl*vl(iquad)) + sqrt(alpha_mlswe(nlayers) * pbl(iquad))
                lam2 = abs(nxr*ur(iquad) + nyr*vr(iquad)) + sqrt(alpha_mlswe(nlayers) * pbr(iquad))
                lamb = max(lam1,lam2)

                dispu = 0.5*lamb*(qbr(3,iquad) - qbl(3,iquad))
                dispv = 0.5*lamb*(qbr(4,iquad) - qbl(4,iquad))
                flux_x = nxl*quu(iquad) + nyl*quv(iquad) - dispu
                flux_y = nxl*qvu(iquad) + nyl*qvv(iquad) - dispv

                dispp = 0.5*lamb*(qbr(2,iquad) - qbl(2,iquad))
                flux = nxl*flux_edge_x(iquad) + nyl*flux_edge_y(iquad) - dispp

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
