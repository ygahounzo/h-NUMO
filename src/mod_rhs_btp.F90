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
    use mod_laplacian_quad, only: btp_create_laplacian, btp_create_laplacian_v2
    use mod_barotropic_terms, only: btp_extract_df
    use mod_input, only: nlayers, method_visc
    use mod_metrics, only: ksiq_x, ksiq_y, ksiq_z, &
                           etaq_x, etaq_y, etaq_z, &
                           zetaq_x, zetaq_y, zetaq_z, &
                           jacq, massinv

    public :: create_rhs_btp

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

        if(method_visc == 1) then
            call btp_create_laplacian_v2(rhs_visc_btp, qprime_df, qb_df)
        else
            call btp_create_laplacian(rhs_visc_btp,qb_df)
        endif 

        rhs(2,:) = rhs(2,:) + rhs_visc_btp(1,:)
        rhs(3,:) = rhs(3,:) + rhs_visc_btp(2,:)

    end subroutine create_rhs_btp

    subroutine create_rhs_btp_volume_qdf(rhs, qb_df, qprime_df)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, &
                                pbprime_df, coriolis_quad, alpha_mlswe
        use mod_variables, only: tau_bot_ave, H_ave, Qu_ave, Quv_ave, Qv_ave, ope_ave, &
                                uvb_ave, btp_mass_flux_ave, ope2_ave
        use mod_variables, only: Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl
        use mod_input, only: nlayers, cd_mlswe, botfr

        implicit none

        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(3,npoin), intent(out) :: rhs

        real :: sc_x, sc_y, Hq, qu, quv, qv
        real :: wq, hi, dhdx, dhdy, coef_fric, tb_u, tb_v, ope, &
                dp, dpp, udp, vdp, ub, vb, Pstress, Pbstress, ubot, vbot, spd
        real :: pp, up, vp, pbq

        integer :: I, Iq, ip

        !speed1 = cd_mlswe/gravity

        rhs = 0.0
        tb_u = 0.0 ; tb_v = 0.0
        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth 
                                                  ! at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth
                                                         ! at which bottom stress is reduced to 0

        do Iq = 1, npoin_q

            dp = 0.0; dpp = 0.0; udp = 0.0; vdp = 0.0
            pp = 0.0; up = 0.0; vp = 0.0; pbq = 0.0

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
                pbq = pbq + hi*pbprime_df(I)
            end do

            wq = wjac(Iq)
            ub = udp/dp; vb = vdp/dp

            if (botfr == 1) then
                ubot = up + ub
                vbot = vp + vb
                spd = (cd_mlswe/gravity)*pp
                tb_u = spd*ubot
                tb_v = spd*vbot
            elseif (botfr == 2) then
                ubot = up + ub
                vbot = vp + vb
                spd = (cd_mlswe/alpha_mlswe(nlayers))*sqrt(ubot**2 + vbot**2)
                tb_u = spd*ubot
                tb_v = spd*vbot
            end if

            sc_x = coriolis_quad(Iq)*vdp + gravity*(tau_wind(1,Iq) - tb_u) - &
                        gravity*dp*grad_zbot_quad(1,Iq)
            sc_y = -coriolis_quad(Iq)*udp + gravity*(tau_wind(2,Iq) - tb_v) - &
                        gravity*dp*grad_zbot_quad(2,Iq)

            ope = 1.0 + (dpp / pbq)
            Hq = (ope**2) * H_bcl(Iq)

            qu = ub * udp + ope * Q_uu_dp(Iq)
            quv = ub * vdp + ope * Q_uv_dp(Iq)
            qv = vb * vdp + ope * Q_vv_dp(Iq)

            H_ave(Iq) = H_ave(Iq) + Hq;  Qu_ave(Iq) = Qu_ave(Iq) + qu
            Qv_ave(Iq) = Qv_ave(Iq) + qv;  Quv_ave(Iq) = Quv_ave(Iq) + quv
            tau_bot_ave(1,Iq) = tau_bot_ave(1,Iq) + tb_u
            tau_bot_ave(2,Iq) = tau_bot_ave(2,Iq) + tb_v
            ope_ave(Iq) = ope_ave(Iq) + ope
            ope2_ave(Iq) = ope2_ave(Iq) + ope**2
            btp_mass_flux_ave(1,Iq) = btp_mass_flux_ave(1,Iq) + udp
            btp_mass_flux_ave(2,Iq) = btp_mass_flux_ave(2,Iq) + vdp
            uvb_ave(1,Iq) = uvb_ave(1,Iq) + ub
            uvb_ave(2,Iq) = uvb_ave(2,Iq) + vb

            do ip = 1,npts
                I = indexq(ip,Iq)
                hi = psih(ip,Iq)
                !Xi derivatives
                dhdx = dpsidx(ip,Iq)
                !Eta derivatives
                dhdy = dpsidy(ip,Iq)

                rhs(1,I) = rhs(1,I) + wq*(dhdx*udp + dhdy*vdp)
                rhs(2,I) = rhs(2,I) + wq*(hi*sc_x + dhdx*(Hq + qu) + quv*dhdy)
                rhs(3,I) = rhs(3,I) + wq*(hi*sc_y + dhdx*quv + dhdy*(Hq + qv))

            end do
        end do

    end subroutine create_rhs_btp_volume_qdf

    subroutine creat_btp_fluxes_qdf(rhs,qb_df_face)

        use mod_basis, only: psiq, ngl
        use mod_grid, only:  npoin, intma, nface,face
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: H_face_ave,ope_face_ave,btp_mass_flux_face_ave, &
                                Qu_face_ave, Qv_face_ave, one_plus_eta_edge_2_ave, &
                                Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, &
                                uvb_face_ave, H_bcl_edge, ope2_face_ave
        use mod_initial, only: alpha_mlswe, pbprime_df_face
        use mod_constants, only: gravity
        use mod_input, only: nlayers

        implicit none

        real, dimension(3, npoin), intent(inout)  :: rhs
        real, dimension(4, 2, ngl, nface), intent(in) :: qb_df_face

        integer :: iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_x, flux_y, flux, nxr, nyr, un
        real, dimension(nq) :: quu, quv, qvu, qvv, H_face_temp, flux_edge_x, flux_edge_y
        real, dimension(nq) :: ul, ur, vl, vr, pbl, pbr, lambq, one_plus_eta_edge
        real :: pU_L, pU_R, lamb, dispu, dispv, pbpert_edge
        real :: qbl(4,nq), qbr(4,nq), lam1, lam2, dispp
        real :: c_minus, c_plus, c_pb_L, c_pb_R, c_pbub_LR, c_pbub_L, c_pbub_R
        real, dimension(nq) :: c_pb_LR

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

                c_minus = sqrt(alpha_mlswe(nlayers) * pbr(iquad))
                c_plus  = sqrt(alpha_mlswe(nlayers) * pbl(iquad))
                c_pb_L = c_minus / (c_minus + c_plus)
                c_pb_R = c_plus / (c_minus + c_plus)
                c_pbub_LR  = 1.0 / (c_minus + c_plus)

                pbpert_edge = c_pb_L * qbl(2,iquad) &
                                            + c_pb_R * qbr(2,iquad) &
                                            + c_pbub_LR * (pU_L + pU_R)

                one_plus_eta_edge(iquad) = 1.0 + (pbpert_edge/pbl(iquad))

                ! Compute mass fluxes at each element face.

                c_pbub_L = c_minus / (c_minus + c_plus)
                c_pbub_R = c_plus / (c_minus + c_plus)
                c_pb_LR(iquad) = c_minus * c_plus / (c_minus + c_plus)

                flux_edge_x(iquad) = c_pbub_L * qbl(3,iquad) &
                                            + c_pbub_R * qbr(3,iquad) &
                                            + c_pb_LR(iquad) * &
                                            (nxl * qbl(2,iquad) + nxr * qbr(2,iquad))

                flux_edge_y(iquad) = c_pbub_L * qbl(4,iquad) &
                                            + c_pbub_R * qbr(4,iquad) &
                                            + c_pb_LR(iquad) * &
                                            (nyl * qbl(2,iquad) + nyr * qbr(2,iquad))

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
            Qu_face_ave(1,:,iface) = Qu_face_ave(1,:,iface) + quu(:)
            Qu_face_ave(2,:,iface) = Qu_face_ave(2,:,iface) + quv(:)
            Qv_face_ave(1,:,iface) = Qv_face_ave(1,:,iface) + qvu(:)
            Qv_face_ave(2,:,iface) = Qv_face_ave(2,:,iface) + qvv(:)
            ope_face_ave(1,:,iface) = ope_face_ave(1,:,iface) + (1.0 + (qbl(2,:)/pbl))
            ope_face_ave(2,:,iface) = ope_face_ave(2,:,iface) + (1.0 + (qbr(2,:)/pbr))
            ope2_face_ave(1,:,iface) = ope2_face_ave(1,:,iface) + (1.0 + (qbl(2,:)/pbl))**2
            ope2_face_ave(2,:,iface) = ope2_face_ave(2,:,iface) + (1.0 + (qbr(2,:)/pbr))**2
            one_plus_eta_edge_2_ave(:,iface) = one_plus_eta_edge_2_ave(:,iface) &
                                                + one_plus_eta_edge(:)**2
            uvb_face_ave(1,1,:,iface) = uvb_face_ave(1,1,:,iface) + ul
            uvb_face_ave(1,2,:,iface) = uvb_face_ave(1,2,:,iface) + ur
            uvb_face_ave(2,1,:,iface) = uvb_face_ave(2,1,:,iface) + vl
            uvb_face_ave(2,2,:,iface) = uvb_face_ave(2,2,:,iface) + vr

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                H_kx = nxl*H_face_temp(iquad)
                H_ky = nyl*H_face_temp(iquad)

                lamb = c_pb_LR(iquad)

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
    
end module mod_rhs_btp
