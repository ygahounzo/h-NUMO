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

        call btp_create_precommunicator(qb_df,qprime_df,4)

        call create_rhs_btp_volume_qdf(rhs, qb_df, qprime_df)

        call create_btp_fluxes_qdf(rhs,qb_df, qprime_df)

        call btp_create_postcommunicator(rhs,4)

        ! Compute RHS viscosity terms

        if(method_visc == 1) then
            call btp_create_laplacian_v2(rhs_visc_btp, qprime_df, qb_df)
        elseif (method_visc > 1) then
            call btp_create_laplacian(rhs_visc_btp,qb_df)
        endif 

        rhs(1,:) = massinv(:)*rhs(1,:)
        rhs(2,:) = massinv(:)*(rhs(2,:) + rhs_visc_btp(1,:))
        rhs(3,:) = massinv(:)*(rhs(3,:) + rhs_visc_btp(2,:))

    end subroutine create_rhs_btp

    subroutine create_rhs_btp_volume_qdf(rhs, qb_df, qprime_df)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, &
                                pbprime_df, coriolis_quad, alpha_mlswe
        use mod_variables, only: tau_bot_ave, H_ave, Qu_ave, Quv_ave, Qv_ave, ope_ave, &
                                uvb_ave, btp_mass_flux_ave, ope2_ave
        use mod_input, only: nlayers, cd_mlswe, botfr, dry_cutoff

        implicit none

        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(3,npoin), intent(out) :: rhs

        real :: sc_x, sc_y, Hq, Qu(2), Qv(2)
        real :: wq, hi, dhdx, dhdy, tb_u, tb_v, ope, &
                dp, dpp, udp, vdp, ub, vb, Pstress, Pbstress, ubot, vbot, spd, pbq
        real, dimension(nlayers) :: pp, up, vp
        real, dimension(nlayers+1) :: pprime
        integer :: I, Iq, ip, k

        !speed1 = cd_mlswe/gravity

        rhs = 0.0
        tb_u = 0.0 ; tb_v = 0.0
        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth 
                                                  ! at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth
                                                         ! at which bottom stress is reduced to 0

        do Iq = 1, npoin_q

            dp = 0.0; dpp = 0.0; udp = 0.0; vdp = 0.0; pbq = 0.0

            do ip = 1,npts
                I = indexq(ip,Iq)
                hi = psih(ip,Iq)

                dp = dp + hi*qb_df(1,I)
                dpp = dpp + hi*qb_df(2,I)
                udp = udp + hi*qb_df(3,I)
                vdp = vdp + hi*qb_df(4,I)
                pbq = pbq + hi*pbprime_df(I)
            end do

            Hq = 0.0 ; pprime(:) = 0.0
            pp = 0.0; up = 0.0; vp = 0.0
            
            do k = 1, nlayers
                do ip = 1,npts
                    I = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    pp(k) = pp(k) + hi*qprime_df(1,I,k)
                    up(k) = up(k) + hi*qprime_df(2,I,k)
                    vp(k) = vp(k) + hi*qprime_df(3,I,k)
                enddo

                pprime(k+1) = pprime(k) + pp(k)
                Hq = Hq + 0.5*alpha_mlswe(k)*(pprime(k+1)**2 - pprime(k)**2)

            enddo

            wq = wjac(Iq)
            ub = udp/dp; vb = vdp/dp

            if (botfr == 1) then
                ubot = up(nlayers) + ub
                vbot = vp(nlayers) + vb
                spd = (cd_mlswe/gravity)*pp(nlayers)
                tb_u = spd*ubot
                tb_v = spd*vbot
            elseif (botfr == 2) then
                ubot = up(nlayers) + ub
                vbot = vp(nlayers) + vb
                spd = (cd_mlswe/alpha_mlswe(nlayers))*sqrt(ubot**2 + vbot**2)
                tb_u = spd*ubot
                tb_v = spd*vbot
            end if

            sc_x = coriolis_quad(Iq)*vdp + gravity*(tau_wind(1,Iq) - tb_u) - &
                        gravity*dp*grad_zbot_quad(1,Iq)
            sc_y = -coriolis_quad(Iq)*udp + gravity*(tau_wind(2,Iq) - tb_v) - &
                        gravity*dp*grad_zbot_quad(2,Iq)

            ope = 1.0 + (dpp / pbq)
            Hq = (ope**2) * Hq

            Qu(1) = ub * udp + ope * sum(up(:)*(pp(:) * up(:)))
            Qu(2) = vb * udp + ope * sum(vp(:)*(pp(:) * up(:)))
            Qv(1) = ub * vdp + ope * sum(up(:)*(pp(:) * vp(:)))
            Qv(2) = vb * vdp + ope * sum(vp(:)*(pp(:) * vp(:)))

            H_ave(Iq) = H_ave(Iq) + Hq;  Qu_ave(Iq) = Qu_ave(Iq) + Qu(1)
            Qv_ave(Iq) = Qv_ave(Iq) + Qv(1);  Quv_ave(Iq) = Quv_ave(Iq) + Qu(2)
            tau_bot_ave(1,Iq) = tau_bot_ave(1,Iq) + tb_u
            tau_bot_ave(2,Iq) = tau_bot_ave(2,Iq) + tb_v
            ope_ave(Iq) = ope_ave(Iq) + ope
            ope2_ave(Iq) = ope2_ave(Iq) + ope**2
            btp_mass_flux_ave(1,Iq) = btp_mass_flux_ave(1,Iq) + udp
            btp_mass_flux_ave(2,Iq) = btp_mass_flux_ave(2,Iq) + vdp
            uvb_ave(1,Iq) = uvb_ave(1,Iq) + ub
            uvb_ave(2,Iq) = uvb_ave(2,Iq) + vb

            Qu(1) = Qu(1) + Hq
            Qv(2) = Qv(2) + Hq

            do ip = 1,npts
                I = indexq(ip,Iq)
                hi = psih(ip,Iq)
                !Xi derivatives
                dhdx = dpsidx(ip,Iq)
                !Eta derivatives
                dhdy = dpsidy(ip,Iq)

                rhs(1,I) = rhs(1,I) + wq*(dhdx*udp + dhdy*vdp)
                rhs(2,I) = rhs(2,I) + wq*(hi*sc_x + dhdx*Qu(1) + dhdy*Qu(2))
                rhs(3,I) = rhs(3,I) + wq*(hi*sc_y + dhdx*Qv(1) + dhdy*Qv(2))

            end do
        end do

    end subroutine create_rhs_btp_volume_qdf

    subroutine create_btp_fluxes_qdf(rhs,qb,qprime_df)

        use mod_basis, only: psiq, ngl
        use mod_grid, only:  npoin, intma, nface, face, face_type
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: H_face_ave,ope_face_ave,btp_mass_flux_face_ave, &
                                Qu_face_ave, Qv_face_ave, one_plus_eta_edge_2_ave, &
                                uvb_face_ave, ope2_face_ave
        use mod_initial, only: alpha_mlswe, pbprime_df_face, pbprime_df
        use mod_constants, only: gravity
        use mod_input, only: nlayers, dry_cutoff

        implicit none

        real, dimension(3, npoin), intent(inout)  :: rhs
        real, dimension(4, npoin), intent(in) :: qb
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        integer :: iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_u, flux_v, flux_pb, nxr, nyr, un
        real, dimension(nq) :: quu, quv, qvu, qvv, H_face_tmp
        real :: pbl, pbr, lambq, one_eta
        real :: pU_L, pU_R, lamb, dispu, dispv, pbpert_edge
        real :: qbl(4), qbr(4), lam1, lam2, dispp
        real :: c_minus, c_plus, c_pb_L, c_pb_R, c_pbub_LR, c_pbub_L, c_pbub_R
        real, dimension(nq) :: c_pb_LR
        real :: Q_uu_q, Q_uv_q, Q_vu_q, Q_vv_q
        real :: ppl, ppr, upl, upr, vpl, vpr, pb_cutoff
        real, dimension(nlayers+1) :: pprime_l, pprime_r
        integer :: itype, k
        real :: pkl, pkr, ukl, ukr, vkl, vkr, H_bcl_ql, H_bcl_qr, ul, ur, vl, vr, H_bcl_q, clam
        real, dimension(2) :: Qu_ql, Qu_qr, Qv_ql, Qv_qr
        real :: flux(3,nq), fxl, fxr, flux_edge_x, flux_edge_y

        do iface = 1, nface

            !Skip boundary faces
            if (face_type(iface) == 2) cycle

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1,nq

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                nxr = -nxl; nyr = -nyl

                qbl = 0.0; qbr = 0.0 ; pbl = 0.0; pbr = 0.0
                do n = 1, ngl
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)
                    I = intma(il,jl,kl,el)

                    hi = psiq(n,iquad)
                    qbl(:) = qbl(:) + hi*qb(:,I)
                    pbl = pbl + hi*pbprime_df(I)
                end do

                if (er > 0) then
                    do n = 1, ngl
                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)
                        I = intma(ir,jr,kr,er)

                        hi = psiq(n,iquad)
                        qbr(:) = qbr(:) + hi*qb(:,I)
                        pbr = pbr + hi*pbprime_df(I)
                    end do
                else
                    qbr(:) = qbl(:)
                    pbr = pbl

                    ! Apply boundary conditions
                    if(er == -4) then
                        un = nxl*qbl(3) + nyl*qbl(4)
                        qbr(3) = qbl(3) - 2.0*un*nxl
                        qbr(4) = qbl(4) - 2.0*un*nyl
                    elseif(er == -2) then
                        qbr(3:4) = -qbl(3:4)
                    end if
                endif 

                pU_L = nxl * qbl(3) + nyl * qbl(4)
                pU_R = nxr * qbr(3) + nyr * qbr(4)

                c_minus = sqrt(alpha_mlswe(nlayers) * pbr)
                c_plus  = sqrt(alpha_mlswe(nlayers) * pbl)
                clam = max(c_minus, c_plus)

                pbpert_edge = 0.5*(qbl(2) + qbr(2)) + (0.5/clam) * (pU_L + pU_R)
                one_eta = 1.0 + (pbpert_edge/pbl)

                ul = qbl(3)/qbl(1); ur = qbr(3)/qbr(1)
                vl = qbl(4)/qbl(1); vr = qbr(4)/qbr(1)

                Qu_ql(1) = ul*qbl(3) ; Qu_ql(2) = vl*qbl(3)
                Qu_qr(1) = ur*qbr(3) ; Qu_qr(2) = vr*qbr(3)

                Qv_ql(1) = ul*qbl(4) ; Qv_ql(2) = vl*qbl(4)
                Qv_qr(1) = ur*qbr(4) ; Qv_qr(2) = vr*qbr(4)

                pprime_l(:) = 0.0; pprime_r(:) = 0.0 ; H_bcl_q = 0.0
                H_bcl_ql = 0.0 ; H_bcl_qr = 0.0
                do k = 1, nlayers
                    ppl = 0.0; ppr = 0.0 ; upl = 0.0; upr = 0.0 ; vpl = 0.0; vpr = 0.0
                    do n = 1, ngl
                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)
                        I = intma(il,jl,kl,el)

                        hi = psiq(n,iquad)
                        ppl = ppl + hi*qprime_df(1,I,k)
                        upl = upl + hi*qprime_df(2,I,k)
                        vpl = vpl + hi*qprime_df(3,I,k)
                    enddo

                    if (er > 0) then
                        do n = 1, ngl
                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)

                            hi = psiq(n,iquad)
                            ppr = ppr + hi*qprime_df(1,I,k)
                            upr = upr + hi*qprime_df(2,I,k)
                            vpr = vpr + hi*qprime_df(3,I,k)
                        enddo
                    else
                        ppr = ppl
                        upr = upl
                        vpr = vpl
                    endif

                    ! Apply boundary conditions
                    if(er == -4) then
                        un = nxl*upl + nyl*vpl
                        upr = upl - 2.0*un*nxl
                            vpr = vpl - 2.0*un*nyl
                    elseif(er == -2) then
                        upr = -upl
                        vpr = -vpl
                    end if

                    Qu_ql(1) = Qu_ql(1) + (upl*(upl*(one_eta*ppl)))
                    Qu_ql(2) = Qu_ql(2) + (vpl*(upl*(one_eta*ppl)))
                    Qu_qr(1) = Qu_qr(1) + (upr*(upr*(one_eta*ppr)))
                    Qu_qr(2) = Qu_qr(2) + (vpr*(upr*(one_eta*ppr)))

                    Qv_ql(1) = Qv_ql(1) + (upl*(vpl*(one_eta*ppl)))
                    Qv_ql(2) = Qv_ql(2) + (vpl*(vpl*(one_eta*ppl)))
                    Qv_qr(1) = Qv_qr(1) + (upr*(vpr*(one_eta*ppr)))
                    Qv_qr(2) = Qv_qr(2) + (vpr*(vpr*(one_eta*ppr)))

                    pprime_l(k+1) = pprime_l(k) + ppl
                    pprime_r(k+1) = pprime_r(k) + ppr

                    H_bcl_ql  = H_bcl_ql + 0.5*alpha_mlswe(k)*(pprime_l(k+1)**2 - pprime_l(k)**2)
                    H_bcl_qr  = H_bcl_qr + 0.5*alpha_mlswe(k)*(pprime_r(k+1)**2 - pprime_r(k)**2)
                    
                end do

                H_bcl_q = 0.5*(H_bcl_ql + H_bcl_qr)

                ! Compute mass fluxes at each element face.

                flux_edge_x = 0.5*(qbl(3) + qbr(3)) + (0.5*clam)*(nxl * qbl(2) + nxr * qbr(2))
                flux_edge_y = 0.5*(qbl(4) + qbr(4)) + (0.5*clam)*(nyl * qbl(2) + nyr * qbr(2))
                ! flux(1,iquad) = nxl*flux_edge_x(iquad) + nyl*flux_edge_y(iquad)

                fxl = nxl*qbl(3) + nyl*qbl(4)
                fxr = nxr*qbr(3) + nyr*qbr(4)
                flux(1,iquad) = 0.5*(fxl - fxr) - 0.5*clam*(qbr(2) - qbl(2))

                ! Compute pressure forcing H_face at each element face.
                H_bcl_q = (one_eta**2) * H_bcl_q

                btp_mass_flux_face_ave(1,iquad,iface) = btp_mass_flux_face_ave(1,iquad,iface) + flux_edge_x
                btp_mass_flux_face_ave(2,iquad,iface) = btp_mass_flux_face_ave(2,iquad,iface) + flux_edge_y
                ! btp_mass_flux_face_ave(1,iquad,iface) = btp_mass_flux_face_ave(1,iquad,iface) + flux(1,iquad)
                ! btp_mass_flux_face_ave(2,iquad,iface) = btp_mass_flux_face_ave(2,iquad,iface) + flux(1,iquad)

                H_face_ave(iquad,iface) = H_face_ave(iquad,iface) + H_bcl_q
                Qu_face_ave(1,iquad,iface) = Qu_face_ave(1,iquad,iface) + 0.5*(Qu_ql(1) + Qu_qr(1))
                Qu_face_ave(2,iquad,iface) = Qu_face_ave(2,iquad,iface) + 0.5*(Qu_ql(2) + Qu_qr(2))
                Qv_face_ave(1,iquad,iface) = Qv_face_ave(1,iquad,iface) + 0.5*(Qv_ql(1) + Qv_qr(1))
                Qv_face_ave(2,iquad,iface) = Qv_face_ave(2,iquad,iface) + 0.5*(Qv_ql(2) + Qv_qr(2))
                ope_face_ave(1,iquad,iface) = ope_face_ave(1,iquad,iface) + (1.0 + (qbl(2)/pbl))
                ope_face_ave(2,iquad,iface) = ope_face_ave(2,iquad,iface) + (1.0 + (qbr(2)/pbr))
                ope2_face_ave(1,iquad,iface) = ope2_face_ave(1,iquad,iface) + (1.0 + (qbl(2)/pbl))**2
                ope2_face_ave(2,iquad,iface) = ope2_face_ave(2,iquad,iface) + (1.0 + (qbr(2)/pbr))**2
                one_plus_eta_edge_2_ave(iquad,iface) = one_plus_eta_edge_2_ave(iquad,iface) &
                                                    + one_eta**2
                uvb_face_ave(1,1,iquad,iface) = uvb_face_ave(1,1,iquad,iface) + ul
                uvb_face_ave(1,2,iquad,iface) = uvb_face_ave(1,2,iquad,iface) + ur
                uvb_face_ave(2,1,iquad,iface) = uvb_face_ave(2,1,iquad,iface) + vl
                uvb_face_ave(2,2,iquad,iface) = uvb_face_ave(2,2,iquad,iface) + vr

                Qu_ql(1) = Qu_ql(1) + (one_eta**2)*H_bcl_ql
                Qu_qr(1) = Qu_qr(1) + (one_eta**2)*H_bcl_qr

                fxl = nxl*Qu_ql(1) + nyl*Qu_ql(2)
                fxr = nxr*Qu_qr(1) + nyr*Qu_qr(2)

                flux(2,iquad) = 0.5*(fxl - fxr) - (0.25*clam)*(qbr(3) - qbl(3))

                Qv_ql(2) = Qv_ql(2) + (one_eta**2)*H_bcl_ql
                Qv_qr(2) = Qv_qr(2) + (one_eta**2)*H_bcl_qr

                fxl = nxl*Qv_ql(1) + nyl*Qv_ql(2)
                fxr = nxr*Qv_qr(1) + nyr*Qv_qr(2)

                flux(3,iquad) = 0.5*(fxl - fxr) - (0.25*clam)*(qbr(4) - qbl(4))

            end do
          

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                flux_pb = flux(1,iquad)
                flux_u = flux(2,iquad)
                flux_v = flux(3,iquad)

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)
                    I = intma(il,jl,kl,el)

                    rhs(1,I) = rhs(1,I) - wq*hi*flux_pb
                    rhs(2,I) = rhs(2,I) - wq*hi*flux_u
                    rhs(3,I) = rhs(3,I) - wq*hi*flux_v

                    if(er > 0) then
                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)
                        I = intma(ir,jr,kr,er)

                        rhs(1,I) = rhs(1,I) + wq*hi*flux_pb
                        rhs(2,I) = rhs(2,I) + wq*hi*flux_u
                        rhs(3,I) = rhs(3,I) + wq*hi*flux_v

                    end if
                end do
            end do
        end do

    end subroutine create_btp_fluxes_qdf
    
end module mod_rhs_btp
