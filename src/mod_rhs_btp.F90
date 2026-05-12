module mod_rhs_btp

    ! Module-level imports: only what create_rhs_btp uses directly.
    ! Each contained subroutine declares its own dependencies via local 'use'.
    use mod_grid,           only: npoin
    use mod_input,          only: nlayers, method_visc
    use mod_metrics,        only: massinv
    use mod_laplacian_quad, only: btp_create_laplacian

    implicit none
    public :: create_rhs_btp

contains

    !===========================================================================
    subroutine create_rhs_btp(rhs_btp, qb_df, qprime_df)
    !===========================================================================
    !  Top-level barotropic RHS driver.
    !
    !  Execution order (GPU data management notes):
    !    1.  MPI pre-comm  — host; qb_df/qprime_df must be on the host here.
    !    2.  GPU upload    — enter data copyin(qb_df,qprime_df) create(rhs)
    !    3.  GPU kernels   — volume integral then face flux
    !    4.  GPU download  — exit data copyout(rhs) delete(qb_df,qprime_df)
    !    5.  MPI post-comm — host; adds neighbour contributions into host rhs
    !    6.  Viscosity + mass-matrix inverse — host
    !
    !  Subroutine arguments (rhs, qb_df, qprime_df) are managed here.
    !  All static module arrays are entered once via openacc_enter_data.
    !===========================================================================
        implicit none

        real, dimension(3,npoin),          intent(inout) :: rhs_btp
        real, dimension(4,npoin),          intent(in)  :: qb_df
        real, dimension(3,npoin,nlayers),  intent(in)  :: qprime_df

        real, dimension(2,npoin) :: rhs_visc_btp

        rhs_visc_btp = 0.0

        ! 1. MPI halo exchange (host).
        !    qb_df / qprime_df must reside on the host at this point.
        !    With standard (non-CUDA-aware) MPI this must precede the upload.
    
        call btp_create_precommunicator(qb_df, qprime_df, 4)

        ! 2. Upload per-call arrays to the device.
        !    rhs  : 'create' (not 'copyin') — volume kernel zeroes it on device.
        !    qb_df, qprime_df : read-only on GPU, 'copyin' only.
        
        !!$acc enter data copyin(qb_df, qprime_df) create(rhs)

        ! 3. GPU kernels.
        !    create_rhs_btp_volume_qdf zeroes rhs on the device first,
        !    then both routines accumulate into it.
    
        call create_rhs_btp_volume_qdf(rhs_btp, qb_df, qprime_df)

        !!$acc update host(rhs_btp)

        call create_btp_fluxes_qdf(rhs_btp, qb_df, qprime_df)

        ! 4. Download rhs to host; free device temporaries.
        !    qb_df / qprime_df were read-only: no copyout needed.
        !!$acc exit data copyout(rhs) delete(qb_df, qprime_df)

        ! 5. MPI post-receive: accumulate neighbour contributions into host rhs.
        call btp_create_postcommunicator(rhs_btp, 4)

        ! 6. Viscosity (host; not yet GPU-ported).
        
        if (method_visc > 0) call btp_create_laplacian(rhs_visc_btp, qb_df)

        ! 7. Mass-matrix inverse scaling (host, after GPU download).

        rhs_btp(1,:) = massinv(:) * rhs_btp(1,:)
        rhs_btp(2,:) = massinv(:) * (rhs_btp(2,:) + rhs_visc_btp(1,:))
        rhs_btp(3,:) = massinv(:) * (rhs_btp(3,:) + rhs_visc_btp(2,:))

    end subroutine create_rhs_btp


    !===========================================================================
    subroutine create_rhs_btp_volume_qdf(rhs_btp, qb_df, qprime_df)
    !===========================================================================
    !  Volume contribution to the barotropic RHS (DG weak form, interior).
    !
    !  GPU parallelism: one gang per quadrature point Iq.  All inner loops
    !  (ip over basis nodes, k over layers) are sequential within the gang.
    !
    !  Data assumptions:
    !    - rhs, qb_df, qprime_df : entered by the caller via '!$acc enter data'
    !    - all module arrays      : entered once via openacc_enter_data
    !    - scalar module variables (botfr, cd_mlswe, gravity, nlayers) :
    !      implicitly firstprivate in the parallel region; no present clause needed.
    !
    !  Race conditions: Iq is unique per node (each node has exactly one owning
    !  Iq), so different gangs scatter to disjoint rhs indices. No atomics needed.
    !===========================================================================
        use mod_grid,      only: npoin_q, npoin
        use mod_basis,     only: npts
        use mod_constants, only: gravity
        use mod_initial,   only: grad_zbot_quad, tau_wind, psih, dpsidx, dpsidy, &
                                 indexq, wjac, pbprime_df, coriolis_quad, alpha_mlswe
        use mod_variables, only: tau_bot_ave, H_ave, Qu_ave, Quv_ave, Qv_ave,    &
                                 ope_ave, uvb_ave, btp_mass_flux_ave, ope2_ave
        use mod_input,     only: nlayers, cd_mlswe, botfr

        implicit none

        real, dimension(4,npoin),         intent(in)  :: qb_df
        real, dimension(3,npoin,nlayers), intent(in)  :: qprime_df
        real, dimension(3,npoin),         intent(inout) :: rhs_btp

        real :: sc_x, sc_y, Hq, Qu1, Qu2, Qv1, Qv2
        real :: wq, hi, dhdx, dhdy, tb_u, tb_v, ope, ope2, grav_dp
        real :: dp, dpp, udp, vdp, ub, vb, ubot, vbot, spd, pbq
        real :: pp_k, up_k, vp_k, pprime_k, pprime_k1
        real :: sum_up2, sum_uv, sum_vp2
        real :: wq_udp, wq_vdp, wq_scx, wq_scy, wq_Qu1, wq_Qu2, wq_Qv1, wq_Qv2
        ! real,    dimension(npts) :: hi_c
        ! integer, dimension(npts) :: I_c
        integer :: Iq, ip, k, I

        !!$acc data present(rhs_btp, qb_df, qprime_df,                                  &
        !!$acc               grad_zbot_quad, tau_wind, psih, dpsidx, dpsidy,        &
        !!$acc               indexq, wjac, pbprime_df, coriolis_quad, alpha_mlswe,  &
        !!$acc               tau_bot_ave, H_ave, Qu_ave, Quv_ave, Qv_ave, ope_ave,  &
        !!$acc               uvb_ave, btp_mass_flux_ave, ope2_ave)

        ! Zero rhs on device. Caller used 'create(rhs)' (not 'copyin'), so
        ! device memory is uninitialised; both kernels accumulate into rhs.
        !!$acc kernels
        rhs_btp = 0.0
        !!$acc end kernels

        !!$acc parallel loop gang                                                    &
        !!$acc   private(dp, dpp, udp, vdp, pbq, hi, Hq, wq, ub, vb,               &
        !!$acc           ubot, vbot, spd, tb_u, tb_v, sc_x, sc_y,                   &
        !!$acc           ope, ope2, grav_dp,                                         &
        !!$acc           pp_k, up_k, vp_k, pprime_k, pprime_k1,                     &
        !!$acc           Qu1, Qu2, Qv1, Qv2,                                        &
        !!$acc           dhdx, dhdy, sum_up2, sum_uv, sum_vp2,                      &
        !!$acc           wq_udp, wq_vdp, wq_scx, wq_scy,                            &
        !!$acc           wq_Qu1, wq_Qu2, wq_Qv1, wq_Qv2,                           &
        !!$acc           ip, k)
        do Iq = 1, npoin_q

            wq = wjac(Iq)

            ! Cache basis values for this quadrature point.
            !!$acc loop seq
            ! do ip = 1, npts
            !     I_c(ip)  = indexq(ip, Iq)
            !     hi_c(ip) = psih(ip, Iq)
            ! end do

            ! Barotropic projection
            dp = 0.0;  dpp = 0.0;  udp = 0.0;  vdp = 0.0;  pbq = 0.0
            !!$acc loop seq
            do ip = 1, npts
                I = indexq(ip, Iq)
                hi = psih(ip, Iq)
                dp  = dp  + hi * qb_df(1, I)
                dpp = dpp + hi * qb_df(2, I)
                udp = udp + hi * qb_df(3, I)
                vdp = vdp + hi * qb_df(4, I)
                pbq = pbq + hi * pbprime_df(I)
            end do

            ub = udp / dp
            vb = vdp / dp

            ! Layer projections, Hq, momentum-flux sums 
            Hq = 0.0;  sum_up2 = 0.0;  sum_uv = 0.0;  sum_vp2 = 0.0
            pprime_k = 0.0     ! pprime_0 = 0 at sea surface

            !!$acc loop seq
            do k = 1, nlayers
                pp_k = 0.0;  up_k = 0.0;  vp_k = 0.0
                !!$acc loop seq
                do ip = 1, npts
                    I = indexq(ip, Iq)
                    hi = psih(ip, Iq)
                    pp_k = pp_k + hi * qprime_df(1, I, k)
                    up_k = up_k + hi * qprime_df(2, I, k)
                    vp_k = vp_k + hi * qprime_df(3, I, k)
                end do

                pprime_k1 = pprime_k + pp_k

                ! Diff-of-squares: pprime_{k+1}^2 - pprime_k^2
                !   = (pprime_{k+1} + pprime_k) * pp_k
                Hq      = Hq      + 0.5 * alpha_mlswe(k) * (pprime_k1 + pprime_k) * pp_k
                sum_up2 = sum_up2 + up_k * up_k * pp_k
                sum_uv  = sum_uv  + up_k * vp_k * pp_k
                sum_vp2 = sum_vp2 + vp_k * vp_k * pp_k

                pprime_k = pprime_k1
            end do

            ! Bottom friction
            tb_u = 0.0;  tb_v = 0.0
            if (botfr /= 0) then
                ubot = up_k + ub
                vbot = vp_k + vb
                if (botfr == 1) then
                    spd = (cd_mlswe / gravity) * pp_k
                else
                    spd = (cd_mlswe / alpha_mlswe(nlayers)) * sqrt(ubot**2 + vbot**2)
                end if
                tb_u = spd * ubot
                tb_v = spd * vbot
            end if

            ! Free-surface factor 
            ope  = 1.0 + (dpp / pbq)
            ope2 = ope * ope
            Hq   = ope2 * Hq

            ! Source terms
            grav_dp = gravity * dp
            sc_x =  coriolis_quad(Iq) * vdp                                  &
                    + gravity * (tau_wind(1,Iq) - tb_u)                       &
                    - grav_dp * grad_zbot_quad(1,Iq)
            sc_y = -coriolis_quad(Iq) * udp                                  &
                    + gravity * (tau_wind(2,Iq) - tb_v)                       &
                    - grav_dp * grad_zbot_quad(2,Iq)

            ! Flux tensors
            Qu1 = ub * udp + ope * sum_up2
            Qu2 = vb * udp + ope * sum_uv
            Qv1 = ub * vdp + ope * sum_uv
            Qv2 = vb * vdp + ope * sum_vp2

            ! Time-average accumulators
            H_ave(Iq)               = H_ave(Iq)               + Hq
            Qu_ave(Iq)              = Qu_ave(Iq)              + Qu1
            Qv_ave(Iq)              = Qv_ave(Iq)              + Qv1
            Quv_ave(Iq)             = Quv_ave(Iq)             + Qu2
            tau_bot_ave(1,Iq)       = tau_bot_ave(1,Iq)       + tb_u
            tau_bot_ave(2,Iq)       = tau_bot_ave(2,Iq)       + tb_v
            ope_ave(Iq)             = ope_ave(Iq)             + ope
            ope2_ave(Iq)            = ope2_ave(Iq)            + ope2
            btp_mass_flux_ave(1,Iq) = btp_mass_flux_ave(1,Iq) + udp
            btp_mass_flux_ave(2,Iq) = btp_mass_flux_ave(2,Iq) + vdp
            uvb_ave(1,Iq)           = uvb_ave(1,Iq)           + ub
            uvb_ave(2,Iq)           = uvb_ave(2,Iq)           + vb

            ! Add hydrostatic pressure to the diagonal flux components.
            Qu1 = Qu1 + Hq
            Qv2 = Qv2 + Hq

            ! ── Weight by quadrature Jacobian ─────────────────────────────────
            wq_udp = wq * udp;   wq_vdp = wq * vdp
            wq_scx = wq * sc_x;  wq_scy = wq * sc_y
            wq_Qu1 = wq * Qu1;   wq_Qu2 = wq * Qu2
            wq_Qv1 = wq * Qv1;   wq_Qv2 = wq * Qv2

            ! Volume RHS
            ! Iq is unique per node
            !!$acc loop seq
            do ip = 1, npts
                I = indexq(ip, Iq)
                hi = psih(ip, Iq)
                dhdx = dpsidx(ip, Iq)
                dhdy = dpsidy(ip, Iq)
                !!$acc atomic update
                rhs_btp(1, I) = rhs_btp(1, I) + dhdx*wq_udp + dhdy*wq_vdp
                !!$acc atomic update
                rhs_btp(2, I) = rhs_btp(2, I) + hi*wq_scx + dhdx*wq_Qu1 + dhdy*wq_Qu2
                !!$acc atomic update
                rhs_btp(3, I) = rhs_btp(3, I) + hi*wq_scy + dhdx*wq_Qv1 + dhdy*wq_Qv2
            end do

        end do
        !!$acc end data

    end subroutine create_rhs_btp_volume_qdf

    !===========================================================================
    subroutine create_btp_fluxes_qdf(rhs_btp, qb, qprime_df)
    !===========================================================================
    !  Face flux contribution to the barotropic RHS (Riemann solver).
    !
    !  GPU parallelism: one gang per face.  The quadrature loop (iquad) and all
    !  inner loops (basis nodes n, layers k) are sequential within each gang.
    !  Running iquad sequentially avoids intra-gang atomics on rhs when different
    !  quadrature points of the same face share boundary nodes.
    !
    !  Key structural change vs. original:
    !    The six runtime-sized gang-private arrays ppl/upl/vpl/ppr/upr/vpr
    !    (dimension nlayers) have been eliminated.  Layer projection is now fused
    !    into the k-loop: each iteration projects from global memory into scalars
    !    pkl/ukl/vkl and pkr/ukr/vkr, immediately accumulates into flux totals,
    !    then discards.  The compiler can keep these in registers.
    !
    !    Scalar variables Qu_ql1/2, Qu_qr1/2, Qv_ql1/2, Qv_qr1/2 replace the
    !    original Qu_ql(2)/Qv_qr(2) fixed-size arrays for the same reason.
    !
    !  Data assumptions:
    !    - rhs, qb, qprime_df : entered by the caller via '!$acc enter data'
    !    - all module arrays  : entered once via openacc_enter_data
    !
    !  Race conditions: rhs(m,I) is shared across faces (different gangs can
    !  reach the same node I).  All rhs scatter writes use '!$acc atomic update'.
    !  Time-average arrays are indexed by (iquad,iface) — each (gang,iquad) pair
    !  is unique; no atomics required for those.
    !===========================================================================
        use mod_basis,     only: psiq, ngl, nq
        use mod_grid,      only: npoin, intma, nface, face, face_type
        use mod_face,      only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: H_face_ave, ope_face_ave, btp_mass_flux_face_ave, &
                                 Qu_face_ave, Qv_face_ave, one_plus_eta_edge_2_ave, &
                                 uvb_face_ave, ope2_face_ave
        use mod_initial,   only: alpha_mlswe, pbprime_df
        use mod_constants, only: gravity
        use mod_input,     only: nlayers

        implicit none

        real, dimension(3,npoin),          intent(inout) :: rhs_btp
        real, dimension(4,npoin),          intent(in)    :: qb
        real, dimension(3,npoin,nlayers),  intent(in)    :: qprime_df

        integer :: iface, iquad, el, er, il, jl, kl, ir, jr, kr, I, n, k
        real    :: wq, hi, nxl, nyl, un
        real    :: pbl, pbr, clam, one_eta, one_eta2, half_clam, quarter_clam
        real    :: pU_L, pU_R, half_clam_dqb2, pbpert_edge
        real    :: qbl(4), qbr(4)
        real    :: c_minus, c_plus
        real    :: ul, ur, vl, vr, opl, opr
        real    :: H_bcl_ql, H_bcl_qr, H_bcl_q, oe2_Hql, oe2_Hqr
        real    :: Qu_ql1, Qu_ql2, Qu_qr1, Qu_qr2
        real    :: Qv_ql1, Qv_ql2, Qv_qr1, Qv_qr2
        real    :: flux_edge_x, flux_edge_y, flux_pb, flux_u, flux_v
        real    :: pprime_lk, pprime_rk, pprime_lk1, pprime_rk1
        real    :: pkl, pkr, ukl, ukr, vkl, vkr
        real    :: ope_ppl_k, ope_ppr_k, uv_cross_l, uv_cross_r
        real    :: wq_flux_pb, wq_flux_u, wq_flux_v
        real,    dimension(ngl) :: hi_c
        integer, dimension(ngl) :: I_l, I_r

        !!$acc data present(rhs_btp, qb, qprime_df,                                     &
        !!$acc               face, face_type, intma, imapl, imapr,                  &
        !!$acc               normal_vector_q, jac_faceq, psiq,                      &
        !!$acc               pbprime_df, alpha_mlswe,                               &
        !!$acc               H_face_ave, ope_face_ave, btp_mass_flux_face_ave,      &
        !!$acc               Qu_face_ave, Qv_face_ave, one_plus_eta_edge_2_ave,     &
        !!$acc               uvb_face_ave, ope2_face_ave)

        ! One gang per face. nface (O(millions)) saturates the GPU without
        ! a vector loop over iquad.  Running iquad sequentially within the gang
        ! avoids: (a) vector-private copies of hi_c/I_l/I_r arrays,
        !         (b) intra-gang atomics on rhs from different iquad points
        !             sharing a boundary node.
        !!$acc parallel loop gang                                                    &
        !!$acc   private(el, er, il, jl, kl, ir, jr, kr, I, n, k,                  &
        !!$acc           I_l, I_r, hi_c,                                            &
        !!$acc           qbl, qbr, pbl, pbr,                                        &
        !!$acc           nxl, nyl, wq, hi, un,                                      &
        !!$acc           clam, one_eta, one_eta2, half_clam, quarter_clam,          &
        !!$acc           pU_L, pU_R, half_clam_dqb2, pbpert_edge,                  &
        !!$acc           c_minus, c_plus, ul, ur, vl, vr, opl, opr,                &
        !!$acc           H_bcl_ql, H_bcl_qr, H_bcl_q, oe2_Hql, oe2_Hqr,           &
        !!$acc           Qu_ql1, Qu_ql2, Qu_qr1, Qu_qr2,                           &
        !!$acc           Qv_ql1, Qv_ql2, Qv_qr1, Qv_qr2,                           &
        !!$acc           flux_edge_x, flux_edge_y, flux_pb, flux_u, flux_v,        &
        !!$acc           pprime_lk, pprime_rk, pprime_lk1, pprime_rk1,             &
        !!$acc           pkl, pkr, ukl, ukr, vkl, vkr,                             &
        !!$acc           ope_ppl_k, ope_ppr_k, uv_cross_l, uv_cross_r,             &
        !!$acc           wq_flux_pb, wq_flux_u, wq_flux_v)
        do iface = 1, nface

            if (face_type(iface) == 2) cycle

            el = face(7, iface)
            er = face(8, iface)

            ! Precompute node indices once per face
            ! I_l / I_r are reused across all iquad iterations within this gang.
            !!$acc loop seq
            do n = 1, ngl
                il = imapl(1,n,1,iface)
                jl = imapl(2,n,1,iface)
                kl = imapl(3,n,1,iface)
                I_l(n) = intma(il, jl, kl, el)
            end do

            if (er > 0) then
                !!$acc loop seq
                do n = 1, ngl
                    ir = imapr(1,n,1,iface)
                    jr = imapr(2,n,1,iface)
                    kr = imapr(3,n,1,iface)
                    I_r(n) = intma(ir, jr, kr, er)
                end do
            end if

            ! Quadrature loop (sequential within gang)
            !!$acc loop seq
            do iquad = 1, nq

                nxl = normal_vector_q(1, iquad, 1, iface)
                nyl = normal_vector_q(2, iquad, 1, iface)

                !!$acc loop seq
                do n = 1, ngl
                    hi_c(n) = psiq(n, iquad)
                end do

                ! Project barotropic LEFT state
                qbl(1) = 0.0;  qbl(2) = 0.0;  qbl(3) = 0.0;  qbl(4) = 0.0
                pbl = 0.0
                !!$acc loop seq
                do n = 1, ngl
                    hi = hi_c(n)
                    I  = I_l(n)
                    qbl(1) = qbl(1) + hi * qb(1,I)
                    qbl(2) = qbl(2) + hi * qb(2,I)
                    qbl(3) = qbl(3) + hi * qb(3,I)
                    qbl(4) = qbl(4) + hi * qb(4,I)
                    pbl    = pbl    + hi * pbprime_df(I)
                end do

                ! Project barotropic RIGHT state or apply barotropic BC
                if (er > 0) then
                    qbr(1) = 0.0;  qbr(2) = 0.0;  qbr(3) = 0.0;  qbr(4) = 0.0
                    pbr = 0.0
                    !!$acc loop seq
                    do n = 1, ngl
                        hi = hi_c(n)
                        I  = I_r(n)
                        qbr(1) = qbr(1) + hi * qb(1,I)
                        qbr(2) = qbr(2) + hi * qb(2,I)
                        qbr(3) = qbr(3) + hi * qb(3,I)
                        qbr(4) = qbr(4) + hi * qb(4,I)
                        pbr    = pbr    + hi * pbprime_df(I)
                    end do
                else
                    ! Boundary: mirror left state, then apply velocity BC.
                    qbr(1) = qbl(1);  qbr(2) = qbl(2)
                    qbr(3) = qbl(3);  qbr(4) = qbl(4)
                    pbr = pbl
                    if (er == -4) then          ! slip wall
                        un     = nxl*qbl(3) + nyl*qbl(4)
                        qbr(3) = qbl(3) - 2.0*un*nxl
                        qbr(4) = qbl(4) - 2.0*un*nyl
                    else if (er == -2) then     ! no-slip wall
                        qbr(3) = -qbl(3);  qbr(4) = -qbl(4)
                    end if
                end if

                ! Wave speeds
                pU_L = nxl*qbl(3) + nyl*qbl(4)
                pU_R = -(nxl*qbr(3) + nyl*qbr(4))

                c_minus      = sqrt(alpha_mlswe(nlayers) * pbr)
                c_plus       = sqrt(alpha_mlswe(nlayers) * pbl)
                clam         = max(c_minus, c_plus)
                half_clam    = 0.5  * clam
                quarter_clam = 0.25 * clam

                pbpert_edge = 0.5*(qbl(2) + qbr(2)) + (0.5/clam)*(pU_L + pU_R)
                one_eta     = 1.0 + pbpert_edge / pbl
                one_eta2    = one_eta * one_eta

                ul = qbl(3) / qbl(1);  ur = qbr(3) / qbr(1)
                vl = qbl(4) / qbl(1);  vr = qbr(4) / qbr(1)

                ! Initialise flux tensors from barotropic state
                Qu_ql1 = ul * qbl(3);  Qu_ql2 = vl * qbl(3)
                Qu_qr1 = ur * qbr(3);  Qu_qr2 = vr * qbr(3)
                Qv_ql1 = ul * qbl(4);  Qv_ql2 = vl * qbl(4)
                Qv_qr1 = ur * qbr(4);  Qv_qr2 = vr * qbr(4)
                H_bcl_ql  = 0.0;  H_bcl_qr  = 0.0
                pprime_lk = 0.0;  pprime_rk = 0.0

                ! Combined layer projection + flux accumulation
                ! Eliminates the six ppl/upl/vpl/ppr/upr/vpr(nlayers) gang-private
                ! arrays: each k iteration projects from global memory into three
                ! scalar pairs, accumulates into running flux totals, then discards.
                ! The compiler can keep pkl/ukl/vkl/pkr/ukr/vkr in registers.
                ! one_eta is available here because it depends only on the
                ! barotropic state, which was projected before this loop.
                !!$acc loop seq
                do k = 1, nlayers

                    ! Project left layer k
                    pkl = 0.0;  ukl = 0.0;  vkl = 0.0
                    !!$acc loop seq
                    do n = 1, ngl
                        hi  = hi_c(n)
                        I   = I_l(n)
                        pkl = pkl + hi * qprime_df(1, I, k)
                        ukl = ukl + hi * qprime_df(2, I, k)
                        vkl = vkl + hi * qprime_df(3, I, k)
                    end do

                    ! Project right layer k, or apply velocity BC for this layer.
                    if (er > 0) then
                        pkr = 0.0;  ukr = 0.0;  vkr = 0.0
                        !!$acc loop seq
                        do n = 1, ngl
                            hi  = hi_c(n)
                            I   = I_r(n)
                            pkr = pkr + hi * qprime_df(1, I, k)
                            ukr = ukr + hi * qprime_df(2, I, k)
                            vkr = vkr + hi * qprime_df(3, I, k)
                        end do
                    else
                        pkr = pkl
                        if (er == -4) then
                            un  = nxl * ukl + nyl * vkl
                            ukr = ukl - 2.0 * un * nxl
                            vkr = vkl - 2.0 * un * nyl
                        else if (er == -2) then
                            ukr = -ukl;  vkr = -vkl
                        else
                            ukr = ukl;  vkr = vkl
                        end if
                    end if

                    ! Flux tensor accumulation for layer k.
                    ope_ppl_k  = one_eta * pkl
                    ope_ppr_k  = one_eta * pkr
                    uv_cross_l = ukl * vkl * ope_ppl_k
                    uv_cross_r = ukr * vkr * ope_ppr_k

                    Qu_ql1 = Qu_ql1 + ukl * ukl * ope_ppl_k
                    Qu_ql2 = Qu_ql2 + uv_cross_l
                    Qv_ql1 = Qv_ql1 + uv_cross_l
                    Qv_ql2 = Qv_ql2 + vkl * vkl * ope_ppl_k

                    Qu_qr1 = Qu_qr1 + ukr * ukr * ope_ppr_k
                    Qu_qr2 = Qu_qr2 + uv_cross_r
                    Qv_qr1 = Qv_qr1 + uv_cross_r
                    Qv_qr2 = Qv_qr2 + vkr * vkr * ope_ppr_k

                    pprime_lk1 = pprime_lk + pkl
                    pprime_rk1 = pprime_rk + pkr
                    H_bcl_ql   = H_bcl_ql  + 0.5*alpha_mlswe(k) * (pprime_lk1 + pprime_lk) * pkl
                    H_bcl_qr   = H_bcl_qr  + 0.5*alpha_mlswe(k) * (pprime_rk1 + pprime_rk) * pkr
                    pprime_lk  = pprime_lk1
                    pprime_rk  = pprime_rk1

                end do  ! k

                H_bcl_q = 0.5 * (H_bcl_ql + H_bcl_qr)

                ! Mass flux
                half_clam_dqb2 = half_clam * (qbl(2) - qbr(2))
                flux_edge_x = 0.5*(qbl(3) + qbr(3)) + nxl * half_clam_dqb2
                flux_edge_y = 0.5*(qbl(4) + qbr(4)) + nyl * half_clam_dqb2
                flux_pb     = nxl*flux_edge_x + nyl*flux_edge_y

                H_bcl_q = one_eta2 * H_bcl_q

                ! Time-average accumulators
                btp_mass_flux_face_ave(1,iquad,iface) = btp_mass_flux_face_ave(1,iquad,iface) + flux_edge_x
                btp_mass_flux_face_ave(2,iquad,iface) = btp_mass_flux_face_ave(2,iquad,iface) + flux_edge_y

                H_face_ave(iquad,iface)    = H_face_ave(iquad,iface)    + H_bcl_q
                Qu_face_ave(1,iquad,iface) = Qu_face_ave(1,iquad,iface) + 0.5*(Qu_ql1 + Qu_qr1)
                Qu_face_ave(2,iquad,iface) = Qu_face_ave(2,iquad,iface) + 0.5*(Qu_ql2 + Qu_qr2)
                Qv_face_ave(1,iquad,iface) = Qv_face_ave(1,iquad,iface) + 0.5*(Qv_ql1 + Qv_qr1)
                Qv_face_ave(2,iquad,iface) = Qv_face_ave(2,iquad,iface) + 0.5*(Qv_ql2 + Qv_qr2)

                opl = 1.0 + qbl(2)/pbl
                opr = 1.0 + qbr(2)/pbr
                ope_face_ave(1,iquad,iface)  = ope_face_ave(1,iquad,iface)  + opl
                ope_face_ave(2,iquad,iface)  = ope_face_ave(2,iquad,iface)  + opr
                ope2_face_ave(1,iquad,iface) = ope2_face_ave(1,iquad,iface) + opl*opl
                ope2_face_ave(2,iquad,iface) = ope2_face_ave(2,iquad,iface) + opr*opr

                one_plus_eta_edge_2_ave(iquad,iface) = one_plus_eta_edge_2_ave(iquad,iface) + one_eta2

                uvb_face_ave(1,1,iquad,iface) = uvb_face_ave(1,1,iquad,iface) + ul
                uvb_face_ave(1,2,iquad,iface) = uvb_face_ave(1,2,iquad,iface) + ur
                uvb_face_ave(2,1,iquad,iface) = uvb_face_ave(2,1,iquad,iface) + vl
                uvb_face_ave(2,2,iquad,iface) = uvb_face_ave(2,2,iquad,iface) + vr

                ! Momentum fluxes
                oe2_Hql = one_eta2 * H_bcl_ql
                oe2_Hqr = one_eta2 * H_bcl_qr

                Qu_ql1 = Qu_ql1 + oe2_Hql
                Qu_qr1 = Qu_qr1 + oe2_Hqr
                flux_u = 0.5 * (nxl*(Qu_ql1 + Qu_qr1) + nyl*(Qu_ql2 + Qu_qr2)) &
                        - quarter_clam * (qbr(3) - qbl(3))

                Qv_ql2 = Qv_ql2 + oe2_Hql
                Qv_qr2 = Qv_qr2 + oe2_Hqr
                flux_v = 0.5 * (nxl*(Qv_ql1 + Qv_qr1) + nyl*(Qv_ql2 + Qv_qr2)) &
                        - quarter_clam * (qbr(4) - qbl(4))

                wq         = jac_faceq(iquad, 1, iface)
                wq_flux_pb = wq * flux_pb
                wq_flux_u  = wq * flux_u
                wq_flux_v  = wq * flux_v

                ! RHS scatter
                ! rhs(m,I) is shared across faces: different gangs (faces) may
                ! write to the same node I simultaneously.
                ! !$acc atomic update serialises conflicting writes.
                !!$acc loop seq
                do n = 1, ngl
                    hi = hi_c(n)
                    I  = I_l(n)
                    !!$acc atomic update
                    rhs_btp(1,I) = rhs_btp(1,I) - hi * wq_flux_pb
                    !!$acc atomic update
                    rhs_btp(2,I) = rhs_btp(2,I) - hi * wq_flux_u
                    !!$acc atomic update
                    rhs_btp(3,I) = rhs_btp(3,I) - hi * wq_flux_v
                end do

                if (er > 0) then
                    !!$acc loop seq
                    do n = 1, ngl
                        hi = hi_c(n)
                        I  = I_r(n)
                        !!$acc atomic update
                        rhs_btp(1,I) = rhs_btp(1,I) + hi * wq_flux_pb
                        !!$acc atomic update
                        rhs_btp(2,I) = rhs_btp(2,I) + hi * wq_flux_u
                        !!$acc atomic update
                        rhs_btp(3,I) = rhs_btp(3,I) + hi * wq_flux_v
                    end do
                end if

            end do  ! iquad

        end do  ! iface
        !!$acc end data

    end subroutine create_btp_fluxes_qdf

end module mod_rhs_btp