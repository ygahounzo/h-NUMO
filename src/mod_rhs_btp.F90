module mod_rhs_btp

   implicit none
   public :: create_rhs_btp

contains

   !===========================================================================
   subroutine create_rhs_btp(G, inp, b, mf, par, btp, init, ref, mpic, mt, tsp, &
      rhs_btp, qb_df, qprime_df)
      !===========================================================================
      !  Top-level barotropic RHS driver.
      !===========================================================================
      use mod_grid,             only: grid
      use mod_input,            only: input
      use mod_basis,            only: basis
      use mod_face,             only: face_CS
      use mod_parallel,         only: parallel_CS
      use mod_variables,        only: btp_CS
      use mod_initial,          only: initial
      use mod_ref,              only: mref
      use mod_mpi_communicator, only: mpi_communicator
      use mod_metrics,          only: metrics
      use mod_tensor,           only: tensor_CS
      use mod_laplacian_quad,   only: btp_create_laplacian

      implicit none

      type(grid),             intent(in)    :: G
      type(input),            intent(in)    :: inp
      type(basis),            intent(in)    :: b
      type(face_CS),          intent(in)    :: mf
      type(parallel_CS),      intent(in)    :: par
      type(btp_CS),           intent(inout) :: btp
      type(initial),          intent(in)    :: init
      type(mref),             intent(inout) :: ref
      type(mpi_communicator), intent(inout) :: mpic
      type(metrics),          intent(in)    :: mt
      type(tensor_CS),        intent(in)    :: tsp

      real, dimension(3,G%npoin),             intent(out) :: rhs_btp
      real, dimension(4,G%npoin),             intent(in)    :: qb_df
      real, dimension(3,G%npoin,inp%nlayers), intent(in)    :: qprime_df

      ! Zero rhs_btp_visc on device unconditionally so the SSPRK GPU kernel
      ! always reads a valid (zero) value even when method_visc == 0.
      !$acc kernels present(btp)
      btp%rhs_btp_visc = 0.0
      !$acc end kernels

      ! 1. MPI halo exchange (host).
      call btp_create_precommunicator(G, inp, b, mf, init, par, ref, mpic, qb_df, qprime_df, 4)

      ! 2. GPU volume kernel — zeros rhs_btp on device, accumulates volume contribution.
      call create_rhs_btp_volume_qdf(G, b, inp, init, btp, tsp, rhs_btp, qb_df, qprime_df)

      ! 3. GPU face kernel — accumulates local face fluxes into device rhs_btp.
      call create_btp_fluxes_qdf(G, b, inp, mf, init, btp, rhs_btp, qb_df, qprime_df)

      ! 4. Post-communicator: mpi_waitall → unpack → upload q_recv to GPU →
      !    GPU create_nbhs_face_df adds MPI boundary fluxes into device rhs_btp
      !    and device btp%*_face_ave.  rhs_btp stays on device until step 5.
      call btp_create_postcommunicator(G, inp, b, mf, par, btp, init, ref, mpic, rhs_btp, 4)

      ! 5. Viscosity (GPU-ported); rhs_btp_visc already zeroed on device above.
      if (inp%method_visc > 0) &
         call btp_create_laplacian(G, inp, b, mf, par, btp, init, ref, mpic, tsp, btp%rhs_btp_visc, qb_df)
      ! rhs_btp stays on device; mass-matrix scaling is fused into the SSPRK update in the caller.

   end subroutine create_rhs_btp


   subroutine create_rhs_btp_volume_qdf(G, b, inp, init, btp, tsp, rhs_btp, qb_df, qprime_df)
      !===========================================================================
      !  Volume contribution to the barotropic RHS (DG weak form, interior).
      !
      !  One gang per element processes all its nqx*nqy quad points sequentially,
      !  accumulates into a gang-private rhs_loc(3,npts), then scatters to global
      !  rhs_btp once with no atomics.
      !  DG nodes are element-local so the final scatter has no race conditions.
      !  Time-average arrays (btp%*_ave) are indexed by Iq which is unique per gang.
      !===========================================================================
      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_constants, only: gravity
      use mod_initial,   only: initial
      use mod_variables, only: btp_CS
      use mod_input,     only: input
      use mod_tensor,    only: tensor_CS

      implicit none

      type(grid),      intent(in)    :: G
      type(basis),     intent(in)    :: b
      type(input),     intent(in)    :: inp
      type(initial),   intent(in)    :: init
      type(btp_CS),    intent(inout) :: btp
      type(tensor_CS), intent(in)    :: tsp

      real, dimension(4,G%npoin),             intent(in)    :: qb_df
      real, dimension(3,G%npoin,inp%nlayers), intent(in)    :: qprime_df
      real, dimension(3,G%npoin),             intent(out) :: rhs_btp

      real :: sc_x, sc_y, Hq, Qu1, Qu2, Qv1, Qv2
      real :: wq, hi, dhdx, dhdy, tb_u, tb_v, ope, ope2, grav_dp
      real :: dp, dpp, udp, vdp, ub, vb, ubot, vbot, spd, pbq
      real :: pp_k, up_k, vp_k, pprime_k, pprime_k1
      real :: sum_up2, sum_uv, sum_vp2
      real :: wq_udp, wq_vdp, wq_scx, wq_scy, wq_Qu1, wq_Qu2, wq_Qv1, wq_Qv2
      real :: rhs_loc(3, b%npts)
      integer :: Iq, ie, iquad, jquad, ip, k, I, ijq
      integer :: npts_l, nqx_l, nqy_l, nelem_l, nlayers_l, botfr_l, npts_q
      real    :: cd_mlswe_l

      npts_l     = b%npts
      nqx_l      = b%nqx
      nqy_l      = b%nqy
      nelem_l    = G%nelem
      nlayers_l  = inp%nlayers
      botfr_l    = inp%botfr
      cd_mlswe_l = real(inp%cd_mlswe)
      npts_q = b%npts_quad

      !$acc data present(init, tsp, btp,                                             &
      !$acc               rhs_btp, qb_df, qprime_df,                                &
      !$acc               init%grad_zbot_quad, init%tau_wind, tsp%psih,             &
      !$acc               tsp%dpsidx, tsp%dpsidy,                                   &
      !$acc               tsp%indexq, tsp%indexq_e, tsp%wjac, init%pbprime_df,      &
      !$acc               init%coriolis_quad, init%alpha_mlswe,                     &
      !$acc               btp%tau_bot_ave, btp%H_ave, btp%Qu_ave, btp%Quv_ave,     &
      !$acc               btp%Qv_ave, btp%ope_ave,                                  &
      !$acc               btp%uvb_ave, btp%btp_mass_flux_ave, btp%ope2_ave)

      ! Zero rhs_btp on device.
      !$acc kernels
      rhs_btp = 0.0
      !$acc end kernels

      !$acc parallel loop gang                                                    &
      !$acc   private(rhs_loc,                                                   &
      !$acc           dp, dpp, udp, vdp, pbq, hi, Hq, wq, ub, vb,               &
      !$acc           ubot, vbot, spd, tb_u, tb_v, sc_x, sc_y,                   &
      !$acc           ope, ope2, grav_dp,                                         &
      !$acc           pp_k, up_k, vp_k, pprime_k, pprime_k1,                     &
      !$acc           Qu1, Qu2, Qv1, Qv2,                                        &
      !$acc           dhdx, dhdy, sum_up2, sum_uv, sum_vp2,                      &
      !$acc           wq_udp, wq_vdp, wq_scx, wq_scy,                            &
      !$acc           wq_Qu1, wq_Qu2, wq_Qv1, wq_Qv2,                           &
      !$acc           Iq, iquad, jquad, ip, k, I, ijq)
      do ie = 1, nelem_l

         ! Zero gang-private accumulator for this element's nodes.
         !$acc loop seq
         do ip = 1, npts_l
            rhs_loc(1, ip) = 0.0
            rhs_loc(2, ip) = 0.0
            rhs_loc(3, ip) = 0.0
         end do

         ! Loop over all nqx*nqy quadrature points of element ie sequentially.
         !$acc loop seq
         do ijq = 1, npts_q

            Iq = tsp%indexq_e(ijq, ie)
            wq = tsp%wjac(Iq)

            ! Barotropic projection
            dp = 0.0;  dpp = 0.0;  udp = 0.0;  vdp = 0.0;  pbq = 0.0
            !$acc loop seq
            do ip = 1, npts_l
               I = tsp%indexq(ip, Iq)
               hi = tsp%psih(ip, Iq)
               dp  = dp  + hi * qb_df(1, I)
               dpp = dpp + hi * qb_df(2, I)
               udp = udp + hi * qb_df(3, I)
               vdp = vdp + hi * qb_df(4, I)
               pbq = pbq + hi * init%pbprime_df(I)
            end do

            ub = udp / dp
            vb = vdp / dp

            ! Layer projections, Hq, momentum-flux sums
            Hq = 0.0;  sum_up2 = 0.0;  sum_uv = 0.0;  sum_vp2 = 0.0
            pprime_k = 0.0

            !$acc loop seq
            do k = 1, nlayers_l
               pp_k = 0.0;  up_k = 0.0;  vp_k = 0.0
               !$acc loop seq
               do ip = 1, npts_l
                  I = tsp%indexq(ip, Iq)
                  hi = tsp%psih(ip, Iq)
                  pp_k = pp_k + hi * qprime_df(1, I, k)
                  up_k = up_k + hi * qprime_df(2, I, k)
                  vp_k = vp_k + hi * qprime_df(3, I, k)
               end do

               pprime_k1 = pprime_k + pp_k

               Hq      = Hq      + 0.5 * init%alpha_mlswe(k) * (pprime_k1 + pprime_k) * pp_k
               sum_up2 = sum_up2 + up_k * up_k * pp_k
               sum_uv  = sum_uv  + up_k * vp_k * pp_k
               sum_vp2 = sum_vp2 + vp_k * vp_k * pp_k

               pprime_k = pprime_k1
            end do

            ! Bottom friction
            tb_u = 0.0;  tb_v = 0.0
            if (botfr_l /= 0) then
               ubot = up_k + ub
               vbot = vp_k + vb
               if (botfr_l == 1) then
                  spd = (cd_mlswe_l / gravity) * pp_k
               else
                  spd = (cd_mlswe_l / init%alpha_mlswe(nlayers_l)) * sqrt(ubot**2 + vbot**2)
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
            sc_x =  init%coriolis_quad(Iq) * vdp                           &
               + gravity * (init%tau_wind(1,Iq) - tb_u)                  &
               - grav_dp * init%grad_zbot_quad(1,Iq)
            sc_y = -init%coriolis_quad(Iq) * udp                           &
               + gravity * (init%tau_wind(2,Iq) - tb_v)                  &
               - grav_dp * init%grad_zbot_quad(2,Iq)

            ! Flux tensors
            Qu1 = ub * udp + ope * sum_up2
            Qu2 = vb * udp + ope * sum_uv
            Qv1 = ub * vdp + ope * sum_uv
            Qv2 = vb * vdp + ope * sum_vp2

            ! Time-average accumulators (Iq unique per gang — no race)
            btp%H_ave(Iq)               = btp%H_ave(Iq)               + Hq
            btp%Qu_ave(Iq)              = btp%Qu_ave(Iq)              + Qu1
            btp%Qv_ave(Iq)              = btp%Qv_ave(Iq)              + Qv1
            btp%Quv_ave(Iq)             = btp%Quv_ave(Iq)             + Qu2
            btp%tau_bot_ave(1,Iq)       = btp%tau_bot_ave(1,Iq)       + tb_u
            btp%tau_bot_ave(2,Iq)       = btp%tau_bot_ave(2,Iq)       + tb_v
            btp%ope_ave(Iq)             = btp%ope_ave(Iq)             + ope
            btp%ope2_ave(Iq)            = btp%ope2_ave(Iq)            + ope2
            btp%btp_mass_flux_ave(1,Iq) = btp%btp_mass_flux_ave(1,Iq) + udp
            btp%btp_mass_flux_ave(2,Iq) = btp%btp_mass_flux_ave(2,Iq) + vdp
            btp%uvb_ave(1,Iq)           = btp%uvb_ave(1,Iq)           + ub
            btp%uvb_ave(2,Iq)           = btp%uvb_ave(2,Iq)           + vb

            ! Add hydrostatic pressure to the diagonal flux components.
            Qu1 = Qu1 + Hq
            Qv2 = Qv2 + Hq

            wq_udp = wq * udp;   wq_vdp = wq * vdp
            wq_scx = wq * sc_x;  wq_scy = wq * sc_y
            wq_Qu1 = wq * Qu1;   wq_Qu2 = wq * Qu2
            wq_Qv1 = wq * Qv1;   wq_Qv2 = wq * Qv2

            ! Accumulate into gang-private array (no race: one gang per element).
            !$acc loop seq
            do ip = 1, npts_l
               hi   = tsp%psih(ip, Iq)
               dhdx = tsp%dpsidx(ip, Iq)
               dhdy = tsp%dpsidy(ip, Iq)
               rhs_loc(1, ip) = rhs_loc(1, ip) + dhdx*wq_udp + dhdy*wq_vdp
               rhs_loc(2, ip) = rhs_loc(2, ip) + hi*wq_scx + dhdx*wq_Qu1 + dhdy*wq_Qu2
               rhs_loc(3, ip) = rhs_loc(3, ip) + hi*wq_scy + dhdx*wq_Qv1 + dhdy*wq_Qv2
            end do

         end do

         ! Scatter local accumulator to global rhs_btp.
         ! No atomics: DG nodes belong to exactly one element.
         Iq = tsp%indexq_e(1, ie)
         !$acc loop seq
         do ip = 1, npts_l
            I = tsp%indexq(ip, Iq)
            rhs_btp(1, I) = rhs_loc(1, ip)
            rhs_btp(2, I) = rhs_loc(2, ip)
            rhs_btp(3, I) = rhs_loc(3, ip)
         end do

      end do
      !$acc end data

   end subroutine create_rhs_btp_volume_qdf


   subroutine create_btp_fluxes_qdf(G, b, inp, mf, init, btp, rhs_btp, qb, qprime_df)
      !===========================================================================
      !  Face flux contribution to the barotropic RHS (Riemann solver).
      !
      !  GPU parallelism: one gang per face.  The quadrature loop (iquad) and all
      !  inner loops (basis nodes n, layers k) are sequential within the gang,
      !  which avoids intra-gang atomics on nodes shared by adjacent iquad points.
      !
      !  Race conditions:
      !    rhs_btp(m,I): multiple faces share boundary nodes — !$acc atomic update.
      !    btp%*_face_ave(iquad,iface): each (iquad,iface) pair is touched by
      !      exactly one gang×sequential-step — no atomics needed.
      !===========================================================================
      use mod_basis,     only: basis
      use mod_grid,      only: grid
      use mod_face,      only: face_CS
      use mod_variables, only: btp_CS
      use mod_initial,   only: initial
      use mod_constants, only: gravity
      use mod_input,     only: input

      implicit none

      type(grid),    intent(in)    :: G
      type(basis),   intent(in)    :: b
      type(input),   intent(in)    :: inp
      type(face_CS), intent(in)    :: mf
      type(initial), intent(in)    :: init
      type(btp_CS),  intent(inout) :: btp

      real, dimension(3,G%npoin),             intent(inout) :: rhs_btp
      real, dimension(4,G%npoin),             intent(in)    :: qb
      real, dimension(3,G%npoin,inp%nlayers), intent(in)    :: qprime_df

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
      real,    dimension(b%ngl) :: hi_c
      integer, dimension(b%ngl) :: I_l, I_r

      ! Capture derived-type scalars before the parallel region so they are
      ! available as firstprivate values inside the GPU kernel.
      integer :: ngl_f, nq_f, nlayers_f, nface_f
      ngl_f     = b%ngl
      nq_f      = b%nq
      nlayers_f = inp%nlayers
      nface_f   = G%nface

      !$acc data present(btp, rhs_btp, qb, qprime_df,                              &
      !$acc               G%face, G%face_type, G%intma,                             &
      !$acc               mf%imapl, mf%imapr, mf%normal_vector_q, mf%jac_faceq,   &
      !$acc               b%psiq,                                                   &
      !$acc               init%pbprime_df, init%alpha_mlswe,                        &
      !$acc               btp%H_face_ave, btp%ope_face_ave, btp%ope2_face_ave,     &
      !$acc               btp%btp_mass_flux_face_ave,                               &
      !$acc               btp%Qu_face_ave, btp%Qv_face_ave,                        &
      !$acc               btp%one_plus_eta_edge_2_ave,                              &
      !$acc               btp%uvb_face_ave)

      !$acc parallel loop gang                                                      &
      !$acc   private(el, er, il, jl, kl, ir, jr, kr, I, n, k,                    &
      !$acc           I_l, I_r, hi_c,                                              &
      !$acc           qbl, qbr, pbl, pbr,                                          &
      !$acc           nxl, nyl, wq, hi, un,                                        &
      !$acc           clam, one_eta, one_eta2, half_clam, quarter_clam,            &
      !$acc           pU_L, pU_R, half_clam_dqb2, pbpert_edge,                    &
      !$acc           c_minus, c_plus, ul, ur, vl, vr, opl, opr,                  &
      !$acc           H_bcl_ql, H_bcl_qr, H_bcl_q, oe2_Hql, oe2_Hqr,             &
      !$acc           Qu_ql1, Qu_ql2, Qu_qr1, Qu_qr2,                             &
      !$acc           Qv_ql1, Qv_ql2, Qv_qr1, Qv_qr2,                             &
      !$acc           flux_edge_x, flux_edge_y, flux_pb, flux_u, flux_v,          &
      !$acc           pprime_lk, pprime_rk, pprime_lk1, pprime_rk1,               &
      !$acc           pkl, pkr, ukl, ukr, vkl, vkr,                               &
      !$acc           ope_ppl_k, ope_ppr_k, uv_cross_l, uv_cross_r,               &
      !$acc           wq_flux_pb, wq_flux_u, wq_flux_v)
      do iface = 1, nface_f

         if (G%face_type(iface) == 2) cycle

         el = G%face(7, iface)
         er = G%face(8, iface)

         ! Precompute node indices once per face; reused across all iquad iterations.
         !$acc loop seq
         do n = 1, ngl_f
            il = mf%imapl(1,n,1,iface)
            jl = mf%imapl(2,n,1,iface)
            kl = mf%imapl(3,n,1,iface)
            I_l(n) = G%intma(il, jl, kl, el)
         end do

         if (er > 0) then
            !$acc loop seq
            do n = 1, ngl_f
               ir = mf%imapr(1,n,1,iface)
               jr = mf%imapr(2,n,1,iface)
               kr = mf%imapr(3,n,1,iface)
               I_r(n) = G%intma(ir, jr, kr, er)
            end do
         end if

         !$acc loop seq
         do iquad = 1, nq_f

            nxl = mf%normal_vector_q(1, iquad, 1, iface)
            nyl = mf%normal_vector_q(2, iquad, 1, iface)

            !$acc loop seq
            do n = 1, ngl_f
               hi_c(n) = b%psiq(n, iquad)
            end do

            ! Project barotropic LEFT state
            qbl(1) = 0.0;  qbl(2) = 0.0;  qbl(3) = 0.0;  qbl(4) = 0.0
            pbl = 0.0
            !$acc loop seq
            do n = 1, ngl_f
               hi = hi_c(n)
               I  = I_l(n)
               qbl(1) = qbl(1) + hi * qb(1,I)
               qbl(2) = qbl(2) + hi * qb(2,I)
               qbl(3) = qbl(3) + hi * qb(3,I)
               qbl(4) = qbl(4) + hi * qb(4,I)
               pbl    = pbl    + hi * init%pbprime_df(I)
            end do

            ! Project barotropic RIGHT state or apply boundary condition.
            if (er > 0) then
               qbr(1) = 0.0;  qbr(2) = 0.0;  qbr(3) = 0.0;  qbr(4) = 0.0
               pbr = 0.0
               !$acc loop seq
               do n = 1, ngl_f
                  hi = hi_c(n)
                  I  = I_r(n)
                  qbr(1) = qbr(1) + hi * qb(1,I)
                  qbr(2) = qbr(2) + hi * qb(2,I)
                  qbr(3) = qbr(3) + hi * qb(3,I)
                  qbr(4) = qbr(4) + hi * qb(4,I)
                  pbr    = pbr    + hi * init%pbprime_df(I)
               end do
            else
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

            c_minus      = sqrt(init%alpha_mlswe(nlayers_f) * pbr)
            c_plus       = sqrt(init%alpha_mlswe(nlayers_f) * pbl)
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

            ! Layer projection + flux accumulation
            !$acc loop seq
            do k = 1, nlayers_f

               pkl = 0.0;  ukl = 0.0;  vkl = 0.0
               !$acc loop seq
               do n = 1, ngl_f
                  hi  = hi_c(n)
                  I   = I_l(n)
                  pkl = pkl + hi * qprime_df(1, I, k)
                  ukl = ukl + hi * qprime_df(2, I, k)
                  vkl = vkl + hi * qprime_df(3, I, k)
               end do

               if (er > 0) then
                  pkr = 0.0;  ukr = 0.0;  vkr = 0.0
                  !$acc loop seq
                  do n = 1, ngl_f
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
               H_bcl_ql   = H_bcl_ql  + 0.5*init%alpha_mlswe(k) * (pprime_lk1 + pprime_lk) * pkl
               H_bcl_qr   = H_bcl_qr  + 0.5*init%alpha_mlswe(k) * (pprime_rk1 + pprime_rk) * pkr
               pprime_lk  = pprime_lk1
               pprime_rk  = pprime_rk1

            end do  ! k

            H_bcl_q = 0.5 * (H_bcl_ql + H_bcl_qr)

            half_clam_dqb2 = half_clam * (qbl(2) - qbr(2))
            flux_edge_x = 0.5*(qbl(3) + qbr(3)) + nxl * half_clam_dqb2
            flux_edge_y = 0.5*(qbl(4) + qbr(4)) + nyl * half_clam_dqb2
            flux_pb     = nxl*flux_edge_x + nyl*flux_edge_y

            H_bcl_q = one_eta2 * H_bcl_q

            ! Time-average accumulators — (iquad,iface) unique per gang×seq step.
            btp%btp_mass_flux_face_ave(1,iquad,iface) = btp%btp_mass_flux_face_ave(1,iquad,iface) + flux_edge_x
            btp%btp_mass_flux_face_ave(2,iquad,iface) = btp%btp_mass_flux_face_ave(2,iquad,iface) + flux_edge_y
            btp%H_face_ave(iquad,iface)    = btp%H_face_ave(iquad,iface)    + H_bcl_q
            btp%Qu_face_ave(1,iquad,iface) = btp%Qu_face_ave(1,iquad,iface) + 0.5*(Qu_ql1 + Qu_qr1)
            btp%Qu_face_ave(2,iquad,iface) = btp%Qu_face_ave(2,iquad,iface) + 0.5*(Qu_ql2 + Qu_qr2)
            btp%Qv_face_ave(1,iquad,iface) = btp%Qv_face_ave(1,iquad,iface) + 0.5*(Qv_ql1 + Qv_qr1)
            btp%Qv_face_ave(2,iquad,iface) = btp%Qv_face_ave(2,iquad,iface) + 0.5*(Qv_ql2 + Qv_qr2)

            opl = 1.0 + qbl(2)/pbl
            opr = 1.0 + qbr(2)/pbr
            btp%ope_face_ave(1,iquad,iface)  = btp%ope_face_ave(1,iquad,iface)  + opl
            btp%ope_face_ave(2,iquad,iface)  = btp%ope_face_ave(2,iquad,iface)  + opr
            btp%ope2_face_ave(1,iquad,iface) = btp%ope2_face_ave(1,iquad,iface) + opl*opl
            btp%ope2_face_ave(2,iquad,iface) = btp%ope2_face_ave(2,iquad,iface) + opr*opr
            btp%one_plus_eta_edge_2_ave(iquad,iface) = btp%one_plus_eta_edge_2_ave(iquad,iface) + one_eta2
            btp%uvb_face_ave(1,1,iquad,iface) = btp%uvb_face_ave(1,1,iquad,iface) + ul
            btp%uvb_face_ave(1,2,iquad,iface) = btp%uvb_face_ave(1,2,iquad,iface) + ur
            btp%uvb_face_ave(2,1,iquad,iface) = btp%uvb_face_ave(2,1,iquad,iface) + vl
            btp%uvb_face_ave(2,2,iquad,iface) = btp%uvb_face_ave(2,2,iquad,iface) + vr

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

            wq         = mf%jac_faceq(iquad, 1, iface)
            wq_flux_pb = wq * flux_pb
            wq_flux_u  = wq * flux_u
            wq_flux_v  = wq * flux_v

            ! RHS scatter — atomics: multiple faces share boundary nodes.
            !$acc loop seq
            do n = 1, ngl_f
               hi = hi_c(n)
               I  = I_l(n)
               !$acc atomic update
               rhs_btp(1,I) = rhs_btp(1,I) - hi * wq_flux_pb
               !$acc atomic update
               rhs_btp(2,I) = rhs_btp(2,I) - hi * wq_flux_u
               !$acc atomic update
               rhs_btp(3,I) = rhs_btp(3,I) - hi * wq_flux_v
            end do

            if (er > 0) then
               !$acc loop seq
               do n = 1, ngl_f
                  hi = hi_c(n)
                  I  = I_r(n)
                  !$acc atomic update
                  rhs_btp(1,I) = rhs_btp(1,I) + hi * wq_flux_pb
                  !$acc atomic update
                  rhs_btp(2,I) = rhs_btp(2,I) + hi * wq_flux_u
                  !$acc atomic update
                  rhs_btp(3,I) = rhs_btp(3,I) + hi * wq_flux_v
               end do
            end if

         end do  ! iquad

      end do  ! iface
      !$acc end data

   end subroutine create_btp_fluxes_qdf

end module mod_rhs_btp
