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

      real, dimension(inp%nvar_btp-1,G%npoin),             intent(out) :: rhs_btp
      real, dimension(inp%nvar_btp,G%npoin),               intent(in)  :: qb_df
      real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers),   intent(in)  :: qprime_df

      ! 1. MPI halo exchange (host).
      call btp_create_precommunicator(G, inp, b, mf, init, par, ref, mpic, qb_df, qprime_df, inp%nvar_btp)

      ! 2. GPU volume kernel — zeros rhs_btp on device, accumulates volume contribution.
#ifdef _OPENACC
      call create_rhs_btp_volume_qdf_qp_gpu(G, b, inp, init, btp, tsp, mt, rhs_btp, qb_df)
#else
      call create_rhs_btp_volume_qdf_qp(G, b, inp, init, btp, tsp, rhs_btp, qb_df)
#endif

      ! 3. GPU face kernel — accumulates local face fluxes into device rhs_btp.
#ifdef _OPENACC
      call create_btp_fluxes_qdf_gpu(G, b, inp, mf, init, btp, rhs_btp, qb_df, qprime_df)
#else
      call create_btp_fluxes_qdf(G, b, inp, mf, init, btp, rhs_btp, qb_df, qprime_df)
#endif

      ! 4. Post-communicator: mpi_waitall → unpack → upload q_recv to GPU →
      !    GPU create_nbhs_face_df adds MPI boundary fluxes into device rhs_btp
      !    and device btp%*_face_ave.  rhs_btp stays on device until step 5.
      call btp_create_postcommunicator(G, inp, b, mf, par, btp, init, ref, mpic, rhs_btp, inp%nvar_btp-1)
      ! rhs_btp stays on device; mass-matrix scaling is fused into the SSPRK update in the caller.

   end subroutine create_rhs_btp


   subroutine create_rhs_btp_volume_qdf(G, b, inp, init, btp, tsp, rhs_btp, qb_df)
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

      real, dimension(4,G%npoin), intent(in)  :: qb_df
      real, dimension(3,G%npoin), intent(out) :: rhs_btp

      real :: sc_x, sc_y, Hq, Qu1, Qu2, Qv1, Qv2
      real :: wq, hi, dhdx, dhdy, tb_u, tb_v, ope, ope2, grav_dp
      real :: dp, dpp, udp, vdp, ub, vb, ubot, vbot, spd, pbq
      real :: pp_k, up_k, vp_k
      real :: sum_up2, sum_uv, sum_vp2
      real :: wq_udp, wq_vdp, wq_scx, wq_scy, wq_Qu1, wq_Qu2, wq_Qv1, wq_Qv2
      real :: rhs_loc(3, b%npts)
      integer :: Iq, ie, iquad, jquad, ip, k, I, ijq
      integer :: npts_l, nqx_l, nqy_l, nelem_l, botfr_l, npts_q
      real    :: cd_mlswe_l, alpha_bot_l

      npts_l     = b%npts
      nqx_l      = b%nqx
      nqy_l      = b%nqy
      nelem_l    = G%nelem
      botfr_l    = inp%botfr
      cd_mlswe_l = real(inp%cd_mlswe)
      alpha_bot_l = real(init%alpha_mlswe(inp%nlayers))
      npts_q = b%npts_quad

      !$acc data present(init, tsp, btp,                                             &
      !$acc               rhs_btp, qb_df,                                            &
      !$acc               init%grad_zbot_quad, init%tau_wind, tsp%psih,             &
      !$acc               tsp%dpsidx, tsp%dpsidy,                                   &
      !$acc               tsp%indexq, tsp%indexq_e, tsp%wjac,                       &
      !$acc               init%coriolis_quad,                                        &
      !$acc               btp%tau_bot_ave, btp%H_ave, btp%Qu_ave,                   &
      !$acc               btp%Qv_ave, btp%ope_ave,                                  &
      !$acc               btp%uvb_ave, btp%btp_mass_flux_ave, btp%ope2_ave,        &
      !$acc               btp%bcl_H, btp%bcl_flux, btp%bcl_btp_flux, btp%pbq)

      ! Zero rhs_btp on device.
      !$acc kernels
      rhs_btp = 0.0
      !$acc end kernels

      !$acc parallel loop gang                                                    &
      !$acc   private(rhs_loc,                                                   &
      !$acc           dp, dpp, udp, vdp, pbq, hi, Hq, wq, ub, vb,               &
      !$acc           ubot, vbot, spd, tb_u, tb_v, sc_x, sc_y,                   &
      !$acc           ope, ope2, grav_dp,                                         &
      !$acc           pp_k, up_k, vp_k,                                          &
      !$acc           Qu1, Qu2, Qv1, Qv2,                                        &
      !$acc           dhdx, dhdy, sum_up2, sum_uv, sum_vp2,                      &
      !$acc           wq_udp, wq_vdp, wq_scx, wq_scy,                            &
      !$acc           wq_Qu1, wq_Qu2, wq_Qv1, wq_Qv2,                           &
      !$acc           Iq, iquad, jquad, ip, I, ijq)                              &
      !$acc   firstprivate(npts_l, botfr_l, cd_mlswe_l, alpha_bot_l, npts_q)
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

            ! Barotropic projection (qb_df changes each substep).
            dp = 0.0;  dpp = 0.0;  udp = 0.0;  vdp = 0.0
            !$acc loop seq
            do ip = 1, npts_l
               I = tsp%indexq(ip, Iq)
               hi = tsp%psih(ip, Iq)
               dp  = dp  + hi * qb_df(1, I)
               dpp = dpp + hi * qb_df(2, I)
               udp = udp + hi * qb_df(3, I)
               vdp = vdp + hi * qb_df(4, I)
            end do

            ub = udp / dp
            vb = vdp / dp

            ! Load BCL layer integrals precomputed once per RK stage.
            Hq      = btp%bcl_H(Iq)
            sum_up2 = btp%bcl_flux(1,1,Iq)
            sum_uv  = btp%bcl_flux(1,2,Iq)
            sum_vp2 = btp%bcl_flux(2,2,Iq)
            pp_k    = btp%bcl_btp_flux(1,Iq)
            up_k    = btp%bcl_btp_flux(2,Iq)
            vp_k    = btp%bcl_btp_flux(3,Iq)
            pbq     = btp%pbq(Iq)

            ! Bottom friction
            tb_u = 0.0;  tb_v = 0.0
            if (botfr_l /= 0) then
               ubot = up_k + ub
               vbot = vp_k + vb
               if (botfr_l == 1) then
                  spd = (cd_mlswe_l / gravity) * pp_k
               else
                  spd = (cd_mlswe_l / alpha_bot_l) * sqrt(ubot**2 + vbot**2)
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
            btp%H_ave(Iq)  = btp%H_ave(Iq)  + Hq
            btp%Qu_ave(1,Iq) = btp%Qu_ave(1,Iq) + Qu1
            btp%Qu_ave(2,Iq) = btp%Qu_ave(2,Iq) + Qu2
            btp%Qv_ave(1,Iq) = btp%Qv_ave(1,Iq) + Qv1
            btp%Qv_ave(2,Iq) = btp%Qv_ave(2,Iq) + Qv2

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

   subroutine create_rhs_btp_volume_qdf_qp_gpu(G, b, inp, init, btp, tsp, mt, rhs_btp, qb_df)
      !===========================================================================
      !  Volume contribution — v1 physics.
      !  BCL layer integrals (H_b, sum_up2, etc.) are precomputed once per RK
      !  stage by btp_bcl_coeffs_qdf into btp%bcl_* arrays; only qb_df
      !  interpolation is done per substep.
      !
      !  Four kernels, no atomics:
      !    A0) one thread per (element, node row m, quad column iqx): qb_df summed
      !        along xi into btp%qb_T, the first pass of the tensor-product
      !        interpolation to quad points;
      !    A) one thread per quad point: the second (eta) pass of the
      !       interpolation from btp%qb_T, then: the physics at that point, combined
      !       with the metric terms and stored as weighted coefficients of the
      !       reference-element basis (psi, dpsi/dxi, dpsi/deta) in
      !       btp%F_vol(Iq,:,:);
      !    B1) one thread per (element, node row m, quad column iqx): the eta
      !        pass of the test-function sum, into btp%qb_T (free once kernel A
      !        has run): for each row r, A_r = sum_iqy psiy*F(r,1) + dpsiy*F(r,3)
      !        and B_r = sum_iqy psiy*F(r,2);
      !    B2) one thread per DOF node: the xi pass,
      !        rhs(r) = sum_iqx psix*A_r + dpsix*B_r, written once per node (DG
      !        nodes belong to exactly one element).
      !  Separate kernels (rather than two vector loops in one gang) avoid the
      !  block barrier between the phases and let kernel B run with far fewer
      !  registers. F_vol has the quad-point index first so that neighbouring
      !  threads in kernel A store to neighbouring addresses.
      !
      !  The physical gradient is dpsi/dx = dpsi/dxi*ksi_x + dpsi/deta*eta_x
      !  + psi*zeta_x (likewise y, z), as in compute_tensor_product. Kernel A
      !  folds the metric terms into the fluxes once per quad point, so kernel B
      !  only needs the 1D basis tables (psiqx/psiqy, dpsiqx/dpsiqy) instead of
      !  the metric terms or the large (npts, npoin_q) tsp%psih/dpsid* arrays.
      !===========================================================================
      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_constants, only: gravity
      use mod_initial,   only: initial
      use mod_variables, only: btp_CS
      use mod_input,     only: input
      use mod_tensor,    only: tensor_CS
      use mod_metrics,   only: metrics
      use mod_initial_mlswe, only: check_btp_thickness

      implicit none

      type(grid),      intent(in)    :: G
      type(basis),     intent(in)    :: b
      type(input),     intent(in)    :: inp
      type(initial),   intent(in)    :: init
      type(btp_CS),    intent(inout) :: btp
      type(tensor_CS), intent(in)    :: tsp
      type(metrics),   intent(in)    :: mt

      real, dimension(inp%nvar_btp,G%npoin),   intent(in)  :: qb_df
      real, dimension(inp%nvar_btp-1,G%npoin), intent(out) :: rhs_btp

      real :: sc_x, sc_y, sc_z, Hq, qu, quv, qvu, qv, quz, qvz, qwu, qwv, qw, H_b
      real :: wq, hi, tb_u, tb_v, tb_w, ope, ope2
      real :: dp, dpp, udp, vdp, wdp, ub, vb, wb, ubot, vbot, wbot, spd, pbq
      real :: pp_k, up_k, vp_k, wp_k
      real :: sum_up2, sum_uv, sum_vu, sum_vp2, sum_uw, sum_vw, sum_ww, sum_wu, sum_wv
      real :: cfx, cfy, cfz
      real :: r1, r2, r3, r4
      real :: h_e, h_n
      real :: kx, ky, kz, ex, ey, ez, zx, zy, zz
      integer :: I, Iq, ip, ie, ijq, Iq0, n, m, iqx, iqy
      integer :: npts_l, npts_q_l, nelem_l, botfr_l, nglx_l, nqx_l, ngly_l, nqy_l, k, k0
      real    :: a1, a2, a3, a4, b1, b2, b3, b4, py, dpy, px, dpx
      real    :: cd_l, alpha_bot_l
      logical :: has_w

      npts_l      = b%npts
      npts_q_l    = b%npts_quad
      nelem_l     = G%nelem
      nglx_l      = b%nglx
      nqx_l       = b%nqx
      ngly_l      = b%ngly
      nqy_l       = b%nqy
      botfr_l     = inp%botfr
      cd_l        = real(inp%cd_mlswe)
      alpha_bot_l = real(init%alpha_mlswe(inp%nlayers))
      has_w       = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex

      !$acc data present(init, tsp, btp, G%intma, G%intma_dg_quad,                  &
      !$acc               rhs_btp, qb_df,                                            &
      !$acc               init%grad_zbot_quad, init%tau_wind, tsp%psih,             &
      !$acc               tsp%dpsidx, tsp%dpsidy, tsp%dpsidz,                       &
      !$acc               tsp%dpsidz_x, tsp%dpsidz_y, tsp%dpsidz_z,                &
      !$acc               tsp%indexq, tsp%indexq_e, tsp%wjac, btp%F_vol,           &
      !$acc               btp%qb_T,                                                  &
      !$acc               b%psiqx, b%psiqy, b%dpsiqx, b%dpsiqy,                    &
      !$acc               mt%ksiq_x, mt%ksiq_y, mt%ksiq_z,                          &
      !$acc               mt%etaq_x, mt%etaq_y, mt%etaq_z,                          &
      !$acc               mt%zetaq_x, mt%zetaq_y, mt%zetaq_z,                       &
      !$acc               init%coriolis_3d_quad,                                     &
      !$acc               btp%tau_bot_ave, btp%H_ave, btp%Qu_ave,                   &
      !$acc               btp%Qv_ave, btp%Qw_ave, btp%ope_ave,                     &
      !$acc               btp%uvb_ave, btp%btp_mass_flux_ave, btp%ope2_ave,        &
      !$acc               btp%bcl_H, btp%bcl_flux, btp%bcl_btp_flux, btp%pbq)

      ! Abort with element/rank diagnostics rather than SIGFPE-ing on ub=udp/dp
      ! below if the barotropic depth has gone negative or NaN (e.g. interface
      ! outcropping near steep topography feeding through into qb_df).
      call check_btp_thickness(b, G, inp, qb_df, 'create_rhs_btp_volume_qdf_qp_gpu')

      ! ---- Kernel A0: interpolate qb_df along xi --------------------------------
      ! qb_T(k,v) = sum_n psiqx(n,iqx) * qb_df(v, node(n,m)) for each element,
      ! node row m and quad column iqx; k runs iqx fastest so that neighbouring
      ! threads store to neighbouring addresses.
      !$acc parallel loop async(1) gang vector collapse(3)                                   &
      !$acc   firstprivate(nelem_l, nglx_l, ngly_l, nqx_l, has_w)                    &
      !$acc   private(Iq0, I, n, k, hi, dp, dpp, udp, vdp, wdp)
      do ie = 1, nelem_l
         do m = 1, ngly_l
            do iqx = 1, nqx_l
               dp = 0.0;  dpp = 0.0;  udp = 0.0;  vdp = 0.0;  wdp = 0.0
               !$acc loop seq
               do n = 1, nglx_l
                  I   = G%intma(n, m, 1, ie)
                  hi  = b%psiqx(n, iqx)
                  dp  = dp  + hi * qb_df(1, I)
                  dpp = dpp + hi * qb_df(2, I)
                  udp = udp + hi * qb_df(3, I)
                  vdp = vdp + hi * qb_df(4, I)
                  if (has_w) wdp = wdp + hi * qb_df(5, I)
               end do
               k = iqx + (m - 1)*nqx_l + (ie - 1)*nqx_l*ngly_l
               btp%qb_T(k, 1) = dp
               btp%qb_T(k, 2) = dpp
               btp%qb_T(k, 3) = udp
               btp%qb_T(k, 4) = vdp
               if (has_w) btp%qb_T(k, 5) = wdp
            end do
         end do
      end do
      !$acc end parallel loop

      ! ---- Kernel A: physics at each quad point (one thread per point) --------
      !$acc parallel loop async(1) gang vector collapse(3)                                   &
      !$acc   firstprivate(nelem_l, botfr_l, cd_l, alpha_bot_l, has_w,             &
      !$acc                ngly_l, nqx_l, nqy_l)                                    &
      !$acc   private(Iq, Iq0, dp, dpp, udp, vdp, wdp, pbq, H_b, Hq, wq, ub, vb, wb, &
      !$acc           ubot, vbot, wbot, spd, tb_u, tb_v, tb_w,                    &
      !$acc           cfx, cfy, cfz, sc_x, sc_y, sc_z, ope, ope2,                 &
      !$acc           pp_k, up_k, vp_k, wp_k,                                      &
      !$acc           sum_up2, sum_uv, sum_vu, sum_vp2, sum_uw, sum_vw, sum_ww,  &
      !$acc           sum_wu, sum_wv,                                              &
      !$acc           qu, quv, qvu, qv, quz, qvz, qwu, qwv, qw,                   &
      !$acc           kx, ky, kz, ex, ey, ez, zx, zy, zz,                         &
      !$acc           hi, k, k0, m)
      do ie = 1, nelem_l
         do iqy = 1, nqy_l
         do iqx = 1, nqx_l

            ! (iqx, iqy): position of this quad point in the element; indexq_e
            ! runs iqx fastest, matching compute_tensor_product.
            Iq  = G%intma_dg_quad(iqx, iqy, 1, ie)

            dp = 0.0;  dpp = 0.0;  udp = 0.0;  vdp = 0.0;  wdp = 0.0

            ! Barotropic state at the quad point: second (eta) pass of the
            ! tensor-product interpolation, from kernel A0's xi sums.
            k0 = iqx + (ie - 1)*nqx_l*ngly_l
            !$acc loop seq
            do m = 1, ngly_l
               k   = k0 + (m - 1)*nqx_l
               hi  = b%psiqy(m, iqy)
               dp  = dp  + hi * btp%qb_T(k, 1)
               dpp = dpp + hi * btp%qb_T(k, 2)
               udp = udp + hi * btp%qb_T(k, 3)
               vdp = vdp + hi * btp%qb_T(k, 4)
               if (has_w) wdp = wdp + hi * btp%qb_T(k, 5)
            end do

            ! Load BCL layer integrals precomputed once per RK stage.
            H_b     = btp%bcl_H(Iq)
            sum_up2 = btp%bcl_flux(1,1,Iq)
            sum_uv  = btp%bcl_flux(1,2,Iq)
            sum_vu  = btp%bcl_flux(2,1,Iq)
            sum_vp2 = btp%bcl_flux(2,2,Iq)
            pp_k    = btp%bcl_btp_flux(1,Iq)
            up_k    = btp%bcl_btp_flux(2,Iq)
            vp_k    = btp%bcl_btp_flux(3,Iq)
            pbq     = btp%pbq(Iq)

            sum_uw = 0.0;  sum_vw = 0.0;  sum_ww = 0.0;  wp_k = 0.0
            sum_wu = 0.0;  sum_wv = 0.0
            if (has_w) then
               sum_uw = btp%bcl_flux(3,1,Iq)
               sum_vw = btp%bcl_flux(3,2,Iq)
               sum_wu = btp%bcl_flux(1,3,Iq)
               sum_wv = btp%bcl_flux(2,3,Iq)
               sum_ww = btp%bcl_flux(3,3,Iq)
               wp_k   = btp%bcl_btp_flux(4,Iq)
            end if

            wq = tsp%wjac(Iq)
            ub = udp / dp;  vb = vdp / dp;  wb = wdp / dp

            ! 3D Coriolis vector (f * r̂); reduces to the f-plane pair (0,0,f) on
            ! cartesian since kvector=(0,0,1) there — see wind_stress_coriolis.
            cfx = init%coriolis_3d_quad(1, Iq)
            cfy = init%coriolis_3d_quad(2, Iq)
            cfz = init%coriolis_3d_quad(3, Iq)

            tb_u = 0.0;  tb_v = 0.0;  tb_w = 0.0
            if (botfr_l == 1) then
               ubot = up_k + ub;  vbot = vp_k + vb;  wbot = wp_k + wb
               spd  = (cd_l / gravity) * pp_k
               tb_u = spd * ubot;  tb_v = spd * vbot;  tb_w = spd * wbot
            elseif (botfr_l == 2) then
               ubot = up_k + ub;  vbot = vp_k + vb;  wbot = wp_k + wb
               spd  = (cd_l / alpha_bot_l) * sqrt(ubot**2 + vbot**2 + wbot**2)
               tb_u = spd * ubot;  tb_v = spd * vbot;  tb_w = spd * wbot
            end if

            ! Source terms (3D Coriolis cross-product -f×p*u; reduces exactly to
            ! the old f-plane form on cartesian since cfx=cfy=0 and wdp=0 there).
            sc_x = -(cfy*wdp - cfz*vdp)                                  &
               + gravity * (init%tau_wind(1, Iq) - tb_u)                 &
               - gravity * dp * init%grad_zbot_quad(1, Iq)
            sc_y = -(cfz*udp - cfx*wdp)                                  &
               + gravity * (init%tau_wind(2, Iq) - tb_v)                 &
               - gravity * dp * init%grad_zbot_quad(2, Iq)

            ope  = 1.0 + (dpp / pbq)
            ope2 = ope * ope
            Hq   = ope2 * H_b

            qu  = ub * udp + ope * sum_up2
            quv = vb * udp + ope * sum_uv
            qvu = ub * vdp + ope * sum_vu
            qv  = vb * vdp + ope * sum_vp2

            ! Time-average accumulators (each Iq is unique per gang — no race).
            btp%H_ave(Iq)    = btp%H_ave(Iq)    + Hq
            btp%Qu_ave(1,Iq) = btp%Qu_ave(1,Iq) + qu
            btp%Qu_ave(2,Iq) = btp%Qu_ave(2,Iq) + quv
            btp%Qv_ave(1,Iq) = btp%Qv_ave(1,Iq) + qvu
            btp%Qv_ave(2,Iq) = btp%Qv_ave(2,Iq) + qv

            btp%tau_bot_ave(1,Iq)       = btp%tau_bot_ave(1,Iq)       + tb_u
            btp%tau_bot_ave(2,Iq)       = btp%tau_bot_ave(2,Iq)       + tb_v
            btp%ope_ave(Iq)             = btp%ope_ave(Iq)             + ope
            btp%ope2_ave(Iq)            = btp%ope2_ave(Iq)            + ope2
            btp%btp_mass_flux_ave(1,Iq) = btp%btp_mass_flux_ave(1,Iq) + udp
            btp%btp_mass_flux_ave(2,Iq) = btp%btp_mass_flux_ave(2,Iq) + vdp
            btp%uvb_ave(1,Iq)           = btp%uvb_ave(1,Iq)           + ub
            btp%uvb_ave(2,Iq)           = btp%uvb_ave(2,Iq)           + vb

            ! Sphere-only: z-momentum flux/source terms and their time averages.
            sc_z = 0.0;  quz = 0.0;  qvz = 0.0;  qwu = 0.0;  qwv = 0.0;  qw = 0.0
            if (has_w) then
               sc_z = -(cfx*vdp - cfy*udp)                                &
                    - gravity * tb_w                                        &
                    - gravity * dp * init%grad_zbot_quad(3, Iq)

               quz = wb * udp + ope * sum_uw
               qvz = wb * vdp + ope * sum_vw
               qwu = ub * wdp + ope * sum_wu
               qwv = vb * wdp + ope * sum_wv
               qw  = wb * wdp + ope * sum_ww

               btp%Qu_ave(3,Iq)            = btp%Qu_ave(3,Iq)            + quz
               btp%Qv_ave(3,Iq)            = btp%Qv_ave(3,Iq)            + qvz
               btp%Qw_ave(1,Iq)            = btp%Qw_ave(1,Iq)            + qwu
               btp%Qw_ave(2,Iq)            = btp%Qw_ave(2,Iq)            + qwv
               btp%Qw_ave(3,Iq)            = btp%Qw_ave(3,Iq)            + qw
               btp%tau_bot_ave(3,Iq)       = btp%tau_bot_ave(3,Iq)       + tb_w
               btp%btp_mass_flux_ave(3,Iq) = btp%btp_mass_flux_ave(3,Iq) + wdp
               btp%uvb_ave(3,Iq)           = btp%uvb_ave(3,Iq)           + wb
            end if

            ! Metric terms at this quad point (zeta_* and the z components are
            ! only used on sphere_hex; zero otherwise).
            kx = mt%ksiq_x(iqx, iqy, 1, ie);  ky = mt%ksiq_y(iqx, iqy, 1, ie)
            ex = mt%etaq_x(iqx, iqy, 1, ie);  ey = mt%etaq_y(iqx, iqy, 1, ie)
            kz = 0.0;  ez = 0.0;  zx = 0.0;  zy = 0.0;  zz = 0.0
            if (has_w) then
               kz = mt%ksiq_z(iqx, iqy, 1, ie);  ez = mt%etaq_z(iqx, iqy, 1, ie)
               zx = mt%zetaq_x(iqx, iqy, 1, ie)
               zy = mt%zetaq_y(iqx, iqy, 1, ie)
               zz = mt%zetaq_z(iqx, iqy, 1, ie)
            end if

            ! Weighted coefficients of psi (j=1), dpsi/dxi (j=2) and dpsi/deta
            ! (j=3) for each equation row at this quad point. With source S and
            ! flux F = (Fx, Fy, Fz): j=1 is S + zeta.F, j=2 is ksi.F, j=3 is eta.F.
            btp%F_vol(Iq,1,1) = wq*(zx*udp + zy*vdp + zz*wdp)
            btp%F_vol(Iq,1,2) = wq*(kx*udp + ky*vdp + kz*wdp)
            btp%F_vol(Iq,1,3) = wq*(ex*udp + ey*vdp + ez*wdp)
            btp%F_vol(Iq,2,1) = wq*(sc_x + zx*(Hq+qu) + zy*quv + zz*quz)
            btp%F_vol(Iq,2,2) = wq*(kx*(Hq+qu) + ky*quv + kz*quz)
            btp%F_vol(Iq,2,3) = wq*(ex*(Hq+qu) + ey*quv + ez*quz)
            btp%F_vol(Iq,3,1) = wq*(sc_y + zx*qvu + zy*(Hq+qv) + zz*qvz)
            btp%F_vol(Iq,3,2) = wq*(kx*qvu + ky*(Hq+qv) + kz*qvz)
            btp%F_vol(Iq,3,3) = wq*(ex*qvu + ey*(Hq+qv) + ez*qvz)
            btp%F_vol(Iq,4,1) = wq*(sc_z + zx*qwu + zy*qwv + zz*(Hq+qw))
            btp%F_vol(Iq,4,2) = wq*(kx*qwu + ky*qwv + kz*(Hq+qw))
            btp%F_vol(Iq,4,3) = wq*(ex*qwu + ey*qwv + ez*(Hq+qw))

         end do
         end do
      end do
      !$acc end parallel loop

      ! ---- Kernel B1: eta pass of the test-function sum ------------------------
      ! For each element, node row m and quad column iqx (same index k as A0):
      !   qb_T(k,2r-1) = sum_iqy psiqy(m,iqy)*F_vol(Iq,r,1) + dpsiqy(m,iqy)*F_vol(Iq,r,3)
      !   qb_T(k,2r)   = sum_iqy psiqy(m,iqy)*F_vol(Iq,r,2)
      !$acc parallel loop async(1) gang vector collapse(3)                                   &
      !$acc   firstprivate(nelem_l, ngly_l, nqx_l, nqy_l, has_w)                     &
      !$acc   private(Iq, ijq, k, iqy, py, dpy, a1, a2, a3, a4, b1, b2, b3, b4)
      do ie = 1, nelem_l
         do m = 1, ngly_l
            do iqx = 1, nqx_l
               a1 = 0.0;  a2 = 0.0;  a3 = 0.0;  a4 = 0.0
               b1 = 0.0;  b2 = 0.0;  b3 = 0.0;  b4 = 0.0
               !$acc loop seq
               do iqy = 1, nqy_l
                  ijq = iqx + (iqy - 1)*nqx_l
                  Iq  = tsp%indexq_e(ijq, ie)
                  py  = b%psiqy(m, iqy)
                  dpy = b%dpsiqy(m, iqy)
                  a1 = a1 + py*btp%F_vol(Iq,1,1) + dpy*btp%F_vol(Iq,1,3)
                  b1 = b1 + py*btp%F_vol(Iq,1,2)
                  a2 = a2 + py*btp%F_vol(Iq,2,1) + dpy*btp%F_vol(Iq,2,3)
                  b2 = b2 + py*btp%F_vol(Iq,2,2)
                  a3 = a3 + py*btp%F_vol(Iq,3,1) + dpy*btp%F_vol(Iq,3,3)
                  b3 = b3 + py*btp%F_vol(Iq,3,2)
                  if (has_w) then
                     a4 = a4 + py*btp%F_vol(Iq,4,1) + dpy*btp%F_vol(Iq,4,3)
                     b4 = b4 + py*btp%F_vol(Iq,4,2)
                  end if
               end do
               k = iqx + (m - 1)*nqx_l + (ie - 1)*nqx_l*ngly_l
               btp%qb_T(k,1) = a1;  btp%qb_T(k,2) = b1
               btp%qb_T(k,3) = a2;  btp%qb_T(k,4) = b2
               btp%qb_T(k,5) = a3;  btp%qb_T(k,6) = b3
               btp%qb_T(k,7) = a4;  btp%qb_T(k,8) = b4
            end do
         end do
      end do
      !$acc end parallel loop

      ! ---- Kernel B2: xi pass, one thread per DOF node -------------------------
      ! No atomics: each DG node belongs to exactly one element.
      !$acc parallel loop async(1) gang vector collapse(3)                                   &
      !$acc   firstprivate(nelem_l, nglx_l, nqx_l, ngly_l, has_w)                    &
      !$acc   private(I, k, k0, iqx, px, dpx, r1, r2, r3, r4)
      do ie = 1, nelem_l
         do m = 1, ngly_l
         do n = 1, nglx_l
            I  = G%intma(n, m, 1, ie)
            k0 = (m - 1)*nqx_l + (ie - 1)*nqx_l*ngly_l
            r1 = 0.0;  r2 = 0.0;  r3 = 0.0;  r4 = 0.0
            !$acc loop seq
            do iqx = 1, nqx_l
               k   = iqx + k0
               px  = b%psiqx(n, iqx)
               dpx = b%dpsiqx(n, iqx)
               r1 = r1 + px*btp%qb_T(k,1) + dpx*btp%qb_T(k,2)
               r2 = r2 + px*btp%qb_T(k,3) + dpx*btp%qb_T(k,4)
               r3 = r3 + px*btp%qb_T(k,5) + dpx*btp%qb_T(k,6)
               if (has_w) r4 = r4 + px*btp%qb_T(k,7) + dpx*btp%qb_T(k,8)
            end do
            rhs_btp(1, I) = r1
            rhs_btp(2, I) = r2
            rhs_btp(3, I) = r3
            if (has_w) rhs_btp(4, I) = r4
         end do
         end do
      end do
      !$acc end parallel loop
      !$acc end data

   end subroutine create_rhs_btp_volume_qdf_qp_gpu

   subroutine create_rhs_btp_volume_qdf_qp(G, b, inp, init, btp, tsp, rhs_btp, qb_df)
      !===========================================================================
      !  CPU version of create_rhs_btp_volume_qdf_qp_gpu (used when built without
      !  OpenACC). Same physics, one pass over quad points with a direct scatter
      !  into rhs_btp: on the CPU the scatter needs no atomics, so this avoids the
      !  F_vol pass and basis recomputation the GPU version uses to remove them.
      !
      !  Volume contribution — v1 physics.
      !  BCL layer integrals (H_b, sum_up2, etc.) are precomputed once per RK
      !  stage by btp_bcl_coeffs_qdf into btp%bcl_* arrays; only qb_df
      !  interpolation is done per substep.
      !  Gang-parallel over do Iq = 1, npoin_q; scatter uses !$acc atomic because
      !  multiple Iq within the same element share DG DOF nodes.
      !===========================================================================
      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_constants, only: gravity
      use mod_initial,   only: initial
      use mod_variables, only: btp_CS
      use mod_input,     only: input
      use mod_tensor,    only: tensor_CS
      use mod_initial_mlswe, only: check_btp_thickness

      implicit none

      type(grid),      intent(in)    :: G
      type(basis),     intent(in)    :: b
      type(input),     intent(in)    :: inp
      type(initial),   intent(in)    :: init
      type(btp_CS),    intent(inout) :: btp
      type(tensor_CS), intent(in)    :: tsp

      real, dimension(inp%nvar_btp,G%npoin),   intent(in)  :: qb_df
      real, dimension(inp%nvar_btp-1,G%npoin), intent(out) :: rhs_btp

      real :: sc_x, sc_y, sc_z, Hq, qu, quv, qvu, qv, quz, qvz, qwu, qwv, qw, H_b
      real :: wq, hi, dhdx, dhdy, dhdz, tb_u, tb_v, tb_w, ope, ope2
      real :: dp, dpp, udp, vdp, wdp, ub, vb, wb, ubot, vbot, wbot, spd, pbq
      real :: pp_k, up_k, vp_k, wp_k
      real :: sum_up2, sum_uv, sum_vu, sum_vp2, sum_uw, sum_vw, sum_ww, sum_wu, sum_wv
      real :: cfx, cfy, cfz
      integer :: I, Iq, ip
      integer :: npts_l, botfr_l
      real    :: cd_l, alpha_bot_l
      logical :: has_w

      npts_l      = b%npts
      botfr_l     = inp%botfr
      cd_l        = real(inp%cd_mlswe)
      alpha_bot_l = real(init%alpha_mlswe(inp%nlayers))
      has_w       = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex

      !$acc data present(init, tsp, btp,                                             &
      !$acc               rhs_btp, qb_df,                                            &
      !$acc               init%grad_zbot_quad, init%tau_wind, tsp%psih,             &
      !$acc               tsp%dpsidx, tsp%dpsidy, tsp%dpsidz,                       &
      !$acc               tsp%dpsidz_x, tsp%dpsidz_y, tsp%dpsidz_z,                &
      !$acc               tsp%indexq, tsp%wjac,                                     &
      !$acc               init%coriolis_3d_quad,                                     &
      !$acc               btp%tau_bot_ave, btp%H_ave, btp%Qu_ave,                   &
      !$acc               btp%Qv_ave, btp%Qw_ave, btp%ope_ave,                     &
      !$acc               btp%uvb_ave, btp%btp_mass_flux_ave, btp%ope2_ave,        &
      !$acc               btp%bcl_H, btp%bcl_flux, btp%bcl_btp_flux, btp%pbq)

      !$acc kernels
      rhs_btp = 0.0
      !$acc end kernels

      ! Abort with element/rank diagnostics rather than SIGFPE-ing on ub=udp/dp
      ! below if the barotropic depth has gone negative or NaN (e.g. interface
      ! outcropping near steep topography feeding through into qb_df).
      call check_btp_thickness(b, G, inp, qb_df, 'create_rhs_btp_volume_qdf_qp')

      !$acc parallel loop gang                                                       &
      !$acc   private(dp, dpp, udp, vdp, wdp, pbq, H_b, Hq, wq, ub, vb, wb,       &
      !$acc           ubot, vbot, wbot, spd, tb_u, tb_v, tb_w,                    &
      !$acc           cfx, cfy, cfz, sc_x, sc_y, sc_z, ope, ope2,                 &
      !$acc           pp_k, up_k, vp_k, wp_k,                                      &
      !$acc           sum_up2, sum_uv, sum_vu, sum_vp2, sum_uw, sum_vw, sum_ww,  &
      !$acc           qu, quv, qvu, qv, quz, qvz, qwu, qwv, qw,                   &
      !$acc           hi, dhdx, dhdy, dhdz, I)                                     &
      !$acc   firstprivate(npts_l, botfr_l, cd_l, alpha_bot_l, has_w)
      do Iq = 1, G%npoin_q

         dp = 0.0;  dpp = 0.0;  udp = 0.0;  vdp = 0.0;  wdp = 0.0

         ! Barotropic interpolation to quad point.
         !$acc loop seq
         do ip = 1, npts_l
            I   = tsp%indexq(ip, Iq)
            hi  = tsp%psih(ip, Iq)
            dp  = dp  + hi * qb_df(1, I)
            dpp = dpp + hi * qb_df(2, I)
            udp = udp + hi * qb_df(3, I)
            vdp = vdp + hi * qb_df(4, I)
            if (has_w) wdp = wdp + hi * qb_df(5, I)
         end do

         ! Load BCL layer integrals precomputed once per RK stage.
         H_b     = btp%bcl_H(Iq)
         sum_up2 = btp%bcl_flux(1,1,Iq)
         sum_uv  = btp%bcl_flux(1,2,Iq)
         sum_vu  = btp%bcl_flux(2,1,Iq)
         sum_vp2 = btp%bcl_flux(2,2,Iq)
         pp_k    = btp%bcl_btp_flux(1,Iq)
         up_k    = btp%bcl_btp_flux(2,Iq)
         vp_k    = btp%bcl_btp_flux(3,Iq)
         pbq     = btp%pbq(Iq)

         sum_uw = 0.0;  sum_vw = 0.0;  sum_ww = 0.0;  wp_k = 0.0
         sum_wu = 0.0;  sum_wv = 0.0
         if (has_w) then
            sum_uw = btp%bcl_flux(3,1,Iq)
            sum_vw = btp%bcl_flux(3,2,Iq)
            sum_wu = btp%bcl_flux(1,3,Iq)
            sum_wv = btp%bcl_flux(2,3,Iq)
            sum_ww = btp%bcl_flux(3,3,Iq)
            wp_k   = btp%bcl_btp_flux(4,Iq)
         end if

         wq = tsp%wjac(Iq)
         ub = udp / dp;  vb = vdp / dp;  wb = wdp / dp

         ! 3D Coriolis vector (f * r̂); reduces to the f-plane pair (0,0,f) on
         ! cartesian since kvector=(0,0,1) there — see wind_stress_coriolis.
         cfx = init%coriolis_3d_quad(1, Iq)
         cfy = init%coriolis_3d_quad(2, Iq)
         cfz = init%coriolis_3d_quad(3, Iq)

         tb_u = 0.0;  tb_v = 0.0;  tb_w = 0.0
         if (botfr_l == 1) then
            ubot = up_k + ub;  vbot = vp_k + vb;  wbot = wp_k + wb
            spd  = (cd_l / gravity) * pp_k
            tb_u = spd * ubot;  tb_v = spd * vbot;  tb_w = spd * wbot
         elseif (botfr_l == 2) then
            ubot = up_k + ub;  vbot = vp_k + vb;  wbot = wp_k + wb
            spd  = (cd_l / alpha_bot_l) * sqrt(ubot**2 + vbot**2 + wbot**2)
            tb_u = spd * ubot;  tb_v = spd * vbot;  tb_w = spd * wbot
         end if

         ! Source terms (3D Coriolis cross-product -f×p*u; reduces exactly to
         ! the old f-plane form on cartesian since cfx=cfy=0 and wdp=0 there).
         sc_x = -(cfy*wdp - cfz*vdp)                                  &
            + gravity * (init%tau_wind(1, Iq) - tb_u)                 &
            - gravity * dp * init%grad_zbot_quad(1, Iq)
         sc_y = -(cfz*udp - cfx*wdp)                                  &
            + gravity * (init%tau_wind(2, Iq) - tb_v)                 &
            - gravity * dp * init%grad_zbot_quad(2, Iq)

         ope  = 1.0 + (dpp / pbq)
         ope2 = ope * ope
         Hq   = ope2 * H_b

         qu  = ub * udp + ope * sum_up2
         quv = vb * udp + ope * sum_uv
         qvu = ub * vdp + ope * sum_vu
         qv  = vb * vdp + ope * sum_vp2

         ! Time-average accumulators (each Iq is unique per gang — no race).
         btp%H_ave(Iq)    = btp%H_ave(Iq)    + Hq
         btp%Qu_ave(1,Iq) = btp%Qu_ave(1,Iq) + qu
         btp%Qu_ave(2,Iq) = btp%Qu_ave(2,Iq) + quv
         btp%Qv_ave(1,Iq) = btp%Qv_ave(1,Iq) + qvu
         btp%Qv_ave(2,Iq) = btp%Qv_ave(2,Iq) + qv

         btp%tau_bot_ave(1,Iq)       = btp%tau_bot_ave(1,Iq)       + tb_u
         btp%tau_bot_ave(2,Iq)       = btp%tau_bot_ave(2,Iq)       + tb_v
         btp%ope_ave(Iq)             = btp%ope_ave(Iq)             + ope
         btp%ope2_ave(Iq)            = btp%ope2_ave(Iq)            + ope2
         btp%btp_mass_flux_ave(1,Iq) = btp%btp_mass_flux_ave(1,Iq) + udp
         btp%btp_mass_flux_ave(2,Iq) = btp%btp_mass_flux_ave(2,Iq) + vdp
         btp%uvb_ave(1,Iq)           = btp%uvb_ave(1,Iq)           + ub
         btp%uvb_ave(2,Iq)           = btp%uvb_ave(2,Iq)           + vb

         ! Sphere-only: z-momentum flux/source terms and their time averages.
         sc_z = 0.0;  quz = 0.0;  qvz = 0.0;  qwu = 0.0;  qwv = 0.0;  qw = 0.0
         if (has_w) then
            sc_z = -(cfx*vdp - cfy*udp)                                &
                 - gravity * tb_w                                        &
                 - gravity * dp * init%grad_zbot_quad(3, Iq)

            quz = wb * udp + ope * sum_uw
            qvz = wb * vdp + ope * sum_vw
            qwu = ub * wdp + ope * sum_wu
            qwv = vb * wdp + ope * sum_wv
            qw  = wb * wdp + ope * sum_ww

            btp%Qu_ave(3,Iq)            = btp%Qu_ave(3,Iq)            + quz
            btp%Qv_ave(3,Iq)            = btp%Qv_ave(3,Iq)            + qvz
            btp%Qw_ave(1,Iq)            = btp%Qw_ave(1,Iq)            + qwu
            btp%Qw_ave(2,Iq)            = btp%Qw_ave(2,Iq)            + qwv
            btp%Qw_ave(3,Iq)            = btp%Qw_ave(3,Iq)            + qw
            btp%tau_bot_ave(3,Iq)       = btp%tau_bot_ave(3,Iq)       + tb_w
            btp%btp_mass_flux_ave(3,Iq) = btp%btp_mass_flux_ave(3,Iq) + wdp
            btp%uvb_ave(3,Iq)           = btp%uvb_ave(3,Iq)           + wb
         end if

         ! Scatter: atomics handle cross-gang races on shared DOF nodes.
         !$acc loop seq
         do ip = 1, npts_l
            I    = tsp%indexq(ip, Iq)
            hi   = tsp%psih(ip, Iq)
            dhdx = tsp%dpsidx(ip, Iq)
            dhdy = tsp%dpsidy(ip, Iq)
            dhdz = 0.0
            if (has_w) then
               dhdx = dhdx + tsp%dpsidz_x(ip, Iq)
               dhdy = dhdy + tsp%dpsidz_y(ip, Iq)
               dhdz = tsp%dpsidz(ip, Iq) + tsp%dpsidz_z(ip, Iq)
            end if
            !$acc atomic update
            rhs_btp(1, I) = rhs_btp(1, I) + wq * (dhdx*udp + dhdy*vdp + dhdz*wdp)
            !$acc atomic update
            rhs_btp(2, I) = rhs_btp(2, I) + wq * (hi*sc_x + dhdx*(Hq+qu) + dhdy*quv + dhdz*quz)
            !$acc atomic update
            rhs_btp(3, I) = rhs_btp(3, I) + wq * (hi*sc_y + dhdx*qvu    + dhdy*(Hq+qv) + dhdz*qvz)
            if (has_w) then
               !$acc atomic update
               rhs_btp(4, I) = rhs_btp(4, I) + wq * (hi*sc_z + dhdx*qwu + dhdy*qwv + dhdz*(Hq+qw))
            end if
         end do

      end do
      !$acc end parallel loop
      !$acc end data

   end subroutine create_rhs_btp_volume_qdf_qp

   subroutine create_btp_fluxes_qdf_gpu(G, b, inp, mf, init, btp, rhs_btp, qb, qprime_df)
      !===========================================================================
      !  Face flux contribution to the barotropic RHS (Riemann solver).
      !
      !  GPU parallelism, two kernels:
      !    1. one thread per (face, quadrature point): Riemann flux, stored
      !       weighted in btp%flux_face_q; the inner loops (basis nodes n,
      !       layers k) are sequential within the thread.
      !    2. one thread per (face, face node): sums the quad-point fluxes
      !       against psiq and scatters them to rhs_btp.  Scattering from
      !       kernel 1 made the nq threads of a face atomically update the
      !       same ngl nodes at once.
      !
      !  Race conditions:
      !    rhs_btp(m,I): corner nodes are shared by two faces of an
      !      element — !$acc atomic update in kernel 2.
      !    btp%*_face_ave(iquad,iface) and flux_face_q: each (iquad,iface)
      !      pair is touched by exactly one thread — no atomics needed.
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

      real, dimension(inp%nvar_btp-1,G%npoin),             intent(inout) :: rhs_btp
      real, dimension(inp%nvar_btp,G%npoin),               intent(in)    :: qb
      real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers),   intent(in)    :: qprime_df

      integer :: iface, iquad, ifq, ifn, el, er, I, n, k
      real    :: wq, hi, nxl, nyl, nzl, un
      real    :: pbl, pbr, clam, one_eta, one_eta2, half_clam, quarter_clam
      real    :: pU_L, pU_R, half_clam_dqb2, pbpert_edge
      real    :: qbl1, qbl2, qbl3, qbl4, qbl5, qbr1, qbr2, qbr3, qbr4, qbr5
      real    :: c_minus, c_plus
      real    :: ul, ur, vl, vr, wl, wr, opl, opr
      real    :: H_bcl_ql, H_bcl_qr, H_bcl_q, oe2_Hql, oe2_Hqr
      real    :: Qu_ql1, Qu_ql2, Qu_ql3, Qu_qr1, Qu_qr2, Qu_qr3
      real    :: Qv_ql1, Qv_ql2, Qv_ql3, Qv_qr1, Qv_qr2, Qv_qr3
      real    :: Qw_ql1, Qw_ql2, Qw_ql3, Qw_qr1, Qw_qr2, Qw_qr3
      real    :: flux_edge_x, flux_edge_y, flux_edge_z, flux_pb, flux_u, flux_v, flux_w
      real    :: pprime_lk, pprime_rk, pprime_lk1, pprime_rk1
      real    :: pkl, pkr, ukl, ukr, vkl, vkr, wkl, wkr
      real    :: ope_ppl_k, ope_ppr_k, uv_cross_l, uv_cross_r
      real    :: wq_flux_pb, wq_flux_u, wq_flux_v, wq_flux_w
      real    :: s_pb, s_u, s_v, s_w
      logical :: has_w

      ! Capture derived-type scalars before the parallel region so they are
      ! available as firstprivate values inside the GPU kernel.
      integer :: ngl_f, nq_f, nlayers_f, nface_f
      real    :: dry_cutoff_f
      ngl_f     = b%ngl
      nq_f      = b%nq
      nlayers_f = inp%nlayers
      nface_f   = G%nface
      dry_cutoff_f = inp%dry_cutoff   ! scalar copy: inp itself is not on the device
      has_w     = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex

      !$acc data present(btp, rhs_btp, qb, qprime_df,                              &
      !$acc               G%face, G%face_type, G%intma,                             &
      !$acc               mf%imapl, mf%imapr, mf%normal_vector_q, mf%jac_faceq,   &
      !$acc               b%psiq,                                                   &
      !$acc               init%pbprime_df, init%alpha_mlswe,                        &
      !$acc               btp%H_face_ave, btp%ope_face_ave, btp%ope2_face_ave,     &
      !$acc               btp%btp_mass_flux_face_ave,                               &
      !$acc               btp%Qu_face_ave, btp%Qv_face_ave, btp%Qw_face_ave,       &
      !$acc               btp%one_plus_eta_edge_2_ave,                              &
      !$acc               btp%uvb_face_ave, btp%flux_face_q, mf%face_I)

      !$acc parallel loop async(1) gang vector                                               &
      !$acc   private(iface, iquad, el, er, I, n, k,                              &
      !$acc           qbl1, qbl2, qbl3, qbl4, qbl5, qbr1, qbr2, qbr3, qbr4, qbr5,  &
      !$acc           pbl, pbr,                                                    &
      !$acc           nxl, nyl, nzl, wq, hi, un,                                   &
      !$acc           clam, one_eta, one_eta2, half_clam, quarter_clam,            &
      !$acc           pU_L, pU_R, half_clam_dqb2, pbpert_edge,                    &
      !$acc           c_minus, c_plus, ul, ur, vl, vr, wl, wr, opl, opr,          &
      !$acc           H_bcl_ql, H_bcl_qr, H_bcl_q, oe2_Hql, oe2_Hqr,             &
      !$acc           Qu_ql1, Qu_ql2, Qu_ql3, Qu_qr1, Qu_qr2, Qu_qr3,            &
      !$acc           Qv_ql1, Qv_ql2, Qv_ql3, Qv_qr1, Qv_qr2, Qv_qr3,            &
      !$acc           Qw_ql1, Qw_ql2, Qw_ql3, Qw_qr1, Qw_qr2, Qw_qr3,            &
      !$acc           flux_edge_x, flux_edge_y, flux_edge_z, flux_pb, flux_u, flux_v, flux_w, &
      !$acc           pprime_lk, pprime_rk, pprime_lk1, pprime_rk1,               &
      !$acc           pkl, pkr, ukl, ukr, vkl, vkr, wkl, wkr,                     &
      !$acc           ope_ppl_k, ope_ppr_k, uv_cross_l, uv_cross_r,               &
      !$acc           wq_flux_pb, wq_flux_u, wq_flux_v, wq_flux_w)              &
      !$acc   firstprivate(has_w)
      do ifq = 1, nface_f*nq_f

         iface = (ifq - 1)/nq_f + 1
         iquad = ifq - (iface - 1)*nq_f

         if (G%face_type(iface) == 2) cycle

         el = G%face(7, iface)
         er = G%face(8, iface)

            nxl = mf%normal_vector_q(1, iquad, 1, iface)
            nyl = mf%normal_vector_q(2, iquad, 1, iface)
            nzl = mf%normal_vector_q(3, iquad, 1, iface)

            ! Project barotropic LEFT state
            qbl1 = 0.0;  qbl2 = 0.0;  qbl3 = 0.0;  qbl4 = 0.0;  qbl5 = 0.0
            pbl = 0.0
            !$acc loop seq
            do n = 1, ngl_f
               hi = b%psiq(n, iquad)
               I  = mf%face_I(n,1,iface)
               qbl1 = qbl1 + hi * qb(1,I)
               qbl2 = qbl2 + hi * qb(2,I)
               qbl3 = qbl3 + hi * qb(3,I)
               qbl4 = qbl4 + hi * qb(4,I)
               if (has_w) qbl5 = qbl5 + hi * qb(5,I)
               pbl    = pbl    + hi * init%pbprime_df(I)
            end do

            ! Project barotropic RIGHT state or apply boundary condition.
            if (er > 0) then
               qbr1 = 0.0;  qbr2 = 0.0;  qbr3 = 0.0;  qbr4 = 0.0;  qbr5 = 0.0
               pbr = 0.0
               !$acc loop seq
               do n = 1, ngl_f
                  hi = b%psiq(n, iquad)
                  I  = mf%face_I(n,2,iface)
                  qbr1 = qbr1 + hi * qb(1,I)
                  qbr2 = qbr2 + hi * qb(2,I)
                  qbr3 = qbr3 + hi * qb(3,I)
                  qbr4 = qbr4 + hi * qb(4,I)
                  if (has_w) qbr5 = qbr5 + hi * qb(5,I)
                  pbr    = pbr    + hi * init%pbprime_df(I)
               end do
            else
               qbr1 = qbl1;  qbr2 = qbl2
               qbr3 = qbl3;  qbr4 = qbl4;  qbr5 = qbl5
               pbr = pbl
               if (er == -4) then          ! slip wall
                  un     = nxl*qbl3 + nyl*qbl4 + nzl*qbl5
                  qbr3 = qbl3 - 2.0*un*nxl
                  qbr4 = qbl4 - 2.0*un*nyl
                  if (has_w) qbr5 = qbl5 - 2.0*un*nzl
               else if (er == -2) then     ! no-slip wall
                  qbr3 = -qbl3;  qbr4 = -qbl4
                  if (has_w) qbr5 = -qbl5
               end if
            end if

            ! Wave speeds
            pU_L = nxl*qbl3 + nyl*qbl4 + nzl*qbl5
            pU_R = -(nxl*qbr3 + nyl*qbr4 + nzl*qbr5)

            c_minus      = sqrt(init%alpha_mlswe(nlayers_f) * pbr)
            c_plus       = sqrt(init%alpha_mlswe(nlayers_f) * pbl)
            clam         = max(c_minus, c_plus)
            half_clam    = 0.5  * clam
            quarter_clam = 0.25 * clam

            pbpert_edge = 0.5*(qbl2 + qbr2) + (0.5/clam)*(pU_L + pU_R)
            one_eta     = 1.0 + pbpert_edge / pbl
            one_eta2    = one_eta * one_eta

            ul = qbl3 / qbl1;  ur = qbr3 / qbr1
            vl = qbl4 / qbl1;  vr = qbr4 / qbr1
            wl = 0.0;  wr = 0.0
            if (has_w) then
               wl = qbl5 / qbl1;  wr = qbr5 / qbr1
            end if

            ! Initialise flux tensors from barotropic state
            Qu_ql1 = ul * qbl3;  Qu_ql2 = vl * qbl3;  Qu_ql3 = 0.0
            Qu_qr1 = ur * qbr3;  Qu_qr2 = vr * qbr3;  Qu_qr3 = 0.0
            Qv_ql1 = ul * qbl4;  Qv_ql2 = vl * qbl4;  Qv_ql3 = 0.0
            Qv_qr1 = ur * qbr4;  Qv_qr2 = vr * qbr4;  Qv_qr3 = 0.0
            Qw_ql1 = 0.0;  Qw_ql2 = 0.0;  Qw_ql3 = 0.0
            Qw_qr1 = 0.0;  Qw_qr2 = 0.0;  Qw_qr3 = 0.0
            if (has_w) then
               Qu_ql3 = wl * qbl3;  Qu_qr3 = wr * qbr3
               Qv_ql3 = wl * qbl4;  Qv_qr3 = wr * qbr4
               Qw_ql1 = ul * qbl5;  Qw_ql2 = vl * qbl5;  Qw_ql3 = wl * qbl5
               Qw_qr1 = ur * qbr5;  Qw_qr2 = vr * qbr5;  Qw_qr3 = wr * qbr5
            end if
            H_bcl_ql  = 0.0;  H_bcl_qr  = 0.0
            pprime_lk = 0.0;  pprime_rk = 0.0

            ! Layer projection + flux accumulation
            !$acc loop seq
            do k = 1, nlayers_f

               pkl = 0.0;  ukl = 0.0;  vkl = 0.0;  wkl = 0.0
               !$acc loop seq
               do n = 1, ngl_f
                  hi  = b%psiq(n, iquad)
                  I   = mf%face_I(n,1,iface)
                  pkl = pkl + hi * qprime_df(1, I, k)
                  ukl = ukl + hi * qprime_df(2, I, k)
                  vkl = vkl + hi * qprime_df(3, I, k)
                  if (has_w) wkl = wkl + hi * qprime_df(4, I, k)
               end do

               if (er > 0) then
                  pkr = 0.0;  ukr = 0.0;  vkr = 0.0;  wkr = 0.0
                  !$acc loop seq
                  do n = 1, ngl_f
                     hi  = b%psiq(n, iquad)
                     I   = mf%face_I(n,2,iface)
                     pkr = pkr + hi * qprime_df(1, I, k)
                     ukr = ukr + hi * qprime_df(2, I, k)
                     vkr = vkr + hi * qprime_df(3, I, k)
                     if (has_w) wkr = wkr + hi * qprime_df(4, I, k)
                  end do
               else
                  pkr = pkl
                  if (er == -4) then
                     un  = nxl * ukl + nyl * vkl + nzl * wkl
                     ukr = ukl - 2.0 * un * nxl
                     vkr = vkl - 2.0 * un * nyl
                     if (has_w) wkr = wkl - 2.0 * un * nzl
                  else if (er == -2) then
                     ukr = -ukl;  vkr = -vkl
                     if (has_w) wkr = -wkl
                  else
                     ukr = ukl;  vkr = vkl;  wkr = wkl
                  end if
               end if

               ! Dry-cell protection, matching the volume/face BCL flux
               ! routines. Thickness is CLAMPED (not zeroed) to a floor so
               ! H_bcl_ql/qr's cumulative pprime bookkeeping stays continuous
               ! across dry_cutoff; velocity is zeroed, which alone fully
               ! suppresses this layer's Qu/Qv/Qw flux terms (each is a
               ! product of two velocity components) regardless of the
               ! floored thickness. Without this, a phantom contribution
               ! would leak into btp%H_face_ave/Qu_face_ave/Qv_face_ave/
               ! Qw_face_ave -- exactly the consistency targets the volume/
               ! face routines' deficit-redistribution logic spreads onto
               ! the wet layers.
               if (pkl < (gravity/init%alpha_mlswe(k)) * dry_cutoff_f) then
                  pkl = (gravity/init%alpha_mlswe(k)) * dry_cutoff_f
                  ukl = 0.0;  vkl = 0.0;  wkl = 0.0
               end if
               if (pkr < (gravity/init%alpha_mlswe(k)) * dry_cutoff_f) then
                  pkr = (gravity/init%alpha_mlswe(k)) * dry_cutoff_f
                  ukr = 0.0;  vkr = 0.0;  wkr = 0.0
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

               if (has_w) then
                  Qu_ql3 = Qu_ql3 + ukl * wkl * ope_ppl_k;  Qu_qr3 = Qu_qr3 + ukr * wkr * ope_ppr_k
                  Qv_ql3 = Qv_ql3 + vkl * wkl * ope_ppl_k;  Qv_qr3 = Qv_qr3 + vkr * wkr * ope_ppr_k
                  Qw_ql1 = Qw_ql1 + ukl * wkl * ope_ppl_k;  Qw_qr1 = Qw_qr1 + ukr * wkr * ope_ppr_k
                  Qw_ql2 = Qw_ql2 + vkl * wkl * ope_ppl_k;  Qw_qr2 = Qw_qr2 + vkr * wkr * ope_ppr_k
                  Qw_ql3 = Qw_ql3 + wkl * wkl * ope_ppl_k;  Qw_qr3 = Qw_qr3 + wkr * wkr * ope_ppr_k
               end if

               pprime_lk1 = pprime_lk + pkl
               pprime_rk1 = pprime_rk + pkr
               H_bcl_ql   = H_bcl_ql  + 0.5*init%alpha_mlswe(k) * (pprime_lk1 + pprime_lk) * pkl
               H_bcl_qr   = H_bcl_qr  + 0.5*init%alpha_mlswe(k) * (pprime_rk1 + pprime_rk) * pkr
               pprime_lk  = pprime_lk1
               pprime_rk  = pprime_rk1

            end do  ! k

            H_bcl_q = 0.5 * (H_bcl_ql + H_bcl_qr)

            half_clam_dqb2 = half_clam * (qbl2 - qbr2)
            flux_edge_x = 0.5*(qbl3 + qbr3) + nxl * half_clam_dqb2
            flux_edge_y = 0.5*(qbl4 + qbr4) + nyl * half_clam_dqb2
            flux_edge_z = 0.0
            if (has_w) flux_edge_z = 0.5*(qbl5 + qbr5) + nzl * half_clam_dqb2
            flux_pb     = nxl*flux_edge_x + nyl*flux_edge_y + nzl*flux_edge_z

            H_bcl_q = one_eta2 * H_bcl_q

            ! Time-average accumulators — (iquad,iface) unique per gang×seq step.
            btp%btp_mass_flux_face_ave(1,iquad,iface) = btp%btp_mass_flux_face_ave(1,iquad,iface) + flux_edge_x
            btp%btp_mass_flux_face_ave(2,iquad,iface) = btp%btp_mass_flux_face_ave(2,iquad,iface) + flux_edge_y
            btp%H_face_ave(iquad,iface)    = btp%H_face_ave(iquad,iface)    + H_bcl_q
            btp%Qu_face_ave(1,iquad,iface) = btp%Qu_face_ave(1,iquad,iface) + 0.5*(Qu_ql1 + Qu_qr1)
            btp%Qu_face_ave(2,iquad,iface) = btp%Qu_face_ave(2,iquad,iface) + 0.5*(Qu_ql2 + Qu_qr2)
            btp%Qv_face_ave(1,iquad,iface) = btp%Qv_face_ave(1,iquad,iface) + 0.5*(Qv_ql1 + Qv_qr1)
            btp%Qv_face_ave(2,iquad,iface) = btp%Qv_face_ave(2,iquad,iface) + 0.5*(Qv_ql2 + Qv_qr2)

            opl = 1.0 + qbl2/pbl
            opr = 1.0 + qbr2/pbr
            btp%ope_face_ave(1,iquad,iface)  = btp%ope_face_ave(1,iquad,iface)  + opl
            btp%ope_face_ave(2,iquad,iface)  = btp%ope_face_ave(2,iquad,iface)  + opr
            btp%ope2_face_ave(1,iquad,iface) = btp%ope2_face_ave(1,iquad,iface) + opl*opl
            btp%ope2_face_ave(2,iquad,iface) = btp%ope2_face_ave(2,iquad,iface) + opr*opr
            btp%one_plus_eta_edge_2_ave(iquad,iface) = btp%one_plus_eta_edge_2_ave(iquad,iface) + one_eta2
            btp%uvb_face_ave(1,1,iquad,iface) = btp%uvb_face_ave(1,1,iquad,iface) + ul
            btp%uvb_face_ave(1,2,iquad,iface) = btp%uvb_face_ave(1,2,iquad,iface) + ur
            btp%uvb_face_ave(2,1,iquad,iface) = btp%uvb_face_ave(2,1,iquad,iface) + vl
            btp%uvb_face_ave(2,2,iquad,iface) = btp%uvb_face_ave(2,2,iquad,iface) + vr

            if (has_w) then
               btp%btp_mass_flux_face_ave(3,iquad,iface) = btp%btp_mass_flux_face_ave(3,iquad,iface) + flux_edge_z
               btp%Qu_face_ave(3,iquad,iface) = btp%Qu_face_ave(3,iquad,iface) + 0.5*(Qu_ql3 + Qu_qr3)
               btp%Qv_face_ave(3,iquad,iface) = btp%Qv_face_ave(3,iquad,iface) + 0.5*(Qv_ql3 + Qv_qr3)
               btp%Qw_face_ave(1,iquad,iface) = btp%Qw_face_ave(1,iquad,iface) + 0.5*(Qw_ql1 + Qw_qr1)
               btp%Qw_face_ave(2,iquad,iface) = btp%Qw_face_ave(2,iquad,iface) + 0.5*(Qw_ql2 + Qw_qr2)
               btp%Qw_face_ave(3,iquad,iface) = btp%Qw_face_ave(3,iquad,iface) + 0.5*(Qw_ql3 + Qw_qr3)
               btp%uvb_face_ave(3,1,iquad,iface) = btp%uvb_face_ave(3,1,iquad,iface) + wl
               btp%uvb_face_ave(3,2,iquad,iface) = btp%uvb_face_ave(3,2,iquad,iface) + wr
            end if

            ! Momentum fluxes
            oe2_Hql = one_eta2 * H_bcl_ql
            oe2_Hqr = one_eta2 * H_bcl_qr

            Qu_ql1 = Qu_ql1 + oe2_Hql
            Qu_qr1 = Qu_qr1 + oe2_Hqr
            Qv_ql2 = Qv_ql2 + oe2_Hql
            Qv_qr2 = Qv_qr2 + oe2_Hqr
            if (has_w) then
               Qw_ql3 = Qw_ql3 + oe2_Hql
               Qw_qr3 = Qw_qr3 + oe2_Hqr
            end if

            flux_u = 0.5 * (nxl*(Qu_ql1 + Qu_qr1) + nyl*(Qu_ql2 + Qu_qr2) + nzl*(Qu_ql3 + Qu_qr3)) &
               - quarter_clam * (qbr3 - qbl3)
            flux_v = 0.5 * (nxl*(Qv_ql1 + Qv_qr1) + nyl*(Qv_ql2 + Qv_qr2) + nzl*(Qv_ql3 + Qv_qr3)) &
               - quarter_clam * (qbr4 - qbl4)
            flux_w = 0.0
            if (has_w) flux_w = 0.5 * (nxl*(Qw_ql1 + Qw_qr1) + nyl*(Qw_ql2 + Qw_qr2) + nzl*(Qw_ql3 + Qw_qr3)) &
               - quarter_clam * (qbr5 - qbl5)

            wq         = mf%jac_faceq(iquad, 1, iface)
            wq_flux_pb = wq * flux_pb
            wq_flux_u  = wq * flux_u
            wq_flux_v  = wq * flux_v
            wq_flux_w  = wq * flux_w

            btp%flux_face_q(ifq,1) = wq_flux_pb
            btp%flux_face_q(ifq,2) = wq_flux_u
            btp%flux_face_q(ifq,3) = wq_flux_v
            btp%flux_face_q(ifq,4) = wq_flux_w

      end do  ! ifq
      !$acc end parallel loop

      ! RHS scatter: one thread per (face, face node).  The left and right
      ! nodes get the same sum with opposite signs.
      !$acc parallel loop async(1) gang vector                                               &
      !$acc   private(iface, n, el, er, I, iquad, ifq, hi, s_pb, s_u, s_v, s_w)    &
      !$acc   firstprivate(has_w)
      do ifn = 1, nface_f*ngl_f

         iface = (ifn - 1)/ngl_f + 1
         n     = ifn - (iface - 1)*ngl_f

         if (G%face_type(iface) == 2) cycle

         el = G%face(7, iface)
         er = G%face(8, iface)

         s_pb = 0.0;  s_u = 0.0;  s_v = 0.0;  s_w = 0.0
         !$acc loop seq
         do iquad = 1, nq_f
            ifq  = iquad + (iface - 1)*nq_f
            hi   = b%psiq(n, iquad)
            s_pb = s_pb + hi * btp%flux_face_q(ifq,1)
            s_u  = s_u  + hi * btp%flux_face_q(ifq,2)
            s_v  = s_v  + hi * btp%flux_face_q(ifq,3)
            if (has_w) s_w = s_w + hi * btp%flux_face_q(ifq,4)
         end do

         I = mf%face_I(n,1,iface)
         !$acc atomic update
         rhs_btp(1,I) = rhs_btp(1,I) - s_pb
         !$acc atomic update
         rhs_btp(2,I) = rhs_btp(2,I) - s_u
         !$acc atomic update
         rhs_btp(3,I) = rhs_btp(3,I) - s_v
         if (has_w) then
            !$acc atomic update
            rhs_btp(4,I) = rhs_btp(4,I) - s_w
         end if

         if (er > 0) then
            I = mf%face_I(n,2,iface)
            !$acc atomic update
            rhs_btp(1,I) = rhs_btp(1,I) + s_pb
            !$acc atomic update
            rhs_btp(2,I) = rhs_btp(2,I) + s_u
            !$acc atomic update
            rhs_btp(3,I) = rhs_btp(3,I) + s_v
            if (has_w) then
               !$acc atomic update
               rhs_btp(4,I) = rhs_btp(4,I) + s_w
            end if
         end if

      end do  ! ifn
      !$acc end parallel loop
      !$acc end data

   end subroutine create_btp_fluxes_qdf_gpu

   subroutine create_btp_fluxes_qdf(G, b, inp, mf, init, btp, rhs_btp, qb, qprime_df)
      !===========================================================================
      !  CPU version of create_btp_fluxes_qdf_gpu (used when built without
      !  OpenACC). Same physics, one face at a time with its node indices cached
      !  once per face; the GPU version instead runs one thread per (face,
      !  quadrature point) and looks the indices up where they are used.
      !
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

      real, dimension(inp%nvar_btp-1,G%npoin),             intent(inout) :: rhs_btp
      real, dimension(inp%nvar_btp,G%npoin),               intent(in)    :: qb
      real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers),   intent(in)    :: qprime_df

      integer :: iface, iquad, el, er, il, jl, kl, ir, jr, kr, I, n, k
      real    :: wq, hi, nxl, nyl, nzl, un
      real    :: pbl, pbr, clam, one_eta, one_eta2, half_clam, quarter_clam
      real    :: pU_L, pU_R, half_clam_dqb2, pbpert_edge
      real    :: qbl(5), qbr(5)
      real    :: c_minus, c_plus
      real    :: ul, ur, vl, vr, wl, wr, opl, opr
      real    :: H_bcl_ql, H_bcl_qr, H_bcl_q, oe2_Hql, oe2_Hqr
      real    :: Qu_ql1, Qu_ql2, Qu_ql3, Qu_qr1, Qu_qr2, Qu_qr3
      real    :: Qv_ql1, Qv_ql2, Qv_ql3, Qv_qr1, Qv_qr2, Qv_qr3
      real    :: Qw_ql1, Qw_ql2, Qw_ql3, Qw_qr1, Qw_qr2, Qw_qr3
      real    :: flux_edge_x, flux_edge_y, flux_edge_z, flux_pb, flux_u, flux_v, flux_w
      real    :: pprime_lk, pprime_rk, pprime_lk1, pprime_rk1
      real    :: pkl, pkr, ukl, ukr, vkl, vkr, wkl, wkr
      real    :: ope_ppl_k, ope_ppr_k, uv_cross_l, uv_cross_r
      real    :: wq_flux_pb, wq_flux_u, wq_flux_v, wq_flux_w
      real,    dimension(b%ngl) :: hi_c
      integer, dimension(b%ngl) :: I_l, I_r
      logical :: has_w

      ! Capture derived-type scalars before the parallel region so they are
      ! available as firstprivate values inside the GPU kernel.
      integer :: ngl_f, nq_f, nlayers_f, nface_f
      ngl_f     = b%ngl
      nq_f      = b%nq
      nlayers_f = inp%nlayers
      nface_f   = G%nface
      has_w     = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex

      !$acc data present(btp, rhs_btp, qb, qprime_df,                              &
      !$acc               G%face, G%face_type, G%intma,                             &
      !$acc               mf%imapl, mf%imapr, mf%normal_vector_q, mf%jac_faceq,   &
      !$acc               b%psiq,                                                   &
      !$acc               init%pbprime_df, init%alpha_mlswe,                        &
      !$acc               btp%H_face_ave, btp%ope_face_ave, btp%ope2_face_ave,     &
      !$acc               btp%btp_mass_flux_face_ave,                               &
      !$acc               btp%Qu_face_ave, btp%Qv_face_ave, btp%Qw_face_ave,       &
      !$acc               btp%one_plus_eta_edge_2_ave,                              &
      !$acc               btp%uvb_face_ave)

      !$acc parallel loop gang                                                      &
      !$acc   private(el, er, il, jl, kl, ir, jr, kr, I, n, k,                    &
      !$acc           I_l, I_r, hi_c,                                              &
      !$acc           qbl, qbr, pbl, pbr,                                          &
      !$acc           nxl, nyl, nzl, wq, hi, un,                                   &
      !$acc           clam, one_eta, one_eta2, half_clam, quarter_clam,            &
      !$acc           pU_L, pU_R, half_clam_dqb2, pbpert_edge,                    &
      !$acc           c_minus, c_plus, ul, ur, vl, vr, wl, wr, opl, opr,          &
      !$acc           H_bcl_ql, H_bcl_qr, H_bcl_q, oe2_Hql, oe2_Hqr,             &
      !$acc           Qu_ql1, Qu_ql2, Qu_ql3, Qu_qr1, Qu_qr2, Qu_qr3,            &
      !$acc           Qv_ql1, Qv_ql2, Qv_ql3, Qv_qr1, Qv_qr2, Qv_qr3,            &
      !$acc           Qw_ql1, Qw_ql2, Qw_ql3, Qw_qr1, Qw_qr2, Qw_qr3,            &
      !$acc           flux_edge_x, flux_edge_y, flux_edge_z, flux_pb, flux_u, flux_v, flux_w, &
      !$acc           pprime_lk, pprime_rk, pprime_lk1, pprime_rk1,               &
      !$acc           pkl, pkr, ukl, ukr, vkl, vkr, wkl, wkr,                     &
      !$acc           ope_ppl_k, ope_ppr_k, uv_cross_l, uv_cross_r,               &
      !$acc           wq_flux_pb, wq_flux_u, wq_flux_v, wq_flux_w)              &
      !$acc   firstprivate(has_w)
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
            nzl = mf%normal_vector_q(3, iquad, 1, iface)

            !$acc loop seq
            do n = 1, ngl_f
               hi_c(n) = b%psiq(n, iquad)
            end do

            ! Project barotropic LEFT state
            qbl(1) = 0.0;  qbl(2) = 0.0;  qbl(3) = 0.0;  qbl(4) = 0.0;  qbl(5) = 0.0
            pbl = 0.0
            !$acc loop seq
            do n = 1, ngl_f
               hi = hi_c(n)
               I  = I_l(n)
               qbl(1) = qbl(1) + hi * qb(1,I)
               qbl(2) = qbl(2) + hi * qb(2,I)
               qbl(3) = qbl(3) + hi * qb(3,I)
               qbl(4) = qbl(4) + hi * qb(4,I)
               if (has_w) qbl(5) = qbl(5) + hi * qb(5,I)
               pbl    = pbl    + hi * init%pbprime_df(I)
            end do

            ! Project barotropic RIGHT state or apply boundary condition.
            if (er > 0) then
               qbr(1) = 0.0;  qbr(2) = 0.0;  qbr(3) = 0.0;  qbr(4) = 0.0;  qbr(5) = 0.0
               pbr = 0.0
               !$acc loop seq
               do n = 1, ngl_f
                  hi = hi_c(n)
                  I  = I_r(n)
                  qbr(1) = qbr(1) + hi * qb(1,I)
                  qbr(2) = qbr(2) + hi * qb(2,I)
                  qbr(3) = qbr(3) + hi * qb(3,I)
                  qbr(4) = qbr(4) + hi * qb(4,I)
                  if (has_w) qbr(5) = qbr(5) + hi * qb(5,I)
                  pbr    = pbr    + hi * init%pbprime_df(I)
               end do
            else
               qbr(1) = qbl(1);  qbr(2) = qbl(2)
               qbr(3) = qbl(3);  qbr(4) = qbl(4);  qbr(5) = qbl(5)
               pbr = pbl
               if (er == -4) then          ! slip wall
                  un     = nxl*qbl(3) + nyl*qbl(4) + nzl*qbl(5)
                  qbr(3) = qbl(3) - 2.0*un*nxl
                  qbr(4) = qbl(4) - 2.0*un*nyl
                  if (has_w) qbr(5) = qbl(5) - 2.0*un*nzl
               else if (er == -2) then     ! no-slip wall
                  qbr(3) = -qbl(3);  qbr(4) = -qbl(4)
                  if (has_w) qbr(5) = -qbl(5)
               end if
            end if

            ! Wave speeds
            pU_L = nxl*qbl(3) + nyl*qbl(4) + nzl*qbl(5)
            pU_R = -(nxl*qbr(3) + nyl*qbr(4) + nzl*qbr(5))

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
            wl = 0.0;  wr = 0.0
            if (has_w) then
               wl = qbl(5) / qbl(1);  wr = qbr(5) / qbr(1)
            end if

            ! Initialise flux tensors from barotropic state
            Qu_ql1 = ul * qbl(3);  Qu_ql2 = vl * qbl(3);  Qu_ql3 = 0.0
            Qu_qr1 = ur * qbr(3);  Qu_qr2 = vr * qbr(3);  Qu_qr3 = 0.0
            Qv_ql1 = ul * qbl(4);  Qv_ql2 = vl * qbl(4);  Qv_ql3 = 0.0
            Qv_qr1 = ur * qbr(4);  Qv_qr2 = vr * qbr(4);  Qv_qr3 = 0.0
            Qw_ql1 = 0.0;  Qw_ql2 = 0.0;  Qw_ql3 = 0.0
            Qw_qr1 = 0.0;  Qw_qr2 = 0.0;  Qw_qr3 = 0.0
            if (has_w) then
               Qu_ql3 = wl * qbl(3);  Qu_qr3 = wr * qbr(3)
               Qv_ql3 = wl * qbl(4);  Qv_qr3 = wr * qbr(4)
               Qw_ql1 = ul * qbl(5);  Qw_ql2 = vl * qbl(5);  Qw_ql3 = wl * qbl(5)
               Qw_qr1 = ur * qbr(5);  Qw_qr2 = vr * qbr(5);  Qw_qr3 = wr * qbr(5)
            end if
            H_bcl_ql  = 0.0;  H_bcl_qr  = 0.0
            pprime_lk = 0.0;  pprime_rk = 0.0

            ! Layer projection + flux accumulation
            !$acc loop seq
            do k = 1, nlayers_f

               pkl = 0.0;  ukl = 0.0;  vkl = 0.0;  wkl = 0.0
               !$acc loop seq
               do n = 1, ngl_f
                  hi  = hi_c(n)
                  I   = I_l(n)
                  pkl = pkl + hi * qprime_df(1, I, k)
                  ukl = ukl + hi * qprime_df(2, I, k)
                  vkl = vkl + hi * qprime_df(3, I, k)
                  if (has_w) wkl = wkl + hi * qprime_df(4, I, k)
               end do

               if (er > 0) then
                  pkr = 0.0;  ukr = 0.0;  vkr = 0.0;  wkr = 0.0
                  !$acc loop seq
                  do n = 1, ngl_f
                     hi  = hi_c(n)
                     I   = I_r(n)
                     pkr = pkr + hi * qprime_df(1, I, k)
                     ukr = ukr + hi * qprime_df(2, I, k)
                     vkr = vkr + hi * qprime_df(3, I, k)
                     if (has_w) wkr = wkr + hi * qprime_df(4, I, k)
                  end do
               else
                  pkr = pkl
                  if (er == -4) then
                     un  = nxl * ukl + nyl * vkl + nzl * wkl
                     ukr = ukl - 2.0 * un * nxl
                     vkr = vkl - 2.0 * un * nyl
                     if (has_w) wkr = wkl - 2.0 * un * nzl
                  else if (er == -2) then
                     ukr = -ukl;  vkr = -vkl
                     if (has_w) wkr = -wkl
                  else
                     ukr = ukl;  vkr = vkl;  wkr = wkl
                  end if
               end if

               ! Dry-cell protection, matching the volume/face BCL flux
               ! routines. Thickness is CLAMPED (not zeroed) to a floor so
               ! H_bcl_ql/qr's cumulative pprime bookkeeping stays continuous
               ! across dry_cutoff; velocity is zeroed, which alone fully
               ! suppresses this layer's Qu/Qv/Qw flux terms (each is a
               ! product of two velocity components) regardless of the
               ! floored thickness. Without this, a phantom contribution
               ! would leak into btp%H_face_ave/Qu_face_ave/Qv_face_ave/
               ! Qw_face_ave -- exactly the consistency targets the volume/
               ! face routines' deficit-redistribution logic spreads onto
               ! the wet layers.
               if (pkl < (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff) then
                  pkl = (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff
                  ukl = 0.0;  vkl = 0.0;  wkl = 0.0
               end if
               if (pkr < (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff) then
                  pkr = (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff
                  ukr = 0.0;  vkr = 0.0;  wkr = 0.0
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

               if (has_w) then
                  Qu_ql3 = Qu_ql3 + ukl * wkl * ope_ppl_k;  Qu_qr3 = Qu_qr3 + ukr * wkr * ope_ppr_k
                  Qv_ql3 = Qv_ql3 + vkl * wkl * ope_ppl_k;  Qv_qr3 = Qv_qr3 + vkr * wkr * ope_ppr_k
                  Qw_ql1 = Qw_ql1 + ukl * wkl * ope_ppl_k;  Qw_qr1 = Qw_qr1 + ukr * wkr * ope_ppr_k
                  Qw_ql2 = Qw_ql2 + vkl * wkl * ope_ppl_k;  Qw_qr2 = Qw_qr2 + vkr * wkr * ope_ppr_k
                  Qw_ql3 = Qw_ql3 + wkl * wkl * ope_ppl_k;  Qw_qr3 = Qw_qr3 + wkr * wkr * ope_ppr_k
               end if

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
            flux_edge_z = 0.0
            if (has_w) flux_edge_z = 0.5*(qbl(5) + qbr(5)) + nzl * half_clam_dqb2
            flux_pb     = nxl*flux_edge_x + nyl*flux_edge_y + nzl*flux_edge_z

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

            if (has_w) then
               btp%btp_mass_flux_face_ave(3,iquad,iface) = btp%btp_mass_flux_face_ave(3,iquad,iface) + flux_edge_z
               btp%Qu_face_ave(3,iquad,iface) = btp%Qu_face_ave(3,iquad,iface) + 0.5*(Qu_ql3 + Qu_qr3)
               btp%Qv_face_ave(3,iquad,iface) = btp%Qv_face_ave(3,iquad,iface) + 0.5*(Qv_ql3 + Qv_qr3)
               btp%Qw_face_ave(1,iquad,iface) = btp%Qw_face_ave(1,iquad,iface) + 0.5*(Qw_ql1 + Qw_qr1)
               btp%Qw_face_ave(2,iquad,iface) = btp%Qw_face_ave(2,iquad,iface) + 0.5*(Qw_ql2 + Qw_qr2)
               btp%Qw_face_ave(3,iquad,iface) = btp%Qw_face_ave(3,iquad,iface) + 0.5*(Qw_ql3 + Qw_qr3)
               btp%uvb_face_ave(3,1,iquad,iface) = btp%uvb_face_ave(3,1,iquad,iface) + wl
               btp%uvb_face_ave(3,2,iquad,iface) = btp%uvb_face_ave(3,2,iquad,iface) + wr
            end if

            ! Momentum fluxes
            oe2_Hql = one_eta2 * H_bcl_ql
            oe2_Hqr = one_eta2 * H_bcl_qr

            Qu_ql1 = Qu_ql1 + oe2_Hql
            Qu_qr1 = Qu_qr1 + oe2_Hqr
            Qv_ql2 = Qv_ql2 + oe2_Hql
            Qv_qr2 = Qv_qr2 + oe2_Hqr
            if (has_w) then
               Qw_ql3 = Qw_ql3 + oe2_Hql
               Qw_qr3 = Qw_qr3 + oe2_Hqr
            end if

            flux_u = 0.5 * (nxl*(Qu_ql1 + Qu_qr1) + nyl*(Qu_ql2 + Qu_qr2) + nzl*(Qu_ql3 + Qu_qr3)) &
               - quarter_clam * (qbr(3) - qbl(3))
            flux_v = 0.5 * (nxl*(Qv_ql1 + Qv_qr1) + nyl*(Qv_ql2 + Qv_qr2) + nzl*(Qv_ql3 + Qv_qr3)) &
               - quarter_clam * (qbr(4) - qbl(4))
            flux_w = 0.0
            if (has_w) flux_w = 0.5 * (nxl*(Qw_ql1 + Qw_qr1) + nyl*(Qw_ql2 + Qw_qr2) + nzl*(Qw_ql3 + Qw_qr3)) &
               - quarter_clam * (qbr(5) - qbl(5))

            wq         = mf%jac_faceq(iquad, 1, iface)
            wq_flux_pb = wq * flux_pb
            wq_flux_u  = wq * flux_u
            wq_flux_v  = wq * flux_v
            wq_flux_w  = wq * flux_w

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
               if (has_w) then
                  !$acc atomic update
                  rhs_btp(4,I) = rhs_btp(4,I) - hi * wq_flux_w
               end if
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
                  if (has_w) then
                     !$acc atomic update
                     rhs_btp(4,I) = rhs_btp(4,I) + hi * wq_flux_w
                  end if
               end do
            end if

         end do  ! iquad

      end do  ! iface
      !$acc end data

   end subroutine create_btp_fluxes_qdf

end module mod_rhs_btp
