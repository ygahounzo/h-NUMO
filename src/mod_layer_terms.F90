! ===============================================================================================
! This module contains the routines for the baroclinic flux terms
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: March 27, 2023
! ===============================================================================================

module mod_layer_terms

    use mod_grid,      only: grid
    use mod_basis,     only: basis
    use mod_input,     only: input
    use mod_initial,   only: initial
    use mod_face,      only: face_CS
    use mod_tensor,    only: tensor_CS
    use mod_metrics,   only: metrics
    use mod_variables, only: btp_CS, bcl_CS
    use mod_constants, only: gravity
    use mod_initial_mlswe,    only: find_dry_elements

    implicit none

    public :: layer_mom_boundary_df,    &
              velocity_df,              &
              evaluate_consistency_face, &
              extract_qprime_df_face,   &
              interpolate_dpp

contains

    subroutine interpolate_dpp(G, inp, b, tsp, dprimeq, dprime_df)

        ! Interpolate layer mass from dofs to quad points

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(tensor_CS), intent(in)    :: tsp

        real, dimension(G%npoin_q, inp%nlayers), intent(inout) :: dprimeq
        real, dimension(G%npoin,   inp%nlayers), intent(inout) :: dprime_df

        integer :: Iq, I, ip
        real    :: hi

        dprimeq = 0.0

        do Iq = 1, G%npoin_q
            do ip = 1, b%npts
                I  = tsp%indexq(ip, Iq)
                hi = tsp%psih(ip, Iq)
                dprimeq(Iq,:) = dprimeq(Iq,:) + dprime_df(I,:) * hi
            end do
        end do

    end subroutine interpolate_dpp

    subroutine evaluate_consistency_face(G, inp, b, mf, btp, bcl, init, &
                                         mass_deficit_mass_face, dprime_df)

        ! Computes the mass deficit face values for consistency correction

        implicit none

        type(grid),    intent(in) :: G
        type(input),   intent(in) :: inp
        type(basis),   intent(in) :: b
        type(face_CS), intent(in) :: mf
        type(btp_CS),  intent(in) :: btp
        type(bcl_CS),  intent(in) :: bcl
        type(initial), intent(in) :: init

        real, intent(out) :: mass_deficit_mass_face(2, 2, b%nq, G%nface, inp%nlayers)
        real, intent(in)  :: dprime_df(G%npoin, inp%nlayers)

        integer :: iface, iquad, k, il, jl, kl, el, er, ir, jr, kr, I, n
        real    :: qprime_l, qprime_r, weights_face_l, weights_face_r, hi
        real    :: pbprime_l, pbprime_r

        mass_deficit_mass_face = 0.0

        do k = 1, inp%nlayers

            do iface = 1, G%nface

                el = G%face(7,iface)
                er = G%face(8,iface)

                do iquad = 1, b%nq

                    qprime_l = 0.0;  qprime_r = 0.0
                    pbprime_l = 0.0; pbprime_r = 0.0

                    do n = 1, b%ngl
                        hi = b%psiq(n,iquad)

                        il = mf%imapl(1,n,1,iface)
                        jl = mf%imapl(2,n,1,iface)
                        kl = mf%imapl(3,n,1,iface)

                        I = G%intma(il,jl,kl,el)
                        qprime_l  = qprime_l  + hi * dprime_df(I,k)
                        pbprime_l = pbprime_l + hi * init%pbprime_df(I)
                    enddo

                    if(er > 0) then
                        do n = 1, b%ngl
                            hi = b%psiq(n,iquad)
                            ir = mf%imapr(1,n,1,iface)
                            jr = mf%imapr(2,n,1,iface)
                            kr = mf%imapr(3,n,1,iface)
                            I = G%intma(ir,jr,kr,er)
                            qprime_r  = qprime_r  + hi * dprime_df(I,k)
                            pbprime_r = pbprime_r + hi * init%pbprime_df(I)
                        enddo
                    else
                        qprime_r  = qprime_l
                        pbprime_r = pbprime_l
                    endif

                    weights_face_l = qprime_l / pbprime_l
                    weights_face_r = qprime_r / pbprime_r

                    mass_deficit_mass_face(1,1,iquad,iface,k) = weights_face_l * &
                        (btp%btp_mass_flux_face_ave(1,iquad,iface) &
                         - bcl%sum_layer_mass_flux_face(1,iquad,iface))
                    mass_deficit_mass_face(2,1,iquad,iface,k) = weights_face_l * &
                        (btp%btp_mass_flux_face_ave(2,iquad,iface) &
                         - bcl%sum_layer_mass_flux_face(2,iquad,iface))

                    mass_deficit_mass_face(1,2,iquad,iface,k) = weights_face_r * &
                        (btp%btp_mass_flux_face_ave(1,iquad,iface) &
                         - bcl%sum_layer_mass_flux_face(1,iquad,iface))
                    mass_deficit_mass_face(2,2,iquad,iface,k) = weights_face_r * &
                        (btp%btp_mass_flux_face_ave(2,iquad,iface) &
                         - bcl%sum_layer_mass_flux_face(2,iquad,iface))
                end do
            end do
        end do

        call bcl_create_communicator(mass_deficit_mass_face, 2, inp%nlayers, b%nq)

    end subroutine evaluate_consistency_face

    subroutine velocity_df(G, inp, q_df, qb_df)

        implicit none

        type(grid),  intent(in)    :: G
        type(input), intent(in)    :: inp

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(inout) :: q_df
        real, dimension(inp%nvar_btp, G%npoin),              intent(in)    :: qb_df

        real    :: ubar, vbar, wbar
        integer :: I, k
        real    :: uv_df(inp%nvar_bcl-1, G%npoin, inp%nlayers)
        logical :: has_w

        has_w = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

        uv_df = 0.0

        do k = 1, inp%nlayers
            uv_df(1,:,k) = q_df(2,:,k) / q_df(1,:,k)
            uv_df(2,:,k) = q_df(3,:,k) / q_df(1,:,k)
            if (has_w) uv_df(3,:,k) = q_df(4,:,k) / q_df(1,:,k)
        end do

        do I = 1, G%npoin
            ubar = 0.0; vbar = 0.0; wbar = 0.0

            do k = 1, inp%nlayers
                ubar = ubar + uv_df(1,I,k) * q_df(1,I,k)
                vbar = vbar + uv_df(2,I,k) * q_df(1,I,k)
                if (has_w) wbar = wbar + uv_df(3,I,k) * q_df(1,I,k)
            end do

            if(qb_df(1,I) > 0.0) then
                ubar = ubar / qb_df(1,I)
                vbar = vbar / qb_df(1,I)
                if (has_w) wbar = wbar / qb_df(1,I)

                do k = 1, inp%nlayers
                    uv_df(1,I,k) = uv_df(1,I,k) - ubar + qb_df(3,I)/qb_df(1,I)
                    uv_df(2,I,k) = uv_df(2,I,k) - vbar + qb_df(4,I)/qb_df(1,I)
                    if (has_w) uv_df(3,I,k) = uv_df(3,I,k) - wbar + qb_df(5,I)/qb_df(1,I)
                end do
            else
                uv_df(:,I,:) = 0.0
            end if
        end do

        do k = 1, inp%nlayers
            q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
            q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
            if (has_w) q_df(4,:,k) = uv_df(3,:,k) * q_df(1,:,k)
        end do

    end subroutine velocity_df

    subroutine extract_velocity(G, inp, b, mt, tsp, init, bcl, uv_df, q_df, qb_df)

      implicit none

      type(grid),    intent(in)  :: G
      type(input),   intent(in)  :: inp
      type(basis),   intent(in)  :: b
      type(metrics), intent(in)  :: mt
      type(tensor_CS), intent(in) :: tsp
      type(initial), intent(in)    :: init
      type(bcl_CS),  intent(inout) :: bcl

      real, dimension(inp%nvar_bcl-1, G%npoin, inp%nlayers), intent(out) :: uv_df
      real, dimension(inp%nvar_bcl,   G%npoin, inp%nlayers), intent(in)  :: q_df
      real, dimension(inp%nvar_btp,   G%npoin),              intent(in)  :: qb_df

      real    :: ubar, vbar, wbar, wjac, wsum, mult
      integer :: I, k, e, n, m
      integer :: nlayers_l, nelem_l, nglx_l, ngly_l
      logical :: has_w
      real, parameter :: eps = 1.0e-20

      real :: dp_avg(inp%nlayers), udp_avg(inp%nlayers), vdp_avg(inp%nlayers)
      real :: wdp_avg(inp%nlayers)
      real :: dp_max(inp%nlayers), dp_min(inp%nlayers)
      real :: dp_cutoff1(inp%nlayers), dp_cutoff2(inp%nlayers), dp_range(inp%nlayers)
      real :: a(inp%nlayers), bc(inp%nlayers), c_td(inp%nlayers)
      real :: r(inp%nlayers, 3)
      real :: u_ave(inp%nlayers), v_ave(inp%nlayers), w_ave(inp%nlayers)
      real :: weight(inp%nlayers)

      ! GPU: classify dry elements (all arrays already on device).
      call find_dry_elements(G, inp, b, tsp, bcl%q_df, init%alpha_mlswe, bcl%dry_flg)

      nlayers_l = inp%nlayers
      nelem_l   = G%nelem
      nglx_l    = b%nglx
      ngly_l    = b%ngly
      has_w     = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

      ! Compute per-layer blending thresholds on CPU (small loop over nlayers).
      do k = 1, nlayers_l
          if (inp%h_cutoff1 < 1.0e-2) then
              dp_cutoff1(k) = (gravity / init%alpha_mlswe(k)) * 10.0
              dp_cutoff2(k) = (gravity / init%alpha_mlswe(k)) * 100.0
          else
              dp_cutoff1(k) = (gravity / init%alpha_mlswe(k)) * inp%h_cutoff1
              dp_cutoff2(k) = (gravity / init%alpha_mlswe(k)) * inp%h_cutoff2
          end if
          dp_range(k) = max(dp_cutoff2(k) - dp_cutoff1(k), eps)
      end do

      ! Copy blending thresholds to device once; shared read-only by both kernels below.
      !$acc data copyin(dp_cutoff1, dp_cutoff2, dp_range)

      ! Part 1: element-parallel velocity blending.
      ! One gang per element; all inner loops are sequential within the gang.
      ! Shared boundary nodes written by multiple gangs with different blended values —
      ! same determinism as the sequential CPU version (last writer wins).
      !$acc parallel loop gang &
      !$acc    present(G%intma, b%wglx, b%wgly, mt%jac, q_df, uv_df) &
      !$acc    firstprivate(nlayers_l, nglx_l, ngly_l, has_w) &
      !$acc    private(dp_avg, udp_avg, vdp_avg, wdp_avg, dp_max, dp_min, &
      !$acc            a, bc, c_td, r, u_ave, v_ave, w_ave, weight, &
      !$acc            wsum, wjac, mult)
      do e = 1, nelem_l

        wsum = 0.0
        !$acc loop seq
        do k = 1, nlayers_l
            dp_avg(k) = 0.0;  udp_avg(k) = 0.0;  vdp_avg(k) = 0.0;  wdp_avg(k) = 0.0
            dp_max(k) = -huge(1.0);  dp_min(k) = huge(1.0)
        end do

        !$acc loop seq
        do m = 1, ngly_l
            !$acc loop seq
            do n = 1, nglx_l
                I    = G%intma(n, m, 1, e)
                wjac = b%wglx(n) * b%wgly(m) * mt%jac(n, m, 1, e)
                wsum = wsum + wjac
                !$acc loop seq
                do k = 1, nlayers_l
                    dp_avg(k)  = dp_avg(k)  + wjac * q_df(1,I,k)
                    udp_avg(k) = udp_avg(k) + wjac * q_df(2,I,k)
                    vdp_avg(k) = vdp_avg(k) + wjac * q_df(3,I,k)
                    if (has_w) wdp_avg(k) = wdp_avg(k) + wjac * q_df(4,I,k)
                    dp_max(k)  = max(dp_max(k), q_df(1,I,k))
                    dp_min(k)  = min(dp_min(k), q_df(1,I,k))
                end do
            end do
        end do

        !$acc loop seq
        do k = 1, nlayers_l
          dp_avg(k)  = dp_avg(k)  / wsum
          udp_avg(k) = udp_avg(k) / wsum
          vdp_avg(k) = vdp_avg(k) / wsum
          if (has_w) wdp_avg(k) = wdp_avg(k) / wsum
        end do

        ! Build tridiagonal system for mass-weighted cell-average velocity.
        !$acc loop seq
        do k = 1, nlayers_l
          weight(k) = (dp_max(k) - dp_cutoff1(k)) / dp_range(k)
          weight(k) = max(min(weight(k), 1.0), 0.0)
          bc(k)   = 1.0
          r(k, 1) = weight(k) * udp_avg(k) / (dp_avg(k) + eps)
          r(k, 2) = weight(k) * vdp_avg(k) / (dp_avg(k) + eps)
          if (has_w) r(k, 3) = weight(k) * wdp_avg(k) / (dp_avg(k) + eps)
        end do
        a(1)            = 0.0
        c_td(1)         = -(1.0 - weight(1))
        a(nlayers_l)    = -(1.0 - weight(nlayers_l))
        c_td(nlayers_l) = 0.0
        !$acc loop seq
        do k = 2, nlayers_l - 1
          a(k)    = -0.5 * (1.0 - weight(k))
          c_td(k) = a(k)
        end do

        ! Forward sweep (Thomas algorithm).
        !$acc loop seq
        do k = 2, nlayers_l
            mult   = a(k) / bc(k-1)
            bc(k)  = bc(k)  - mult * c_td(k-1)
            r(k,1) = r(k,1) - mult * r(k-1,1)
            r(k,2) = r(k,2) - mult * r(k-1,2)
            if (has_w) r(k,3) = r(k,3) - mult * r(k-1,3)
        end do
        ! Back substitution.
        u_ave(nlayers_l) = r(nlayers_l,1) / bc(nlayers_l)
        v_ave(nlayers_l) = r(nlayers_l,2) / bc(nlayers_l)
        if (has_w) w_ave(nlayers_l) = r(nlayers_l,3) / bc(nlayers_l)
        !$acc loop seq
        do k = nlayers_l-1, 1, -1
          u_ave(k) = (r(k,1) - c_td(k) * u_ave(k+1)) / bc(k)
          v_ave(k) = (r(k,2) - c_td(k) * v_ave(k+1)) / bc(k)
          if (has_w) w_ave(k) = (r(k,3) - c_td(k) * w_ave(k+1)) / bc(k)
        end do

        ! Recompute weight using dp_min for the pointwise blending.
        !$acc loop seq
        do k = 1, nlayers_l
          weight(k) = (dp_min(k) - dp_cutoff1(k)) / dp_range(k)
          weight(k) = max(min(weight(k), 1.0), 0.0)
        end do

        ! Write blended velocities at all element nodes.
        !$acc loop seq
        do m = 1, ngly_l
          !$acc loop seq
          do n = 1, nglx_l
            I = G%intma(n, m, 1, e)
            !$acc loop seq
            do k = 1, nlayers_l
              uv_df(1,I,k) = weight(k) * q_df(2,I,k) / (q_df(1,I,k) + eps) &
                            + (1.0 - weight(k)) * u_ave(k)
              uv_df(2,I,k) = weight(k) * q_df(3,I,k) / (q_df(1,I,k) + eps) &
                            + (1.0 - weight(k)) * v_ave(k)
              if (has_w) uv_df(3,I,k) = weight(k) * q_df(4,I,k) / (q_df(1,I,k) + eps) &
                                       + (1.0 - weight(k)) * w_ave(k)
            end do
          end do
        end do

      end do
      !$acc end parallel loop

      ! Part 2: element-parallel barotropic consistency correction.
      ! One gang per element; uses dry_flg(e,k) to skip fully-dry layers, matching
      ! the original CPU semantics. Shared boundary nodes may be written by multiple
      ! gangs with different corrections — same non-determinism as the sequential
      ! CPU version (last writer wins).
      !$acc parallel loop gang &
      !$acc    present(G%intma, uv_df, q_df, qb_df, bcl%dry_flg) &
      !$acc    firstprivate(nlayers_l, nglx_l, ngly_l, has_w) &
      !$acc    private(ubar, vbar, wbar)
      do e = 1, nelem_l
        !$acc loop seq
        do m = 1, ngly_l
          !$acc loop seq
          do n = 1, nglx_l
            I = G%intma(n, m, 1, e)
            ubar = 0.0;  vbar = 0.0;  wbar = 0.0
            !$acc loop seq
            do k = 1, nlayers_l
              if (bcl%dry_flg(e, k) == 2) cycle
              ubar = ubar + uv_df(1,I,k) * q_df(1,I,k)
              vbar = vbar + uv_df(2,I,k) * q_df(1,I,k)
              if (has_w) wbar = wbar + uv_df(3,I,k) * q_df(1,I,k)
            end do
            if (qb_df(1,I) > 0.0) then
              ubar = ubar / qb_df(1,I)
              vbar = vbar / qb_df(1,I)
              if (has_w) wbar = wbar / qb_df(1,I)
              !$acc loop seq
              do k = 1, nlayers_l
                if (bcl%dry_flg(e, k) == 2) cycle
                uv_df(1,I,k) = uv_df(1,I,k) - (ubar - qb_df(3,I)/qb_df(1,I))
                uv_df(2,I,k) = uv_df(2,I,k) - (vbar - qb_df(4,I)/qb_df(1,I))
                if (has_w) uv_df(3,I,k) = uv_df(3,I,k) - (wbar - qb_df(5,I)/qb_df(1,I))
              end do
            else
              !$acc loop seq
              do k = 1, nlayers_l
                uv_df(1,I,k) = 0.0
                uv_df(2,I,k) = 0.0
                if (has_w) uv_df(3,I,k) = 0.0
              end do
            end if
          end do
        end do
      end do
      !$acc end parallel loop

      !$acc end data

    end subroutine extract_velocity

    subroutine extract_qprime_df_face(G, inp, b, mt, tsp, init, bcl, qprime_df, q_df, qb_df)

      implicit none

      type(grid),      intent(in) :: G
      type(input),     intent(in) :: inp
      type(basis),     intent(in) :: b
      type(metrics),   intent(in) :: mt
      type(tensor_CS), intent(in) :: tsp
      type(initial),   intent(in)    :: init
      type(bcl_CS),    intent(inout) :: bcl

      real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(out) :: qprime_df
      real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: q_df
      real, dimension(inp%nvar_btp, G%npoin),              intent(in)  :: qb_df

      integer :: k, I, npoin_l, nlayers_l
      real    :: ope
      real    :: uv_df(inp%nvar_bcl-1, G%npoin, inp%nlayers)
      logical :: has_w

      npoin_l   = G%npoin
      nlayers_l = inp%nlayers
      has_w     = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

      !$acc data create(uv_df)

      !$acc kernels present(qprime_df)
      qprime_df = 0.0
      !$acc end kernels

      call extract_velocity(G, inp, b, mt, tsp, init, bcl, uv_df, q_df, qb_df)

      ! Gang over nodes: accumulate ope = sum_k(h_k)/H0 sequentially, then
      ! write qprime.  No cross-node dependency — no atomics needed.
      !$acc parallel loop gang &
      !$acc    present(qprime_df, q_df, qb_df, uv_df, init%pbprime_df) &
      !$acc    firstprivate(npoin_l, nlayers_l, has_w) private(ope)
      do I = 1, npoin_l
          ope = 0.0
          !$acc loop seq
          do k = 1, nlayers_l
              ope = ope + q_df(1,I,k)
          end do
          ope = ope / init%pbprime_df(I)

          !$acc loop seq
          do k = 1, nlayers_l
              qprime_df(1,I,k) = q_df(1,I,k) / ope
              qprime_df(2,I,k) = uv_df(1,I,k) - qb_df(3,I)/qb_df(1,I)
              qprime_df(3,I,k) = uv_df(2,I,k) - qb_df(4,I)/qb_df(1,I)
              if (has_w) qprime_df(4,I,k) = uv_df(3,I,k) - qb_df(5,I)/qb_df(1,I)
          end do
      end do
      !$acc end parallel loop

      !$acc end data

    end subroutine extract_qprime_df_face

    subroutine layer_mom_boundary_df(G, inp, b, mf, init, q)

        implicit none

        type(grid),    intent(in)    :: G
        type(input),   intent(in)    :: inp
        type(basis),   intent(in)    :: b
        type(face_CS), intent(in)    :: mf
        type(initial), intent(in)    :: init

        real, intent(inout) :: q(inp%nvar_bcl, G%npoin, inp%nlayers)

        integer :: iface, n, il, jl, kl, el, er, I, k
        integer :: nface_l, ngl_l, nlayers_l
        real    :: nx, ny, nz, upnl_k
        logical :: has_w

        nface_l   = G%nface
        ngl_l     = b%ngl
        nlayers_l = inp%nlayers
        has_w     = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

        ! Sphere tangency constraint: remove radial momentum component at all
        ! nodes per layer so momentum stays tangent to the sphere surface.
        if (has_w) then
            !$acc parallel loop gang collapse(2) &
            !$acc    present(init%kvector, q) private(nx, ny, nz, upnl_k)
            do k = 1, nlayers_l
                do I = 1, G%npoin
                    nx     = init%kvector(1,I)
                    ny     = init%kvector(2,I)
                    nz     = init%kvector(3,I)
                    upnl_k = q(2,I,k)*nx + q(3,I,k)*ny + q(4,I,k)*nz
                    q(2,I,k) = q(2,I,k) - upnl_k*nx
                    q(3,I,k) = q(3,I,k) - upnl_k*ny
                    q(4,I,k) = q(4,I,k) - upnl_k*nz
                end do
            end do
            !$acc end parallel loop
        end if

        ! Gang over faces; seq over face nodes and layers.
        ! Atomics guard writes at corner nodes shared between two wall faces.
        !$acc data copy(q)

        !$acc parallel loop gang &
        !$acc    present(G%face, G%intma, mf%imapl, mf%normal_vector, q) &
        !$acc    firstprivate(nface_l, ngl_l, nlayers_l, has_w) &
        !$acc    private(il, jl, kl, el, er, I, nx, ny, nz, upnl_k)
        do iface = 1, nface_l

            el = G%face(7,iface)
            er = G%face(8,iface)

            if (er == -4) then
                !$acc loop seq
                do n = 1, ngl_l
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    nx = mf%normal_vector(1,n,1,iface)
                    ny = mf%normal_vector(2,n,1,iface)
                    nz = 0.0
                    if (has_w) nz = mf%normal_vector(3,n,1,iface)
                    !$acc loop seq
                    do k = 1, nlayers_l
                        upnl_k = q(2,I,k)*nx + q(3,I,k)*ny
                        if (has_w) upnl_k = upnl_k + q(4,I,k)*nz
                        !$acc atomic update
                        q(2,I,k) = q(2,I,k) - upnl_k*nx
                        !$acc atomic update
                        q(3,I,k) = q(3,I,k) - upnl_k*ny
                        if (has_w) then
                            !$acc atomic update
                            q(4,I,k) = q(4,I,k) - upnl_k*nz
                        end if
                    end do
                end do

            elseif (er == -2) then
                !$acc loop seq
                do n = 1, ngl_l
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    !$acc loop seq
                    do k = 1, nlayers_l
                        !$acc atomic write
                        q(2,I,k) = 0.0
                        !$acc atomic write
                        q(3,I,k) = 0.0
                        if (has_w) then
                            !$acc atomic write
                            q(4,I,k) = 0.0
                        end if
                    end do
                end do
            end if
        end do
        !$acc end parallel loop

        !$acc end data

        ! Solid body rotation: restore momentum direction from the fixed initial
        ! (t=0) analytic wind field, held constant for the whole run.
        if (has_w .and. trim(inp%test_case) == 'solid_body') then
            !$acc parallel loop gang collapse(2) present(q, init%q_df)
            do k = 1, nlayers_l
                do I = 1, G%npoin
                    q(2,I,k) = (init%q_df(2,I,k)/init%q_df(1,I,k)) * q(1,I,k)
                    q(3,I,k) = (init%q_df(3,I,k)/init%q_df(1,I,k)) * q(1,I,k)
                    q(4,I,k) = (init%q_df(4,I,k)/init%q_df(1,I,k)) * q(1,I,k)
                end do
            end do
            !$acc end parallel loop
        end if

    end subroutine layer_mom_boundary_df

end module mod_layer_terms
