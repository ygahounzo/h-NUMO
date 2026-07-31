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

    ! Sized by nc = inp%nvar_bcl-1 (2 momentum components flat, 3 on sphere)
    ! and looped/sectioned over ivar — no has_w branching, matching the
    ! convention in rhs_layer_shear_stress. qb_df's barotropic momentum
    ! components are contiguous at qb_df(3:2+nc,I) (u,v,[w]); confirmed
    ! against extract_velocity/extract_qprime_df_face's identical indexing.
    subroutine velocity_df(G, inp, q_df, qb_df)

        implicit none

        type(grid),  intent(in)    :: G
        type(input), intent(in)    :: inp

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(inout) :: q_df
        real, dimension(inp%nvar_btp, G%npoin),              intent(in)    :: qb_df

        real    :: bar(inp%nvar_bcl-1)
        integer :: I, k, ivar, nc
        real    :: uv_df(inp%nvar_bcl-1, G%npoin, inp%nlayers)

        nc = inp%nvar_bcl - 1

        uv_df = 0.0

        do k = 1, inp%nlayers
            do ivar = 1, nc
                uv_df(ivar,:,k) = q_df(ivar+1,:,k) / q_df(1,:,k)
            end do
        end do

        do I = 1, G%npoin
            bar = 0.0

            do k = 1, inp%nlayers
                bar(:) = bar(:) + uv_df(:,I,k) * q_df(1,I,k)
            end do

            if (qb_df(1,I) > 0.0) then
                bar(:) = bar(:) / qb_df(1,I)

                do k = 1, inp%nlayers
                    uv_df(:,I,k) = uv_df(:,I,k) - bar(:) + qb_df(3:2+nc,I)/qb_df(1,I)
                end do
            else
                uv_df(:,I,:) = 0.0
            end if
        end do

        do k = 1, inp%nlayers
            do ivar = 1, nc
                q_df(ivar+1,:,k) = uv_df(ivar,:,k) * q_df(1,:,k)
            end do
        end do

    end subroutine velocity_df

    ! Sized by nc = inp%nvar_bcl-1 (2 momentum components flat, 3 on sphere);
    ! has_w branching replaced by array-section slicing (q_df(2:nc+1,...),
    ! qb_df(3:2+nc,...)) since I/k are always fixed scalars at each point of
    ! use, matching the convention in rhs_layer_shear_stress/velocity_df.
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

      real    :: wjac, wsum, mult, btp_threshold
      integer :: I, k, e, n, m, nc
      integer :: nlayers_l, nelem_l, nglx_l, ngly_l
      real, parameter :: eps = 1.0e-20

      real :: dp_avg(inp%nlayers)
      real :: mdp_avg(inp%nvar_bcl-1, inp%nlayers)
      real :: dp_max(inp%nlayers), dp_min(inp%nlayers)
      real :: dp_cutoff1(inp%nlayers), dp_cutoff2(inp%nlayers), dp_range(inp%nlayers)
      real :: a(inp%nlayers), bc(inp%nlayers), c_td(inp%nlayers)
      real :: r(inp%nlayers, inp%nvar_bcl-1)
      real :: vel_ave(inp%nvar_bcl-1, inp%nlayers)
      real :: weight(inp%nlayers)
      real :: bar(inp%nvar_bcl-1)

      ! GPU: classify dry elements (all arrays already on device).
      call find_dry_elements(G, inp, b, tsp, bcl%q_df, init%alpha_mlswe, bcl%dry_flg)

      nlayers_l = inp%nlayers
      nelem_l   = G%nelem
      nglx_l    = b%nglx
      ngly_l    = b%ngly
      nc        = inp%nvar_bcl - 1

      ! Compute per-layer blending thresholds on CPU (small loop over nlayers).
      btp_threshold = 0.0
      do k = 1, nlayers_l
          if (inp%h_cutoff1 < 1.0e-2) then
              dp_cutoff1(k) = (gravity / init%alpha_mlswe(k)) * 10.0
              dp_cutoff2(k) = (gravity / init%alpha_mlswe(k)) * 100.0
          else
              dp_cutoff1(k) = (gravity / init%alpha_mlswe(k)) * inp%h_cutoff1
              dp_cutoff2(k) = (gravity / init%alpha_mlswe(k)) * inp%h_cutoff2
          end if
          dp_range(k) = max(dp_cutoff2(k) - dp_cutoff1(k), eps)
          ! BTP threshold: sum of per-layer floors (same convention as btp_poslimiter).
          btp_threshold = btp_threshold + (gravity / init%alpha_mlswe(k)) * inp%dry_cutoff
      end do

      ! Copy blending thresholds to device once; shared read-only by both kernels below.
      !$acc data copyin(dp_cutoff1, dp_cutoff2, dp_range)

      ! Part 1: element-parallel velocity blending.
      ! One gang per element; all inner loops are sequential within the gang.
      ! Shared boundary nodes written by multiple gangs with different blended values —
      ! same determinism as the sequential CPU version (last writer wins).
      !$acc parallel loop gang &
      !$acc    present(G%intma, b%wglx, b%wgly, mt%jac, q_df, uv_df) &
      !$acc    firstprivate(nlayers_l, nglx_l, ngly_l, nc) &
      !$acc    private(dp_avg, mdp_avg, dp_max, dp_min, &
      !$acc            a, bc, c_td, r, vel_ave, weight, &
      !$acc            wsum, wjac, mult)
      do e = 1, nelem_l

        wsum = 0.0
        !$acc loop seq
        do k = 1, nlayers_l
            dp_avg(k) = 0.0;  mdp_avg(:,k) = 0.0
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
                    dp_avg(k)    = dp_avg(k)    + wjac * q_df(1,I,k)
                    mdp_avg(:,k) = mdp_avg(:,k) + wjac * q_df(2:nc+1,I,k)
                    dp_max(k)    = max(dp_max(k), q_df(1,I,k))
                    dp_min(k)    = min(dp_min(k), q_df(1,I,k))
                end do
            end do
        end do

        !$acc loop seq
        do k = 1, nlayers_l
          dp_avg(k)    = dp_avg(k)    / wsum
          mdp_avg(:,k) = mdp_avg(:,k) / wsum
        end do

        ! Build tridiagonal system for mass-weighted cell-average velocity.
        !$acc loop seq
        do k = 1, nlayers_l
          weight(k) = (dp_max(k) - dp_cutoff1(k)) / dp_range(k)
          weight(k) = max(min(weight(k), 1.0), 0.0)
          bc(k)   = 1.0
          r(k, :) = weight(k) * mdp_avg(:,k) / (dp_avg(k) + eps)
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
            r(k,:) = r(k,:) - mult * r(k-1,:)
        end do
        ! Back substitution.
        vel_ave(:,nlayers_l) = r(nlayers_l,:) / bc(nlayers_l)
        !$acc loop seq
        do k = nlayers_l-1, 1, -1
          vel_ave(:,k) = (r(k,:) - c_td(k) * vel_ave(:,k+1)) / bc(k)
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
              uv_df(:,I,k) = weight(k) * q_df(2:nc+1,I,k) / (q_df(1,I,k) + eps) &
                            + (1.0 - weight(k)) * vel_ave(:,k)
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
      !$acc    firstprivate(nlayers_l, nglx_l, ngly_l, nc, btp_threshold) &
      !$acc    private(bar)
      do e = 1, nelem_l
        !$acc loop seq
        do m = 1, ngly_l
          !$acc loop seq
          do n = 1, nglx_l
            I = G%intma(n, m, 1, e)
            bar = 0.0
            !$acc loop seq
            do k = 1, nlayers_l
              if (bcl%dry_flg(e, k) == 2) cycle
              bar(:) = bar(:) + uv_df(:,I,k) * q_df(1,I,k)
            end do
            if (qb_df(1,I) > btp_threshold) then
              bar(:) = bar(:) / qb_df(1,I)
              !$acc loop seq
              do k = 1, nlayers_l
                if (bcl%dry_flg(e, k) == 2) cycle
                uv_df(:,I,k) = uv_df(:,I,k) - (bar(:) - qb_df(3:2+nc,I)/qb_df(1,I))
              end do
            else
              !$acc loop seq
              do k = 1, nlayers_l
                uv_df(:,I,k) = 0.0
              end do
            end if
          end do
        end do
      end do
      !$acc end parallel loop

      !$acc end data

    end subroutine extract_velocity

    ! Sized by nc = inp%nvar_bcl-1 (2 momentum components flat, 3 on sphere);
    ! has_w branching replaced by array-section slicing since I/k are always
    ! fixed scalars at each point of use, matching the convention in
    ! rhs_layer_shear_stress/velocity_df/extract_velocity.
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

      integer :: k, I, e, n, m, nelem_l, nlayers_l, nglx_l, ngly_l, nc
      real    :: ope, btp_threshold, max_qprime_spd, dp_threshold
      real    :: uv_df(inp%nvar_bcl-1, G%npoin, inp%nlayers)
      real    :: vel_btp(inp%nvar_bcl-1)

      nelem_l   = G%nelem
      nlayers_l = inp%nlayers
      nglx_l    = b%nglx
      ngly_l    = b%ngly
      nc        = inp%nvar_bcl - 1

      ! BTP threshold: same convention as btp_poslimiter.
      btp_threshold = 0.0
      do k = 1, nlayers_l
          btp_threshold = btp_threshold + (gravity / init%alpha_mlswe(k)) * inp%dry_cutoff
      end do

      ! Per-layer pressure threshold for the qprime positivity limiter.
      ! Approximate as btp_threshold/nlayers (exact when all layers have equal buoyancy).
      dp_threshold   = btp_threshold / real(nlayers_l)
      max_qprime_spd = 100.0

      !$acc data create(uv_df)

      !$acc kernels present(qprime_df)
      qprime_df = 0.0
      !$acc end kernels

      call extract_velocity(G, inp, b, mt, tsp, init, bcl, uv_df, q_df, qb_df)

      ! Element-parallel: one gang per element; uses dry_flg(e,k) to skip fully-dry
      ! layers (left at the qprime_df=0.0 init above), and guards the qb_df(1,I)
      ! divisions the same way extract_velocity's Part 2 barotropic correction does.
      ! Shared boundary nodes may be written by multiple gangs — same
      ! non-determinism as extract_velocity.
      !$acc parallel loop gang &
      !$acc    present(G%intma, qprime_df, q_df, qb_df, uv_df, init%pbprime_df, bcl%dry_flg) &
      !$acc    firstprivate(nelem_l, nlayers_l, nglx_l, ngly_l, nc, btp_threshold) &
      !$acc    private(ope, vel_btp)
      do e = 1, nelem_l
        !$acc loop seq
        do m = 1, ngly_l
          !$acc loop seq
          do n = 1, nglx_l
            I = G%intma(n, m, 1, e)

            ope = 0.0
            !$acc loop seq
            do k = 1, nlayers_l
                if (bcl%dry_flg(e, k) == 2) cycle
                ope = ope + q_df(1,I,k)
            end do
            ope = ope / init%pbprime_df(I)

            if (qb_df(1,I) > btp_threshold) then
              vel_btp(:) = qb_df(3:2+nc,I) / qb_df(1,I)
            else
              vel_btp(:) = 0.0
            end if

            !$acc loop seq
            do k = 1, nlayers_l
                if (bcl%dry_flg(e, k) == 2) cycle
                if (ope > 0.0) qprime_df(1,I,k) = q_df(1,I,k) / ope
                ! When BTP column is at minimum, ope is tiny and qprime_df(1) inflates
                ! by 1/ope. Zero the velocity deviation to prevent quadratic flux blow-up
                ! in btp_bcl_coeffs_qdf (flux += pp_k * up_k^2 with inflated pp_k).
                if (qb_df(1,I) > btp_threshold) then
                    qprime_df(2:nc+1,I,k) = uv_df(:,I,k) - vel_btp(:)
                else
                    qprime_df(2:nc+1,I,k) = 0.0
                end if
            end do
          end do
        end do
      end do
      !$acc end parallel loop

      ! Positivity limiter on qprime_df:
      !   - Near-dry layers (qprime_df(1) <= dp_threshold): zero velocity deviation.
      !   - All other nodes: smooth tanh saturation — u_soft = max_spd * tanh(u/max_spd).
      !     Unlike a hard clip, tanh is C-infinity smooth: no kink, no checkerboard noise.
      !     For |u| << max_spd it is nearly linear; for |u| >> max_spd it saturates to ±max_spd.
      ! !$acc parallel loop gang &
      ! !$acc    present(G%intma, qprime_df) &
      ! !$acc    firstprivate(nelem_l, nlayers_l, nglx_l, ngly_l, has_w, dp_threshold, max_qprime_spd)
      ! do e = 1, nelem_l
      !   !$acc loop seq
      !   do m = 1, ngly_l
      !     !$acc loop seq
      !     do n = 1, nglx_l
      !       I = G%intma(n, m, 1, e)
      !       !$acc loop seq
      !       do k = 1, nlayers_l
      !         if (qprime_df(1,I,k) <= dp_threshold) then
      !           qprime_df(2,I,k) = 0.0
      !           qprime_df(3,I,k) = 0.0
      !           if (has_w) qprime_df(4,I,k) = 0.0
      !         else
      !           if (qprime_df(2,I,k) >  max_qprime_spd) qprime_df(2,I,k) =  max_qprime_spd
      !           if (qprime_df(2,I,k) < -max_qprime_spd) qprime_df(2,I,k) = -max_qprime_spd
      !           if (qprime_df(3,I,k) >  max_qprime_spd) qprime_df(3,I,k) =  max_qprime_spd
      !           if (qprime_df(3,I,k) < -max_qprime_spd) qprime_df(3,I,k) = -max_qprime_spd
      !           if (has_w) then
      !             if (qprime_df(4,I,k) >  max_qprime_spd) qprime_df(4,I,k) =  max_qprime_spd
      !             if (qprime_df(4,I,k) < -max_qprime_spd) qprime_df(4,I,k) = -max_qprime_spd
      !           end if
      !         end if
      !       end do
      !     end do
      !   end do
      ! end do
      ! !$acc end parallel loop

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

        integer :: iface, n, il, jl, kl, el, er, I, k, ivar
        integer :: nface_l, ngl_l, nlayers_l, nc
        real    :: nx, ny, nz, upnl_k
        real    :: nvec(inp%nvar_bcl-1)
        logical :: has_w

        nface_l   = G%nface
        ngl_l     = b%ngl
        nlayers_l = inp%nlayers
        nc        = inp%nvar_bcl - 1
        has_w     = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex; still needed below -- the tangency projection and solid_body restore are sphere-only physics, not just wider components

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
        !$acc    firstprivate(nface_l, ngl_l, nlayers_l, nc) &
        !$acc    private(il, jl, kl, el, er, I, nvec, upnl_k)
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
                    ! Only ever reads components 1:nc of normal_vector (2 flat,
                    ! 3 sphere), same as the has_w-gated nz read this replaces.
                    nvec(1:nc) = mf%normal_vector(1:nc,n,1,iface)
                    !$acc loop seq
                    do k = 1, nlayers_l
                        ! Reflection v' = v - (v.n)n is dimension-agnostic
                        upnl_k = 0.0
                        !$acc loop seq
                        do ivar = 1, nc
                            upnl_k = upnl_k + q(ivar+1,I,k)*nvec(ivar)
                        end do
                        !$acc loop seq
                        do ivar = 1, nc
                            !$acc atomic update
                            q(ivar+1,I,k) = q(ivar+1,I,k) - upnl_k*nvec(ivar)
                        end do
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
                        !$acc loop seq
                        do ivar = 1, nc
                            !$acc atomic write
                            q(ivar+1,I,k) = 0.0
                        end do
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
