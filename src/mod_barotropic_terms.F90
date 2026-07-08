! ================================================================================================
! This module contains the routines for the barotropic flux terms
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: March 27, 2023
! ================================================================================================

module mod_barotropic_terms

    implicit none

    public :: &
                btp_mom_boundary_df, &
                btp_bcl_coeffs_qdf

    contains

    subroutine btp_mom_boundary_df(G, b, mf, init, inp, qb)

        use mod_basis,    only: basis
        use mod_grid,     only: grid
        use mod_face,     only: face_CS
        use mod_initial,  only: initial
        use mod_input,    only: input

        implicit none

        type(grid),    intent(in)    :: G
        type(basis),   intent(in)    :: b
        type(face_CS), intent(in)    :: mf
        type(initial), intent(in)    :: init
        type(input),   intent(in)    :: inp

        real, intent(inout) :: qb(inp%nvar_btp,G%npoin)

        integer :: iface, il, jl, el, er, I, kl, n
        real :: nx, ny, nz, unl
        logical :: has_w

        has_w = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex

        ! Sphere tangency constraint: remove radial momentum component at all nodes
        ! so momentum stays tangent to the sphere surface (kvector = r_hat).
        if (has_w) then
            !$acc parallel loop present(init%kvector, qb) private(nx, ny, nz, unl)
            do I = 1, G%npoin
                nx  = init%kvector(1,I)
                ny  = init%kvector(2,I)
                nz  = init%kvector(3,I)
                unl = qb(3,I)*nx + qb(4,I)*ny + qb(5,I)*nz
                qb(3,I) = qb(3,I) - unl*nx
                qb(4,I) = qb(4,I) - unl*ny
                qb(5,I) = qb(5,I) - unl*nz
            end do
            !$acc end parallel loop
        end if

        !$acc parallel loop present(G%face, G%intma, mf%imapl, mf%normal_vector, qb) &
        !$acc    private(il, jl, kl, el, er, I, nx, ny, nz, unl) firstprivate(has_w)
        do iface = 1, G%nface

            el = G%face(7,iface)
            er = G%face(8,iface)

            if (er == -4) then
                !$acc loop seq
                do n = 1, b%ngl
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    nx = mf%normal_vector(1,n,1,iface)
                    ny = mf%normal_vector(2,n,1,iface)
                    nz = 0.0
                    if (has_w) nz = mf%normal_vector(3,n,1,iface)
                    unl = qb(3,I)*nx + qb(4,I)*ny
                    if (has_w) unl = unl + qb(5,I)*nz
                    !$acc atomic update
                    qb(3,I) = qb(3,I) - unl*nx
                    !$acc atomic update
                    qb(4,I) = qb(4,I) - unl*ny
                    if (has_w) then
                        !$acc atomic update
                        qb(5,I) = qb(5,I) - unl*nz
                    end if
                end do

            elseif (er == -2) then
                !$acc loop seq
                do n = 1, b%ngl
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    !$acc atomic write
                    qb(3,I) = 0.0
                    !$acc atomic write
                    qb(4,I) = 0.0
                    if (has_w) then
                        !$acc atomic write
                        qb(5,I) = 0.0
                    end if
                end do
            end if
        end do
        !$acc end parallel loop

        ! Solid body rotation: restore momentum direction from the fixed initial
        ! (t=0) analytic wind field, held constant for the whole run.
        if (has_w .and. trim(inp%test_case) == 'solid_body') then
            !$acc parallel loop present(qb, init%qb_df)
            do I = 1, G%npoin
                qb(3,I) = (init%qb_df(3,I)/init%qb_df(1,I)) * qb(1,I)
                qb(4,I) = (init%qb_df(4,I)/init%qb_df(1,I)) * qb(1,I)
                qb(5,I) = (init%qb_df(5,I)/init%qb_df(1,I)) * qb(1,I)
            end do
            !$acc end parallel loop
        end if

    end subroutine btp_mom_boundary_df

    subroutine btp_bcl_coeffs_qdf(G, inp, b, tsp, bcl, btp, qprime_df, alpha_mlswe, pbprime_df)

        ! Compute baroclinic coefficients in the advective barotropic momentum fluxes
        ! and in the barotropic pressure forcing.  Also precomputes the BCL layer
        ! integrals at quad points (btp%bcl_*) so they need not be recomputed each
        ! BTP substep — qprime_df is fixed during BTP subcycling.

        use mod_grid,      only: grid
        use mod_input,     only: input
        use mod_basis,     only: basis
        use mod_constants, only: gravity
        use mod_tensor,    only: tensor_CS
        use mod_variables, only: bcl_CS, btp_CS

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(tensor_CS), intent(in)    :: tsp
        type(bcl_CS),    intent(inout) :: bcl
        type(btp_CS),    intent(inout) :: btp

        real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers), intent(in) :: qprime_df
        real, dimension(inp%nlayers),           intent(in) :: alpha_mlswe
        real, dimension(G%npoin),               intent(in) :: pbprime_df

        integer :: Iq, ip, k, npoin_l, nlayers_l, npts_l, npoin_q_l, I
        real :: dhdx, dhdy, dhdz, dpp
        real :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz
        real :: btg(9), pbps
        real :: pprime_k, pprime_k1
        real :: H_b, H_bclq, flux(3,3), pbq, hi, dry_k
        real :: pp_k, up_k, vp_k, wp_k
        logical :: has_w

        npoin_l   = G%npoin
        npoin_q_l = G%npoin_q
        nlayers_l = inp%nlayers
        npts_l    = b%npts
        has_w     = (inp%nvar_bcl == 4) ! w (vertical velocity) is only carried on sphere_hex

        ! Fused single kernel: gang over nodes, seq over layers and neighbours.
        ! Eliminates the graduv temporary and the sequential k-loop with its
        ! 3*nlayers kernel launches.  btp sums accumulated privately per node.
        !$acc parallel loop gang &
        !$acc    present(bcl%dpprime_visc, bcl%dpp_graduvw, bcl%dpp_uvp, &
        !$acc            btp%btp_dpp_graduvw, &
        !$acc            btp%pbprime_visc, qprime_df, &
        !$acc            tsp%index_df, tsp%dpsidx_df, tsp%dpsidy_df, &
        !$acc            tsp%dpsidz_df, tsp%dpsidz_df_x, tsp%dpsidz_df_y, &
        !$acc            tsp%dpsidz_df_z) &
        !$acc    firstprivate(npoin_l, nlayers_l, npts_l, has_w) &
        !$acc    private(dhdx, dhdy, dhdz, dpp, &
        !$acc            du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, &
        !$acc            dw_dx, dw_dy, dw_dz, &
        !$acc            btg, pbps, I)
        do Iq = 1, npoin_l
            btg = 0.0; pbps = 0.0
            !$acc loop seq
            do k = 1, nlayers_l
                dpp = qprime_df(1,Iq,k)
                bcl%dpprime_visc(Iq,k) = dpp

                du_dx = 0.0; du_dy = 0.0; du_dz = 0.0
                dv_dx = 0.0; dv_dy = 0.0; dv_dz = 0.0
                dw_dx = 0.0; dw_dy = 0.0; dw_dz = 0.0
                !$acc loop seq
                do ip = 1, npts_l
                    I    = tsp%index_df(ip,Iq)
                    dhdx = tsp%dpsidx_df(ip,Iq)
                    dhdy = tsp%dpsidy_df(ip,Iq)
                    dhdz = 0.0
                    if (has_w) then
                        dhdx = dhdx + tsp%dpsidz_df_x(ip,Iq)
                        dhdy = dhdy + tsp%dpsidz_df_y(ip,Iq)
                        dhdz = tsp%dpsidz_df(ip,Iq) + tsp%dpsidz_df_z(ip,Iq)
                    end if
                    du_dx = du_dx + dhdx*qprime_df(2,I,k)
                    du_dy = du_dy + dhdy*qprime_df(2,I,k)
                    dv_dx = dv_dx + dhdx*qprime_df(3,I,k)
                    dv_dy = dv_dy + dhdy*qprime_df(3,I,k)
                    if (has_w) then
                        du_dz = du_dz + dhdz*qprime_df(2,I,k)
                        dv_dz = dv_dz + dhdz*qprime_df(3,I,k)
                        dw_dx = dw_dx + dhdx*qprime_df(4,I,k)
                        dw_dy = dw_dy + dhdy*qprime_df(4,I,k)
                        dw_dz = dw_dz + dhdz*qprime_df(4,I,k)
                    end if
                end do

                ! x/y gradient components — sphere-corrected metric (no size change)
                bcl%dpp_graduvw(1,Iq,k) = dpp*du_dx
                bcl%dpp_graduvw(2,Iq,k) = dpp*du_dy
                bcl%dpp_graduvw(3,Iq,k) = dpp*dv_dx
                bcl%dpp_graduvw(4,Iq,k) = dpp*dv_dy

                bcl%dpp_uvp(1,Iq,k) = dpp*qprime_df(2,Iq,k)
                bcl%dpp_uvp(2,Iq,k) = dpp*qprime_df(3,Iq,k)

                btg(1:4) = btg(1:4) + bcl%dpp_graduvw(1:4,Iq,k)
                pbps = pbps + dpp

                ! Sphere-only: z-gradient of u/v (for future 3D laplacian) then full 3D grad(w).
                if (has_w) then
                    bcl%dpp_graduvw(5,Iq,k) = dpp*du_dz
                    bcl%dpp_graduvw(6,Iq,k) = dpp*dv_dz
                    bcl%dpp_graduvw(7,Iq,k) = dpp*dw_dx
                    bcl%dpp_graduvw(8,Iq,k) = dpp*dw_dy
                    bcl%dpp_graduvw(9,Iq,k) = dpp*dw_dz

                    bcl%dpp_uvp(3,Iq,k) = dpp*qprime_df(4,Iq,k)

                    btg(5:9) = btg(5:9) + bcl%dpp_graduvw(5:9,Iq,k)
                end if
            end do
            btp%pbprime_visc(Iq) = pbps
            if (has_w) then
                btp%btp_dpp_graduvw(:,Iq) = btg(1:9)
            else
                btp%btp_dpp_graduvw(:,Iq) = btg(1:4)
            end if
        end do
        !$acc end parallel loop

        ! Precompute BCL layer integrals at quad points.
        ! qprime_df is fixed during BTP subcycling, so this runs once per RK stage.
        !$acc parallel loop gang                                                   &
        !$acc    private(pp_k, up_k, vp_k, wp_k, pprime_k, pprime_k1,            &
        !$acc            H_b, H_bclq, flux, pbq, hi, I)                          &
        !$acc    firstprivate(npoin_q_l, nlayers_l, npts_l, has_w)                 &
        !$acc    present(tsp%indexq, tsp%psih, qprime_df, alpha_mlswe, pbprime_df, &
        !$acc            btp%bcl_H, btp%bcl_flux, btp%bcl_btp_flux, btp%pbq)
        do Iq = 1, npoin_q_l
           H_b = 0.0
           pbq = 0.0;  pprime_k = 0.0

           flux = 0.0

           ! Interpolate pbprime_df to quad point (static — same for all stages).
           !$acc loop seq
           do ip = 1, npts_l
              I   = tsp%indexq(ip, Iq)
              pbq = pbq + tsp%psih(ip, Iq) * pbprime_df(I)
           end do

           ! Layer loop: accumulate H_b and the velocity-product flux tensor.
           !$acc loop seq
           do k = 1, nlayers_l
              pp_k = 0.0;  up_k = 0.0;  vp_k = 0.0;  wp_k = 0.0

              !$acc loop seq
              do ip = 1, npts_l
                 I    = tsp%indexq(ip, Iq)
                 hi   = tsp%psih(ip, Iq)
                 pp_k = pp_k + hi * qprime_df(1, I, k)
                 up_k = up_k + hi * qprime_df(2, I, k)
                 vp_k = vp_k + hi * qprime_df(3, I, k)
                 if (has_w) wp_k = wp_k + hi * qprime_df(4, I, k)
              end do

              ! Dry-cell protection, matching the volume/face BCL flux
              ! routines. Thickness is CLAMPED (not zeroed) to a floor so
              ! H_bclq's cumulative pprime_k bookkeeping stays continuous
              ! across dry_cutoff; velocity is zeroed, which alone fully
              ! suppresses this layer's flux terms (every flux(i,j) entry
              ! below is a product of two velocity components) regardless of
              ! the floored thickness. Without this, a phantom contribution
              ! would leak into bcl_H/bcl_flux, then into
              ! btp%H_ave/Qu_ave/Qv_ave (accumulated every BTP substep in
              ! mod_rhs_btp.F90) -- exactly the consistency targets the
              ! volume/face routines' deficit-redistribution logic spreads
              ! onto the wet layers.
              if (pp_k < (gravity/alpha_mlswe(k)) * inp%dry_cutoff) then
                 pp_k = (gravity/alpha_mlswe(k)) * inp%dry_cutoff
                 up_k = 0.0;  vp_k = 0.0;  wp_k = 0.0
              end if

              pprime_k1 = pprime_k + pp_k
              H_bclq = 0.5 * alpha_mlswe(k) * (pprime_k1**2 - pprime_k**2)
              H_b = H_b + H_bclq

              ! flux(i,j) = sum_k vel_i,k * (pp_k * vel_j,k), for velocity components
              ! (1=u,2=v,3=w). The 3rd row/col only exists (and is only computed) on sphere.
              flux(1,1) = flux(1,1) + up_k*(pp_k*up_k)
              flux(2,1) = flux(2,1) + vp_k*(pp_k*up_k)
              flux(1,2) = flux(1,2) + up_k*(pp_k*vp_k)
              flux(2,2) = flux(2,2) + vp_k*(pp_k*vp_k)
              if (has_w) then
                 flux(3,1) = flux(3,1) + wp_k*(pp_k*up_k)
                 flux(3,2) = flux(3,2) + wp_k*(pp_k*vp_k)
                 flux(1,3) = flux(1,3) + up_k*(pp_k*wp_k)
                 flux(2,3) = flux(2,3) + vp_k*(pp_k*wp_k)
                 flux(3,3) = flux(3,3) + wp_k*(pp_k*wp_k)
              end if

              pprime_k = pprime_k1
           end do
           ! After the k loop, pp_k/up_k/vp_k/wp_k hold the bottom-layer values.

           btp%bcl_H(Iq) = H_b
           btp%pbq(Iq)   = pbq
           if (has_w) then
              btp%bcl_flux(:,:,Iq)    = flux(1:3,1:3)
              btp%bcl_btp_flux(1,Iq)  = pp_k
              btp%bcl_btp_flux(2,Iq)  = up_k
              btp%bcl_btp_flux(3,Iq)  = vp_k
              btp%bcl_btp_flux(4,Iq)  = wp_k
           else
              btp%bcl_flux(:,:,Iq)    = flux(1:2,1:2)
              btp%bcl_btp_flux(1,Iq)  = pp_k
              btp%bcl_btp_flux(2,Iq)  = up_k
              btp%bcl_btp_flux(3,Iq)  = vp_k
           end if
        end do
        !$acc end parallel loop

    end subroutine btp_bcl_coeffs_qdf

end module mod_barotropic_terms
