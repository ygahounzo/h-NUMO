! =======================================================================
!>@brief This subroutine builds the RHS vector for the DG method for MLSWE
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: April 03, 2023
! =======================================================================
module mod_create_rhs_mlswe

    use mod_constants,        only: gravity
    use mod_grid,             only: grid
    use mod_basis,            only: basis
    use mod_input,            only: input
    use mod_initial,          only: initial
    use mod_metrics,          only: metrics
    use mod_face,             only: face_CS
    use mod_tensor,           only: tensor_CS
    use mod_variables,        only: btp_CS, bcl_CS
    use mod_parallel,         only: parallel_CS
    use mod_ref,              only: mref
    use mod_mpi_communicator, only: mpi_communicator

    implicit none

    public :: rhs_layer_shear_stress, layer_mass_rhs, &
              bcl_rhs, bcl_apply_implicit_vertical_viscosity

contains

    subroutine layer_mass_rhs(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, dp_advec, qprime_df)

        implicit none

        type(grid),             intent(in)    :: G
        type(input),            intent(in)    :: inp
        type(basis),            intent(in)    :: b
        type(face_CS),          intent(in)    :: mf
        type(parallel_CS),      intent(in)    :: par
        type(btp_CS),           intent(inout) :: btp
        type(bcl_CS),           intent(inout) :: bcl
        type(initial),          intent(in)    :: init
        type(mref),             intent(inout) :: ref
        type(mpi_communicator), intent(inout) :: mpic
        type(tensor_CS),        intent(in)    :: tsp
        type(metrics),          intent(in)    :: mt

        real, dimension(G%npoin, inp%nlayers),    intent(out) :: dp_advec
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: qprime_df

        integer :: k

        call bcl_create_precommunicator(G, inp, b, mf, par, ref, mpic, qprime_df)
        if (trim(inp%geometry_type) == 'sphere_hex' .or. trim(inp%geometry_type) == 'sphere_ico') then
            call create_layers_volume_mass_sphere(G, inp, b, btp, bcl, init, tsp, dp_advec, qprime_df)
            call create_layer_mass_flux_sphere(G, inp, b, mf, btp, bcl, init, dp_advec, qprime_df)
        else
            call create_layers_volume_mass(G, inp, b, btp, bcl, init, tsp, dp_advec, qprime_df)
            call create_layer_mass_flux(G, inp, b, mf, btp, bcl, init, dp_advec, qprime_df)
        end if
        call bcl_create_postcommunicator_continuity(G, inp, b, mf, par, btp, init, ref, mpic, dp_advec)

        do k = 1, inp%nlayers
            dp_advec(:,k) = mt%massinv(:) * dp_advec(:,k)
        end do

    end subroutine layer_mass_rhs

    subroutine bcl_rhs(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs, qprime_df, q_df)

        implicit none

        type(grid),             intent(in)    :: G
        type(input),            intent(in)    :: inp
        type(basis),            intent(in)    :: b
        type(face_CS),          intent(in)    :: mf
        type(parallel_CS),      intent(in)    :: par
        type(btp_CS),           intent(inout) :: btp
        type(bcl_CS),           intent(inout) :: bcl
        type(initial),          intent(in)    :: init
        type(mref),             intent(inout) :: ref
        type(mpi_communicator), intent(inout) :: mpic
        type(tensor_CS),        intent(in)    :: tsp
        type(metrics),          intent(in)    :: mt

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(out) :: rhs
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: q_df, qprime_df

        if (trim(inp%geometry_type) == 'sphere_hex' .or. trim(inp%geometry_type) == 'sphere_ico') then
            call bcl_rhs_sphere(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs, qprime_df, q_df)
            return
        end if

        !$acc data present(qprime_df, rhs) copyin(q_df)
        call bcl_create_precommunicator(G, inp, b, mf, par, ref, mpic, qprime_df)
        call create_rhs_dynamics_volume_bcl_qp(G, inp, b, btp, bcl, init, tsp, rhs, qprime_df, q_df)
        call Apply_bcl_fluxes(G, inp, b, mf, btp, bcl, init, rhs, qprime_df)
        call bcl_create_postcommunicator(G, inp, b, mf, par, btp, init, ref, mpic, rhs)
        !$acc end data
        ! rhs stays on device; create_rhs_bcl owns the download after the combining loop.

    end subroutine bcl_rhs

    ! sphere_hex entry point for the BCL RHS (mass,u,v,w).
    subroutine bcl_rhs_sphere(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs, qprime_df, q_df)

        implicit none

        type(grid),             intent(in)    :: G
        type(input),            intent(in)    :: inp
        type(basis),            intent(in)    :: b
        type(face_CS),          intent(in)    :: mf
        type(parallel_CS),      intent(in)    :: par
        type(btp_CS),           intent(inout) :: btp
        type(bcl_CS),           intent(inout) :: bcl
        type(initial),          intent(in)    :: init
        type(mref),             intent(inout) :: ref
        type(mpi_communicator), intent(inout) :: mpic
        type(tensor_CS),        intent(in)    :: tsp
        type(metrics),          intent(in)    :: mt

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(out) :: rhs
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: q_df
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: qprime_df

        !$acc data present(qprime_df, rhs) copyin(q_df)

        call bcl_create_precommunicator(G, inp, b, mf, par, ref, mpic, qprime_df)

        call create_rhs_dynamics_volume_bcl_sphere(G, inp, b, btp, bcl, init, tsp, rhs, qprime_df, q_df)

        call Apply_bcl_fluxes_sphere(G, inp, b, mf, btp, bcl, init, rhs, qprime_df)

        call bcl_create_postcommunicator_sphere(G, inp, b, mf, par, btp, init, ref, mpic, rhs)

        !$acc end data

    end subroutine bcl_rhs_sphere

    ! Inter-layer vertical viscosity via a per-quadrature-point tridiagonal
    ! (Thomas algorithm) solve, identical formula to the do_ad_visc block
    ! already used inside create_rhs_dynamics_volume_bcl_qp/_sphere -- kept
    ! standalone so it can be applied once per full timestep by
    ! bcl_apply_implicit_vertical_viscosity rather than once per RK stage.
    ! Sized by nc = inp%nvar_bcl-1 (2 momentum components flat, 3 on sphere)
    ! and looped over with ivar, matching the generic inp%nvar_bcl-loop
    ! convention already used elsewhere (pack_and_send_df_bcl,
    ! unpack_data_dg_general_bcl) -- no has_w branching, and flat runs never
    ! touch a 3rd component at all. init%coriolis_quad(Iq) is the scalar
    ! Coriolis magnitude, already computed identically for flat and sphere
    ! geometry (mod_initial_mlswe.F90) -- equal to sqrt(cfx**2+cfy**2+cfz**2)
    ! on sphere since coriolis_3d_quad = f*r_hat and |r_hat|=1, so no
    ! separate sphere formula is needed here.
    subroutine rhs_layer_shear_stress(G, inp, b, init, tsp, rhs_stress, q_df)

        implicit none

        type(grid),      intent(in) :: G
        type(input),     intent(in) :: inp
        type(basis),     intent(in) :: b
        type(initial),   intent(in) :: init
        type(tensor_CS), intent(in) :: tsp

        real, intent(in)  :: q_df(inp%nvar_bcl, G%npoin, inp%nlayers)
        real, intent(out) :: rhs_stress(inp%nvar_bcl-1, G%npoin, inp%nlayers)

        real, dimension(inp%nvar_bcl-1, inp%nlayers+1) :: tau
        real, dimension(inp%nlayers)                   :: a, bc, c, dp
        real, dimension(inp%nvar_bcl-1, inp%nlayers)   :: mdp, r, uv
        real, dimension(inp%nvar_bcl-1) :: tau_q
        real    :: wq, hi, coeff, mult, coeff1
        integer :: k, Iq, I, ip, ivar, nc

        nc = inp%nvar_bcl - 1

        rhs_stress = 0.0

        do Iq = 1, G%npoin_q

            dp = 0.0; mdp = 0.0
            do ip = 1, b%npts
                I  = tsp%indexq(ip,Iq)
                hi = tsp%psih(ip,Iq)
                dp(:) = dp(:) + hi * q_df(1,I,:)
                do ivar = 1, nc
                    mdp(ivar,:) = mdp(ivar,:) + hi * q_df(ivar+1,I,:)
                end do
            end do

            coeff  = max(sqrt(0.5*init%coriolis_quad(Iq)*inp%ad_mlswe)/init%alpha_mlswe(1), &
                         inp%ad_mlswe/(init%alpha_mlswe(1) * inp%max_shear_dz))
            coeff1 = gravity * inp%dt * coeff

            do k = 1, inp%nlayers
                a(k)   = -coeff
                bc(k)  = dp(k) + 2.0*coeff1
                c(k)   = -coeff1
                r(:,k) = mdp(:,k) / dp(k)
            end do

            bc(1)           = dp(1) + coeff1
            bc(inp%nlayers) = dp(inp%nlayers) + coeff1
            a(1)            = 0.0
            c(inp%nlayers)  = 0.0

            do k = 2, inp%nlayers
                mult   = a(k) / bc(k-1)
                bc(k)  = bc(k) - mult*c(k-1)
                r(:,k) = r(:,k) - mult*r(:,k-1)
            end do

            r(:,inp%nlayers)  = r(:,inp%nlayers) / bc(inp%nlayers)
            uv(:,inp%nlayers) = r(:,inp%nlayers)

            do k = inp%nlayers-1, 1, -1
                r(:,k)  = (r(:,k) - c(k)*r(:,k+1)) / bc(k)
                uv(:,k) = r(:,k)
            end do

            ! Full zero (not just index 1) -- matches the live
            ! create_rhs_dynamics_volume_bcl_qp/_sphere pattern; the original
            ! index-1-only zeroing here left tau(:,nlayers+1) formally
            ! undefined (only correct in practice because this project's
            ! build uses -finit-real=zero).
            tau = 0.0
            do k = 2, inp%nlayers
                tau(:,k) = coeff*(uv(:,k-1) - uv(:,k))
            end do

            wq = tsp%wjac(Iq)

            do k = 1, inp%nlayers
                tau_q(:) = gravity*(tau(:,k) - tau(:,k+1))
                do ip = 1, b%npts
                    I  = tsp%indexq(ip,Iq)
                    hi = tsp%psih(ip,Iq)
                    rhs_stress(:,I,k) = rhs_stress(:,I,k) + wq*hi*tau_q(:)
                end do
            end do

        end do

    end subroutine rhs_layer_shear_stress

    ! Applies the ad_mlswe inter-layer vertical viscosity term implicitly,
    ! once per full timestep (Thetis Eq. 43 pattern), instead of embedding it
    ! inside the per-RK-stage RHS (create_rhs_dynamics_volume_bcl_qp/_sphere
    ! no longer compute it). No-op when inp%ad_mlswe<=0, so it is safe to
    ! call unconditionally after every integrator's explicit stages finish.
    subroutine bcl_apply_implicit_vertical_viscosity(G, inp, b, mf, init, tsp, mt, bcl, q_df, qb_df)

        use mod_layer_terms,   only: layer_mom_boundary_df, extract_velocity
        use mod_initial_mlswe, only: poslimiter

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(face_CS),   intent(in)    :: mf
        type(initial),   intent(in)    :: init
        type(tensor_CS), intent(in)    :: tsp
        type(metrics),   intent(in)    :: mt
        type(bcl_CS),    intent(inout) :: bcl

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(inout) :: q_df
        real, dimension(inp%nvar_btp, G%npoin),               intent(in)    :: qb_df

        real, dimension(inp%nvar_bcl-1, G%npoin, inp%nlayers) :: rhs_stress
        integer :: k, ivar

        if (inp%ad_mlswe <= 0.0) return

        call rhs_layer_shear_stress(G, inp, b, init, tsp, rhs_stress, q_df)

        do k = 1, inp%nlayers
            do ivar = 1, inp%nvar_bcl-1
                q_df(ivar+1,:,k) = q_df(ivar+1,:,k) + inp%dt*mt%massinv(:)*rhs_stress(ivar,:,k)
            end do
        end do

        call layer_mom_boundary_df(G, inp, b, mf, init, q_df)
        call poslimiter(b, G, inp, mt, q_df, init%alpha_mlswe)
        call extract_velocity(G, inp, b, mt, tsp, init, bcl, bcl%uv_df, q_df, qb_df)

        do k = 1, inp%nlayers
            do ivar = 1, inp%nvar_bcl-1
                q_df(ivar+1,:,k) = bcl%uv_df(ivar,:,k) * q_df(1,:,k)
            end do
        end do

        call poslimiter(b, G, inp, mt, q_df, init%alpha_mlswe)

    end subroutine bcl_apply_implicit_vertical_viscosity

    subroutine create_rhs_dynamics_volume_bcl_qp(G, inp, b, btp, bcl, init, tsp, rhs, qprime_df, q_df)

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(btp_CS),    intent(in)    :: btp
        type(bcl_CS),    intent(inout) :: bcl
        type(initial),   intent(in)    :: init
        type(tensor_CS), intent(in)    :: tsp

        real, dimension(3, G%npoin, inp%nlayers), intent(out) :: rhs
        real, dimension(3, G%npoin, inp%nlayers), intent(in)  :: q_df, qprime_df

        real :: wq, hi, dhdx, dhdy, tau_wind_u, tau_wind_v, temp1
        real :: Hq, source_x, source_y, Pstress, Pbstress
        integer :: k, I, Iq, ip
        real, dimension(inp%nlayers+1) :: pprime_temp, z
        real :: tempbot, weight, acceleration, pbq
        real, dimension(3) :: qp, qb
        real, dimension(inp%nlayers) :: temp_uu, temp_vv, H_tmp, u_udp, v_vdp, udp, vdp, dp
        real, dimension(2,inp%nlayers) :: u_vdp
        real :: p_tmp(inp%nlayers+1), u, v, weightq, one_over_sumuq, one_over_sumvq
        real :: uu_dp_deficitq, uv_dp_deficitq, vu_dp_deficitq, vv_dp_deficitq, gradz(2,inp%nlayers+1)
        real, parameter :: eps1 = 1.0e-20
        real :: flux(3,3)

        ! Scalar captures — inp/b/G are not on device; string comparison must stay on host.
        integer :: npoin_q_l, npts_l, nlayers_l

        npoin_q_l      = G%npoin_q
        npts_l         = b%npts
        nlayers_l      = inp%nlayers

        Pstress  = (gravity/init%alpha_mlswe(1))           * 50.0
        Pbstress = (gravity/init%alpha_mlswe(inp%nlayers)) * 10.0
        bcl%sum_layer_mass_flux = 0.0

        !$acc data present(tsp%indexq, tsp%psih, tsp%dpsidx, tsp%dpsidy, tsp%wjac,   &
        !$acc              qprime_df, rhs,                                              &
        !$acc              init%alpha_mlswe, init%zbot_df, init%pbprime_df,            &
        !$acc              init%grad_zbot_quad, init%tau_wind, init%coriolis_quad,     &
        !$acc              btp%ope_ave, btp%uvb_ave, btp%ope2_ave, btp%ope2_ave_df,   &
        !$acc              btp%H_ave, btp%btp_mass_flux_ave, btp%Qu_ave,               &
        !$acc              btp%Qv_ave, btp%tau_bot_ave)

        !$acc kernels present(rhs)
        rhs = 0.0
        !$acc end kernels

        !$acc parallel loop gang                                                        &
        !$acc   private(pprime_temp, z, qp, qb, temp_uu, temp_vv, H_tmp,              &
        !$acc           u_udp, v_vdp, u_vdp, udp, vdp, dp, p_tmp, gradz, flux,        &
        !$acc           wq, hi, dhdx, dhdy, tau_wind_u, tau_wind_v, temp1,             &
        !$acc           Hq, source_x, source_y, pbq, tempbot, weight, acceleration,   &
        !$acc           u, v, weightq, one_over_sumuq, one_over_sumvq,                 &
        !$acc           uu_dp_deficitq, uv_dp_deficitq, vu_dp_deficitq, vv_dp_deficitq, &
        !$acc           k, I, ip)                                 &
        !$acc   firstprivate(Pstress, Pbstress, npoin_q_l, npts_l, nlayers_l)
        do Iq = 1, npoin_q_l

            p_tmp(1)       = 0.0
            pprime_temp(:) = 0.0
            temp_uu = 0.0; temp_vv = 0.0

            !$acc loop seq
            do k = 1, nlayers_l
                qp = 0.0; qb = 0.0
                !$acc loop seq
                do ip = 1, npts_l
                    I  = tsp%indexq(ip,Iq)
                    hi = tsp%psih(ip,Iq)
                    qp(1) = qp(1) + hi*qprime_df(1,I,k)
                    qp(2) = qp(2) + hi*qprime_df(2,I,k)
                    qp(3) = qp(3) + hi*qprime_df(3,I,k)
                end do

                qb(1) = btp%ope_ave(Iq)
                qb(2) = btp%uvb_ave(1,Iq)
                qb(3) = btp%uvb_ave(2,Iq)

                p_tmp(k+1) = p_tmp(k) + sqrt(btp%ope2_ave(Iq)) * qp(1)
                H_tmp(k)   = 0.5*init%alpha_mlswe(k) * (p_tmp(k+1)**2 - p_tmp(k)**2)

                dp(k) = qp(1) * qb(1)
                u     = qp(2) + qb(2)
                v     = qp(3) + qb(3)

                udp(k)     = u * dp(k)
                vdp(k)     = v * dp(k)
                u_udp(k)   = u * udp(k)
                v_vdp(k)   = v * vdp(k)
                u_vdp(1,k) = v * udp(k)
                u_vdp(2,k) = u * vdp(k)

                temp_uu(k) = abs(udp(k)) + eps1
                temp_vv(k) = abs(vdp(k)) + eps1

                pprime_temp(k+1) = pprime_temp(k) + qp(1)
            end do

            gradz(:,:) = 0.0; pbq = 0.0
            !$acc loop seq
            do ip = 1, npts_l
                I = tsp%indexq(ip,Iq)
                z(nlayers_l+1) = init%zbot_df(I)
                !$acc loop seq
                do k = nlayers_l, 1, -1
                    z(k) = z(k+1) + (init%alpha_mlswe(k)/gravity) * &
                           (sqrt(btp%ope2_ave_df(I))*qprime_df(1,I,k))
                    gradz(1,k) = gradz(1,k) + tsp%dpsidx(ip,Iq) * z(k)
                    gradz(2,k) = gradz(2,k) + tsp%dpsidy(ip,Iq) * z(k)
                end do
                gradz(1,nlayers_l+1) = init%grad_zbot_quad(1,Iq)
                gradz(2,nlayers_l+1) = init%grad_zbot_quad(2,Iq)
                pbq = pbq + tsp%psih(ip,Iq) * init%pbprime_df(I)
            end do

            uu_dp_deficitq = btp%Qu_ave(1,Iq) - sum(u_udp(:))
            uv_dp_deficitq = btp%Qu_ave(2,Iq) - sum(u_vdp(1,:))
            vu_dp_deficitq = btp%Qv_ave(1,Iq) - sum(u_vdp(2,:))
            vv_dp_deficitq = btp%Qv_ave(2,Iq) - sum(v_vdp(:))

            one_over_sumuq = 1.0/sum(temp_uu(:))
            one_over_sumvq = 1.0/sum(temp_vv(:))

            wq = tsp%wjac(Iq)

            !$acc loop seq
            do k = 1, nlayers_l

                weightq    = temp_uu(k) * one_over_sumuq
                u_udp(k)   = u_udp(k)   + weightq * uu_dp_deficitq
                u_vdp(1,k) = u_vdp(1,k) + weightq * uv_dp_deficitq

                weightq    = temp_vv(k) * one_over_sumvq
                u_vdp(2,k) = u_vdp(2,k) + weightq * vu_dp_deficitq
                v_vdp(k)   = v_vdp(k)   + weightq * vv_dp_deficitq

                Hq = H_tmp(k)
                weight = 1.0
                acceleration = sum(H_tmp(:))
                if (acceleration > 0.0) weight = btp%H_ave(Iq) / acceleration
                Hq = Hq * weight

                flux(1,1) = udp(k) + (dp(k)/sum(dp(:))) * (btp%btp_mass_flux_ave(1,Iq) - sum(udp(:)))
                flux(2,1) = vdp(k) + (dp(k)/sum(dp(:))) * (btp%btp_mass_flux_ave(2,Iq) - sum(vdp(:)))

                flux(1,2) = u_udp(k) + Hq
                flux(2,2) = u_vdp(1,k)
                flux(1,3) = u_vdp(2,k)
                flux(2,3) = v_vdp(k) + Hq

                temp1      = (min(pprime_temp(k+1), Pstress) - min(pprime_temp(k), Pstress)) / Pstress
                tau_wind_u = temp1 * init%tau_wind(1,Iq)
                tau_wind_v = temp1 * init%tau_wind(2,Iq)

                tempbot = min(Pbstress, pbq-pprime_temp(k+1)) - min(Pbstress, pbq-pprime_temp(k))
                tempbot = tempbot / Pbstress

                source_x = gravity*(tau_wind_u - tempbot*btp%tau_bot_ave(1,Iq) + &
                            p_tmp(k)*gradz(1,k) - p_tmp(k+1)*gradz(1,k+1))
                source_y = gravity*(tau_wind_v - tempbot*btp%tau_bot_ave(2,Iq) + &
                            p_tmp(k)*gradz(2,k) - p_tmp(k+1)*gradz(2,k+1))

                ! ad_mlswe inter-layer vertical viscosity is no longer added
                ! here — applied once per full timestep by
                ! bcl_apply_implicit_vertical_viscosity instead.

                !$acc loop seq
                do ip = 1, npts_l
                    I    = tsp%indexq(ip,Iq)
                    hi   = tsp%psih(ip,Iq)
                    dhdx = tsp%dpsidx(ip,Iq)
                    dhdy = tsp%dpsidy(ip,Iq)
                    !$acc atomic update
                    rhs(1,I,k) = rhs(1,I,k) + wq*(dhdx*flux(1,1) + dhdy*flux(2,1))
                    !$acc atomic update
                    rhs(2,I,k) = rhs(2,I,k) + wq*(hi*source_x + dhdx*flux(1,2) + dhdy*flux(2,2))
                    !$acc atomic update
                    rhs(3,I,k) = rhs(3,I,k) + wq*(hi*source_y + dhdx*flux(1,3) + dhdy*flux(2,3))
                end do
            end do
        end do
        !$acc end parallel loop

        !$acc end data

    end subroutine create_rhs_dynamics_volume_bcl_qp

    ! sphere_hex counterpart of create_rhs_dynamics_volume_bcl_qp: carries the
    ! w-momentum row (mass,u,v,w) and the full 3D pressure gradient. Forked
    ! rather than has_w-gated, per the reference (hnumo-sphere) design (the
    ! bcl RHS pipeline there is itself split cartesian/sphere). Coriolis is
    ! never added here — handled by the semi-implicit rotation in each
    ! ti_*_bcl integrator instead.
    subroutine create_rhs_dynamics_volume_bcl_sphere(G, inp, b, btp, bcl, init, tsp, rhs, qprime_df, q_df)

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(btp_CS),    intent(in)    :: btp
        type(bcl_CS),    intent(inout) :: bcl
        type(initial),   intent(in)    :: init
        type(tensor_CS), intent(in)    :: tsp

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(out) :: rhs
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: q_df
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: qprime_df

        real :: wq, hi, dhdx, dhdy, dhdz, tau_wind_u, tau_wind_v, temp1
        real :: Hq, source_x, source_y, source_z, Pstress, Pbstress
        integer :: k, I, Iq, ip
        real :: pprime_temp(inp%nlayers+1), z(inp%nlayers+1), p_tmp(inp%nlayers+1)
        ! grad_lapz(d,k): gradient of ∇²z_k in direction d at quadrature point Iq.
        ! Used by Chen (2025) APE; k=1 and k=nlayers+1 stay 0 (no APE at boundaries).
        real :: grad_lapz(3, inp%nlayers+1)
        real :: tempbot, weight, acceleration, pbq
        real :: qp(4), qb(4)
        real :: dp(inp%nlayers), dpp(inp%nlayers), udp(inp%nlayers), vdp(inp%nlayers), wdp(inp%nlayers)
        real :: H_tmp(inp%nlayers)
        real :: temp_uu(inp%nlayers), temp_vv(inp%nlayers), temp_ww(inp%nlayers)
        real :: u_udp(3,inp%nlayers), v_vdp(3,inp%nlayers), w_wdp(3,inp%nlayers)
        real :: uu_dp_deficitq(3), vv_dp_deficitq(3), ww_dp_deficitq(3)
        real :: gradz(3, inp%nlayers+1)
        real :: u, v, w, weightq, weight_dp
        real :: one_over_sumuq, one_over_sumvq, one_over_sumwq
        real :: flux(3,4)
        real :: drho_over_rho
        real, parameter :: eps1 = 1.0e-20
        logical :: is_dry

        integer :: npoin_q_l, npts_l, nlayers_l
        real    :: Pstress_l, Pbstress_l

        npoin_q_l = G%npoin_q
        npts_l    = b%npts
        nlayers_l = inp%nlayers
        Pstress_l  = (gravity/init%alpha_mlswe(1))           * 50.0
        Pbstress_l = (gravity/init%alpha_mlswe(inp%nlayers)) * 10.0

        bcl%sum_layer_mass_flux = 0.0

        !$acc data present(tsp%indexq, tsp%psih, tsp%wjac,                              &
        !$acc              tsp%dpsidx, tsp%dpsidy, tsp%dpsidz,                           &
        !$acc              tsp%dpsidz_x, tsp%dpsidz_y, tsp%dpsidz_z,                    &
        !$acc              qprime_df, rhs,                                                &
        !$acc              init%alpha_mlswe, init%zbot_df, init%pbprime_df,              &
        !$acc              init%grad_zbot_quad, init%tau_wind,                            &
        !$acc              btp%ope_ave, btp%uvb_ave, btp%ope2_ave, btp%ope2_ave_df,     &
        !$acc              btp%H_ave, btp%btp_mass_flux_ave,                             &
        !$acc              btp%Qu_ave, btp%Qv_ave, btp%Qw_ave, btp%tau_bot_ave,         &
        !$acc              bcl%lap_z_df)

        !$acc kernels present(rhs)
        rhs = 0.0
        !$acc end kernels

        !$acc parallel loop gang                                                          &
        !$acc   private(pprime_temp, z, p_tmp, grad_lapz, qp, qb, gradz, flux,           &
        !$acc           dp, dpp, udp, vdp, wdp, H_tmp,                                  &
        !$acc           temp_uu, temp_vv, temp_ww,                                       &
        !$acc           u_udp, v_vdp, w_wdp,                                             &
        !$acc           uu_dp_deficitq, vv_dp_deficitq, ww_dp_deficitq,                 &
        !$acc           u, v, w, wq, hi, dhdx, dhdy, dhdz,              &
        !$acc           pbq, Hq, source_x, source_y, source_z,                          &
        !$acc           temp1, tau_wind_u, tau_wind_v, tempbot,                          &
        !$acc           weight, acceleration, weight_dp, weightq,                        &
        !$acc           one_over_sumuq, one_over_sumvq, one_over_sumwq, k, I, ip, is_dry, drho_over_rho)  &
        !$acc   firstprivate(Pstress_l, Pbstress_l, npoin_q_l, npts_l, nlayers_l)
        do Iq = 1, npoin_q_l

            p_tmp(1)       = 0.0
            pprime_temp(:) = 0.0
            temp_uu = 0.0;  temp_vv = 0.0;  temp_ww = 0.0

            ! ---------------------------------------------------------------
            !  Layer loop 1: interpolate primed state, compute velocities
            !  and momentum-flux tensor entries.
            ! ---------------------------------------------------------------
            !$acc loop seq
            do k = 1, nlayers_l
                qp = 0.0
                !$acc loop seq
                do ip = 1, npts_l
                    I  = tsp%indexq(ip,Iq)
                    hi = tsp%psih(ip,Iq)
                    qp(1) = qp(1) + hi*qprime_df(1,I,k)
                    qp(2) = qp(2) + hi*qprime_df(2,I,k)
                    qp(3) = qp(3) + hi*qprime_df(3,I,k)
                    qp(4) = qp(4) + hi*qprime_df(4,I,k)
                end do

                qb(1) = btp%ope_ave(Iq)
                qb(2) = btp%uvb_ave(1,Iq)
                qb(3) = btp%uvb_ave(2,Iq)
                qb(4) = btp%uvb_ave(3,Iq)

                ! Dry-cell protection. Thickness (qp(1)) is CLAMPED, not
                ! zeroed, to a floor so the H_tmp/p_tmp cumulative pressure
                ! bookkeeping stays continuous across the threshold (a hard
                ! zero would make H_tmp jump discontinuously as this layer
                ! crosses dry_cutoff). Velocity/momentum are ZEROED outright
                ! -- flux terms are velocity products, so this alone fully
                ! suppresses a thin-and-fast layer's flux contribution
                ! regardless of the (now-floored) thickness.
                !
                ! is_dry is decided ONCE here and reused below (rather than
                ! re-derived from dp(k) after clamping) to avoid dividing by
                ! qb(1) to land dp(k) exactly on the threshold -- qb(1) can
                ! be arbitrarily small, and a fresh dp(k)-based check after a
                ! non-divided clamp could disagree with this one whenever
                ! qb(1) isn't exactly 1.
                is_dry = (qp(1) * qb(1) < (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff)
                if (is_dry) qp(1) = (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff

                p_tmp(k+1) = p_tmp(k) + sqrt(btp%ope2_ave(Iq)) * qp(1)
                H_tmp(k)   = 0.5*init%alpha_mlswe(k) * (p_tmp(k+1)**2 - p_tmp(k)**2)

                dp(k)  = qp(1) * qb(1)
                dpp(k) = qp(1)

                if (is_dry) then
                    udp(k) = 0.0;  vdp(k) = 0.0;  wdp(k) = 0.0
                    u_udp(1,k) = 0.0;  u_udp(2,k) = 0.0;  u_udp(3,k) = 0.0
                    v_vdp(1,k) = 0.0;  v_vdp(2,k) = 0.0;  v_vdp(3,k) = 0.0
                    w_wdp(1,k) = 0.0;  w_wdp(2,k) = 0.0;  w_wdp(3,k) = 0.0
                    temp_uu(k) = eps1;  temp_vv(k) = eps1;  temp_ww(k) = eps1
                else
                    u = qp(2) + qb(2)
                    v = qp(3) + qb(3)
                    w = qp(4) + qb(4)

                    udp(k) = u*dp(k);  vdp(k) = v*dp(k);  wdp(k) = w*dp(k)

                    u_udp(1,k) = u*udp(k);  u_udp(2,k) = v*udp(k);  u_udp(3,k) = w*udp(k)
                    v_vdp(1,k) = u*vdp(k);  v_vdp(2,k) = v*vdp(k);  v_vdp(3,k) = w*vdp(k)
                    w_wdp(1,k) = u*wdp(k);  w_wdp(2,k) = v*wdp(k);  w_wdp(3,k) = w*wdp(k)

                    temp_uu(k) = abs(udp(k)) + eps1
                    temp_vv(k) = abs(vdp(k)) + eps1
                    temp_ww(k) = abs(wdp(k)) + eps1
                end if

                pprime_temp(k+1) = pprime_temp(k) + qp(1)
            end do

            ! ---------------------------------------------------------------
            !  3D gradient of layer interface elevations z_k.
            !  Sphere curvature corrections: dh/dx_i += dpsidz * dz/dx_i.
            ! ---------------------------------------------------------------
            gradz(:,:) = 0.0;  grad_lapz(:,:) = 0.0;  pbq = 0.0
            !$acc loop seq
            do ip = 1, npts_l
                I = tsp%indexq(ip,Iq)
                z(nlayers_l+1) = init%zbot_df(I)
                !$acc loop seq
                do k = nlayers_l, 1, -1
                    z(k) = z(k+1) + (init%alpha_mlswe(k)/gravity) * &
                           (sqrt(btp%ope2_ave_df(I))*qprime_df(1,I,k))
                    gradz(1,k) = gradz(1,k) + (tsp%dpsidx(ip,Iq)+tsp%dpsidz_x(ip,Iq))*z(k)
                    gradz(2,k) = gradz(2,k) + (tsp%dpsidy(ip,Iq)+tsp%dpsidz_y(ip,Iq))*z(k)
                    gradz(3,k) = gradz(3,k) + (tsp%dpsidz(ip,Iq)+tsp%dpsidz_z(ip,Iq))*z(k)
                end do
                ! Chen (2025) APE: (Δρ/ρ)×∇(∇²z_k) at this quad point for internal
                ! interfaces.  Density ratio (ρ_k−ρ_{k-1})/ρ_k = 1 − α_k/α_{k-1}
                ! is folded in here; boundary slots (k=1, k=nlayers+1) stay 0.
                if (inp%c_APE > 0.0) then
                    !$acc loop seq
                    do k = 2, nlayers_l
                        drho_over_rho = 1.0 - init%alpha_mlswe(k)/init%alpha_mlswe(k-1)
                        grad_lapz(1,k) = grad_lapz(1,k) + &
                            drho_over_rho*(tsp%dpsidx(ip,Iq)+tsp%dpsidz_x(ip,Iq))*bcl%lap_z_df(I,k)
                        grad_lapz(2,k) = grad_lapz(2,k) + &
                            drho_over_rho*(tsp%dpsidy(ip,Iq)+tsp%dpsidz_y(ip,Iq))*bcl%lap_z_df(I,k)
                        grad_lapz(3,k) = grad_lapz(3,k) + &
                            drho_over_rho*(tsp%dpsidz(ip,Iq)+tsp%dpsidz_z(ip,Iq))*bcl%lap_z_df(I,k)
                    end do
                end if
                gradz(1,nlayers_l+1) = init%grad_zbot_quad(1,Iq)
                gradz(2,nlayers_l+1) = init%grad_zbot_quad(2,Iq)
                gradz(3,nlayers_l+1) = init%grad_zbot_quad(3,Iq)
                pbq = pbq + tsp%psih(ip,Iq)*init%pbprime_df(I)
            end do

            ! ---------------------------------------------------------------
            !  Consistency deficit corrections for all three flux tensors.
            ! ---------------------------------------------------------------
            uu_dp_deficitq(1) = btp%Qu_ave(1,Iq) - sum(u_udp(1,:))
            uu_dp_deficitq(2) = btp%Qu_ave(2,Iq) - sum(u_udp(2,:))
            uu_dp_deficitq(3) = btp%Qu_ave(3,Iq) - sum(u_udp(3,:))

            vv_dp_deficitq(1) = btp%Qv_ave(1,Iq) - sum(v_vdp(1,:))
            vv_dp_deficitq(2) = btp%Qv_ave(2,Iq) - sum(v_vdp(2,:))
            vv_dp_deficitq(3) = btp%Qv_ave(3,Iq) - sum(v_vdp(3,:))

            ww_dp_deficitq(1) = btp%Qw_ave(1,Iq) - sum(w_wdp(1,:))
            ww_dp_deficitq(2) = btp%Qw_ave(2,Iq) - sum(w_wdp(2,:))
            ww_dp_deficitq(3) = btp%Qw_ave(3,Iq) - sum(w_wdp(3,:))

            one_over_sumuq = 1.0/sum(temp_uu(:))
            one_over_sumvq = 1.0/sum(temp_vv(:))
            one_over_sumwq = 1.0/sum(temp_ww(:))

            wq = tsp%wjac(Iq)

            ! ---------------------------------------------------------------
            !  Layer loop 2: assemble fluxes and accumulate into global rhs.
            ! ---------------------------------------------------------------
            !$acc loop seq
            do k = 1, nlayers_l

                ! Distribute deficit to each layer proportionally -- but only
                ! among wet layers; a (near-)dry layer shouldn't absorb any of
                ! the barotropic-consistency deficit correction.
                if (dp(k) > (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff) then
                    weightq = temp_uu(k)*one_over_sumuq
                    u_udp(1,k) = u_udp(1,k) + weightq*uu_dp_deficitq(1)
                    u_udp(2,k) = u_udp(2,k) + weightq*uu_dp_deficitq(2)
                    u_udp(3,k) = u_udp(3,k) + weightq*uu_dp_deficitq(3)

                    weightq = temp_vv(k)*one_over_sumvq
                    v_vdp(1,k) = v_vdp(1,k) + weightq*vv_dp_deficitq(1)
                    v_vdp(2,k) = v_vdp(2,k) + weightq*vv_dp_deficitq(2)
                    v_vdp(3,k) = v_vdp(3,k) + weightq*vv_dp_deficitq(3)

                    weightq = temp_ww(k)*one_over_sumwq
                    w_wdp(1,k) = w_wdp(1,k) + weightq*ww_dp_deficitq(1)
                    w_wdp(2,k) = w_wdp(2,k) + weightq*ww_dp_deficitq(2)
                    w_wdp(3,k) = w_wdp(3,k) + weightq*ww_dp_deficitq(3)
                end if

                ! Scale Hq so that sum over layers matches BTP-averaged H
                ! (skipped when this layer's own H_tmp is zero, i.e. dry).
                Hq = H_tmp(k)
                if (Hq > 0.0) then
                    weight = 1.0
                    acceleration = sum(H_tmp(:))
                    if (acceleration > 0.0) weight = btp%H_ave(Iq) / acceleration
                    Hq = Hq * weight
                end if

                ! Mass flux: consistency correction distributes BTP mass
                ! surplus among wet layers only; a dry layer's flux stays at
                ! its own (zero) value rather than absorbing a share of the
                ! surplus/deficit meant for layers that actually carry water.
                if (dp(k) > (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff) then
                    weight_dp  = abs(dp(k)) / (sum(abs(dp(:))) + eps1)
                    flux(1,1) = udp(k) + weight_dp*(btp%btp_mass_flux_ave(1,Iq) - sum(udp(:)))
                    flux(2,1) = vdp(k) + weight_dp*(btp%btp_mass_flux_ave(2,Iq) - sum(vdp(:)))
                    flux(3,1) = wdp(k) + weight_dp*(btp%btp_mass_flux_ave(3,Iq) - sum(wdp(:)))
                else
                    flux(1,1) = udp(k)
                    flux(2,1) = vdp(k)
                    flux(3,1) = wdp(k)
                end if

                ! 3×4 momentum flux tensor: [space dir, momentum component].
                flux(1,2) = u_udp(1,k) + Hq;  flux(2,2) = u_udp(2,k);       flux(3,2) = u_udp(3,k)
                flux(1,3) = v_vdp(1,k);        flux(2,3) = v_vdp(2,k) + Hq;  flux(3,3) = v_vdp(3,k)
                flux(1,4) = w_wdp(1,k);        flux(2,4) = w_wdp(2,k);        flux(3,4) = w_wdp(3,k) + Hq

                ! Wind penetration fraction and bottom drag fraction for this layer.
                temp1 = (min(pprime_temp(k+1), Pstress_l) - min(pprime_temp(k), Pstress_l)) / Pstress_l
                tau_wind_u = temp1 * init%tau_wind(1,Iq)
                tau_wind_v = temp1 * init%tau_wind(2,Iq)

                tempbot = (min(Pbstress_l, pbq-pprime_temp(k+1)) - &
                           min(Pbstress_l, pbq-pprime_temp(k))) / Pbstress_l

                ! 3D source: pressure gradient + wind/drag. Coriolis is never
                ! added here — handled by the semi-implicit rotation in each
                ! ti_*_bcl integrator instead.
                source_x = gravity*(p_tmp(k)*gradz(1,k) - p_tmp(k+1)*gradz(1,k+1)) &
                           + gravity*(tau_wind_u - tempbot*btp%tau_bot_ave(1,Iq))
                source_y = gravity*(p_tmp(k)*gradz(2,k) - p_tmp(k+1)*gradz(2,k+1)) &
                           + gravity*(tau_wind_v - tempbot*btp%tau_bot_ave(2,Iq))
                source_z = gravity*(p_tmp(k)*gradz(3,k) - p_tmp(k+1)*gradz(3,k+1)) &
                           - gravity*tempbot*btp%tau_bot_ave(3,Iq)

                ! Chen (2025) APE stabilization: force on layer k uses interface z_k (top of layer k).
                ! Chen Eq. 49: Φ̃^k = Φ^k − ν_k Δ(P^k/ρ^k)  →  momentum source = +ν_k ∇(∇²P^k/ρ^k)
                ! Dominant baroclinic term: +ν_k (Δρ/ρ)_k g ∇(∇²z_k).
                ! grad_lapz(:,k) carries (Δρ/ρ)_k × ∇(∇²z_k) from the accumulation above.
                ! grad_lapz(:,1) = 0 (free surface, ν₁=0 per Chen) → no force on layer 1.
                if (inp%c_APE > 0.0) then
                    source_x = source_x + gravity*inp%c_APE*dp(k)*grad_lapz(1,k)
                    source_y = source_y + gravity*inp%c_APE*dp(k)*grad_lapz(2,k)
                    source_z = source_z + gravity*inp%c_APE*dp(k)*grad_lapz(3,k)
                end if

                ! ad_mlswe inter-layer vertical viscosity is no longer added
                ! here — applied once per full timestep by
                ! bcl_apply_implicit_vertical_viscosity instead.

                !$acc loop seq
                do ip = 1, npts_l
                    I    = tsp%indexq(ip,Iq)
                    hi   = tsp%psih(ip,Iq)
                    dhdx = tsp%dpsidx(ip,Iq) + tsp%dpsidz_x(ip,Iq)
                    dhdy = tsp%dpsidy(ip,Iq) + tsp%dpsidz_y(ip,Iq)
                    dhdz = tsp%dpsidz(ip,Iq) + tsp%dpsidz_z(ip,Iq)
                    !$acc atomic update
                    rhs(1,I,k) = rhs(1,I,k) + wq*(dhdx*flux(1,1) + dhdy*flux(2,1) + dhdz*flux(3,1))
                    !$acc atomic update
                    rhs(2,I,k) = rhs(2,I,k) + wq*(hi*source_x + dhdx*flux(1,2) + dhdy*flux(2,2) + dhdz*flux(3,2))
                    !$acc atomic update
                    rhs(3,I,k) = rhs(3,I,k) + wq*(hi*source_y + dhdx*flux(1,3) + dhdy*flux(2,3) + dhdz*flux(3,3))
                    !$acc atomic update
                    rhs(4,I,k) = rhs(4,I,k) + wq*(hi*source_z + dhdx*flux(1,4) + dhdy*flux(2,4) + dhdz*flux(3,4))
                end do

            end do
        end do
        !$acc end parallel loop

        !$acc end data

    end subroutine create_rhs_dynamics_volume_bcl_sphere

    subroutine Apply_bcl_fluxes(G, inp, b, mf, btp, bcl, init, rhs, qprime_df)

        implicit none

        type(grid),    intent(in)    :: G
        type(input),   intent(in)    :: inp
        type(basis),   intent(in)    :: b
        type(face_CS), intent(in)    :: mf
        type(btp_CS),  intent(in)    :: btp
        type(bcl_CS),  intent(inout) :: bcl
        type(initial), intent(in)    :: init

        real, dimension(3, G%npoin, inp%nlayers), intent(inout) :: rhs
        real, dimension(3, G%npoin, inp%nlayers), intent(in)    :: qprime_df

        real, dimension(inp%nlayers)     :: alpha_over_g, g_over_alpha
        real, dimension(2,inp%nlayers+1) :: p_face, z_face
        real, dimension(inp%nlayers+1)   :: p_edge_plus, p_edge_minus, p2l, p2r
        real, dimension(inp%nlayers+1)   :: z_edge_plus, z_edge_minus
        real, dimension(3,inp%nlayers)   :: ql, qr
        real, dimension(3)               :: qbl, qbr
        integer :: iface, k, iquad, ktemp, I
        real :: z_intersect_top, z_intersect_bot, dz_intersect, H_r_plus, H_r_minus, acceleration
        real :: p_intersect_bot, p_intersect_top, one_plus_eta_edge
        real :: H_corr1, p_inc1, H_corr2, p_inc2, weight, ope_l, ope_r
        integer :: el, er
        real :: ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv, flux, un
        real, dimension(inp%nlayers) :: udpl, udpr, vdpl, vdpr
        real :: uu_dp_flux_deficit(2), vv_dp_flux_deficit(2), dp_deficit(2)
        real, parameter :: eps1 = 1.0e-20
        integer :: il, jl, ir, jr, kl, kr, n
        real :: wq, hi, flux_xl, flux_xr, flux_yl, flux_yr, hl, hr
        real, dimension(2,inp%nlayers,b%nq) :: H_face, udp_flux, vdp_flux, dp_flux
        real :: dp_lr(2,inp%nlayers)
        integer :: ngl_f, nq_f, nlayers_f, nface_f

        ngl_f     = b%ngl
        nq_f      = b%nq
        nlayers_f = inp%nlayers
        nface_f   = G%nface

        !$acc data present(G%face, G%face_type, G%intma,                                              &
        !$acc              mf%imapl, mf%imapr, mf%normal_vector_q, mf%jac_faceq,                     &
        !$acc              b%psiq,                                                                     &
        !$acc              init%alpha_mlswe, init%zbot_face,                                          &
        !$acc              btp%ope_face_ave, btp%uvb_face_ave, btp%ope2_face_ave,                     &
        !$acc              btp%one_plus_eta_edge_2_ave, btp%H_face_ave,                               &
        !$acc              btp%btp_mass_flux_face_ave, btp%Qu_face_ave, btp%Qv_face_ave,              &
        !$acc              qprime_df, rhs)

        !$acc parallel loop gang                                                                       &
        !$acc   private(alpha_over_g, g_over_alpha,                                                   &
        !$acc           ql, qr, qbl, qbr,                                                             &
        !$acc           p_face, z_face,                                                               &
        !$acc           p_edge_plus, p_edge_minus, p2l, p2r,                                          &
        !$acc           z_edge_plus, z_edge_minus,                                                    &
        !$acc           H_face, udp_flux, vdp_flux, dp_flux, dp_lr,                                  &
        !$acc           udpl, udpr, vdpl, vdpr,                                                       &
        !$acc           uu_dp_flux_deficit, vv_dp_flux_deficit, dp_deficit,                           &
        !$acc           el, er, k, iquad, ktemp, I, il, jl, ir, jr, kl, kr, n,                       &
        !$acc           nxl, nyl, wq, hi, uu, vv, un, flux,                                          &
        !$acc           flux_xl, flux_xr, flux_yl, flux_yr, hl, hr,                                  &
        !$acc           dpl, dpr, ul, ur, vl, vr, ope_l, ope_r, one_plus_eta_edge,                   &
        !$acc           weight, acceleration,                                                          &
        !$acc           H_r_plus, H_r_minus, z_intersect_top, z_intersect_bot, dz_intersect,          &
        !$acc           p_intersect_bot, p_intersect_top,                                             &
        !$acc           H_corr1, p_inc1, H_corr2, p_inc2)
        do iface = 1, nface_f

            if (G%face_type(iface) == 2) cycle

            el = G%face(7,iface)
            er = G%face(8,iface)

            !$acc loop seq
            do k = 1, nlayers_f
                alpha_over_g(k) = init%alpha_mlswe(k) / gravity
                g_over_alpha(k) = gravity / init%alpha_mlswe(k)
            end do

            !$acc loop seq
            do iquad = 1, nq_f

            qbl(1) = btp%ope_face_ave(1,iquad,iface)
            qbl(2) = btp%uvb_face_ave(1,1,iquad,iface)
            qbl(3) = btp%uvb_face_ave(2,1,iquad,iface)
            qbr(1) = btp%ope_face_ave(2,iquad,iface)
            qbr(2) = btp%uvb_face_ave(1,2,iquad,iface)
            qbr(3) = btp%uvb_face_ave(2,2,iquad,iface)

            nxl = mf%normal_vector_q(1,iquad,1,iface)
            nyl = mf%normal_vector_q(2,iquad,1,iface)

            ql = 0.0; qr = 0.0
            !$acc loop seq
            do k = 1, nlayers_f

                !$acc loop seq
                do n = 1, ngl_f
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    hi = b%psiq(n,iquad)
                    ql(1,k) = ql(1,k) + hi*qprime_df(1,I,k)
                    ql(2,k) = ql(2,k) + hi*qprime_df(2,I,k)
                    ql(3,k) = ql(3,k) + hi*qprime_df(3,I,k)
                end do

                if (er > 0) then
                    !$acc loop seq
                    do n = 1, ngl_f
                        ir = mf%imapr(1,n,1,iface)
                        jr = mf%imapr(2,n,1,iface)
                        kr = mf%imapr(3,n,1,iface)
                        I  = G%intma(ir,jr,kr,er)
                        hi = b%psiq(n,iquad)
                        qr(1,k) = qr(1,k) + hi*qprime_df(1,I,k)
                        qr(2,k) = qr(2,k) + hi*qprime_df(2,I,k)
                        qr(3,k) = qr(3,k) + hi*qprime_df(3,I,k)
                    end do
                else
                    qr(:,k) = ql(:,k)
                    if (er == -4) then
                        un = ql(2,k)*nxl + ql(3,k)*nyl
                        qr(2,k) = ql(2,k) - 2.0*un*nxl
                        qr(3,k) = ql(3,k) - 2.0*un*nyl
                    elseif (er == -2) then
                        qr(2,k) = -ql(2,k)
                        qr(3,k) = -ql(3,k)
                    end if
                end if

                dpl = qbl(1) * ql(1,k)
                dpr = qbr(1) * qr(1,k)
                ul  = ql(2,k) + qbl(2)
                ur  = qr(2,k) + qbr(2)
                vl  = ql(3,k) + qbl(3)
                vr  = qr(3,k) + qbr(3)

                dp_lr(1,k) = dpl
                dp_lr(2,k) = dpr

                uu = 0.5*(ul+ur)
                vv = 0.5*(vl+vr)
                udpl(k) = ul*dpl
                udpr(k) = ur*dpr
                vdpl(k) = vl*dpl
                vdpr(k) = vr*dpr

                un = uu*nxl + vv*nyl
                if(un > 0.0) then
                    dp_flux(1,k,iquad)  = uu * dpl
                    udp_flux(1,k,iquad) = uu * (ul*dpl)
                    vdp_flux(1,k,iquad) = uu * (vl*dpl)
                    dp_flux(2,k,iquad)  = vv * dpl
                    udp_flux(2,k,iquad) = vv * (ul*dpl)
                    vdp_flux(2,k,iquad) = vv * (vl*dpl)
                else
                    dp_flux(1,k,iquad)  = uu * dpr
                    udp_flux(1,k,iquad) = uu * (ur*dpr)
                    vdp_flux(1,k,iquad) = uu * (vr*dpr)
                    dp_flux(2,k,iquad)  = vv * dpr
                    udp_flux(2,k,iquad) = vv * (ur*dpr)
                    vdp_flux(2,k,iquad) = vv * (vr*dpr)
                end if

            end do

            uu_dp_flux_deficit(1) = btp%Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,:,iquad))
            uu_dp_flux_deficit(2) = btp%Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,:,iquad))
            vv_dp_flux_deficit(1) = btp%Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,:,iquad))
            vv_dp_flux_deficit(2) = btp%Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,:,iquad))

            dp_deficit(1) = btp%btp_mass_flux_face_ave(1,iquad,iface) - sum(dp_flux(1,:,iquad))
            dp_deficit(2) = btp%btp_mass_flux_face_ave(2,iquad,iface) - sum(dp_flux(2,:,iquad))

            !$acc loop seq
            do k = 1, nlayers_f

                weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:))+eps1))
                if ((dp_deficit(1)*nxl + dp_deficit(2)*nyl) < 0.0) &
                    weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:))+eps1))
                dp_flux(1,k,iquad) = dp_flux(1,k,iquad) + weight * dp_deficit(1)
                dp_flux(2,k,iquad) = dp_flux(2,k,iquad) + weight * dp_deficit(2)

                weight = abs(udpl(k)) / (sum(abs(udpl(:))+eps1))
                if ((uu_dp_flux_deficit(1)*nxl + uu_dp_flux_deficit(2)*nyl) < 0.0) &
                    weight = abs(udpr(k)) / (sum(abs(udpr(:))+eps1))
                udp_flux(1,k,iquad) = udp_flux(1,k,iquad) + weight * uu_dp_flux_deficit(1)
                udp_flux(2,k,iquad) = udp_flux(2,k,iquad) + weight * uu_dp_flux_deficit(2)

                weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
                if ((vv_dp_flux_deficit(1)*nxl + vv_dp_flux_deficit(2)*nyl) < 0.0) &
                    weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
                vdp_flux(1,k,iquad) = vdp_flux(1,k,iquad) + weight * vv_dp_flux_deficit(1)
                vdp_flux(2,k,iquad) = vdp_flux(2,k,iquad) + weight * vv_dp_flux_deficit(2)
            end do

            z_face = 0.0; p_face = 0.0
            z_edge_plus = 0.0; z_edge_minus = 0.0
            p_edge_plus = 0.0; p_edge_minus = 0.0

            ope_l = sqrt(btp%ope2_face_ave(1,iquad,iface))
            ope_r = sqrt(btp%ope2_face_ave(2,iquad,iface))
            p_face(1,1) = 0.0
            p_face(2,1) = 0.0
            !$acc loop seq
            do k = 1, nlayers_f
                p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,k)
                p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,k)
            end do

            one_plus_eta_edge = sqrt(btp%one_plus_eta_edge_2_ave(iquad,iface))
            z_face(1,nlayers_f+1)    = init%zbot_face(1,iquad,iface)
            z_face(2,nlayers_f+1)    = init%zbot_face(2,iquad,iface)
            z_edge_plus(nlayers_f+1) = init%zbot_face(1,iquad,iface)
            z_edge_minus(nlayers_f+1)= init%zbot_face(2,iquad,iface)
            !$acc loop seq
            do k = nlayers_f, 1, -1
                z_face(1,k)    = z_face(1,k+1)    + alpha_over_g(k) * (ope_l * ql(1,k))
                z_face(2,k)    = z_face(2,k+1)    + alpha_over_g(k) * (ope_r * qr(1,k))
                z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k) * (one_plus_eta_edge * ql(1,k))
                z_edge_minus(k)= z_edge_minus(k+1)+ alpha_over_g(k) * (one_plus_eta_edge * qr(1,k))
            end do

            p_edge_plus(2)  = one_plus_eta_edge * ql(1,1)
            p_edge_minus(2) = one_plus_eta_edge * qr(1,1)
            !$acc loop seq
            do k = 2, nlayers_f
                p_edge_plus(k+1)  = p_edge_plus(k)  + one_plus_eta_edge * ql(1,k)
                p_edge_minus(k+1) = p_edge_minus(k) + one_plus_eta_edge * qr(1,k)
            end do

            !$acc loop seq
            do k = 1, nlayers_f

                H_r_plus = 0.5*init%alpha_mlswe(k)*(p_edge_plus(k+1)**2 - p_edge_plus(k)**2)

                H_r_minus = 0.0
                !$acc loop seq
                do ktemp = 1, nlayers_f
                    z_intersect_top = min(z_edge_minus(ktemp),   z_edge_plus(k))
                    z_intersect_bot = max(z_edge_minus(ktemp+1), z_edge_plus(k+1))
                    dz_intersect    = z_intersect_top - z_intersect_bot
                    if (dz_intersect > 0.0) then
                        p_intersect_bot = p_edge_minus(ktemp+1) &
                                - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_minus(ktemp+1))
                        p_intersect_top = p_edge_minus(ktemp+1) &
                                - g_over_alpha(ktemp)*(z_intersect_top - z_edge_minus(ktemp+1))
                        H_r_minus = H_r_minus + &
                                0.5*init%alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
                    end if
                end do
                H_face(1,k,iquad) = 0.5*(H_r_plus + H_r_minus)

                H_r_minus = 0.5*init%alpha_mlswe(k)*(p_edge_minus(k+1)**2 - p_edge_minus(k)**2)

                H_r_plus = 0.0
                !$acc loop seq
                do ktemp = 1, nlayers_f
                    z_intersect_top = min(z_edge_plus(ktemp),   z_edge_minus(k))
                    z_intersect_bot = max(z_edge_plus(ktemp+1), z_edge_minus(k+1))
                    dz_intersect    = z_intersect_top - z_intersect_bot
                    if (dz_intersect > 0.0) then
                        p_intersect_bot = p_edge_plus(ktemp+1) &
                            - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_plus(ktemp+1))
                        p_intersect_top = p_edge_plus(ktemp+1) &
                            - g_over_alpha(ktemp)*(z_intersect_top - z_edge_plus(ktemp+1))
                        H_r_plus = H_r_plus &
                            + 0.5*init%alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
                    end if
                end do
                H_face(2,k,iquad) = 0.5*(H_r_plus + H_r_minus)
            end do

            if(er == -4) then
                p2l = 0.0; p2r = 0.0
                !$acc loop seq
                do k = 1, nlayers_f
                    p2l(k+1) = p_face(1,k+1)
                    H_face(1,k,iquad) = 0.5*init%alpha_mlswe(k)*(p2l(k+1)**2 - p2l(k)**2)
                    p2r(k+1) = p_face(2,k+1)
                    H_face(2,k,iquad) = 0.5*init%alpha_mlswe(k)*(p2r(k+1)**2 - p2r(k)**2)
                end do
            end if

            if(er /= -4) then
                !$acc loop seq
                do k = 1, nlayers_f-1
                    p_inc1  = p_face(1,k+1) + g_over_alpha(k)*(z_face(1,k+1) - z_edge_plus(k+1))
                    H_corr1 = 0.5*init%alpha_mlswe(k)*(p_inc1**2 - p_face(1,k+1)**2)
                    H_face(1,k,iquad)   = H_face(1,k,iquad)   - H_corr1
                    H_face(1,k+1,iquad) = H_face(1,k+1,iquad) + H_corr1

                    p_inc2  = p_face(2,k+1) + g_over_alpha(k)*(z_face(2,k+1) - z_edge_minus(k+1))
                    H_corr2 = 0.5*init%alpha_mlswe(k)*(p_inc2**2 - p_face(2,k+1)**2)
                    H_face(2,k,iquad)   = H_face(2,k,iquad)   - H_corr2
                    H_face(2,k+1,iquad) = H_face(2,k+1,iquad) + H_corr2
                end do
            end if

            weight = 1.0
            acceleration = sum(H_face(1,:,iquad))
            if(abs(acceleration) > eps1) weight = btp%H_face_ave(iquad,iface) / acceleration
            H_face(1,:,iquad) = H_face(1,:,iquad) * weight

            weight = 1.0
            acceleration = sum(H_face(2,:,iquad))
            if(abs(acceleration) > eps1) weight = btp%H_face_ave(iquad,iface) / acceleration
            H_face(2,:,iquad) = H_face(2,:,iquad) * weight

            wq = mf%jac_faceq(iquad,1,iface)

            !$acc loop seq
            do k = 1, nlayers_f

                flux = nxl*dp_flux(1,k,iquad) + nyl*dp_flux(2,k,iquad)

                hl = H_face(1,k,iquad)
                hr = H_face(2,k,iquad)

                flux_xl = nxl*(udp_flux(1,k,iquad) + hl) + nyl*udp_flux(2,k,iquad)
                flux_xr = nxl*(udp_flux(1,k,iquad) + hr) + nyl*udp_flux(2,k,iquad)
                flux_yl = nxl*vdp_flux(1,k,iquad) + nyl*(vdp_flux(2,k,iquad) + hl)
                flux_yr = nxl*vdp_flux(1,k,iquad) + nyl*(vdp_flux(2,k,iquad) + hr)

                !$acc loop seq
                do n = 1, ngl_f
                    hi = b%psiq(n,iquad)
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    !$acc atomic update
                    rhs(1,I,k) = rhs(1,I,k) - wq*hi*flux
                    !$acc atomic update
                    rhs(2,I,k) = rhs(2,I,k) - wq*hi*flux_xl
                    !$acc atomic update
                    rhs(3,I,k) = rhs(3,I,k) - wq*hi*flux_yl
                end do

                if(er > 0) then
                    !$acc loop seq
                    do n = 1, ngl_f
                        hi = b%psiq(n,iquad)
                        ir = mf%imapr(1,n,1,iface)
                        jr = mf%imapr(2,n,1,iface)
                        kr = mf%imapr(3,n,1,iface)
                        I  = G%intma(ir,jr,kr,er)
                        !$acc atomic update
                        rhs(1,I,k) = rhs(1,I,k) + wq*hi*flux
                        !$acc atomic update
                        rhs(2,I,k) = rhs(2,I,k) + wq*hi*flux_xr
                        !$acc atomic update
                        rhs(3,I,k) = rhs(3,I,k) + wq*hi*flux_yr
                    end do
                end if
            end do

            end do  ! iquad

        end do  ! iface
        !$acc end data

    end subroutine Apply_bcl_fluxes

    ! sphere_hex counterpart of Apply_bcl_fluxes: 4-wide (mass,u,v,w) local
    ! face Riemann fluxes with nzl normal and two-sided H_face(1,:,:)/H_face(2,:,:)
    ! interface-elevation correction, mirroring create_btp_fluxes_qdf's structure.
    subroutine Apply_bcl_fluxes_sphere(G, inp, b, mf, btp, bcl, init, rhs, qprime_df)

        implicit none

        type(grid),    intent(in)    :: G
        type(input),   intent(in)    :: inp
        type(basis),   intent(in)    :: b
        type(face_CS), intent(in)    :: mf
        type(btp_CS),  intent(in)    :: btp
        type(bcl_CS),  intent(inout) :: bcl
        type(initial), intent(in)    :: init

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(inout) :: rhs
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)    :: qprime_df

        real, dimension(inp%nlayers)     :: alpha_over_g, g_over_alpha
        real, dimension(2,inp%nlayers+1) :: p_face, z_face
        real, dimension(inp%nlayers+1)   :: p_edge_plus, p_edge_minus, p2l, p2r
        real, dimension(inp%nlayers+1)   :: z_edge_plus, z_edge_minus
        real, dimension(inp%nvar_bcl,inp%nlayers)   :: ql, qr
        real, dimension(inp%nvar_bcl)               :: qbl, qbr
        integer :: iface, k, iquad, ktemp, I
        real :: z_intersect_top, z_intersect_bot, dz_intersect, H_r_plus, H_r_minus, acceleration
        real :: p_intersect_bot, p_intersect_top, one_plus_eta_edge
        real :: H_corr1, p_inc1, H_corr2, p_inc2, weight, ope_l, ope_r
        integer :: el, er
        real :: ul, ur, vl, vr, wl, wr, dpl, dpr, nxl, nyl, nzl, uu, vv, ww, un
        real, dimension(inp%nlayers) :: udpl, udpr, vdpl, vdpr, wdpl, wdpr
        real :: uu_dp_flux_deficit(3), vv_dp_flux_deficit(3), ww_dp_flux_deficit(3), dp_deficit(3)
        real, parameter :: eps1 = 1.0e-20
        logical :: is_dry_l, is_dry_r
        integer :: il, jl, ir, jr, kl, kr, n
        real :: wq, hi, flux_xl, flux_xr, flux_yl, flux_yr, flux_zl, flux_zr, hl, hr
        real, dimension(2,inp%nlayers,b%nq) :: H_face
        real, dimension(3,inp%nlayers,b%nq) :: udp_flux, vdp_flux, dp_flux, wdp_flux
        real :: dp_lr(2,inp%nlayers)
        integer :: ngl_f, nq_f, nlayers_f, nface_f

        ngl_f     = b%ngl
        nq_f      = b%nq
        nlayers_f = inp%nlayers
        nface_f   = G%nface

        !$acc data present(G%face, G%face_type, G%intma,                                              &
        !$acc              mf%imapl, mf%imapr, mf%normal_vector_q, mf%jac_faceq,                     &
        !$acc              b%psiq,                                                                     &
        !$acc              init%alpha_mlswe, init%zbot_face,                                          &
        !$acc              btp%ope_face_ave, btp%uvb_face_ave, btp%ope2_face_ave,                     &
        !$acc              btp%one_plus_eta_edge_2_ave, btp%H_face_ave,                               &
        !$acc              btp%btp_mass_flux_face_ave, btp%Qu_face_ave, btp%Qv_face_ave,              &
        !$acc              btp%Qw_face_ave,                                                           &
        !$acc              qprime_df, rhs)

        !$acc parallel loop gang                                                                       &
        !$acc   private(alpha_over_g, g_over_alpha,                                                   &
        !$acc           ql, qr, qbl, qbr,                                                             &
        !$acc           p_face, z_face,                                                               &
        !$acc           p_edge_plus, p_edge_minus, p2l, p2r,                                          &
        !$acc           z_edge_plus, z_edge_minus,                                                    &
        !$acc           H_face, udp_flux, vdp_flux, dp_flux, wdp_flux, dp_lr,                        &
        !$acc           udpl, udpr, vdpl, vdpr, wdpl, wdpr,                                          &
        !$acc           uu_dp_flux_deficit, vv_dp_flux_deficit, ww_dp_flux_deficit, dp_deficit,      &
        !$acc           el, er, k, iquad, ktemp, I, il, jl, ir, jr, kl, kr, n,                       &
        !$acc           nxl, nyl, nzl, wq, hi, uu, vv, ww, un,                                       &
        !$acc           flux_xl, flux_xr, flux_yl, flux_yr, flux_zl, flux_zr, hl, hr,                &
        !$acc           dpl, dpr, ul, ur, vl, vr, wl, wr, ope_l, ope_r, one_plus_eta_edge,           &
        !$acc           weight, acceleration, is_dry_l, is_dry_r,                                     &
        !$acc           H_r_plus, H_r_minus, z_intersect_top, z_intersect_bot, dz_intersect,          &
        !$acc           p_intersect_bot, p_intersect_top,                                             &
        !$acc           H_corr1, p_inc1, H_corr2, p_inc2)
        do iface = 1, nface_f

            if (G%face_type(iface) == 2) cycle

            el = G%face(7,iface)
            er = G%face(8,iface)

            !$acc loop seq
            do k = 1, nlayers_f
                alpha_over_g(k) = init%alpha_mlswe(k) / gravity
                g_over_alpha(k) = gravity / init%alpha_mlswe(k)
            end do

            !$acc loop seq
            do iquad = 1, nq_f

            qbl(1) = btp%ope_face_ave(1,iquad,iface)
            qbl(2) = btp%uvb_face_ave(1,1,iquad,iface)
            qbl(3) = btp%uvb_face_ave(2,1,iquad,iface)
            qbl(4) = btp%uvb_face_ave(3,1,iquad,iface)
            qbr(1) = btp%ope_face_ave(2,iquad,iface)
            qbr(2) = btp%uvb_face_ave(1,2,iquad,iface)
            qbr(3) = btp%uvb_face_ave(2,2,iquad,iface)
            qbr(4) = btp%uvb_face_ave(3,2,iquad,iface)

            nxl = mf%normal_vector_q(1,iquad,1,iface)
            nyl = mf%normal_vector_q(2,iquad,1,iface)
            nzl = mf%normal_vector_q(3,iquad,1,iface)

            ql = 0.0; qr = 0.0
            !$acc loop seq
            do k = 1, nlayers_f

                !$acc loop seq
                do n = 1, ngl_f
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    hi = b%psiq(n,iquad)
                    ql(1,k) = ql(1,k) + hi*qprime_df(1,I,k)
                    ql(2,k) = ql(2,k) + hi*qprime_df(2,I,k)
                    ql(3,k) = ql(3,k) + hi*qprime_df(3,I,k)
                    ql(4,k) = ql(4,k) + hi*qprime_df(4,I,k)
                end do

                if (er > 0) then
                    !$acc loop seq
                    do n = 1, ngl_f
                        ir = mf%imapr(1,n,1,iface)
                        jr = mf%imapr(2,n,1,iface)
                        kr = mf%imapr(3,n,1,iface)
                        I  = G%intma(ir,jr,kr,er)
                        hi = b%psiq(n,iquad)
                        qr(1,k) = qr(1,k) + hi*qprime_df(1,I,k)
                        qr(2,k) = qr(2,k) + hi*qprime_df(2,I,k)
                        qr(3,k) = qr(3,k) + hi*qprime_df(3,I,k)
                        qr(4,k) = qr(4,k) + hi*qprime_df(4,I,k)
                    end do
                else
                    qr(:,k) = ql(:,k)
                    if (er == -4) then
                        un = ql(2,k)*nxl + ql(3,k)*nyl + ql(4,k)*nzl
                        qr(2,k) = ql(2,k) - 2.0*un*nxl
                        qr(3,k) = ql(3,k) - 2.0*un*nyl
                        qr(4,k) = ql(4,k) - 2.0*un*nzl
                    elseif (er == -2) then
                        qr(2,k) = -ql(2,k)
                        qr(3,k) = -ql(3,k)
                        qr(4,k) = -ql(4,k)
                    end if
                end if

                ! Dry-cell protection. Thickness is CLAMPED, not zeroed, to a
                ! floor so the H_face pressure bookkeeping (built below from
                ! z_face/p_face, which derive from ql(1,k)/qr(1,k)) stays
                ! continuous across dry_cutoff -- a hard zero would make it
                ! jump discontinuously as a layer crosses the threshold. The
                ! velocity perturbation IS zeroed: it keeps a dry layer's
                ! potentially noisy velocity estimate out of the uu/vv/ww
                ! averages shared with the other side, while the actual flux
                ! suppression for this layer comes from is_dry_l/is_dry_r
                ! below.
                !
                ! is_dry_l/is_dry_r are decided ONCE here and reused below
                ! (rather than re-derived from dpl/dpr after clamping) to
                ! avoid dividing by qbl(1)/qbr(1) to land dpl/dpr exactly on
                ! the threshold -- those can be arbitrarily small, and a
                ! fresh dpl/dpr-based check after a non-divided clamp could
                ! disagree with this one whenever qbl(1)/qbr(1) isn't
                ! exactly 1.
                is_dry_l = (ql(1,k) * qbl(1) < g_over_alpha(k) * inp%dry_cutoff)
                if (is_dry_l) then
                    ql(1,k) = g_over_alpha(k) * inp%dry_cutoff
                    ql(2,k) = 0.0;  ql(3,k) = 0.0;  ql(4,k) = 0.0
                end if
                is_dry_r = (qr(1,k) * qbr(1) < g_over_alpha(k) * inp%dry_cutoff)
                if (is_dry_r) then
                    qr(1,k) = g_over_alpha(k) * inp%dry_cutoff
                    qr(2,k) = 0.0;  qr(3,k) = 0.0;  qr(4,k) = 0.0
                end if

                dpl = qbl(1) * ql(1,k)
                dpr = qbr(1) * qr(1,k)
                ul  = ql(2,k) + qbl(2)
                ur  = qr(2,k) + qbr(2)
                vl  = ql(3,k) + qbl(3)
                vr  = qr(3,k) + qbr(3)
                wl  = ql(4,k) + qbl(4)
                wr  = qr(4,k) + qbr(4)

                dp_lr(1,k) = dpl
                dp_lr(2,k) = dpr

                uu = 0.5*(ul+ur)
                vv = 0.5*(vl+vr)
                ww = 0.5*(wl+wr)

                if (.not. is_dry_l) then
                    udpl(k) = ul*dpl
                    vdpl(k) = vl*dpl
                    wdpl(k) = wl*dpl
                else
                    udpl(k) = 0.0
                    vdpl(k) = 0.0
                    wdpl(k) = 0.0
                end if

                if (.not. is_dry_r) then
                    udpr(k) = ur*dpr
                    vdpr(k) = vr*dpr
                    wdpr(k) = wr*dpr
                else
                    udpr(k) = 0.0
                    vdpr(k) = 0.0
                    wdpr(k) = 0.0
                end if

                un = uu*nxl + vv*nyl + ww*nzl
                if(un > 0.0) then
                    dp_flux(1,k,iquad)  = uu * dpl
                    udp_flux(1,k,iquad) = uu * udpl(k)
                    vdp_flux(1,k,iquad) = uu * vdpl(k)
                    wdp_flux(1,k,iquad) = uu * wdpl(k)
                    dp_flux(2,k,iquad)  = vv * dpl
                    udp_flux(2,k,iquad) = vv * udpl(k)
                    vdp_flux(2,k,iquad) = vv * vdpl(k)
                    wdp_flux(2,k,iquad) = vv * wdpl(k)
                    dp_flux(3,k,iquad)  = ww * dpl
                    udp_flux(3,k,iquad) = ww * udpl(k)
                    vdp_flux(3,k,iquad) = ww * vdpl(k)
                    wdp_flux(3,k,iquad) = ww * wdpl(k)
                else
                    dp_flux(1,k,iquad)  = uu * dpr
                    udp_flux(1,k,iquad) = uu * udpr(k)
                    vdp_flux(1,k,iquad) = uu * vdpr(k)
                    wdp_flux(1,k,iquad) = uu * wdpr(k)
                    dp_flux(2,k,iquad)  = vv * dpr
                    udp_flux(2,k,iquad) = vv * udpr(k)
                    vdp_flux(2,k,iquad) = vv * vdpr(k)
                    wdp_flux(2,k,iquad) = vv * wdpr(k)
                    dp_flux(3,k,iquad)  = ww * dpr
                    udp_flux(3,k,iquad) = ww * udpr(k)
                    vdp_flux(3,k,iquad) = ww * vdpr(k)
                    wdp_flux(3,k,iquad) = ww * wdpr(k)
                end if

            end do

            uu_dp_flux_deficit(1) = btp%Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,:,iquad))
            uu_dp_flux_deficit(2) = btp%Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,:,iquad))
            uu_dp_flux_deficit(3) = btp%Qu_face_ave(3,iquad,iface) - sum(udp_flux(3,:,iquad))
            vv_dp_flux_deficit(1) = btp%Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,:,iquad))
            vv_dp_flux_deficit(2) = btp%Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,:,iquad))
            vv_dp_flux_deficit(3) = btp%Qv_face_ave(3,iquad,iface) - sum(vdp_flux(3,:,iquad))
            ww_dp_flux_deficit(1) = btp%Qw_face_ave(1,iquad,iface) - sum(wdp_flux(1,:,iquad))
            ww_dp_flux_deficit(2) = btp%Qw_face_ave(2,iquad,iface) - sum(wdp_flux(2,:,iquad))
            ww_dp_flux_deficit(3) = btp%Qw_face_ave(3,iquad,iface) - sum(wdp_flux(3,:,iquad))

            dp_deficit(1) = btp%btp_mass_flux_face_ave(1,iquad,iface) - sum(dp_flux(1,:,iquad))
            dp_deficit(2) = btp%btp_mass_flux_face_ave(2,iquad,iface) - sum(dp_flux(2,:,iquad))
            dp_deficit(3) = btp%btp_mass_flux_face_ave(3,iquad,iface) - sum(dp_flux(3,:,iquad))

            !$acc loop seq
            do k = 1, nlayers_f

                weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:))+eps1))
                if ((dp_deficit(1)*nxl + dp_deficit(2)*nyl + dp_deficit(3)*nzl) < 0.0) &
                    weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:))+eps1))
                dp_flux(1,k,iquad) = dp_flux(1,k,iquad) + weight * dp_deficit(1)
                dp_flux(2,k,iquad) = dp_flux(2,k,iquad) + weight * dp_deficit(2)
                dp_flux(3,k,iquad) = dp_flux(3,k,iquad) + weight * dp_deficit(3)

                weight = abs(udpl(k)) / (sum(abs(udpl(:))+eps1))
                if ((uu_dp_flux_deficit(1)*nxl + uu_dp_flux_deficit(2)*nyl + uu_dp_flux_deficit(3)*nzl) < 0.0) &
                    weight = abs(udpr(k)) / (sum(abs(udpr(:))+eps1))
                udp_flux(1,k,iquad) = udp_flux(1,k,iquad) + weight * uu_dp_flux_deficit(1)
                udp_flux(2,k,iquad) = udp_flux(2,k,iquad) + weight * uu_dp_flux_deficit(2)
                udp_flux(3,k,iquad) = udp_flux(3,k,iquad) + weight * uu_dp_flux_deficit(3)

                weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
                if ((vv_dp_flux_deficit(1)*nxl + vv_dp_flux_deficit(2)*nyl + vv_dp_flux_deficit(3)*nzl) < 0.0) &
                    weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
                vdp_flux(1,k,iquad) = vdp_flux(1,k,iquad) + weight * vv_dp_flux_deficit(1)
                vdp_flux(2,k,iquad) = vdp_flux(2,k,iquad) + weight * vv_dp_flux_deficit(2)
                vdp_flux(3,k,iquad) = vdp_flux(3,k,iquad) + weight * vv_dp_flux_deficit(3)

                weight = abs(wdpl(k)) / (sum(abs(wdpl(:))+eps1))
                if ((ww_dp_flux_deficit(1)*nxl + ww_dp_flux_deficit(2)*nyl + ww_dp_flux_deficit(3)*nzl) < 0.0) &
                    weight = abs(wdpr(k)) / (sum(abs(wdpr(:))+eps1))
                wdp_flux(1,k,iquad) = wdp_flux(1,k,iquad) + weight * ww_dp_flux_deficit(1)
                wdp_flux(2,k,iquad) = wdp_flux(2,k,iquad) + weight * ww_dp_flux_deficit(2)
                wdp_flux(3,k,iquad) = wdp_flux(3,k,iquad) + weight * ww_dp_flux_deficit(3)
            end do

            z_face = 0.0; p_face = 0.0
            z_edge_plus = 0.0; z_edge_minus = 0.0
            p_edge_plus = 0.0; p_edge_minus = 0.0

            ope_l = sqrt(btp%ope2_face_ave(1,iquad,iface))
            ope_r = sqrt(btp%ope2_face_ave(2,iquad,iface))
            p_face(1,1) = 0.0
            p_face(2,1) = 0.0
            !$acc loop seq
            do k = 1, nlayers_f
                p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,k)
                p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,k)
            end do

            one_plus_eta_edge = sqrt(btp%one_plus_eta_edge_2_ave(iquad,iface))
            z_face(1,nlayers_f+1)    = init%zbot_face(1,iquad,iface)
            z_face(2,nlayers_f+1)    = init%zbot_face(2,iquad,iface)
            z_edge_plus(nlayers_f+1) = init%zbot_face(1,iquad,iface)
            z_edge_minus(nlayers_f+1)= init%zbot_face(2,iquad,iface)
            !$acc loop seq
            do k = nlayers_f, 1, -1
                z_face(1,k)    = z_face(1,k+1)    + alpha_over_g(k) * (ope_l * ql(1,k))
                z_face(2,k)    = z_face(2,k+1)    + alpha_over_g(k) * (ope_r * qr(1,k))
                z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k) * (one_plus_eta_edge * ql(1,k))
                z_edge_minus(k)= z_edge_minus(k+1)+ alpha_over_g(k) * (one_plus_eta_edge * qr(1,k))
            end do

            p_edge_plus(2)  = one_plus_eta_edge * ql(1,1)
            p_edge_minus(2) = one_plus_eta_edge * qr(1,1)
            !$acc loop seq
            do k = 2, nlayers_f
                p_edge_plus(k+1)  = p_edge_plus(k)  + one_plus_eta_edge * ql(1,k)
                p_edge_minus(k+1) = p_edge_minus(k) + one_plus_eta_edge * qr(1,k)
            end do

            !$acc loop seq
            do k = 1, nlayers_f

                H_r_plus = 0.5*init%alpha_mlswe(k)*(p_edge_plus(k+1)**2 - p_edge_plus(k)**2)

                H_r_minus = 0.0
                !$acc loop seq
                do ktemp = 1, nlayers_f
                    z_intersect_top = min(z_edge_minus(ktemp),   z_edge_plus(k))
                    z_intersect_bot = max(z_edge_minus(ktemp+1), z_edge_plus(k+1))
                    dz_intersect    = z_intersect_top - z_intersect_bot
                    if (dz_intersect > 0.0) then
                        p_intersect_bot = p_edge_minus(ktemp+1) &
                                - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_minus(ktemp+1))
                        p_intersect_top = p_edge_minus(ktemp+1) &
                                - g_over_alpha(ktemp)*(z_intersect_top - z_edge_minus(ktemp+1))
                        H_r_minus = H_r_minus + &
                                0.5*init%alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
                    end if
                end do
                H_face(1,k,iquad) = 0.5*(H_r_plus + H_r_minus)

                H_r_minus = 0.5*init%alpha_mlswe(k)*(p_edge_minus(k+1)**2 - p_edge_minus(k)**2)

                H_r_plus = 0.0
                !$acc loop seq
                do ktemp = 1, nlayers_f
                    z_intersect_top = min(z_edge_plus(ktemp),   z_edge_minus(k))
                    z_intersect_bot = max(z_edge_plus(ktemp+1), z_edge_minus(k+1))
                    dz_intersect    = z_intersect_top - z_intersect_bot
                    if (dz_intersect > 0.0) then
                        p_intersect_bot = p_edge_plus(ktemp+1) &
                            - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_plus(ktemp+1))
                        p_intersect_top = p_edge_plus(ktemp+1) &
                            - g_over_alpha(ktemp)*(z_intersect_top - z_edge_plus(ktemp+1))
                        H_r_plus = H_r_plus &
                            + 0.5*init%alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
                    end if
                end do
                H_face(2,k,iquad) = 0.5*(H_r_plus + H_r_minus)
            end do

            if(er == -4) then
                p2l = 0.0; p2r = 0.0
                !$acc loop seq
                do k = 1, nlayers_f
                    p2l(k+1) = p_face(1,k+1)
                    H_face(1,k,iquad) = 0.5*init%alpha_mlswe(k)*(p2l(k+1)**2 - p2l(k)**2)
                    p2r(k+1) = p_face(2,k+1)
                    H_face(2,k,iquad) = 0.5*init%alpha_mlswe(k)*(p2r(k+1)**2 - p2r(k)**2)
                end do
            end if

            if(er > 0) then
                !$acc loop seq
                do k = 1, nlayers_f-1
                    p_inc1  = p_face(1,k+1) + g_over_alpha(k)*(z_face(1,k+1) - z_edge_plus(k+1))
                    H_corr1 = 0.5*init%alpha_mlswe(k)*(p_inc1**2 - p_face(1,k+1)**2)
                    H_face(1,k,iquad)   = H_face(1,k,iquad)   - H_corr1
                    H_face(1,k+1,iquad) = H_face(1,k+1,iquad) + H_corr1

                    p_inc2  = p_face(2,k+1) + g_over_alpha(k)*(z_face(2,k+1) - z_edge_minus(k+1))
                    H_corr2 = 0.5*init%alpha_mlswe(k)*(p_inc2**2 - p_face(2,k+1)**2)
                    H_face(2,k,iquad)   = H_face(2,k,iquad)   - H_corr2
                    H_face(2,k+1,iquad) = H_face(2,k+1,iquad) + H_corr2
                end do
            end if

            ! Rescale so the layer sum matches the BTP-averaged H, skipping
            ! layers whose H_face is already zero (dry) so a poorly
            ! conditioned weight can't corrupt an already-correct zero.
            acceleration = sum(H_face(1,:,iquad))
            if (abs(acceleration) > eps1) then
                weight = btp%H_face_ave(iquad,iface) / acceleration
                do k = 1, nlayers_f
                    if (H_face(1,k,iquad) /= 0.0) H_face(1,k,iquad) = H_face(1,k,iquad) * weight
                end do
            end if

            acceleration = sum(H_face(2,:,iquad))
            if (abs(acceleration) > eps1) then
                weight = btp%H_face_ave(iquad,iface) / acceleration
                do k = 1, nlayers_f
                    if (H_face(2,k,iquad) /= 0.0) H_face(2,k,iquad) = H_face(2,k,iquad) * weight
                end do
            end if

            wq = mf%jac_faceq(iquad,1,iface)

            !$acc loop seq
            do k = 1, nlayers_f

                hl = H_face(1,k,iquad)
                hr = H_face(2,k,iquad)

                flux_xl = nxl*(udp_flux(1,k,iquad) + hl) + nyl*udp_flux(2,k,iquad) + nzl*udp_flux(3,k,iquad)
                flux_xr = nxl*(udp_flux(1,k,iquad) + hr) + nyl*udp_flux(2,k,iquad) + nzl*udp_flux(3,k,iquad)
                flux_yl = nxl*vdp_flux(1,k,iquad) + nyl*(vdp_flux(2,k,iquad) + hl) + nzl*vdp_flux(3,k,iquad)
                flux_yr = nxl*vdp_flux(1,k,iquad) + nyl*(vdp_flux(2,k,iquad) + hr) + nzl*vdp_flux(3,k,iquad)
                flux_zl = nxl*wdp_flux(1,k,iquad) + nyl*wdp_flux(2,k,iquad) + nzl*(wdp_flux(3,k,iquad) + hl)
                flux_zr = nxl*wdp_flux(1,k,iquad) + nyl*wdp_flux(2,k,iquad) + nzl*(wdp_flux(3,k,iquad) + hr)

                !$acc loop seq
                do n = 1, ngl_f
                    hi = b%psiq(n,iquad)
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    !$acc atomic update
                    rhs(1,I,k) = rhs(1,I,k) - wq*hi*(nxl*dp_flux(1,k,iquad) + nyl*dp_flux(2,k,iquad) + nzl*dp_flux(3,k,iquad))
                    !$acc atomic update
                    rhs(2,I,k) = rhs(2,I,k) - wq*hi*flux_xl
                    !$acc atomic update
                    rhs(3,I,k) = rhs(3,I,k) - wq*hi*flux_yl
                    !$acc atomic update
                    rhs(4,I,k) = rhs(4,I,k) - wq*hi*flux_zl
                end do

                if(er > 0) then
                    !$acc loop seq
                    do n = 1, ngl_f
                        hi = b%psiq(n,iquad)
                        ir = mf%imapr(1,n,1,iface)
                        jr = mf%imapr(2,n,1,iface)
                        kr = mf%imapr(3,n,1,iface)
                        I  = G%intma(ir,jr,kr,er)
                        !$acc atomic update
                        rhs(1,I,k) = rhs(1,I,k) + wq*hi*(nxl*dp_flux(1,k,iquad) + nyl*dp_flux(2,k,iquad) + nzl*dp_flux(3,k,iquad))
                        !$acc atomic update
                        rhs(2,I,k) = rhs(2,I,k) + wq*hi*flux_xr
                        !$acc atomic update
                        rhs(3,I,k) = rhs(3,I,k) + wq*hi*flux_yr
                        !$acc atomic update
                        rhs(4,I,k) = rhs(4,I,k) + wq*hi*flux_zr
                    end do
                end if
            end do

            end do  ! iquad

        end do  ! iface
        !$acc end data

    end subroutine Apply_bcl_fluxes_sphere

    subroutine create_layers_volume_mass(G, inp, b, btp, bcl, init, tsp, dp_advec, qprime_df)

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(btp_CS),    intent(in)    :: btp
        type(bcl_CS),    intent(inout) :: bcl
        type(initial),   intent(in)    :: init
        type(tensor_CS), intent(in)    :: tsp

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: qprime_df
        real, dimension(G%npoin, inp%nlayers),    intent(out) :: dp_advec

        real :: wq, hi, dp_temp
        integer :: k, I, Iq, ip
        real, dimension(inp%nvar_bcl) :: qp, qb
        real, parameter :: eps = 1.0e-10
        real, dimension(inp%nlayers) :: dpp, udp, vdp
        real :: flux(2)

        dp_advec = 0.0
        bcl%sum_layer_mass_flux = 0.0

        do concurrent(Iq = 1:G%npoin_q)

            qb(1) = btp%ope_ave(Iq)
            qb(2) = btp%uvb_ave(1,Iq)
            qb(3) = btp%uvb_ave(2,Iq)
            wq    = tsp%wjac(Iq)

            do k = 1, inp%nlayers
                qp = 0.0
                do ip = 1, b%npts
                    I  = tsp%indexq(ip,Iq)
                    hi = tsp%psih(ip,Iq)
                    qp(:) = qp(:) + hi*qprime_df(:,I,k)
                end do
                dpp(k)  = qp(1) * qb(1)
                dp_temp = dpp(k)
                udp(k)  = (qp(2) + qb(2)) * dp_temp
                vdp(k)  = (qp(3) + qb(3)) * dp_temp
            end do

            do k = 1, inp%nlayers
                flux(1) = udp(k) + (dpp(k)/sum(dpp(:))) * (btp%btp_mass_flux_ave(1,Iq) - sum(udp(:)))
                flux(2) = vdp(k) + (dpp(k)/sum(dpp(:))) * (btp%btp_mass_flux_ave(2,Iq) - sum(vdp(:)))
                do ip = 1, b%npts
                    I = tsp%indexq(ip,Iq)
                    dp_advec(I,k) = dp_advec(I,k) + wq*(tsp%dpsidx(ip,Iq)*flux(1) + tsp%dpsidy(ip,Iq)*flux(2))
                end do
            end do
        end do

    end subroutine create_layers_volume_mass

    subroutine create_layer_mass_flux(G, inp, b, mf, btp, bcl, init, dp_advec, qprime_df)

        implicit none

        type(grid),    intent(in)    :: G
        type(input),   intent(in)    :: inp
        type(basis),   intent(in)    :: b
        type(face_CS), intent(in)    :: mf
        type(btp_CS),  intent(in)    :: btp
        type(bcl_CS),  intent(inout) :: bcl
        type(initial), intent(in)    :: init

        real, dimension(G%npoin, inp%nlayers),    intent(inout) :: dp_advec
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)    :: qprime_df

        integer :: k, iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n
        real :: wq, nxl, nyl, hi
        real, dimension(b%nq, inp%nlayers) :: flux_edge_u, flux_edge_v, flux_dp
        real :: dpl, dpr, uu, vv, flux, un, ul, ur, vl, vr, weight, dp_deficit(2)
        real :: dp_lr(2,inp%nlayers)
        real, dimension(inp%nvar_bcl) :: ql, qr, qbl, qbr
        real, parameter :: eps = 1.0e-10
        real :: flux_u, flux_v

        bcl%sum_layer_mass_flux_face = 0.0

        do concurrent(iface = 1:G%nface, iquad = 1:b%nq)

            if (G%face_type(iface) == 2) cycle

            el = G%face(7,iface)
            er = G%face(8,iface)

            qbl(1) = btp%ope_face_ave(1,iquad,iface)
            qbl(2) = btp%uvb_face_ave(1,1,iquad,iface)
            qbl(3) = btp%uvb_face_ave(2,1,iquad,iface)
            qbr(1) = btp%ope_face_ave(2,iquad,iface)
            qbr(2) = btp%uvb_face_ave(1,2,iquad,iface)
            qbr(3) = btp%uvb_face_ave(2,2,iquad,iface)

            nxl = mf%normal_vector_q(1,iquad,1,iface)
            nyl = mf%normal_vector_q(2,iquad,1,iface)

            do k = 1, inp%nlayers

                ql = 0.0; qr = 0.0
                do n = 1, b%ngl
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    hi = b%psiq(n,iquad)
                    ql(:) = ql(:) + hi*qprime_df(:,I,k)
                end do

                if (er > 0) then
                    do n = 1, b%ngl
                        ir = mf%imapr(1,n,1,iface)
                        jr = mf%imapr(2,n,1,iface)
                        kr = mf%imapr(3,n,1,iface)
                        I  = G%intma(ir,jr,kr,er)
                        hi = b%psiq(n,iquad)
                        qr(:) = qr(:) + hi*qprime_df(:,I,k)
                    end do
                else
                    qr(:) = ql(:)
                    if (er == -4) then
                        un = ql(2)*nxl + ql(3)*nyl
                        qr(2) = ql(2) - 2.0*un*nxl
                        qr(3) = ql(3) - 2.0*un*nyl
                    elseif (er == -2) then
                        qr(2) = -ql(2)
                        qr(3) = -ql(3)
                    end if
                end if

                dpl = qbl(1) * ql(1)
                dpr = qbr(1) * qr(1)

                dp_lr(1,k) = dpl
                dp_lr(2,k) = dpr

                ul = ql(2) + qbl(2)
                ur = qr(2) + qbr(2)
                vl = ql(3) + qbl(3)
                vr = qr(3) + qbr(3)

                uu = 0.5*(ul + ur)
                vv = 0.5*(vl + vr)

                un = uu*nxl + vv*nyl
                if(un > 0.0) then
                    flux_edge_u(iquad,k) = uu * dpl
                    flux_edge_v(iquad,k) = vv * dpl
                else
                    flux_edge_u(iquad,k) = uu * dpr
                    flux_edge_v(iquad,k) = vv * dpr
                end if

                flux_dp(iquad,k) = nxl*flux_edge_u(iquad,k) + nyl*flux_edge_v(iquad,k)

            end do

            dp_deficit(1) = btp%btp_mass_flux_face_ave(1,iquad,iface) - sum(flux_edge_u(iquad,:))
            dp_deficit(2) = btp%btp_mass_flux_face_ave(2,iquad,iface) - sum(flux_edge_v(iquad,:))

            wq = mf%jac_faceq(iquad,1,iface)

            do k = 1, inp%nlayers

                nxl = mf%normal_vector_q(1,iquad,1,iface)
                nyl = mf%normal_vector_q(2,iquad,1,iface)

                weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:))+eps))
                if ((dp_deficit(1)*nxl + dp_deficit(2)*nyl) < 0.0) &
                    weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:))+eps))
                flux_u = flux_edge_u(iquad,k) + weight*dp_deficit(1)
                flux_v = flux_edge_v(iquad,k) + weight*dp_deficit(2)

                flux = nxl*flux_u + nyl*flux_v

                do n = 1, b%ngl
                    hi = b%psiq(n,iquad)
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    dp_advec(I,k) = dp_advec(I,k) - wq*hi*flux
                end do

                if(er > 0) then
                    do n = 1, b%ngl
                        hi = b%psiq(n,iquad)
                        ir = mf%imapr(1,n,1,iface)
                        jr = mf%imapr(2,n,1,iface)
                        kr = mf%imapr(3,n,1,iface)
                        I  = G%intma(ir,jr,kr,er)
                        dp_advec(I,k) = dp_advec(I,k) + wq*hi*flux
                    end do
                end if

            end do
        end do

    end subroutine create_layer_mass_flux

    ! Sphere counterpart of create_layers_volume_mass: 4-wide (mass,u,v,w)
    ! interpolation and a sphere-corrected 3-component test-function gradient
    ! (dpsidx/y/z + curvature terms dpsidz_x/y/z), mirroring the mass-specific
    ! subset of create_rhs_dynamics_volume_bcl_sphere's layer loop 1 and its
    ! mass-flux consistency correction (mod_create_rhs_mlswe.F90, layer loop 1
    ! and the flux(*,1) block) -- no Coriolis/pressure-gradient/viscosity/APE,
    ! those are momentum-only.
    subroutine create_layers_volume_mass_sphere(G, inp, b, btp, bcl, init, tsp, dp_advec, qprime_df)

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(btp_CS),    intent(in)    :: btp
        type(bcl_CS),    intent(inout) :: bcl
        type(initial),   intent(in)    :: init
        type(tensor_CS), intent(in)    :: tsp

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: qprime_df
        real, dimension(G%npoin, inp%nlayers),    intent(out) :: dp_advec

        real :: wq, hi, dhdx, dhdy, dhdz
        integer :: k, I, Iq, ip
        real, dimension(4) :: qp, qb
        real, parameter :: eps1 = 1.0e-20
        real, dimension(inp%nlayers) :: dp, udp, vdp, wdp
        real :: flux(3), weight_dp
        logical :: is_dry

        dp_advec = 0.0

        do concurrent(Iq = 1:G%npoin_q)

            qb(1) = btp%ope_ave(Iq)
            qb(2) = btp%uvb_ave(1,Iq)
            qb(3) = btp%uvb_ave(2,Iq)
            qb(4) = btp%uvb_ave(3,Iq)
            wq    = tsp%wjac(Iq)

            do k = 1, inp%nlayers
                qp = 0.0
                do ip = 1, b%npts
                    I  = tsp%indexq(ip,Iq)
                    hi = tsp%psih(ip,Iq)
                    qp(1) = qp(1) + hi*qprime_df(1,I,k)
                    qp(2) = qp(2) + hi*qprime_df(2,I,k)
                    qp(3) = qp(3) + hi*qprime_df(3,I,k)
                    qp(4) = qp(4) + hi*qprime_df(4,I,k)
                end do

                ! Dry-cell protection, same convention as create_rhs_dynamics_volume_bcl_sphere.
                is_dry = (qp(1) * qb(1) < (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff)
                if (is_dry) qp(1) = (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff

                dp(k) = qp(1) * qb(1)

                if (is_dry) then
                    udp(k) = 0.0;  vdp(k) = 0.0;  wdp(k) = 0.0
                else
                    udp(k) = (qp(2) + qb(2)) * dp(k)
                    vdp(k) = (qp(3) + qb(3)) * dp(k)
                    wdp(k) = (qp(4) + qb(4)) * dp(k)
                end if
            end do

            do k = 1, inp%nlayers
                ! Mass flux: consistency correction distributes BTP mass
                ! surplus among wet layers only (same weighting as the
                ! combined sphere routine's flux(*,1) block).
                if (dp(k) > (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff) then
                    weight_dp = abs(dp(k)) / (sum(abs(dp(:))) + eps1)
                    flux(1) = udp(k) + weight_dp*(btp%btp_mass_flux_ave(1,Iq) - sum(udp(:)))
                    flux(2) = vdp(k) + weight_dp*(btp%btp_mass_flux_ave(2,Iq) - sum(vdp(:)))
                    flux(3) = wdp(k) + weight_dp*(btp%btp_mass_flux_ave(3,Iq) - sum(wdp(:)))
                else
                    flux(1) = udp(k)
                    flux(2) = vdp(k)
                    flux(3) = wdp(k)
                end if

                do ip = 1, b%npts
                    I = tsp%indexq(ip,Iq)
                    dhdx = tsp%dpsidx(ip,Iq) + tsp%dpsidz_x(ip,Iq)
                    dhdy = tsp%dpsidy(ip,Iq) + tsp%dpsidz_y(ip,Iq)
                    dhdz = tsp%dpsidz(ip,Iq) + tsp%dpsidz_z(ip,Iq)
                    dp_advec(I,k) = dp_advec(I,k) + wq*(dhdx*flux(1) + dhdy*flux(2) + dhdz*flux(3))
                end do
            end do
        end do

    end subroutine create_layers_volume_mass_sphere

    ! Sphere counterpart of create_layer_mass_flux: 4-wide interpolation,
    ! 3-component face normal, 3D wall reflection, mirroring the mass-specific
    ! subset of Apply_bcl_fluxes_sphere (ql/qr interpolation, dry-cell
    ! clamping, upwind dp_flux selection, dp_deficit consistency correction,
    ! final scatter) -- no momentum flux tensor (udp_flux/vdp_flux/wdp_flux)
    ! or its consistency correction, those are momentum-only.
    subroutine create_layer_mass_flux_sphere(G, inp, b, mf, btp, bcl, init, dp_advec, qprime_df)

        implicit none

        type(grid),    intent(in)    :: G
        type(input),   intent(in)    :: inp
        type(basis),   intent(in)    :: b
        type(face_CS), intent(in)    :: mf
        type(btp_CS),  intent(in)    :: btp
        type(bcl_CS),  intent(inout) :: bcl
        type(initial), intent(in)    :: init

        real, dimension(G%npoin, inp%nlayers),    intent(inout) :: dp_advec
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)    :: qprime_df

        integer :: k, iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n
        real :: wq, nxl, nyl, nzl, hi
        real, dimension(3,inp%nlayers,b%nq) :: dp_flux
        real :: dpl, dpr, uu, vv, ww, un, ul, ur, vl, vr, wl, wr
        real :: dp_lr(2,inp%nlayers)
        real, dimension(inp%nvar_bcl,inp%nlayers) :: ql, qr
        real, dimension(inp%nvar_bcl) :: qbl, qbr
        real, parameter :: eps1 = 1.0e-20
        real :: dp_deficit(3), weight
        real :: g_over_alpha(inp%nlayers)
        logical :: is_dry_l, is_dry_r
        integer :: ngl_f, nq_f, nlayers_f, nface_f

        ngl_f     = b%ngl
        nq_f      = b%nq
        nlayers_f = inp%nlayers
        nface_f   = G%nface

        do concurrent(iface = 1:nface_f, iquad = 1:nq_f)

            if (G%face_type(iface) == 2) cycle

            el = G%face(7,iface)
            er = G%face(8,iface)

            do k = 1, nlayers_f
                g_over_alpha(k) = gravity / init%alpha_mlswe(k)
            end do

            qbl(1) = btp%ope_face_ave(1,iquad,iface)
            qbl(2) = btp%uvb_face_ave(1,1,iquad,iface)
            qbl(3) = btp%uvb_face_ave(2,1,iquad,iface)
            qbl(4) = btp%uvb_face_ave(3,1,iquad,iface)
            qbr(1) = btp%ope_face_ave(2,iquad,iface)
            qbr(2) = btp%uvb_face_ave(1,2,iquad,iface)
            qbr(3) = btp%uvb_face_ave(2,2,iquad,iface)
            qbr(4) = btp%uvb_face_ave(3,2,iquad,iface)

            nxl = mf%normal_vector_q(1,iquad,1,iface)
            nyl = mf%normal_vector_q(2,iquad,1,iface)
            nzl = mf%normal_vector_q(3,iquad,1,iface)

            ql = 0.0; qr = 0.0
            do k = 1, nlayers_f

                do n = 1, ngl_f
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    hi = b%psiq(n,iquad)
                    ql(1,k) = ql(1,k) + hi*qprime_df(1,I,k)
                    ql(2,k) = ql(2,k) + hi*qprime_df(2,I,k)
                    ql(3,k) = ql(3,k) + hi*qprime_df(3,I,k)
                    ql(4,k) = ql(4,k) + hi*qprime_df(4,I,k)
                end do

                if (er > 0) then
                    do n = 1, ngl_f
                        ir = mf%imapr(1,n,1,iface)
                        jr = mf%imapr(2,n,1,iface)
                        kr = mf%imapr(3,n,1,iface)
                        I  = G%intma(ir,jr,kr,er)
                        hi = b%psiq(n,iquad)
                        qr(1,k) = qr(1,k) + hi*qprime_df(1,I,k)
                        qr(2,k) = qr(2,k) + hi*qprime_df(2,I,k)
                        qr(3,k) = qr(3,k) + hi*qprime_df(3,I,k)
                        qr(4,k) = qr(4,k) + hi*qprime_df(4,I,k)
                    end do
                else
                    qr(:,k) = ql(:,k)
                    if (er == -4) then
                        un = ql(2,k)*nxl + ql(3,k)*nyl + ql(4,k)*nzl
                        qr(2,k) = ql(2,k) - 2.0*un*nxl
                        qr(3,k) = ql(3,k) - 2.0*un*nyl
                        qr(4,k) = ql(4,k) - 2.0*un*nzl
                    elseif (er == -2) then
                        qr(2,k) = -ql(2,k)
                        qr(3,k) = -ql(3,k)
                        qr(4,k) = -ql(4,k)
                    end if
                end if

                is_dry_l = (ql(1,k) * qbl(1) < g_over_alpha(k) * inp%dry_cutoff)
                if (is_dry_l) then
                    ql(1,k) = g_over_alpha(k) * inp%dry_cutoff
                    ql(2,k) = 0.0;  ql(3,k) = 0.0;  ql(4,k) = 0.0
                end if
                is_dry_r = (qr(1,k) * qbr(1) < g_over_alpha(k) * inp%dry_cutoff)
                if (is_dry_r) then
                    qr(1,k) = g_over_alpha(k) * inp%dry_cutoff
                    qr(2,k) = 0.0;  qr(3,k) = 0.0;  qr(4,k) = 0.0
                end if

                dpl = qbl(1) * ql(1,k)
                dpr = qbr(1) * qr(1,k)
                ul  = ql(2,k) + qbl(2)
                ur  = qr(2,k) + qbr(2)
                vl  = ql(3,k) + qbl(3)
                vr  = qr(3,k) + qbr(3)
                wl  = ql(4,k) + qbl(4)
                wr  = qr(4,k) + qbr(4)

                dp_lr(1,k) = dpl
                dp_lr(2,k) = dpr

                uu = 0.5*(ul+ur)
                vv = 0.5*(vl+vr)
                ww = 0.5*(wl+wr)

                un = uu*nxl + vv*nyl + ww*nzl
                if (un > 0.0) then
                    dp_flux(1,k,iquad) = uu * dpl
                    dp_flux(2,k,iquad) = vv * dpl
                    dp_flux(3,k,iquad) = ww * dpl
                else
                    dp_flux(1,k,iquad) = uu * dpr
                    dp_flux(2,k,iquad) = vv * dpr
                    dp_flux(3,k,iquad) = ww * dpr
                end if

            end do

            dp_deficit(1) = btp%btp_mass_flux_face_ave(1,iquad,iface) - sum(dp_flux(1,:,iquad))
            dp_deficit(2) = btp%btp_mass_flux_face_ave(2,iquad,iface) - sum(dp_flux(2,:,iquad))
            dp_deficit(3) = btp%btp_mass_flux_face_ave(3,iquad,iface) - sum(dp_flux(3,:,iquad))

            wq = mf%jac_faceq(iquad,1,iface)

            do k = 1, nlayers_f

                weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:))+eps1))
                if ((dp_deficit(1)*nxl + dp_deficit(2)*nyl + dp_deficit(3)*nzl) < 0.0) &
                    weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:))+eps1))
                dp_flux(1,k,iquad) = dp_flux(1,k,iquad) + weight * dp_deficit(1)
                dp_flux(2,k,iquad) = dp_flux(2,k,iquad) + weight * dp_deficit(2)
                dp_flux(3,k,iquad) = dp_flux(3,k,iquad) + weight * dp_deficit(3)

                do n = 1, ngl_f
                    hi = b%psiq(n,iquad)
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    dp_advec(I,k) = dp_advec(I,k) - wq*hi*(nxl*dp_flux(1,k,iquad) + nyl*dp_flux(2,k,iquad) + nzl*dp_flux(3,k,iquad))
                end do

                if (er > 0) then
                    do n = 1, ngl_f
                        hi = b%psiq(n,iquad)
                        ir = mf%imapr(1,n,1,iface)
                        jr = mf%imapr(2,n,1,iface)
                        kr = mf%imapr(3,n,1,iface)
                        I  = G%intma(ir,jr,kr,er)
                        dp_advec(I,k) = dp_advec(I,k) + wq*hi*(nxl*dp_flux(1,k,iquad) + nyl*dp_flux(2,k,iquad) + nzl*dp_flux(3,k,iquad))
                    end do
                end if

            end do

        end do

    end subroutine create_layer_mass_flux_sphere

end module mod_create_rhs_mlswe
