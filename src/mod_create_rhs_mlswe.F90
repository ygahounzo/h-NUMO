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

    public :: layer_momentum_rhs, &
              interpolate_layer_from_quad_to_node, rhs_layer_shear_stress, layer_mass_rhs, &
              bcl_rhs

contains

    subroutine layer_momentum_rhs(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs_mom, qprime_df, q_df)

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

        real, dimension(2, G%npoin, inp%nlayers), intent(out) :: rhs_mom
        real, dimension(3, G%npoin, inp%nlayers), intent(in)  :: q_df, qprime_df

        call bcl_create_precommunicator(G, inp, b, mf, par, ref, mpic, qprime_df)
        call create_rhs_dynamics_volume_layers(G, inp, b, btp, init, tsp, rhs_mom, qprime_df, q_df)
        call Apply_layers_fluxes(G, inp, b, mf, btp, init, rhs_mom, qprime_df)
        call bcl_create_postcommunicator_momentum(G, inp, b, mf, par, btp, init, ref, mpic, rhs_mom)

    end subroutine layer_momentum_rhs

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
        real, dimension(3, G%npoin, inp%nlayers), intent(in)  :: qprime_df

        integer :: k

        call bcl_create_precommunicator(G, inp, b, mf, par, ref, mpic, qprime_df)
        call create_layers_volume_mass(G, inp, b, btp, bcl, init, tsp, dp_advec, qprime_df)
        call create_layer_mass_flux(G, inp, b, mf, btp, bcl, init, dp_advec, qprime_df)
        call bcl_create_postcommunicator_continuity(G, inp, b, mf, par, btp, ref, mpic, dp_advec)

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

        real, dimension(3, G%npoin, inp%nlayers), intent(out) :: rhs
        real, dimension(3, G%npoin, inp%nlayers), intent(in)  :: q_df, qprime_df

        !$acc data present(qprime_df, rhs) copyin(q_df)
        call bcl_create_precommunicator(G, inp, b, mf, par, ref, mpic, qprime_df)
        call create_rhs_dynamics_volume_bcl_qp(G, inp, b, btp, bcl, init, tsp, rhs, qprime_df, q_df)
        call Apply_bcl_fluxes(G, inp, b, mf, btp, bcl, init, rhs, qprime_df)
        call bcl_create_postcommunicator(G, inp, b, mf, par, btp, init, ref, mpic, rhs)
        !$acc end data
        ! rhs stays on device; create_rhs_bcl owns the download after the combining loop.

    end subroutine bcl_rhs

    subroutine interpolate_layer_from_quad_to_node(G, inp, b, tsp, mt, q_df, q)

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(tensor_CS), intent(in)    :: tsp
        type(metrics),   intent(in)    :: mt

        real, dimension(3, G%npoin,   inp%nlayers), intent(inout) :: q_df
        real, dimension(3, G%npoin_q, inp%nlayers), intent(in)    :: q

        real    :: wq, hi
        integer :: Iq, ip, k, I
        real    :: var_u(inp%nlayers), var_v(inp%nlayers)

        q_df(2:3,:,:) = 0.0

        do Iq = 1, G%npoin_q

            wq    = tsp%wjac(Iq)
            var_u = q(2,Iq,:)
            var_v = q(3,Iq,:)

            do ip = 1, b%npts
                I  = tsp%indexq(ip,Iq)
                hi = tsp%psih(ip,Iq)
                q_df(2,I,:) = q_df(2,I,:) + wq*var_u(:)*hi
                q_df(3,I,:) = q_df(3,I,:) + wq*var_v(:)*hi
            end do
        end do

        do k = 1, inp%nlayers
            q_df(2,:,k) = mt%massinv(:) * q_df(2,:,k)
            q_df(3,:,k) = mt%massinv(:) * q_df(3,:,k)
        end do

    end subroutine interpolate_layer_from_quad_to_node

    subroutine rhs_layer_shear_stress(G, inp, b, init, tsp, rhs_stress, q_df)

        implicit none

        type(grid),      intent(in) :: G
        type(input),     intent(in) :: inp
        type(basis),     intent(in) :: b
        type(initial),   intent(in) :: init
        type(tensor_CS), intent(in) :: tsp

        real, intent(in)  :: q_df(3, G%npoin, inp%nlayers)
        real, intent(out) :: rhs_stress(2, G%npoin, inp%nlayers)

        real, dimension(inp%nlayers+1) :: tau_u, tau_v
        real, dimension(inp%nlayers)   :: a, bc, c, dp, udp, vdp
        real, dimension(2,inp%nlayers) :: r, uv
        real    :: wq, hi, tau_u_q, tau_v_q, coeff, mult, coeff1
        integer :: k, Iq, I, ip

        rhs_stress = 0.0

        do Iq = 1, G%npoin_q

            dp = 0.0; udp = 0.0; vdp = 0.0
            do ip = 1, b%npts
                I  = tsp%indexq(ip,Iq)
                hi = tsp%psih(ip,Iq)
                dp(:)  = dp(:)  + hi * q_df(1,I,:)
                udp(:) = udp(:) + hi * q_df(2,I,:)
                vdp(:) = vdp(:) + hi * q_df(3,I,:)
            end do

            coeff  = max(sqrt(0.5*init%coriolis_quad(Iq)*inp%ad_mlswe)/init%alpha_mlswe(1), &
                         inp%ad_mlswe/(init%alpha_mlswe(1) * inp%max_shear_dz))
            coeff1 = gravity * inp%dt * coeff

            do k = 1, inp%nlayers
                a(k)   = -coeff
                bc(k)  = dp(k) + 2.0*coeff1
                c(k)   = -coeff1
                r(1,k) = udp(k)/dp(k)
                r(2,k) = vdp(k)/dp(k)
            end do

            bc(1)           = dp(1) + coeff1
            bc(inp%nlayers) = dp(inp%nlayers) + coeff1
            a(1)            = 0.0
            c(inp%nlayers)  = 0.0

            do k = 2, inp%nlayers
                mult   = a(k) / bc(k-1)
                bc(k)  = bc(k) - mult*c(k-1)
                r(1,k) = r(1,k) - mult*r(1,k-1)
                r(2,k) = r(2,k) - mult*r(2,k-1)
            end do

            r(1,inp%nlayers)  = r(1,inp%nlayers)  / bc(inp%nlayers)
            r(2,inp%nlayers)  = r(2,inp%nlayers)  / bc(inp%nlayers)
            uv(1,inp%nlayers) = r(1,inp%nlayers)
            uv(2,inp%nlayers) = r(2,inp%nlayers)

            do k = inp%nlayers-1, 1, -1
                r(1,k)  = (r(1,k) - c(k)*r(1,k+1)) / bc(k)
                r(2,k)  = (r(2,k) - c(k)*r(2,k+1)) / bc(k)
                uv(1,k) = r(1,k)
                uv(2,k) = r(2,k)
            end do

            tau_u(1) = 0.0; tau_v(1) = 0.0
            do k = 2, inp%nlayers
                tau_u(k) = coeff*(uv(1,k-1) - uv(1,k))
                tau_v(k) = coeff*(uv(2,k-1) - uv(2,k))
            end do

            wq = tsp%wjac(Iq)

            do k = 1, inp%nlayers
                tau_u_q = gravity*(tau_u(k) - tau_u(k+1))
                tau_v_q = gravity*(tau_v(k) - tau_v(k+1))
                do ip = 1, b%npts
                    I  = tsp%indexq(ip,Iq)
                    hi = tsp%psih(ip,Iq)
                    rhs_stress(1,I,k) = rhs_stress(1,I,k) + wq*hi*tau_u_q
                    rhs_stress(2,I,k) = rhs_stress(2,I,k) + wq*hi*tau_v_q
                end do
            end do

        end do

    end subroutine rhs_layer_shear_stress

    subroutine create_rhs_dynamics_volume_layers(G, inp, b, btp, init, tsp, rhs_mom, qprime_df, q_df)

        implicit none

        type(grid),      intent(in) :: G
        type(input),     intent(in) :: inp
        type(basis),     intent(in) :: b
        type(btp_CS),    intent(in) :: btp
        type(initial),   intent(in) :: init
        type(tensor_CS), intent(in) :: tsp

        real, dimension(2, G%npoin, inp%nlayers), intent(out) :: rhs_mom
        real, dimension(3, G%npoin, inp%nlayers), intent(in)  :: q_df, qprime_df

        real :: wq, hi, dhdx, dhdy, tau_wind_u, tau_wind_v, temp1
        real :: Hq, source_x, source_y, Pstress, Pbstress
        integer :: k, I, Iq, ip
        real, dimension(inp%nlayers+1) :: pprime_temp, z
        real :: tempbot, weight, acceleration, pbq
        real, dimension(3) :: qp, qb
        real, dimension(inp%nlayers) :: temp_uu, temp_vv, H_tmp, u_udp, v_vdp
        real, dimension(2,inp%nlayers) :: u_vdp
        real :: p_tmp(inp%nlayers+1), u, v, dp, weightq, one_over_sumuq, one_over_sumvq
        real :: uu_dp_deficitq, uv_dp_deficitq, vv_dp_deficitq, gradz(2,inp%nlayers+1)
        real, parameter :: eps1 = 1.0e-20
        real :: flux(2,2)

        rhs_mom = 0.0
        pprime_temp = 0.0

        Pstress  = (gravity/init%alpha_mlswe(1))            * 50.0
        Pbstress = (gravity/init%alpha_mlswe(inp%nlayers))  * 10.0

        do concurrent(Iq = 1:G%npoin_q)

            p_tmp(:)       = 0.0
            pprime_temp(:) = 0.0
            temp_uu = 0.0; temp_vv = 0.0

            do k = 1, inp%nlayers
                qp = 0.0; qb = 0.0
                do ip = 1, b%npts
                    I  = tsp%indexq(ip,Iq)
                    hi = tsp%psih(ip,Iq)
                    qp(:) = qp(:) + hi*qprime_df(:,I,k)
                end do

                qb(1) = btp%ope_ave(Iq)
                qb(2) = btp%uvb_ave(1,Iq)
                qb(3) = btp%uvb_ave(2,Iq)

                p_tmp(k+1) = p_tmp(k) + sqrt(btp%ope2_ave(Iq)) * qp(1)
                H_tmp(k)   = 0.5*init%alpha_mlswe(k) * (p_tmp(k+1)**2 - p_tmp(k)**2)

                dp = qp(1) * qb(1)
                u  = qp(2) + qb(2)
                v  = qp(3) + qb(3)

                u_udp(k)   = dp * u*u
                v_vdp(k)   = dp * v*v
                u_vdp(1,k) = u * v * dp
                u_vdp(2,k) = v * u * dp

                temp_uu(k) = abs(dp*u) + eps1
                temp_vv(k) = abs(dp*v) + eps1

                pprime_temp(k+1) = pprime_temp(k) + qp(1)
            end do

            gradz(:,:) = 0.0; pbq = 0.0
            do ip = 1, b%npts
                I = tsp%indexq(ip,Iq)
                z(inp%nlayers+1) = init%zbot_df(I)
                do k = inp%nlayers, 1, -1
                    z(k) = z(k+1) + (init%alpha_mlswe(k)/gravity) * &
                           (sqrt(btp%ope2_ave_df(I))*qprime_df(1,I,k))
                    gradz(1,k) = gradz(1,k) + tsp%dpsidx(ip,Iq) * z(k)
                    gradz(2,k) = gradz(2,k) + tsp%dpsidy(ip,Iq) * z(k)
                end do
                gradz(1,inp%nlayers+1) = init%grad_zbot_quad(1,Iq)
                gradz(2,inp%nlayers+1) = init%grad_zbot_quad(2,Iq)
                pbq = pbq + tsp%psih(ip,Iq) * init%pbprime_df(I)
            end do

            uu_dp_deficitq = btp%Qu_ave(Iq)  - sum(u_udp(:))
            uv_dp_deficitq = btp%Quv_ave(Iq) - sum(u_vdp(1,:))
            vv_dp_deficitq = btp%Qv_ave(Iq)  - sum(v_vdp(:))

            one_over_sumuq = 1.0/sum(temp_uu(:))
            one_over_sumvq = 1.0/sum(temp_vv(:))

            wq = tsp%wjac(Iq)

            do k = 1, inp%nlayers

                weightq    = temp_uu(k) * one_over_sumuq
                flux(1,1)  = u_udp(k)   + weightq * uu_dp_deficitq
                flux(2,1)  = u_vdp(1,k) + weightq * uv_dp_deficitq

                weightq    = temp_vv(k) * one_over_sumvq
                flux(1,2)  = u_vdp(2,k) + weightq * uv_dp_deficitq
                flux(2,2)  = v_vdp(k)   + weightq * vv_dp_deficitq

                Hq = H_tmp(k)
                weight = 1.0
                acceleration = sum(H_tmp(:))
                if(acceleration > 0.0) weight = btp%H_ave(Iq) / acceleration
                Hq = Hq * weight

                flux(1,1) = flux(1,1) + Hq
                flux(2,2) = flux(2,2) + Hq

                temp1      = (min(pprime_temp(k+1), Pstress) - min(pprime_temp(k), Pstress)) / Pstress
                tau_wind_u = temp1 * init%tau_wind(1,Iq)
                tau_wind_v = temp1 * init%tau_wind(2,Iq)

                tempbot = min(Pbstress, pbq-pprime_temp(k+1)) - min(Pbstress, pbq-pprime_temp(k))
                tempbot = tempbot / Pbstress

                source_x = gravity*(tau_wind_u - tempbot*btp%tau_bot_ave(1,Iq) + &
                            p_tmp(k)*gradz(1,k) - p_tmp(k+1)*gradz(1,k+1))
                source_y = gravity*(tau_wind_v - tempbot*btp%tau_bot_ave(2,Iq) + &
                            p_tmp(k)*gradz(2,k) - p_tmp(k+1)*gradz(2,k+1))

                do ip = 1, b%npts
                    I    = tsp%indexq(ip,Iq)
                    hi   = tsp%psih(ip,Iq)
                    dhdx = tsp%dpsidx(ip,Iq)
                    dhdy = tsp%dpsidy(ip,Iq)
                    rhs_mom(1,I,k) = rhs_mom(1,I,k) + wq*(hi*source_x + dhdx*flux(1,1) + dhdy*flux(1,2))
                    rhs_mom(2,I,k) = rhs_mom(2,I,k) + wq*(hi*source_y + dhdx*flux(2,1) + dhdy*flux(2,2))
                end do
            end do
        end do

    end subroutine create_rhs_dynamics_volume_layers

    subroutine create_rhs_dynamics_volume_bcl(G, inp, b, btp, bcl, init, tsp, rhs, qprime_df, q_df)
       !===========================================================================
       !  Volume contribution to the BCL RHS (DG weak form, interior).
       !
       !  One element at a time: accumulates into a local rhs_loc(3,npts,nlayers)
       !  then scatters to global rhs once with no race conditions.
       !===========================================================================

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
        integer :: k, I, Iq, ip, ie, ijq
        real, dimension(inp%nlayers+1) :: pprime_temp, z
        real :: tempbot, weight, acceleration, pbq
        real, dimension(3) :: qp, qb
        real, dimension(inp%nlayers) :: temp_uu, temp_vv, H_tmp, u_udp, v_vdp, udp, vdp, dp
        real, dimension(2,inp%nlayers) :: u_vdp
        real :: p_tmp(inp%nlayers+1), u, v, weightq, one_over_sumuq, one_over_sumvq
        real :: uu_dp_deficitq, uv_dp_deficitq, vv_dp_deficitq, gradz(2,inp%nlayers+1)
        real, parameter :: eps1 = 1.0e-20
        real :: flux(3,3)
        real, dimension(inp%nlayers)   :: a_visc, bc_visc, c_visc
        real, dimension(2,inp%nlayers) :: r_visc, uv_visc
        real, dimension(inp%nlayers+1) :: tau_u, tau_v
        real :: coeff, coeff1, mult
        real :: rhs_loc(3, b%npts, inp%nlayers)
        integer :: npts_l, npts_q, nelem_l, nlayers_l
        real    :: Pstress_l, Pbstress_l, ad_mlswe_l, max_shear_dz_l, dt_l
        logical :: is_rk3_or_lsrk3, do_ad_visc

        npts_l       = b%npts
        npts_q       = b%npts_quad
        nelem_l      = G%nelem
        nlayers_l    = inp%nlayers
        Pstress_l    = (gravity/init%alpha_mlswe(1))           * 50.0
        Pbstress_l   = (gravity/init%alpha_mlswe(inp%nlayers)) * 10.0
        ad_mlswe_l   = real(inp%ad_mlswe)
        max_shear_dz_l = inp%max_shear_dz
        dt_l           = inp%dt
        is_rk3_or_lsrk3 = (trim(inp%bcl_time_method) == 'rk3' .or. &
                       trim(inp%bcl_time_method) == 'lsrk3')
        do_ad_visc  = (inp%ad_mlswe > 0 .and. is_rk3_or_lsrk3)

        Pstress  = Pstress_l
        Pbstress = Pbstress_l

        bcl%sum_layer_mass_flux = 0.0

        !$acc data present(tsp%indexq, tsp%indexq_e, tsp%psih, tsp%dpsidx, tsp%dpsidy, tsp%wjac,    &
        !$acc              btp%ope_ave, btp%ope2_ave, btp%uvb_ave, btp%btp_mass_flux_ave,             &
        !$acc              btp%H_ave, btp%Qu_ave, btp%Quv_ave, btp%Qv_ave, btp%tau_bot_ave,          &
        !$acc              btp%ope2_ave_df,                                                            &
        !$acc              init%alpha_mlswe, init%zbot_df, init%grad_zbot_quad,                       &
        !$acc              init%pbprime_df, init%tau_wind, init%coriolis_quad,                        &
        !$acc              qprime_df, q_df, rhs)

        !$acc kernels
        rhs = 0.0
        !$acc end kernels

        !$acc parallel loop gang                                                                       &
        !$acc   private(rhs_loc, p_tmp, pprime_temp, z, gradz,                                        &
        !$acc           qp, qb,                                                                        &
        !$acc           temp_uu, temp_vv, H_tmp, u_udp, v_vdp, udp, vdp, dp, u_vdp,                  &
        !$acc           flux, a_visc, bc_visc, c_visc, r_visc, uv_visc, tau_u, tau_v,                 &
        !$acc           wq, hi, dhdx, dhdy, pbq, Hq, I, Iq, ijq, ip, k,                              &
        !$acc           source_x, source_y, temp1, tau_wind_u, tau_wind_v,                             &
        !$acc           tempbot, weight, acceleration, weightq,                                         &
        !$acc           one_over_sumuq, one_over_sumvq,                                                &
        !$acc           uu_dp_deficitq, uv_dp_deficitq, vv_dp_deficitq,                               &
        !$acc           u, v, coeff, coeff1, mult)
        do ie = 1, nelem_l

            !$acc loop seq
            do k = 1, nlayers_l
                do ip = 1, npts_l
                    rhs_loc(1,ip,k) = 0.0
                    rhs_loc(2,ip,k) = 0.0
                    rhs_loc(3,ip,k) = 0.0
                end do
            end do

            !$acc loop seq
            do ijq = 1, npts_q

                Iq = tsp%indexq_e(ijq, ie)
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

                if (do_ad_visc) then
                    coeff  = max(sqrt(0.5*init%coriolis_quad(Iq)*ad_mlswe_l)/init%alpha_mlswe(1), &
                                 ad_mlswe_l/(init%alpha_mlswe(1) * max_shear_dz_l))
                    coeff1 = gravity * dt_l * coeff

                    !$acc loop seq
                    do k = 1, nlayers_l
                        a_visc(k)   = -coeff
                        bc_visc(k)  = dp(k) + 2.0*coeff1
                        c_visc(k)   = -coeff1
                        r_visc(1,k) = udp(k)/dp(k)
                        r_visc(2,k) = vdp(k)/dp(k)
                    end do

                    bc_visc(1)          = dp(1) + coeff1
                    bc_visc(nlayers_l)  = dp(nlayers_l) + coeff1
                    a_visc(1)           = 0.0
                    c_visc(nlayers_l)   = 0.0

                    !$acc loop seq
                    do k = 2, nlayers_l
                        mult        = a_visc(k) / bc_visc(k-1)
                        bc_visc(k)  = bc_visc(k) - mult*c_visc(k-1)
                        r_visc(1,k) = r_visc(1,k) - mult*r_visc(1,k-1)
                        r_visc(2,k) = r_visc(2,k) - mult*r_visc(2,k-1)
                    end do

                    uv_visc(1,nlayers_l) = r_visc(1,nlayers_l) / bc_visc(nlayers_l)
                    uv_visc(2,nlayers_l) = r_visc(2,nlayers_l) / bc_visc(nlayers_l)
                    !$acc loop seq
                    do k = nlayers_l-1, 1, -1
                        uv_visc(1,k) = (r_visc(1,k) - c_visc(k)*uv_visc(1,k+1)) / bc_visc(k)
                        uv_visc(2,k) = (r_visc(2,k) - c_visc(k)*uv_visc(2,k+1)) / bc_visc(k)
                    end do

                    tau_u(1) = 0.0;  tau_v(1) = 0.0
                    !$acc loop seq
                    do k = 2, nlayers_l
                        tau_u(k) = coeff*(uv_visc(1,k-1) - uv_visc(1,k))
                        tau_v(k) = coeff*(uv_visc(2,k-1) - uv_visc(2,k))
                    end do
                    tau_u(nlayers_l+1) = 0.0;  tau_v(nlayers_l+1) = 0.0
                end if

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

                uu_dp_deficitq = btp%Qu_ave(Iq)  - sum(u_udp(:))
                uv_dp_deficitq = btp%Quv_ave(Iq) - sum(u_vdp(1,:))
                vv_dp_deficitq = btp%Qv_ave(Iq)  - sum(v_vdp(:))

                one_over_sumuq = 1.0/sum(temp_uu(:))
                one_over_sumvq = 1.0/sum(temp_vv(:))

                wq = tsp%wjac(Iq)

                !$acc loop seq
                do k = 1, nlayers_l

                    weightq    = temp_uu(k) * one_over_sumuq
                    u_udp(k)   = u_udp(k)   + weightq * uu_dp_deficitq
                    u_vdp(1,k) = u_vdp(1,k) + weightq * uv_dp_deficitq

                    weightq    = temp_vv(k) * one_over_sumvq
                    u_vdp(2,k) = u_vdp(2,k) + weightq * uv_dp_deficitq
                    v_vdp(k)   = v_vdp(k)   + weightq * vv_dp_deficitq

                    Hq = H_tmp(k)
                    weight = 1.0
                    acceleration = sum(H_tmp(:))
                    if(acceleration > 0.0) weight = btp%H_ave(Iq) / acceleration
                    Hq = Hq * weight

                    flux(1,1) = udp(k) + (dp(k)/sum(dp(:))) * (btp%btp_mass_flux_ave(1,Iq) - sum(udp(:)))
                    flux(2,1) = vdp(k) + (dp(k)/sum(dp(:))) * (btp%btp_mass_flux_ave(2,Iq) - sum(vdp(:)))

                    flux(1,2) = u_udp(k) + Hq
                    flux(2,2) = u_vdp(1,k)
                    flux(1,3) = u_vdp(2,k)
                    flux(2,3) = v_vdp(k) + Hq

                    temp1      = (min(pprime_temp(k+1), Pstress_l) - min(pprime_temp(k), Pstress_l)) / Pstress_l
                    tau_wind_u = temp1 * init%tau_wind(1,Iq)
                    tau_wind_v = temp1 * init%tau_wind(2,Iq)

                    tempbot = min(Pbstress_l, pbq-pprime_temp(k+1)) - min(Pbstress_l, pbq-pprime_temp(k))
                    tempbot = tempbot / Pbstress_l

                    source_x = gravity*(tau_wind_u - tempbot*btp%tau_bot_ave(1,Iq) + &
                                p_tmp(k)*gradz(1,k) - p_tmp(k+1)*gradz(1,k+1))
                    source_y = gravity*(tau_wind_v - tempbot*btp%tau_bot_ave(2,Iq) + &
                                p_tmp(k)*gradz(2,k) - p_tmp(k+1)*gradz(2,k+1))

                    if (is_rk3_or_lsrk3) then
                        source_x = source_x + init%coriolis_quad(Iq) * vdp(k)
                        source_y = source_y - init%coriolis_quad(Iq) * udp(k)
                        if (do_ad_visc) then
                            source_x = source_x + gravity*(tau_u(k) - tau_u(k+1))
                            source_y = source_y + gravity*(tau_v(k) - tau_v(k+1))
                        end if
                    end if

                    !$acc loop seq
                    do ip = 1, npts_l
                        I    = tsp%indexq(ip,Iq)
                        hi   = tsp%psih(ip,Iq)
                        dhdx = tsp%dpsidx(ip,Iq)
                        dhdy = tsp%dpsidy(ip,Iq)
                        rhs_loc(1,ip,k) = rhs_loc(1,ip,k) + wq*(dhdx*flux(1,1) + dhdy*flux(2,1))
                        rhs_loc(2,ip,k) = rhs_loc(2,ip,k) + wq*(hi*source_x + dhdx*flux(1,2) + dhdy*flux(2,2))
                        rhs_loc(3,ip,k) = rhs_loc(3,ip,k) + wq*(hi*source_y + dhdx*flux(1,3) + dhdy*flux(2,3))
                    end do
                end do

            end do  ! ijq

            ! Scatter gang-private rhs_loc to global rhs (no atomics: DG DOFs are element-exclusive).
            Iq = tsp%indexq_e(1, ie)
            !$acc loop seq
            do k = 1, nlayers_l
                !$acc loop seq
                do ip = 1, npts_l
                    I = tsp%indexq(ip, Iq)
                    rhs(1,I,k) = rhs_loc(1,ip,k)
                    rhs(2,I,k) = rhs_loc(2,ip,k)
                    rhs(3,I,k) = rhs_loc(3,ip,k)
                end do
            end do

        end do  ! ie
        !$acc end data

    end subroutine create_rhs_dynamics_volume_bcl

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
        real :: uu_dp_deficitq, uv_dp_deficitq, vv_dp_deficitq, gradz(2,inp%nlayers+1)
        real, parameter :: eps1 = 1.0e-20
        real :: flux(3,3)
        real, dimension(inp%nlayers)   :: a_visc, bc_visc, c_visc
        real, dimension(2,inp%nlayers) :: r_visc, uv_visc
        real, dimension(inp%nlayers+1) :: tau_u, tau_v
        real :: coeff, coeff1, mult

        ! Scalar captures — inp/b/G are not on device; string comparison must stay on host.
        integer :: npoin_q_l, npts_l, nlayers_l
        real    :: ad_mlswe_l, max_shear_dz_l, dt_l
        logical :: is_rk3_or_lsrk3, do_ad_visc

        npoin_q_l      = G%npoin_q
        npts_l         = b%npts
        nlayers_l      = inp%nlayers
        ad_mlswe_l     = real(inp%ad_mlswe)
        max_shear_dz_l = real(inp%max_shear_dz)
        dt_l           = inp%dt
        is_rk3_or_lsrk3 = (trim(inp%bcl_time_method) == 'rk3' .or. &
                            trim(inp%bcl_time_method) == 'lsrk3')
        do_ad_visc  = (inp%ad_mlswe > 0 .and. is_rk3_or_lsrk3)

        Pstress  = (gravity/init%alpha_mlswe(1))           * 50.0
        Pbstress = (gravity/init%alpha_mlswe(inp%nlayers)) * 10.0
        bcl%sum_layer_mass_flux = 0.0

        !$acc data present(tsp%indexq, tsp%psih, tsp%dpsidx, tsp%dpsidy, tsp%wjac,   &
        !$acc              qprime_df, rhs,                                              &
        !$acc              init%alpha_mlswe, init%zbot_df, init%pbprime_df,            &
        !$acc              init%grad_zbot_quad, init%tau_wind, init%coriolis_quad,     &
        !$acc              btp%ope_ave, btp%uvb_ave, btp%ope2_ave, btp%ope2_ave_df,   &
        !$acc              btp%H_ave, btp%btp_mass_flux_ave, btp%Qu_ave,               &
        !$acc              btp%Quv_ave, btp%Qv_ave, btp%tau_bot_ave)

        !$acc kernels present(rhs)
        rhs = 0.0
        !$acc end kernels

        !$acc parallel loop gang                                                        &
        !$acc   private(pprime_temp, z, qp, qb, temp_uu, temp_vv, H_tmp,              &
        !$acc           u_udp, v_vdp, u_vdp, udp, vdp, dp, p_tmp, gradz, flux,        &
        !$acc           a_visc, bc_visc, c_visc, r_visc, uv_visc, tau_u, tau_v,       &
        !$acc           wq, hi, dhdx, dhdy, tau_wind_u, tau_wind_v, temp1,             &
        !$acc           Hq, source_x, source_y, pbq, tempbot, weight, acceleration,   &
        !$acc           u, v, weightq, one_over_sumuq, one_over_sumvq,                 &
        !$acc           uu_dp_deficitq, uv_dp_deficitq, vv_dp_deficitq,               &
        !$acc           coeff, coeff1, mult, k, I, ip)                                 &
        !$acc   firstprivate(Pstress, Pbstress, npoin_q_l, npts_l, nlayers_l,         &
        !$acc                ad_mlswe_l, max_shear_dz_l, dt_l,                        &
        !$acc                is_rk3_or_lsrk3, do_ad_visc)
        do Iq = 1, npoin_q_l

            p_tmp(1)       = 0.0
            pprime_temp(:) = 0.0
            temp_uu = 0.0; temp_vv = 0.0
            tau_u   = 0.0; tau_v   = 0.0

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

            if (do_ad_visc) then
                coeff  = max(sqrt(0.5*init%coriolis_quad(Iq)*ad_mlswe_l)/init%alpha_mlswe(1), &
                             ad_mlswe_l/(init%alpha_mlswe(1) * max_shear_dz_l))
                coeff1 = gravity * dt_l * coeff

                !$acc loop seq
                do k = 1, nlayers_l
                    a_visc(k)   = -coeff
                    bc_visc(k)  = dp(k) + 2.0*coeff1
                    c_visc(k)   = -coeff1
                    r_visc(1,k) = udp(k)/dp(k)
                    r_visc(2,k) = vdp(k)/dp(k)
                end do

                bc_visc(1)         = dp(1) + coeff1
                bc_visc(nlayers_l) = dp(nlayers_l) + coeff1
                a_visc(1)          = 0.0
                c_visc(nlayers_l)  = 0.0

                !$acc loop seq
                do k = 2, nlayers_l
                    mult        = a_visc(k) / bc_visc(k-1)
                    bc_visc(k)  = bc_visc(k) - mult*c_visc(k-1)
                    r_visc(1,k) = r_visc(1,k) - mult*r_visc(1,k-1)
                    r_visc(2,k) = r_visc(2,k) - mult*r_visc(2,k-1)
                end do

                uv_visc(1,nlayers_l) = r_visc(1,nlayers_l) / bc_visc(nlayers_l)
                uv_visc(2,nlayers_l) = r_visc(2,nlayers_l) / bc_visc(nlayers_l)
                !$acc loop seq
                do k = nlayers_l-1, 1, -1
                    uv_visc(1,k) = (r_visc(1,k) - c_visc(k)*uv_visc(1,k+1)) / bc_visc(k)
                    uv_visc(2,k) = (r_visc(2,k) - c_visc(k)*uv_visc(2,k+1)) / bc_visc(k)
                end do

                !$acc loop seq
                do k = 2, nlayers_l
                    tau_u(k) = coeff*(uv_visc(1,k-1) - uv_visc(1,k))
                    tau_v(k) = coeff*(uv_visc(2,k-1) - uv_visc(2,k))
                end do
            end if

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

            uu_dp_deficitq = btp%Qu_ave(Iq)  - sum(u_udp(:))
            uv_dp_deficitq = btp%Quv_ave(Iq) - sum(u_vdp(1,:))
            vv_dp_deficitq = btp%Qv_ave(Iq)  - sum(v_vdp(:))

            one_over_sumuq = 1.0/sum(temp_uu(:))
            one_over_sumvq = 1.0/sum(temp_vv(:))

            wq = tsp%wjac(Iq)

            !$acc loop seq
            do k = 1, nlayers_l

                weightq    = temp_uu(k) * one_over_sumuq
                u_udp(k)   = u_udp(k)   + weightq * uu_dp_deficitq
                u_vdp(1,k) = u_vdp(1,k) + weightq * uv_dp_deficitq

                weightq    = temp_vv(k) * one_over_sumvq
                u_vdp(2,k) = u_vdp(2,k) + weightq * uv_dp_deficitq
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

                if (is_rk3_or_lsrk3) then
                    source_x = source_x + init%coriolis_quad(Iq) * vdp(k)
                    source_y = source_y - init%coriolis_quad(Iq) * udp(k)

                    if (do_ad_visc) then
                        source_x = source_x + gravity*(tau_u(k) - tau_u(k+1))
                        source_y = source_y + gravity*(tau_v(k) - tau_v(k+1))
                    end if
                end if

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

    subroutine Apply_layers_fluxes(G, inp, b, mf, btp, init, rhs_mom, qprime_df)

        implicit none

        type(grid),    intent(in)    :: G
        type(input),   intent(in)    :: inp
        type(basis),   intent(in)    :: b
        type(face_CS), intent(in)    :: mf
        type(btp_CS),  intent(in)    :: btp
        type(initial), intent(in)    :: init

        real, dimension(2, G%npoin, inp%nlayers), intent(inout) :: rhs_mom
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
        real :: ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv, un
        real, dimension(inp%nlayers) :: udpl, udpr, vdpl, vdpr
        real :: uu_dp_flux_deficit(2), vv_dp_flux_deficit(2)
        real, parameter :: eps1 = 1.0e-20
        integer :: il, jl, ir, jr, kl, kr, n
        real :: wq, hi, flux_xl, flux_xr, flux_yl, flux_yr, hl, hr
        real, dimension(2,inp%nlayers,b%nq) :: H_face, udp_flux, vdp_flux

        do k = 1, inp%nlayers
            alpha_over_g(k) = init%alpha_mlswe(k) / gravity
            g_over_alpha(k) = gravity / init%alpha_mlswe(k)
        end do

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

            ql = 0.0; qr = 0.0

            do k = 1, inp%nlayers

                do n = 1, b%ngl
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    hi = b%psiq(n,iquad)
                    ql(:,k) = ql(:,k) + hi*qprime_df(:,I,k)
                end do

                if (er > 0) then
                    do n = 1, b%ngl
                        ir = mf%imapr(1,n,1,iface)
                        jr = mf%imapr(2,n,1,iface)
                        kr = mf%imapr(3,n,1,iface)
                        I  = G%intma(ir,jr,kr,er)
                        hi = b%psiq(n,iquad)
                        qr(:,k) = qr(:,k) + hi*qprime_df(:,I,k)
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

                uu = 0.5*(ul+ur)
                vv = 0.5*(vl+vr)
                udpl(k) = ul*dpl
                udpr(k) = ur*dpr
                vdpl(k) = vl*dpl
                vdpr(k) = vr*dpr

                un = uu*nxl + vv*nyl
                if(un > 0.0) then
                    udp_flux(1,k,iquad) = uu * (ul*dpl)
                    vdp_flux(1,k,iquad) = uu * (vl*dpl)
                    udp_flux(2,k,iquad) = vv * (ul*dpl)
                    vdp_flux(2,k,iquad) = vv * (vl*dpl)
                else
                    udp_flux(1,k,iquad) = uu * (ur*dpr)
                    vdp_flux(1,k,iquad) = uu * (vr*dpr)
                    udp_flux(2,k,iquad) = vv * (ur*dpr)
                    vdp_flux(2,k,iquad) = vv * (vr*dpr)
                end if

            end do

            uu_dp_flux_deficit(1) = btp%Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,:,iquad))
            uu_dp_flux_deficit(2) = btp%Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,:,iquad))
            vv_dp_flux_deficit(1) = btp%Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,:,iquad))
            vv_dp_flux_deficit(2) = btp%Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,:,iquad))

            do k = 1, inp%nlayers

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
            do k = 1, inp%nlayers
                p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,k)
                p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,k)
            end do

            one_plus_eta_edge = sqrt(btp%one_plus_eta_edge_2_ave(iquad,iface))
            z_face(1,inp%nlayers+1)    = init%zbot_face(1,iquad,iface)
            z_face(2,inp%nlayers+1)    = init%zbot_face(2,iquad,iface)
            z_edge_plus(inp%nlayers+1) = init%zbot_face(1,iquad,iface)
            z_edge_minus(inp%nlayers+1)= init%zbot_face(2,iquad,iface)
            do k = inp%nlayers, 1, -1
                z_face(1,k)    = z_face(1,k+1)    + alpha_over_g(k) * (ope_l * ql(1,k))
                z_face(2,k)    = z_face(2,k+1)    + alpha_over_g(k) * (ope_r * qr(1,k))
                z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k) * (one_plus_eta_edge * ql(1,k))
                z_edge_minus(k)= z_edge_minus(k+1)+ alpha_over_g(k) * (one_plus_eta_edge * qr(1,k))
            end do

            p_edge_plus(2)  = one_plus_eta_edge * ql(1,1)
            p_edge_minus(2) = one_plus_eta_edge * qr(1,1)
            do k = 2, inp%nlayers
                p_edge_plus(k+1)  = p_edge_plus(k)  + one_plus_eta_edge * ql(1,k)
                p_edge_minus(k+1) = p_edge_minus(k) + one_plus_eta_edge * qr(1,k)
            end do

            do k = 1, inp%nlayers

                H_r_plus = 0.5*init%alpha_mlswe(k)*(p_edge_plus(k+1)**2 - p_edge_plus(k)**2)

                H_r_minus = 0.0
                do ktemp = 1, inp%nlayers
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
                do ktemp = 1, inp%nlayers
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
                do k = 1, inp%nlayers
                    p2l(k+1) = p_face(1,k+1)
                    H_face(1,k,iquad) = 0.5*init%alpha_mlswe(k)*(p2l(k+1)**2 - p2l(k)**2)
                    p2r(k+1) = p_face(2,k+1)
                    H_face(2,k,iquad) = 0.5*init%alpha_mlswe(k)*(p2r(k+1)**2 - p2r(k)**2)
                end do
            end if

            if(er /= -4) then
                do k = 1, inp%nlayers-1
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

            do k = 1, inp%nlayers

                hl = H_face(1,k,iquad)
                hr = H_face(2,k,iquad)

                flux_xl = nxl*(udp_flux(1,k,iquad) + hl) + nyl*udp_flux(2,k,iquad)
                flux_xr = nxl*(udp_flux(1,k,iquad) + hr) + nyl*udp_flux(2,k,iquad)
                flux_yl = nxl*vdp_flux(1,k,iquad) + nyl*(vdp_flux(2,k,iquad) + hl)
                flux_yr = nxl*vdp_flux(1,k,iquad) + nyl*(vdp_flux(2,k,iquad) + hr)

                do n = 1, b%ngl
                    hi = b%psiq(n,iquad)
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    rhs_mom(1,I,k) = rhs_mom(1,I,k) - wq*hi*flux_xl
                    rhs_mom(2,I,k) = rhs_mom(2,I,k) - wq*hi*flux_yl

                    if(er > 0) then
                        ir = mf%imapr(1,n,1,iface)
                        jr = mf%imapr(2,n,1,iface)
                        kr = mf%imapr(3,n,1,iface)
                        I  = G%intma(ir,jr,kr,er)
                        rhs_mom(1,I,k) = rhs_mom(1,I,k) + wq*hi*flux_xr
                        rhs_mom(2,I,k) = rhs_mom(2,I,k) + wq*hi*flux_yr
                    end if
                end do
            end do
        end do

    end subroutine Apply_layers_fluxes

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

    subroutine create_layers_volume_mass(G, inp, b, btp, bcl, init, tsp, dp_advec, qprime_df)

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(btp_CS),    intent(in)    :: btp
        type(bcl_CS),    intent(inout) :: bcl
        type(initial),   intent(in)    :: init
        type(tensor_CS), intent(in)    :: tsp

        real, dimension(3, G%npoin, inp%nlayers), intent(in)  :: qprime_df
        real, dimension(G%npoin, inp%nlayers),    intent(out) :: dp_advec

        real :: wq, hi, dp_temp
        integer :: k, I, Iq, ip
        real, dimension(3) :: qp, qb
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
        real, dimension(3, G%npoin, inp%nlayers), intent(in)    :: qprime_df

        integer :: k, iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n
        real :: wq, nxl, nyl, hi
        real, dimension(b%nq, inp%nlayers) :: flux_edge_u, flux_edge_v, flux_dp
        real :: dpl, dpr, uu, vv, flux, un, ul, ur, vl, vr, weight, dp_deficit(2)
        real :: dp_lr(2,inp%nlayers)
        real, dimension(3) :: ql, qr, qbl, qbr
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

end module mod_create_rhs_mlswe
