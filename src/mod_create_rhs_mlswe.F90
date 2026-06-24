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

    public :: rhs_layer_shear_stress, bcl_rhs

contains

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

        call bcl_create_precommunicator(G, inp, b, mf, par, ref, mpic, qprime_df)
        call create_rhs_dynamics_volume_bcl(G, inp, b, btp, bcl, init, tsp, mt, rhs, qprime_df, q_df)
        call Apply_bcl_fluxes(G, inp, b, mf, btp, bcl, init, rhs, qprime_df)
        call bcl_create_postcommunicator(G, inp, b, mf, par, btp, init, ref, mpic, rhs)

    end subroutine bcl_rhs

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

    subroutine create_rhs_dynamics_volume_bcl(G, inp, b, btp, bcl, init, tsp, mt, rhs, qprime_df, q_df)

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(btp_CS),    intent(in)    :: btp
        type(bcl_CS),    intent(inout) :: bcl
        type(initial),   intent(in)    :: init
        type(tensor_CS), intent(in)    :: tsp
        type(metrics),   intent(in)    :: mt

        real, dimension(3, G%npoin, inp%nlayers), intent(out) :: rhs
        real, dimension(3, G%npoin, inp%nlayers), intent(in)  :: q_df, qprime_df

        real :: wq, hi, dhdx, dhdy, bot_layer, tau_wind_u, tau_wind_v, temp1
        real :: Hq, source_x, source_y, Pstress, Pbstress
        integer :: k, I, Iq, ip, ie, ijq
        integer :: iq0, jq0, n0, n_L, n_R, m_L, m_R
        integer :: ip_xl, ip_xr, ip_yl, ip_yr
        real, dimension(inp%nlayers+1) :: pprime_temp
        real :: tempbot, weight, acceleration, pbq, ramp
        real :: Dxi_z, Deta_z, e_x, e_y, nx_m, ny_m
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
        real, dimension(b%npts, inp%nlayers+1) :: z_elem
        logical :: wet_elem(b%npts, inp%nlayers)
        logical :: semi_dry(inp%nlayers)
        real :: rhs_loc(3, b%npts, inp%nlayers)

        rhs = 0.0
        bot_layer = 0.0
        bcl%sum_layer_mass_flux = 0.0

        Pstress  = (gravity/init%alpha_mlswe(1))            * 50.0
        Pbstress = (gravity/init%alpha_mlswe(inp%nlayers))  * 10.0

        do ie = 1, G%nelem

            rhs_loc = 0.0

            ! Pre-compute interface heights and wet flags at all DOF nodes of element ie.
            ! For DG each DOF belongs to one element; any Iq within the element
            ! gives the same ip->I mapping via tsp%indexq.
            Iq = tsp%indexq_e(1, ie)
            do ip = 1, b%npts
                I = tsp%indexq(ip, Iq)
                z_elem(ip, inp%nlayers+1) = init%zbot_df(I)
                do k = inp%nlayers, 1, -1
                    z_elem(ip, k) = z_elem(ip, k+1) + (init%alpha_mlswe(k)/gravity) * &
                                    (sqrt(btp%ope2_ave_df(I)) * qprime_df(1, I, k))
                    wet_elem(ip, k) = (bcl%dry_flg(ie, k) == 0)
                end do
            end do

            do ijq = 1, b%npts_quad

                Iq = tsp%indexq_e(ijq, ie)

                p_tmp(1)       = 0.0
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

                    if (qp(1) * qb(1) < (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff) then
                        qp(1) = 0.0
                        p_tmp(k+1) = p_tmp(k)
                        H_tmp(k)   = 0.0
                        dp(k) = 0.0
                        udp(k) = 0.0; vdp(k) = 0.0
                        u_udp(k) = 0.0; v_vdp(k) = 0.0
                        u_vdp(1,k) = 0.0; u_vdp(2,k) = 0.0
                        temp_uu(k) = eps1; temp_vv(k) = eps1
                    else
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
                    end if

                    pprime_temp(k+1) = pprime_temp(k) + qp(1)
                end do

                ! DG gradient of interface heights using pre-computed z_elem.
                gradz(:,:) = 0.0; pbq = 0.0
                do ip = 1, b%npts
                    I = tsp%indexq(ip, Iq)
                    do k = 1, inp%nlayers
                        gradz(1,k) = gradz(1,k) + tsp%dpsidx(ip,Iq) * z_elem(ip,k)
                        gradz(2,k) = gradz(2,k) + tsp%dpsidy(ip,Iq) * z_elem(ip,k)
                    end do
                    pbq = pbq + tsp%psih(ip,Iq) * init%pbprime_df(I)
                end do
                gradz(1,inp%nlayers+1) = init%grad_zbot_quad(1,Iq)
                gradz(2,inp%nlayers+1) = init%grad_zbot_quad(2,Iq)

                ! FD override for semi-dry layers: replace the DG gradz with a
                ! nearest-neighbor FD that uses only wet DOF nodes, suppressing
                ! the spurious pressure gradient across the dry/wet interface.
                ! if (bcl%dry_flg(ie, k) == 1) then
                !     iq0 = mod(ijq-1, b%nqx) + 1
                !     jq0 = (ijq-1) / b%nqx + 1

                !     ! Bracket quad point in ξ: n_L = largest n with xglx(n) <= xnqx(iq0)
                !     n_L = 1
                !     do n0 = 1, b%nglx - 1
                !         if (b%xglx(n0) <= b%xnqx(iq0)) n_L = n0
                !     end do
                !     n_R = min(n_L + 1, b%nglx)

                !     ! Bracket quad point in η: m_L = largest m with xgly(m) <= xnqy(jq0)
                !     m_L = 1
                !     do n0 = 1, b%ngly - 1
                !         if (b%xgly(n0) <= b%xnqy(jq0)) m_L = n0
                !     end do
                !     m_R = min(m_L + 1, b%ngly)

                !     e_x  = mt%ksiq_x(iq0, jq0, 1, ie)
                !     e_y  = mt%ksiq_y(iq0, jq0, 1, ie)
                !     nx_m = mt%etaq_x(iq0, jq0, 1, ie)
                !     ny_m = mt%etaq_y(iq0, jq0, 1, ie)

                !     do k = 1, inp%nlayers
                !         if (.not. semi_dry(k)) cycle

                !         ! ξ FD: left/right DOF nodes in η row m_L
                !         ip_xl = (m_L-1)*b%nglx + n_L
                !         ip_xr = (m_L-1)*b%nglx + n_R
                !         if (wet_elem(ip_xl,k) .and. wet_elem(ip_xr,k)) then
                !             Dxi_z = (z_elem(ip_xr,k) - z_elem(ip_xl,k)) / &
                !                     (b%xglx(n_R) - b%xglx(n_L))
                !         else
                !             Dxi_z = 0.0
                !         end if

                !         ! η FD: bottom/top DOF nodes in ξ column n_L
                !         ip_yl = (m_L-1)*b%nglx + n_L
                !         ip_yr =  m_L   *b%nglx + n_L
                !         if (wet_elem(ip_yl,k) .and. wet_elem(ip_yr,k)) then
                !             Deta_z = (z_elem(ip_yr,k) - z_elem(ip_yl,k)) / &
                !                      (b%xgly(m_R) - b%xgly(m_L))
                !         else
                !             Deta_z = 0.0
                !         end if

                !         gradz(1,k) = Dxi_z * e_x + Deta_z * nx_m
                !         gradz(2,k) = Dxi_z * e_y + Deta_z * ny_m
                !     end do
                ! end if

                uu_dp_deficitq = btp%Qu_ave(Iq)  - sum(u_udp(:))
                uv_dp_deficitq = btp%Quv_ave(Iq) - sum(u_vdp(1,:))
                vv_dp_deficitq = btp%Qv_ave(Iq)  - sum(v_vdp(:))

                one_over_sumuq = 1.0/sum(temp_uu(:))
                one_over_sumvq = 1.0/sum(temp_vv(:))

                ! ramp = 1.0
                ! do k = 1, inp%nlayers
                !     if (dp(k) > (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff) then
                !         ramp = min(ramp, dp(k) / (100.0 * (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff))
                !     else
                !         ramp = 0.0
                !     end if
                ! end do

                wq = tsp%wjac(Iq)

                do k = 1, inp%nlayers

                    ! Skip entirely if every layer of this element is fully dry.
                    ! if (bcl%dry_flg(ie,k) == 2) cycle

                    if (dp(k) > (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff) then
                        weightq    = temp_uu(k) * one_over_sumuq
                        u_udp(k)   = u_udp(k)   + weightq * uu_dp_deficitq
                        u_vdp(1,k) = u_vdp(1,k) + weightq * uv_dp_deficitq

                        weightq    = temp_vv(k) * one_over_sumvq
                        u_vdp(2,k) = u_vdp(2,k) + weightq * uv_dp_deficitq
                        v_vdp(k)   = v_vdp(k)   + weightq * vv_dp_deficitq
                    end if

                    ! ramp = 1.0
                    ! if (dp(k) > (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff) then
                    !     ramp = min(ramp, dp(k) / (100.0 * (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff))
                    ! else
                    !     ramp = 0.0
                    ! end if
                    Hq = H_tmp(k)
                    if (Hq > 0.0) then
                        weight = 1.0
                        acceleration = sum(H_tmp(:))
                        ! if(acceleration > 0.0) weight = 1.0 + ramp * (btp%H_ave(Iq) / acceleration - 1.0)
                        if(acceleration > 0.0) weight = btp%H_ave(Iq) / acceleration
                        Hq = Hq * weight
                    end if

                    if (dp(k) > (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff) then
                        flux(1,1) = udp(k) + (dp(k)/sum(dp(:))) * (btp%btp_mass_flux_ave(1,Iq) - sum(udp(:)))
                        flux(2,1) = vdp(k) + (dp(k)/sum(dp(:))) * (btp%btp_mass_flux_ave(2,Iq) - sum(vdp(:)))
                    else
                        flux(1,1) = udp(k)
                        flux(2,1) = vdp(k)
                    end if

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

                    if (trim(inp%bcl_time_method) == 'rk3' .or. trim(inp%bcl_time_method) == 'lsrk3') then
                        source_x = source_x + init%coriolis_quad(Iq) * vdp(k)
                        source_y = source_y - init%coriolis_quad(Iq) * udp(k)
                    endif

                    do ip = 1, b%npts
                        hi   = tsp%psih(ip,Iq)
                        dhdx = tsp%dpsidx(ip,Iq)
                        dhdy = tsp%dpsidy(ip,Iq)
                        rhs_loc(1,ip,k) = rhs_loc(1,ip,k) + wq*(dhdx*flux(1,1) + dhdy*flux(2,1))
                        rhs_loc(2,ip,k) = rhs_loc(2,ip,k) + wq*(hi*source_x + dhdx*flux(1,2) + dhdy*flux(2,2))
                        rhs_loc(3,ip,k) = rhs_loc(3,ip,k) + wq*(hi*source_y + dhdx*flux(1,3) + dhdy*flux(2,3))
                    end do
                end do
            end do

            ! Scatter element-local RHS to global array.
            Iq = tsp%indexq_e(1, ie)
            do ip = 1, b%npts
                I = tsp%indexq(ip, Iq)
                do k = 1, inp%nlayers
                    rhs(1,I,k) = rhs_loc(1,ip,k)
                    rhs(2,I,k) = rhs_loc(2,ip,k)
                    rhs(3,I,k) = rhs_loc(3,ip,k)
                end do
            end do
        end do

    end subroutine create_rhs_dynamics_volume_bcl

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
        real :: H_corr1, p_inc1, H_corr2, p_inc2, weight, ope_l, ope_r, ramp
        integer :: el, er
        real :: ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv, flux, un
        real, dimension(inp%nlayers) :: udpl, udpr, vdpl, vdpr
        real :: uu_dp_flux_deficit(2), vv_dp_flux_deficit(2), dp_deficit(2)
        real, parameter :: eps1 = 1.0e-20
        integer :: il, jl, ir, jr, kl, kr, n
        real :: wq, hi, flux_xl, flux_xr, flux_yl, flux_yr, hl, hr
        real, dimension(2,inp%nlayers,b%nq) :: H_face, udp_flux, vdp_flux, dp_flux
        real :: dp_lr(2,inp%nlayers)

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

                ! Zero pressure-gradient contribution from layers below dry threshold.
                ! Using zero (not a floor) prevents overflow when BTP ope is large.
                if (ql(1,k) * qbl(1) < g_over_alpha(k) * inp%dry_cutoff) then
                    ql(1,k) = 0.0
                    ql(2,k) = 0.0; ql(3,k) = 0.0
                end if
                if (qr(1,k) * qbr(1) < g_over_alpha(k) * inp%dry_cutoff) then
                    qr(1,k) = 0.0
                    qr(2,k) = 0.0; qr(3,k) = 0.0
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

                if (dpl > g_over_alpha(k) * inp%dry_cutoff) then
                    udpl(k) = ul*dpl
                    vdpl(k) = vl*dpl
                else
                    udpl(k) = 0.0
                    vdpl(k) = 0.0
                end if

                if (dpr > g_over_alpha(k) * inp%dry_cutoff) then
                    udpr(k) = ur*dpr
                    vdpr(k) = vr*dpr
                else
                    udpr(k) = 0.0
                    vdpr(k) = 0.0
                end if

                un = uu*nxl + vv*nyl
                if(un > 0.0) then
                    dp_flux(1,k,iquad)  = uu * dpl
                    udp_flux(1,k,iquad) = uu * udpl(k)
                    vdp_flux(1,k,iquad) = uu * vdpl(k)
                    dp_flux(2,k,iquad)  = vv * dpl
                    udp_flux(2,k,iquad) = vv * udpl(k)
                    vdp_flux(2,k,iquad) = vv * vdpl(k)
                else
                    dp_flux(1,k,iquad)  = uu * dpr
                    udp_flux(1,k,iquad) = uu * udpr(k)
                    vdp_flux(1,k,iquad) = uu * vdpr(k)
                    dp_flux(2,k,iquad)  = vv * dpr
                    udp_flux(2,k,iquad) = vv * udpr(k)
                    vdp_flux(2,k,iquad) = vv * vdpr(k)
                end if

            end do

            uu_dp_flux_deficit(1) = btp%Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,:,iquad))
            uu_dp_flux_deficit(2) = btp%Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,:,iquad))
            vv_dp_flux_deficit(1) = btp%Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,:,iquad))
            vv_dp_flux_deficit(2) = btp%Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,:,iquad))

            dp_deficit(1) = btp%btp_mass_flux_face_ave(1,iquad,iface) - sum(dp_flux(1,:,iquad))
            dp_deficit(2) = btp%btp_mass_flux_face_ave(2,iquad,iface) - sum(dp_flux(2,:,iquad))

            do k = 1, inp%nlayers

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

            ! Zero H_face for vanishing layers (correction loop can leave
            ! non-physical values when dry layer k borders a wet layer).
            ! do k = 1, inp%nlayers
            !     if (dp_lr(1,k) < g_over_alpha(k) * inp%dry_cutoff) H_face(1,k,iquad) = 0.0
            !     if (dp_lr(2,k) < g_over_alpha(k) * inp%dry_cutoff) H_face(2,k,iquad) = 0.0
            ! end do

            ! if (dp_lr(1,k) < g_over_alpha(k) * inp%dry_cutoff) H_face(1,k,iquad) = 0.0
            ! if (dp_lr(2,k) < g_over_alpha(k) * inp%dry_cutoff) H_face(2,k,iquad) = 0.0

            ! Weight correction: skip entirely when any layer is at the dry threshold.
            ! For fully-wet faces, apply the full BTP-BCL scaling.
            ! ramp = 1.0
            ! do k = 1, inp%nlayers
            !     if (dp_lr(1,k) > g_over_alpha(k) * inp%dry_cutoff) then
            !         ramp = min(ramp, dp_lr(1,k) / (100.0 * (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff))
            !     else
            !         ramp = 0.0
            !     end if
            ! end do

            ! if (dp_lr(1,k) > g_over_alpha(k) * inp%dry_cutoff) then
            !     ramp = min(ramp, dp_lr(1,k) / (100.0 * (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff))
            ! else
            !     ramp = 0.0
            ! end if
            acceleration = sum(H_face(1,:,iquad))
            if (abs(acceleration) > eps1) then
                ! weight = 1.0 + ramp * (btp%H_face_ave(iquad,iface) / acceleration - 1.0)
                weight = btp%H_face_ave(iquad,iface) / acceleration
                do k = 1, inp%nlayers
                    if (H_face(1,k,iquad) /= 0.0) H_face(1,k,iquad) = H_face(1,k,iquad) * weight
                end do
            end if

            ramp = 1.0
            do k = 1, inp%nlayers
                if (dp_lr(2,k) > g_over_alpha(k) * inp%dry_cutoff) then
                    ramp = min(ramp, dp_lr(2,k) / (100.0 * (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff))
                else
                    ramp = 0.0
                end if
            end do

            ! if (dp_lr(2,k) > g_over_alpha(k) * inp%dry_cutoff) then
            !     ramp = min(ramp, dp_lr(2,k) / (100.0 * (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff))
            ! else
            !     ramp = 0.0
            ! end if
            acceleration = sum(H_face(2,:,iquad))
            if (abs(acceleration) > eps1) then
                ! weight = 1.0 + ramp * (btp%H_face_ave(iquad,iface) / acceleration - 1.0)
                weight = btp%H_face_ave(iquad,iface) / acceleration
                do k = 1, inp%nlayers
                    if (H_face(2,k,iquad) /= 0.0) H_face(2,k,iquad) = H_face(2,k,iquad) * weight
                end do
            end if

            wq = mf%jac_faceq(iquad,1,iface)

            do k = 1, inp%nlayers

                ! Skip faces where both neighboring elements are fully dry in all layers.
                ! if (bcl%dry_flg(el,k) == 2 .and. &
                !     (er <= 0 .or. bcl%dry_flg(er,k) == 2)) cycle

                ! if (bcl%dry_flg(el,k) == 2 .and. (er <= 0 .or. bcl%dry_flg(er,k) == 2)) cycle

                flux = nxl*dp_flux(1,k,iquad) + nyl*dp_flux(2,k,iquad)

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
                    rhs(1,I,k) = rhs(1,I,k) - wq*hi*flux
                    rhs(2,I,k) = rhs(2,I,k) - wq*hi*flux_xl
                    rhs(3,I,k) = rhs(3,I,k) - wq*hi*flux_yl
                end do

                if(er > 0) then
                    do n = 1, b%ngl
                        hi = b%psiq(n,iquad)
                        ir = mf%imapr(1,n,1,iface)
                        jr = mf%imapr(2,n,1,iface)
                        kr = mf%imapr(3,n,1,iface)
                        I  = G%intma(ir,jr,kr,er)
                        rhs(1,I,k) = rhs(1,I,k) + wq*hi*flux
                        rhs(2,I,k) = rhs(2,I,k) + wq*hi*flux_xr
                        rhs(3,I,k) = rhs(3,I,k) + wq*hi*flux_yr
                    end do
                end if
            end do
        end do

    end subroutine Apply_bcl_fluxes

end module mod_create_rhs_mlswe
