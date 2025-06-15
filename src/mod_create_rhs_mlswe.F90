! =======================================================================
!>@brief This subroutine builds the RHS vector for the DG method for MLSWE
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: April 03, 2023
! =======================================================================
module mod_create_rhs_mlswe

    use mod_constants, only: gravity
    use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad
    use mod_basis, only: nglx, ngly, nglz, npts, dpsiqx, dpsiqy, dpsiqz, nqx, nqy, nqz, &
                            psiqx, psiqy, psiqz
    use mod_grid, only: intma, coord
    use mod_initial, only: q_ref, coriolis_constant, kvector, nvar
    use mod_input, only: nlayers
    use mod_metrics, only: ksiq_x, ksiq_y, ksiq_z, &
                           etaq_x, etaq_y, etaq_z, &
                           zetaq_x, zetaq_y, zetaq_z, &
                           jacq, massinv

    public :: layer_momentum_rhs, btp_mass_advection_rhs, &
              interpolate_layer_from_quad_to_node, rhs_layer_shear_stress, layer_mass_rhs, &
              consistency_mass_rhs
              
    contains

    subroutine layer_momentum_rhs(rhs_mom, rhs_visc, qprime_df, q_df, qprime_df_face)

        use mod_input, only: nlayers
        use mod_metrics, only: massinv
        use mod_grid, only: npoin, npoin_q, intma, intma_dg_quad, nface
        use mod_basis, only: nq, ngl

        implicit none
        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom
        real, dimension(2,npoin,nlayers), intent(in) :: rhs_visc
        real, dimension(3,npoin,nlayers), intent(in) :: q_df, qprime_df
        real, dimension(3,2,ngl,nface,nlayers), intent(in) :: qprime_df_face

        integer :: k,I

        call create_rhs_dynamics_volume_layers(rhs_mom, qprime_df, q_df)
        call Apply_layers_fluxes(rhs_mom, qprime_df_face)

        do k = 1, nlayers
            rhs_mom(1,:,k) = massinv(:)*rhs_mom(1,:,k) + rhs_visc(1,:,k)
            rhs_mom(2,:,k) = massinv(:)*rhs_mom(2,:,k) + rhs_visc(2,:,k)
        end do

    end subroutine layer_momentum_rhs

    subroutine layer_mass_rhs(dp_advec, qprime_df, qprime_df_face)

        use mod_metrics, only: massinv
        use mod_input, only: nlayers
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq, ngl

        implicit none

        real, dimension(npoin, nlayers), intent(out) :: dp_advec
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(3,2,ngl,nface,nlayers), intent(in) :: qprime_df_face
        
        integer :: k

        ! Compute the mass advection term for the degree of freedom for dp in each layer
        call create_layers_volume_mass(dp_advec, qprime_df)

        ! Compute the mass flux term 
        call create_layer_mass_flux(dp_advec, qprime_df_face)

        do k = 1, nlayers
            dp_advec(:,k) = massinv(:)*dp_advec(:,k)
        end do
        
    end subroutine layer_mass_rhs

    subroutine consistency_mass_rhs(dp_advec, qprime, flux_deficit_mass_face)

        use mod_metrics, only: massinv
        use mod_input, only: nlayers
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq

        implicit none

        real, dimension(npoin, nlayers), intent(out) :: dp_advec
        real, dimension(3, npoin_q, nlayers), intent(in) :: qprime
        real, dimension(2,2,nq,nface,nlayers), intent(in)  :: flux_deficit_mass_face
        
        integer :: k

        ! Compute the mass advection term for the degree of freedom for dp in each layer
        call create_consistency_volume_mass(dp_advec, qprime)

        ! Compute the mass flux term 
        call create_consistency_mass_flux(dp_advec, flux_deficit_mass_face)

    end subroutine consistency_mass_rhs

    subroutine interpolate_layer_from_quad_to_node(q_df,q)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_metrics, only: massinv
        use mod_input, only: nlayers
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac

        implicit none
        
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: q

        real :: wq, var_u(nlayers), var_v(nlayers), hi
        integer :: iquad, jquad, kquad, e, k, m, n, l, I, Iq, ip

        q_df(2:3,:,:) = 0.0
        
        ! do k = 1, nlayers

        do Iq = 1,npoin_q
            
            wq = wjac(Iq)
            var_u = q(2,Iq,:)
            var_v = q(3,Iq,:)

            do ip = 1, npts

                I = indexq(ip,Iq)
                hi = psih(ip,Iq)

                q_df(2,I,:) = q_df(2,I,:) + wq*var_u(:)*hi
                q_df(3,I,:) = q_df(3,I,:) + wq*var_v(:)*hi
            end do
        end do

        do k = 1, nlayers
            q_df(2,:,k) = massinv(:)*q_df(2,:,k)
            q_df(3,:,k) = massinv(:)*q_df(3,:,k)  
        end do
        
    end subroutine interpolate_layer_from_quad_to_node

    subroutine rhs_layer_shear_stress(rhs_stress,q)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_input, only: nlayers, ad_mlswe, max_shear_dz
        use mod_initial, only: coriolis_quad, alpha_mlswe
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac

        implicit none
        
        real, intent(in) :: q(2,npoin_q,nlayers)
        real, intent(out) :: rhs_stress(2,npoin,nlayers)
        
        real :: tau_u(npoin_q,nlayers+1), tau_v(npoin_q,nlayers+1), coeff(npoin_q)
        real :: tau_u_q(npoin_q, nlayers), tau_v_q(npoin_q, nlayers)
        real :: wq, hi
        integer :: k, e, iquad, jquad, kquad, l, m, n, Iq, I, ip
        
        tau_u = 0.0
        tau_v = 0.0
        rhs_stress = 0.0
        
        ! coeff = sqrt(0.5d0*coriolis_quad(:)*ad_mlswe)/alpha_mlswe(1)

        !do Iq = 1,npoin_q
        coeff(:) = max(sqrt(0.5*coriolis_quad(:)*ad_mlswe)/alpha_mlswe(1), ad_mlswe/(alpha_mlswe(1) * max_shear_dz))
        !end do
        
        do k = 2, nlayers   ! loop over interfaces between layers
            tau_u(:,k) = coeff*(q(1,:,k-1) - q(1,:,k))
            tau_v(:,k) = coeff*(q(2,:,k-1) - q(2,:,k))
        end do

        do k = 1, nlayers
            tau_u_q(:,k) = gravity*(tau_u(:,k) - tau_u(:,k+1))
            tau_v_q(:,k) = gravity*(tau_v(:,k) - tau_v(:,k+1))
        end do
        
        ! do k=1,nlayers
            
        !     tau_u_q = gravity*(tau_u(:,k) - tau_u(:,k+1))
        !     tau_v_q = gravity*(tau_v(:,k) - tau_v(:,k+1))

        do Iq = 1,npoin_q

            wq = wjac(Iq)
            do ip = 1, npts

                I = indexq(ip,Iq)
                hi = psih(ip,Iq)

                rhs_stress(1,I,:) = rhs_stress(1,I,:) + wq*hi*tau_u_q(Iq,:)
                rhs_stress(2,I,:) = rhs_stress(2,I,:) + wq*hi*tau_v_q(Iq,:)
            end do
        end do

        do k = 1, nlayers
            rhs_stress(1,:,k) = massinv(:)*rhs_stress(1,:,k)
            rhs_stress(2,:,k) = massinv(:)*rhs_stress(2,:,k)
        end do
        
    end subroutine rhs_layer_shear_stress

    subroutine create_rhs_dynamics_volume_layers(rhs_mom, qprime_df,q_df)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts, dpsiqx
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, alpha_mlswe, tau_wind, pbprime, zbot_df
        use mod_variables, only: tau_wind_int, tau_bot_int, tau_bot_ave, ope_ave, uvb_ave, H_ave, &
                                    Qu_ave, Qv_ave, Quv_ave, uvb_ave, ope_ave_df

        implicit none

        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom
        real, dimension(3,npoin,nlayers), intent(in) :: q_df, qprime_df

        real :: wq, hi, dhdx, dhdy, bot_layer, tau_wind_u, tau_wind_v, temp1
        real :: Hq, var_uu, var_uv, var_vu, var_vv, source_x, source_y, Pstress
        integer :: k, I, Iq, ip
        real, dimension(nlayers+1) :: pprime_temp
        real :: temp_dp, temp_u, temp_v, Ptop_k, Pbot_k, tempbot, Pbstress
        real :: weight, acceleration
        real, dimension(3) :: qp, qb
        real, dimension(nlayers) :: temp_uu, temp_vv, H_tmp, u_udp, v_vdp
        real, dimension(2,nlayers) :: u_vdp
        real :: p_tmp(nlayers+1), u, v, dp, weightq, one_over_sumuq, one_over_sumvq
        real :: uu_dp_deficitq, uv_dp_deficitq, vv_dp_deficitq, gradz(2,nlayers+1)
        real :: z_elv(npoin,nlayers+1)
        real, parameter :: eps1 = 1.0e-20 !  Parameter used to prevent division by zero.

        rhs_mom = 0.0
        bot_layer = 0.0
        pprime_temp = 0.0

        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth at which bottom stress is reduced to 0

        ! Find layer interfaces
        z_elv(:,nlayers+1) = zbot_df(:)
        do k = nlayers,1,-1
            z_elv(:,k) = z_elv(:,k+1) + (alpha_mlswe(k)/gravity) * (qprime_df(1,:,k)*ope_ave_df(:))
        end do

        do Iq = 1, npoin_q

            p_tmp(1) = 0.0
            temp_uu = 0.0; temp_vv = 0.0
            do k = 1,nlayers
                qp = 0.0; qb = 0.0
                do ip = 1,npts
                    I = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    qp(:) = qp(:) + hi*qprime_df(:,I,k)
                    !qb(:) = qb(:) + hi*uvb_ope_ave_df(:,I)

                    temp_uu(k) = temp_uu(k) + hi*q_df(2,I,k)
                    temp_vv(k) = temp_vv(k) + hi*q_df(3,I,k)
                enddo

                qb(1) = ope_ave(Iq)
                qb(2) = uvb_ave(1,Iq)
                qb(3) = uvb_ave(2,Iq)

                p_tmp(k+1) = p_tmp(k) + qp(1) * qb(1)
                H_tmp(k) = 0.5*alpha_mlswe(k) * (p_tmp(k+1)**2 - p_tmp(k)**2)

                dp = qp(1) * qb(1)
                u = qp(2) + qb(2)
                v = qp(3) + qb(3)

                u_udp(k) = dp * u**2
                v_vdp(k) = dp * v**2
                u_vdp(1,k) = u * v * dp
                u_vdp(2,k) = v * u * dp

                temp_uu(k) = abs(temp_uu(k)) + eps1
                temp_vv(k) = abs(temp_vv(k)) + eps1
            end do

            gradz = 0.0
            do ip = 1,npts
                I = indexq(ip,Iq)
                gradz(1,:) = gradz(1,:) + dpsidx(ip,Iq)*z_elv(I,:)
                gradz(2,:) = gradz(2,:) + dpsidy(ip,Iq)*z_elv(I,:)
            enddo

            ! Consistency
            uu_dp_deficitq = Qu_ave(Iq) - sum(u_udp(:))
            uv_dp_deficitq = Quv_ave(Iq) - sum(u_vdp(1,:))
            vv_dp_deficitq = Qv_ave(Iq) - sum(v_vdp(:))

            one_over_sumuq = 1.0/sum(temp_uu(:))
            one_over_sumvq = 1.0/sum(temp_vv(:))

            wq = wjac(Iq)
            pprime_temp = 0.0

            do k = 1,nlayers

                pprime_temp(k+1) = pprime_temp(k) + qp(k)

                weightq = temp_uu(k) * one_over_sumuq
                u_udp(k) = u_udp(k) + weightq * uu_dp_deficitq
                u_vdp(1,k) = u_vdp(1,k) + weightq * uv_dp_deficitq

                weightq = temp_vv(k) * one_over_sumvq
                u_vdp(2,k) = u_vdp(2,k) + weightq * uv_dp_deficitq
                v_vdp(k) = v_vdp(k) + weightq * vv_dp_deficitq

                Hq = H_tmp(k)
                
                ! Adjust the values of  H_r, at quadrature points, so that the vertical sum of
                ! H_r  over all layers equals the time average of the barotropic forcing  H  over all barotropic substeps of the baroclinic
                ! time interval.
                ! The difference between the time-averaged  H  and the vertical sum of
                ! H_r  must be distributed over the layers via some sort of
                ! weighting scheme.
                ! Weight according to the current value of H_r.
                !   That is, the weight for layer r is
                !   H_r / (sum of H_s over all layers s).
                !   The adjusted  H_r  is then
                !   (H_r)_adjusted  =   H_r + [ H_r/(sum H_s)] * [ H_ave - sum(H_s) ]
                !                       =   H_r +   H_r * H_ave/(sum H_s)  -  H_r
                !                       =   H_r * H_ave/(sum H_s)
                !   Therefore, at each quadrature point and cell edge, multiply
                !   the current value of  H_r  by the layer-independent ratio
                !   H_ave/(sum H_s),  which should be approximately equal to  1.

                weight = 1.0
                acceleration = sum(H_tmp(:))
                if(acceleration > 0.0) then
                    weight = H_ave(Iq) / acceleration
                end if
                Hq = Hq * weight

                var_uu = u_udp(k)
                var_uv = u_vdp(1,k)
                var_vu = u_vdp(2,k)
                var_vv = v_vdp(k)

                temp1 = (min(pprime_temp(k+1), Pstress) - min(pprime_temp(k), Pstress))/ Pstress
                tau_wind_u = temp1*tau_wind(1,Iq)
                tau_wind_v = temp1*tau_wind(2,Iq)

                tempbot = min(Pbstress,pbprime(Iq)-pprime_temp(k+1)) - min(Pbstress, pbprime(Iq) - pprime_temp(k))
                tempbot = tempbot / Pbstress

                source_x = gravity*(tau_wind_u - tempbot*tau_bot_ave(1,Iq) + &
                            p_tmp(k) * gradz(1,k) - p_tmp(k+1) * gradz(1,k+1))
                source_y = gravity*(tau_wind_v - tempbot*tau_bot_ave(2,Iq) + &
                            p_tmp(k) * gradz(2,k) - p_tmp(k+1) * gradz(2,k+1))
                do ip = 1, npts

                    I = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    !Xi derivatives
                    dhdx = dpsidx(ip,Iq)
                    !Eta derivatives
                    dhdy = dpsidy(ip,Iq)

                    rhs_mom(1,I,k) = rhs_mom(1,I,k) + wq*(hi*source_x + dhdx*(Hq + var_uu) + var_uv*dhdy)
                    rhs_mom(2,I,k) = rhs_mom(2,I,k) + wq*(hi*source_y + var_vu*dhdx + dhdy*(Hq +var_vv))

                end do
            end do
        end do !k

    end subroutine create_rhs_dynamics_volume_layers

    subroutine Apply_layers_fluxes(rhs_mom, qprime_df_face)

        ! This routine computes the layer momentum advection flux terms using upwind flux
        ! This routine computes the layer momentum pressure terms

        use mod_constants, only : gravity
        use mod_initial, only : alpha_mlswe, pbprime_face, zbot_face
        use mod_grid, only : nface, npoin, npoin_q, face, intma
        use mod_basis, only : nq, psiq, ngl
        use mod_input, only : nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: ope_face_ave, H_face_ave, one_plus_eta_edge_2_ave, &
                                uvb_face_ave, Quv_face_ave, Qu_face_ave, Qv_face_ave

        implicit none

        real, dimension(2, npoin, nlayers), intent(inout) :: rhs_mom
        real, dimension(3,2,ngl,nface,nlayers), intent(in) :: qprime_df_face

        real, dimension(nlayers) :: alpha_over_g, g_over_alpha
        real, dimension(2,nlayers+1) :: p_face, z_face
        real, dimension(nlayers+1) :: p_edge_plus, p_edge_minus, p2l, p2r, z_edge_plus, z_edge_minus
        real, dimension(3,nq,nlayers) :: ql, qr
        real, dimension(3,nq) :: qbl, qbr
        integer :: iface, ilr, k, iquad, ktemp, I
        real :: z_intersect_top,z_intersect_bot, dz_intersect, H_r_plus, H_r_minus, acceleration
        real :: p_intersect_bot, p_intersect_top, one_plus_eta_edge
        real :: H_corr,p_inc, weight, H_corr1,p_inc1, H_corr2,p_inc2, temp, ope_l, ope_r
        integer :: Iq, el, er
        real ::  ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv
        real, dimension(nq,nlayers) :: udpl, udpr, vdpl, vdpr
        real :: one_over_sum_l, one_over_sum_r, uu_dp_flux_deficit, uv_dp_flux_deficit
        real :: vu_dp_flux_deficit, vv_dp_flux_deficit
        real, parameter :: eps1 = 1.0e-20 !  Parameter used to prevent division by zero.
        integer :: il, jl, ir, jr, kl, kr, jquad, n, m
        real :: wq, hi, hlx_k, hly_k, hrx_k, hry_k, flux_x, flux_y, hx_k, hy_k
        real, dimension(2,nq,nlayers) :: H_face, udp_flux, vdp_flux

        do k=1,nlayers
            alpha_over_g(k) = alpha_mlswe(k)/gravity
            g_over_alpha(k) = gravity/alpha_mlswe(k)
        enddo

        ! Compute H_r at the element face
        do iface = 1, nface

            ! Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            qbl(1,:) = ope_face_ave(1,:,iface)
            qbl(2,:) = uvb_face_ave(1,1,:,iface)
            qbl(3,:) = uvb_face_ave(2,1,:,iface)
            qbr(1,:) = ope_face_ave(2,:,iface)
            qbr(2,:) = uvb_face_ave(1,2,:,iface)
            qbr(3,:) = uvb_face_ave(2,2,:,iface)

            ql = 0.0; qr = 0.0
            do iquad = 1,nq

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                do k = 1,nlayers

                    do n = 1, ngl
                        hi = psiq(n,iquad)
                        ql(:,iquad,k) = ql(:,iquad,k) + hi*qprime_df_face(:,1,n,iface,k)
                        qr(:,iquad,k) = qr(:,iquad,k) + hi*qprime_df_face(:,2,n,iface,k)
                    enddo

                    ! Left side of the edge
                    dpl = qbl(1,iquad) * ql(1,iquad,k)
                    dpr = qbr(1,iquad) * qr(1,iquad,k)
                    ul = ql(2,iquad,k)+qbl(2,iquad)
                    ur = qr(2,iquad,k)+qbr(2,iquad)
                    vl = ql(3,iquad,k)+qbl(3,iquad)
                    vr = qr(3,iquad,k)+qbr(3,iquad)

                    uu = 0.5*(ul+ur)
                    vv = 0.5*(vl+vr)
                    udpl(iquad,k) = ul*dpl
                    udpr(iquad,k) = ur*dpr
                    vdpl(iquad,k) = vl*dpl
                    vdpr(iquad,k) = vr*dpr

                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    if(uu*nxl > 0.0) then
                        udp_flux(1,iquad,k) = uu * (ul*dpl)
                        vdp_flux(1,iquad,k) = uu * (vl*dpl)
                    else
                        udp_flux(1,iquad,k) = uu * (ur*dpr)
                        vdp_flux(1,iquad,k) = uu * (vr*dpr)
                    endif
                    if(vv*nyl > 0.0) then
                        udp_flux(2,iquad,k) = vv * (ul*dpl)
                        vdp_flux(2,iquad,k) = vv * (vl*dpl)
                    else
                        udp_flux(2,iquad,k) = vv * (ur*dpr)
                        vdp_flux(2,iquad,k) = vv * (vr*dpr)
                    endif

                enddo

                uu_dp_flux_deficit = Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,iquad,:))
                uv_dp_flux_deficit = Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,iquad,:))
                vu_dp_flux_deficit = Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,iquad,:))
                vv_dp_flux_deficit = Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,iquad,:))

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                ! Adjust the fluxes for the u-momentum equation
                one_over_sum_l = 1.0 / sum(abs(udpl(iquad,:))+eps1)
                one_over_sum_r = 1.0 / sum(abs(udpr(iquad,:))+eps1)
                !x-direction
                if(uu_dp_flux_deficit*nxl > 0.0) then
                    do k = 1,nlayers
                        weight = abs(udpl(iquad,k)) * one_over_sum_l
                        udp_flux(1,iquad,k) = udp_flux(1,iquad,k) + weight * uu_dp_flux_deficit
                    end do
                else
                    do k = 1,nlayers
                        weight = abs(udpr(iquad,k)) * one_over_sum_r
                        udp_flux(1,iquad,k) = udp_flux(1,iquad,k) + weight * uu_dp_flux_deficit
                    end do
                end if
                !y-direction
                if(uv_dp_flux_deficit*nyl > 0.0) then
                    do k = 1,nlayers
                        weight = abs(udpl(iquad,k)) * one_over_sum_l
                        udp_flux(2,iquad,k) = udp_flux(2,iquad,k) + weight * uv_dp_flux_deficit
                    end do
                else
                    do k = 1,nlayers
                        weight = abs(udpr(iquad,k)) * one_over_sum_r
                        udp_flux(2,iquad,k) = udp_flux(2,iquad,k) + weight * uv_dp_flux_deficit
                    end do
                end if
                ! Adjust the fluxes for the v-momentum equation
                one_over_sum_l = 1.0 / sum(abs(vdpl(iquad,:))+eps1)
                one_over_sum_r = 1.0 / sum(abs(vdpr(iquad,:))+eps1)
                !x-direction
                if(vu_dp_flux_deficit*nxl > 0.0) then
                    do k = 1,nlayers
                        weight = abs(vdpl(iquad,k)) * one_over_sum_l
                        vdp_flux(1,iquad,k) = vdp_flux(1,iquad,k) + weight * vu_dp_flux_deficit
                    end do
                else
                    do k = 1,nlayers
                        weight = abs(vdpr(iquad,k)) * one_over_sum_r
                        vdp_flux(1,iquad,k) = vdp_flux(1,iquad,k) + weight * vu_dp_flux_deficit
                    end do
                end if
                !y-direction
                if(vv_dp_flux_deficit*nyl > 0.0) then
                    do k = 1,nlayers
                        weight = abs(vdpl(iquad,k)) * one_over_sum_l
                        vdp_flux(2,iquad,k) = vdp_flux(2,iquad,k) + weight * vv_dp_flux_deficit
                    end do
                else
                    do k = 1,nlayers
                        weight = abs(vdpr(iquad,k)) * one_over_sum_r
                        vdp_flux(2,iquad,k) = vdp_flux(2,iquad,k) + weight * vv_dp_flux_deficit
                    end do
                end if

                z_face = 0.0 ; p_face = 0.0
                z_edge_plus = 0.0 ; z_edge_minus = 0.0
                p_edge_plus = 0.0; p_edge_minus = 0.0

                !Store Left Side Variables
                ope_l = ope_face_ave(1,iquad,iface)
                ope_r = ope_face_ave(2,iquad,iface)
                p_face(1,1) = 0.0
                p_face(2,1) = 0.0
                do k=1,nlayers
                    p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,iquad,k)
                    p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,iquad,k)
                end do

                one_plus_eta_edge = one_plus_eta_edge_2_ave(iquad,iface)
                z_face(1,nlayers+1) = zbot_face(1,iquad,iface)
                z_face(2,nlayers+1) = zbot_face(2,iquad,iface)
                z_edge_plus(nlayers+1) = zbot_face(1,iquad,iface)
                z_edge_minus(nlayers+1) = zbot_face(2,iquad,iface)
                do k=nlayers,1,-1
                    z_face(1,k) = z_face(1,k+1) + alpha_over_g(k) * (ope_l * ql(1,iquad,k))
                    z_face(2,k) = z_face(2,k+1) + alpha_over_g(k) * (ope_r * qr(1,iquad,k))
                    z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k) * (one_plus_eta_edge * ql(1,iquad,k))
                    z_edge_minus(k) = z_edge_minus(k+1) + alpha_over_g(k) * (one_plus_eta_edge * qr(1,iquad,k))
                end do

                p_edge_plus(2) = one_plus_eta_edge * ql(1,iquad,1)
                p_edge_minus(2) = one_plus_eta_edge * qr(1,iquad,1)
                do k = 2,nlayers
                    p_edge_plus(k+1) = p_edge_plus(k) + one_plus_eta_edge * ql(1,iquad,k)
                    p_edge_minus(k+1) = p_edge_minus(k) + one_plus_eta_edge * qr(1,iquad,k)
                end do

                do k = 1, nlayers

                    ! Computation from + side for layer k
                    H_r_plus = 0.5*alpha_mlswe(k)*(p_edge_plus(k+1)**2 - p_edge_plus(k)**2)

                    ! Computation from - side for layer k
                    H_r_minus = 0.0
                    do ktemp = 1, nlayers

                        z_intersect_top = min(z_edge_minus(ktemp), z_edge_plus(k))
                        z_intersect_bot = max(z_edge_minus(ktemp+1), z_edge_plus(k+1))
                        dz_intersect = z_intersect_top - z_intersect_bot

                        if (dz_intersect > 0.0) then
                            p_intersect_bot = p_edge_minus(ktemp+1) - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_minus(ktemp+1))
                            p_intersect_top = p_edge_minus(ktemp+1) - g_over_alpha(ktemp)*(z_intersect_top - z_edge_minus(ktemp+1))
                            H_r_minus = H_r_minus + 0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)

                        end if
                    end do
                    H_face(1,iquad,k) = 0.5*(H_r_plus + H_r_minus) !computation of H_r for the left side
                    ! Computation from - side for layer k
                    H_r_minus = 0.5*alpha_mlswe(k)*(p_edge_minus(k+1)**2 - p_edge_minus(k)**2)

                    ! Computation from + side for layer k
                    H_r_plus = 0.0
                    do ktemp = 1, nlayers

                        z_intersect_top = min(z_edge_plus(ktemp), z_edge_minus(k))
                        z_intersect_bot = max(z_edge_plus(ktemp+1), z_edge_minus(k+1))
                        dz_intersect = z_intersect_top - z_intersect_bot

                        if (dz_intersect > 0.0) then
                            p_intersect_bot = p_edge_plus(ktemp+1) - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_plus(ktemp+1))
                            p_intersect_top = p_edge_plus(ktemp+1) - g_over_alpha(ktemp)*(z_intersect_top - z_edge_plus(ktemp+1))
                            H_r_plus = H_r_plus + 0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
                        end if
                    end do
                    H_face(2,iquad,k) = 0.5*(H_r_plus + H_r_minus) ! computation of H_r for the right side
                end do !k

                ! Wall Boundary conditions
                if(er == -4) then
                    p2l = 0.0 ; p2r = 0.0
                    do k = 1,nlayers
                        p2l(k+1) = p_face(1,k+1)
                        H_face(1,iquad,k) = 0.5*alpha_mlswe(k)*(p2l(k+1)**2 - p2l(k)**2)
                        p2r(k+1) = p_face(2,k+1)
                        H_face(2,iquad,k) = 0.5*alpha_mlswe(k)*(p2r(k+1)**2 - p2r(k)**2)

                    end do
                end if

                if(er /= -4) then
                    do k = 1, nlayers-1          ! interface at the bottom of layer k
                        ! Corrections at the left side of a face.
                        p_inc1 = g_over_alpha(k)*(z_face(1,k+1) - z_edge_plus(k+1))
                        H_corr1 = 0.5 * alpha_mlswe(k) * ((p_face(1,k+1) + p_inc1)**2 - p_face(1,k+1)**2)
                        H_face(1,iquad,k) = H_face(1,iquad,k) - H_corr1
                        H_face(1,iquad,k+1) = H_face(1,iquad,k+1) + H_corr1

                        ! Corrections at the right side of a face.
                        p_inc2 = g_over_alpha(k)*(z_face(2,k+1) - z_edge_minus(k+1))
                        H_corr2 = 0.5 * alpha_mlswe(k) * ((p_face(2,k+1) + p_inc2)**2 - p_face(2,k+1)**2)
                        H_face(2,iquad,k) = H_face(2,iquad,k) - H_corr2
                        H_face(2,iquad,k+1) = H_face(2,iquad,k+1) + H_corr2

                    end do
                end if

                ! Adjust the values of  H_k, at element faces, so that the vertical sum of
                ! H_k  over all layers equals the time average of the barotropic forcing  H  over all barotropic substeps of the baroclinic
                ! time interval.

                ! The difference between the time-averaged  H  and the vertical sum of
                ! H_k  must be distributed over the layers via some sort of
                ! weighting scheme.
                ! Weight according to the current value of H_k.
                !   That is, the weight for layer k is
                !   H_k / (sum of H_s over all layers s).
                !   The adjusted  H_k  is then
                !   (H_k)_adjusted  =   H_k + [ H_k/(sum H_s)] * [ H_ave - sum(H_s) ]
                !                       =   H_k +   H_k * H_ave/(sum H_s)  -  H_k
                !                       =   H_k * H_ave/(sum H_s)
                !   Therefore, at each quadrature point and cell edge, multiply
                !   the current value of  H_r  by the layer-independent ratio
                !   H_ave/(sum H_s),  which should be approximately equal to  1.

                ! Left side of face
                weight = 1.0
                acceleration = sum(H_face(1,iquad,:))
                if(acceleration > 0.0) then
                    weight = H_face_ave(iquad,iface) / acceleration
                end if
                H_face(1,iquad,:) = H_face(1,iquad,:) * weight

                ! Right side of face
                weight = 1.0
                acceleration = sum(H_face(2,iquad,:))
                if(acceleration > 0.0) then
                    weight = H_face_ave(iquad,iface) / acceleration
                end if
                H_face(2,iquad,:) = H_face(2,iquad,:) * weight

            end do ! iquad
            
            ! Do Gauss-Lobatto Integration
            do k = 1,nlayers

                do iquad = 1, nq

                    wq = jac_faceq(iquad,1,iface)
                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    hlx_k = nxl*H_face(1,iquad,k)
                    hrx_k = nxl*H_face(2,iquad,k)
                    hly_k = nyl*H_face(1,iquad,k)
                    hry_k = nyl*H_face(2,iquad,k)
                    flux_x = nxl*udp_flux(1,iquad,k) + nyl*udp_flux(2,iquad,k)
                    flux_y = nxl*vdp_flux(1,iquad,k) + nyl*vdp_flux(2,iquad,k)

                    do n = 1, ngl

                        hi = psiq(n,iquad)
                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)
                        I = intma(il,jl,kl,el)

                        rhs_mom(1,I,k) = rhs_mom(1,I,k) - wq*hi*(hlx_k + flux_x)
                        rhs_mom(2,I,k) = rhs_mom(2,I,k) - wq*hi*(hly_k + flux_y)

                        if(er > 0) then

                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)

                            rhs_mom(1,I,k) = rhs_mom(1,I,k) + wq*hi*(hrx_k + flux_x)
                            rhs_mom(2,I,k) = rhs_mom(2,I,k) + wq*hi*(hry_k + flux_y)

                        end if
                    end do
                end do
            end do
        end do ! iface

    end subroutine Apply_layers_fluxes

    subroutine create_layers_volume_mass(dp_advec, qprime_df)

        use mod_grid, only : npoin_q, npoin, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac
        use mod_variables, only: uvb_ave, ope_ave, sum_layer_mass_flux

        implicit none

        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(npoin,nlayers), intent(out) :: dp_advec

        real :: wq, hi, dhde, dhdn
        integer :: k, I, Iq, ip
        real :: dp_temp, udp, vdp, dp, u, v, opeq
        real, dimension(3) :: qp, qb

        dp_advec = 0.0
        sum_layer_mass_flux = 0.0

        do Iq = 1, npoin_q

            qb(1) = ope_ave(Iq)
            qb(2) = uvb_ave(1,Iq)
            qb(3) = uvb_ave(2,Iq)
            wq = wjac(Iq)

            do k=1,nlayers

                qp = 0.0
                do ip = 1,npts
                    I = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    qp(:) = qp(:) + hi*qprime_df(:,I,k)
                    !qb(:) = qb(:) + hi*uvb_ope_ave_df(:,I)
                enddo

                dp_temp = qp(1) * qb(1)
                udp = (qp(2)+qb(2)) * dp_temp
                vdp = (qp(3)+qb(3)) * dp_temp

                sum_layer_mass_flux(1,Iq) = sum_layer_mass_flux(1,Iq) + udp
                sum_layer_mass_flux(2,Iq) = sum_layer_mass_flux(2,Iq) + vdp

                do ip = 1, npts

                    I = indexq(ip,Iq)
                    dp_advec(I,k) = dp_advec(I,k) + wq*(dpsidx(ip,Iq)*udp + dpsidy(ip,Iq)*vdp)

                end do
            end do

        end do !k

    end subroutine create_layers_volume_mass

    subroutine create_consistency_volume_mass(dp_advec, dprime_df)

        use mod_grid, only : npoin_q, npoin, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, pbprime
        use mod_variables, only: sum_layer_mass_flux, btp_mass_flux_ave

        implicit none 

        real, dimension(npoin,nlayers), intent(in) :: dprime_df
        real, dimension(npoin,nlayers), intent(out) :: dp_advec

        integer :: k, I, Iq, ip
        real :: weight, udp, vdp, wq, hi, dp

        dp_advec = 0.0

        do k=1,nlayers

            do Iq = 1, npoin_q
                dp = 0.0
                do ip = 1,npts
                    I = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    dp = dp + hi*dprime_df(I,k)
                enddo
                weight = dp/pbprime(Iq)

                udp = weight * (btp_mass_flux_ave(1,Iq) - sum_layer_mass_flux(1,Iq))
                vdp = weight * (btp_mass_flux_ave(2,Iq) - sum_layer_mass_flux(2,Iq)) 
                wq = wjac(Iq)

                do ip = 1, npts
                    I = indexq(ip,Iq)
                    dp_advec(I,k) = dp_advec(I,k) + wq*(dpsidx(ip,Iq)*udp + dpsidy(ip,Iq)*vdp)
                end do
            end do
        end do !k

    end subroutine create_consistency_volume_mass

    subroutine create_layer_mass_flux(dp_advec, qprime_df_face)

        use mod_basis, only: nq, psiq, ngl
        use mod_grid, only:  npoin_q, intma,  nface, face
        use mod_input, only: nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: sum_layer_mass_flux_face, ope_face_ave, uvb_face_ave

        implicit none

        real, dimension(npoin, nlayers), intent(inout) :: dp_advec
        real, dimension(3, 2, ngl, nface, nlayers), intent(in) :: qprime_df_face

        integer :: k, iface, iquad, el, er, il, jl, ir, jr, I
        integer :: kl, kr, jquad, n, m
        real :: wq, nxl, nyl, hi
        real, dimension(nq) :: ul,ur,vl,vr, flux_edge_u, flux_edge_v
        real :: dpl, dpr, uu, vv, flux
        real, dimension(3,nq) :: ql, qr, qbl, qbr

        sum_layer_mass_flux_face = 0.0

        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)
            
            qbl(1,:) = ope_face_ave(1,:,iface)
            qbl(2,:) = uvb_face_ave(1,1,:,iface)
            qbl(3,:) = uvb_face_ave(2,1,:,iface)
            qbr(1,:) = ope_face_ave(2,:,iface)
            qbr(2,:) = uvb_face_ave(1,2,:,iface)
            qbr(3,:) = uvb_face_ave(2,2,:,iface)

            do k = 1, nlayers

                do iquad = 1,nq

                    ! In the following computation of fluxes at element faces,
                    ! flux_edge iface a numerical approximation to the mass flux at element faces (each face has nq quadrature points)
                    ! Here we are using centered fluxes, so we need to compute the fluxes at the left and right edges of each face

                    ql = 0.0; qr = 0.0
                    do n = 1, ngl
                        hi = psiq(n,iquad)
                        ql(:,iquad) = ql(:,iquad) + hi*qprime_df_face(:,1,n,iface,k)
                        qr(:,iquad) = qr(:,iquad) + hi*qprime_df_face(:,2,n,iface,k)
                        !qbl(:,iquad) = qbl(:,iquad) + hi*uvb_ope_face_ave_df(:,1,n,iface)
                        !qbr(:,iquad) = qbr(:,iquad) + hi*uvb_ope_face_ave_df(:,2,n,iface)
                    enddo

                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    uu = 0.5*((ql(2,iquad)+qbl(2,iquad)) + (qr(2,iquad)+qbr(2,iquad)))
                    vv = 0.5*((ql(3,iquad)+qbl(3,iquad)) + (qr(3,iquad)+qbr(3,iquad)))

                    dpl = qbl(1,iquad) * ql(1,iquad)
                    dpr = qbr(1,iquad) * qr(1,iquad)

                    if(uu*nxl > 0.0) then
                        flux_edge_u(iquad) = uu * dpl
                    else
                        flux_edge_u(iquad) = uu * dpr
                    endif

                    if(vv*nyl > 0.0) then
                        flux_edge_v(iquad) = vv * dpl
                    else
                        flux_edge_v(iquad) = vv * dpr
                    endif
                end do

                sum_layer_mass_flux_face(1,:,iface) = sum_layer_mass_flux_face(1,:,iface) + flux_edge_u(:)
                sum_layer_mass_flux_face(2,:,iface) = sum_layer_mass_flux_face(2,:,iface) + flux_edge_v(:)

                do iquad = 1, nq

                    wq = jac_faceq(iquad,1,iface)
                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    flux = nxl*flux_edge_u(iquad) + nyl*flux_edge_v(iquad)

                    do n = 1, ngl

                        hi = psiq(n,iquad)
                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)
                        I = intma(il,jl,kl,el)

                        dp_advec(I,k) = dp_advec(I,k) - wq*hi*flux

                        if(er > 0) then

                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)

                            dp_advec(I,k) = dp_advec(I,k) + wq*hi*flux

                        end if
                    end do
                end do
            end do
        end do

    end subroutine create_layer_mass_flux

    subroutine create_consistency_mass_flux(dp_advec, flux_deficit_mass_face)

        use mod_basis, only: nq, psiq, ngl
        use mod_grid, only:  npoin_q, intma,  nface, face
        use mod_input, only: nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: uvb_face_ave, ope_face_ave, sum_layer_mass_flux_face
        use mod_initial, only: pbprime_face

        implicit none

        real, dimension(npoin, nlayers), intent(inout) :: dp_advec
        real, dimension(2,2,nq,nface,nlayers), intent(in)  :: flux_deficit_mass_face
    
        integer :: k, iface, iquad, el, er, il, jl, ir, jr, I
        integer :: kl, kr, jquad, n, m
        real :: wq, nxl, nyl, hi, flux
        real, dimension(nq) :: flux_edge_u, flux_edge_v
    
        do k = 1, nlayers
            do iface = 1, nface

                el = face(7,iface)
                er = face(8,iface)

                do iquad = 1,nq 

                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    if(flux_deficit_mass_face(1,1,iquad,iface,k)*nxl > 0.0) then 

                        flux_edge_u(iquad) = flux_deficit_mass_face(1,1,iquad,iface,k)
                    else 
                        flux_edge_u(iquad) = flux_deficit_mass_face(1,2,iquad,iface,k)
                    end if 

                    if(flux_deficit_mass_face(2,1,iquad,iface,k)*nyl > 0.0) then 
                        flux_edge_v(iquad) = flux_deficit_mass_face(2,1,iquad,iface,k)
                    else 
                        flux_edge_v(iquad) = flux_deficit_mass_face(2,2,iquad,iface,k)
                    end if 
                enddo 

                do iquad = 1, nq

                    wq = jac_faceq(iquad,1,iface)

                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    flux = nxl*flux_edge_u(iquad) + nyl*flux_edge_v(iquad)

                    do n = 1, ngl

                        hi = psiq(n,iquad)
                        
                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)
                        I = intma(il,jl,kl,el)

                        dp_advec(I,k) = dp_advec(I,k) - wq*hi*flux

                        if(er > 0) then

                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)

                            dp_advec(I,k) = dp_advec(I,k) + wq*hi*flux

                        end if
                    end do
                end do
            end do
        end do
    
    end subroutine create_consistency_mass_flux

end module mod_create_rhs_mlswe
