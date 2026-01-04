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

    public :: layer_momentum_rhs, &
              interpolate_layer_from_quad_to_node, rhs_layer_shear_stress, layer_mass_rhs, &
              bcl_rhs
              
    contains

    subroutine layer_momentum_rhs(rhs_mom, rhs_visc, qprime_df, q_df)

        use mod_input, only: nlayers
        use mod_metrics, only: massinv
        use mod_grid, only: npoin, npoin_q, intma, intma_dg_quad, nface
        use mod_basis, only: nq, ngl

        implicit none
        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom
        real, dimension(2,npoin,nlayers), intent(in) :: rhs_visc
        real, dimension(3,npoin,nlayers), intent(in) :: q_df, qprime_df

        integer :: k,I

        call bcl_create_precommunicator(qprime_df)

        call create_rhs_dynamics_volume_layers(rhs_mom, qprime_df, q_df)

        call Apply_layers_fluxes(rhs_mom, qprime_df)

        call bcl_create_postcommunicator_momentum(rhs_mom)

        do k = 1, nlayers
            rhs_mom(1,:,k) = massinv(:)*rhs_mom(1,:,k) + rhs_visc(1,:,k)
            rhs_mom(2,:,k) = massinv(:)*rhs_mom(2,:,k) + rhs_visc(2,:,k)
        end do

    end subroutine layer_momentum_rhs

    subroutine layer_mass_rhs(dp_advec, qprime_df)

        use mod_metrics, only: massinv
        use mod_input, only: nlayers
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq, ngl

        implicit none

        real, dimension(npoin, nlayers), intent(out) :: dp_advec
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        
        integer :: k

        call bcl_create_precommunicator(qprime_df)
        ! Compute the mass advection term for the degree of freedom for dp in each layer
        call create_layers_volume_mass(dp_advec, qprime_df)

        ! Compute the mass flux term 
        call create_layer_mass_flux(dp_advec, qprime_df)

        call bcl_create_postcommunicator_continuity(dp_advec)

        do k = 1, nlayers
            dp_advec(:,k) = massinv(:)*dp_advec(:,k)
        end do
        
    end subroutine layer_mass_rhs

    subroutine bcl_rhs(rhs, rhs_visc, qprime_df, q_df)

        use mod_input, only: nlayers
        use mod_metrics, only: massinv
        use mod_grid, only: npoin, npoin_q, intma, intma_dg_quad, nface
        use mod_basis, only: nq, ngl

        implicit none
        real, dimension(3,npoin,nlayers), intent(out) :: rhs
        real, dimension(2,npoin,nlayers), intent(in) :: rhs_visc
        real, dimension(3,npoin,nlayers), intent(in) :: q_df, qprime_df

        integer :: k,I

        call bcl_create_precommunicator(qprime_df)

        call create_rhs_dynamics_volume_bcl(rhs, qprime_df, q_df)
        call Apply_bcl_fluxes(rhs, qprime_df)

        call bcl_create_postcommunicator(rhs)

        do k = 1, nlayers
            rhs(1,:,k) = massinv(:)*rhs(1,:,k)
            rhs(2,:,k) = massinv(:)*rhs(2,:,k) + rhs_visc(1,:,k)
            rhs(3,:,k) = massinv(:)*rhs(3,:,k) + rhs_visc(2,:,k)
        end do

    end subroutine bcl_rhs

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

    subroutine rhs_layer_shear_stress(rhs_stress,q_df)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_input, only: nlayers, ad_mlswe, max_shear_dz, dt
        use mod_initial, only: coriolis_quad, alpha_mlswe
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac

        implicit none
        
        real, intent(in) :: q_df(3,npoin,nlayers)
        real, intent(out) :: rhs_stress(2,npoin,nlayers)
        
        real, dimension(nlayers+1) :: tau_u, tau_v
        real, dimension(nlayers)   :: a,b,c, dp, udp, vdp
        real, dimension(2,nlayers) :: r, uv
        real :: wq, hi, tau_u_q, tau_v_q, coeff, mult, coeff1
        integer :: k, e, iquad, jquad, kquad, l, m, n, Iq, I, ip

        ! Use shear stress between layers to compute  u  and  v  in each layer. 
        ! In order to do this, implement the shear stress implicitly and
        ! solve a linear system in the vertical direction, at each
        ! horizontal location.   
        ! This computation assumes that the shear stress at an interface
        ! between layers has the form 
        ! tau_u = rho * A_D * (u_upper - u_lower) / Delta z
        ! tau_v = rho * A_D * (v_upper - v_lower) / Delta z
        ! Here,  Delta z  is the thickness of the Ekman layer 
        ! (i.e., sqrt(2*A_D/f)) or the value of  max_shear_dz,  whichever 
        ! is less.
        
        rhs_stress = 0.0
        
        ! Compute the shear stress
        do Iq = 1,npoin_q

            ! Projection from nodal to quadrature points
            dp = 0.0 ; udp = 0.0 ; vdp = 0.0
            do ip = 1,npts
                I = indexq(ip,Iq)
                hi = psih(ip,Iq)
                dp(:) = dp(:) + hi*q_df(1,I,:)
                udp(:) = udp(:) + hi*q_df(2,I,:)
                vdp(:) = vdp(:) + hi*q_df(3,I,:)
            enddo

            ! The linear system is structured as follows.
            ! a  refers to the subdiagonal.
            ! b  refers to the diagonal.
            ! c  refers to the superdiagonal.
            ! r(1,:)  initially contains the right side for the system
            ! for  u, and at the end it contains the solution, which consists of 
            ! values of  u. 
            ! r(2,:)  initially contains the right side for the system
            ! for  v, and at the end it contains the solution, which consists of
            ! values of  v.

            coeff = max(sqrt(0.5*coriolis_quad(Iq)*ad_mlswe)/alpha_mlswe(1), &
                                  ad_mlswe/(alpha_mlswe(1) * max_shear_dz))

            coeff1 = gravity*dt*coeff

            ! Set up the system
            do k = 1, nlayers
                a(k) = -coeff
                b(k) = dp(k) + 2.0*coeff1
                c(k) = -coeff1
                r(1,k) = udp(k)/dp(k)
                r(2,k) = vdp(k)/dp(k)
            end do
            
            b(1) = dp(1) + coeff1
            b(nlayers) = dp(nlayers) + coeff1
            a(1) = 0.0
            c(nlayers) = 0.0
            
            ! Solve the system
            do k = 2, nlayers
                mult = a(k) / b(k-1)
                b(k) = b(k) - mult*c(k-1)
                r(1,k) = r(1,k) - mult*r(1,k-1)
                r(2,k) = r(2,k) - mult*r(2,k-1)
            end do
            
            r(1,nlayers) = r(1,nlayers) / b(nlayers)
            r(2,nlayers) = r(2,nlayers) / b(nlayers)
            uv(1,nlayers) = r(1,nlayers)
            uv(2,nlayers) = r(2,nlayers)
            
            do k = nlayers-1,1,-1
                r(1,k) = (r(1,k) - c(k)*r(1,k+1)) / b(k)
                r(2,k) = (r(2,k) - c(k)*r(2,k+1)) / b(k)
                uv(1,k) = r(1,k)
                uv(2,k) = r(2,k)
            end do

            ! tau_u(k)  and  tau_v(k) contain shear stresses along interface k,
            ! which is the bottom of layer k.

            tau_u(1) = 0.0 ; tau_v(1) = 0.0

            do k=2,nlayers
                tau_u(k) = coeff*(uv(1,k-1) - uv(1,k))
                tau_v(k) = coeff*(uv(2,k-1) - uv(2,k))
            end do

            ! Do Gauss-Lobatto Integration
            wq = wjac(Iq)

            do k = 1, nlayers

                tau_u_q = gravity*(tau_u(k) - tau_u(k+1))
                tau_v_q = gravity*(tau_v(k) - tau_v(k+1))
                
                do ip = 1, npts

                    I = indexq(ip,Iq)
                    hi = psih(ip,Iq)

                    rhs_stress(1,I,k) = rhs_stress(1,I,k) + wq*hi*tau_u_q
                    rhs_stress(2,I,k) = rhs_stress(2,I,k) + wq*hi*tau_v_q
                end do
            end do

        end do ! Iq

        ! Apply the mass inverse to the shear stress
        ! do k = 1, nlayers
        !     rhs_stress(1,:,k) = massinv(:)*rhs_stress(1,:,k)
        !     rhs_stress(2,:,k) = massinv(:)*rhs_stress(2,:,k)
        ! end do
        
    end subroutine rhs_layer_shear_stress

    subroutine create_rhs_dynamics_volume_layers(rhs_mom, qprime_df,q_df)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts, dpsiqx
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, alpha_mlswe, &
                                tau_wind, pbprime_df, zbot_df
        use mod_variables, only: tau_bot_ave, ope_ave, uvb_ave, H_ave, &
                                    Qu_ave, Qv_ave, Quv_ave, uvb_ave, ope2_ave_df, ope2_ave

        implicit none

        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom
        real, dimension(3,npoin,nlayers), intent(in) :: q_df, qprime_df

        real :: wq, hi, dhdx, dhdy, bot_layer, tau_wind_u, tau_wind_v, temp1
        real :: Hq, var_uu, var_uv, var_vu, var_vv, source_x, source_y, Pstress
        integer :: k, I, Iq, ip
        real, dimension(nlayers+1) :: pprime_temp
        real :: temp_dp, temp_u, temp_v, Ptop_k, Pbot_k, tempbot, Pbstress
        real :: weight, acceleration, pbq
        real, dimension(3) :: qp, qb
        real, dimension(nlayers) :: temp_uu, temp_vv, H_tmp, u_udp, v_vdp
        real, dimension(2,nlayers) :: u_vdp
        real :: p_tmp(nlayers+1), u, v, dp, weightq, one_over_sumuq, one_over_sumvq
        real :: uu_dp_deficitq, uv_dp_deficitq, vv_dp_deficitq, gradz(2,nlayers+1)
        real :: z_elv(npoin,nlayers+1)
        real, parameter :: eps1 = 1.0e-20, eps = 1.0e-10 !  Parameter used to prevent division by zero.

        rhs_mom = 0.0
        bot_layer = 0.0
        pprime_temp = 0.0

        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth 
                                                  ! at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth 
                                                         ! at which bottom stress is reduced to 0

        ! Find layer interfaces
        z_elv(:,nlayers+1) = zbot_df(:)
        do k = nlayers,1,-1
            z_elv(:,k) = z_elv(:,k+1) + (alpha_mlswe(k)/gravity) * &
                                        (sqrt(ope2_ave_df(:))*qprime_df(1,:,k))
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

                    temp_uu(k) = temp_uu(k) + hi*q_df(2,I,k)
                    temp_vv(k) = temp_vv(k) + hi*q_df(3,I,k)
                enddo

                if (qp(1) <= (gravity/alpha_mlswe(k))*eps) then
                    qp(1) = (gravity/alpha_mlswe(k))*eps
                    qp(2:3) = 0.0
                    temp_uu(k) = 0.0 ; temp_vv(k) = 0.0
                end if

                qb(1) = ope_ave(Iq)
                qb(2) = uvb_ave(1,Iq)
                qb(3) = uvb_ave(2,Iq)

                p_tmp(k+1) = p_tmp(k) + sqrt(ope2_ave(Iq)) * qp(1)
                H_tmp(k) = 0.5*alpha_mlswe(k) * (p_tmp(k+1)**2 - p_tmp(k)**2)

                dp = qp(1) * qb(1)
                u = qp(2) + qb(2)
                v = qp(3) + qb(3)

                u_udp(k) = dp * u*u
                v_vdp(k) = dp * v*v
                u_vdp(1,k) = u * v * dp
                u_vdp(2,k) = v * u * dp

                temp_uu(k) = abs(temp_uu(k)) + eps1
                temp_vv(k) = abs(temp_vv(k)) + eps1
            end do

            gradz = 0.0 ; pbq = 0.0
            do ip = 1,npts
                I = indexq(ip,Iq)
                gradz(1,:) = gradz(1,:) + dpsidx(ip,Iq)*z_elv(I,:)
                gradz(2,:) = gradz(2,:) + dpsidy(ip,Iq)*z_elv(I,:)
                pbq = pbq + psih(ip,Iq)*pbprime_df(I)
            enddo

            ! Consistency
            uu_dp_deficitq = Qu_ave(1,Iq) - sum(u_udp(:))
            uv_dp_deficitq = Qu_ave(2,Iq) - sum(u_vdp(1,:))
            vv_dp_deficitq = Qv_ave(2,Iq) - sum(v_vdp(:))

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
                ! H_r  over all layers equals the time average of the barotropic forcing  H  
                ! over all barotropic substeps of the baroclinic
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

                tempbot = min(Pbstress,pbq-pprime_temp(k+1)) - min(Pbstress, pbq &
                            - pprime_temp(k))
                tempbot = tempbot / Pbstress

                source_x = gravity*(tau_wind_u - tempbot*tau_bot_ave(1,Iq) + &
                            p_tmp(k) * gradz(1,k) - p_tmp(k+1) * gradz(1,k+1))
                source_y = gravity*(tau_wind_v - tempbot*tau_bot_ave(2,Iq) + &
                            p_tmp(k) * gradz(2,k) - p_tmp(k+1) * gradz(2,k+1))

                ! Do Gauss-Lobatto Integration
                do ip = 1, npts

                    I = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    !Xi derivatives
                    dhdx = dpsidx(ip,Iq)
                    !Eta derivatives
                    dhdy = dpsidy(ip,Iq)

                    rhs_mom(1,I,k) = rhs_mom(1,I,k) + wq*(hi*source_x &
                                                    + dhdx*(Hq + var_uu) + var_uv*dhdy)
                    rhs_mom(2,I,k) = rhs_mom(2,I,k) + wq*(hi*source_y &
                                                    + var_vu*dhdx + dhdy*(Hq +var_vv))

                end do
            end do
        end do !k

    end subroutine create_rhs_dynamics_volume_layers

    subroutine Apply_layers_fluxes(rhs_mom, qprime_df)

        ! This routine computes the layer momentum advection flux terms using upwind flux
        ! This routine computes the layer momentum pressure terms

        use mod_constants, only : gravity
        use mod_initial, only : alpha_mlswe, zbot_face
        use mod_grid, only : nface, npoin, npoin_q, face, intma, face_type
        use mod_basis, only : nq, psiq, ngl
        use mod_input, only : nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: ope_face_ave, H_face_ave, one_plus_eta_edge_2_ave, &
                                uvb_face_ave, Quv_face_ave, Qu_face_ave, Qv_face_ave, ope2_face_ave

        implicit none

        real, dimension(2, npoin, nlayers), intent(inout) :: rhs_mom
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

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
        real ::  ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv, hi, wq
        real, dimension(nq,nlayers) :: udpl, udpr, vdpl, vdpr
        real :: uu_dp_flux_deficit(2), vv_dp_flux_deficit(2)
        real, parameter :: eps1 = 1.0e-20, eps = 1.0e-10 !  Parameter used to prevent division by zero.
        integer :: il, jl, ir,jr, kl, kr, n, m, hl_k, hr_k, flux_x, flux_y, un
        real, dimension(2,nq,nlayers) :: H_face, udp_flux, vdp_flux

        do k=1,nlayers
            alpha_over_g(k) = alpha_mlswe(k)/gravity
            g_over_alpha(k) = gravity/alpha_mlswe(k)
        enddo

        ! Compute H_r at the element face
        do iface = 1, nface

            if (face_type(iface) == 2) cycle

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

                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)
                        I = intma(il,jl,kl,el)
                        hi = psiq(n,iquad)

                        ql(:,iquad,k) = ql(:,iquad,k) + hi*qprime_df(:,I,k)
                    enddo

                    if (er > 0) then
                        do n = 1, ngl

                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)
                            hi = psiq(n,iquad)

                            qr(:,iquad,k) = qr(:,iquad,k) + hi*qprime_df(:,I,k)
                        enddo
                    else 
                        qr(:,iquad,k) = ql(:,iquad,k)

                    endif 

                    if (ql(1,iquad,k) <= (gravity/alpha_mlswe(k))*eps) then
                        ql(1,iquad,k) = (gravity/alpha_mlswe(k))*eps
                        ql(2:3,iquad,k) = 0.0
                    end if
                    if (qr(1,iquad,k) <= (gravity/alpha_mlswe(k))*eps) then
                        qr(1,iquad,k) = (gravity/alpha_mlswe(k))*eps
                        qr(2:3,iquad,k) = 0.0
                    end if

                    ! Apply wall boundary conditions
                    if (er == -4) then
                        un = ql(2,iquad,k)*nxl + ql(3,iquad,k)*nyl
                        qr(2,iquad,k) = ql(2,iquad,k) - 2.0*un*nxl
                        qr(3,iquad,k) = ql(3,iquad,k) - 2.0*un*nyl
                    elseif (er == -2) then
                        qr(2,iquad,k) = -ql(2,iquad,k)
                        qr(3,iquad,k) = -ql(3,iquad,k)
                    endif

                    ! Left side of the edge
                    dpl = qbl(1,iquad) * ql(1,iquad,k)
                    dpr = qbr(1,iquad) * qr(1,iquad,k)
                    ul = ql(2,iquad,k) + qbl(2,iquad)
                    ur = qr(2,iquad,k) + qbr(2,iquad)
                    vl = ql(3,iquad,k) + qbl(3,iquad)
                    vr = qr(3,iquad,k) + qbr(3,iquad)

                    uu = 0.5*(ul+ur)
                    vv = 0.5*(vl+vr)
                    udpl(iquad,k) = ul*dpl
                    udpr(iquad,k) = ur*dpr
                    vdpl(iquad,k) = vl*dpl
                    vdpr(iquad,k) = vr*dpr

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

                uu_dp_flux_deficit(1) = Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,iquad,:))
                uu_dp_flux_deficit(2) = Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,iquad,:))
                vv_dp_flux_deficit(1) = Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,iquad,:))
                vv_dp_flux_deficit(2) = Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,iquad,:))

                ! Adjust the fluxes for the u-momentum equation
                do k = 1,nlayers

                    ! Adjust the fluxes for the u-momentum equation
                    !x-direction
                    weight = abs(udpl(iquad,k)) / sum(abs(udpl(iquad,:))+eps1)
                    if(uu_dp_flux_deficit(1)*nxl < 0.0) &
                        weight = abs(udpr(iquad,k)) / sum(abs(udpr(iquad,:))+eps1)

                    udp_flux(1,iquad,k) = udp_flux(1,iquad,k) + weight * uu_dp_flux_deficit(1)

                    !y-direction
                    weight = abs(udpl(iquad,k)) / sum(abs(udpl(iquad,:))+eps1)
                    if(uu_dp_flux_deficit(2)*nyl < 0.0) &
                        weight = abs(udpr(iquad,k)) / sum(abs(udpr(iquad,:))+eps1)

                    udp_flux(2,iquad,k) = udp_flux(2,iquad,k) + weight * uu_dp_flux_deficit(2)
                
                    ! Adjust the fluxes for the v-momentum equation
                    !x-direction
                    weight = abs(vdpl(iquad,k)) / sum(abs(vdpl(iquad,:))+eps1)
                    if(vv_dp_flux_deficit(1)*nxl < 0.0) &
                        weight = abs(vdpr(iquad,k)) / sum(abs(vdpr(iquad,:))+eps1)

                    vdp_flux(1,iquad,k) = vdp_flux(1,iquad,k) + weight * vv_dp_flux_deficit(1)
  
                    !y-direction
                    weight = abs(vdpl(iquad,k)) / sum(abs(vdpl(iquad,:))+eps1)
                    if(vv_dp_flux_deficit(2)*nyl < 0.0) &
                        weight = abs(vdpr(iquad,k)) / sum(abs(vdpr(iquad,:))+eps1)

                    vdp_flux(2,iquad,k) = vdp_flux(2,iquad,k) + weight * vv_dp_flux_deficit(2)
                end do

                z_face = 0.0 ; p_face = 0.0
                z_edge_plus = 0.0 ; z_edge_minus = 0.0
                p_edge_plus = 0.0; p_edge_minus = 0.0

                !Store Left Side Variables
                ope_l = sqrt(ope2_face_ave(1,iquad,iface))
                ope_r = sqrt(ope2_face_ave(2,iquad,iface))
                p_face(1,1) = 0.0 ; p_face(2,1) = 0.0
                do k=1,nlayers
                    p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,iquad,k)
                    p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,iquad,k)
                end do

                one_plus_eta_edge = sqrt(one_plus_eta_edge_2_ave(iquad,iface))
                z_face(1,nlayers+1) = zbot_face(1,iquad,iface)
                z_face(2,nlayers+1) = zbot_face(2,iquad,iface)
                z_edge_plus(nlayers+1) = zbot_face(1,iquad,iface)
                z_edge_minus(nlayers+1) = zbot_face(2,iquad,iface)
                do k=nlayers,1,-1
                    z_face(1,k) = z_face(1,k+1) + alpha_over_g(k) * (ope_l * ql(1,iquad,k))
                    z_face(2,k) = z_face(2,k+1) + alpha_over_g(k) * (ope_r * qr(1,iquad,k))
                    z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k) * &
                                                        (one_plus_eta_edge * ql(1,iquad,k))
                    z_edge_minus(k) = z_edge_minus(k+1) + alpha_over_g(k) * &
                                                        (one_plus_eta_edge * qr(1,iquad,k))
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
                            p_intersect_bot = p_edge_minus(ktemp+1) &
                                    - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_minus(ktemp+1))
                            p_intersect_top = p_edge_minus(ktemp+1) &
                                    - g_over_alpha(ktemp)*(z_intersect_top - z_edge_minus(ktemp+1))
                            H_r_minus = H_r_minus + &
                                    0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)

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
                            p_intersect_bot = p_edge_plus(ktemp+1) &
                                - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_plus(ktemp+1))
                            p_intersect_top = p_edge_plus(ktemp+1) &
                                - g_over_alpha(ktemp)*(z_intersect_top - z_edge_plus(ktemp+1))
                            H_r_plus = H_r_plus &
                                + 0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
                        end if
                    end do
                    H_face(2,iquad,k) = 0.5*(H_r_plus + H_r_minus) ! H_r for the right side
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
                        H_corr1 = 0.5 * alpha_mlswe(k) * ((p_face(1,k+1) &
                                    + p_inc1)**2 - p_face(1,k+1)**2)
                        H_face(1,iquad,k) = H_face(1,iquad,k) - H_corr1
                        H_face(1,iquad,k+1) = H_face(1,iquad,k+1) + H_corr1

                        ! Corrections at the right side of a face.
                        p_inc2 = g_over_alpha(k)*(z_face(2,k+1) - z_edge_minus(k+1))
                        H_corr2 = 0.5 * alpha_mlswe(k) * ((p_face(2,k+1) &
                                    + p_inc2)**2 - p_face(2,k+1)**2)
                        H_face(2,iquad,k) = H_face(2,iquad,k) - H_corr2
                        H_face(2,iquad,k+1) = H_face(2,iquad,k+1) + H_corr2

                    end do
                end if

                ! Adjust the values of  H_k, at element faces, so that the vertical sum of
                ! H_k  over all layers equals the time average of the barotropic forcing  H 
                ! over all barotropic substeps of the baroclinic
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
                if(acceleration > 0.0) weight = H_face_ave(iquad,iface) / acceleration
                H_face(1,iquad,:) = H_face(1,iquad,:) * weight

                ! Right side of face
                weight = 1.0
                acceleration = sum(H_face(2,iquad,:))
                if(acceleration > 0.0) weight = H_face_ave(iquad,iface) / acceleration
                H_face(2,iquad,:) = H_face(2,iquad,:) * weight

            end do ! iquad
            
            ! Do Gauss-Lobatto Integration
            do k = 1,nlayers

                do iquad = 1, nq

                    wq = jac_faceq(iquad,1,iface)
                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    hl_k = H_face(1,iquad,k)
                    hr_k = H_face(2,iquad,k)
                    flux_x = nxl*udp_flux(1,iquad,k) + nyl*udp_flux(2,iquad,k)
                    flux_y = nxl*vdp_flux(1,iquad,k) + nyl*vdp_flux(2,iquad,k)

                    do n = 1, ngl

                        hi = psiq(n,iquad)
                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)
                        I = intma(il,jl,kl,el)

                        rhs_mom(1,I,k) = rhs_mom(1,I,k) - wq*hi*(nxl*hl_k + flux_x)
                        rhs_mom(2,I,k) = rhs_mom(2,I,k) - wq*hi*(nyl*hl_k + flux_y)

                        if(er > 0) then

                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)

                            rhs_mom(1,I,k) = rhs_mom(1,I,k) + wq*hi*(nxl*hr_k + flux_x)
                            rhs_mom(2,I,k) = rhs_mom(2,I,k) + wq*hi*(nyl*hr_k + flux_y)

                        end if
                    end do
                end do
            end do
        end do ! iface

    end subroutine Apply_layers_fluxes

    subroutine create_rhs_dynamics_volume_bcl(rhs, qprime_df, q_df)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts, dpsiqx
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, alpha_mlswe, &
                                tau_wind, pbprime_df, zbot_df, coriolis_quad
        use mod_variables, only: tau_bot_ave, ope_ave, uvb_ave, H_ave, &
                                    Qu_ave, Qv_ave, Quv_ave, uvb_ave, &
                                    ope2_ave_df, ope2_ave, sum_layer_mass_flux, &
                                    btp_mass_flux_ave

        implicit none

        real, dimension(3,npoin,nlayers), intent(out) :: rhs
        real, dimension(3,npoin,nlayers), intent(in) :: q_df, qprime_df

        real :: wq, hi, dhdx, dhdy, bot_layer, tau_wind_u, tau_wind_v, temp1
        real :: Hq, var_uu, var_uv, var_vu, var_vv, source_x, source_y, Pstress
        integer :: k, I, Iq, ip
        real, dimension(nlayers+1) :: pprime_temp
        real :: temp_dp, temp_u, temp_v, Ptop_k, Pbot_k, tempbot, Pbstress
        real :: weight, acceleration, pbq
        real, dimension(3) :: qp, qb
        real, dimension(nlayers) :: temp_uu, temp_vv, H_tmp, u_udp, v_vdp, udp, vdp, dp, dpp
        real, dimension(2,nlayers) :: u_vdp
        real :: p_tmp(nlayers+1), u, v, weightq, one_over_sumuq, one_over_sumvq
        real :: uu_dp_deficitq(2), vv_dp_deficitq(2), gradz(2,nlayers+1)
        real :: z_elv(npoin,nlayers+1)
        real, parameter :: eps1 = 1.0e-20, eps = 1.0e-10 !  Parameter used to prevent division by zero.
        real :: flux(2,3), udp_cstcy, vdp_cstcy, weight_dp, sum_layer_mass(2)

        rhs = 0.0
        bot_layer = 0.0
        pprime_temp = 0.0
        sum_layer_mass_flux = 0.0

        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth 
                                                  ! at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth 
                                                         ! at which bottom stress is reduced to 0

        ! Find layer interfaces
        z_elv(:,nlayers+1) = zbot_df(:)
        do k = nlayers,1,-1
            z_elv(:,k) = z_elv(:,k+1) + (alpha_mlswe(k)/gravity) * &
                                        (sqrt(ope2_ave_df(:))*qprime_df(1,:,k))
        end do

        do Iq = 1, npoin_q

            p_tmp(1) = 0.0
            temp_uu = 0.0; temp_vv = 0.0
            sum_layer_mass = 0.0
            do k = 1,nlayers
                qp = 0.0; qb = 0.0
                do ip = 1,npts
                    I = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    qp(:) = qp(:) + hi*qprime_df(:,I,k)

                    temp_uu(k) = temp_uu(k) + hi*q_df(2,I,k)
                    temp_vv(k) = temp_vv(k) + hi*q_df(3,I,k)
                enddo

                if (qp(1) <= (gravity/alpha_mlswe(k))*eps) then
                    qp(1) = (gravity/alpha_mlswe(k))*eps
                    qp(2:3) = 0.0
                    temp_uu(k) = 0.0 ; temp_vv(k) = 0.0
                end if

                qb(1) = ope_ave(Iq)
                qb(2) = uvb_ave(1,Iq)
                qb(3) = uvb_ave(2,Iq)

                p_tmp(k+1) = p_tmp(k) + sqrt(ope2_ave(Iq)) * qp(1)
                H_tmp(k) = 0.5*alpha_mlswe(k) * (p_tmp(k+1)**2 - p_tmp(k)**2)

                dpp(k) = qp(1)
                dp(k) = qp(1) * qb(1)
                u = qp(2) + qb(2)
                v = qp(3) + qb(3)

                udp(k) = u*dp(k)
                vdp(k) = v*dp(k)

                u_udp(k) = u * udp(k)
                v_vdp(k) = v * vdp(k)
                u_vdp(1,k) = v * udp(k)
                u_vdp(2,k) = u * vdp(k)

                temp_uu(k) = abs(temp_uu(k)) + eps1
                temp_vv(k) = abs(temp_vv(k)) + eps1

                sum_layer_mass(1) = sum_layer_mass(1) + udp(k)
                sum_layer_mass(2) = sum_layer_mass(2) + vdp(k)
            end do

            gradz = 0.0 ; pbq = 0.0
            do ip = 1,npts
                I = indexq(ip,Iq)
                gradz(1,:) = gradz(1,:) + dpsidx(ip,Iq)*z_elv(I,:)
                gradz(2,:) = gradz(2,:) + dpsidy(ip,Iq)*z_elv(I,:)
                pbq = pbq + psih(ip,Iq)*pbprime_df(I)
            enddo

            ! Consistency
            uu_dp_deficitq(1) = Qu_ave(1,Iq) - sum(u_udp(:))
            uu_dp_deficitq(2) = Qu_ave(2,Iq) - sum(u_vdp(1,:))
            vv_dp_deficitq(1) = Qv_ave(1,Iq) - sum(u_vdp(2,:))
            vv_dp_deficitq(2) = Qv_ave(2,Iq) - sum(v_vdp(:))

            one_over_sumuq = 1.0/sum(temp_uu(:))
            one_over_sumvq = 1.0/sum(temp_vv(:))

            wq = wjac(Iq)
            pprime_temp = 0.0

            do k = 1,nlayers

                weight_dp = dpp(k) / pbq
                udp_cstcy = weight_dp * (btp_mass_flux_ave(1,Iq) - sum_layer_mass(1))
                vdp_cstcy = weight_dp * (btp_mass_flux_ave(2,Iq) - sum_layer_mass(2))

                pprime_temp(k+1) = pprime_temp(k) + qp(k)

                weightq = temp_uu(k) * one_over_sumuq
                u_udp(k) = u_udp(k) + weightq * uu_dp_deficitq(1)
                u_vdp(1,k) = u_vdp(1,k) + weightq * uu_dp_deficitq(2)

                weightq = temp_vv(k) * one_over_sumvq
                u_vdp(2,k) = u_vdp(2,k) + weightq * vv_dp_deficitq(1)
                v_vdp(k) = v_vdp(k) + weightq * vv_dp_deficitq(2)

                Hq = H_tmp(k)
                
                ! Adjust the values of  H_r, at quadrature points, so that the vertical sum of
                ! H_r  over all layers equals the time average of the barotropic forcing  H  
                ! over all barotropic substeps of the baroclinic
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

                flux(1,1) = udp(k) + udp_cstcy
                flux(2,1) = vdp(k) + vdp_cstcy

                flux(1,2) = u_udp(k) + Hq
                flux(2,2) = u_vdp(1,k)

                flux(1,3) = u_vdp(2,k)
                flux(2,3) = v_vdp(k) + Hq

                temp1 = (min(pprime_temp(k+1), Pstress) - min(pprime_temp(k), Pstress))/ Pstress
                tau_wind_u = temp1*tau_wind(1,Iq)
                tau_wind_v = temp1*tau_wind(2,Iq)

                tempbot = min(Pbstress,pbq-pprime_temp(k+1)) - min(Pbstress, pbq &
                            - pprime_temp(k))
                tempbot = tempbot / Pbstress

                source_x = gravity*(tau_wind_u - tempbot*tau_bot_ave(1,Iq) + &
                            p_tmp(k) * gradz(1,k) - p_tmp(k+1) * gradz(1,k+1))
                source_y = gravity*(tau_wind_v - tempbot*tau_bot_ave(2,Iq) + &
                            p_tmp(k) * gradz(2,k) - p_tmp(k+1) * gradz(2,k+1))

                source_x = source_x + coriolis_quad(Iq)*vdp(k)
                source_y = source_y - coriolis_quad(Iq)*udp(k)

                ! Do Gauss-Lobatto Integration
                do ip = 1, npts

                    I = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    !Xi derivatives
                    dhdx = dpsidx(ip,Iq)
                    !Eta derivatives
                    dhdy = dpsidy(ip,Iq)

                    rhs(1,I,k) = rhs(1,I,k) + wq*(dhdx*flux(1,1) + dhdy*flux(2,1))

                    rhs(2,I,k) = rhs(2,I,k) + wq*(hi*source_x &
                                                    + dhdx*flux(1,2) + flux(2,2)*dhdy)
                    rhs(3,I,k) = rhs(3,I,k) + wq*(hi*source_y &
                                                    + flux(1,3)*dhdx + dhdy*flux(2,3))

                end do
            end do
        end do !k

    end subroutine create_rhs_dynamics_volume_bcl

    subroutine Apply_bcl_fluxes(rhs, qprime_df)

        ! This routine computes the layer momentum advection flux terms using upwind flux
        ! This routine computes the layer momentum pressure terms

        use mod_constants, only : gravity
        use mod_initial, only : alpha_mlswe, zbot_face, pbprime_df
        use mod_grid, only : nface, npoin, npoin_q, face, intma, face_type
        use mod_basis, only : nq, psiq, ngl
        use mod_input, only : nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: ope_face_ave, H_face_ave, one_plus_eta_edge_2_ave, &
                                uvb_face_ave, Quv_face_ave, Qu_face_ave, Qv_face_ave, ope2_face_ave
        use mod_variables, only: sum_layer_mass_flux_face, btp_mass_flux_face_ave

        implicit none

        real, dimension(3, npoin, nlayers), intent(inout) :: rhs
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

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
        real :: uu_dp_flux_deficit(2), vv_dp_flux_deficit(2)
        real, parameter :: eps1 = 1.0e-20, eps = 1.0e-10 !  Parameter used to prevent division by zero.
        integer :: il, jl, ir, jr, kl, kr, jquad, n, m
        real :: wq, hi, hl_k, hr_k, flux, flux_u, flux_v, un
        real, dimension(2,nq,nlayers) :: H_face, dp_flux, udp_flux, vdp_flux
        real :: dp_deficit_l(2), dp_deficit_r(2), dp_deficit(2,nlayers)
        real :: pbprime_l, pbprime_r, wght_l(nlayers), wght_r(nlayers)

        do k=1,nlayers
            alpha_over_g(k) = alpha_mlswe(k)/gravity
            g_over_alpha(k) = gravity/alpha_mlswe(k)
        enddo

        ! Compute H_r at the element face
        do iface = 1, nface

            if (face_type(iface) == 2) cycle

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

                    pbprime_l = 0.0 ; pbprime_r = 0.0
                    do n = 1, ngl

                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)
                        I = intma(il,jl,kl,el)
                        hi = psiq(n,iquad)

                        ql(:,iquad,k) = ql(:,iquad,k) + hi*qprime_df(:,I,k)
                        pbprime_l = pbprime_l + hi*pbprime_df(I)
                    enddo

                    if (er > 0) then
                        do n = 1, ngl

                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)
                            hi = psiq(n,iquad)

                            qr(:,iquad,k) = qr(:,iquad,k) + hi*qprime_df(:,I,k)
                            pbprime_r = pbprime_r + hi*pbprime_df(I)
                        enddo
                    else 
                        qr(:,iquad,k) = ql(:,iquad,k)
                        pbprime_r = pbprime_l
                    endif

                    ! Apply wall boundary conditions
                    if (er == -4) then
                        un = ql(2,iquad,k)*nxl + ql(3,iquad,k)*nyl
                        qr(2,iquad,k) = ql(2,iquad,k) - 2.0*un*nxl
                        qr(3,iquad,k) = ql(3,iquad,k) - 2.0*un*nyl
                    elseif (er == -2) then
                        qr(2,iquad,k) = -ql(2,iquad,k)
                        qr(3,iquad,k) = -ql(3,iquad,k)
                    endif

                    ! Left side of the edge
                    dpl = qbl(1,iquad) * ql(1,iquad,k)
                    dpr = qbr(1,iquad) * qr(1,iquad,k)
                    ul = ql(2,iquad,k) + qbl(2,iquad)
                    ur = qr(2,iquad,k) + qbr(2,iquad)
                    vl = ql(3,iquad,k) + qbl(3,iquad)
                    vr = qr(3,iquad,k) + qbr(3,iquad)

                    uu = 0.5*(ul + ur)
                    vv = 0.5*(vl + vr)
                    udpl(iquad,k) = ul*dpl
                    udpr(iquad,k) = ur*dpr
                    vdpl(iquad,k) = vl*dpl
                    vdpr(iquad,k) = vr*dpr

                    if(uu*nxl > 0.0) then
                        dp_flux(1,iquad,k) = uu*dpl
                        udp_flux(1,iquad,k) = uu * (ul*dpl)
                        vdp_flux(1,iquad,k) = uu * (vl*dpl)
                    else
                        dp_flux(1,iquad,k) = uu*dpr
                        udp_flux(1,iquad,k) = uu * (ur*dpr)
                        vdp_flux(1,iquad,k) = uu * (vr*dpr)
                    endif
                    if(vv*nyl > 0.0) then
                        dp_flux(2,iquad,k)  = vv * dpl
                        udp_flux(2,iquad,k) = vv * (ul*dpl)
                        vdp_flux(2,iquad,k) = vv * (vl*dpl)
                    else
                        dp_flux(2,iquad,k)  = vv * dpr
                        udp_flux(2,iquad,k) = vv * (ur*dpr)
                        vdp_flux(2,iquad,k) = vv * (vr*dpr)
                    endif

                    wght_l(k) = ql(1,iquad,k) / pbprime_l
                    wght_r(k) = qr(1,iquad,k) / pbprime_r

                enddo

                ! Decficit for mass equation
                dp_deficit_l(1) = btp_mass_flux_face_ave(1,iquad,iface) &
                                                        - sum(dp_flux(1,iquad,:))
                dp_deficit_l(2) = btp_mass_flux_face_ave(2,iquad,iface) &
                                                    - sum(dp_flux(2,iquad,:))

                dp_deficit_r(1) = btp_mass_flux_face_ave(1,iquad,iface) &
                                                    - sum(dp_flux(1,iquad,:))
                dp_deficit_r(2) = btp_mass_flux_face_ave(2,iquad,iface) &
                                                        - sum(dp_flux(2,iquad,:))

                ! Adjustments for consistency
                uu_dp_flux_deficit(1) = Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,iquad,:))
                uu_dp_flux_deficit(2) = Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,iquad,:))
                vv_dp_flux_deficit(1) = Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,iquad,:))
                vv_dp_flux_deficit(2) = Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,iquad,:))

                do k = 1,nlayers

                    ! Adjust the fluxes for the mass equation
                    weight = wght_l(k)
                    if(dp_deficit_l(1)*nxl <= 0.0) weight = wght_r(k)
                    dp_flux(1,iquad,k) = dp_flux(1,iquad,k) + weight* dp_deficit_l(1)

                    weight = wght_l(k)
                    if(dp_deficit_l(2)*nyl <= 0.0) weight = wght_r(k)
                    dp_flux(2,iquad,k) = dp_flux(2,iquad,k) + weight*dp_deficit_l(2)

                    ! Adjust the fluxes for the u-momentum equation
                    !x-direction
                    weight = abs(udpl(iquad,k)) / sum(abs(udpl(iquad,:))+eps1)
                    if(uu_dp_flux_deficit(1)*nxl < 0.0) &
                        weight = abs(udpr(iquad,k)) / sum(abs(udpr(iquad,:))+eps1)

                    udp_flux(1,iquad,k) = udp_flux(1,iquad,k) + weight * uu_dp_flux_deficit(1)

                    !y-direction
                    weight = abs(udpl(iquad,k)) / sum(abs(udpl(iquad,:))+eps1)
                    if(uu_dp_flux_deficit(2)*nyl < 0.0) &
                        weight = abs(udpr(iquad,k)) / sum(abs(udpr(iquad,:))+eps1)

                    udp_flux(2,iquad,k) = udp_flux(2,iquad,k) + weight * uu_dp_flux_deficit(2)
                
                    ! Adjust the fluxes for the v-momentum equation
                    !x-direction
                    weight = abs(vdpl(iquad,k)) / sum(abs(vdpl(iquad,:))+eps1)
                    if(vv_dp_flux_deficit(1)*nxl < 0.0) &
                        weight = abs(vdpr(iquad,k)) / sum(abs(vdpr(iquad,:))+eps1)

                    vdp_flux(1,iquad,k) = vdp_flux(1,iquad,k) + weight * vv_dp_flux_deficit(1)
  
                    !y-direction
                    weight = abs(vdpl(iquad,k)) / sum(abs(vdpl(iquad,:))+eps1)
                    if(vv_dp_flux_deficit(2)*nyl < 0.0) &
                        weight = abs(vdpr(iquad,k)) / sum(abs(vdpr(iquad,:))+eps1)

                    vdp_flux(2,iquad,k) = vdp_flux(2,iquad,k) + weight * vv_dp_flux_deficit(2)
                end do

                z_face = 0.0 ; p_face = 0.0
                z_edge_plus = 0.0 ; z_edge_minus = 0.0
                p_edge_plus = 0.0; p_edge_minus = 0.0

                !Store Left Side Variables
                ope_l = sqrt(ope2_face_ave(1,iquad,iface))
                ope_r = sqrt(ope2_face_ave(2,iquad,iface))
                p_face(1,1) = 0.0
                p_face(2,1) = 0.0
                do k=1,nlayers
                    p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,iquad,k)
                    p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,iquad,k)
                end do

                one_plus_eta_edge = sqrt(one_plus_eta_edge_2_ave(iquad,iface))
                z_face(1,nlayers+1) = zbot_face(1,iquad,iface)
                z_face(2,nlayers+1) = zbot_face(2,iquad,iface)
                z_edge_plus(nlayers+1) = zbot_face(1,iquad,iface)
                z_edge_minus(nlayers+1) = zbot_face(2,iquad,iface)
                do k=nlayers,1,-1
                    z_face(1,k) = z_face(1,k+1) + alpha_over_g(k) * (ope_l * ql(1,iquad,k))
                    z_face(2,k) = z_face(2,k+1) + alpha_over_g(k) * (ope_r * qr(1,iquad,k))
                    z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k) * &
                                                        (one_plus_eta_edge * ql(1,iquad,k))
                    z_edge_minus(k) = z_edge_minus(k+1) + alpha_over_g(k) * &
                                                        (one_plus_eta_edge * qr(1,iquad,k))
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
                            p_intersect_bot = p_edge_minus(ktemp+1) &
                                    - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_minus(ktemp+1))
                            p_intersect_top = p_edge_minus(ktemp+1) &
                                    - g_over_alpha(ktemp)*(z_intersect_top - z_edge_minus(ktemp+1))
                            H_r_minus = H_r_minus + &
                                    0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)

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
                            p_intersect_bot = p_edge_plus(ktemp+1) &
                                - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_plus(ktemp+1))
                            p_intersect_top = p_edge_plus(ktemp+1) &
                                - g_over_alpha(ktemp)*(z_intersect_top - z_edge_plus(ktemp+1))
                            H_r_plus = H_r_plus &
                                + 0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
                        end if
                    end do
                    H_face(2,iquad,k) = 0.5*(H_r_plus + H_r_minus) ! H_r for the right side
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
                        H_corr1 = 0.5 * alpha_mlswe(k) * ((p_face(1,k+1) &
                                    + p_inc1)**2 - p_face(1,k+1)**2)
                        H_face(1,iquad,k) = H_face(1,iquad,k) - H_corr1
                        H_face(1,iquad,k+1) = H_face(1,iquad,k+1) + H_corr1

                        ! Corrections at the right side of a face.
                        p_inc2 = g_over_alpha(k)*(z_face(2,k+1) - z_edge_minus(k+1))
                        H_corr2 = 0.5 * alpha_mlswe(k) * ((p_face(2,k+1) &
                                    + p_inc2)**2 - p_face(2,k+1)**2)
                        H_face(2,iquad,k) = H_face(2,iquad,k) - H_corr2
                        H_face(2,iquad,k+1) = H_face(2,iquad,k+1) + H_corr2

                    end do
                end if

                ! Adjust the values of  H_k, at element faces, so that the vertical sum of
                ! H_k  over all layers equals the time average of the barotropic forcing  H 
                ! over all barotropic substeps of the baroclinic
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

                    hl_k = H_face(1,iquad,k)
                    hr_k = H_face(2,iquad,k)

                    flux = nxl*dp_flux(1,iquad,k) + nyl*dp_flux(2,iquad,k)
                    flux_u = nxl*udp_flux(1,iquad,k) + nyl*udp_flux(2,iquad,k)
                    flux_v = nxl*vdp_flux(1,iquad,k) + nyl*vdp_flux(2,iquad,k)

                    do n = 1, ngl

                        hi = psiq(n,iquad)
                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)
                        I = intma(il,jl,kl,el)

                        rhs(1,I,k) = rhs(1,I,k) - wq*hi*flux
                        rhs(2,I,k) = rhs(2,I,k) - wq*hi*(nxl*hl_k + flux_u)
                        rhs(3,I,k) = rhs(3,I,k) - wq*hi*(nyl*hl_k + flux_v)

                        if(er > 0) then

                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)

                            rhs(1,I,k) = rhs(1,I,k) + wq*hi*flux
                            rhs(2,I,k) = rhs(2,I,k) + wq*hi*(nxl*hr_k + flux_u)
                            rhs(3,I,k) = rhs(3,I,k) + wq*hi*(nyl*hr_k + flux_v)

                        end if
                    end do
                end do
            end do
        end do ! iface

    end subroutine Apply_bcl_fluxes

    subroutine create_layers_volume_mass(dp_advec, qprime_df)

        use mod_grid, only : npoin_q, npoin, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, alpha_mlswe, pbprime_df
        use mod_variables, only: uvb_ave, ope_ave, sum_layer_mass_flux, btp_mass_flux_ave

        implicit none

        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(npoin,nlayers), intent(out) :: dp_advec

        real :: wq, hi, dhde, dhdn
        integer :: k, I, Iq, ip
        real, dimension(nlayers) :: dp, udp, vdp
        real, dimension(3) :: qp, qb
        real, parameter :: eps = 1.0e-10
        real :: pbq, weight_dp, flux(2)
        real, dimension(2) :: sum_layer_mass

        dp_advec = 0.0
        sum_layer_mass_flux = 0.0

        do Iq = 1, npoin_q

            qb(1) = ope_ave(Iq)
            qb(2) = uvb_ave(1,Iq)
            qb(3) = uvb_ave(2,Iq)
            wq = wjac(Iq)

            sum_layer_mass = 0.0
            do k=1,nlayers

                qp = 0.0 ; pbq = 0.0
                do ip = 1,npts
                    I = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    qp(:) = qp(:) + hi*qprime_df(:,I,k)
                    pbq = pbq + hi*pbprime_df(I)
                enddo

                if (qp(1) <= (gravity/alpha_mlswe(k))*eps) then
                    qp(1) = (gravity/alpha_mlswe(k))*eps
                    qp(2:3) = 0.0
                end if

                dp(k) = qp(1) * qb(1)
                udp(k) = (qp(2)+qb(2)) * dp(k)
                vdp(k) = (qp(3)+qb(3)) * dp(k)

                sum_layer_mass(1) = sum_layer_mass(1) + udp(k)
                sum_layer_mass(2) = sum_layer_mass(2) + vdp(k)
            end do !k

            do k = 1,nlayers

                weight_dp = dp(k) / pbq
                flux(1) = udp(k) + weight_dp * (btp_mass_flux_ave(1,Iq) - sum_layer_mass(1))
                flux(2) = vdp(k) + weight_dp * (btp_mass_flux_ave(2,Iq) - sum_layer_mass(2))

                do ip = 1, npts

                    I = indexq(ip,Iq)
                    dp_advec(I,k) = dp_advec(I,k) + wq*(dpsidx(ip,Iq)*flux(1) + dpsidy(ip,Iq)*flux(2))

                end do
            end do

        end do !k

    end subroutine create_layers_volume_mass

    subroutine create_layer_mass_flux(rhs, qprime_df)

        ! This routine computes the layer momentum advection flux terms using upwind flux
        ! This routine computes the layer momentum pressure terms

        use mod_basis, only: nq, psiq, ngl
        use mod_grid, only:  npoin_q, intma,  nface, face, face_type
        use mod_input, only: nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: sum_layer_mass_flux_face, ope_face_ave, uvb_face_ave, btp_mass_flux_face_ave
        use mod_constants, only : gravity
        use mod_initial, only : alpha_mlswe, pbprime_df

        implicit none

        real, dimension(npoin, nlayers), intent(inout) :: rhs
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        real, dimension(3,nq,nlayers) :: ql, qr
        real, dimension(3,nq) :: qbl, qbr
        integer :: iface, ilr, k, iquad, I, Iq, el, er
        integer :: il, jl, ir, jr, kl, kr, jquad, n, m
        real ::  ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv
        real :: wq, hi, flux, un
        real, dimension(2,nq,nlayers) :: dp_flux
        real :: dp_deficit_l(2), dp_deficit_r(2), dp_deficit(2,nlayers)
        real :: pbprime_l, pbprime_r, wght_l(nlayers), wght_r(nlayers)

        ! Compute H_r at the element face
        do iface = 1, nface

            if (face_type(iface) == 2) cycle

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

                    pbprime_l = 0.0 ; pbprime_r = 0.0
                    do n = 1, ngl

                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)
                        I = intma(il,jl,kl,el)
                        hi = psiq(n,iquad)

                        ql(:,iquad,k) = ql(:,iquad,k) + hi*qprime_df(:,I,k)
                        pbprime_l = pbprime_l + hi*pbprime_df(I)
                    enddo

                    if (er > 0) then
                        do n = 1, ngl

                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)
                            hi = psiq(n,iquad)

                            qr(:,iquad,k) = qr(:,iquad,k) + hi*qprime_df(:,I,k)
                            pbprime_r = pbprime_r + hi*pbprime_df(I)
                        enddo
                    else 
                        qr(:,iquad,k) = ql(:,iquad,k)
                        pbprime_r = pbprime_l
                    endif

                    ! Apply wall boundary conditions
                    if (er == -4) then
                        un = ql(2,iquad,k)*nxl + ql(3,iquad,k)*nyl
                        qr(2,iquad,k) = ql(2,iquad,k) - 2.0*un*nxl
                        qr(3,iquad,k) = ql(3,iquad,k) - 2.0*un*nyl
                    elseif (er == -2) then
                        qr(2,iquad,k) = -ql(2,iquad,k)
                        qr(3,iquad,k) = -ql(3,iquad,k)
                    endif

                    ! Left side of the edge
                    dpl = qbl(1,iquad) * ql(1,iquad,k)
                    dpr = qbr(1,iquad) * qr(1,iquad,k)
                    ul = ql(2,iquad,k) + qbl(2,iquad)
                    ur = qr(2,iquad,k) + qbr(2,iquad)
                    vl = ql(3,iquad,k) + qbl(3,iquad)
                    vr = qr(3,iquad,k) + qbr(3,iquad)

                    uu = 0.5*(ul+ur) ; vv = 0.5*(vl+vr)
                    dp_flux(1,iquad,k) = uu*dpl
                    if(uu*nxl < 0.0) dp_flux(1,iquad,k) = uu*dpr
                    
                    dp_flux(2,iquad,k)  = vv * dpl
                    if(vv*nyl < 0.0) dp_flux(2,iquad,k)  = vv * dpr

                    wght_l(k) = dpl / pbprime_l
                    wght_r(k) = dpr / pbprime_r

                enddo

                ! Decficit for mass equation

                dp_deficit_l(1) = btp_mass_flux_face_ave(1,iquad,iface) &
                                                        - sum(dp_flux(1,iquad,:))
                dp_deficit_l(2) = btp_mass_flux_face_ave(2,iquad,iface) &
                                                    - sum(dp_flux(2,iquad,:))

                dp_deficit_r(1) = btp_mass_flux_face_ave(1,iquad,iface) &
                                                    - sum(dp_flux(1,iquad,:))
                dp_deficit_r(2) = btp_mass_flux_face_ave(2,iquad,iface) &
                                                    - sum(dp_flux(2,iquad,:))

                ! Adjustments for consistency
                do k = 1,nlayers
                    ! Adjust the fluxes for the mass equation
                    dp_flux(1,iquad,k) = dp_flux(1,iquad,k) + wght_l(k)* dp_deficit_l(1)
                    if(dp_deficit_l(1)*nxl < 0.0) &
                        dp_flux(1,iquad,k) =  dp_flux(1,iquad,k) + wght_r(k)* dp_deficit_r(1)

                    dp_flux(2,iquad,k) = dp_flux(2,iquad,k) + wght_l(k)* dp_deficit_l(2)
                    if(dp_deficit_l(2)*nyl < 0.0) &
                        dp_flux(2,iquad,k) =  dp_flux(2,iquad,k) + wght_r(k)* dp_deficit_r(2)
                end do

            end do ! iquad
            
            ! Do Gauss-Lobatto Integration
            do k = 1,nlayers

                do iquad = 1, nq

                    wq = jac_faceq(iquad,1,iface)
                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    flux = nxl*dp_flux(1,iquad,k) + nyl*dp_flux(2,iquad,k)

                    do n = 1, ngl

                        hi = psiq(n,iquad)
                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)
                        I = intma(il,jl,kl,el)

                        rhs(I,k) = rhs(I,k) - wq*hi*flux

                        if(er > 0) then

                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)

                            rhs(I,k) = rhs(I,k) + wq*hi*flux
                        
                        end if
                    end do
                end do
            end do
        end do ! iface

    end subroutine create_layer_mass_flux

end module mod_create_rhs_mlswe
