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
        use mod_input, only: nlayers, dry_cutoff
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, alpha_mlswe, &
                                tau_wind, pbprime_df, zbot_df, grad_zbot_quad
        use mod_variables, only: tau_bot_ave, ope_ave, uvb_ave, H_ave, &
                                    Qu_ave, Qv_ave, Quv_ave, uvb_ave, ope2_ave_df, ope2_ave

        implicit none

        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom
        real, dimension(3,npoin,nlayers), intent(in) :: q_df, qprime_df

        real :: wq, hi, dhdx, dhdy, bot_layer, tau_wind_u, tau_wind_v, temp1
        real :: Hq, var_uu, var_uv, var_vu, var_vv, source_x, source_y, Pstress
        integer :: k, I, Iq, ip
        real, dimension(nlayers+1) :: pprime_temp, z
        real :: temp_dp, temp_u, temp_v, Ptop_k, Pbot_k, tempbot, Pbstress
        real :: weight, acceleration, pbq
        real, dimension(3) :: qp, qb
        real, dimension(nlayers) :: temp_uu, temp_vv, H_tmp, u_udp, v_vdp
        real, dimension(2,nlayers) :: u_vdp
        real :: p_tmp(nlayers+1), u, v, dp, weightq, one_over_sumuq, one_over_sumvq
        real :: uu_dp_deficitq, uv_dp_deficitq, vv_dp_deficitq, gradz(2,nlayers+1)
        real :: z_elv(npoin,nlayers+1)
        real, parameter :: eps1 = 1.0e-10 !  Parameter used to prevent division by zero.
        real :: flux(2,2)

        rhs_mom = 0.0
        bot_layer = 0.0
        pprime_temp = 0.0

        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth 
                                                  ! at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth 
                                                         ! at which bottom stress is reduced to 0

        ! Find layer interfaces
        ! z_elv(:,nlayers+1) = zbot_df(:)
        ! do k = nlayers,1,-1
        !     z_elv(:,k) = z_elv(:,k+1) + (alpha_mlswe(k)/gravity) * &
        !                                 (sqrt(ope2_ave_df(:))*qprime_df(1,:,k))
        ! end do

        do Iq = 1, npoin_q

            ! init
            p_tmp(1) = 0.0

            qb(1) = ope_ave(Iq)
            qb(2) = uvb_ave(1,Iq)
            qb(3) = uvb_ave(2,Iq)

            !---------------------------------------
            ! Build layer quantities
            !---------------------------------------
            do k = 1, nlayers

                qp = 0.0
                do ip = 1, npts
                    I  = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    qp(:) = qp(:) + hi*qprime_df(:,I,k)
                end do

                dp = qp(1) * qb(1)
                u  = qp(2) + qb(2)
                v  = qp(3) + qb(3)

                if (qp(1) <= (gravity/alpha_mlswe(k)) * dry_cutoff) then
                    qp(1) = (gravity/alpha_mlswe(k)) * dry_cutoff
                    dp = qp(1) * qb(1)
                    u = 0.0
                    v = 0.0
                end if

                u_udp(k)    = u*(u*dp)
                v_vdp(k)    = v*(v*dp)
                u_vdp(1,k)  = v*(u*dp)
                u_vdp(2,k)  = u*(v*dp)

                p_tmp(k+1) = p_tmp(k) + sqrt(ope2_ave(Iq)) * qp(1)
                H_tmp(k)   = 0.5*alpha_mlswe(k) * (p_tmp(k+1)**2 - p_tmp(k)**2)

                temp_uu(k)  = abs(dp*u) + eps1
                temp_vv(k)  = abs(dp*v) + eps1

            end do

            !---------------------------------------
            ! gradz and pbq
            !---------------------------------------

            gradz(:,:) = 0.0 ; pbq = 0.0
            do ip = 1, npts
                I = indexq(ip,Iq)
                ! z(:) = 0.0
                z(nlayers+1) = zbot_df(I)
                do k = nlayers,1,-1
                    z(k) = z(k+1) + (alpha_mlswe(k)/gravity)*(sqrt(ope2_ave_df(I))*qprime_df(1,I,k))

                    gradz(1,k) = gradz(1,k) + dpsidx(ip,Iq) * z(k)
                    gradz(2,k) = gradz(2,k) + dpsidy(ip,Iq) * z(k)

                    ! if (abs(z(k) - z(k+1)) <= dry_cutoff) then
                    !     gradz(1,k) = 0.0
                    !     gradz(2,k) = 0.0
                    ! end if
                enddo 

                gradz(1,nlayers+1) = grad_zbot_quad(1,Iq)
                gradz(2,nlayers+1) = grad_zbot_quad(2,Iq)

                pbq = pbq + psih(ip,Iq)  * pbprime_df(I)
            end do

            !---------------------------------------
            ! Consistency deficits
            !---------------------------------------
            uu_dp_deficitq = Qu_ave(Iq)  - sum(u_udp(:))
            uv_dp_deficitq = Quv_ave(Iq) - sum(u_vdp(1,:))
            vv_dp_deficitq = Qv_ave(Iq)  - sum(v_vdp(:))

            one_over_sumuq = 1.0 / sum(temp_uu(:))
            one_over_sumvq = 1.0 / sum(temp_vv(:))

            wq = wjac(Iq)

            acceleration = sum(H_tmp(:))
            weight = 1.0
            if (acceleration > 0.0) weight = H_ave(Iq) / acceleration

            !---------------------------------------
            ! Main layer loop + RHS update
            !---------------------------------------
            do k = 1, nlayers

                weightq   = temp_uu(k) * one_over_sumuq
                flux(1,1) = u_udp(k)    + weightq * uu_dp_deficitq
                flux(2,1) = u_vdp(1,k)  + weightq * uv_dp_deficitq

                weightq   = temp_vv(k) * one_over_sumvq
                flux(1,2) = u_vdp(2,k)  + weightq * uv_dp_deficitq
                flux(2,2) = v_vdp(k)    + weightq * vv_dp_deficitq

                Hq = H_tmp(k) * weight
                flux(1,1) = flux(1,1) + Hq
                flux(2,2) = flux(2,2) + Hq

                temp1      = (min(p_tmp(k+1), Pstress) - min(p_tmp(k), Pstress)) / Pstress
                tau_wind_u = temp1 * tau_wind(1,Iq)
                tau_wind_v = temp1 * tau_wind(2,Iq)

                tempbot = min(Pbstress, pbq - p_tmp(k+1)) - min(Pbstress, pbq - p_tmp(k))
                tempbot = tempbot / Pbstress

                source_x = gravity*(tau_wind_u - tempbot*tau_bot_ave(1,Iq) + p_tmp(k)*gradz(1,k) - p_tmp(k+1)*gradz(1,k+1))
                source_y = gravity*(tau_wind_v - tempbot*tau_bot_ave(2,Iq) + p_tmp(k)*gradz(2,k) - p_tmp(k+1)*gradz(2,k+1))

                do ip = 1, npts
                    I  = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    dhdx = dpsidx(ip,Iq)
                    dhdy = dpsidy(ip,Iq)

                    rhs_mom(1,I,k) = rhs_mom(1,I,k) + wq*(hi*source_x + dhdx*flux(1,1) + dhdy*flux(1,2))
                    rhs_mom(2,I,k) = rhs_mom(2,I,k) + wq*(hi*source_y + dhdx*flux(2,1) + dhdy*flux(2,2))
                end do

            end do

        end do

    end subroutine create_rhs_dynamics_volume_layers

    subroutine create_rhs_dynamics_volume_bcl(rhs, qprime_df, q_df)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts, dpsiqx
        use mod_input, only: nlayers, dry_cutoff
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, alpha_mlswe, &
                                tau_wind, pbprime_df, zbot_df, grad_zbot_quad
        use mod_variables, only: tau_bot_ave, ope_ave, uvb_ave, H_ave, &
                                    Qu_ave, Qv_ave, Quv_ave, uvb_ave, &
                                    ope2_ave_df, ope2_ave, sum_layer_mass_flux, btp_mass_flux_ave

        implicit none

        real, dimension(3,npoin,nlayers), intent(out) :: rhs
        real, dimension(3,npoin,nlayers), intent(in) :: q_df, qprime_df

        real :: wq, hi, dhdx, dhdy, bot_layer, tau_wind_u, tau_wind_v, temp1
        real :: Hq, var_uu, var_uv, var_vu, var_vv, source_x, source_y, Pstress
        integer :: k, I, Iq, ip
        real, dimension(nlayers+1) :: pprime_temp, z
        real :: temp_dp, temp_u, temp_v, Ptop_k, Pbot_k, tempbot, Pbstress
        real :: weight, acceleration, pbq
        real, dimension(3) :: qp, qb
        real, dimension(nlayers) :: temp_uu, temp_vv, H_tmp, u_udp, v_vdp, udp, vdp,  dp
        real, dimension(2,nlayers) :: u_vdp
        real :: p_tmp(nlayers+1), u, v,weightq, one_over_sumuq, one_over_sumvq
        real :: uu_dp_deficitq, uv_dp_deficitq, vv_dp_deficitq, gradz(2,nlayers+1)
        real :: z_elv(npoin,nlayers+1)
        real, parameter :: eps1 = 1.0e-10 !  Parameter used to prevent division by zero.
        real :: flux(3,3), cst_dp(nlayers)

        rhs = 0.0
        bot_layer = 0.0
        pprime_temp = 0.0
        sum_layer_mass_flux = 0.0

        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth 
                                                  ! at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth 
                                                         ! at which bottom stress is reduced to 0

        do Iq = 1, npoin_q

            qb(1) = ope_ave(Iq)
            qb(2) = uvb_ave(1,Iq)
            qb(3) = uvb_ave(2,Iq)

            p_tmp(1) = 0.0
            cst_dp(:) = 1.0

            ! ---- build layer quantities ----
            do k = 1, nlayers

                qp(:) = 0.0
                do ip = 1, npts
                    I  = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    qp(1) = qp(1) + hi * qprime_df(1, I, k)   ! (1)=p'
                    qp(2) = qp(2) + hi * qprime_df(2, I, k)   ! (2)=u'
                    qp(3) = qp(3) + hi * qprime_df(3, I, k)   ! (3)=v'
                end do

                dp(k) = qp(1) * qb(1)
                u     = qp(2) + qb(2)
                v     = qp(3) + qb(3)

                if (qp(1) < (gravity/alpha_mlswe(k)) * dry_cutoff) then
                    qp(1) = (gravity/alpha_mlswe(k)) * dry_cutoff
                    dp(k) = qp(1) * qb(1)
                    u = 0.0
                    v = 0.0
                    cst_dp(k) = 0.0
                end if

                udp(k) = u * dp(k)
                vdp(k) = v * dp(k)

                u_udp(k) = u * udp(k)
                v_vdp(k) = v * vdp(k)

                u_vdp(1,k) = v * udp(k)
                u_vdp(2,k) = u * vdp(k)

                p_tmp(k+1) = p_tmp(k) + sqrt(ope2_ave(Iq)) * qp(1)
                H_tmp(k)   = 0.5*alpha_mlswe(k) * (p_tmp(k+1)**2 - p_tmp(k)**2)

                temp_uu(k) = abs(udp(k)) !+ eps1
                temp_vv(k) = abs(vdp(k)) !+ eps1

            end do

            ! ---- gradz and pbq ----
            gradz(:,:) = 0.0 ; pbq = 0.0
            do ip = 1, npts
                I = indexq(ip,Iq)
                z(nlayers+1) = zbot_df(I)
                do k = nlayers,1,-1
                    z(k) = z(k+1) + (alpha_mlswe(k)/gravity)*(sqrt(ope2_ave_df(I))*qprime_df(1,I,k))

                    gradz(1,k) = gradz(1,k) + dpsidx(ip,Iq) * z(k)
                    gradz(2,k) = gradz(2,k) + dpsidy(ip,Iq) * z(k)
                end do

                gradz(1,nlayers+1) = grad_zbot_quad(1,Iq)
                gradz(2,nlayers+1) = grad_zbot_quad(2,Iq)

                pbq = pbq + psih(ip,Iq)  * pbprime_df(I)
            end do

            ! ---- deficits (consistency) ----
            uu_dp_deficitq = Qu_ave(Iq)  - sum(u_udp(:))
            uv_dp_deficitq = Quv_ave(Iq) - sum(u_vdp(1,:))
            vv_dp_deficitq = Qv_ave(Iq)  - sum(v_vdp(:))

            one_over_sumuq = 1.0 / sum(temp_uu(:) + eps1)
            one_over_sumvq = 1.0 / sum(temp_vv(:) + eps1)

            wq = wjac(Iq)

            do k = 1, nlayers

                

                weightq  = temp_uu(k) * one_over_sumuq
                
                u_udp(k)  = u_udp(k)     + weightq * uu_dp_deficitq
                u_vdp(1,k)= u_vdp(1,k)   + weightq * uv_dp_deficitq

                weightq  = temp_vv(k) * one_over_sumvq
                u_vdp(2,k)= u_vdp(2,k)   + weightq * uv_dp_deficitq
                v_vdp(k)  = v_vdp(k)     + weightq * vv_dp_deficitq

                Hq = H_tmp(k)

                acceleration = sum(H_tmp(:))
                weight = 1.0
                if (acceleration > 0.0) weight = H_ave(Iq) / acceleration
                Hq = Hq * weight

                flux(1,1) = udp(k) + cst_dp(k) * (dp(k)/sum(dp(:))) * (btp_mass_flux_ave(1,Iq) - sum(udp(:)))
                flux(2,1) = vdp(k) + cst_dp(k) * (dp(k)/sum(dp(:))) * (btp_mass_flux_ave(2,Iq) - sum(vdp(:)))

                flux(1,2) = u_udp(k) + Hq
                flux(2,2) = u_vdp(1,k)
                flux(1,3) = u_vdp(2,k)
                flux(2,3) = v_vdp(k) + Hq

                temp1      = (min(p_tmp(k+1), Pstress) - min(p_tmp(k), Pstress)) / Pstress
                tau_wind_u = temp1 * tau_wind(1,Iq)
                tau_wind_v = temp1 * tau_wind(2,Iq)

                tempbot = min(Pbstress, pbq - p_tmp(k+1)) - min(Pbstress, pbq - p_tmp(k))
                tempbot = tempbot / Pbstress

                source_x = gravity*( tau_wind_u - tempbot*tau_bot_ave(1,Iq) + p_tmp(k)*gradz(1,k) - p_tmp(k+1)*gradz(1,k+1) )
                source_y = gravity*( tau_wind_v - tempbot*tau_bot_ave(2,Iq) + p_tmp(k)*gradz(2,k) - p_tmp(k+1)*gradz(2,k+1) )

                do ip = 1, npts
                    I  = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    dhdx = dpsidx(ip,Iq)
                    dhdy = dpsidy(ip,Iq)

                    rhs(1,I,k) = rhs(1,I,k) + wq*(dhdx*flux(1,1) + dhdy*flux(2,1))
                    rhs(2,I,k) = rhs(2,I,k) + wq*(hi*source_x + dhdx*flux(1,2) + dhdy*flux(2,2))
                    rhs(3,I,k) = rhs(3,I,k) + wq*(hi*source_y + dhdx*flux(1,3) + dhdy*flux(2,3))
                end do

            end do

        end do !Iq

    end subroutine create_rhs_dynamics_volume_bcl

    subroutine Apply_layers_fluxes(rhs_mom, qprime_df)

        ! This routine computes the layer momentum advection flux terms using upwind flux
        ! This routine computes the layer momentum pressure terms

        use mod_constants, only : gravity
        use mod_initial, only : alpha_mlswe, zbot_face
        use mod_grid, only : nface, npoin, npoin_q, face, intma, face_type
        use mod_basis, only : nq, psiq, ngl
        use mod_input, only : nlayers, dry_cutoff
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: ope_face_ave, H_face_ave, one_plus_eta_edge_2_ave, &
                                uvb_face_ave, Quv_face_ave, Qu_face_ave, Qv_face_ave, ope2_face_ave

        implicit none

        real, dimension(2, npoin, nlayers), intent(inout) :: rhs_mom
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        real, dimension(nlayers) :: alpha_over_g, g_over_alpha
        real, dimension(2,nlayers+1) :: p_face, z_face
        real, dimension(nlayers+1) :: p_edge_plus, p_edge_minus, p2l, p2r, z_edge_plus, z_edge_minus
        real, dimension(3,nlayers) :: ql, qr
        real, dimension(3,nq) :: qbl, qbr
        integer :: iface, ilr, k, iquad, ktemp, I
        real :: z_intersect_top,z_intersect_bot, dz_intersect, H_r_plus, H_r_minus, acceleration
        real :: p_intersect_bot, p_intersect_top, one_plus_eta_edge
        real :: H_corr,p_inc, weight, H_corr1,p_inc1, H_corr2,p_inc2, temp, ope_l, ope_r
        integer :: Iq, el, er
        real ::  ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv
        real, dimension(nlayers) :: udpl, udpr, vdpl, vdpr
        real :: uu_dp_flux_deficit(2), vv_dp_flux_deficit(2)
        real, parameter :: eps1 = 1.0e-20 !  Parameter used to prevent division by zero.
        integer :: il, jl, ir, jr, kl, kr, jquad, n, m
        real :: wq, hi, hlx_k, hly_k, hrx_k, hry_k, flux_x, flux_y, hx_k, hy_k, un
        real, dimension(2,nlayers,nq) :: H_face, udp_flux, vdp_flux
        real, dimension(nlayers,nq) :: flux_ul, flux_ur, flux_vl, flux_vr
        real :: flux_xl, flux_xr, flux_yl, flux_yr

        do k=1,nlayers
            alpha_over_g(k) = alpha_mlswe(k)/gravity
            g_over_alpha(k) = gravity/alpha_mlswe(k)
        enddo

        ! Compute H_r at the element face
        do concurrent(iface = 1:nface, iquad = 1:nq)

            if (face_type(iface) == 2) cycle

            ! Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            qbl(1,iquad) = ope_face_ave(1,iquad,iface)
            qbl(2,iquad) = uvb_face_ave(1,1,iquad,iface)
            qbl(3,iquad) = uvb_face_ave(2,1,iquad,iface)
            qbr(1,iquad) = ope_face_ave(2,iquad,iface)
            qbr(2,iquad) = uvb_face_ave(1,2,iquad,iface)
            qbr(3,iquad) = uvb_face_ave(2,2,iquad,iface)

            nxl = normal_vector_q(1,iquad,1,iface)
            nyl = normal_vector_q(2,iquad,1,iface)

            ql = 0.0; qr = 0.0

            do k = 1,nlayers

                do n = 1, ngl

                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)
                    I = intma(il,jl,kl,el)
                    hi = psiq(n,iquad)

                    ql(1,k) = ql(1,k) + hi*qprime_df(1,I,k)
                    ql(2,k) = ql(2,k) + hi*qprime_df(2,I,k)
                    ql(3,k) = ql(3,k) + hi*qprime_df(3,I,k)
                enddo

                if (er > 0) then
                    do n = 1, ngl

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)
                        I = intma(ir,jr,kr,er)
                        hi = psiq(n,iquad)

                        qr(1,k) = qr(1,k) + hi*qprime_df(1,I,k)
                        qr(2,k) = qr(2,k) + hi*qprime_df(2,I,k)
                        qr(3,k) = qr(3,k) + hi*qprime_df(3,I,k)
                    enddo
                else 
                    qr(:,k) = ql(:,k)
                    ! Apply wall boundary conditions
                    if (er == -4) then
                        un = ql(2,k)*nxl + ql(3,k)*nyl
                        qr(2,k) = ql(2,k) - 2.0*un*nxl
                        qr(3,k) = ql(3,k) - 2.0*un*nyl
                    elseif (er == -2) then
                        qr(2,k) = -ql(2,k)
                        qr(3,k) = -ql(3,k)
                    endif

                endif 

                ! Compute the fluxes
                dpl = qbl(1,iquad) * ql(1,k)
                dpr = qbr(1,iquad) * qr(1,k)
                ul = ql(2,k) + qbl(2,iquad)
                ur = qr(2,k) + qbr(2,iquad)
                vl = ql(3,k) + qbl(3,iquad)
                vr = qr(3,k) + qbr(3,iquad)

                if (ql(1,k) <= (gravity/alpha_mlswe(k)) * dry_cutoff) then
                    ql(1,k) = (gravity/alpha_mlswe(k)) * dry_cutoff
                    dpl = qbl(1,iquad) * ql(1,k)
                    ul = 0.0
                    vl = 0.0
                end if

                if (qr(1,k) <= (gravity/alpha_mlswe(k)) * dry_cutoff) then
                    qr(1,k) = (gravity/alpha_mlswe(k)) * dry_cutoff
                    dpr = qbr(1,iquad) * qr(1,k)
                    ur = 0.0
                    vr = 0.0
                end if

                uu = 0.5*(ul+ur)
                vv = 0.5*(vl+vr)
                udpl(k) = ul*dpl
                udpr(k) = ur*dpr
                vdpl(k) = vl*dpl
                vdpr(k) = vr*dpr

                if(uu*nxl > 0.0) then
                    udp_flux(1,k,iquad) = uu * (ul*dpl)
                    vdp_flux(1,k,iquad) = uu * (vl*dpl)
                else
                    udp_flux(1,k,iquad) = uu * (ur*dpr)
                    vdp_flux(1,k,iquad) = uu * (vr*dpr)
                endif
                if(vv*nyl > 0.0) then
                    udp_flux(2,k,iquad) = vv * (ul*dpl)
                    vdp_flux(2,k,iquad) = vv * (vl*dpl)
                else
                    udp_flux(2,k,iquad) = vv * (ur*dpr)
                    vdp_flux(2,k,iquad) = vv * (vr*dpr)
                endif

            enddo

            uu_dp_flux_deficit(1) = Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,:,iquad))
            uu_dp_flux_deficit(2) = Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,:,iquad))
            vv_dp_flux_deficit(1) = Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,:,iquad))
            vv_dp_flux_deficit(2) = Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,:,iquad))

            do k = 1,nlayers

                ! Adjust the fluxes for the u-momentum equation
                !x-direction
                weight = abs(udpl(k)) / (sum(abs(udpl(:))+eps1))
                if(uu_dp_flux_deficit(1)*nxl < 0.0) &
                    weight = abs(udpr(k)) / (sum(abs(udpr(:))+eps1))
                udp_flux(1,k,iquad) = udp_flux(1,k,iquad) + weight * uu_dp_flux_deficit(1)

                !y-direction
                weight = abs(udpl(k)) / (sum(abs(udpl(:))+eps1))
                if(uu_dp_flux_deficit(2)*nyl < 0.0) &
                    weight = abs(udpr(k)) / (sum(abs(udpr(:))+eps1))
                udp_flux(2,k,iquad) = udp_flux(2,k,iquad) + weight * uu_dp_flux_deficit(2)

                ! Adjust the fluxes for the v-momentum equation
                !x-direction
                weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
                if(vv_dp_flux_deficit(1)*nxl < 0.0) &
                    weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
                vdp_flux(1,k,iquad) = vdp_flux(1,k,iquad) + weight * vv_dp_flux_deficit(1)

                !y-direction
                weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
                if(vv_dp_flux_deficit(2)*nyl < 0.0) &
                    weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
                vdp_flux(2,k,iquad) = vdp_flux(2,k,iquad) + weight * vv_dp_flux_deficit(2)
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
                p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,k)
                p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,k)
            end do

            one_plus_eta_edge = sqrt(one_plus_eta_edge_2_ave(iquad,iface))
            z_face(1,nlayers+1) = zbot_face(1,iquad,iface)
            z_face(2,nlayers+1) = zbot_face(2,iquad,iface)
            z_edge_plus(nlayers+1) = zbot_face(1,iquad,iface)
            z_edge_minus(nlayers+1) = zbot_face(2,iquad,iface)
            do k=nlayers,1,-1
                z_face(1,k) = z_face(1,k+1) + alpha_over_g(k) * (ope_l * ql(1,k))
                z_face(2,k) = z_face(2,k+1) + alpha_over_g(k) * (ope_r * qr(1,k))
                z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k) * &
                                                    (one_plus_eta_edge * ql(1,k))
                z_edge_minus(k) = z_edge_minus(k+1) + alpha_over_g(k) * &
                                                    (one_plus_eta_edge * qr(1,k))
            end do

            p_edge_plus(2) = one_plus_eta_edge * ql(1,1)
            p_edge_minus(2) = one_plus_eta_edge * qr(1,1)
            do k = 2,nlayers
                p_edge_plus(k+1) = p_edge_plus(k) + one_plus_eta_edge * ql(1,k)
                p_edge_minus(k+1) = p_edge_minus(k) + one_plus_eta_edge * qr(1,k)
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
                H_face(1,k,iquad) = 0.5*(H_r_plus + H_r_minus) !computation of H_r for the left side
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
                H_face(2,k,iquad) = 0.5*(H_r_plus + H_r_minus) ! H_r for the right side
            end do !k

            ! Wall Boundary conditions
            if(er == -4) then
                p2l = 0.0 ; p2r = 0.0
                do k = 1,nlayers
                    p2l(k+1) = p_face(1,k+1)
                    H_face(1,k,iquad) = 0.5*alpha_mlswe(k)*(p2l(k+1)**2 - p2l(k)**2)
                    p2r(k+1) = p_face(2,k+1)
                    H_face(2,k,iquad) = 0.5*alpha_mlswe(k)*(p2r(k+1)**2 - p2r(k)**2)

                end do
            end if

            if(er /= -4) then
                do k = 1, nlayers-1          ! interface at the bottom of layer k
                    ! Corrections at the left side of a face.
                    p_inc1 = g_over_alpha(k)*(z_face(1,k+1) - z_edge_plus(k+1))
                    H_corr1 = 0.5 * alpha_mlswe(k) * ((p_face(1,k+1) &
                                + p_inc1)**2 - p_face(1,k+1)**2)
                    H_face(1,k,iquad) = H_face(1,k,iquad) - H_corr1
                    H_face(1,k+1,iquad) = H_face(1,k+1,iquad) + H_corr1

                    ! Corrections at the right side of a face.
                    p_inc2 = g_over_alpha(k)*(z_face(2,k+1) - z_edge_minus(k+1))
                    H_corr2 = 0.5 * alpha_mlswe(k) * ((p_face(2,k+1) &
                                + p_inc2)**2 - p_face(2,k+1)**2)
                    H_face(2,k,iquad) = H_face(2,k,iquad) - H_corr2
                    H_face(2,k+1,iquad) = H_face(2,k+1,iquad) + H_corr2

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
            acceleration = sum(H_face(1,:,iquad))
            if(acceleration > 0.0) then
                weight = H_face_ave(iquad,iface) / acceleration
            end if
            H_face(1,:,iquad) = H_face(1,:,iquad) * weight

            ! Right side of face
            weight = 1.0
            acceleration = sum(H_face(2,:,iquad))
            if(acceleration > 0.0) then
                weight = H_face_ave(iquad,iface) / acceleration
            end if
            H_face(2,:,iquad) = H_face(2,:,iquad) * weight
            
            ! Do Gauss-Lobatto Integration
            wq = jac_faceq(iquad,1,iface)

            do k = 1,nlayers

                flux_xl = nxl*(udp_flux(1,k,iquad) + H_face(1,k,iquad)) + nyl*udp_flux(2,k,iquad)
                flux_xr = nxl*(udp_flux(1,k,iquad) + H_face(2,k,iquad)) + nyl*udp_flux(2,k,iquad)

                flux_yl = nxl*vdp_flux(1,k,iquad) + nyl*(vdp_flux(2,k,iquad) + H_face(1,k,iquad))
                flux_yr = nxl*vdp_flux(1,k,iquad) + nyl*(vdp_flux(2,k,iquad) + H_face(2,k,iquad))

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)
                    I = intma(il,jl,kl,el)

                    rhs_mom(1,I,k) = rhs_mom(1,I,k) - wq*hi*flux_xl
                    rhs_mom(2,I,k) = rhs_mom(2,I,k) - wq*hi*flux_yl

                enddo

                if(er > 0) then

                    do n = 1, ngl

                        hi = psiq(n,iquad)
                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)
                        I = intma(ir,jr,kr,er)

                        rhs_mom(1,I,k) = rhs_mom(1,I,k) + wq*hi*flux_xr
                        rhs_mom(2,I,k) = rhs_mom(2,I,k) + wq*hi*flux_yr

                    enddo
                end if
            end do
        end do ! iface

    end subroutine Apply_layers_fluxes

    subroutine Apply_bcl_fluxes(rhs, qprime_df)

        ! This routine computes the layer momentum advection flux terms using upwind flux
        ! This routine computes the layer momentum pressure terms

        use mod_constants, only : gravity
        use mod_initial, only : alpha_mlswe, zbot_face
        use mod_grid, only : nface, npoin, npoin_q, face, intma, face_type
        use mod_basis, only : nq, psiq, ngl
        use mod_input, only : nlayers, dry_cutoff
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
        real, dimension(3, nlayers) :: ql, qr
        real, dimension(3) :: qbl, qbr
        integer :: iface, ilr, k, iquad, ktemp, I
        real :: z_intersect_top,z_intersect_bot, dz_intersect, H_r_plus, H_r_minus, acceleration
        real :: p_intersect_bot, p_intersect_top, one_plus_eta_edge
        real :: H_corr,p_inc, weight, H_corr1,p_inc1, H_corr2,p_inc2, temp, ope_l, ope_r
        integer :: Iq, el, er
        real ::  ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv
        real, dimension(nlayers) :: udpl, udpr, vdpl, vdpr
        real :: uu_dp_flux_deficit(2), vv_dp_flux_deficit(2)
        real, parameter :: eps1 = 1.0e-20 !  Parameter used to prevent division by zero.
        integer :: il, jl, ir, jr, kl, kr, jquad, n, m
        real :: wq, hi, hlx_k, hly_k, hrx_k, hry_k, flux_x, flux_y, hx_k, hy_k, flux, un
        real, dimension(2,nq,nlayers) :: H_face, udp_flux, vdp_flux, dp_flux
        real :: dp_lr(2,nlayers), dp_deficit(2)
        real, dimension(nq,nlayers) :: flux_ul, flux_ur, flux_vl, flux_vr, flux_dp
        real :: flux_xl, flux_xr, flux_yl, flux_yr

        do k=1,nlayers
            alpha_over_g(k) = alpha_mlswe(k)/gravity
            g_over_alpha(k) = gravity/alpha_mlswe(k)
        enddo

        ! Compute H_r at the element face
        do concurrent(iface = 1:nface, iquad = 1:nq)

            if (face_type(iface) == 2) cycle

            ! Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)


            qbl(1) = ope_face_ave(1,iquad,iface)
            qbl(2) = uvb_face_ave(1,1,iquad,iface)
            qbl(3) = uvb_face_ave(2,1,iquad,iface)
            qbr(1) = ope_face_ave(2,iquad,iface)
            qbr(2) = uvb_face_ave(1,2,iquad,iface)
            qbr(3) = uvb_face_ave(2,2,iquad,iface)

            nxl = normal_vector_q(1,iquad,1,iface)
            nyl = normal_vector_q(2,iquad,1,iface)

            ql = 0.0; qr = 0.0
            do k = 1,nlayers

                do n = 1, ngl

                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)
                    I = intma(il,jl,kl,el)
                    hi = psiq(n,iquad)

                    ql(:,k) = ql(:,k) + hi*qprime_df(:,I,k)
                enddo

                if (er > 0) then
                    do n = 1, ngl

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)
                        I = intma(ir,jr,kr,er)
                        hi = psiq(n,iquad)

                        qr(:,k) = qr(:,k) + hi*qprime_df(:,I,k)
                    enddo
                else 
                    qr(:,k) = ql(:,k)
                    ! Apply wall boundary conditions
                    if (er == -4) then
                        un = ql(2,k)*nxl + ql(3,k)*nyl
                        qr(2,k) = ql(2,k) - 2.0*un*nxl
                        qr(3,k) = ql(3,k) - 2.0*un*nyl
                    elseif (er == -2) then
                        qr(2,k) = -ql(2,k)
                        qr(3,k) = -ql(3,k)
                    endif
                endif

                ! Left side of the edge
                dpl = qbl(1) * ql(1,k)
                dpr = qbr(1) * qr(1,k)
                ul = ql(2,k) + qbl(2)
                ur = qr(2,k) + qbr(2)
                vl = ql(3,k) + qbl(3)
                vr = qr(3,k) + qbr(3)

                if (ql(1,k) <= (gravity/alpha_mlswe(k)) * dry_cutoff) then
                    ql(1,k) = (gravity/alpha_mlswe(k)) * dry_cutoff
                    dpl = qbl(1) * ql(1,k)
                    ul = 0.0
                    vl = 0.0
                end if

                if (qr(1,k) <= (gravity/alpha_mlswe(k)) * dry_cutoff) then
                    qr(1,k) = (gravity/alpha_mlswe(k)) * dry_cutoff
                    dpr = qbr(1) * qr(1,k)
                    ur = 0.0
                    vr = 0.0
                end if

                ! if ((ql(1,k) <= (gravity/alpha_mlswe(k)) * dry_cutoff) .or. (qr(1,k) <= (gravity/alpha_mlswe(k)) * dry_cutoff)) then
                !     ql(1,k) = (gravity/alpha_mlswe(k)) * dry_cutoff
                !     dpl = qbl(1) * ql(1,k)
                !     ul = 0.0
                !     vl = 0.0

                !     qr(1,k) = (gravity/alpha_mlswe(k)) * dry_cutoff
                !     dpr = qbr(1) * qr(1,k)
                !     ur = 0.0
                !     vr = 0.0
                ! end if

                dp_lr(1,k) = dpl
                dp_lr(2,k) = dpr

                uu = 0.5*(ul+ur)
                vv = 0.5*(vl+vr)
                udpl(k) = ul*dpl
                udpr(k) = ur*dpr
                vdpl(k) = vl*dpl
                vdpr(k) = vr*dpr

                if(uu*nxl > 0.0) then
                    dp_flux(1,iquad,k) = uu * dpl
                    udp_flux(1,iquad,k) = uu * (ul*dpl)
                    vdp_flux(1,iquad,k) = uu * (vl*dpl)
                else
                    dp_flux(1,iquad,k) = uu * dpr
                    udp_flux(1,iquad,k) = uu * (ur*dpr)
                    vdp_flux(1,iquad,k) = uu * (vr*dpr)
                endif
                if(vv*nyl > 0.0) then
                    dp_flux(2,iquad,k) = vv * dpl
                    udp_flux(2,iquad,k) = vv * (ul*dpl)
                    vdp_flux(2,iquad,k) = vv * (vl*dpl)
                else
                    dp_flux(2,iquad,k) = vv * dpr
                    udp_flux(2,iquad,k) = vv * (ur*dpr)
                    vdp_flux(2,iquad,k) = vv * (vr*dpr)
                endif

            enddo

            uu_dp_flux_deficit(1) = Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,iquad,:))
            uu_dp_flux_deficit(2) = Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,iquad,:))
            vv_dp_flux_deficit(1) = Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,iquad,:))
            vv_dp_flux_deficit(2) = Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,iquad,:))

            dp_deficit(1) = btp_mass_flux_face_ave(1,iquad,iface) - sum(dp_flux(1,iquad,:))
            dp_deficit(2) = btp_mass_flux_face_ave(2,iquad,iface) - sum(dp_flux(2,iquad,:))

            do k = 1,nlayers

                weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:))+eps1))
                if (dp_deficit(1)*nxl < 0.0) &
                    weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:))+eps1))
                dp_flux(1,iquad,k) = dp_flux(1,iquad,k) + weight * dp_deficit(1)

                weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:))+eps1))
                if (dp_deficit(2)*nyl < 0.0) &
                    weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:))+eps1))
                dp_flux(2,iquad,k) = dp_flux(2,iquad,k) + weight * dp_deficit(2)

                ! Adjust the fluxes for the u-momentum equation
                !x-direction
                weight = abs(udpl(k)) / (sum(abs(udpl(:))+eps1))
                if(uu_dp_flux_deficit(1)*nxl < 0.0) &
                    weight = abs(udpr(k)) / (sum(abs(udpr(:))+eps1))
                udp_flux(1,iquad,k) = udp_flux(1,iquad,k) + weight * uu_dp_flux_deficit(1)

                !y-direction
                weight = abs(udpl(k)) / (sum(abs(udpl(:))+eps1))
                if(uu_dp_flux_deficit(2)*nyl < 0.0) &
                    weight = abs(udpr(k)) / (sum(abs(udpr(:))+eps1))
                udp_flux(2,iquad,k) = udp_flux(2,iquad,k) + weight * uu_dp_flux_deficit(2)

                ! Adjust the fluxes for the v-momentum equation
                !x-direction
                weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
                if(vv_dp_flux_deficit(1)*nxl < 0.0) &
                    weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
                vdp_flux(1,iquad,k) = vdp_flux(1,iquad,k) + weight * vv_dp_flux_deficit(1)

                !y-direction
                weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
                if(vv_dp_flux_deficit(2)*nyl < 0.0) &
                    weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
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
                p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,k)
                p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,k)
            end do

            one_plus_eta_edge = sqrt(one_plus_eta_edge_2_ave(iquad,iface))
            z_face(1,nlayers+1) = zbot_face(1,iquad,iface)
            z_face(2,nlayers+1) = zbot_face(2,iquad,iface)
            z_edge_plus(nlayers+1) = zbot_face(1,iquad,iface)
            z_edge_minus(nlayers+1) = zbot_face(2,iquad,iface)
            do k=nlayers,1,-1
                z_face(1,k) = z_face(1,k+1) + alpha_over_g(k) * (ope_l * ql(1,k))
                z_face(2,k) = z_face(2,k+1) + alpha_over_g(k) * (ope_r * qr(1,k))
                z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k) * &
                                                    (one_plus_eta_edge * ql(1,k))
                z_edge_minus(k) = z_edge_minus(k+1) + alpha_over_g(k) * &
                                                    (one_plus_eta_edge * qr(1,k))
            end do

            p_edge_plus(2) = one_plus_eta_edge * ql(1,1)
            p_edge_minus(2) = one_plus_eta_edge * qr(1,1)
            do k = 2,nlayers
                p_edge_plus(k+1) = p_edge_plus(k) + one_plus_eta_edge * ql(1,k)
                p_edge_minus(k+1) = p_edge_minus(k) + one_plus_eta_edge * qr(1,k)
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
        
            ! Do Gauss-Lobatto Integration

            wq = jac_faceq(iquad,1,iface)

            do k = 1,nlayers

                flux = nxl*dp_flux(1,iquad,k) + nyl*dp_flux(2,iquad,k)
                flux_xl = nxl*(udp_flux(1,iquad,k) + H_face(1,iquad,k)) + nyl*udp_flux(2,iquad,k)
                flux_yl = nxl*vdp_flux(1,iquad,k) + nyl*(vdp_flux(2,iquad,k) + H_face(1,iquad,k))

                flux_xr = nxl*(udp_flux(1,iquad,k) + H_face(2,iquad,k)) + nyl*udp_flux(2,iquad,k)
                flux_yr = nxl*vdp_flux(1,iquad,k) + nyl*(vdp_flux(2,iquad,k) + H_face(2,iquad,k))

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)
                    I = intma(il,jl,kl,el)

                    rhs(1,I,k) = rhs(1,I,k) - wq*hi*flux
                    rhs(2,I,k) = rhs(2,I,k) - wq*hi*flux_xl
                    rhs(3,I,k) = rhs(3,I,k) - wq*hi*flux_yl

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)
                        I = intma(ir,jr,kr,er)

                        rhs(1,I,k) = rhs(1,I,k) + wq*hi*flux
                        rhs(2,I,k) = rhs(2,I,k) + wq*hi*flux_xr
                        rhs(3,I,k) = rhs(3,I,k) + wq*hi*flux_yr

                    end if
                end do
            end do
        end do ! iface, iquad

    end subroutine Apply_bcl_fluxes

    subroutine create_layers_volume_mass(dp_advec, qprime_df)

        use mod_grid, only : npoin_q, npoin, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_input, only: nlayers, dry_cutoff
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, alpha_mlswe
        use mod_variables, only: uvb_ave, ope_ave, sum_layer_mass_flux, btp_mass_flux_ave

        implicit none

        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(npoin,nlayers), intent(out) :: dp_advec

        real :: wq, hi, dhde, dhdn
        integer :: k, I, Iq, ip
        real :: dp_temp, dp, u, v, opeq, flux(2), sum_dpp, sum_udp, sum_vdp
        real, dimension(3) :: qp, qb
        real, parameter :: eps = 1.0e-10
        real, dimension(nlayers) :: dpp, udp, vdp

        dp_advec = 0.0
        sum_layer_mass_flux = 0.0

        do concurrent (Iq = 1:npoin_q)

            qb(1) = ope_ave(Iq)
            qb(2) = uvb_ave(1,Iq)
            qb(3) = uvb_ave(2,Iq)
            wq    = wjac(Iq)

            ! 1) Build layer arrays (private)
            do k = 1, nlayers

                qp = 0.0
                do ip = 1, npts
                    I  = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    qp(1) = qp(1) + hi * qprime_df(1, I, k)
                    qp(2) = qp(2) + hi * qprime_df(2, I, k)
                    qp(3) = qp(3) + hi * qprime_df(3, I, k)
                end do

                dp_temp   = qp(1) * qb(1)
                dpp(k)  = dp_temp
                udp(k)  = (qp(2) + qb(2)) * dp_temp
                vdp(k)  = (qp(3) + qb(3)) * dp_temp

                if (qp(1) <= (gravity/alpha_mlswe(k)) * dry_cutoff) then
                    qp(1) = (gravity/alpha_mlswe(k)) * dry_cutoff
                    dp_temp = qp(1) * qb(1)
                    dpp(k)  = dp_temp
                    udp(k)  = 0.0
                    vdp(k)  = 0.0
                end if

            end do

            ! 2) Precompute sums once
            sum_dpp = sum(dpp(:))
            sum_udp = sum(udp(:))
            sum_vdp = sum(vdp(:))

            ! 3) Flux and advection update
            do k = 1, nlayers
                flux(1) = udp(k) + (dpp(k)/sum_dpp) * (btp_mass_flux_ave(1,Iq) - sum_udp)
                flux(2) = vdp(k) + (dpp(k)/sum_dpp) * (btp_mass_flux_ave(2,Iq) - sum_vdp)

                do ip = 1, npts
                    I = indexq(ip,Iq)
                    dp_advec(I,k) = dp_advec(I,k) + wq * ( dpsidx(ip,Iq)*flux(1) + dpsidy(ip,Iq)*flux(2) )
                end do
            end do

        end do ! Iq

    end subroutine create_layers_volume_mass

    subroutine create_layer_mass_flux(dp_advec, qprime_df)

        use mod_basis, only: nq, psiq, ngl
        use mod_grid, only:  npoin_q, intma,  nface, face, face_type
        use mod_input, only: nlayers, dry_cutoff
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: sum_layer_mass_flux_face, ope_face_ave, uvb_face_ave, btp_mass_flux_face_ave
        use mod_constants, only : gravity
        use mod_initial, only : alpha_mlswe

        implicit none

        real, dimension(npoin, nlayers), intent(inout) :: dp_advec
        real, dimension(3, npoin, nlayers), intent(in) :: qprime_df

        integer :: k, iface, iquad, el, er, il, jl, ir, jr, I
        integer :: kl, kr, jquad, n, m
        real :: wq, nxl, nyl, hi
        real, dimension(nlayers,nq) :: flux_edge_u, flux_edge_v, flux_dp
        real :: dpl, dpr, uu, vv, flux, un, ul, ur, vl, vr, weight, dp_lr(2,nlayers), dp_deficit(2)
        real, dimension(3) :: ql, qr
        real, dimension(3,nq) :: qbl, qbr
        real, parameter :: eps = 1.0e-10

        do concurrent (iface = 1:nface, iquad = 1:nq)

            if (face_type(iface) == 2) cycle

            el = face(7,iface)
            er = face(8,iface)

            qbl(1,iquad) = ope_face_ave(1,iquad,iface)
            qbl(2,iquad) = uvb_face_ave(1,1,iquad,iface)
            qbl(3,iquad) = uvb_face_ave(2,1,iquad,iface)

            qbr(1,iquad) = ope_face_ave(2,iquad,iface)
            qbr(2,iquad) = uvb_face_ave(1,2,iquad,iface)
            qbr(3,iquad) = uvb_face_ave(2,2,iquad,iface)

            !---------------------------------------
            ! Build flux_edge_u/v for each iquad,k
            !---------------------------------------

            nxl = normal_vector_q(1,iquad,1,iface)
            nyl = normal_vector_q(2,iquad,1,iface)

            do k = 1, nlayers

                ql = 0.0 ; qr = 0.0
                do n = 1, ngl
                    il = imapl(1,n,1,iface); jl = imapl(2,n,1,iface); kl = imapl(3,n,1,iface)
                    I  = intma(il,jl,kl,el)
                    hi = psiq(n,iquad)
                    ql(1) = ql(1) + hi*qprime_df(1,I,k)
                    ql(2) = ql(2) + hi*qprime_df(2,I,k)
                    ql(3) = ql(3) + hi*qprime_df(3,I,k)
                end do

                if (er > 0) then
                    do n = 1, ngl
                        ir = imapr(1,n,1,iface); jr = imapr(2,n,1,iface); kr = imapr(3,n,1,iface)
                        I  = intma(ir,jr,kr,er)
                        hi = psiq(n,iquad)
                        qr(1) = qr(1) + hi*qprime_df(1,I,k)
                        qr(2) = qr(2) + hi*qprime_df(2,I,k)
                        qr(3) = qr(3) + hi*qprime_df(3,I,k)
                    end do
                else
                    qr(:) = ql(:)
                    if (er == -4) then
                        un = ql(2)*nxl + ql(3)*nyl
                        qr(2) = ql(2) - 2.0*un*nxl
                        qr(3) = ql(3) - 2.0*un*nyl
                    else if (er == -2) then
                        qr(2) = -ql(2)
                        qr(3) = -ql(3)
                    end if
                end if

                dpl = qbl(1,iquad) * ql(1)
                dpr = qbr(1,iquad) * qr(1)

                dp_lr(1,k) = dpl
                dp_lr(2,k) = dpr

                ul = ql(2) + qbl(2,iquad)
                ur = qr(2) + qbr(2,iquad)
                vl = ql(3) + qbl(3,iquad)
                vr = qr(3) + qbr(3,iquad)

                if (ql(1) <= (gravity/alpha_mlswe(k)) * dry_cutoff) then
                    ql(1) = (gravity/alpha_mlswe(k)) * dry_cutoff
                    dpl = qbl(1,iquad) * ql(1)
                    dp_lr(1,k) = dpl
                    ul = 0.0
                    vl = 0.0
                end if

                if (qr(1) <= (gravity/alpha_mlswe(k)) * dry_cutoff) then
                    qr(1) = (gravity/alpha_mlswe(k)) * dry_cutoff
                    dpr = qbr(1,iquad) * qr(1)
                    dp_lr(2,k) = dpr
                    ur = 0.0
                    vr = 0.0
                end if

                ! if ((ql(1) <= (gravity/alpha_mlswe(k)) * dry_cutoff) .or. (qr(1) <= (gravity/alpha_mlswe(k)) * dry_cutoff)) then
                !     ql(1) = (gravity/alpha_mlswe(k)) * dry_cutoff
                !     dpl = qbl(1,iquad) * ql(1)
                !     dp_lr(1,k) = dpl
                !     ul = 0.0
                !     vl = 0.0

                !     qr(1) = (gravity/alpha_mlswe(k)) * dry_cutoff
                !     dpr = qbr(1,iquad) * qr(1)
                !     dp_lr(2,k) = dpr
                !     ur = 0.0
                !     vr = 0.0
                ! end if

                uu = 0.5*(ul + ur)
                vv = 0.5*(vl + vr)

                if (uu*nxl > 0.0) then
                    flux_edge_u(k,iquad) = uu * dpl
                else
                    flux_edge_u(k,iquad) = uu * dpr
                end if

                if (vv*nyl > 0.0) then
                    flux_edge_v(k,iquad) = vv * dpl
                else
                    flux_edge_v(k,iquad) = vv * dpr
                end if

            end do ! k

            dp_deficit(1) = btp_mass_flux_face_ave(1,iquad,iface) - sum(flux_edge_u(:,iquad))
            dp_deficit(2) = btp_mass_flux_face_ave(2,iquad,iface) - sum(flux_edge_v(:,iquad))

            do k = 1, nlayers

                weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:)) + eps))
                if (dp_deficit(1)*nxl < 0.0) weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:)) + eps))
                flux_edge_u(k,iquad) = flux_edge_u(k,iquad) + weight*dp_deficit(1)

                weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:)) + eps))
                if (dp_deficit(1)*nxl < 0.0) weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:)) + eps))
                flux_edge_v(k,iquad) = flux_edge_v(k,iquad) + weight*dp_deficit(2)

            end do ! k

            !---------------------------------------
            ! Update dp_advec (RACE across faces unless disjoint I)
            !---------------------------------------

            wq  = jac_faceq(iquad,1,iface)
 
            do k = 1, nlayers

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                flux = nxl*flux_edge_u(k,iquad) + nyl*flux_edge_v(k,iquad)

                do n = 1, ngl
                    hi = psiq(n,iquad)
                    il = imapl(1,n,1,iface); jl = imapl(2,n,1,iface); kl = imapl(3,n,1,iface)
                    I  = intma(il,jl,kl,el)
                    dp_advec(I,k) = dp_advec(I,k) - wq*hi*flux
                end do

                if (er > 0) then
                    do n = 1, ngl
                        hi = psiq(n,iquad)
                        ir = imapr(1,n,1,iface); jr = imapr(2,n,1,iface); kr = imapr(3,n,1,iface)
                        I  = intma(ir,jr,kr,er)
                        dp_advec(I,k) = dp_advec(I,k) + wq*hi*flux
                    end do
                end if

            end do
        end do ! iface, iquad

    end subroutine create_layer_mass_flux

end module mod_create_rhs_mlswe
