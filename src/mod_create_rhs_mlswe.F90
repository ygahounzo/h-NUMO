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

    use mod_metrics, only: &
        ksiq_x, ksiq_y, ksiq_z, &
        etaq_x, etaq_y, etaq_z, &
        zetaq_x, zetaq_y, zetaq_z, &
        jacq, massinv

    public :: layer_momentum_rhs, btp_mass_advection_rhs, &
              interpolate_layer_from_quad_to_node, rhs_layer_shear_stress, layer_mass_rhs, &
              consistency_mass_rhs
              
contains


    subroutine layer_momentum_rhs(rhs_mom, rhs_visc, qprime, qprime_face)

        use mod_input, only: nlayers
        use mod_metrics, only: massinv
        use mod_grid, only: npoin, npoin_q, intma, intma_dg_quad, nface
        use mod_basis, only: nq

        implicit none
        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom
        real, dimension(2,npoin,nlayers), intent(in) :: rhs_visc
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime 
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face 

        integer :: k,I

        call create_rhs_dynamics_volume_layers(rhs_mom, qprime)

        call Apply_layers_fluxes(rhs_mom, qprime_face)

        do k = 1, nlayers

            rhs_mom(1,:,k) = massinv(:)*rhs_mom(1,:,k) + rhs_visc(1,:,k)
            rhs_mom(2,:,k) = massinv(:)*rhs_mom(2,:,k) + rhs_visc(2,:,k)
        end do

    end subroutine layer_momentum_rhs

    subroutine layer_mass_rhs(dp_advec, qprime, qprime_face)

        use mod_metrics, only: massinv
        use mod_input, only: nlayers
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq

        implicit none

        real, dimension(npoin, nlayers), intent(out) :: dp_advec
        real, dimension(3, npoin_q, nlayers), intent(in) :: qprime
        real, dimension(3, 2, nq, nface, nlayers), intent(in) :: qprime_face
        
        integer :: k

        ! Compute the mass advection term for the degree of freedom for dp in each layer

        call create_layers_volume_mass(dp_advec, qprime)

        ! Compute the mass flux term 

        call create_layer_mass_flux(dp_advec, qprime_face)

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

        !do k = 1, nlayers

        !    dp_advec(:,k) = massinv(:)*dp_advec(:,k)
        !end do
        
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

    subroutine create_rhs_dynamics_volume_layers(rhs_mom, qprime)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts, dpsiqx
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, alpha_mlswe, tau_wind

        use mod_variables, only: H_r, u_udp_temp, v_vdp_temp, p, u_vdp_temp, grad_z, &
                    tau_wind_int, tau_bot_int, tau_bot_ave, ope_ave, uvb_ave

        implicit none 

        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime 

        real :: wq, hi, dhdx, dhdy, bot_layer, tau_wind_u, tau_wind_v, temp1
        real :: Hq, var_uu, var_uv, var_vu, var_vv, source_x, source_y, Pstress

        integer :: k, I, Iq, ip
        real, dimension(npoin_q,nlayers+1) :: pprime_temp 
        real :: temp_dp, temp_u, temp_v, Ptop_k, Pbot_k, tempbot, Pbstress

        rhs_mom = 0.0
        bot_layer = 0.0

        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth at which bottom stress is reduced to 0
        do k = 1,nlayers
        
            pprime_temp(:, k+1) = pprime_temp(:, k) + qprime(1,:,k)

            if(k == nlayers) bot_layer = 1.0

            do Iq = 1, npoin_q

                wq = wjac(Iq)

                Hq = H_r(Iq,k)

                !temp_dp = qprime(1,Iq,k) * ope_ave(Iq)
                !temp_u = qprime(2,Iq,k) + uvb_ave(1,Iq)
                !temp_v = qprime(3,Iq,k) + uvb_ave(2,Iq)

                !var_uu = temp_dp * temp_u**2  
                !var_uv = temp_u * temp_v * temp_dp
                !var_vu = temp_v * temp_u * temp_dp
                !var_vv = temp_dp * temp_v**2

                var_uu = u_udp_temp(Iq,k)   
                var_uv = u_vdp_temp(1,Iq,k)
                var_vu = u_vdp_temp(2,Iq,k)
                var_vv = v_vdp_temp(Iq,k)  

                temp1 = (min(pprime_temp(Iq,k+1), Pstress) - min(pprime_temp(Iq,k), Pstress))/ Pstress
                tau_wind_u = temp1*tau_wind(1,Iq)
                tau_wind_v = temp1*tau_wind(2,Iq)
                
                Ptop_k = pprime_temp(Iq,k)
                Pbot_k = pprime_temp(Iq,k+1)
                tempbot = max(pprime_temp(Iq,nlayers+1)-Pbstress,Pbot_k)-max(pprime_temp(Iq,nlayers+1)-Pbstress,min(Ptop_k,Pbot_k))
                
                tempbot = 1.0 !tempbot/qprime(1,Iq,k)

                !source_x = gravity*(tau_wind_u - bot_layer*tau_bot_ave(1,Iq) + &
                !            p(Iq,k) * grad_z(1,Iq,k) - p(Iq,k+1) * grad_z(1,Iq,k+1))
                !source_y = gravity*(tau_wind_v - bot_layer*tau_bot_ave(2,Iq) + &
                !            p(Iq,k) * grad_z(2,Iq,k) - p(Iq,k+1) * grad_z(2,Iq,k+1))                                                                                                                      
                source_x = gravity*(tau_wind_u - tempbot*tau_bot_ave(1,Iq) + &
                            p(Iq,k) * grad_z(1,Iq,k) - p(Iq,k+1) * grad_z(1,Iq,k+1))
                source_y = gravity*(tau_wind_v - tempbot*tau_bot_ave(2,Iq) + &
                            p(Iq,k) * grad_z(2,Iq,k) - p(Iq,k+1) * grad_z(2,Iq,k+1))                                                                                                                      
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

    subroutine create_rhs_dynamics_volume_layers_v2(rhs_mom, qprime)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts, dpsiqx
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, alpha_mlswe, tau_wind

        use mod_variables, only: H_r,u_udp_temp, v_vdp_temp, p, u_vdp_temp, grad_z, &
                    tau_wind_int, tau_bot_int, tau_bot_ave, ope_ave, uvb_ave

        implicit none

        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime

        real :: wq, hi, dhdx, dhdy, bot_layer, tau_wind_u, tau_wind_v, temp1
        real :: Hq, var_uu, var_uv, var_vu, var_vv, source_x, source_y, Pstress

        integer :: k, I, Iq, ip
        real, dimension(npoin_q,nlayers+1) :: pprime_temp
        real :: temp_dp, temp_u, temp_v, Ptop_k, Pbot_k, tempbot, Pbstress

        rhs_mom = 0.0
        bot_layer = 0.0

        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth at which wind stress is reduced to 0
        Pbstress = (gravity/alpha_mlswe(nlayers)) * 10.0 ! pressure corresponding to 10m depth at which bottom stress is reduced to 0
        do k = 1,nlayers

            pprime_temp(:, k+1) = pprime_temp(:, k) + qprime(1,:,k)

            if(k == nlayers) bot_layer = 1.0

            do Iq = 1, npoin_q

                wq = wjac(Iq)

                Hq = H_r(Iq,k)

                temp_dp = qprime(1,Iq,k) * ope_ave(Iq)
                temp_u = qprime(2,Iq,k) + uvb_ave(1,Iq)
                temp_v = qprime(3,Iq,k) + uvb_ave(2,Iq)

                var_uu = temp_dp * temp_u**2
                var_uv = temp_u * temp_v * temp_dp
                var_vu = temp_v * temp_u * temp_dp
                var_vv = temp_dp * temp_v**2

                temp1 = (min(pprime_temp(Iq,k+1), Pstress) - min(pprime_temp(Iq,k), Pstress))/ Pstress
                tau_wind_u = temp1*tau_wind(1,Iq)
                tau_wind_v = temp1*tau_wind(2,Iq)

                Ptop_k = pprime_temp(Iq,k)
                Pbot_k = pprime_temp(Iq,k+1)
                tempbot = max(pprime_temp(Iq,nlayers+1)-Pbstress,Pbot_k)-max(pprime_temp(Iq,nlayers+1)-Pbstress,min(Ptop_k,Pbot_k))

                tempbot = 1.0 !tempbot/qprime(1,Iq,k)

                !source_x = gravity*(tau_wind_u - bot_layer*tau_bot_ave(1,Iq) + &
                !            p(Iq,k) * grad_z(1,Iq,k) - p(Iq,k+1) * grad_z(1,Iq,k+1))
                !source_y = gravity*(tau_wind_v - bot_layer*tau_bot_ave(2,Iq) + &
                !            p(Iq,k) * grad_z(2,Iq,k) - p(Iq,k+1) * grad_z(2,Iq,k+1))
                source_x = gravity*(tau_wind_u - tempbot*tau_bot_ave(1,Iq) + &
                            p(Iq,k) * grad_z(1,Iq,k) - p(Iq,k+1) * grad_z(1,Iq,k+1))
                source_y = gravity*(tau_wind_v - tempbot*tau_bot_ave(2,Iq) + &
                            p(Iq,k) * grad_z(2,Iq,k) - p(Iq,k+1) * grad_z(2,Iq,k+1))
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

    end subroutine create_rhs_dynamics_volume_layers_v2

    subroutine Apply_layers_fluxes(rhs_mom, qprime_face)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq, psiq, ngl
        use mod_grid, only:  npoin, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_input, only: nlayers, dt, dt_btp
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe, coeff_mass_pbpert_LR

        use mod_variables, only: udp_flux_edge, vdp_flux_edge, H_r_face, u_edge, v_edge, uvb_face_ave, ope_face_ave

        implicit none

        real, dimension(2, npoin, nlayers), intent(inout) :: rhs_mom
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face 
    
        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, nq_i, nq_j
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j, n, m
        real :: wq, hi, nxl, nyl, uu, vv
        real :: hlx_k, hly_k, hrx_k, hry_k, flux_x, flux_y, hx_k, hy_k
        real, dimension(nq) :: h_l, h_r
        real, dimension(2,nq) :: udp_flux, vdp_flux
        real :: utemp, vtemp, dptemp, udp_left_temp, udp_right_temp, vdp_left_temp, vdp_right_temp
    
        do k = 1, nlayers

            do iface = 1, nface

                !Store Left Side Variables

                el = face(7,iface)
                er = face(8,iface)

                h_l = H_r_face(1,:,iface,k)
                h_r = H_r_face(2,:,iface,k)

                udp_flux(1,:) = udp_flux_edge(1,:,iface,k)
                udp_flux(2,:) = udp_flux_edge(2,:,iface,k)

                vdp_flux(1,:) = vdp_flux_edge(1,:,iface,k)
                vdp_flux(2,:) = vdp_flux_edge(2,:,iface,k)

                do iquad = 1, nq

                    wq = jac_faceq(iquad,1,iface)

                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    hlx_k = nxl*h_l(iquad)
                    hrx_k = nxl*h_r(iquad)

                    hly_k = nyl*h_l(iquad)  
                    hry_k = nyl*h_r(iquad)
        
                    flux_x = nxl*udp_flux(1,iquad) + nyl*udp_flux(2,iquad)
                    flux_y = nxl*vdp_flux(1,iquad) + nyl*vdp_flux(2,iquad) 

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
        end do
    
    end subroutine Apply_layers_fluxes

    subroutine Apply_layers_fluxes_v2(rhs_mom, qprime_face)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq, psiq, ngl
        use mod_grid, only:  npoin, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_input, only: nlayers, dt, dt_btp
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe, coeff_mass_pbpert_LR

        use mod_variables, only: udp_flux_edge, vdp_flux_edge, H_r_face, u_edge, v_edge, uvb_face_ave, ope_face_ave

        implicit none

        real, dimension(2, npoin, nlayers), intent(inout) :: rhs_mom
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face

        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, nq_i, nq_j
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j, n, m
        real :: wq, hi, nxl, nyl, uu, vv
        real :: hlx_k, hly_k, hrx_k, hry_k, flux_x, flux_y, hx_k, hy_k
        real, dimension(nq) :: h_l, h_r
        real, dimension(2,nq) :: udp_flux, vdp_flux
        real :: utemp, vtemp, dptemp, udp_left_temp, udp_right_temp, vdp_left_temp, vdp_right_temp

        do k = 1, nlayers

            do iface = 1, nface

                !Store Left Side Variables

                el = face(7,iface)
                er = face(8,iface)

                do iquad = 1, nq

                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    utemp = qprime_face(2,1,iquad,iface,k) + uvb_face_ave(1,1,iquad,iface)
                    vtemp = qprime_face(3,1,iquad,iface,k) + uvb_face_ave(2,1,iquad,iface)
                    dptemp = qprime_face(1,1,iquad,iface,k) * ope_face_ave(1,iquad,iface)

                    udp_left_temp = utemp * dptemp
                    vdp_left_temp = vtemp * dptemp

                    ! Right side of the edge
                    utemp = qprime_face(2,2,iquad,iface,k) + uvb_face_ave(1,2,iquad,iface)
                    vtemp = qprime_face(3,2,iquad,iface,k) + uvb_face_ave(2,2,iquad,iface)
                    dptemp = qprime_face(1,2,iquad,iface,k) * ope_face_ave(2,iquad,iface)

                    udp_right_temp = utemp * dptemp
                    vdp_right_temp = vtemp * dptemp

                    uu = 0.5*(u_edge(1,iquad,iface,k) + u_edge(2,iquad,iface,k))
                    vv = 0.5*(v_edge(1,iquad,iface,k) + v_edge(2,iquad,iface,k))

                    if(uu*nxl > 0.0) then
                        udp_flux(1,iquad) = uu * udp_left_temp
                        vdp_flux(1,iquad) = uu * vdp_left_temp
                    else
                        udp_flux(1,iquad) = uu * udp_right_temp
                        vdp_flux(1,iquad) = uu * vdp_right_temp
                    endif

                    if(vv*nyl > 0.0) then
                        udp_flux(2,iquad) = vv * udp_left_temp
                        vdp_flux(2,iquad) = vv * vdp_left_temp
                    else
                        udp_flux(2,iquad) = vv * udp_right_temp
                        vdp_flux(2,iquad) = vv * vdp_right_temp
                    endif
                end do

                h_l = H_r_face(1,:,iface,k)
                h_r = H_r_face(2,:,iface,k)

                do iquad = 1, nq

                    wq = jac_faceq(iquad,1,iface)

                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    hlx_k = nxl*h_l(iquad)
                    hrx_k = nxl*h_r(iquad)

                    hly_k = nyl*h_l(iquad)
                    hry_k = nyl*h_r(iquad)

                    flux_x = nxl*udp_flux(1,iquad) + nyl*udp_flux(2,iquad)
                    flux_y = nxl*vdp_flux(1,iquad) + nyl*vdp_flux(2,iquad)

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
        end do

    end subroutine Apply_layers_fluxes_v2

    subroutine create_layers_volume_mass(dp_advec, qprime)

        use mod_grid, only : npoin_q, npoin, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac
        use mod_variables, only: ope_ave, uvb_ave, sum_layer_mass_flux

        implicit none 

        real, dimension(3, npoin_q, nlayers), intent(in) :: qprime
        real, dimension(npoin,nlayers), intent(out) :: dp_advec


        real :: wq, hi, dhde, dhdn

        integer :: k, I, Iq, ip
        real :: dp_temp, udp, vdp

        dp_advec = 0.0

        sum_layer_mass_flux = 0.0

        do k=1,nlayers

            do Iq = 1, npoin_q

                dp_temp = qprime(1,Iq,k) * ope_ave(Iq)
                udp = (qprime(2,Iq,k) + uvb_ave(1,Iq)) * dp_temp
                vdp = (qprime(3,Iq,k) + uvb_ave(2,Iq)) * dp_temp

                sum_layer_mass_flux(1,Iq) = sum_layer_mass_flux(1,Iq) + udp
                sum_layer_mass_flux(2,Iq) = sum_layer_mass_flux(2,Iq) + vdp
                    
                wq = wjac(Iq)

                do ip = 1, npts

                    I = indexq(ip,Iq)

                    dp_advec(I,k) = dp_advec(I,k) + wq*(dpsidx(ip,Iq)*udp + dpsidy(ip,Iq)*vdp)

                end do
            end do

        end do !k

    end subroutine create_layers_volume_mass

    subroutine create_consistency_volume_mass(dp_advec, qprime)

        use mod_grid, only : npoin_q, npoin, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, pbprime
        use mod_variables, only: sum_layer_mass_flux, btp_mass_flux_ave

        implicit none 

        real, dimension(3, npoin_q, nlayers), intent(in) :: qprime
        real, dimension(npoin,nlayers), intent(out) :: dp_advec


        real :: wq, hi, dhde, dhdn

        integer :: k, I, Iq, ip
        real :: weight, udp, vdp

        dp_advec = 0.0

        do k=1,nlayers

            do Iq = 1, npoin_q

                weight = qprime(1,Iq,k) / pbprime(Iq)

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

    subroutine create_layer_mass_flux(dp_advec, qprime_face)

        use mod_basis, only: nq, psiq, ngl
        use mod_grid, only:  npoin_q, intma,  nface, face
        use mod_input, only: nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_variables, only: uvb_face_ave, ope_face_ave, sum_layer_mass_flux_face, u_edge, v_edge

        implicit none

        real, dimension(npoin, nlayers), intent(inout) :: dp_advec
        real, dimension(3, 2, nq, nface, nlayers), intent(in) :: qprime_face
    
        integer :: k, iface, iquad, el, er, il, jl, ir, jr, I
        integer :: kl, kr, jquad, n, m
        real :: wq, nxl, nyl, hi
        real, dimension(nq) :: ul,ur,vl,vr, flux_edge_u, flux_edge_v
        real :: dpl, dpr, uu, vv, flux

        sum_layer_mass_flux_face = 0.0
    
        do k = 1, nlayers
            do iface = 1, nface

                do iquad = 1,nq

                    ! In the following computation of fluxes at element faces,
                    ! flux_edge iface a numerical approximation to the mass flux at element faces (each face has nq quadrature points)
                    ! Here we are using centered fluxes, so we need to compute the fluxes at the left and right edges of each face

                    ul(iquad) = qprime_face(2,1,iquad,iface,k) + uvb_face_ave(1,1,iquad,iface)
                    ur(iquad) = qprime_face(2,2,iquad,iface,k) + uvb_face_ave(1,2,iquad,iface)
                    vl(iquad) = qprime_face(3,1,iquad,iface,k) + uvb_face_ave(2,1,iquad,iface)
                    vr(iquad) = qprime_face(3,2,iquad,iface,k) + uvb_face_ave(2,2,iquad,iface)

                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    uu = 0.5*(ul(iquad) + ur(iquad))
                    vv = 0.5*(vl(iquad) + vr(iquad))
                    
                    dpl = ope_face_ave(1,iquad,iface) * qprime_face(1,1,iquad,iface,k)
                    dpr = ope_face_ave(2,iquad,iface) * qprime_face(1,2,iquad,iface,k)

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

                u_edge(1,:,iface,k) = ul
                u_edge(2,:,iface,k) = ur
                v_edge(1,:,iface,k) = vl
                v_edge(2,:,iface,k) = vr

                sum_layer_mass_flux_face(1,:,iface) = sum_layer_mass_flux_face(1,:,iface) + flux_edge_u(:)
                sum_layer_mass_flux_face(2,:,iface) = sum_layer_mass_flux_face(2,:,iface) + flux_edge_v(:)

                !Store Left Side Variables
                el = face(7,iface)
                er = face(8,iface)

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
