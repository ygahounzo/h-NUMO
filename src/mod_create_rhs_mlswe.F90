
module mod_create_rhs_mlswe


    use mod_constants, only: gravity

    use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad

    use mod_basis, only: nglx, ngly, nglz, npts, dpsiqx, dpsiqy, dpsiqz, nqx, nqy, nqz, &
        psiqx, psiqy, psiqz

    use mod_grid, only: intma, coord

    use mod_initial, only: q_ref, coriolis_constant, kvector, nvar

    use mod_barotropic_terms, only: massinv_rhs

    use mod_input, only: nlayers, mass_exact


    use mod_metrics, only: &
        ksiq_x, ksiq_y, ksiq_z, &
        etaq_x, etaq_y, etaq_z, &
        zetaq_x, zetaq_y, zetaq_z, &
        jacq, massinv

    public :: layer_momentum_rhs, layer_mass_advection_rhs, create_rhs_btp_momentum, btp_mass_advection_rhs, &
              interpolate_layer_from_quad_to_node, rhs_layer_shear_stress, interpolate_layer_from_quad_to_node_1d, &
              create_rhs_btp_momentum_new, create_rhs_btp_momentum_new1, create_rhs_btp_momentum1, create_rhs_btp_momentum_new2, &
              layer_mass_advection_rhs1, create_rhs_btp_dynamics_volume_new1, Apply_btp_fluxes_new, Apply_btp_fluxes_new_v4
              
contains


    subroutine layer_momentum_rhs(rhs_mom, H_r,p,grad_z,u_udp_temp,u_vdp_temp,v_vdp_temp,udp_flux_edge,vdp_flux_edge,H_r_face,tau_wind_int,tau_bot_int,rhs_visc,q_face, disp)

        use mod_input, only: nlayers
        use mod_metrics, only: massinv
        use mod_grid, only: npoin, npoin_q, intma, intma_dg_quad, nface
        use mod_basis, only: nq

        implicit none

        real, dimension(npoin_q, nlayers), intent(in) :: H_r
        real, dimension(npoin_q, nlayers+1), intent(in) :: p
        real, dimension(2, npoin_q, nlayers+1), intent(in) :: grad_z,tau_wind_int, tau_bot_int
        real, dimension(npoin_q, nlayers), intent(in) :: u_udp_temp, v_vdp_temp
        real, dimension(2, npoin_q, nlayers), intent(in) :: u_vdp_temp
        real, dimension(2, nq, nface, nlayers), intent(in) :: udp_flux_edge, vdp_flux_edge, H_r_face
        real, dimension(2,npoin,nlayers), intent(in) :: rhs_visc
        real, dimension(3, 2, nq, nface, nlayers), intent(in) :: q_face
        real, dimension(nq, nface, nlayers), intent(in) :: disp

        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom

        real, dimension(2,npoin_q) :: qb2

        integer :: k,I

        call create_rhs_dynamics_volume_layers(rhs_mom, H_r, p, grad_z, u_udp_temp, u_vdp_temp, v_vdp_temp, tau_wind_int, tau_bot_int)

        call Apply_layers_fluxes(rhs_mom, udp_flux_edge, vdp_flux_edge, H_r_face, q_face, disp)

        do k = 1, nlayers

            rhs_mom(1,:,k) = massinv(:)*rhs_mom(1,:,k)
            rhs_mom(2,:,k) = massinv(:)*rhs_mom(2,:,k)
        end do

    end subroutine layer_momentum_rhs

    subroutine layer_mass_advection_rhs(dp_advec, uvdp_temp, flux_edge)

        use mod_metrics, only: massinv
        use mod_input, only: nlayers
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq

        implicit none

        real, dimension(npoin, nlayers), intent(out) :: dp_advec
        real, dimension(2, npoin_q, nlayers), intent(in) :: uvdp_temp
        real, dimension(2, nq, nface, nlayers), intent(in) :: flux_edge
        
        integer :: k,I

        ! Compute the mass advection term for the degree of freedom for dp in each layer

        call create_layers_volume_mass(dp_advec, uvdp_temp)

        ! Compute the mass flux term 

        call create_layer_mass_flux(dp_advec, flux_edge)

        do k = 1, nlayers
            dp_advec(:,k) = massinv(:)*dp_advec(:,k)
        end do
        
    end subroutine layer_mass_advection_rhs

    subroutine layer_mass_advection_rhs1(dp_advec, uvdp_temp, flux_edge, q_face, disp)

        use mod_metrics, only: massinv
        use mod_input, only: nlayers
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq

        implicit none

        real, dimension(npoin, nlayers), intent(out) :: dp_advec
        real, dimension(2, npoin_q, nlayers), intent(in) :: uvdp_temp
        real, dimension(2, nq, nface, nlayers), intent(in) :: flux_edge
        real, dimension(3, 2, nq, nface, nlayers), intent(in) :: q_face
        real, dimension(nq, nface, nlayers), intent(in) :: disp
        
        integer :: k,I

        ! Compute the mass advection term for the degree of freedom for dp in each layer

        call create_layers_volume_mass(dp_advec, uvdp_temp)

        ! Compute the mass flux term 

        call create_layer_mass_flux1(dp_advec, flux_edge, q_face, disp)

        do k = 1, nlayers

            dp_advec(:,k) = massinv(:)*dp_advec(:,k)
        end do
        
    end subroutine layer_mass_advection_rhs1

    subroutine create_rhs_btp_momentum(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc,qb_face)


        use mod_metrics, only: massinv
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq

        implicit none

        real, dimension(2, npoin), intent(out) :: rhs_mom
        real, dimension(npoin_q), intent(in) :: H, Quu, Qvv, Quv
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(2, nq, nface), intent(in) :: Qu_face, Qv_face
        real, dimension(nq, nface), intent(in) :: H_face
        real, dimension(2, npoin_q), intent(in) :: tau_bot
        real, dimension(2, npoin), intent(in) :: rhs_visc
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face

        real, dimension(2,npoin_q) :: qb2
        

        integer :: I

        call create_rhs_btp_dynamics_volume(rhs_mom, H, qb, Quu, Quv, Qvv,tau_bot)

        call Apply_btp_fluxes(rhs_mom,H_face,Qu_face,Qv_face, qb_face)

    end subroutine create_rhs_btp_momentum

    subroutine create_rhs_btp_momentum1(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc,qb_face)


        use mod_metrics, only: massinv
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq

        implicit none

        real, dimension(2, npoin), intent(out) :: rhs_mom
        real, dimension(npoin_q), intent(in) :: H, Quu, Qvv, Quv
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(2, nq, nface), intent(in) :: Qu_face, Qv_face
        real, dimension(nq, nface), intent(in) :: H_face
        real, dimension(2, npoin_q), intent(in) :: tau_bot
        real, dimension(2, npoin), intent(in) :: rhs_visc
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face

        ! real, dimension(2,npoin_q) :: qb2
        

        integer :: I

        call create_rhs_btp_dynamics_volume1(rhs_mom, H, qb, Quu, Quv, Qvv,tau_bot)

        call Apply_btp_fluxes(rhs_mom,H_face,Qu_face,Qv_face, qb_face)

    end subroutine create_rhs_btp_momentum1

    subroutine create_rhs_btp_momentum_new(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc,btp_mass_flux, flux_edge,qb_face)


        use mod_metrics, only: massinv
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq

        implicit none

        real, dimension(3, npoin), intent(out) :: rhs_mom
        real, dimension(npoin_q), intent(in) :: H, Quu, Qvv, Quv
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(2, nq, nface), intent(in) :: Qu_face, Qv_face, flux_edge
        real, dimension(nq, nface), intent(in) :: H_face
        real, dimension(2, npoin_q), intent(in) :: tau_bot, btp_mass_flux
        real, dimension(2, npoin), intent(in) :: rhs_visc
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face


        call create_rhs_btp_dynamics_volume_new(rhs_mom, H, qb, Quu, Quv, Qvv,tau_bot, btp_mass_flux)

        call Apply_btp_fluxes_new(rhs_mom,H_face,Qu_face,Qv_face, qb_face, flux_edge)

    end subroutine create_rhs_btp_momentum_new

    subroutine create_rhs_btp_momentum_new1(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,rhs_visc,btp_mass_flux, flux_edge,qb_face)


        use mod_metrics, only: massinv
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq

        implicit none

        real, dimension(3, npoin), intent(out) :: rhs_mom
        real, dimension(npoin_q), intent(in) :: H, Quu, Qvv, Quv
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(2, nq, nface), intent(in) :: Qu_face, Qv_face, flux_edge
        real, dimension(nq, nface), intent(in) :: H_face
        real, dimension(2, npoin_q), intent(in) :: tau_bot, btp_mass_flux
        real, dimension(2, npoin), intent(in) :: rhs_visc
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face


        call create_rhs_btp_dynamics_volume_new1(rhs_mom, H, qb, Quu, Quv, Qvv,tau_bot, btp_mass_flux)

        call Apply_btp_fluxes_new(rhs_mom,H_face,Qu_face,Qv_face, qb_face, flux_edge)

    end subroutine create_rhs_btp_momentum_new1

    subroutine create_rhs_btp_momentum_new2(rhs_mom,Quu,Qvv,Quv,H,qb,H_face,Qu_face,Qv_face,tau_bot,btp_mass_flux, flux_edge, grad_uvdp,grad_uvdp_face, qb_face)


        use mod_metrics, only: massinv
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq

        implicit none

        real, dimension(3, npoin), intent(out) :: rhs_mom
        real, dimension(npoin_q), intent(in) :: H, Quu, Qvv, Quv
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(2, nq, nface), intent(in) :: Qu_face, Qv_face, flux_edge
        real, dimension(nq, nface), intent(in) :: H_face
        real, dimension(2, npoin_q), intent(in) :: tau_bot, btp_mass_flux
        real, dimension(4,2,nq,nface), intent(in) :: grad_uvdp_face
        real, dimension(2,2,npoin_q), intent(in) :: grad_uvdp
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face

        call create_rhs_btp_dynamics_volume_new_v2(rhs_mom, H, qb, Quu, Quv, Qvv,tau_bot, btp_mass_flux, grad_uvdp)

        call Apply_btp_fluxes_new_v2(rhs_mom,H_face,Qu_face,Qv_face, flux_edge, grad_uvdp_face, qb_face)

    end subroutine create_rhs_btp_momentum_new2

    subroutine btp_mass_advection_rhs(pb_advec, btp_mass_flux, flux_edge)

        use mod_metrics, only: massinv 
        use mod_grid, only: npoin_q, npoin, nface
        use mod_basis, only: nq

        implicit none

        real, dimension(2, npoin_q), intent(in) :: btp_mass_flux
        real, dimension(2, nq, nface), intent(in) :: flux_edge
        real, dimension(npoin), intent(out) :: pb_advec

        call create_rhs_btp_dynamics_volume_mass(pb_advec, btp_mass_flux)
        call  Apply_btp_flux_mass(pb_advec, flux_edge)

    end subroutine btp_mass_advection_rhs


    subroutine interpolate_layer_from_quad_to_node(q_df,q)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_metrics, only: massinv
        use mod_input, only: nlayers
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac

        implicit none
        
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(5,npoin_q,nlayers), intent(in) :: q

        real :: wq, var_u(nlayers), var_v(nlayers), hi
        integer :: iquad, jquad, kquad, e, k, m, n, l, I, Iq, ip

        q_df(2:3,:,:) = 0.0
        
        ! do k = 1, nlayers

        do Iq = 1,npoin_q
            
            wq = wjac(Iq)
            var_u = q(4,Iq,:)
            var_v = q(5,Iq,:)

            do ip = 1, npts

                I = indexq(Iq,ip)

                hi = psih(Iq,ip)

                q_df(2,I,:) = q_df(2,I,:) + wq*var_u(:)*hi
                q_df(3,I,:) = q_df(3,I,:) + wq*var_v(:)*hi

            end do
        end do

        do k = 1, nlayers

            if(mass_exact) then 
                call massinv_rhs(q_df(2:3,:,k),2) 
            else
                q_df(2,:,k) = massinv(:)*q_df(2,:,k)
                q_df(3,:,k) = massinv(:)*q_df(3,:,k)
            end if    

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

                I = indexq(Iq,ip)

                hi = psih(Iq,ip)

                rhs_stress(1,I,:) = rhs_stress(1,I,:) + wq*hi*tau_u_q(Iq,:)
                rhs_stress(2,I,:) = rhs_stress(2,I,:) + wq*hi*tau_v_q(Iq,:)

            end do

        end do

        do k = 1, nlayers
            if(mass_exact) then 
                call massinv_rhs(rhs_stress(:,:,k),2)
            else
                rhs_stress(1,:,k) = massinv(:)*rhs_stress(1,:,k)
                rhs_stress(2,:,k) = massinv(:)*rhs_stress(2,:,k)
            end if
        end do
        
    end subroutine rhs_layer_shear_stress

    subroutine create_rhs_dynamics_volume_layers(rhs_mom, H_r, p, grad_z, u_udp_temp, u_vdp_temp, v_vdp_temp, tau_wind_int, tau_bot_int)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts, dpsiqx
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac

        implicit none 

        real, dimension(npoin_q, nlayers), intent(in) :: H_r
        real, dimension(npoin_q, nlayers+1), intent(in) :: p
        real, dimension(2, npoin_q, nlayers+1), intent(in) :: grad_z,tau_wind_int, tau_bot_int
        real, dimension(npoin_q, nlayers), intent(in) :: u_udp_temp, v_vdp_temp
        real, dimension(2, npoin_q, nlayers), intent(in) :: u_vdp_temp

        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom

        real, dimension(npoin_q,nlayers) :: tau_u_wind, tau_u_bot, tau_v_wind, tau_v_bot, pz_x, pz_y


        real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, hi, &
                h_e, h_n, dhdx, dhdy, h_c
        real, dimension(nlayers) :: u_mom, v_mom, Hq, var_uu, var_uv, var_vu, var_vv, u_cori, v_cori

        integer :: k, e, iquad, jquad, kquad, l, n, m, I, Iq, ip

        rhs_mom = 0.0

        do k=1,nlayers

            tau_u_wind(:,k) = gravity*tau_wind_int(1,:,k)
            tau_u_bot(:,k) = -gravity*tau_bot_int(1,:,k)

            tau_v_wind(:,k) = gravity*tau_wind_int(2,:,k)
            tau_v_bot(:,k) = -gravity*tau_bot_int(2,:,k)

            pz_x(:,k) = gravity*(p(:,k) * grad_z(1,:,k) - p(:,k+1) * grad_z(1,:,k+1))
            pz_y(:,k) = gravity*(p(:,k) * grad_z(2,:,k) - p(:,k+1) * grad_z(2,:,k+1))
        end do

        do Iq = 1, npoin_q

            wq = wjac(Iq)

            Hq = H_r(Iq,:)

            var_uu = u_udp_temp(Iq,:) 
            var_uv = u_vdp_temp(1,Iq,:)
            var_vu = u_vdp_temp(2,Iq,:)
            var_vv = v_vdp_temp(Iq,:)

            do ip = 1, npts

                I = indexq(Iq,ip)

                hi = psih(Iq,ip)

                rhs_mom(1,I,:) = rhs_mom(1,I,:) + wq*hi*(tau_u_wind(Iq,:) + tau_u_bot(Iq,:) + pz_x(Iq,:))
                rhs_mom(2,I,:) = rhs_mom(2,I,:) + wq*hi*(tau_v_wind(Iq,:) + tau_v_bot(Iq,:) + pz_y(Iq,:))

                !Xi derivatives
                dhdx = dpsidx(Iq,ip)
                !Eta derivatives
                dhdy = dpsidy(Iq,ip)

                rhs_mom(1,I,:) = rhs_mom(1,I,:) + wq*dhdx*Hq(:) 
                rhs_mom(2,I,:) = rhs_mom(2,I,:) + wq*dhdy*Hq(:)

                u_mom(:) = var_uu(:)*dhdx + var_uv(:)*dhdy
                v_mom(:) = var_vu(:)*dhdx + var_vv(:)*dhdy

                rhs_mom(1,I,:) = rhs_mom(1,I,:) + wq*u_mom(:)
                rhs_mom(2,I,:) = rhs_mom(2,I,:) + wq*v_mom(:)

            end do
        end do

        ! end do !k

    end subroutine create_rhs_dynamics_volume_layers


    subroutine Apply_layers_fluxes(rhs_mom,udp_flux_edge, vdp_flux_edge, H_r_face, q_face, disp)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq, psiq, ngl
        use mod_grid, only:  npoin, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_input, only: nlayers, dt, dt_btp
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe, coeff_mass_pbpert_LR

        implicit none

        real, dimension(2, npoin, nlayers), intent(inout) :: rhs_mom
        real, dimension(2, nq, nface, nlayers), intent(in) :: udp_flux_edge, vdp_flux_edge, H_r_face
        real, dimension(3, 2, nq, nface, nlayers), intent(in) :: q_face
        real, dimension(nq, nface, nlayers), intent(in) :: disp
    
        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, nq_i, nq_j
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j, n, m
        real :: wq, hi, nxl, nyl, gprime, h1, h2, clam_temp, lam
        real, dimension(nlayers) :: hlx_k, hly_k, hrx_k, hry_k, flux_x, flux_y, hx_k, hy_k
        real, dimension(nq,nlayers) :: h_l, h_r, var_uudp, var_uvdp, var_vudp, var_vvdp
        real, dimension(nlayers) :: unl, unr, Hl, Hr, dispu, dispv, ul, vl, ur, vr
        real :: h1_l, h1_r, h2_l, h2_r, claml1, clamr1, claml2, clamr2, clam1, clam2
        real :: Ucl, Ucr, hgl, hgr, claml, clamr, clam

        !gprime = gravity *alpha_mlswe(2)*(1.0/alpha_mlswe(2) - 1.0/alpha_mlswe(1))
    
        ! do k = 1, nlayers

        do iface = 1, nface

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            h_l = H_r_face(1,:,iface,:)
            h_r = H_r_face(2,:,iface,:)

            var_uudp = udp_flux_edge(1,:,iface,:)
            var_uvdp = udp_flux_edge(2,:,iface,:)

            var_vudp = vdp_flux_edge(1,:,iface,:)
            var_vvdp = vdp_flux_edge(2,:,iface,:)   

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)
                
                clam1 = disp(iquad,iface,1)
                clam2 = disp(iquad,iface,2)
                
                dispu(1) = clam1*(q_face(2,2,iquad,iface,1) - q_face(2,1,iquad,iface,1))
                dispu(2) = clam2*(q_face(2,2,iquad,iface,2) - q_face(2,1,iquad,iface,2))
                dispv(1) = clam1*(q_face(3,2,iquad,iface,1) - q_face(3,1,iquad,iface,1))
                dispv(2) = clam2*(q_face(3,2,iquad,iface,2) - q_face(3,1,iquad,iface,2))

                hlx_k = nxl*h_l(iquad,:)
                hrx_k = nxl*h_r(iquad,:)
                
                ! hx_k = 0.5*nxl*(h_l(iquad,:) + h_r(iquad,:))

                hly_k = nyl*h_l(iquad,:)  
                hry_k = nyl*h_r(iquad,:)

                ! hy_k = 0.5*nyl*(h_l(iquad,:) + h_r(iquad,:))
    
                flux_x = nxl*var_uudp(iquad,:) + nyl*var_uvdp(iquad,:) !- dispu(:)
                flux_y = nxl*var_vudp(iquad,:) + nyl*var_vvdp(iquad,:) !- dispv(:)

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    rhs_mom(1,I,:) = rhs_mom(1,I,:) - wq*hi*(hlx_k + flux_x)
                    rhs_mom(2,I,:) = rhs_mom(2,I,:) - wq*hi*(hly_k + flux_y)

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        rhs_mom(1,I,:) = rhs_mom(1,I,:) + wq*hi*(hrx_k + flux_x)
                        rhs_mom(2,I,:) = rhs_mom(2,I,:) + wq*hi*(hry_k + flux_y)

                    end if

                end do
            end do

        end do
        ! end do
    
    end subroutine Apply_layers_fluxes


    subroutine create_layers_volume_mass(dp_advec, uvdp_temp)

        use mod_grid, only : npoin_q, npoin, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_input, only: nlayers
        use mod_constants, only: gravity
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac


        implicit none 

        real, dimension(2, npoin_q, nlayers), intent(in) :: uvdp_temp
        real, dimension(npoin,nlayers), intent(out) :: dp_advec


        real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, hi, dhde, dhdn

        integer :: k, e, iquad, jquad, kquad, l, n, m, I, Iq, ip
        real, dimension(nlayers) :: var_u, var_v

        dp_advec = 0.0

        ! do k=1,nlayers

        do Iq = 1, npoin_q
                
            wq = wjac(Iq)

            var_u = uvdp_temp(1,Iq,:) 
            var_v = uvdp_temp(2,Iq,:)

            do ip = 1, npts

                I = indexq(Iq,ip)

                dp_advec(I,:) = dp_advec(I,:) + wq*(dpsidx(Iq,ip)*var_u + dpsidy(Iq,ip)*var_v)

            end do
        end do

        ! end do !k

    end subroutine create_layers_volume_mass

    subroutine create_layer_mass_flux(dp_advec, flux_edge)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq, psiq, ngl
        use mod_grid, only:  npoin_q, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_input, only: nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq

        implicit none

        real, dimension(npoin, nlayers), intent(inout) :: dp_advec
        real, dimension(2, nq, nface, nlayers), intent(in) :: flux_edge
    
        integer :: k, iface, iquad, el, er, il, jl, ir, jr, I
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j, n, m
        real :: wq, nxl, nyl, hi
        real, dimension(nq,nlayers) :: var_x, var_y
        real, dimension(nlayers) :: qx_k, qy_k, flux
    
        ! do k = 1, nlayers
        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            var_x = flux_edge(1,:,iface,:)
            var_y = flux_edge(2,:,iface,:)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                qx_k = var_x(iquad,:);  qy_k = var_y(iquad,:)

                flux = nxl*qx_k + nyl*qy_k 

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    dp_advec(I,:) = dp_advec(I,:) - wq*hi*flux

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        dp_advec(I,:) = dp_advec(I,:) + wq*hi*flux

                    end if

                end do
            end do

        end do
        ! end do
    
    end subroutine create_layer_mass_flux

    subroutine create_layer_mass_flux1(dp_advec, flux_edge, q_face, disp)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq, psiq, ngl
        use mod_grid, only:  npoin_q, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_input, only: nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe

        implicit none

        real, dimension(npoin, nlayers), intent(inout) :: dp_advec
        real, dimension(2, nq, nface, nlayers), intent(in) :: flux_edge
        real, dimension(3, 2, nq, nface, nlayers), intent(in) :: q_face
        real, dimension(nq, nface, nlayers), intent(in) :: disp
    
        integer :: k, iface, iquad, el, er, il, jl, ir, jr, I
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j, n, m
        real :: wq, nxl, nyl, hi
        real, dimension(nq,nlayers) :: var_x, var_y
        real, dimension(nlayers) :: qx_k, qy_k, flux, ul, vl, ur, vr, unl, unr, disp_temp
        real :: h1_l, h1_r, h2_l, h2_r, claml1, clamr1, claml2, clamr2, clam1, clam2
        real :: claml, clamr, gprime, clam, Ucl, Ucr, hgl, hgr

        !gprime = gravity *alpha_mlswe(2)*(1.0/alpha_mlswe(2) - 1.0/alpha_mlswe(1))
    
        ! do k = 1, nlayers
        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            var_x = flux_edge(1,:,iface,:)
            var_y = flux_edge(2,:,iface,:)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                clam1 = disp(iquad,iface,1)
                clam2 = disp(iquad,iface,2)
                
                disp_temp(1) = clam1*(q_face(1,2,iquad,iface,1) - q_face(1,1,iquad,iface,1))
                disp_temp(2) = clam2*(q_face(1,2,iquad,iface,2) - q_face(1,1,iquad,iface,2))

                qx_k = var_x(iquad,:);  qy_k = var_y(iquad,:)

                flux = nxl*qx_k + nyl*qy_k !- disp_temp(:)

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    dp_advec(I,:) = dp_advec(I,:) - wq*hi*flux

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        dp_advec(I,:) = dp_advec(I,:) + wq*hi*flux

                    end if

                end do
            end do

        end do
        ! end do
    
    end subroutine create_layer_mass_flux1


    subroutine create_rhs_btp_dynamics_volume(rhs_mom, H, qb, Quu, Quv, Qvv, tau_bot)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, coriolis_quad


        implicit none 

        real, dimension(npoin_q), intent(in) :: H, Quu, Qvv, Quv
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(2, npoin_q), intent(in) :: tau_bot

        real, dimension(2,npoin), intent(out) :: rhs_mom

        real, dimension(npoin_q) :: stress_u, stress_v, pb_x, pb_y

        real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, u_cori, v_cori, Hq, var_uu, var_uv, var_vu, var_vv, hi, &
                h_e, h_n, dhdx, dhdy, u_mom, v_mom, h_c

        integer :: k, e, iquad, jquad, kquad, l, n, m, I, Iq, ip

        rhs_mom = 0.0

        stress_u(:) = gravity*(tau_wind(1,:) - tau_bot(1,:))
        stress_v(:) = gravity*(tau_wind(2,:) - tau_bot(2,:))

        pb_x(:) = gravity*qb(1,:)*grad_zbot_quad(1,:)
        pb_y(:) = gravity*qb(1,:)*grad_zbot_quad(2,:)

        do Iq = 1, npoin_q

            wq = wjac(Iq)

            do ip = 1,npts

                I = indexq(Iq,ip)

                hi = psih(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*hi*(stress_u(Iq) - pb_x(Iq))
                rhs_mom(2,I) = rhs_mom(2,I) + wq*hi*(stress_v(Iq) - pb_y(Iq))

                !Xi derivatives
                dhdx = dpsidx(Iq,ip)
                !Eta derivatives
                dhdy = dpsidy(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*dhdx*H(Iq)
                rhs_mom(2,I) = rhs_mom(2,I) + wq*dhdy*H(Iq)

                u_mom = Quu(Iq)*dhdx + Quv(Iq)*dhdy
                v_mom = Quv(Iq)*dhdx + Qvv(Iq)*dhdy

                rhs_mom(1,I) = rhs_mom(1,I) + wq*u_mom
                rhs_mom(2,I) = rhs_mom(2,I) + wq*v_mom

            end do
        end do

    end subroutine create_rhs_btp_dynamics_volume

    subroutine create_rhs_btp_dynamics_volume1(rhs_mom, H, qb, Quu, Quv, Qvv, tau_bot)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, coriolis_quad


        implicit none 

        real, dimension(npoin_q), intent(in) :: H, Quu, Qvv, Quv
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(2, npoin_q), intent(in) :: tau_bot

        real, dimension(2,npoin), intent(out) :: rhs_mom

        real, dimension(npoin_q) :: stress_u, stress_v, pb_x, pb_y

        real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, u_cori, v_cori, Hq, var_uu, var_uv, var_vu, var_vv, hi, &
                h_e, h_n, dhdx, dhdy, u_mom, v_mom, h_c

        integer :: k, e, iquad, jquad, kquad, l, n, m, I, Iq, ip

        rhs_mom = 0.0

        stress_u(:) = gravity*(tau_wind(1,:) - tau_bot(1,:))
        stress_v(:) = gravity*(tau_wind(2,:) - tau_bot(2,:))

        pb_x(:) = gravity*qb(1,:)*grad_zbot_quad(1,:)
        pb_y(:) = gravity*qb(1,:)*grad_zbot_quad(2,:)

        do Iq = 1, npoin_q

            wq = wjac(Iq)
            u_cori = coriolis_quad(Iq)*qb(4,Iq)
            v_cori = -coriolis_quad(Iq)*qb(3,Iq)

            do ip = 1,npts

                I = indexq(Iq,ip)

                hi = psih(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*hi*(stress_u(Iq) - pb_x(Iq) + u_cori)
                rhs_mom(2,I) = rhs_mom(2,I) + wq*hi*(stress_v(Iq) - pb_y(Iq) + v_cori)

                !Xi derivatives
                dhdx = dpsidx(Iq,ip)
                !Eta derivatives
                dhdy = dpsidy(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*dhdx*H(Iq)
                rhs_mom(2,I) = rhs_mom(2,I) + wq*dhdy*H(Iq)

                u_mom = Quu(Iq)*dhdx + Quv(Iq)*dhdy
                v_mom = Quv(Iq)*dhdx + Qvv(Iq)*dhdy

                rhs_mom(1,I) = rhs_mom(1,I) + wq*u_mom
                rhs_mom(2,I) = rhs_mom(2,I) + wq*v_mom

            end do
        end do

    end subroutine create_rhs_btp_dynamics_volume1

    subroutine create_rhs_btp_dynamics_volume_new(rhs_mom, H, qb, Quu, Quv, Qvv, tau_bot, btp_mass_flux)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, coriolis_quad


        implicit none 

        real, dimension(npoin_q), intent(in) :: H, Quu, Qvv, Quv
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(2, npoin_q), intent(in) :: tau_bot, btp_mass_flux

        real, dimension(3,npoin), intent(out) :: rhs_mom

        real, dimension(npoin_q) :: stress_u, stress_v, pb_x, pb_y

        real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, u_cori, v_cori, Hq, var_uu, var_uv, var_vu, var_vv, hi, &
                h_e, h_n, dhdx, dhdy, u_mom, v_mom, h_c, var_u, var_v

        integer :: k, e, iquad, jquad, kquad, l, n, m, I, Iq, ip

        rhs_mom = 0.0

        stress_u(:) = gravity*(tau_wind(1,:) - tau_bot(1,:))
        stress_v(:) = gravity*(tau_wind(2,:) - tau_bot(2,:))

        pb_x(:) = gravity*qb(1,:)*grad_zbot_quad(1,:)
        pb_y(:) = gravity*qb(1,:)*grad_zbot_quad(2,:)

        do Iq = 1, npoin_q

            wq = wjac(Iq)
            var_u = btp_mass_flux(1,Iq)
            var_v = btp_mass_flux(2,Iq)

            do ip = 1,npts

                I = indexq(Iq,ip)

                hi = psih(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*hi*(stress_u(Iq) - pb_x(Iq))
                rhs_mom(2,I) = rhs_mom(2,I) + wq*hi*(stress_v(Iq) - pb_y(Iq))

                !Xi derivatives
                dhdx = dpsidx(Iq,ip)
                !Eta derivatives
                dhdy = dpsidy(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*dhdx*H(Iq)
                rhs_mom(2,I) = rhs_mom(2,I) + wq*dhdy*H(Iq)

                u_mom = Quu(Iq)*dhdx + Quv(Iq)*dhdy
                v_mom = Quv(Iq)*dhdx + Qvv(Iq)*dhdy

                rhs_mom(1,I) = rhs_mom(1,I) + wq*u_mom
                rhs_mom(2,I) = rhs_mom(2,I) + wq*v_mom

                rhs_mom(3,I) = rhs_mom(3,I) + wq*(dhdx*var_u + dhdy*var_v)

            end do
        end do

    end subroutine create_rhs_btp_dynamics_volume_new

    subroutine create_rhs_btp_dynamics_volume_new1(rhs_mom, H, qb, Quu, Quv, Qvv, tau_bot, btp_mass_flux)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, coriolis_quad


        implicit none 

        real, dimension(npoin_q), intent(in) :: H, Quu, Qvv, Quv
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(2, npoin_q), intent(in) :: tau_bot, btp_mass_flux

        real, dimension(3,npoin), intent(out) :: rhs_mom

        real, dimension(npoin_q) :: stress_u, stress_v, pb_x, pb_y

        real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, u_cori, v_cori, Hq, var_uu, var_uv, var_vu, var_vv, hi, &
                h_e, h_n, dhdx, dhdy, u_mom, v_mom, h_c, var_u, var_v

        integer :: k, e, iquad, jquad, kquad, l, n, m, I, Iq, ip

        rhs_mom = 0.0

        stress_u(:) = gravity*(tau_wind(1,:) - tau_bot(1,:))
        stress_v(:) = gravity*(tau_wind(2,:) - tau_bot(2,:))

        pb_x(:) = gravity*qb(1,:)*grad_zbot_quad(1,:)
        pb_y(:) = gravity*qb(1,:)*grad_zbot_quad(2,:)

        do Iq = 1, npoin_q

            wq = wjac(Iq)
            var_u = btp_mass_flux(1,Iq)
            var_v = btp_mass_flux(2,Iq)

            u_cori = coriolis_quad(Iq)*qb(4,Iq)
            v_cori = -coriolis_quad(Iq)*qb(3,Iq)

            do ip = 1,npts

                I = indexq(Iq,ip)

                hi = psih(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*hi*(stress_u(Iq) - pb_x(Iq) + u_cori)
                rhs_mom(2,I) = rhs_mom(2,I) + wq*hi*(stress_v(Iq) - pb_y(Iq) + v_cori)

                !Xi derivatives
                dhdx = dpsidx(Iq,ip)
                !Eta derivatives
                dhdy = dpsidy(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*dhdx*H(Iq)
                rhs_mom(2,I) = rhs_mom(2,I) + wq*dhdy*H(Iq)

                u_mom = Quu(Iq)*dhdx + Quv(Iq)*dhdy
                v_mom = Quv(Iq)*dhdx + Qvv(Iq)*dhdy

                rhs_mom(1,I) = rhs_mom(1,I) + wq*u_mom
                rhs_mom(2,I) = rhs_mom(2,I) + wq*v_mom

                rhs_mom(3,I) = rhs_mom(3,I) + wq*(dhdx*var_u + dhdy*var_v)

            end do
        end do

    end subroutine create_rhs_btp_dynamics_volume_new1

    subroutine create_rhs_btp_dynamics_volume_new_v2(rhs_mom, H, qb, Quu, Quv, Qvv, tau_bot, btp_mass_flux, grad_uvdp)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, coriolis_quad
        use mod_input, only: visc_mlswe


        implicit none 

        real, dimension(npoin_q), intent(in) :: H, Quu, Qvv, Quv
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(2, npoin_q), intent(in) :: tau_bot, btp_mass_flux
        real, dimension(2,2,npoin_q), intent(in) :: grad_uvdp

        real, dimension(3,npoin), intent(out) :: rhs_mom

        real, dimension(npoin_q) :: stress_u, stress_v, pb_x, pb_y

        real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, u_cori, v_cori, Hq, var_uu, var_uv, var_vu, var_vv, hi, &
                h_e, h_n, dhdx, dhdy, u_mom, v_mom, h_c, var_u, var_v

        integer :: k, e, iquad, jquad, kquad, l, n, m, I, Iq, ip

        rhs_mom = 0.0

        stress_u(:) = gravity*(tau_wind(1,:) - tau_bot(1,:))
        stress_v(:) = gravity*(tau_wind(2,:) - tau_bot(2,:))

        pb_x(:) = gravity*qb(1,:)*grad_zbot_quad(1,:)
        pb_y(:) = gravity*qb(1,:)*grad_zbot_quad(2,:)

        do Iq = 1, npoin_q

            wq = wjac(Iq)
            var_u = btp_mass_flux(1,Iq)
            var_v = btp_mass_flux(2,Iq)

            u_cori = coriolis_quad(Iq)*qb(4,Iq)
            v_cori = -coriolis_quad(Iq)*qb(3,Iq)

            do ip = 1,npts

                I = indexq(Iq,ip)

                hi = psih(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*hi*(stress_u(Iq) - pb_x(Iq) + u_cori)
                rhs_mom(2,I) = rhs_mom(2,I) + wq*hi*(stress_v(Iq) - pb_y(Iq) + v_cori)

                !Xi derivatives
                dhdx = dpsidx(Iq,ip)
                !Eta derivatives
                dhdy = dpsidy(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*dhdx*H(Iq)
                rhs_mom(2,I) = rhs_mom(2,I) + wq*dhdy*H(Iq)

                u_mom = Quu(Iq)*dhdx + Quv(Iq)*dhdy
                v_mom = Quv(Iq)*dhdx + Qvv(Iq)*dhdy

                rhs_mom(1,I) = rhs_mom(1,I) + wq*u_mom
                rhs_mom(2,I) = rhs_mom(2,I) + wq*v_mom

                u_mom = dhdx*grad_uvdp(1,1,Iq) + dhdy*grad_uvdp(1,2,Iq)
                v_mom = dhdx*grad_uvdp(2,1,Iq) + dhdy*grad_uvdp(2,2,Iq)

                rhs_mom(1,I) = rhs_mom(1,I) - wq*visc_mlswe*u_mom
                rhs_mom(2,I) = rhs_mom(2,I) - wq*visc_mlswe*v_mom

                rhs_mom(3,I) = rhs_mom(3,I) + wq*(dhdx*var_u + dhdy*var_v)

            end do
        end do

    end subroutine create_rhs_btp_dynamics_volume_new_v2

    subroutine create_rhs_btp_dynamics_volume_new2(rhs_mom, H, qb, Quu, Quv, Qvv, tau_bot, btp_mass_flux,qb2)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts
        use mod_constants, only: gravity
        use mod_initial, only: grad_zbot_quad, tau_wind, psih, dpsidx,dpsidy, indexq, wjac, coriolis_quad


        implicit none 

        real, dimension(npoin_q), intent(in) :: H, Quu, Qvv, Quv
        real, dimension(4,npoin_q), intent(in) :: qb,qb2
        real, dimension(2, npoin_q), intent(in) :: tau_bot, btp_mass_flux

        real, dimension(3,npoin), intent(out) :: rhs_mom

        real, dimension(npoin_q) :: stress_u, stress_v, pb_x, pb_y

        real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, u_cori, v_cori, Hq, var_uu, var_uv, var_vu, var_vv, hi, &
                h_e, h_n, dhdx, dhdy, u_mom, v_mom, h_c, var_u, var_v

        integer :: k, e, iquad, jquad, kquad, l, n, m, I, Iq, ip

        rhs_mom = 0.0

        stress_u(:) = gravity*(tau_wind(1,:) - tau_bot(1,:))
        stress_v(:) = gravity*(tau_wind(2,:) - tau_bot(2,:))

        pb_x(:) = gravity*qb(1,:)*grad_zbot_quad(1,:)
        pb_y(:) = gravity*qb(1,:)*grad_zbot_quad(2,:)

        do Iq = 1, npoin_q

            wq = wjac(Iq)
            var_u = btp_mass_flux(1,Iq)
            var_v = btp_mass_flux(2,Iq)

            u_cori = coriolis_quad(Iq)*qb2(4,Iq)
            v_cori = -coriolis_quad(Iq)*qb2(3,Iq)

            do ip = 1,npts

                I = indexq(Iq,ip)

                hi = psih(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*hi*(stress_u(Iq) - pb_x(Iq) + u_cori)
                rhs_mom(2,I) = rhs_mom(2,I) + wq*hi*(stress_v(Iq) - pb_y(Iq) + v_cori)

                !Xi derivatives
                dhdx = dpsidx(Iq,ip)
                !Eta derivatives
                dhdy = dpsidy(Iq,ip)

                rhs_mom(1,I) = rhs_mom(1,I) + wq*dhdx*H(Iq)
                rhs_mom(2,I) = rhs_mom(2,I) + wq*dhdy*H(Iq)

                u_mom = Quu(Iq)*dhdx + Quv(Iq)*dhdy
                v_mom = Quv(Iq)*dhdx + Qvv(Iq)*dhdy

                rhs_mom(1,I) = rhs_mom(1,I) + wq*u_mom
                rhs_mom(2,I) = rhs_mom(2,I) + wq*v_mom

                rhs_mom(3,I) = rhs_mom(3,I) + wq*(dhdx*var_u + dhdy*var_v)

            end do
        end do

    end subroutine create_rhs_btp_dynamics_volume_new2


    subroutine Apply_btp_fluxes(rhs_mom,H_face,Qu_face,Qv_face,qb_face)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq, psiq, ngl
        use mod_grid, only:  npoin, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_input, only: nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe, coeff_mass_pbpert_LR, pbprime_face

        implicit none

        real, dimension(2, npoin), intent(inout)  :: rhs_mom
        real, dimension(2, nq, nface), intent(in) :: Qu_face, Qv_face
        real, dimension(nq, nface), intent(in) :: H_face
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face
    
        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, nq_i, nq_j
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_x, flux_y
        real, dimension(nq) :: quu, quv, qvu, qvv, Hface_l
        integer :: n, m
        real :: unl, unr, hbl, hbr, lambl, lambr, lamb, dispu, dispv, ubl, vbl, ubr, vbr

        quu = 0.0
        quv = 0.0
        qvu = 0.0
        qvv = 0.0
        Hface_l = 0.0

        do iface = 1, nface

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            Hface_l(:) = H_face(:,iface)

            quu(:) = Qu_face(1,:,iface)
            quv(:) = Qu_face(2,:,iface)

            qvu(:) = Qv_face(1,:,iface)
            qvv(:) = Qv_face(2,:,iface)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                H_kx = nxl*Hface_l(iquad)
                H_ky = nyl*Hface_l(iquad) 

                lamb = coeff_mass_pbpert_LR(iquad,iface)

                dispu = 0.5*lamb*(qb_face(3,2,iquad,iface) - qb_face(3,1,iquad,iface))
                dispv = 0.5*lamb*(qb_face(4,2,iquad,iface) - qb_face(4,1,iquad,iface))

                flux_x = nxl*quu(iquad) + nyl*quv(iquad) - dispu
                flux_y = nxl*qvu(iquad) + nyl*qvv(iquad) - dispv

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    rhs_mom(1,I) = rhs_mom(1,I) - wq*hi*(H_kx + flux_x)
                    rhs_mom(2,I) = rhs_mom(2,I) - wq*hi*(H_ky + flux_y)

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        rhs_mom(1,I) = rhs_mom(1,I) + wq*hi*(H_kx + flux_x)
                        rhs_mom(2,I) = rhs_mom(2,I) + wq*hi*(H_ky + flux_y)

                    end if

                end do
            end do

        end do

        rhs_mom(1,:) = massinv(:)*rhs_mom(1,:)
        rhs_mom(2,:) = massinv(:)*rhs_mom(2,:)
    
    end subroutine Apply_btp_fluxes

    subroutine Apply_btp_fluxes_new(rhs_mom,H_face,Qu_face,Qv_face,qb_face, flux_edge)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq, psiq, ngl
        use mod_grid, only:  npoin, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_input, only: nlayers
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe, coeff_mass_pbpert_LR, pbprime_face

        implicit none

        real, dimension(3, npoin), intent(inout)  :: rhs_mom
        real, dimension(2, nq, nface), intent(in) :: Qu_face, Qv_face, flux_edge
        real, dimension(nq, nface), intent(in) :: H_face
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face
    
        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, nq_i, nq_j
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_x, flux_y, flux
        real, dimension(nq) :: quu, quv, qvu, qvv, Hface_l, var_x, var_y
        integer :: n, m
        real :: unl, unr, hbl, hbr, lambl, lambr, lamb, dispu, dispv, ubl, vbl, ubr, vbr

        quu = 0.0
        quv = 0.0
        qvu = 0.0
        qvv = 0.0
        Hface_l = 0.0

        do iface = 1, nface

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            Hface_l(:) = H_face(:,iface)

            quu(:) = Qu_face(1,:,iface)
            quv(:) = Qu_face(2,:,iface)

            qvu(:) = Qv_face(1,:,iface)
            qvv(:) = Qv_face(2,:,iface)

            var_x = flux_edge(1,:,iface)
            var_y = flux_edge(2,:,iface)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                H_kx = nxl*Hface_l(iquad)
                H_ky = nyl*Hface_l(iquad) 

                lamb = coeff_mass_pbpert_LR(iquad,iface)

                dispu = 0.5*lamb*(qb_face(3,2,iquad,iface) - qb_face(3,1,iquad,iface))
                dispv = 0.5*lamb*(qb_face(4,2,iquad,iface) - qb_face(4,1,iquad,iface))

                flux_x = nxl*quu(iquad) + nyl*quv(iquad) - dispu
                flux_y = nxl*qvu(iquad) + nyl*qvv(iquad) - dispv

                flux = nxl*var_x(iquad) + nyl*var_y(iquad)

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    rhs_mom(1,I) = rhs_mom(1,I) - wq*hi*(H_kx + flux_x)
                    rhs_mom(2,I) = rhs_mom(2,I) - wq*hi*(H_ky + flux_y)
                    rhs_mom(3,I) = rhs_mom(3,I) - wq*hi*flux

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        rhs_mom(1,I) = rhs_mom(1,I) + wq*hi*(H_kx + flux_x)
                        rhs_mom(2,I) = rhs_mom(2,I) + wq*hi*(H_ky + flux_y)
                        rhs_mom(3,I) = rhs_mom(3,I) + wq*hi*flux

                    end if

                end do
            end do

        end do

        rhs_mom(1,:) = massinv(:)*rhs_mom(1,:)
        rhs_mom(2,:) = massinv(:)*rhs_mom(2,:)
        rhs_mom(3,:) = massinv(:)*rhs_mom(3,:)
    
    end subroutine Apply_btp_fluxes_new 

    subroutine Apply_btp_fluxes_new_v2(rhs_mom,H_face,Qu_face,Qv_face, flux_edge, grad_uvdp_face, qb_face)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq, psiq, ngl
        use mod_grid, only:  npoin, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_input, only: nlayers, visc_mlswe
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe, coeff_mass_pbpert_LR, pbprime_face

        implicit none

        real, dimension(3, npoin), intent(inout)  :: rhs_mom
        real, dimension(2, nq, nface), intent(in) :: Qu_face, Qv_face, flux_edge
        real, dimension(nq, nface), intent(in) :: H_face
        real, dimension(4,2,nq,nface), intent(in) :: grad_uvdp_face
        real, dimension(4, 2, nq, nface), intent(in) :: qb_face
    
        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, nq_i, nq_j
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j
        real :: wq, hi, nxl, nyl, H_kx, H_ky, flux_x, flux_y, flux, flux_qu, flux_qv
        real, dimension(nq) :: quu, quv, qvu, qvv, Hface_l, var_x, var_y
        integer :: n, m
        real :: unl, unr, hbl, hbr, lambl, lambr, lamb, dispu, dispv, ubl, vbl, ubr, vbr
        real, dimension(2) :: qu_mean, qv_mean

        quu = 0.0
        quv = 0.0
        qvu = 0.0
        qvv = 0.0
        Hface_l = 0.0

        do iface = 1, nface

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            Hface_l(:) = H_face(:,iface)

            quu(:) = Qu_face(1,:,iface)
            quv(:) = Qu_face(2,:,iface)

            qvu(:) = Qv_face(1,:,iface)
            qvv(:) = Qv_face(2,:,iface)

            var_x = flux_edge(1,:,iface)
            var_y = flux_edge(2,:,iface)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                H_kx = nxl*Hface_l(iquad)
                H_ky = nyl*Hface_l(iquad)

                lamb = coeff_mass_pbpert_LR(iquad,iface)

                dispu = 0.5*lamb*(qb_face(3,2,iquad,iface) - qb_face(3,1,iquad,iface))
                dispv = 0.5*lamb*(qb_face(4,2,iquad,iface) - qb_face(4,1,iquad,iface))

                flux_x = nxl*quu(iquad) + nyl*quv(iquad) - dispu
                flux_y = nxl*qvu(iquad) + nyl*qvv(iquad) - dispv

                flux = nxl*var_x(iquad) + nyl*var_y(iquad)

                qu_mean(:) = 0.5*(grad_uvdp_face(1:2,1,iquad,iface) + grad_uvdp_face(1:2,2,iquad,iface))
                qv_mean(:) = 0.5*(grad_uvdp_face(3:4,1,iquad,iface) + grad_uvdp_face(3:4,2,iquad,iface))

                flux_qu = nxl*qu_mean(1) + nyl*qu_mean(2)
                flux_qv = nxl*qv_mean(1) + nyl*qv_mean(2)

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    rhs_mom(1,I) = rhs_mom(1,I) - wq*hi*(H_kx + flux_x - visc_mlswe*flux_qu)
                    rhs_mom(2,I) = rhs_mom(2,I) - wq*hi*(H_ky + flux_y - visc_mlswe*flux_qv)
                    rhs_mom(3,I) = rhs_mom(3,I) - wq*hi*flux

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        rhs_mom(1,I) = rhs_mom(1,I) + wq*hi*(H_kx + flux_x - visc_mlswe*flux_qu)
                        rhs_mom(2,I) = rhs_mom(2,I) + wq*hi*(H_ky + flux_y - visc_mlswe*flux_qv)
                        rhs_mom(3,I) = rhs_mom(3,I) + wq*hi*flux

                    end if

                end do
            end do

        end do

        rhs_mom(1,:) = massinv(:)*rhs_mom(1,:)
        rhs_mom(2,:) = massinv(:)*rhs_mom(2,:)
        rhs_mom(3,:) = massinv(:)*rhs_mom(3,:)
    
    end subroutine Apply_btp_fluxes_new_v2

    subroutine create_rhs_btp_dynamics_volume_mass(pb_advec, btp_mass_flux)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma

        use mod_basis, only: npts

        use mod_constants, only: gravity

        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac


        implicit none 

        real, dimension(2, npoin_q), intent(in) :: btp_mass_flux
        real, dimension(npoin), intent(out) :: pb_advec


        real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, hi, dhde, dhdn
        real :: h_e, h_c, h_n, var_u, var_v

        integer :: k, e, iquad, jquad, kquad, l, n, m, I, Iq, ip

        pb_advec = 0.0

        do Iq = 1,npoin_q
            
            wq = wjac(Iq)

            var_u = btp_mass_flux(1,Iq)
            var_v = btp_mass_flux(2,Iq)

            do ip = 1,npts

                I = indexq(Iq,ip)

                pb_advec(I) = pb_advec(I) + wq*(dpsidx(Iq,ip)*var_u + dpsidy(Iq,ip)*var_v)

            end do
        end do

    end subroutine create_rhs_btp_dynamics_volume_mass

    subroutine Apply_btp_flux_mass(pb_advec, flux_edge)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq, psiq, ngl
        use mod_grid, only:  npoin_q, intma, mod_grid_get_face_nq, nface,face, mod_grid_get_face_ngl
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq

        implicit none

        real, dimension(npoin_q), intent(inout) :: pb_advec
        real, dimension(2, nq, nface), intent(in) :: flux_edge
    
        integer :: k, iface, iquad, el, er, il, jl, ir, jr, I
        integer :: plane_ij, kl, kr, jquad, ngl_i, ngl_j, n, m
        real :: wq, nxl, nyl, flux
        real, dimension(nq) :: var_x, var_y
        real :: hi
    
        do iface = 1, nface

            !Store Left Side Variables

            el = face(7,iface)
            er = face(8,iface)

            var_x = flux_edge(1,:,iface)
            var_y = flux_edge(2,:,iface)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                flux = nxl*var_x(iquad) + nyl*var_y(iquad)

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    pb_advec(I) = pb_advec(I) - wq*hi*flux

                    if(er > 0) then

                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        pb_advec(I) = pb_advec(I) + wq*hi*flux

                    end if !er > 0

                end do !n
            end do !iquad

        end do

        pb_advec(:) = massinv(:)*pb_advec(:)
    
    end subroutine Apply_btp_flux_mass

end module mod_create_rhs_mlswe
