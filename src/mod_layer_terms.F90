! ===========================================================================================================================
! This module contains the routines for the baroclinic flux terms
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: March 27, 2023
! ==========================================================================================================================

module mod_layer_terms
    
    use mod_grid, only: npoin_q, nface, intma_dg_quad, face
    use mod_basis, only: nqx, nqy, nqz, nq
    use mod_input, only: nlayers, flux_type
        
    implicit none

    public :: layer_mass_advection_terms, evaluate_dp, evaluate_dp_face, &
            compute_momentum_edge_values, layer_momentum_advec_terms, layer_pressure_terms, layer_windbot_stress_terms, &
            shear_stress_system, consistency_mass_terms1, consistency_mass_terms, layer_mom_boundary_df,  &
            filter_mlswe, layer_mom_boundary, evaluate_mom, velocity_df, evaluate_mom_face, &
            interpolate_mom, velocity, velocity_face, evaluate_mom_face_all, layer_momentum_advec_terms_upwind

    contains

    subroutine layer_mass_advection_terms(sum_layer_mass_flux,sum_layer_mass_flux_face,u_edge,v_edge,uvdp_temp,flux_edge, &
                        qprime,qprime_face)

        ! This routines compute the layer mass flux terms using upwind or centered flux

        use mod_grid, only: npoin_q, nface, intma_dg_quad, face
        use mod_basis, only: nqx, nqy, nqz, nq
        use mod_input, only: nlayers, bcl_flux
        use mod_face, only: imapl_q, imapr_q, normal_vector_q
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe
        use mod_variables, only: uvb_ave,ope_ave, uvb_face_ave,ope_face_ave

        implicit none

        real, dimension(3,npoin_q, nlayers),    intent(in)    :: qprime
        real, dimension(3,2,nq,nface,nlayers),  intent(in)    :: qprime_face
        real, dimension(2,npoin_q),             intent(out)   :: sum_layer_mass_flux
        real, dimension(2,nq,nface),          intent(out)     :: sum_layer_mass_flux_face
        real, dimension(2,nq, nface, nlayers),    intent(out) :: u_edge, v_edge
        real, dimension(2,npoin_q,nlayers),       intent(out) :: uvdp_temp
        real, dimension(2,nq, nface, nlayers),  intent(out)   :: flux_edge

        real, dimension(nq) :: dp_left, dp_right, dp_temp(npoin_q)
        integer :: k, iface, iquad
        real :: nxl, nyl, h_l, h_r, claml, clamr, clam, h1_l, h1_r, h2_l, h2_r, Ucl, Ucr, hgl, hgr, gprime, h1, h2, hg
        real, dimension(nlayers) :: unl, unr, u, v
        real :: clam1, clam2, s1_l, s1_r, s2_l, s2_r, alph, uu, vv, dpl, dpr

        alph = alpha_mlswe(1)/alpha_mlswe(2)

        gprime = gravity *alpha_mlswe(2)*(1.0/alpha_mlswe(2) - 1.0/alpha_mlswe(1))

        sum_layer_mass_flux = 0.0
        sum_layer_mass_flux_face = 0.0
        flux_edge = 0.0

        do k = 1, nlayers

            dp_temp = qprime(1,:,k) * ope_ave(:)
            uvdp_temp(1,:,k) = (qprime(2,:,k) + uvb_ave(1,:)) * dp_temp
            uvdp_temp(2,:,k) = (qprime(3,:,k) + uvb_ave(2,:)) * dp_temp

            sum_layer_mass_flux(1,:) = sum_layer_mass_flux(1,:) + uvdp_temp(1,:,k)
            sum_layer_mass_flux(2,:) = sum_layer_mass_flux(2,:) + uvdp_temp(2,:,k)

        end do

        if(flux_type == 'upwind') then 

            do iface = 1,nface

                do k = 1,nlayers

                    do iquad = 1,nq

                        ! In the following computation of fluxes at element faces,
                        ! flux_edge iface a numerical approximation to the mass flux at element faces (each face has nq quadrature points)
                        ! Here we are using centered fluxes, so we need to compute the fluxes at the left and right edges of each face

                        u_edge(1,iquad,iface,k) = qprime_face(2,1,iquad,iface,k) + uvb_face_ave(1,1,iquad,iface)
                        u_edge(2,iquad,iface,k) = qprime_face(2,2,iquad,iface,k) + uvb_face_ave(1,2,iquad,iface)
                        v_edge(1,iquad,iface,k) = qprime_face(3,1,iquad,iface,k) + uvb_face_ave(2,1,iquad,iface)
                        v_edge(2,iquad,iface,k) = qprime_face(3,2,iquad,iface,k) + uvb_face_ave(2,2,iquad,iface)

                        nxl = normal_vector_q(1,iquad,1,iface)
                        nyl = normal_vector_q(2,iquad,1,iface)

                        uu = 0.5*(u_edge(1,iquad,iface,k) + u_edge(2,iquad,iface,k))
                        vv = 0.5*(v_edge(1,iquad,iface,k) + v_edge(2,iquad,iface,k))
                        
                        dpl = ope_face_ave(1,iquad,iface) * qprime_face(1,1,iquad,iface,k)
                        dpr = ope_face_ave(2,iquad,iface) * qprime_face(1,2,iquad,iface,k)

                        if(uu*nxl > 0.0) then 
                            flux_edge(1,iquad,iface,k) = uu * dpl
                        else 
                            flux_edge(1,iquad,iface,k) = uu * dpr 
                        endif 

                        if(vv*nyl > 0.0) then 
                            flux_edge(2,iquad,iface,k) = vv * dpl
                        else 
                            flux_edge(2,iquad,iface,k) = vv * dpr 
                        endif 

                        sum_layer_mass_flux_face(1,iquad,iface) = sum_layer_mass_flux_face(1,iquad,iface) + flux_edge(1,iquad,iface,k)
                        sum_layer_mass_flux_face(2,iquad,iface) = sum_layer_mass_flux_face(2,iquad,iface) + flux_edge(2,iquad,iface,k)
                    end do 

                end do

            end do
        else 
            do iface = 1,nface

                do k = 1,nlayers

                    ! In the following computation of fluxes at element faces,
                    ! flux_edge iface a numerical approximation to the mass flux at element faces (each face has nq quadrature points)
                    ! Here we are using centered fluxes, so we need to compute the fluxes at the left and right edges of each face

                    u_edge(1,:,iface,k) = qprime_face(2,1,:,iface,k) + uvb_face_ave(1,1,:,iface)
                    u_edge(2,:,iface,k) = qprime_face(2,2,:,iface,k) + uvb_face_ave(1,2,:,iface)
                    v_edge(1,:,iface,k) = qprime_face(3,1,:,iface,k) + uvb_face_ave(2,1,:,iface)
                    v_edge(2,:,iface,k) = qprime_face(3,2,:,iface,k) + uvb_face_ave(2,2,:,iface)

                    dp_left = ope_face_ave(1,:,iface) * qprime_face(1,1,:,iface,k)
                    dp_right = ope_face_ave(2,:,iface) * qprime_face(1,2,:,iface,k)

                    flux_edge(1,:,iface,k) = 0.5*(u_edge(1,:,iface,k) * dp_left + u_edge(2,:,iface,k) * dp_right)
                    flux_edge(2,:,iface,k) = 0.5*(v_edge(1,:,iface,k) * dp_left + v_edge(2,:,iface,k) * dp_right)

                    sum_layer_mass_flux_face(1,:,iface) = sum_layer_mass_flux_face(1,:,iface) + flux_edge(1,:,iface,k)
                    sum_layer_mass_flux_face(2,:,iface) = sum_layer_mass_flux_face(2,:,iface) + flux_edge(2,:,iface,k)

                end do

            end do
        end if 

    end subroutine layer_mass_advection_terms

    subroutine evaluate_dp(q,qprime,q_df,pbprime)
    
        ! This routine interpolate layer mass from dofs to quad points 

        use mod_basis, only: npts
        use mod_grid, only: npoin, npoin_q, intma, intma_dg_quad
        use mod_input, only: nlayers

        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac

        implicit none
        real, dimension(3,npoin_q,nlayers), intent(inout) :: q
        real, dimension(3,npoin_q,nlayers), intent(inout) :: qprime
        real, dimension(3,npoin,nlayers), intent(in) :: q_df
        real, dimension(npoin_q), intent(in) :: pbprime

        integer :: k, Iq,I, ip
        real :: hi
        real, dimension(npoin_q) :: pb_temp, one_plus_eta_temp


        qprime(1,:,:) = 0.0
        q(1,:,:) = 0.0
        pb_temp = 0.0

        ! do k = 1, nlayers

        do Iq = 1,npoin_q   
            do ip = 1,npts

                I = indexq(Iq,ip)
                hi = psih(Iq,ip)

                q(1,Iq,:) = q(1,Iq,:) + q_df(1,I,:)*hi

            end do
            pb_temp(Iq) = pb_temp(Iq) + sum(q(1,Iq,:))

        end do
        ! end do
      
        one_plus_eta_temp(:) = pb_temp(:) / pbprime(:)
        do k = 1, nlayers
            qprime(1,:,k) = q(1,:,k) / one_plus_eta_temp(:)
        end do

    end subroutine evaluate_dp

    subroutine evaluate_dp_face(q_face, qprime_face,q,qprime)

        ! This routines extracts the layer mass face values 

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq
        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq, nface, face
        !use mod_initial, only: pbprime, pbprime_face
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q

        implicit none
    
        ! Input arguments
        real, intent(in) :: q(3,npoin_q,nlayers)
        ! Output arguments
        real, intent(inout) :: q_face(3,2,nq,nface,nlayers)
        real, intent(inout) :: qprime_face(3,2,nq,nface,nlayers)
        real, intent(in) :: qprime(3,npoin_q,nlayers)
    
        ! Local variables
        real :: pb_temp_face(2,nq, nface)
        integer :: iface, iquad, jquad, k, il, jl, kl, el, er, ilocl, ilocr, ir, jr, kr, I
        real :: one_plus_eta_temp(2)

        ! do k = 1,nlayers
        qprime_face(1,:,:,:,:) = 0.0
        q_face(1,:,:,:,:) = 0.0

        do iface = 1, nface

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                il = imapl_q(1,iquad,1,iface)
                jl = imapl_q(2,iquad,1,iface)
                kl = imapl_q(3,iquad,1,iface)

                I = intma_dg_quad(il,jl,kl,el)

                q_face(1,1,iquad,iface,:) = q(1,I,:)
                qprime_face(1,1,iquad,iface,:) = qprime(1,I,:)

                if(er > 0) then

                    ir = imapr_q(1,iquad,1,iface)
                    jr = imapr_q(2,iquad,1,iface)
                    kr = imapr_q(3,iquad,1,iface)

                    I = intma_dg_quad(ir,jr,kr,er)

                    q_face(1,2,iquad,iface,:) = q(1,I,:)
                    qprime_face(1,2,iquad,iface,:) = qprime(1,I,:)

                else

                    q_face(1,2,iquad,iface,:) = q_face(1,1,iquad,iface,:)
                    qprime_face(1,2,iquad,iface,:) = qprime_face(1,1,iquad,iface,:)

                end if
            end do
        end do
        ! end do
    
    end subroutine evaluate_dp_face
    
    subroutine consistency_mass_terms(flux_adjustment, flux_adjust_edge, qprime, &
        sum_layer_mass_flux, qprime_face, flux_deficit_mass_face)

        ! This routines computes the consistency flux terms 

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, nq, ngl
        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq,nface,face, npoin, intma, nelem
        use mod_initial, only: pbprime, pbprime_face
        use mod_input, only: nlayers
        use mod_variables, only: btp_mass_flux_ave, ope_face_ave
        use mod_face, only: normal_vector_q

        implicit none

        real, dimension(2,npoin_q, nlayers), intent(out)      :: flux_adjustment
        real, dimension(2,nq, nface, nlayers), intent(out)    :: flux_adjust_edge

        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real,dimension(2,npoin_q), intent(in)          :: sum_layer_mass_flux
        real, dimension(2,2,nq,nface), intent(in) :: flux_deficit_mass_face
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face

        real, dimension(npoin_q) :: weights
        !real, dimension(nq) :: weights_face_l, weights_face_r, weight
        integer :: iface, k, iquad
        real :: weights_face_l, weights_face_r, fu, fv, nxl, nyl


        flux_adjust_edge = 0.0

        ! Consistency terms for interior element points
        do k = 1, nlayers

            weights(:) = qprime(1,:,k) / pbprime(:)

            flux_adjustment(1,:,k) = weights(:) * (btp_mass_flux_ave(1,:) - sum_layer_mass_flux(1,:))
            flux_adjustment(2,:,k) = weights(:) * (btp_mass_flux_ave(2,:) - sum_layer_mass_flux(2,:))   

        end do

        ! Face terms 

        do iface = 1,nface

            do iquad = 1,nq 
                do k = 1,nlayers

                    weights_face_l = qprime_face(1,1,iquad,iface,k) / pbprime_face(1,iquad,iface)
                    weights_face_r = qprime_face(1,2,iquad,iface,k) / pbprime_face(2,iquad,iface)

                    !fu = 0.5*(flux_deficit_mass_face(1,1,iquad,iface) + flux_deficit_mass_face(1,2,iquad,iface))
                    !fv = 0.5*(flux_deficit_mass_face(2,1,iquad,iface) + flux_deficit_mass_face(2,2,iquad,iface))

                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    if(flux_deficit_mass_face(1,1,iquad,iface)*nxl > 0.0) then 

                        flux_adjust_edge(1,iquad,iface,k) = weights_face_l*flux_deficit_mass_face(1,1,iquad,iface)
                    else 
                        flux_adjust_edge(1,iquad,iface,k) = weights_face_r*flux_deficit_mass_face(1,2,iquad,iface)
                    end if 

                    if(flux_deficit_mass_face(2,1,iquad,iface)*nyl > 0.0) then 
                        flux_adjust_edge(2,iquad,iface,k) = weights_face_l*flux_deficit_mass_face(2,1,iquad,iface)
                    else 
                        flux_adjust_edge(2,iquad,iface,k) = weights_face_r*flux_deficit_mass_face(2,2,iquad,iface)
                    end if 
                end do 
            end do 
        end do

    end subroutine consistency_mass_terms

    subroutine consistency_mass_terms1(flux_adjustment, flux_adjust_edge, q_df, qprime, &
        sum_layer_mass_flux, flux_deficit_mass_face, q_face)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, nq, ngl
        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq,nface,face, npoin, intma
        use mod_initial, only: pbprime, pbprime_face, pbprime_edge
        use mod_input, only: nlayers
        use mod_face, only: normal_vector_q
        use mod_variables, only: ope_face_ave, btp_mass_flux_ave

        implicit none

        real, dimension(2,npoin_q, nlayers), intent(out)      :: flux_adjustment
        real, dimension(2,nq, nface, nlayers), intent(out)    :: flux_adjust_edge

        real, dimension(3,npoin,nlayers), intent(in)   :: q_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real,dimension(2,npoin_q), intent(in)          :: sum_layer_mass_flux
        real, dimension(2,2,nq,nface), intent(in)        :: flux_deficit_mass_face
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: q_face

        real, dimension(2,2) :: flux_adjustment_face
        real, dimension(nlayers)                  :: dp_avail_left, dp_avail_right
        real :: sum_left, sum_right, pb_edge, p_left, p_right, temp, one_plus_eta_edge
        integer :: el, er, il, jl, iquad, k, iface
        real, dimension(npoin_q) :: weights
        real :: weights_face, qml, qmr, weights_face_r, weights_face_l, nxl, nyl
        integer :: I, ir, jr, kr, kl, ier,m, Iq

        flux_adjustment = 0.0
        flux_adjust_edge = 0.0
        flux_adjustment_face = 0.0

        ! Consistency terms for interior cell points
        do k = 1, nlayers
            weights(:) = qprime(1,:,k) / pbprime(:)
            flux_adjustment(1,:,k) = weights(:) * (btp_mass_flux_ave(1,:) - sum_layer_mass_flux(1,:))
            flux_adjustment(2,:,k) = weights(:) * (btp_mass_flux_ave(2,:) - sum_layer_mass_flux(2,:))
        end do

        ! Consistency terms for face points
        do iface = 1,nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)
                
            do iquad = 1,nq

                one_plus_eta_edge = 0.5*(ope_face_ave(2,iquad,iface) + ope_face_ave(1,iquad,iface))
                pb_edge = pbprime_edge(iquad,iface) * one_plus_eta_edge

                dp_avail_left = 0.0
                dp_avail_right = 0.0

                sum_left = 0.0
                sum_right = 0.0

                p_left = 0.0 
                p_right = 0.0
                    
                do k=1,nlayers
                        
                    temp = max(pb_edge - p_left, 0.0)

                    qml = q_face(1,1,iquad,iface,k) 
                    dp_avail_left(k)  = min(temp, qml)
                    sum_left  = sum_left + dp_avail_left(k)
                    p_left  = p_left  + qml

                    temp = max(pb_edge - p_right, 0.0)

                    qmr = q_face(1,2,iquad,iface,k)
                    dp_avail_right(k) = min(temp, qmr)
                    sum_right = sum_right + dp_avail_right(k)
                    p_right = p_right + qmr

                end do

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                do k = 1,nlayers

                    if (sum_left > 0.0) then
                        weights_face_l = dp_avail_left(k) / sum_left
                        flux_adjustment_face(1,1) = weights_face_l * flux_deficit_mass_face(1,1,iquad,iface)
                        flux_adjustment_face(2,1) = weights_face_l * flux_deficit_mass_face(2,1,iquad,iface)
                    end if
                    
                    if (sum_right > 0.0) then
                        weights_face_r = dp_avail_right(k) / sum_right
                        flux_adjustment_face(1,2) = weights_face_r * flux_deficit_mass_face(1,2,iquad,iface)
                        flux_adjustment_face(2,2) = weights_face_r * flux_deficit_mass_face(2,2,iquad,iface)
                    end if

                    if (flux_deficit_mass_face(1,1,iquad,iface)*nxl > 0.0) then
                        flux_adjust_edge(1,iquad,iface,k) = flux_adjustment_face(1,1)
                    else
                        flux_adjust_edge(1,iquad,iface,k) = flux_adjustment_face(1,2)
                    end if

                    if (flux_deficit_mass_face(2,1,iquad,iface)*nyl > 0.0) then
                        flux_adjust_edge(2,iquad,iface,k) = flux_adjustment_face(2,1)
                    else
                        flux_adjust_edge(2,iquad,iface,k) = flux_adjustment_face(2,2)
                    end if

                end do

            end do
        end do

    end subroutine consistency_mass_terms1

    subroutine compute_momentum_edge_values(udp_left, vdp_left, udp_right, vdp_right, qprime_face)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, nq
        use mod_grid, only:  nface
        use mod_input, only: nlayers, bcl_flux
        use mod_face, only: imapl_q, imapr_q, normal_vector_q
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe
        use mod_variables, only: uvb_face_ave, ope_face_ave
  
        implicit none
        
        ! Input variables
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face 
        
        ! Output variables
        real, dimension(nq,nface,nlayers), intent(out) :: udp_left, vdp_left, udp_right, vdp_right
        
        ! Local variables
        integer :: k, iface, iquad
        real :: one_plus_eta_left, dpprime_left
        real, dimension(nq) :: u_left, v_left, u_right, v_right, dp_left, dp_right
        real :: one_plus_eta_right, dpprime_right, h_l, h_r, claml, clamr, clam, nxl, nyl
        real :: h1_l, h1_r, h2_l, h2_r, Ucl, Ucr, hgl, hgr, gprime, h1, h2, UU, hg, clam1, clam2
        real, dimension(nlayers) :: unl, unr, ul, ur, vl, vr, u, v

        gprime = gravity *alpha_mlswe(2)*(1.0/alpha_mlswe(2) - 1.0/alpha_mlswe(1))
        
        ! Compute the momentum values at each cell edge
        
        do iface = 1, nface
            do k = 1, nlayers

                ! Left side of the edge

                u_left = qprime_face(2,1,:,iface,k) + uvb_face_ave(1,1,:,iface)
                v_left = qprime_face(3,1,:,iface,k) + uvb_face_ave(2,1,:,iface)
                dp_left = qprime_face(1,1,:,iface,k) * ope_face_ave(1,:,iface)
                
                udp_left(:,iface,k) = u_left * dp_left
                vdp_left(:,iface,k) = v_left * dp_left
                
                ! Right side of the edge
                u_right = qprime_face(2,2,:,iface,k) + uvb_face_ave(1,2,:,iface)
                v_right = qprime_face(3,2,:,iface,k) + uvb_face_ave(2,2,:,iface)
                dp_right = qprime_face(1,2,:,iface,k) * ope_face_ave(2,:,iface)
                
                udp_right(:,iface,k) = u_right * dp_right
                vdp_right(:,iface,k) = v_right * dp_right
            end do
        end do
        
    end subroutine compute_momentum_edge_values

    subroutine layer_momentum_advec_terms(u_udp_temp, u_vdp_temp, v_vdp_temp, udp_flux_edge, vdp_flux_edge, &
        q, qprime, u_edge, v_edge, udp_left, vdp_left, udp_right, vdp_right)

        ! This routine computes the layer momentum advection flux terms using centered flux 

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, nq
        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq,nface
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q, normal_vector_q
        use mod_variables, only: Quv_face_ave, Qu_face_ave, Qv_face_ave, ope_ave, Qu_ave, Qv_ave, Quv_ave, uvb_ave

        implicit none

        real, dimension(3, npoin_q, nlayers), intent(in)    :: q 
        real, dimension(3, npoin_q, nlayers), intent(in)    :: qprime
        real, dimension(2,nq, nface, nlayers), intent(in)     :: u_edge, v_edge
        real, dimension(nq, nface, nlayers), intent(in)     :: udp_left, vdp_left, udp_right, vdp_right

        real, dimension(npoin_q, nlayers), intent(out)      :: u_udp_temp, v_vdp_temp
        real, dimension(2, npoin_q, nlayers), intent(out)   :: u_vdp_temp
        real, dimension(2, nq, nface, nlayers), intent(out) :: udp_flux_edge, vdp_flux_edge

        real, dimension(nlayers) :: udp_abs_temp, vdp_abs_temp
        real, dimension(npoin_q) :: uu_dp_flux_deficitq, uv_dp_flux_deficitq, vv_dp_flux_deficitq, one_over_sumuq, one_over_sumvq, weightq
        real :: sum_uu, sum_uv, sum_vv, uu_dp_flux_deficit, uv_dp_flux_deficit, vv_dp_flux_deficit
        real :: sumu, sumv, one_over_sumu, one_over_sumv, weight
        real, dimension(npoin_q) :: sumu1, sumv1, sumuv1
        integer :: k, iface, I, iquad , Iq, el, er, il, jl, kl, mom_consistency
        real :: vu_dp_flux_deficit, sum_vu, u_min, u_max, v_min, v_max
        real, parameter :: eps1 = 1.0e-12 !  Parameter used to prevent division by zero.
        real :: nx, ny, un(nlayers), uu, vv, nxl, nyl
        real, dimension(npoin_q) :: temp_u, temp_v, temp_dp
        real, dimension(npoin_q,nlayers) :: temp_uu, temp_vv

        mom_consistency = 1

        do k = 1, nlayers

            temp_dp = qprime(1,:,k) * ope_ave(:)

            temp_u = qprime(2,:,k) + uvb_ave(1,:)
            u_udp_temp(:, k) = temp_dp * temp_u**2 
            
            temp_v = qprime(3,:,k) + uvb_ave(2,:)
            v_vdp_temp(:, k) = temp_dp * temp_v**2 

            u_vdp_temp(1,:,k) = temp_u * temp_v * temp_dp
            u_vdp_temp(2,:,k) = temp_v * temp_u * temp_dp

        end do

        ! In the following computation of momentum fluxes at cell edges, 
        ! udp_flux_edge  iface a numerical approximation to the flux of udp
        ! at cells edge in layer k, and vdp_flux_edge  iface a numerical approximation to the flux of vdp.

        ! *****   Flux, with available momentum   *****

        do k = 1,nlayers
            do iface = 1,nface

                do iquad = 1, nq

                    udp_flux_edge(1,iquad,iface,k) = 0.5*(udp_left(iquad,iface,k) * u_edge(1,iquad,iface,k) + udp_right(iquad,iface,k) * u_edge(2,iquad,iface,k))
                    udp_flux_edge(2,iquad,iface,k) = 0.5*(udp_left(iquad,iface,k) * v_edge(1,iquad,iface,k) + udp_right(iquad,iface,k) * v_edge(2,iquad,iface,k))

                    vdp_flux_edge(1,iquad,iface,k) = 0.5*(vdp_left(iquad,iface,k) * u_edge(1,iquad,iface,k) + vdp_right(iquad,iface,k) * u_edge(2,iquad,iface,k))
                    vdp_flux_edge(2,iquad,iface,k) = 0.5*(vdp_left(iquad,iface,k) * v_edge(1,iquad,iface,k) + vdp_right(iquad,iface,k) * v_edge(2,iquad,iface,k))
                end do

            end do
        end do

        ! *****   End computation of fluxes at cell edges   *****

        if(mom_consistency == 1) then 

            uu_dp_flux_deficitq = Qu_ave(:) - sum(u_udp_temp(:,:), dim=2)
            uv_dp_flux_deficitq = Quv_ave(:) - sum(u_vdp_temp(1,:,:), dim=2)
            vv_dp_flux_deficitq = Qv_ave(:) - sum(v_vdp_temp(:,:), dim=2)

            temp_uu(:,:) = abs(q(2,:,:)) 
            temp_vv(:,:) = abs(q(3,:,:))
            
            one_over_sumuq = 1.0 / (sum(temp_uu(:,:), dim=2) + eps1)
            one_over_sumvq = 1.0 / (sum(temp_vv(:,:), dim=2) + eps1)
            
            do k=1,nlayers
        
                weightq = temp_uu(:,k) * one_over_sumuq
                u_udp_temp(:,k) = u_udp_temp(:,k) + weightq * uu_dp_flux_deficitq
                u_vdp_temp(1,:,k) = u_vdp_temp(1,:,k) + weightq * uv_dp_flux_deficitq
        
                weightq = temp_vv(:,k) * one_over_sumvq
                u_vdp_temp(2,:,k) = u_vdp_temp(2,:,k) + weightq * uv_dp_flux_deficitq
                v_vdp_temp(:,k) = v_vdp_temp(:,k) + weightq * vv_dp_flux_deficitq
            enddo

            do iface = 1,nface   ! loop over cell edges

                !Store Left Side Variables
                el = face(7,iface)
                er = face(8,iface)

                do iquad = 1,nq

                    sum_uu = 0.0
                    sum_uv = 0.0
                    sum_vu = 0.0
                    sum_vv = 0.0

                    do k = 1,nlayers

                        sum_uu = sum_uu + udp_flux_edge(1,iquad,iface,k)
                        sum_uv = sum_uv + udp_flux_edge(2,iquad,iface,k)
                        sum_vu = sum_vu + vdp_flux_edge(1,iquad,iface,k)
                        sum_vv = sum_vv + vdp_flux_edge(2,iquad,iface,k)

                    end do

                    uu_dp_flux_deficit = Qu_face_ave(1,iquad,iface) - sum_uu
                    uv_dp_flux_deficit = Qu_face_ave(2,iquad,iface) - sum_uv
                    vu_dp_flux_deficit = Qv_face_ave(1,iquad,iface) - sum_vu
                    vv_dp_flux_deficit = Qv_face_ave(2,iquad,iface) - sum_vv

                    ! Adjust the fluxes for the u-momentum equation

                    sumu = 0.0

                    do k = 1,nlayers
                        udp_abs_temp(k) = 0.5*(abs(udp_left(iquad,iface,k)) + abs(udp_right(iquad,iface,k))) !+ eps1
                        sumu = sumu + udp_abs_temp(k)
                    end do
                    one_over_sumu = 1.0 / (sumu +eps1)

                    do k = 1,nlayers
                        weight = udp_abs_temp(k) * one_over_sumu
                        udp_flux_edge(1,iquad,iface,k) = udp_flux_edge(1,iquad,iface,k) + weight * uu_dp_flux_deficit
                        udp_flux_edge(2,iquad,iface,k) = udp_flux_edge(2,iquad,iface,k) + weight * uv_dp_flux_deficit
                    end do

                    ! Adjust the fluxes for the v-momentum equation

                    sumv = 0.0

                    do k = 1,nlayers
                        vdp_abs_temp(k) = 0.5*(abs(vdp_left(iquad,iface,k)) + abs(vdp_right(iquad,iface,k))) !+ eps1
                        sumv = sumv + vdp_abs_temp(k)
                    end do
                    one_over_sumv = 1.0 / (sumv + eps1)

                    do k = 1,nlayers
                        weight = vdp_abs_temp(k) * one_over_sumv
                        vdp_flux_edge(1,iquad,iface,k) = vdp_flux_edge(1,iquad,iface,k) + weight * vu_dp_flux_deficit
                        vdp_flux_edge(2,iquad,iface,k) = vdp_flux_edge(2,iquad,iface,k) + weight * vv_dp_flux_deficit
                    end do

                end do   ! end loop over quadrature points
            end do   ! end loop over cell edges
        end if
          
    end subroutine layer_momentum_advec_terms

    subroutine layer_momentum_advec_terms_upwind(u_udp_temp, u_vdp_temp, v_vdp_temp, udp_flux_edge, vdp_flux_edge, &
        q, qprime, u_edge, v_edge, udp_left, vdp_left, udp_right, vdp_right)

        ! This routine computes the layer momentum advection flux terms using upwind flux 

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, nq
        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq,nface
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q, normal_vector_q
        use mod_variables, only: Quv_face_ave, Qu_face_ave, Qv_face_ave, ope_ave, Qu_ave, Qv_ave, Quv_ave, uvb_ave

        implicit none

        real, dimension(3, npoin_q, nlayers), intent(in)    :: q 
        real, dimension(3, npoin_q, nlayers), intent(in)    :: qprime
        real, dimension(2,nq, nface, nlayers), intent(in)     :: u_edge, v_edge
        real, dimension(nq, nface, nlayers), intent(in)     :: udp_left, vdp_left, udp_right, vdp_right

        real, dimension(npoin_q, nlayers), intent(out)      :: u_udp_temp, v_vdp_temp
        real, dimension(2, npoin_q, nlayers), intent(out)   :: u_vdp_temp
        real, dimension(2, nq, nface, nlayers), intent(out) :: udp_flux_edge, vdp_flux_edge

        real, dimension(nlayers) :: udp_abs_temp, vdp_abs_temp
        real, dimension(npoin_q) :: uu_dp_flux_deficitq, uv_dp_flux_deficitq, vv_dp_flux_deficitq, one_over_sumuq, one_over_sumvq, weightq
        real :: sum_uu, sum_uv, sum_vv, uu_dp_flux_deficit, uv_dp_flux_deficit, vv_dp_flux_deficit
        real :: sumu, sumv, one_over_sumu, one_over_sumv, weight
        real, dimension(npoin_q) :: sumu1, sumv1, sumuv1
        integer :: k, iface, I, iquad , Iq, el, er, il, jl, kl, mom_consistency
        real :: vu_dp_flux_deficit, sum_vu, u_min, u_max, v_min, v_max
        real, parameter :: eps1 = 1.0e-12 !  Parameter used to prevent division by zero.
        real :: un(nlayers), uu, vv, nxl, nyl
        real, dimension(npoin_q) :: temp_u, temp_v, temp_dp
        real, dimension(npoin_q,nlayers) :: temp_uu, temp_vv

        mom_consistency = 1

        do k = 1, nlayers

            temp_dp = qprime(1,:,k) * ope_ave(:)

            temp_u = qprime(2,:,k) + uvb_ave(1,:)
            u_udp_temp(:, k) = temp_dp * temp_u**2 
            
            temp_v = qprime(3,:,k) + uvb_ave(2,:)
            v_vdp_temp(:, k) = temp_dp * temp_v**2 

            u_vdp_temp(1,:,k) = temp_u * temp_v * temp_dp
            u_vdp_temp(2,:,k) = temp_v * temp_u * temp_dp

        end do

        ! In the following computation of momentum fluxes at cell edges, 
        ! udp_flux_edge  iface a numerical approximation to the flux of udp
        ! at cells edge in layer k, and vdp_flux_edge  iface a numerical approximation to the flux of vdp.

        ! *****   Flux, with available momentum   *****

        do k = 1,nlayers
            do iface = 1,nface

                do iquad = 1, nq

                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    uu = 0.5*(u_edge(1,iquad,iface,k) + u_edge(2,iquad,iface,k))
                    vv = 0.5*(v_edge(1,iquad,iface,k) + v_edge(2,iquad,iface,k))

                    if(uu*nxl > 0.0) then 
                        udp_flux_edge(1,iquad,iface,k) = uu * udp_left(iquad,iface,k)
                        vdp_flux_edge(1,iquad,iface,k) = uu * vdp_left(iquad,iface,k)
                    else 
                        udp_flux_edge(1,iquad,iface,k) = uu * udp_right(iquad,iface,k)
                        vdp_flux_edge(1,iquad,iface,k) = uu * vdp_right(iquad,iface,k)
                    endif 

                    if(vv*nyl > 0.0) then 
                        udp_flux_edge(2,iquad,iface,k) = vv * udp_left(iquad,iface,k)
                        vdp_flux_edge(2,iquad,iface,k) = vv * vdp_left(iquad,iface,k)
                    else 
                        udp_flux_edge(2,iquad,iface,k) = vv * udp_right(iquad,iface,k)
                        vdp_flux_edge(2,iquad,iface,k) = vv * vdp_right(iquad,iface,k)
                    endif 
                end do

            end do
        end do

        ! *****   End computation of fluxes at cell edges   *****

        if(mom_consistency == 1) then 

            uu_dp_flux_deficitq = Qu_ave(:) - sum(u_udp_temp(:,:), dim=2)
            uv_dp_flux_deficitq = Quv_ave(:) - sum(u_vdp_temp(1,:,:), dim=2)
            vv_dp_flux_deficitq = Qv_ave(:) - sum(v_vdp_temp(:,:), dim=2)

            one_over_sumuq = 0.0
            one_over_sumvq = 0.0

            do k = 1,nlayers
                temp_uu(:,k) = abs(q(2,:,k)) + eps1
                temp_vv(:,k) = abs(q(3,:,k)) + eps1
                
                one_over_sumuq = one_over_sumuq + temp_uu(:,k)
                one_over_sumvq = one_over_sumvq + temp_vv(:,k)
            end do 

            one_over_sumuq = 1.0/one_over_sumuq
            one_over_sumvq = 1.0/one_over_sumvq
            
            do k=1,nlayers
        
                weightq = temp_uu(:,k) * one_over_sumuq
                u_udp_temp(:,k) = u_udp_temp(:,k) + weightq * uu_dp_flux_deficitq
                u_vdp_temp(1,:,k) = u_vdp_temp(1,:,k) + weightq * uv_dp_flux_deficitq
        
                weightq = temp_vv(:,k) * one_over_sumvq
                u_vdp_temp(2,:,k) = u_vdp_temp(2,:,k) + weightq * uv_dp_flux_deficitq
                v_vdp_temp(:,k) = v_vdp_temp(:,k) + weightq * vv_dp_flux_deficitq
            enddo

            do iface = 1,nface   ! loop over cell edges

                !Store Left Side Variables
                el = face(7,iface)
                er = face(8,iface)

                do iquad = 1,nq

                    nxl = normal_vector_q(1,iquad,1,iface)
                    nyl = normal_vector_q(2,iquad,1,iface)

                    sum_uu = 0.0
                    sum_uv = 0.0
                    sum_vu = 0.0
                    sum_vv = 0.0

                    do k = 1,nlayers

                        sum_uu = sum_uu + udp_flux_edge(1,iquad,iface,k)
                        sum_uv = sum_uv + udp_flux_edge(2,iquad,iface,k)
                        sum_vu = sum_vu + vdp_flux_edge(1,iquad,iface,k)
                        sum_vv = sum_vv + vdp_flux_edge(2,iquad,iface,k)

                    end do

                    uu_dp_flux_deficit = Qu_face_ave(1,iquad,iface) - sum_uu
                    uv_dp_flux_deficit = Qu_face_ave(2,iquad,iface) - sum_uv
                    vu_dp_flux_deficit = Qv_face_ave(1,iquad,iface) - sum_vu
                    vv_dp_flux_deficit = Qv_face_ave(2,iquad,iface) - sum_vv

                    ! Adjust the fluxes for the u-momentum equation
                    !x-direction

                    if(uu_dp_flux_deficit*nxl > 0.0) then 
                        sumu = 0.0

                        do k = 1,nlayers
                            udp_abs_temp(k) = abs(udp_left(iquad,iface,k)) + eps1
                            sumu = sumu + udp_abs_temp(k)
                        end do
                        one_over_sumu = 1.0 / sumu
                    else 
                        sumu = 0.0

                        do k = 1,nlayers
                            udp_abs_temp(k) = abs(udp_right(iquad,iface,k)) + eps1
                            sumu = sumu + udp_abs_temp(k)
                        end do
                        one_over_sumu = 1.0 / sumu
                    end if 

                    do k = 1,nlayers
                        weight = udp_abs_temp(k) * one_over_sumu
                        udp_flux_edge(1,iquad,iface,k) = udp_flux_edge(1,iquad,iface,k) + weight * uu_dp_flux_deficit
                    end do

                    !y-direction

                    if(uv_dp_flux_deficit*nyl > 0.0) then 
                        sumu = 0.0

                        do k = 1,nlayers
                            udp_abs_temp(k) = abs(udp_left(iquad,iface,k)) + eps1
                            sumu = sumu + udp_abs_temp(k)
                        end do
                        one_over_sumu = 1.0 / sumu
                    else 
                        sumu = 0.0

                        do k = 1,nlayers
                            udp_abs_temp(k) = abs(udp_right(iquad,iface,k)) + eps1
                            sumu = sumu + udp_abs_temp(k)
                        end do
                        one_over_sumu = 1.0 / sumu
                    end if 

                    do k = 1,nlayers
                        weight = udp_abs_temp(k) * one_over_sumu
                        udp_flux_edge(2,iquad,iface,k) = udp_flux_edge(2,iquad,iface,k) + weight * uv_dp_flux_deficit
                    end do

                    ! Adjust the fluxes for the v-momentum equation
                    !x-direction

                    if (vu_dp_flux_deficit*nxl > 0.0) then 
                        sumv = 0.0

                        do k = 1,nlayers
                            vdp_abs_temp(k) = abs(vdp_left(iquad,iface,k)) + eps1
                            sumv = sumv + vdp_abs_temp(k)
                        end do
                        one_over_sumv = 1.0 / sumv
                    else 
                        sumv = 0.0

                        do k = 1,nlayers
                            vdp_abs_temp(k) = abs(vdp_right(iquad,iface,k)) + eps1
                            sumv = sumv + vdp_abs_temp(k)
                        end do
                        one_over_sumv = 1.0 / sumv
                    end if 

                    do k = 1,nlayers
                        weight = vdp_abs_temp(k) * one_over_sumv
                        vdp_flux_edge(1,iquad,iface,k) = vdp_flux_edge(1,iquad,iface,k) + weight * vu_dp_flux_deficit
                    end do

                    !y-direction

                    if (vv_dp_flux_deficit*nyl > 0.0) then 
                        sumv = 0.0

                        do k = 1,nlayers
                            vdp_abs_temp(k) = abs(vdp_left(iquad,iface,k)) + eps1
                            sumv = sumv + vdp_abs_temp(k)
                        end do
                        one_over_sumv = 1.0 / sumv
                    else 
                        sumv = 0.0

                        do k = 1,nlayers
                            vdp_abs_temp(k) = abs(vdp_right(iquad,iface,k)) + eps1
                            sumv = sumv + vdp_abs_temp(k)
                        end do
                        one_over_sumv = 1.0 / sumv
                    end if 

                    do k = 1,nlayers
                        weight = vdp_abs_temp(k) * one_over_sumv
                        vdp_flux_edge(2,iquad,iface,k) = vdp_flux_edge(2,iquad,iface,k) + weight * vv_dp_flux_deficit
                    end do

                end do   ! end loop over quadrature points
            end do   ! end loop over cell edges
        end if
          
    end subroutine layer_momentum_advec_terms_upwind

    subroutine layer_pressure_terms(H_r, H_r_face, p, z_elev, qprime, qprime_face, zbot_face, qprime_df ,grad_z)

        ! This routine computes the layer momentum pressure terms

        use mod_constants, only : gravity
        use mod_initial, only : alpha_mlswe, zbot_df, pbprime, pbprime_face, zbot
        use mod_grid, only : nface, npoin, npoin_q, face, intma_dg_quad
        use mod_basis, only : nq
        use mod_input, only : nlayers, adjust_H_vertical_sum
        use mod_face, only: imapl_q, imapr_q
        use mod_Tensorproduct, only: interpolate_layer_from_quad_to_node_1d, compute_gradient_quad
        use mod_variables, only: ope_ave, H_ave, ope_face_ave, H_face_ave, one_plus_eta_edge_2_ave, ope_ave_df, ope2_ave

        implicit none

        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face
        real, dimension(2,nq,nface), intent(in) :: zbot_face
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(npoin_q,nlayers), intent(out) :: H_r
        real, dimension(2,nq,nface,nlayers), intent(out) :: H_r_face
        real, dimension(npoin_q,nlayers+1), intent(out) :: p
        real, dimension(npoin,nlayers+1), intent(out) :: z_elev
        real, dimension(2,npoin_q,nlayers+1), intent(out) :: grad_z

        ! local variables
        real, dimension(nlayers) :: alpha_over_g, g_over_alpha
        real, dimension(nq,nface,nlayers+1) :: z_edge_plus, z_edge_minus
        real, dimension(2,nq,nface,nlayers+1) :: p_face, z_face
        real, dimension(nq,nface,nlayers+1) :: p_edge_plus, p_edge_minus
        real, dimension(nq,nlayers) :: dp_plus, dp_minus, dp_temp
        real, dimension(npoin_q,nlayers) :: dpprime_H
        real, dimension(2,nq,nface,nlayers) :: dpprime_H_face
        integer :: iface, ilr, k, iquad, ktemp, I
        real :: z_intersect_top,z_intersect_bot, dz_intersect, H_r_plus, H_r_minus, acceleration
        real :: p_intersect_bot, p_intersect_top
        real, dimension(nq) :: one_plus_eta_edge, one_plus_eta_cell_face
        real, dimension(nlayers+1) :: p2l, p2r
        real :: H_corr,p_inc, weight, H_corr1,p_inc1, H_corr2,p_inc2, temp
        integer :: Iq, el, er, jl, kl, il, ilocr, ilocl, ir, jr, kr
        real :: ope_df(npoin), z_r(npoin_q,nlayers+1), z_temp(npoin)
        real, dimension(nq) :: one_plus_eta_edge2,one_plus_eta_edge1

        real, dimension(npoin_q,nlayers+1) :: zt
        real, dimension(2,npoin_q) :: z_xy

        real, dimension(2, 2,nq,nface,nlayers+1) :: com_face
        
        z_elev = 0.0
        zt = 0.0

        dpprime_H(:,:) = qprime(1,:,:)
        
        dpprime_H_face(:,:,:,:) = qprime_face(1,:,:,:,:)
        z_face = 0.0
        p_face = 0.0
        
        do k=1,nlayers
            alpha_over_g(k) = alpha_mlswe(k)/gravity
            g_over_alpha(k) = gravity/alpha_mlswe(k)
        enddo

        do iface = 1,nface

            !Store Left Side Variables
            one_plus_eta_cell_face = ope_face_ave(1,:,iface)
            p_face(1,:,iface,1) = 0.0

            do k=1,nlayers
                dp_temp(:,k) = one_plus_eta_cell_face * dpprime_H_face(1,:,iface,k)
                p_face(1,:,iface,k+1) = p_face(1,:,iface,k) + dp_temp(:,k)
            end do
            
            z_face(1,:,iface,nlayers+1) = zbot_face(1,:,iface)
            do k=nlayers,1,-1
                z_face(1,:,iface,k) = z_face(1,:,iface,k+1) + alpha_over_g(k) * dp_temp(:,k)
            end do

            !Store Right Side Variables
            one_plus_eta_cell_face = ope_face_ave(2,:,iface)
            p_face(2,:,iface,1) = 0.0

            do k=1,nlayers
                dp_temp(:,k) = one_plus_eta_cell_face * dpprime_H_face(2,:,iface,k)
                p_face(2,:,iface,k+1) = p_face(2,:,iface,k) + dp_temp(:,k)
            end do
            
            z_face(2,:,iface,nlayers+1) = zbot_face(2,:,iface)
            do k=nlayers,1,-1
                z_face(2,:,iface,k) = z_face(2,:,iface,k+1) + alpha_over_g(k) * dp_temp(:,k)
            end do 


            ! one_plus_eta_edge = sqrt(one_plus_eta_edge_2_ave(:,iface))
            one_plus_eta_edge = one_plus_eta_edge_2_ave(:,iface)

            do k = 1,nlayers
                dp_plus(:,k) = one_plus_eta_edge * dpprime_H_face(1,:,iface,k)
                dp_minus(:,k) = one_plus_eta_edge * dpprime_H_face(2,:,iface,k)
            end do

            z_edge_plus(:,iface,nlayers+1) = zbot_face(1,:,iface)
            z_edge_minus(:,iface,nlayers+1) = zbot_face(2,:,iface)

            do k = nlayers,1,-1
                z_edge_plus(:,iface,k) = z_edge_plus(:,iface,k+1) + alpha_over_g(k) * dp_plus(:,k)
                z_edge_minus(:,iface,k) = z_edge_minus(:,iface,k+1) + alpha_over_g(k) * dp_minus(:,k)
            end do

            p_edge_plus(:,iface,2) = dp_plus(:,1)
            p_edge_minus(:,iface,2) = dp_minus(:,1)

            do k = 2,nlayers
                p_edge_plus(:,iface,k+1) = p_edge_plus(:,iface,k) + dp_plus(:,k)
                p_edge_minus(:,iface,k+1) = p_edge_minus(:,iface,k) + dp_minus(:,k)
            end do

        end do  ! end loop over cell faces
        
        ! Find elevation at the top of layer k

        z_elev(:,nlayers+1) = zbot_df(:)
            
        do k = nlayers,1,-1
            z_elev(:,k) = z_elev(:,k+1) + alpha_over_g(k) * (qprime_df(1,:,k) * ope_ave_df(:))
        end do

        ! Compute the gradient of elevation z

        do k = 1, nlayers+1
            call compute_gradient_quad(grad_z(:,:,k),z_elev(:,k))
        end do

        ! Find pressure at the bottom of layer k and layer height
        p(:,1) = 0.0
        do k = 1,nlayers
            p(:,k+1) = p(:,k) + dpprime_H(:,k) * ope_ave(:)
            H_r(:,k) = 0.5*alpha_mlswe(k) * (p(:,k+1)**2 - p(:,k)**2)
        end do

        ! Compute H_r at the element face

        do iface = 1, nface

            ! Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq
                do k = 1, nlayers
        
                    ! ====== Begin computation of H_r for the right side ====
        
                    ! Computation from + side for layer k
                    H_r_plus = 0.5*alpha_mlswe(k)*(p_edge_plus(iquad, iface, k+1)**2 - p_edge_plus(iquad, iface, k)**2)
        
                    ! Computation from - side for layer k
                    H_r_minus = 0.0
    
                    do ktemp = 1, nlayers
        
                        z_intersect_top = min(z_edge_minus(iquad, iface, ktemp), z_edge_plus(iquad, iface, k))
                        z_intersect_bot = max(z_edge_minus(iquad, iface, ktemp+1), z_edge_plus(iquad, iface, k+1))
                        dz_intersect = z_intersect_top - z_intersect_bot
        
                        if (dz_intersect > 0.0) then
        
                            p_intersect_bot = p_edge_minus(iquad,iface,ktemp+1) - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_minus(iquad,iface,ktemp+1))
        
                            p_intersect_top = p_edge_minus(iquad,iface,ktemp+1) - g_over_alpha(ktemp)*(z_intersect_top - z_edge_minus(iquad,iface,ktemp+1))
        
                            H_r_minus = H_r_minus + 0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
        
                        end if
        
                    end do
        
                    H_r_face(1, iquad, iface, k) = 0.5*(H_r_plus + H_r_minus)
        
                    ! ====== End computation of H_r for the left side =====
        
                    ! ====== Begin computation of H_r for the right side ====
        
                    ! Computation from - side for layer k

                    H_r_minus = 0.5*alpha_mlswe(k)*(p_edge_minus(iquad, iface, k+1)**2 - p_edge_minus(iquad, iface, k)**2)
        
                    ! Computation from + side for layer k
                    H_r_plus = 0.0
        
                    do ktemp = 1, nlayers
        
                        z_intersect_top = min(z_edge_plus(iquad, iface, ktemp), z_edge_minus(iquad, iface, k))
                        z_intersect_bot = max(z_edge_plus(iquad, iface, ktemp+1), z_edge_minus(iquad, iface, k+1))
                        dz_intersect = z_intersect_top - z_intersect_bot
        
                        if (dz_intersect > 0.0) then
        
                            p_intersect_bot = p_edge_plus(iquad, iface, ktemp+1) - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_plus(iquad, iface, ktemp+1))
        
                            p_intersect_top = p_edge_plus(iquad, iface, ktemp+1) - g_over_alpha(ktemp)*(z_intersect_top - z_edge_plus(iquad, iface, ktemp+1))
        
                            H_r_plus = H_r_plus + 0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
        
                        end if
        
                    end do
        
                    H_r_face(2, iquad, iface, k) = 0.5*(H_r_plus + H_r_minus)
        
                    ! ===== End computation of H_r for the right side =====
        
                end do

                ! Wall Boundary conditions

                if(er == -4) then 

                    p2l = 0.0
                    p2r = 0.0

                    do k = 1,nlayers

                        p2l(k+1) = p_face(1,iquad,iface,k+1)
                        H_r_face(1,iquad,iface,k) = 0.5*alpha_mlswe(k)*(p2l(k+1)**2 - p2l(k)**2)

                        p2r(k+1) = p_face(2,iquad,iface,k+1)
                        H_r_face(2,iquad,iface,k) = 0.5*alpha_mlswe(k)*(p2r(k+1)**2 - p2r(k)**2)

                    end do

                end if
            end do
        end do

        ! Adjust the values of  H_r, at quadrature points, so that the vertical sum of
        ! H_r  over all layers equals the time average of the barotropic forcing  H  over all barotropic substeps of the baroclinic 
        ! time interval.

        if(adjust_H_vertical_sum == 1) then 
            do I = 1, npoin_q

                acceleration = 0.0
                if (pbprime(I) > 0.0) then
                    acceleration = (H_ave(I) - sum(H_r(I,:))) / pbprime(I)
                end if
                H_r(I,:) = H_r(I,:) + dpprime_H(I,:) * acceleration
            end do

        elseif(adjust_H_vertical_sum == 2) then 
            do I = 1, npoin_q
                weight = 1.0
                acceleration = sum(H_r(I,:))
                if(acceleration > 0.0) then 
                    weight = H_ave(I) / acceleration
                end if
                H_r(I,:) = H_r(I,:) * weight
            end do
        end if
        
        do iface = 1, nface
            el = face(7,iface)
            er = face(8,iface)
            do iquad = 1, nq
                if(er /= -4) then
                    do k = 1, nlayers-1          ! interface at the bottom of layer k
            
                        ! Corrections at the left side of a face.
                        p_inc1 = g_over_alpha(k)*(z_face(1, iquad, iface, k+1) - z_edge_plus(iquad, iface, k+1))
                        H_corr1 = 0.5 * alpha_mlswe(k) * ((p_face(1, iquad, iface, k+1) + p_inc1)**2 - p_face(1, iquad, iface, k+1)**2)
                        H_r_face(1, iquad, iface, k) = H_r_face(1, iquad, iface, k) - H_corr1
                        H_r_face(1, iquad, iface, k+1) = H_r_face(1, iquad, iface, k+1) + H_corr1
            
                        ! Corrections at the right side of a face.
                        p_inc2 = g_over_alpha(k)*(z_face(2, iquad, iface, k+1) - z_edge_minus(iquad, iface, k+1))
                        H_corr2 = 0.5 * alpha_mlswe(k) * ((p_face(2, iquad, iface, k+1) + p_inc2)**2 - p_face(2, iquad, iface, k+1)**2)
                        H_r_face(2, iquad, iface, k) = H_r_face(2, iquad, iface, k) - H_corr2
                        H_r_face(2, iquad, iface, k+1) = H_r_face(2, iquad, iface, k+1) + H_corr2
            
                    end do
                end if
            
                ! Adjust the vertical sums of  H_r,  for the sake of consistency
                ! between the layer equations and the barotropic equations.

                ! Adjust the values of  H_r, at element faces, so that the vertical sum of
                ! H_r  over all layers equals the time average of the barotropic forcing  H  over all barotropic substeps of the baroclinic 
                ! time interval.

                ! The difference between the time-averaged  H  and the vertical sum of
                ! H_r  must be distributed over the layers via some sort of 
                ! weighting scheme. We conider the following approaches:

                ! (1) Weight according to layer thickness.  
                !       That is, the weight for layer r is  
                !       (Delta p)_r / (sum of (Delta p)_s over all layers s) =  (Delta p)_r / p_b
                !                                                          =  (Delta p')_r / p'_b
                !       The adjusted  H_r  is then
                !       (H_r)_adjusted = H_r + [(Delta p')_r / p'_b] * [ H_ave - sum(H_s) ] 
                !
                ! (2) Weight according to the current value of H_r.
                !       That is, the weight for layer r is  
                !       H_r / (sum of H_s over all layers s).  
                !       The adjusted  H_r  is then
                !       (H_r)_adjusted  =   H_r + [ H_r/(sum H_s)] * [ H_ave - sum(H_s) ]
                !                       =   H_r +   H_r * H_ave/(sum H_s)  -  H_r
                !                       =   H_r * H_ave/(sum H_s)
                !       Therefore, at each quadrature point and cell edge, multiply
                !       the current value of  H_r  by the layer-independent ratio
                !       H_ave/(sum H_s),  which should be approximately equal to  1.  

                if(adjust_H_vertical_sum == 1) then 
                    ! Left side of face
                    acceleration = 0.0
                    if (pbprime_face(1,iquad,iface) > 0.0) then
                        acceleration = (H_face_ave(iquad,iface) - sum(H_r_face(1,iquad,iface,:))) / pbprime_face(1,iquad,iface)
                    end if
                    H_r_face(1,iquad,iface,:) = H_r_face(1,iquad,iface,:) + dpprime_H_face(1,iquad,iface,:)*acceleration

                    ! Right side of face
                    acceleration = 0.0
                    if (pbprime_face(2,iquad,iface) > 0.0) then
                        acceleration = (H_face_ave(iquad,iface) - sum(H_r_face(2,iquad,iface,:))) / pbprime_face(2,iquad,iface)
                    end if
                    H_r_face(2,iquad,iface,:) = H_r_face(2,iquad,iface,:) + dpprime_H_face(2,iquad,iface,:)*acceleration

                elseif(adjust_H_vertical_sum == 2) then
                    ! Left side of face
                    weight = 1.0
                    acceleration = sum(H_r_face(1,iquad,iface,:))
                    if(acceleration > 0.0) then 
                        weight = H_face_ave(iquad,iface) / acceleration
                    end if
                    H_r_face(1,iquad,iface,:) = H_r_face(1,iquad,iface,:) * weight
                    ! Right side of face
                    weight = 1.0
                    acceleration = sum(H_r_face(2,iquad,iface,:))
                    if(acceleration > 0.0) then 
                        weight = H_face_ave(iquad,iface) / acceleration
                    end if
                    H_r_face(2,iquad,iface,:) = H_r_face(2,iquad,iface,:) * weight
                end if

            end do
        end do

    end subroutine layer_pressure_terms

    subroutine layer_windbot_stress_terms(tau_wind_int, tau_bot_int,qprime)

        use mod_initial, only : pbprime, tau_wind, alpha_mlswe
        use mod_grid, only :  npoin_q
        use mod_input, only: nlayers, dp_tau_wind, dp_tau_bot
        use mod_mpi_utilities, only : irank
        use mod_constants, only : gravity
        use mod_variables, only: tau_bot_ave, tau_wind_ave

        implicit none
    
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime 
        real, dimension(2,npoin_q,nlayers), intent(out) :: tau_wind_int, tau_bot_int ! output arrays
    
        real, dimension(npoin_q,nlayers+1) :: pprime_temp 
    
        integer :: k, i, Iq ! loop indices
        real, dimension(npoin_q) :: temp ! temporary variable
        character :: num
        real :: ptop, pbot, htop, Pstress, temp1

        pprime_temp = 0.0
        tau_wind_int = 0.0
        tau_bot_int = 0.0
        temp = 0.0

        do k = 1, nlayers
            pprime_temp(:, k+1) = pprime_temp(:, k) + qprime(1,:,k)
        end do

        ! write(num, '(i1)') irank
        ! open(unit=100, file='tau_wind_int'//trim(num)//'.dat', status='unknown')
        
        Pstress = (gravity/alpha_mlswe(1)) * 50.0 ! pressure corresponding to 50m depth at which wind stress is reduced to 0
    
        do k = 1, nlayers
            do Iq = 1,npoin_q 
                temp1 = (min(pprime_temp(Iq,k+1), Pstress) - min(pprime_temp(Iq,k), Pstress))/ Pstress
                tau_wind_int(1,Iq,k) = temp1*tau_wind(1,Iq)
                tau_wind_int(2,Iq,k) = temp1*tau_wind(2,Iq)
            end do
        end do

        tau_bot_int(1,:,nlayers) = tau_bot_ave(1,:)
        tau_bot_int(2,:,nlayers) = tau_bot_ave(2,:)
    
    end subroutine layer_windbot_stress_terms

    subroutine velocity_df(q_df, qb_df, flag_pred)

        use mod_grid, only : npoin, intma, face, nface 
        use mod_basis, only : ngl
        use mod_input, only: nlayers, dp_cutoff1, dp_cutoff2, icase, ifilter
        use mod_face, only: imapl

        implicit none
        
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(4,npoin), intent(in) :: qb_df
        integer, intent(in) :: flag_pred
        
        real, dimension(nlayers) :: a, b, c
        real, dimension(nlayers,2) :: r
        real, dimension(nlayers) :: weight
        real :: eps = 1.0e-12
        real :: ubar, vbar, mult, dpbar
        integer :: I, k, el, er, n , iface, il, jl, kl, Iq
        real :: uv_df(2,npoin,nlayers)

        uv_df = 0.0

        do k = 1,nlayers
            uv_df(1,:,k) = q_df(2,:,k) / q_df(1,:,k)
            uv_df(2,:,k) = q_df(3,:,k) / q_df(1,:,k)
        end do

        do I = 1, npoin
            ubar = 0.0
            vbar = 0.0
            
            do k = 1, nlayers
                ubar = ubar + uv_df(1,I,k) * q_df(1,I,k)
                vbar = vbar + uv_df(2,I,k) * q_df(1,I,k)
            end do

            if(qb_df(1,I) > 0.0) then
                
                ubar = ubar / qb_df(1,I)
                vbar = vbar / qb_df(1,I)

                do k = 1, nlayers

                    uv_df(1,I,k) = uv_df(1,I,k) - ubar + qb_df(3,I)/qb_df(1,I)
                    uv_df(2,I,k) = uv_df(2,I,k) - vbar + qb_df(4,I)/qb_df(1,I)
                end do

            else
                uv_df(:,I,:) = 0.0
            end if
        end do

        if(flag_pred == 0 .and. ifilter > 0) then 
            do k = 1,nlayers
                call filter_mlswe(uv_df(:,:,k),2)
            end do
            call layer_mom_boundary_df(uv_df(:,:,:))
        end if

        do k = 1,nlayers
            q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
            q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
        end do

    end subroutine velocity_df

    subroutine velocity(uv,q, qb)

        use mod_grid, only : npoin_q
        use mod_input, only: nlayers

        implicit none
        
        real, dimension(3,npoin_q,nlayers), intent(in) :: q
        real, dimension(4,npoin_q), intent(in) :: qb

        real, dimension(2,npoin_q,nlayers), intent(out) :: uv
    
        real :: eps = 1.0e-4
        real :: ubar, vbar, weight
        integer :: Iq, k

        do k = 1,nlayers
            uv(1,:,k) = q(2,:,k) / q(1,:,k)
            uv(2,:,k) = q(3,:,k) / q(1,:,k)
        end do

        do Iq = 1, npoin_q
            ubar = 0.0
            vbar = 0.0
            
            do k = 1, nlayers
                ubar = ubar + uv(1,Iq,k) * q(1,Iq,k)
                vbar = vbar + uv(2,Iq,k) * q(1,Iq,k)
            end do

            if(qb(1,Iq) > 0.0) then
                
                ubar = ubar / qb(1,Iq)
                vbar = vbar / qb(1,Iq)

                do k = 1, nlayers
                    uv(1,Iq,k) = uv(1,Iq,k) - ubar + qb(3,Iq)/qb(1,Iq)
                    uv(2,Iq,k) = uv(2,Iq,k) - vbar + qb(4,Iq)/qb(1,Iq)
                end do

            else
                uv(:,Iq,:) = 0.0
            end if
        end do

    end subroutine velocity

    subroutine interpolate_mom(q_df,q, qb,flag_pred)

        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad
        use mod_input, only: nlayers, ifilter
        use mod_create_rhs_mlswe, only: interpolate_layer_from_quad_to_node

        use mod_initial, only: psih, indexq

        implicit none
      
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(3,npoin_q,nlayers), intent(inout) :: q
        real, dimension(4,npoin_q), intent(in) :: qb
        integer, intent(in) :: flag_pred

        integer :: k, e, kquad,jquad, iquad, Iq, m, n, l,I, ip
        real :: hi
        real, dimension(2,npoin_q,nlayers) :: uv
        
        ! Interpolate momentum dofs to quad points
        call evaluate_mom(q,q_df)

        ! Extract velocity 
        call velocity(uv,q, qb)

        ! Compute the new mometum value based on the velocity and layer pressure 
        do k = 1,nlayers
            q(2,:,k) = uv(1,:,k) * q(1,:,k)
            q(3,:,k) = uv(2,:,k) * q(1,:,k)
        end do 

        ! Interpolate from quad to nodal space

        call interpolate_layer_from_quad_to_node(q_df, q)

    end subroutine interpolate_mom

    subroutine evaluate_mom(q,q_df)

        ! Interpolate from dofs to quadrature points 

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad
        use mod_input, only: nlayers

        use mod_initial, only: psih, indexq

        implicit none

        real, intent(inout) :: q(3,npoin_q,nlayers)
        real, intent(in) :: q_df(3,npoin,nlayers)

        integer :: k, Iq, m, n, l,I, ip
        real :: hi

        q(2:3,:,:) = 0.0

        do Iq = 1,npoin_q
            do ip = 1,npts

                I = indexq(Iq,ip)
                hi = psih(Iq,ip)
                
                q(2,Iq,:) = q(2,Iq,:) + q_df(2,I,:)*hi
                q(3,Iq,:) = q(3,Iq,:) + q_df(3,I,:)*hi
            end do
        end do

    end subroutine evaluate_mom
    
    subroutine evaluate_mom_face_all(q_face, uv_face, q, uv)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq
        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq, nface,face
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q, normal_vector_q

        implicit none

        real, dimension(3, 2, nq, nface, nlayers), intent(inout) :: q_face
        real, dimension(2, 2, nq, nface, nlayers), intent(out) :: uv_face
        real, dimension(3, npoin_q, nlayers), intent(in) :: q
        real, dimension(2, npoin_q, nlayers), intent(in) :: uv
    
        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, kl, kr, jquad
        real :: nx, ny, un(nlayers)
    
        q_face(2:3,:,:,:,:) = 0.0
        uv_face = 0.0
    
        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                il = imapl_q(1,iquad,1,iface)
                jl = imapl_q(2,iquad,1,iface)
                kl = imapl_q(3,iquad,1,iface)

                I = intma_dg_quad(il,jl,kl,el)

                q_face(2:3,1,iquad,iface,:) = q(2:3,I,:)
                uv_face(:,1,iquad,iface,:) = uv(:,I,:)

                if(er > 0) then

                    ir = imapr_q(1,iquad,1,iface)
                    jr = imapr_q(2,iquad,1,iface)
                    kr = imapr_q(3,iquad,1,iface)

                    I = intma_dg_quad(ir,jr,kr,er)

                    q_face(2:3,2,iquad,iface,:) = q(2:3,I,:)
                    uv_face(:,2,iquad,iface,:) = uv(:,I,:)

                else 
                    q_face(2:3,2,iquad,iface,:) = q_face(2:3,1,iquad,iface,:)
                    uv_face(:,2,iquad,iface,:) = uv_face(:,1,iquad,iface,:)
                    
                    if(er == -4) then

                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)

                        un = q(2,I,:)*nx + q(3,I,:)*ny

                        q_face(2,2,iquad,iface,:) = q(2,I,:) - 2.0*un*nx
                        q_face(3,2,iquad,iface,:) = q(3,I,:) - 2.0*un*ny

                        un = uv(1,I,:)*nx + uv(2,I,:)*ny

                        uv_face(1,2,iquad,iface,:) = uv(1,I,:) - 2.0*un*nx
                        uv_face(2,2,iquad,iface,:) = uv(2,I,:) - 2.0*un*ny

                    elseif(er == -2) then 
                        q_face(2:3,2,iquad,iface,:) = -q_face(2:3,1,iquad,iface,:)

                        uv_face(:,2,iquad,iface,:) = -uv_face(:,1,iquad,iface,:)
                        
                    end if
                end if

            end do
        end do
    
    end subroutine evaluate_mom_face_all

    subroutine evaluate_mom_face(q_face, q)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq
        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq, nface,face
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q, normal_vector_q

        implicit none

        real, dimension(3, 2, nq, nface, nlayers), intent(inout) :: q_face
        real, dimension(3, npoin_q, nlayers), intent(in) :: q
    
        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, kl, kr, jquad
        real :: nx, ny, un(nlayers)
    
        q_face(2:3,:,:,:,:) = 0.0
    
        ! do k = 1, nlayers
        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                il = imapl_q(1,iquad,1,iface)
                jl = imapl_q(2,iquad,1,iface)
                kl = imapl_q(3,iquad,1,iface)

                I = intma_dg_quad(il,jl,kl,el)

                q_face(2:3,1,iquad,iface,:) = q(2:3,I,:)

                if(er > 0) then

                    ir = imapr_q(1,iquad,1,iface)
                    jr = imapr_q(2,iquad,1,iface)
                    kr = imapr_q(3,iquad,1,iface)

                    I = intma_dg_quad(ir,jr,kr,er)

                    q_face(2:3,2,iquad,iface,:) = q(2:3,I,:)

                else 
                    q_face(2:3,2,iquad,iface,:) = q_face(2:3,1,iquad,iface,:)
                    if(er == -4) then

                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)

                        un = q(2,I,:)*nx + q(3,I,:)*ny

                        q_face(2,2,iquad,iface,:) = q(2,I,:) - 2.0*un*nx
                        q_face(3,2,iquad,iface,:) = q(3,I,:) - 2.0*un*ny
                    elseif(er == -2) then 
                        q_face(2:3,2,iquad,iface,:) = -q_face(2:3,1,iquad,iface,:)
                        
                    end if
                end if

            end do
        end do
        ! end do
    
    end subroutine evaluate_mom_face

    subroutine velocity_face(uv_face, uv)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq
        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq, nface,face
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q, normal_vector_q

        implicit none

        real, dimension(2, 2, nq, nface, nlayers), intent(out) :: uv_face
        real, dimension(2, npoin_q, nlayers), intent(in) :: uv
    
        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, kl, kr, jquad
        real :: nx, ny, un(nlayers)
    
        uv_face = 0.0
    
        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                il = imapl_q(1,iquad,1,iface)
                jl = imapl_q(2,iquad,1,iface)
                kl = imapl_q(3,iquad,1,iface)

                I = intma_dg_quad(il,jl,kl,el)

                uv_face(:,1,iquad,iface,:) = uv(:,I,:)

                if(er > 0) then

                    ir = imapr_q(1,iquad,1,iface)
                    jr = imapr_q(2,iquad,1,iface)
                    kr = imapr_q(3,iquad,1,iface)

                    I = intma_dg_quad(ir,jr,kr,er)

                    uv_face(:,2,iquad,iface,:) = uv(:,I,:)

                else 
                    uv_face(:,2,iquad,iface,:) = uv_face(:,1,iquad,iface,:)
                    if(er == -4) then

                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)

                        un = uv(1,I,:)*nx + uv(2,I,:)*ny

                        uv_face(1,2,iquad,iface,:) = uv(1,I,:) - 2.0*un*nx
                        uv_face(2,2,iquad,iface,:) = uv(2,I,:) - 2.0*un*ny
                    elseif(er == -2) then 
                        uv_face(:,2,iquad,iface,:) = -uv_face(:,1,iquad,iface,:)
                        
                    end if
                end if

            end do
        end do
    
    end subroutine velocity_face

    subroutine layer_mom_boundary_df(q)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl,nq
        use mod_grid, only:  npoin, intma, nface, face,mod_grid_get_face_nq
        use mod_initial, only: pbprime_face
        use mod_face, only: imapl, imapr, normal_vector
        use mod_input, only: nlayers, mlswe_bc_strong

        implicit none

        real, intent(inout) :: q(2,npoin,nlayers)

        integer :: iface, n, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr,k
        real :: nx, ny, upnl(nlayers)


        do iface = 1, nface                 

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            if(er == -4) then
        
                do n = 1, ngl

                    il=imapl(1,n,1,iface)
                    jl=imapl(2,n,1,iface)
                    kl=imapl(3,n,1,iface)
                    I=intma(il,jl,kl,el)

                    nx = normal_vector(1,n,1,iface)
                    ny = normal_vector(2,n,1,iface)

                    upnl = q(1,I,:)*nx + q(2,I,:)*ny

                    q(1,I,:) = q(1,I,:) - upnl*nx
                    q(2,I,:) = q(2,I,:) - upnl*ny

                end do
            end if
    
        end do
      
    end subroutine layer_mom_boundary_df 

    subroutine layer_mom_boundary(q)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl,nq
        use mod_grid, only:  npoin_q, intma_dg_quad, nface, face,mod_grid_get_face_nq
        use mod_initial, only: pbprime_face
        use mod_face, only: imapl_q, imapr_q, normal_vector_q
        use mod_input, only: nlayers

        implicit none

        real, intent(inout) :: q(3,npoin_q,nlayers)

        integer :: iface, n, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr,k
        real :: nx, ny, upnl(nlayers)

      
        ! do k = 1,nlayers
        do iface = 1, nface                 

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            if(er == -4) then
        
                do n = 1, nq

                    il=imapl_q(1,n,1,iface)
                    jl=imapl_q(2,n,1,iface)
                    kl=imapl_q(3,n,1,iface)
                    I=intma_dg_quad(il,jl,kl,el)

                    nx = normal_vector_q(1,n,1,iface)
                    ny = normal_vector_q(2,n,1,iface)

                    upnl = q(2,I,:)*nx + q(3,I,:)*ny

                    q(2,I,:) = q(2,I,:) - upnl*nx
                    q(3,I,:) = q(3,I,:) - upnl*ny

                end do
            elseif(er == -2) then
        
                do n = 1, nq

                    il=imapl_q(1,n,1,iface)
                    jl=imapl_q(2,n,1,iface)
                    kl=imapl_q(3,n,1,iface)
                    I=intma_dg_quad(il,jl,kl,el)

                    q(2,I,:) = 0.0
                    q(3,I,:) = 0.0

                end do
            end if
    
        end do

        ! end do
      
    end subroutine layer_mom_boundary

    
    subroutine shear_stress_system(uv,q)

        use mod_constants, only : gravity
        use mod_initial, only : alpha_mlswe, coriolis_quad
        use mod_grid, only : nface, npoin, npoin_q
        use mod_input, only : dt, ad_mlswe, nlayers, max_shear_dz
  
        implicit none
        
        real, dimension(3,npoin_q,nlayers), intent(in) :: q
        real, dimension(2,npoin_q,nlayers), intent(out) :: uv
        real, dimension(npoin_q) :: coeff
        real, dimension(nlayers) :: a,b,c
        real, dimension(nlayers,2) :: r
        real :: mult
        integer :: Iq,k

        do Iq = 1,npoin_q
            coeff(Iq) = gravity*dt*max(sqrt(0.5*coriolis_quad(Iq)*ad_mlswe)/alpha_mlswe(1), ad_mlswe/(alpha_mlswe(1) * max_shear_dz))
        end do
        
        do Iq = 1, npoin_q
            do k = 1, nlayers
                a(k) = -coeff(Iq)
                b(k) = q(1,Iq,k) + 2.0*coeff(Iq)
                c(k) = -coeff(Iq)
                r(k,1) = q(2,Iq,k)/q(1,Iq,k)
                r(k,2) = q(3,Iq,k)/q(1,Iq,k)
            end do
            
            b(1) = q(1,Iq,1) + coeff(Iq)
            b(nlayers) = q(1,Iq,nlayers) + coeff(Iq)
            a(1) = 0.0
            c(nlayers) = 0.0
            
            do k = 2, nlayers
                mult = a(k) / b(k-1)
                b(k) = b(k) - mult*c(k-1)
                r(k,1) = r(k,1) - mult*r(k-1,1)
                r(k,2) = r(k,2) - mult*r(k-1,2)
            end do
            
            r(nlayers,1) = r(nlayers,1) / b(nlayers)
            r(nlayers,2) = r(nlayers,2) / b(nlayers)
            uv(1,Iq,nlayers) = r(nlayers,1)
            uv(2,Iq,nlayers) = r(nlayers,2)
            
            do k = nlayers-1,1,-1
                r(k,1) = (r(k,1) - c(k)*r(k+1,1)) / b(k)
                r(k,2) = (r(k,2) - c(k)*r(k+1,2)) / b(k)
                uv(1,Iq,k) = r(k,1)
                uv(2,Iq,k) = r(k,2)
            end do
        end do
        
    end subroutine shear_stress_system

    subroutine filter_mlswe(q,nvarb)

        use mod_basis, only: nglx, ngly, nglz, npts, ngl, f, fx, fy, fz

        use mod_filter, only: b, b_data

        use mod_grid, only: intma, npoin, nelem

        use mod_initial, only: nvar, nvart

        use mod_gradient, only: compute_local_gradient_filter_v3

        use mod_metrics, only: jac, massinv
    
        implicit none
    
        !global arrays
        integer, intent(in) :: nvarb
        real, intent(inout) :: q(nvarb,npoin)
    
        !local arrays
        real, dimension(nvarb,npts) :: qq, fqf, qq_i, qq_ij, qq_ijk
        real, dimension(npts) :: jac_e
        real :: wq
        integer :: inode(npts)
        integer :: i, j, k, l, m, e, ip, ip1, ii, ie
        integer :: ndim, ndim2
        real r_k, u_k, v_k, w_k, t_k
    
        !Store dimensions of MxM object
        ndim  = nvarb
        ndim2 = nvarb
        
        !Initialize
        b=0.0

        !loop thru the elements
        do e=1,nelem
        
            !Store Element Variables
            ii=0
            do j=1,ngly
                do i=1,nglx
                    ip=intma(i,j,1,e)
                    ii=ii+1
                    inode(ii)=ip
                    do m=1,nvarb
                        qq(m,ii) = q(m,ip)
                    end do
                    jac_e(ii)=jac(i,j,1,e)
                end do !i
            end do !j
    
            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_filter_v3(fqf,qq,nglx,ngly,nglz,ndim)

            !Do Numerical Integration
            do i=1,npts
                ip=inode(i)

                !Gauss-Lobatto Weight and Jacobian
                wq=jac_e(i)
                
                !Store RHS
                do m=1,ndim
                    b(m,ip)=b(m,ip) + wq*fqf(m,i)
                end do !m
            
            end do !i
    
        end do !e
    
        !Apply Mass Matrix Inverse

        do m=1,ndim
            b(m,:) = b(m,:)*massinv(:)
        end do !m
        !q=b
        q(1:ndim2,:) = b(1:ndim2,:)
    
    end subroutine filter_mlswe
    
end module mod_layer_terms