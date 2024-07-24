! ===========================================================================================================================
! This module contains the routines for the barotropic flux terms
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: March 27, 2023
! ==========================================================================================================================

module mod_barotropic_terms

    
    use mod_grid, only: npoin_q, nface, intma_dg_quad
    use mod_basis, only: nqx, nqy, nqz, nq
    use mod_input, only: nlayers
        
    implicit none

    public :: btp_evaluate_mom, &
                btp_evaluate_mom_face, btp_evaluate_pb, btp_evaluate_pb_face, btp_mass_advection_terms, &
                btp_bcl_coeffs, btp_mom_boundary_df, evaluate_quprime, &
                compute_btp_terms, btp_evaluate_mom_dp_face, btp_evaluate_mom_dp, &
                evaluate_quprime2, evaluate_visc_terms, restart_mlswe_variales, restart_mlswe, &
                compute_gradient_uv, compute_btp_mom_terms, compute_gradient_uv_v1, btp_bcl_grad_coeffs, &
                compute_btp_volume_terms, compute_btp_face_terms, btp_interpolate_avg, btp_extract_df, &
                btp_interpolate_avg_v1, btp_extract_face, btp_extract_face_v1

    contains

    subroutine compute_btp_terms(qb,qb_face, qprime, qb_df)

        ! This routine computes the flux terms for barotropic system

        use mod_grid, only: npoin_q, nface, npoin, face, intma_dg_quad
        use mod_basis, only: nq
        use mod_face, only: imapl_q, imapr_q, normal_vector_q
        use mod_input, only: nlayers, cd_mlswe, botfr
        use mod_initial, only: alpha_mlswe
        use mod_initial, only: coeff_pbpert_L, coeff_pbub_LR, coeff_pbpert_R, one_over_pbprime_edge, &
                                one_over_pbprime, one_over_pbprime_face, one_over_pbprime_df, &
                                coeff_mass_pbub_L, coeff_mass_pbub_R, coeff_mass_pbpert_LR
        use mod_constants, only: gravity

        use mod_variables, only: Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, H_bcl_edge
        use mod_variables, only: Qu_face, Qv_face, one_plus_eta_face, flux_edge, Quu, Qvv, Quv, H, one_plus_eta, &
                                    one_plus_eta_df, tau_bot, btp_mass_flux, H_face, one_plus_eta_edge_2

        implicit none

        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(4,2,nq,nface), intent(in) :: qb_face
        real, dimension(4,npoin_q), intent(in) :: qb

        integer :: iface, iquad, Iq
        real, dimension(nq) :: ul, ur, vl, vr, upl, upr, vpl, vpr
        real :: pbub_left, pbvb_left, pbpert_left
        real :: pbub_right, pbvb_right, pbpert_right
        real :: pU_L, pU_R, nxl, nyl, nxr, nyr, rho
        real :: pbpert_edge, eta_edge, speed1
        integer :: ir, jr, kr, kl, ier, jquad, el, er, jl, il
        real, dimension(npoin_q) :: ubot, vbot, speed
        real, dimension(nq,nface) :: one_plus_eta_edge
        real, dimension(npoin_q) :: ub, vb 

        ub = qb(3,:)/qb(1,:)
        vb = qb(4,:)/qb(1,:)

        tau_bot = 0.0

        ! Bottom friction terms

        if(botfr == 1) then 

            ubot = qprime(1,:,nlayers)*(qprime(2,:,nlayers) + ub)
            vbot = qprime(1,:,nlayers)*(qprime(3,:,nlayers) + vb)

            speed1 = cd_mlswe/gravity
        
            tau_bot(1,:) = speed1 * ubot
            tau_bot(2,:) = speed1 * vbot
        elseif(botfr == 2) then 

            rho = 1.0/alpha_mlswe(nlayers)

            ubot = qprime(2,:,nlayers) + ub
            vbot = qprime(3,:,nlayers) + vb
            speed = rho*cd_mlswe*sqrt(ubot**2 + vbot**2)
            
            tau_bot(1,:) = speed * ubot
            tau_bot(2,:) = speed * vbot
        end if 
        ! Compute (1 + eta) and (1 + eta)**2 quadrature points of each grid cell.

        one_plus_eta(:) = 1.0 + qb(2,:) * one_over_pbprime(:)

        one_plus_eta_df(:) = 1.0 + qb_df(2,:) * one_over_pbprime_df(:)

        btp_mass_flux = qb(3:4,:)

        ! Compute pressure forcing H

        H(:) = (one_plus_eta(:)**2) * H_bcl(:)

        ! Compute  Q_u  and  Q_v  at the quadrature points in each cell.
        
        Quu(:) = ub * qb(3,:) + one_plus_eta(:) * Q_uu_dp(:)
        Quv(:) = ub * qb(4,:) + one_plus_eta(:) * Q_uv_dp(:)
        Qvv(:) = vb * qb(4,:) + one_plus_eta(:) * Q_vv_dp(:)

        ! Compute eta

        do iface = 1, nface

            el = face(7,iface)
            er = face(8,iface)

            ! Use a Riemann problem to determine interpolated values of 
            ! pbpert = p'_b*eta at element faces.
            ! Here, eta iface the relative perturbation in bottom pressure, 
            ! relative to the global rest state that iface used to define the 
            ! barotropic-baroclinic splitting.

            ! Regard the interpolation process as Lax-Friedrichs interpolation,
            ! or a generalization thereof.
            ! Then compute the relative perturbations in bottom pressure 
            ! at each cell edge.  That is, 
            ! eta_edge(i) = pbpert_edge i / (p'_b at edge i)
            ! Also compute  (1 + eta_edge)**2 . 

            ! Compute pbpert_edge, eta_edge, and one_plus_eta_edge_2

            do iquad = 1, nq
                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)
        
                nxr = -nxl
                nyr = -nyl
        
                pbub_left = qb_face(3, 1, iquad, iface)
                pbvb_left = qb_face(4, 1, iquad, iface)
                pbpert_left = qb_face(2, 1, iquad, iface)
        
                pbub_right = qb_face(3, 2, iquad, iface)
                pbvb_right = qb_face(4, 2, iquad, iface)
                pbpert_right = qb_face(2, 2, iquad, iface)
        
                pU_L = nxl * pbub_left + nyl * pbvb_left
                pU_R = nxr * pbub_right + nyr * pbvb_right
        
                pbpert_edge = coeff_pbpert_L(iquad, iface) * pbpert_left &
                                            + coeff_pbpert_R(iquad, iface) * pbpert_right &
                                            + coeff_pbub_LR(iquad, iface) * (pU_L + pU_R)
        
                eta_edge = pbpert_edge * one_over_pbprime_edge(iquad, iface)
                one_plus_eta_edge(iquad, iface) = 1.0 + eta_edge
                one_plus_eta_edge_2(iquad, iface) = (1.0 + eta_edge)

                ! Compute mass fluxes at each element face.

                flux_edge(1,iquad,iface) = coeff_mass_pbub_L(iquad,iface) * pbub_left &
                                            + coeff_mass_pbub_R(iquad,iface) * pbub_right &
                                            + coeff_mass_pbpert_LR(iquad,iface) * (nxl * pbpert_left + nxr * pbpert_right)

                flux_edge(2,iquad,iface) = coeff_mass_pbub_L(iquad,iface) * pbvb_left &
                                            + coeff_mass_pbub_R(iquad,iface) * pbvb_right &
                                            + coeff_mass_pbpert_LR(iquad,iface) * (nyl * pbpert_left + nyr * pbpert_right)
            end do

            ! Compute (1 + eta) and (1 + eta)**2 at each element faces( 1 = left & 2=right side).

            one_plus_eta_face(1,:,iface) = 1.0 + qb_face(2,1,:,iface) * one_over_pbprime_face(1,:,iface)
            one_plus_eta_face(2,:,iface) = 1.0 + qb_face(2,2,:,iface) * one_over_pbprime_face(2,:,iface)

            ! Compute pressure forcing H_face at each element face.

            H_face(:,iface) = (one_plus_eta_edge_2(:,iface)**2) * H_bcl_edge(:,iface)

            if(er == -4) then 

                do iquad = 1, nq

                    il = imapl_q(1,iquad,1,iface)
                    jl = imapl_q(2,iquad,1,iface)
                    kl = imapl_q(3,iquad,1,iface)

                    Iq = intma_dg_quad(il,jl,kl,el)

                    H_face(iquad,iface) = H(Iq)
                    
                end do
            end if

            ! Compute momentum flux terms at each element face 

            ul = qb_face(3,1,:,iface)/qb_face(1,1,:,iface); ur = qb_face(3,2,:,iface)/qb_face(1,2,:,iface)
            vl = qb_face(4,1,:,iface)/qb_face(1,1,:,iface); vr = qb_face(4,2,:,iface)/qb_face(1,2,:,iface)

            upl = qb_face(3, 1, :, iface); upr = qb_face(3, 2, :, iface)
            vpl = qb_face(4, 1, :, iface); vpr = qb_face(4, 2, :, iface)

            Qu_face(1,:,iface) = 0.5*(ul*upl + ur*upr) + one_plus_eta_edge(:,iface) * Q_uu_dp_edge(:, iface)
            Qu_face(2,:,iface) = 0.5*(vl*upl + vr*upr) + one_plus_eta_edge(:,iface) * Q_uv_dp_edge(:, iface)

            Qv_face(1,:,iface) = 0.5*(ul*vpl + ur*vpr) + one_plus_eta_edge(:,iface) * Q_uv_dp_edge(:, iface)
            Qv_face(2,:,iface) = 0.5*(vl*vpl + vr*vpr) + one_plus_eta_edge(:,iface) * Q_vv_dp_edge(:, iface)

        end do

    end subroutine compute_btp_terms

    subroutine compute_btp_volume_terms(qb, qprime, qb_df)

        ! This routine computes the flux terms for barotropic system

        use mod_grid, only: npoin_q, nface, npoin, face, intma_dg_quad
        use mod_basis, only: nq
        use mod_face, only: imapl_q, imapr_q, normal_vector_q
        use mod_input, only: nlayers, cd_mlswe, botfr
        use mod_initial, only: alpha_mlswe
        use mod_initial, only: one_over_pbprime, one_over_pbprime_df

        use mod_constants, only: gravity

        use mod_variables, only: Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl
        use mod_variables, only: Quu, Qvv, Quv, H, one_plus_eta, tau_bot, btp_mass_flux, one_plus_eta_df

        implicit none

        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(4,npoin_q), intent(in) :: qb

        integer :: iface, iquad, Iq
        real, dimension(nq) :: ul, ur, vl, vr, upl, upr, vpl, vpr
        real :: pbub_left, pbvb_left, pbpert_left
        real :: pbub_right, pbvb_right, pbpert_right
        real :: pU_L, pU_R, nxl, nyl, nxr, nyr, rho
        real :: pbpert_edge, eta_edge, speed1
        integer :: ir, jr, kr, kl, ier, jquad, el, er, jl, il
        real, dimension(npoin_q) :: ubot, vbot, speed
        real, dimension(nq,nface) :: one_plus_eta_edge
        real, dimension(npoin_q) :: ub, vb 

        ub = qb(3,:)/qb(1,:)
        vb = qb(4,:)/qb(1,:)

        tau_bot = 0.0

        ! Bottom friction terms

        if(botfr == 1) then 

            ubot = qprime(1,:,nlayers)*(qprime(2,:,nlayers) + ub)
            vbot = qprime(1,:,nlayers)*(qprime(3,:,nlayers) + vb)

            speed1 = cd_mlswe/gravity
        
            tau_bot(1,:) = speed1 * ubot
            tau_bot(2,:) = speed1 * vbot
        elseif(botfr == 2) then 

            rho = 1.0/alpha_mlswe(nlayers)

            ubot = qprime(2,:,nlayers) + ub
            vbot = qprime(3,:,nlayers) + vb
            speed = rho*cd_mlswe*sqrt(ubot**2 + vbot**2)
            
            tau_bot(1,:) = speed * ubot
            tau_bot(2,:) = speed * vbot
        end if 
        ! Compute (1 + eta) and (1 + eta)**2 quadrature points of each grid cell.

        one_plus_eta(:) = 1.0 + qb(2,:) * one_over_pbprime(:)

        one_plus_eta_df(:) = 1.0 + qb_df(2,:) * one_over_pbprime_df(:)

        btp_mass_flux = qb(3:4,:)

        ! Compute pressure forcing H

        H(:) = (one_plus_eta(:)**2) * H_bcl(:)

        ! Compute  Q_u  and  Q_v  at the quadrature points in each cell.
        
        Quu(:) = ub * qb(3,:) + one_plus_eta(:) * Q_uu_dp(:)
        Quv(:) = ub * qb(4,:) + one_plus_eta(:) * Q_uv_dp(:)
        Qvv(:) = vb * qb(4,:) + one_plus_eta(:) * Q_vv_dp(:)

    end subroutine compute_btp_volume_terms

    subroutine compute_btp_face_terms(qb_face)

        ! This routine computes the flux terms for barotropic system

        use mod_grid, only: npoin_q, nface, npoin, face, intma_dg_quad
        use mod_basis, only: nq
        use mod_face, only: imapl_q, imapr_q, normal_vector_q
        use mod_input, only: nlayers, cd_mlswe, botfr
        use mod_initial, only: alpha_mlswe
        use mod_initial, only: coeff_pbpert_L, coeff_pbub_LR, coeff_pbpert_R, one_over_pbprime_edge, &
                                one_over_pbprime, one_over_pbprime_face, one_over_pbprime_df, &
                                coeff_mass_pbub_L, coeff_mass_pbub_R, coeff_mass_pbpert_LR
        use mod_constants, only: gravity

        use mod_variables, only: Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, H_bcl_edge
        use mod_variables, only: Qu_face, Qv_face, one_plus_eta_face, flux_edge, &
                                    one_plus_eta_df, H_face, one_plus_eta_edge_2

        implicit none
        real, dimension(4,2,nq,nface), intent(in) :: qb_face

        integer :: iface, iquad, Iq
        real, dimension(nq) :: ul, ur, vl, vr, upl, upr, vpl, vpr
        real :: pbub_left, pbvb_left, pbpert_left
        real :: pbub_right, pbvb_right, pbpert_right
        real :: pU_L, pU_R, nxl, nyl, nxr, nyr, rho
        real :: pbpert_edge, eta_edge, speed1
        integer :: ir, jr, kr, kl, ier, jquad, el, er, jl, il
        real, dimension(nq,nface) :: one_plus_eta_edge

        ! Compute eta

        do iface = 1, nface

            el = face(7,iface)
            er = face(8,iface)

            ! Use a Riemann problem to determine interpolated values of 
            ! pbpert = p'_b*eta at element faces.
            ! Here, eta iface the relative perturbation in bottom pressure, 
            ! relative to the global rest state that iface used to define the 
            ! barotropic-baroclinic splitting.

            ! Regard the interpolation process as Lax-Friedrichs interpolation,
            ! or a generalization thereof.
            ! Then compute the relative perturbations in bottom pressure 
            ! at each cell edge.  That is, 
            ! eta_edge(i) = pbpert_edge i / (p'_b at edge i)
            ! Also compute  (1 + eta_edge)**2 . 

            ! Compute pbpert_edge, eta_edge, and one_plus_eta_edge_2

            do iquad = 1, nq
                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)
        
                nxr = -nxl
                nyr = -nyl
        
                pbub_left = qb_face(3, 1, iquad, iface)
                pbvb_left = qb_face(4, 1, iquad, iface)
                pbpert_left = qb_face(2, 1, iquad, iface)
        
                pbub_right = qb_face(3, 2, iquad, iface)
                pbvb_right = qb_face(4, 2, iquad, iface)
                pbpert_right = qb_face(2, 2, iquad, iface)
        
                pU_L = nxl * pbub_left + nyl * pbvb_left
                pU_R = nxr * pbub_right + nyr * pbvb_right
        
                pbpert_edge = coeff_pbpert_L(iquad, iface) * pbpert_left &
                                            + coeff_pbpert_R(iquad, iface) * pbpert_right &
                                            + coeff_pbub_LR(iquad, iface) * (pU_L + pU_R)
        
                eta_edge = pbpert_edge * one_over_pbprime_edge(iquad, iface)
                one_plus_eta_edge(iquad, iface) = 1.0 + eta_edge
                one_plus_eta_edge_2(iquad, iface) = (1.0 + eta_edge)

                ! Compute mass fluxes at each element face.

                flux_edge(1,iquad,iface) = coeff_mass_pbub_L(iquad,iface) * pbub_left &
                                            + coeff_mass_pbub_R(iquad,iface) * pbub_right &
                                            + coeff_mass_pbpert_LR(iquad,iface) * (nxl * pbpert_left + nxr * pbpert_right)

                flux_edge(2,iquad,iface) = coeff_mass_pbub_L(iquad,iface) * pbvb_left &
                                            + coeff_mass_pbub_R(iquad,iface) * pbvb_right &
                                            + coeff_mass_pbpert_LR(iquad,iface) * (nyl * pbpert_left + nyr * pbpert_right)
            end do

            ! Compute (1 + eta) and (1 + eta)**2 at each element faces( 1 = left & 2=right side).

            one_plus_eta_face(1,:,iface) = 1.0 + qb_face(2,1,:,iface) * one_over_pbprime_face(1,:,iface)
            one_plus_eta_face(2,:,iface) = 1.0 + qb_face(2,2,:,iface) * one_over_pbprime_face(2,:,iface)

            ! Compute pressure forcing H_face at each element face.

            H_face(:,iface) = (one_plus_eta_edge_2(:,iface)**2) * H_bcl_edge(:,iface)

            ! Compute momentum flux terms at each element face 

            ul = qb_face(3,1,:,iface)/qb_face(1,1,:,iface); ur = qb_face(3,2,:,iface)/qb_face(1,2,:,iface)
            vl = qb_face(4,1,:,iface)/qb_face(1,1,:,iface); vr = qb_face(4,2,:,iface)/qb_face(1,2,:,iface)

            upl = qb_face(3, 1, :, iface); upr = qb_face(3, 2, :, iface)
            vpl = qb_face(4, 1, :, iface); vpr = qb_face(4, 2, :, iface)

            Qu_face(1,:,iface) = 0.5*(ul*upl + ur*upr) + one_plus_eta_edge(:,iface) * Q_uu_dp_edge(:, iface)
            Qu_face(2,:,iface) = 0.5*(vl*upl + vr*upr) + one_plus_eta_edge(:,iface) * Q_uv_dp_edge(:, iface)

            Qv_face(1,:,iface) = 0.5*(ul*vpl + ur*vpr) + one_plus_eta_edge(:,iface) * Q_uv_dp_edge(:, iface)
            Qv_face(2,:,iface) = 0.5*(vl*vpl + vr*vpr) + one_plus_eta_edge(:,iface) * Q_vv_dp_edge(:, iface)

        end do

    end subroutine compute_btp_face_terms

    subroutine btp_interpolate_avg(uvdp_ave)

        use mod_variables, only: ope_ave, ope_face_ave, uvb_face_ave, uvb_ave, uvb_ave_df, ope_ave_df, btp_mass_flux_ave, &
                                Qu_ave, Quv_ave, Qv_ave, tau_bot_ave, tau_bot_ave_temp, Qu_ave_temp, Quv_ave_temp, Qv_ave_temp

        use mod_face, only: imapl_q, imapr_q, normal_vector_q

        ! Interpolate btp mass pb from nodal to quad points

        use mod_basis, only: ngl, nq, psiq, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad, face, nface

        use mod_initial, only: psih, indexq

        implicit none

        real, dimension(2,npoin), intent(in) :: uvdp_ave

        integer :: e, iquad, jquad, I, Iq, m, n, l, kquad, ip
        integer :: iface, il, jl, ir, jr, el, er, ilocl, ilocr,kl,kr
        real :: un, nx, ny, hn
        
        uvb_ave  = 0.0
        ope_ave = 0.0
        btp_mass_flux_ave = 0.0
        !Qu_ave = 0.0
        !Quv_ave = 0.0
        !Qv_ave = 0.0
        !tau_bot_ave = 0.0

        do Iq = 1,npoin_q
            do ip = 1,npts

                I = indexq(ip,Iq)
                hn = psih(ip,Iq)

                uvb_ave(1,Iq) = uvb_ave(1,Iq) + hn*uvb_ave_df(1,I)
                uvb_ave(2,Iq) = uvb_ave(2,Iq) + hn*uvb_ave_df(2,I)

                ope_ave(Iq) = ope_ave(Iq) + hn*ope_ave_df(I)

                !Qu_ave(Iq) = Qu_ave(Iq) + hn*Qu_ave_temp(I)
                !Quv_ave(Iq) = Quv_ave(Iq) + hn*Quv_ave_temp(I)
                !Qv_ave(Iq) = Qv_ave(Iq) + hn*Qv_ave_temp(I)

                !tau_bot_ave(1,Iq) = tau_bot_ave(1,Iq) + hn*tau_bot_ave_temp(1,I)
                !tau_bot_ave(2,Iq) = tau_bot_ave(2,Iq) + hn*tau_bot_ave_temp(2,I)

                btp_mass_flux_ave(1,Iq) = btp_mass_flux_ave(1,Iq) + hn*uvdp_ave(1,I)
                btp_mass_flux_ave(2,Iq) = btp_mass_flux_ave(2,Iq) + hn*uvdp_ave(2,I)

            end do

        end do

        do iface = 1, nface                  !i specifies the grid cell

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                I=intma_dg_quad(il,jl,kl,el)

                uvb_face_ave(:,1,iquad,iface) = uvb_ave(:,I)
                !ope_face_ave(1,iquad,iface) = ope_ave(I)

                if(er > 0) then
                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    I=intma_dg_quad(ir,jr,kr,er)

                    uvb_face_ave(:,2,iquad,iface) = uvb_ave(:,I)
                    !ope_face_ave(2,iquad,iface) = ope_ave(I)
                else    
                    uvb_face_ave(:,2,iquad,iface) = uvb_face_ave(:,1,iquad,iface)
                    !ope_face_ave(2,iquad,iface) = ope_face_ave(1,iquad,iface)

                    if(er == -4) then

                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)

                        un = nx*uvb_ave(1,I) + ny*uvb_ave(2,I)

                        uvb_face_ave(1,2,iquad,iface) = uvb_ave(1,I) - 2.0*un*nx
                        uvb_face_ave(2,2,iquad,iface) = uvb_ave(2,I) - 2.0*un*ny

                    elseif(er == -2) then 
                        uvb_face_ave(:,2,iquad,iface) = -uvb_face_ave(:,1,iquad,iface)

                    end if
                end if

            end do
      
        end do
        
    end subroutine btp_interpolate_avg


    subroutine btp_interpolate_avg_v1()

        use mod_variables, only: uvb_face_ave, uvb_ave, uvb_ave_df

        use mod_face, only: imapl_q, imapr_q, normal_vector_q

        ! Interpolate btp mass pb from nodal to quad points

        use mod_basis, only: ngl, nq, psiq, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad, face, nface

        use mod_initial, only: psih, indexq

        implicit none

        integer :: e, iquad, jquad, I, Iq, m, n, l, kquad, ip
        integer :: iface, il, jl, ir, jr, el, er, ilocl, ilocr,kl,kr
        real :: un, nx, ny, hn
        
        uvb_ave  = 0.0
        uvb_face_ave = 0.0

        do Iq = 1,npoin_q
            do ip = 1,npts

                I = indexq(ip,Iq)
                hn = psih(ip,Iq)

                uvb_ave(1,Iq) = uvb_ave(1,Iq) + hn*uvb_ave_df(1,I)
                uvb_ave(2,Iq) = uvb_ave(2,Iq) + hn*uvb_ave_df(2,I)

            end do

        end do

        do iface = 1, nface                  !i specifies the grid cell

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                I=intma_dg_quad(il,jl,kl,el)

                uvb_face_ave(:,1,iquad,iface) = uvb_ave(:,I)

                if(er > 0) then
                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    I=intma_dg_quad(ir,jr,kr,er)

                    uvb_face_ave(:,2,iquad,iface) = uvb_ave(:,I)
                else    
                    uvb_face_ave(:,2,iquad,iface) = uvb_face_ave(:,1,iquad,iface)

                    if(er == -4) then

                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)

                        un = nx*uvb_ave(1,I) + ny*uvb_ave(2,I)

                        uvb_face_ave(1,2,iquad,iface) = uvb_ave(1,I) - 2.0*un*nx
                        uvb_face_ave(2,2,iquad,iface) = uvb_ave(2,I) - 2.0*un*ny

                    elseif(er == -2) then 
                        uvb_face_ave(:,2,iquad,iface) = -uvb_face_ave(:,1,iquad,iface)

                    end if
                end if

            end do
      
        end do

        call create_communicator_quad(uvb_face_ave,2)
        
    end subroutine btp_interpolate_avg_v1

    subroutine compute_btp_mom_terms(qb,qb_face,qprime,qb_df)

        use mod_grid, only: npoin_q, nface, npoin, face, intma_dg_quad
        use mod_basis, only: nq
        use mod_face, only: imapl_q, imapr_q, normal_vector_q
        use mod_input, only: nlayers, cd_mlswe, botfr
        use mod_initial, only: alpha_mlswe
        use mod_initial, only: coeff_pbpert_L, coeff_pbub_LR, coeff_pbpert_R, one_over_pbprime_edge, &
                                one_over_pbprime, one_over_pbprime_face, one_over_pbprime_df, &
                                coeff_mass_pbub_L, coeff_mass_pbub_R, coeff_mass_pbpert_LR
        use mod_constants, only: gravity

        use mod_variables, only: Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, H_bcl_edge
        use mod_variables, only: Qu_face, Qv_face, one_plus_eta_face, Quu, Qvv, Quv, H, one_plus_eta, &
                                    one_plus_eta_df, tau_bot, H_face, one_plus_eta_edge_2

        implicit none

        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(4,2,nq,nface), intent(in) :: qb_face
        real, dimension(4,npoin_q), intent(in) :: qb

        integer :: iface, iquad, Iq
        real, dimension(nq) :: ul, ur, vl, vr, upl, upr, vpl, vpr
        real :: pbub_left, pbvb_left, pbpert_left
        real :: pbub_right, pbvb_right, pbpert_right
        real :: pU_L, pU_R, nxl, nyl, nxr, nyr, rho
        real :: pbpert_edge, eta_edge, speed1
        integer :: ir, jr, kr, kl, ier, jquad, el, er, jl, il
        real, dimension(npoin_q) :: ubot, vbot, speed
        real, dimension(nq,nface) :: one_plus_eta_edge
        real, dimension(npoin_q) :: ub, vb

        ub = qb(3,:)/qb(1,:)
        vb = qb(4,:)/qb(1,:)


        ! Stress terms

        if(botfr == 1) then 

            ubot = qprime(1,:,nlayers)*(qprime(2,:,nlayers) + ub)
            vbot = qprime(1,:,nlayers)*(qprime(3,:,nlayers) + vb)
            speed1 = cd_mlswe/gravity
        
            tau_bot(1,:) = speed1 * ubot
            tau_bot(2,:) = speed1 * vbot
        elseif(botfr == 2) then 

            rho = 1.0/alpha_mlswe(nlayers)

            ubot = qprime(2,:,nlayers) + ub
            vbot = qprime(3,:,nlayers) + vb
            speed = rho*cd_mlswe*sqrt(ubot**2 + vbot**2)
            
            tau_bot(1,:) = speed * ubot
            tau_bot(2,:) = speed * vbot
        end if 
        ! Compute (1 + eta) and (1 + eta)**2 quadrature points of each grid cell.

        one_plus_eta(:) = 1.0 + qb(2,:) * one_over_pbprime(:)

        one_plus_eta_df(:) = 1.0 + qb_df(2,:) * one_over_pbprime_df(:)

        ! Compute H

        H(:) = (one_plus_eta(:)**2) * H_bcl(:)

        ! Compute  Q_u  and  Q_v  at the quadrature points in each cell.
        
        Quu(:) = ub * qb(3,:) + one_plus_eta(:) * Q_uu_dp(:)
        Quv(:) = ub * qb(4,:) + one_plus_eta(:) * Q_uv_dp(:)
        Qvv(:) = vb * qb(4,:) + one_plus_eta(:) * Q_vv_dp(:)

        ! Compute eta

        do iface = 1, nface

            el = face(7,iface)
            er = face(8,iface)

            ! Use a Riemann problem to determine interpolated values of 
            ! pbpert = p'_b*eta at element faces.
            ! Here, eta iface the relative perturbation in bottom pressure, 
            ! relative to the global rest state that iface used to define the 
            ! barotropic-baroclinic splitting.

            ! Regard the interpolation process as Lax-Friedrichs interpolation,
            ! or a generalization thereof.
            ! Then compute the relative perturbations in bottom pressure 
            ! at each cell edge.  That is, 
            ! eta_edge(i) = pbpert_edge i / (p'_b at edge i)
            ! Also compute  (1 + eta_edge)**2 . 

            ! Compute pbpert_edge, eta_edge, and one_plus_eta_edge_2

            do iquad = 1, nq
                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)
        
                nxr = -nxl
                nyr = -nyl
        
                pbub_left = qb_face(3, 1, iquad, iface)
                pbvb_left = qb_face(4, 1, iquad, iface)
                pbpert_left = qb_face(2, 1, iquad, iface)
        
                pbub_right = qb_face(3, 2, iquad, iface)
                pbvb_right = qb_face(4, 2, iquad, iface)
                pbpert_right = qb_face(2, 2, iquad, iface)
        
                pU_L = nxl * pbub_left + nyl * pbvb_left
                pU_R = nxr * pbub_right + nyr * pbvb_right
        
                pbpert_edge = coeff_pbpert_L(iquad, iface) * pbpert_left &
                                            + coeff_pbpert_R(iquad, iface) * pbpert_right &
                                            + coeff_pbub_LR(iquad, iface) * (pU_L + pU_R)
        
                eta_edge = pbpert_edge * one_over_pbprime_edge(iquad, iface)
                one_plus_eta_edge(iquad, iface) = 1.0 + eta_edge
                one_plus_eta_edge_2(iquad, iface) = (1.0 + eta_edge)

            end do

            ! Compute (1 + eta) and (1 + eta)**2 at each element faces( 1 = left & 2=right side).

            one_plus_eta_face(1,:,iface) = 1.0 + qb_face(2,1,:,iface) * one_over_pbprime_face(1,:,iface)
            one_plus_eta_face(2,:,iface) = 1.0 + qb_face(2,2,:,iface) * one_over_pbprime_face(2,:,iface)

            ! Compute H_face at each element face.

            H_face(:,iface) = (one_plus_eta_edge_2(:,iface)**2) * H_bcl_edge(:,iface)

            if(er ==-4) then 

                do iquad = 1, nq

                    il = imapl_q(1,iquad,1,iface)
                    jl = imapl_q(2,iquad,1,iface)
                    kl = imapl_q(3,iquad,1,iface)

                    Iq = intma_dg_quad(il,jl,kl,el)

                    H_face(iquad,iface) = H(Iq)
                    
                end do
            end if

            ul = qb_face(3,1,:,iface)/qb_face(1,1,:,iface); ur = qb_face(3,2,:,iface)/qb_face(1,2,:,iface)
            vl = qb_face(4,1,:,iface)/qb_face(1,1,:,iface); vr = qb_face(4,2,:,iface)/qb_face(1,2,:,iface)

            upl = qb_face(3, 1, :, iface); upr = qb_face(3, 2, :, iface)
            vpl = qb_face(4, 1, :, iface); vpr = qb_face(4, 2, :, iface)

            Qu_face(1,:,iface) = 0.5*(ul*upl + ur*upr) + one_plus_eta_edge(:,iface) * Q_uu_dp_edge(:, iface)
            Qu_face(2,:,iface) = 0.5*(vl*upl + vr*upr) + one_plus_eta_edge(:,iface) * Q_uv_dp_edge(:, iface)

            Qv_face(1,:,iface) = 0.5*(ul*vpl + ur*vpr) + one_plus_eta_edge(:,iface) * Q_uv_dp_edge(:, iface)
            Qv_face(2,:,iface) = 0.5*(vl*vpl + vr*vpr) + one_plus_eta_edge(:,iface) * Q_vv_dp_edge(:, iface)

        end do

    end subroutine compute_btp_mom_terms

    subroutine btp_mass_advection_terms(qb_face)

        use mod_basis, only: nq
        use mod_grid, only: nface, mod_grid_get_face_nq, face
        use mod_face, only: normal_vector_q
        use mod_initial, only: coeff_mass_pbub_L, coeff_mass_pbub_R, coeff_mass_pbpert_LR
        use mod_face, only: imapl_q, imapr_q
        use mod_variables, only: flux_edge

        implicit none

        real, intent(in) :: qb_face(4,2,nq,nface)

        integer :: iface, iquad, nq_i, nq_j, plane_ij, ilocl, il, jl ,kl, Iq, el, er
        real :: nxl, nyl, nxr, nyr, pbub_left, pbvb_left, pbpert_left, pbub_right, pbvb_right, pbpert_right
        real :: unpl

        ! =============== Begin computation of fluxes at boundaries ================

        do iface = 1,nface

            ilocl = face(5,iface)
            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                nxr = -nxl
                nyr = -nyl

                pbub_left = qb_face(3,1,iquad,iface)
                pbvb_left = qb_face(4,1,iquad,iface)
                pbpert_left = qb_face(2,1,iquad,iface)

                pbub_right = qb_face(3,2,iquad,iface)
                pbvb_right = qb_face(4,2,iquad,iface)          
                pbpert_right = qb_face(2,2,iquad,iface)

                flux_edge(1,iquad,iface) = coeff_mass_pbub_L(iquad,iface) * pbub_left &
                                            + coeff_mass_pbub_R(iquad,iface) * pbub_right &
                                            + coeff_mass_pbpert_LR(iquad,iface) * (nxl * pbpert_left + nxr * pbpert_right)

                flux_edge(2,iquad,iface) = coeff_mass_pbub_L(iquad,iface) * pbvb_left &
                                            + coeff_mass_pbub_R(iquad,iface) * pbvb_right &
                                            + coeff_mass_pbpert_LR(iquad,iface) * (nyl * pbpert_left + nyr * pbpert_right)

            end do 
        enddo

    end subroutine btp_mass_advection_terms


    subroutine btp_evaluate_pb(qb, qb_df)

        ! Interpolate btp mass pb from nodal to quad points

        use mod_initial, only: pbprime
        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad

        use mod_initial, only: psih, indexq

        implicit none
        

        real, intent(inout) :: qb(4,npoin_q) 
        real, intent(in) :: qb_df(4,npoin) 

        integer :: e, iquad, jquad, I, Iq, m, n, l, kquad, ip
        real :: hn
        
        qb(1:2,:) = 0.0

        do Iq = 1,npoin_q
            do ip = 1,npts

                I = indexq(ip,Iq)
                hn = psih(ip,Iq)

                qb(2,Iq) = qb(2,Iq) + hn*qb_df(2,I)

            end do

            qb(1,Iq) = qb(2,Iq) + pbprime(Iq)

        end do
        
    end subroutine btp_evaluate_pb


    subroutine btp_evaluate_mom(qb,qb_df)

        ! Interpolate btp momentum ub*pb, vb*pb from nodal to quad points


        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad

        use mod_initial, only: psih, indexq
    
        implicit none
        
        real, intent(out) :: qb(4,npoin_q)
        real, intent(in) :: qb_df(4,npoin)

        integer :: e, iquad, jquad, I, Iq, m, n, l, kquad, ip
        real :: hi

        ! Compute values of certain barotropic functions at quadrature points
        ! in each cell and endpoints of each cell, at barotropic
        ! time level m+1.
        
        qb(3:4,:) = 0.0

        ! Evaluate pbub and pbvb.

        do Iq = 1,npoin_q
            do ip = 1,npts
                
                I = indexq(ip,Iq)
                hi = psih(ip,Iq)
                
                qb(3,Iq) = qb(3,Iq) + hi*qb_df(3,I)
                qb(4,Iq) = qb(4,Iq) + hi*qb_df(4,I)
                
            end do

            if (qb(1,Iq) <= 0.0) then
                write(*,*) 'Nonpositive depth in barotropic.m'
                stop
            end if
        end do
    
    end subroutine btp_evaluate_mom

    subroutine btp_evaluate_mom_dp(qb,qb_df)

        ! Interpolate btp pb, ubpb, vbpb from nodal to quad points

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad

        use mod_initial, only: psih, indexq, pbprime
    
        implicit none
        
        real, intent(out) :: qb(4,npoin_q)
        real, intent(in) :: qb_df(4,npoin)

        integer :: e, iquad, jquad, I, Iq, m, n, l, kquad, ip
        real :: hi

        ! Compute values of certain barotropic functions at quadrature points
        ! in each cell and endpoints of each cell, at barotropic
        ! time level m+1.
        
        qb = 0.0

        ! Evaluate pbub and pbvb.

        do Iq = 1,npoin_q
            do ip = 1,npts
                
                I = indexq(ip,Iq)
                hi = psih(ip,Iq)

                qb(2,Iq) = qb(2,Iq) + hi*qb_df(2,I)
                qb(3,Iq) = qb(3,Iq) + hi*qb_df(3,I)
                qb(4,Iq) = qb(4,Iq) + hi*qb_df(4,I)
                
            end do

            qb(1,Iq) = qb(2,Iq) + pbprime(Iq)
        end do

        if (any(qb(1,:) <= 0.0)) then
            write(*,*) 'Nonpositive depth in barotropic.m'
            stop
        end if
    
    end subroutine btp_evaluate_mom_dp
    
    subroutine btp_evaluate_pb_face(qb_face, qb)

        ! Extract btp pb (qb(1,:)) face values to qb_face(1,:,:,:)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl, nq
        use mod_grid, only:  npoin_q, intma_dg_quad, nface, face,mod_grid_get_face_nq
        use mod_initial, only: pbprime_face
        use mod_face, only: imapl_q, imapr_q
  
        implicit none
        
        ! Input variables
        real, intent(inout) :: qb_face(4,2,nq,nface)
        real, intent(in) :: qb(4,npoin_q)
        
        ! Local variables
        integer :: il, jl, ir, jr, iquad, jquad, el, er, ilocl, ilocr, I,kl,kr, iface, nq_i, nq_j, plane_ij

        qb_face(1:2,:,:,:) = 0.0
        
        ! Evaluate  pbpert = p'_b*eta  and  pb.
        do iface = 1, nface        
          
            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)
          
            do iquad = 1,nq   
                
                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                I=intma_dg_quad(il,jl,kl,el)

                qb_face(1,1,iquad,iface) = qb(1,I)
                qb_face(2,1,iquad,iface) = qb(2,I)

                if(er > 0) then
                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    I=intma_dg_quad(ir,jr,kr,er)

                    qb_face(1,2,iquad,iface) = qb(1,I)
                    qb_face(2,2,iquad,iface) = qb(2,I)
                else    
                    qb_face(1,2,iquad,iface) = qb_face(1,1,iquad,iface)
                    qb_face(2,2,iquad,iface) = qb_face(2,1,iquad,iface)
                end if
            
            enddo
          
        enddo
      
    end subroutine btp_evaluate_pb_face

    subroutine btp_evaluate_mom_face(qb_face, qb)

        ! Extract btp ubpb, vbpb (qb(3:4,:)) face values to qb_face(3:4,:,:,:)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl,nq
        use mod_grid, only:  npoin_q, intma_dg_quad, nface, face,mod_grid_get_face_nq
        use mod_initial, only: pbprime_face
        use mod_face, only: imapl_q, imapr_q, normal_vector_q

        implicit none

        real, intent(inout) :: qb_face(4,2,nq,nface)
        real, intent(in) :: qb(4,npoin_q)

        integer :: iface, ilr, iquad,jquad, m, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr,plane_ij,nq_i,nq_j
        real :: un, nx, ny

      
        ! Compute values of certain barotropic functions at quadrature points
        ! in each cell and endpoints of each cell, at barotropic
        ! time level m+1.
        ! The values at endpoints should be interpreted as one-sided limits.
      
        ! In the case of  pbub,  pbvb,  and  pbpert,  use degrees of freedom 
        ! and basis functions.  In the case of  ub  and  vb,  use division.
      
        qb_face(3:4,:,:,:) = 0.0
      
        ! Evaluate  pbub  and  pbvb.
      
        do iface = 1, nface                  !i specifies the grid cell

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                ! Left
                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                I=intma_dg_quad(il,jl,kl,el)

                qb_face(3:4,1,iquad,iface) = qb(3:4,I)

                if(er > 0) then
                    ! Right
                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    I=intma_dg_quad(ir,jr,kr,er)

                    qb_face(3:4,2,iquad,iface) = qb(3:4,I)
                else

                    qb_face(3:4,2,iquad,iface) = qb_face(3:4,1,iquad,iface)

                    if(er == -4) then

                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)

                        un = nx*qb(3,I) + ny*qb(4,I)

                        qb_face(3,2,iquad,iface) = qb(3,I) - 2.0*un*nx
                        qb_face(4,2,iquad,iface) = qb(4,I) - 2.0*un*ny

                    elseif(er == -2) then 
                        qb_face(3:4,2,iquad,iface) = -qb_face(3:4,1,iquad,iface)

                    end if

                end if

            end do
      
        end do
      
    end subroutine btp_evaluate_mom_face

    subroutine btp_evaluate_mom_dp_face(qb_face, qb)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl,nq
        use mod_grid, only:  npoin_q, intma_dg_quad, nface, face,mod_grid_get_face_nq
        use mod_initial, only: pbprime_face
        use mod_face, only: imapl_q, imapr_q, normal_vector_q

        implicit none

        real, intent(inout) :: qb_face(4,2,nq,nface)
        real, intent(in) :: qb(4,npoin_q)

        integer :: iface, ilr, iquad,jquad, m, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr,plane_ij,nq_i,nq_j
        real :: un, nx, ny

      
        ! Compute values of certain barotropic functions at quadrature points
        ! in each cell and endpoints of each cell, at barotropic
        ! time level m+1.
        ! The values at endpoints should be interpreted as one-sided limits.
      
        ! In the case of  pbub,  pbvb,  and  pbpert,  use degrees of freedom 
        ! and basis functions.  In the case of  ub  and  vb,  use division.
      
        qb_face = 0.0
      
        ! Evaluate  pbub  and  pbvb.
      
        do iface = 1, nface                  !i specifies the grid cell

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                ! Left
                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                I=intma_dg_quad(il,jl,kl,el)

                qb_face(1:4,1,iquad,iface) = qb(1:4,I)

                if(er > 0) then
                    ! Right
                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    I=intma_dg_quad(ir,jr,kr,er)

                    qb_face(1:4,2,iquad,iface) = qb(1:4,I)
                else

                    qb_face(1:4,2,iquad,iface) = qb_face(1:4,1,iquad,iface)

                    if(er == -4) then

                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)

                        un = nx*qb(3,I) + ny*qb(4,I)

                        qb_face(3,2,iquad,iface) = qb(3,I) - 2.0*un*nx
                        qb_face(4,2,iquad,iface) = qb(4,I) - 2.0*un*ny

                    elseif(er == -2) then 
                        qb_face(3:4,2,iquad,iface) = -qb_face(3:4,1,iquad,iface)
                        
                    end if

                end if

            end do
      
        end do
      
    end subroutine btp_evaluate_mom_dp_face  

    subroutine btp_extract_df(qb_df_face, qb_df)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl, nq, psiq
        use mod_grid, only:  npoin, intma, nface, face
        use mod_face, only: imapl, imapr, normal_vector

        implicit none

        real, intent(inout) :: qb_df_face(4,2,ngl,nface)
        real, intent(in) :: qb_df(4,npoin)

        integer :: iface, ilr, iquad,jquad, n, m, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr
        real :: un, nx, ny, hi

      
        ! Compute values of certain barotropic functions at quadrature points
        ! in each cell and endpoints of each cell, at barotropic
        ! time level m+1.
        ! The values at endpoints should be interpreted as one-sided limits.
      
        ! In the case of  pbub,  pbvb,  and  pbpert,  use degrees of freedom 
        ! and basis functions.  In the case of  ub  and  vb,  use division.
      
        qb_df_face = 0.0
      
        ! Evaluate  pbub  and  pbvb.
      
        do iface = 1, nface                  !i specifies the grid cell

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do n = 1,ngl

                ! Left
                il=imapl(1,n,1,iface)
                jl=imapl(2,n,1,iface)
                kl=imapl(3,n,1,iface)
                I=intma(il,jl,kl,el)

                qb_df_face(1:4,1,n,iface) = qb_df(1:4,I)

                if(er > 0) then

                    ! Right
                    ir=imapr(1,n,1,iface)
                    jr=imapr(2,n,1,iface)
                    kr=imapr(3,n,1,iface)
                    I=intma(ir,jr,kr,er)

                    qb_df_face(1:4,2,n,iface) = qb_df(1:4,I)
                else

                    qb_df_face(1:4,2,n,iface) = qb_df_face(1:4,1,n,iface)

                    if(er == -4) then

                        nx = normal_vector(1,n,1,iface)
                        ny = normal_vector(2,n,1,iface)

                        un = nx*qb_df_face(3,1,n,iface) + ny*qb_df_face(4,1,n,iface)

                        qb_df_face(3,2,n,iface) = qb_df_face(3,1,n,iface) - 2.0*un*nx
                        qb_df_face(4,2,n,iface) = qb_df_face(4,1,n,iface) - 2.0*un*ny

                    elseif(er == -2) then 
                        qb_df_face(3:4,2,n,iface) = -qb_df_face(3:4,1,n,iface)
                        
                    end if

                end if
        
            end do
        end do 
      
    end subroutine btp_extract_df 

    subroutine btp_extract_face(qb_face, qb_df)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl, nq, psiq
        use mod_grid, only:  npoin, intma, nface, face
        use mod_face, only: imapl, imapr, normal_vector_q

        implicit none

        real, intent(out) :: qb_face(4,2,nq,nface)
        real, intent(in) :: qb_df(4,npoin)

        integer :: iface, ilr, iquad,jquad, n, m, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr
        real :: un, nx, ny, hi

      
        ! Compute values of certain barotropic functions at quadrature points
        ! in each cell and endpoints of each cell, at barotropic
        ! time level m+1.
        ! The values at endpoints should be interpreted as one-sided limits.
      
        ! In the case of  pbub,  pbvb,  and  pbpert,  use degrees of freedom 
        ! and basis functions.  In the case of  ub  and  vb,  use division.
      
        qb_face = 0.0
      
        ! Evaluate  pbub  and  pbvb.
      
        do iface = 1, nface                  !i specifies the grid cell

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1,nq

                ! Left
                do n = 1,ngl

                    hi = psiq(n,iquad)

                    il=imapl(1,n,1,iface)
                    jl=imapl(2,n,1,iface)
                    kl=imapl(3,n,1,iface)
                    I=intma(il,jl,kl,el)

                    qb_face(1:4,1,iquad,iface) = qb_face(1:4,1,iquad,iface) + hi*qb_df(1:4,I)
                end do 

                ! Right
                if(er > 0) then

                    do n = 1,ngl

                        hi = psiq(n,iquad)

                        ir=imapr(1,n,1,iface)
                        jr=imapr(2,n,1,iface)
                        kr=imapr(3,n,1,iface)
                        I=intma(ir,jr,kr,er)

                        qb_face(1:4,2,iquad,iface) = qb_face(1:4,2,iquad,iface) + hi*qb_df(1:4,I)
                    end do 
                else

                    qb_face(1:4,2,iquad,iface) = qb_face(1:4,1,iquad,iface)

                    if(er == -4) then

                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)

                        un = nx*qb_face(3,1,iquad,iface) + ny*qb_face(4,1,iquad,iface)

                        qb_face(3,2,iquad,iface) = qb_face(3,1,iquad,iface) - 2.0*un*nx
                        qb_face(4,2,iquad,iface) = qb_face(4,1,iquad,iface) - 2.0*un*ny

                    elseif(er == -2) then 
                        qb_face(3:4,2,iquad,iface) = -qb_face(3:4,1,iquad,iface)
                        
                    end if

                end if
            end do
        end do 
      
    end subroutine btp_extract_face 

    subroutine btp_extract_face_v1(qb_face, qb_df_face)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl, nq, psiq
        use mod_grid, only:  npoin, intma, nface, face
        use mod_face, only: imapl, imapr, normal_vector_q

        implicit none

        real, intent(out) :: qb_face(4,2,nq,nface)
        real, intent(in) :: qb_df_face(4,2,ngl,nface)

        integer :: iface, ilr, iquad,jquad, n, m, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr
        real :: un, nx, ny, hi

      
        ! Compute values of certain barotropic functions at quadrature points
        ! in each cell and endpoints of each cell, at barotropic
        ! time level m+1.
        ! The values at endpoints should be interpreted as one-sided limits.
      
        ! In the case of  pbub,  pbvb,  and  pbpert,  use degrees of freedom 
        ! and basis functions.  In the case of  ub  and  vb,  use division.
      
        qb_face = 0.0
      
        ! Evaluate  pbub  and  pbvb.
      
        do iface = 1, nface                  !i specifies the grid cell

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1,nq

                ! Left
                do n = 1,ngl

                    hi = psiq(n,iquad)

                    il=imapl(1,n,1,iface)
                    jl=imapl(2,n,1,iface)
                    kl=imapl(3,n,1,iface)
                    I=intma(il,jl,kl,el)

                    qb_face(1:4,1,iquad,iface) = qb_face(1:4,1,iquad,iface) + hi*qb_df_face(1:4,1,n,iface)
                end do 

                ! Right
                if(er > 0) then

                    do n = 1,ngl

                        hi = psiq(n,iquad)

                        ir=imapr(1,n,1,iface)
                        jr=imapr(2,n,1,iface)
                        kr=imapr(3,n,1,iface)
                        I=intma(ir,jr,kr,er)

                        qb_face(1:4,2,iquad,iface) = qb_face(1:4,2,iquad,iface) + hi*qb_df_face(1:4,2,n,iface)
                    end do 

                end if
            end do
        end do 
      
    end subroutine btp_extract_face_v1

    subroutine btp_mom_boundary_df(qb)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl,nq
        use mod_grid, only:  npoin, intma, nface, face,mod_grid_get_face_nq
        use mod_initial, only: pbprime_face
        use mod_face, only: imapl, imapr, normal_vector
        use mod_input, only: mlswe_bc_strong

        implicit none

        real, intent(inout) :: qb(2,npoin)

        integer :: iface, ilr, iquad,jquad, m, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr,plane_ij,nq_i,nq_j
        integer :: bcflag
        real :: nx, ny, unl, upnl

      
	    do iface = 1, nface                  !i specifies the grid cell

		!Store Left Side Variables
		el = face(7,iface)
		er = face(8,iface)

		if(er == -4) then

		    do iquad = 1, ngl

		        il=imapl(1,iquad,1,iface)
		        jl=imapl(2,iquad,1,iface)
		        kl=imapl(3,iquad,1,iface)
		        I=intma(il,jl,kl,el)

		        nx = normal_vector(1,iquad,1,iface)
		        ny = normal_vector(2,iquad,1,iface)

		        unl = qb(1,I)*nx + qb(2,I)*ny

		        qb(1,I) = qb(1,I) - unl*nx
		        qb(2,I) = qb(2,I) - unl*ny
		        
		    end do

		end if

	    end do
      
    end subroutine btp_mom_boundary_df
 

    subroutine btp_bcl_coeffs(qprime,qprime_face)
        
        ! Compute baroclinic coefficients in the advective barotropic momentum fluxes and in the barotropic pressure forcing. 

        use mod_grid, only: npoin_q, nface, npoin
        use mod_input, only: nlayers, dpprime_visc_min
        use mod_basis, only: nqx, nqy, nqz, nq
        use mod_initial, only: alpha_mlswe

        use mod_Tensorproduct, only: compute_gradient_quad, interpolate_layer_from_quad_to_node_1d
        use mod_variables, only: Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, H_bcl_edge, &
                                    Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge
        
        implicit none

        real, intent(in)  :: qprime(3,npoin_q,nlayers), qprime_face(3,2,nq,nface,nlayers)
        
        integer :: k, iface, Iq, iquad
        real, dimension(nq) :: left_uudp, right_uudp, left_uvdp
        real, dimension(nq) :: right_uvdp, right_vvdp, left_vvdp

        real, dimension(nq,nlayers+1) :: pprime_l, pprime_l2, pprime_r, pprime_r2
        real, dimension(npoin_q,nlayers+1) :: pprime, pprime_2
        real, dimension(nq) :: left_dp, right_dp
        
        ! Compute Q_uu_dp, Q_uv_dp, Q_vv_dp at quadrature points
        Q_uu_dp = 0.0
        Q_uv_dp = 0.0
        Q_vv_dp = 0.0
        H_bcl = 0.0
        H_bcl_edge = 0.0

        ! Compute Q_up_up_quad and Q_up_vp_quad at quadrature points.

        pprime_2(:,1) = 0.0
        pprime(:,1) = 0.0

        do k = 1, nlayers
            Q_uu_dp(:) = Q_uu_dp(:) + qprime(2,:,k)*(qprime(2,:,k) * qprime(1,:,k))
            Q_uv_dp(:) = Q_uv_dp(:) + qprime(3,:,k)*(qprime(2,:,k) * qprime(1,:,k))
            Q_vv_dp(:) = Q_vv_dp(:) + qprime(3,:,k)*(qprime(3,:,k) * qprime(1,:,k))

            ! Pressure term in barotropic momentum equation

            pprime(:,k+1) = pprime(:,k) + qprime(1,:,k)
            pprime_2(:,k+1) = pprime(:,k+1)**2
            H_bcl(:) = H_bcl(:) + 0.5*alpha_mlswe(k)*(pprime_2(:,k+1) - pprime_2(:,k))

        end do
        
        do iface = 1, nface
                
            left_uudp = 0.0
            right_uudp = 0.0
            left_uvdp = 0.0
            right_uvdp = 0.0
            right_vvdp = 0.0
            left_vvdp = 0.0

            left_dp = 0.0
            right_dp = 0.0

            pprime_l(:,1) = 0.0
            pprime_l2(:,1) = 0.0
            pprime_r(:,1) = 0.0
            pprime_r2(:,1) = 0.0

            do k = 1, nlayers

                left_uudp = left_uudp + qprime_face(2,1,:,iface,k)*qprime_face(2,1,:,iface,k)*qprime_face(1,1,:,iface,k)
                left_uvdp = left_uvdp + qprime_face(3,1,:,iface,k)*qprime_face(2,1,:,iface,k)*qprime_face(1,1,:,iface,k)
                left_vvdp = left_vvdp + qprime_face(3,1,:,iface,k)*qprime_face(3,1,:,iface,k)*qprime_face(1,1,:,iface,k)
                right_uudp = right_uudp + qprime_face(2,2,:,iface,k)*qprime_face(2,2,:,iface,k)*qprime_face(1,2,:,iface,k)
                right_uvdp = right_uvdp + qprime_face(3,2,:,iface,k)*qprime_face(2,2,:,iface,k)*qprime_face(1,2,:,iface,k)
                right_vvdp = right_vvdp + qprime_face(3,2,:,iface,k)*qprime_face(3,2,:,iface,k)*qprime_face(1,2,:,iface,k)

                pprime_l(:,k+1) = pprime_l(:,k) + qprime_face(1,1,:,iface,k)
                pprime_l2(:,k+1) = pprime_l(:,k+1)**2
                left_dp = left_dp + 0.5*alpha_mlswe(k) *(pprime_l2(:,k+1) - pprime_l2(:,k))

                pprime_r(:,k+1) = pprime_r(:,k) + qprime_face(1,2,:,iface,k)
                pprime_r2(:,k+1) = pprime_r(:,k+1)**2
                right_dp = right_dp + 0.5*alpha_mlswe(k) *(pprime_r2(:,k+1) - pprime_r2(:,k))

            end do

            Q_uu_dp_edge(:,iface) = 0.5*(left_uudp + right_uudp)
            !Q_uu_dp_edge(2,:,iface) = right_uudp
            Q_uv_dp_edge(:,iface) = 0.5*(left_uvdp + right_uvdp)
            !Q_uv_dp_edge(2,:,iface) = right_uvdp
            Q_vv_dp_edge(:,iface) = 0.5*(left_vvdp + right_vvdp)
            !Q_vv_dp_edge(2,:,iface) = right_vvdp

            H_bcl_edge(:,iface) = 0.5*(left_dp + right_dp)
            !H_bcl_edge(2,:,iface) = right_dp

        end do

    end subroutine btp_bcl_coeffs


    subroutine btp_bcl_grad_coeffs(qprime_df)

        ! Compute baroclinic coefficients that appear in the weak forms of the horizontal viscosity terms in the 
        ! layer momentum equations and in the vertically-summed barotropic momentum equations. 
        
        use mod_grid, only: npoin_q, nface, npoin, face, intma
        use mod_input, only: nlayers, dpprime_visc_min
        use mod_basis, only: nqx, nqy, nqz, nq, ngl
        use mod_initial, only: alpha_mlswe
        use mod_variables, only: dpprime_visc,pbprime_visc,btp_dpp_graduv, btp_dpp_uvp, &
                dpp_uvp, dpp_graduv, graduv_dpp_face, btp_graduv_dpp_face
        use mod_face, only: imapl_q, imapr_q, normal_vector_q, imapl, imapr, normal_vector

        implicit none

        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        real, dimension(4, npoin) :: graduv

        integer :: iface, il, jl, kl, ir, jr, kr, iel, ier, iquad, Iq, c_jump,k, ivar
        real, dimension(2,npoin) :: Uk
        real :: un(nlayers), nx, ny

        btp_dpp_graduv = 0.0
        pbprime_visc = 0.0

        do k = 1,nlayers

            Uk = qprime_df(2:3,:,k)
            
            call compute_gradient_uv_v1(graduv, Uk)
            
            dpp_graduv(1,:,k) = dpprime_visc(:,k)*graduv(1,:)
            dpp_graduv(2,:,k) = dpprime_visc(:,k)*graduv(2,:)
            dpp_graduv(3,:,k) = dpprime_visc(:,k)*graduv(3,:)
            dpp_graduv(4,:,k) = dpprime_visc(:,k)*graduv(4,:)

            dpp_uvp(1,:,k) =  dpprime_visc(:,k)*qprime_df(2,:,k)
            dpp_uvp(2,:,k) =  dpprime_visc(:,k)*qprime_df(3,:,k)

            btp_dpp_graduv(:,:) = btp_dpp_graduv(:,:) + dpp_graduv(:,:,k)
            pbprime_visc(:) =  pbprime_visc(:) + dpprime_visc(:,k)

        end do

        do iface=1,nface
            !-------------------------------------
            !Store Left and Right Side Variables
            !-------------------------------------
            iel=face(7,iface)
            ier=face(8,iface)
    
            !----------------------------Left Element
            do iquad = 1,ngl
                !Get Pointers
                il=imapl(1,iquad,1,iface)
                jl=imapl(2,iquad,1,iface)
                kl=imapl(3,iquad,1,iface)
                Iq = intma(il,jl,kl,iel)

                !Variables
                graduv_dpp_face(1,1,iquad,iface,:) = dpp_graduv(1,Iq,:)
                graduv_dpp_face(2,1,iquad,iface,:) = dpp_graduv(2,Iq,:)
                graduv_dpp_face(3,1,iquad,iface,:) = dpp_graduv(3,Iq,:)
                graduv_dpp_face(4,1,iquad,iface,:) = dpp_graduv(4,Iq,:)

                graduv_dpp_face(5,1,iquad,iface,:) = dpprime_visc(Iq,:)

                if (ier > 0 ) then

                    !Get Pointers
                    ir=imapr(1,iquad,1,iface)
                    jr=imapr(2,iquad,1,iface)
                    kr=imapr(3,iquad,1,iface)
                    Iq=intma(ir,jr,kr,ier)

                    !Variables
                    graduv_dpp_face(1,2,iquad,iface,:) = dpp_graduv(1,Iq,:)
                    graduv_dpp_face(2,2,iquad,iface,:) = dpp_graduv(2,Iq,:)
                    graduv_dpp_face(3,2,iquad,iface,:) = dpp_graduv(3,Iq,:)
                    graduv_dpp_face(4,2,iquad,iface,:) = dpp_graduv(4,Iq,:)

                    graduv_dpp_face(5,2,iquad,iface,:) = dpprime_visc(Iq,:)

                else
                    !default values
                    graduv_dpp_face(1,2,iquad,iface,:) = graduv_dpp_face(1,1,iquad,iface,:)
                    graduv_dpp_face(2,2,iquad,iface,:) = graduv_dpp_face(2,1,iquad,iface,:)

                    graduv_dpp_face(3,2,iquad,iface,:) = graduv_dpp_face(3,1,iquad,iface,:)
                    graduv_dpp_face(4,2,iquad,iface,:) = graduv_dpp_face(4,1,iquad,iface,:)

                    graduv_dpp_face(5,2,iquad,iface,:) = graduv_dpp_face(5,1,iquad,iface,:)

                    if(ier == -4) then 
                        nx = normal_vector(1,iquad,1,iface)
                        ny = normal_vector(2,iquad,1,iface)

                        un = dpp_graduv(1,Iq,:)*nx + dpp_graduv(2,Iq,:)*ny

                        graduv_dpp_face(1,2,iquad,iface,:) = dpp_graduv(1,Iq,:) - 2.0*un*nx
                        graduv_dpp_face(2,2,iquad,iface,:) = dpp_graduv(2,Iq,:) - 2.0*un*ny

                        un = dpp_graduv(3,Iq,:)*nx + dpp_graduv(4,Iq,:)*ny

                        graduv_dpp_face(3,2,iquad,iface,:) = dpp_graduv(3,Iq,:) - 2.0*un*nx
                        graduv_dpp_face(4,2,iquad,iface,:) = dpp_graduv(4,Iq,:) - 2.0*un*ny

                    end if 
                end if
            end do 

        end do

        call bcl_create_communicator(graduv_dpp_face,5,nlayers,ngl)

        btp_graduv_dpp_face = 0.0

        do iface=1,nface
            do iquad = 1,ngl
                
                do k = 1,nlayers

                    btp_graduv_dpp_face(:,1,iquad,iface) = btp_graduv_dpp_face(:,1,iquad,iface) + graduv_dpp_face(:,1,iquad,iface,k)

                    btp_graduv_dpp_face(:,2,iquad,iface) = btp_graduv_dpp_face(:,2,iquad,iface) + graduv_dpp_face(:,2,iquad,iface,k)

                end do 

            end do 
        end do 

    end subroutine btp_bcl_grad_coeffs

    subroutine evaluate_quprime(quprime, qvprime, qp, qp_face)

        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq, nface,face, intma, npoin
        use mod_basis, only: npts, nq, ngl, psiq
        use mod_metrics, only: massinv
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq

        implicit none

        real, dimension(2, npoin_q), intent(out) :: quprime
        real, dimension(2, npoin_q), intent(out) :: qvprime

        real, dimension(2, npoin_q), intent(in) :: qp
        real, dimension(2, 2, nq, nface), intent(in) :: qp_face

        real, dimension(2, npoin) :: qu_df, qv_df

        real :: nxl, nyl, dhdx, dhdy, wq, var_u, var_v, hi, um, vm
        integer :: Iq, ip, I, iface, el, er, iquad, ir, jr, kr, il, jl, kl, n

        quprime = 0.0
        qvprime = 0.0
        qu_df = 0.0
        qv_df = 0.0

        do Iq = 1,npoin_q
            
            wq = wjac(Iq)
            var_u = qp(1,Iq)
            var_v = qp(2,Iq)

            do ip = 1, npts

                I = indexq(ip,Iq)

                !Xi derivatives
                dhdx = dpsidx(ip,Iq)
                !Eta derivatives
                dhdy = dpsidy(ip,Iq)

                qu_df(1,I) = qu_df(1,I) - wq*var_u*dhdx
                qu_df(2,I) = qu_df(2,I) - wq*var_u*dhdy

                qv_df(1,I) = qv_df(1,I) - wq*var_v*dhdx
                qv_df(2,I) = qv_df(2,I) - wq*var_v*dhdy

            end do
        end do

        do iface = 1,nface 

            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                um = 0.5*(qp_face(1,1,iquad,iface) + qp_face(1,2,iquad,iface))
                vm = 0.5*(qp_face(2,1,iquad,iface) + qp_face(2,2,iquad,iface))

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    qu_df(1,I) = qu_df(1,I) + wq*hi*nxl*um
                    qu_df(2,I) = qu_df(2,I) + wq*hi*nyl*um

                    qv_df(1,I) = qv_df(1,I) + wq*hi*nxl*vm
                    qv_df(2,I) = qv_df(2,I) + wq*hi*nyl*vm

                    if(er > 0) then
                            
                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        qu_df(1,I) = qu_df(1,I) - wq*hi*nxl*um
                        qu_df(2,I) = qu_df(2,I) - wq*hi*nyl*um

                        qv_df(1,I) = qv_df(1,I) - wq*hi*nxl*vm
                        qv_df(2,I) = qv_df(2,I) - wq*hi*nyl*vm

                    end if

                end do
            end do 
        end do

        qu_df(1,:) = massinv(:)*qu_df(1,:)
        qu_df(2,:) = massinv(:)*qu_df(2,:)
        qv_df(1,:) = massinv(:)*qv_df(1,:)
        qv_df(2,:) = massinv(:)*qv_df(2,:)

        do Iq = 1,npoin_q
            do ip = 1,npts
                
                I = indexq(ip,Iq)
                hi = psih(ip,Iq)
                
                quprime(1,Iq) = quprime(1,Iq) + hi*qu_df(1,I)
                quprime(2,Iq) = quprime(2,Iq) + hi*qu_df(2,I)
                
                qvprime(1,Iq) = qvprime(1,Iq) + hi*qv_df(1,I)
                qvprime(2,Iq) = qvprime(2,Iq) + hi*qv_df(2,I)
                
            end do
        end do


    end subroutine evaluate_quprime

    subroutine evaluate_visc_terms(grad_qpvisc,grad_qpvisc_face, qpvisc, qpvisc_face)

        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq, nface,face, intma, npoin
        use mod_basis, only: npts, nq, ngl, psiq
        use mod_metrics, only: massinv
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac
        use mod_face, only: imapl_q, imapr_q, normal_vector_q, jac_faceq

        implicit none

        real, dimension(2,npoin_q), intent(in) :: qpvisc
        real, dimension(2,2,nq,nface), intent(out) :: qpvisc_face
        real, dimension(2,2,npoin_q), intent(out) :: grad_qpvisc
        real, dimension(4,2,nq,nface), intent(out) :: grad_qpvisc_face

        integer :: Iq, iface, iel, ier, iquad, ir, jr, kr, il, jl, kl

        qpvisc_face = 0.0

        do iface=1,nface
            !-------------------------------------
            !Store Left and Right Side Variables
            !-------------------------------------
            iel=face(7,iface)
            ier=face(8,iface)

            !----------------------------Left Element
            do iquad = 1,nq
                !Get Pointers
                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                Iq = intma_dg_quad(il,jl,kl,iel)

                !Variables
                qpvisc_face(:,1,iquad,iface) = qpvisc(:,Iq)

                if (ier > 0 ) then

                    !Get Pointers
                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    Iq=intma_dg_quad(ir,jr,kr,ier)

                    !Variables
                    qpvisc_face(:,2,iquad,iface) = qpvisc(:,Iq)

                else
                    !default values
                    qpvisc_face(:,2,iquad,iface) = qpvisc_face(:,1,iquad,iface)
                end if
            end do 
        end do

        call create_communicator_quad(qpvisc_face,2)
        call evaluate_quprime2(grad_qpvisc, qpvisc,qpvisc_face)

        do iface=1,nface
            !-------------------------------------
            !Store Left and Right Side Variables
            !-------------------------------------
            iel=face(7,iface)
            ier=face(8,iface)

            !----------------------------Left Element
            do iquad = 1,nq
                !Get Pointers
                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                Iq = intma_dg_quad(il,jl,kl,iel)

                !Variables

                grad_qpvisc_face(1:2,1,iquad,iface) = grad_qpvisc(1,:,Iq)
                grad_qpvisc_face(3:4,1,iquad,iface) = grad_qpvisc(2,:,Iq)

                if (ier > 0 ) then

                    !Get Pointers
                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    Iq=intma_dg_quad(ir,jr,kr,ier)

                    !Variables
                    grad_qpvisc_face(1:2,2,iquad,iface) = grad_qpvisc(1,:,Iq)
                    grad_qpvisc_face(3:4,2,iquad,iface) = grad_qpvisc(2,:,Iq)
          
                else
                    !default values
                    grad_qpvisc_face(1:2,2,iquad,iface) = grad_qpvisc_face(1:2,1,iquad,iface)
                    grad_qpvisc_face(3:4,2,iquad,iface) = grad_qpvisc_face(3:4,1,iquad,iface)
                end if
            end do 
        end do

        call create_communicator_quad(grad_qpvisc_face,4)

    end subroutine evaluate_visc_terms


    subroutine evaluate_quprime2(quvprime, qp, qp_face)

        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq, nface,face, intma, npoin
        use mod_basis, only: npts, nq, ngl, psiq
        use mod_metrics, only: massinv
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq

        implicit none

        real, dimension(2,2, npoin_q), intent(out) :: quvprime

        real, dimension(2, npoin_q), intent(in) :: qp
        real, dimension(2, 2, nq, nface), intent(in) :: qp_face

        real, dimension(2, npoin) :: qu_df, qv_df

        real :: nxl, nyl, dhdx, dhdy, wq, var_u, var_v, hi, um, vm
        integer :: Iq, ip, I, iface, el, er, iquad, ir, jr, kr, il, jl, kl, n

        quvprime = 0.0
        qu_df = 0.0
        qv_df = 0.0

        do Iq = 1,npoin_q
            
            wq = wjac(Iq)
            var_u = qp(1,Iq)
            var_v = qp(2,Iq)

            do ip = 1, npts

                I = indexq(ip,Iq)

                !Xi derivatives
                dhdx = dpsidx(ip,Iq)
                !Eta derivatives
                dhdy = dpsidy(ip,Iq)

                qu_df(1,I) = qu_df(1,I) - wq*var_u*dhdx
                qu_df(2,I) = qu_df(2,I) - wq*var_u*dhdy

                qv_df(1,I) = qv_df(1,I) - wq*var_v*dhdx
                qv_df(2,I) = qv_df(2,I) - wq*var_v*dhdy

            end do
        end do

        do iface = 1,nface 

            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                um = 0.5*(qp_face(1,1,iquad,iface) + qp_face(1,2,iquad,iface))
                vm = 0.5*(qp_face(2,1,iquad,iface) + qp_face(2,2,iquad,iface))

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    qu_df(1,I) = qu_df(1,I) + wq*hi*nxl*um
                    qu_df(2,I) = qu_df(2,I) + wq*hi*nyl*um

                    qv_df(1,I) = qv_df(1,I) + wq*hi*nxl*vm
                    qv_df(2,I) = qv_df(2,I) + wq*hi*nyl*vm

                    if(er > 0) then
                            
                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        qu_df(1,I) = qu_df(1,I) - wq*hi*nxl*um
                        qu_df(2,I) = qu_df(2,I) - wq*hi*nyl*um

                        qv_df(1,I) = qv_df(1,I) - wq*hi*nxl*vm
                        qv_df(2,I) = qv_df(2,I) - wq*hi*nyl*vm

                    end if

                end do
            end do 
        end do

        qu_df(1,:) = massinv(:)*qu_df(1,:)
        qu_df(2,:) = massinv(:)*qu_df(2,:)
        qv_df(1,:) = massinv(:)*qv_df(1,:)
        qv_df(2,:) = massinv(:)*qv_df(2,:)

        do Iq = 1,npoin_q
            do ip = 1,npts
                
                I = indexq(ip,Iq)
                hi = psih(ip,Iq)
                
                quvprime(1,1,Iq) = quvprime(1,1,Iq) + hi*qu_df(1,I)
                quvprime(1,2,Iq) = quvprime(1,2,Iq) + hi*qu_df(2,I)
                
                quvprime(2,1,Iq) = quvprime(2,1,Iq) + hi*qv_df(1,I)
                quvprime(2,2,Iq) = quvprime(2,2,Iq) + hi*qv_df(2,I)
                
            end do
        end do


    end subroutine evaluate_quprime2


    subroutine restart_mlswe_variales(qr,qbr,qprimer,q,qb,qprime, q_df,qb_df,qprime_df,q_face,qb_face,qprime_face)

        use mod_grid, only : npoin_q, npoin, intma_dg_quad, intma,nface,face
        use mod_basis, only: npts
        use mod_input, only: nlayers
        use mod_initial, only: psih, indexq, pbprime_df
        use mod_face, only: imapl_q, imapr_q


        real, dimension(3,npoin,nlayers), intent(in) :: qprimer
        real, dimension(3,npoin,nlayers), intent(in) :: qr
        real, dimension(3,npoin), intent(in) :: qbr

        real, dimension(5,npoin_q,nlayers), intent(out) :: q
        real, dimension(5,npoin,nlayers), intent(out) :: q_df
        real, dimension(6,npoin_q), intent(out) :: qb
        real, dimension(6,npoin), intent(out) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(out) :: qprime
        real, dimension(3,npoin,nlayers), intent(out) :: qprime_df
        real, dimension(5,2,nq,nface,nlayers), intent(out) :: q_face
        real, dimension(3,2,nq,nface,nlayers), intent(out) :: qprime_face
        real, dimension(6,2,nq,nface), intent(out) :: qb_face

        integer :: I, Iq, ip, iface, el, er, iquad, ir, jr, kr, il, jl, kl, n,k, ilocl, ilocr
        real :: hi

        q = 0.0
        qb = 0.0
        qprime = 0.0
        q_face = 0.0
        qprime_face = 0.0
        qb_face = 0.0

        ! Barotropic variables at the dofs (nodal points)
        qb_df(1,:) = qbr(1,:) + pbprime_df(:)
        qb_df(2:4,:) = qbr(1:3,:)
        qb_df(5,:) = qbr(2,:)*qbr(1,:)
        qb_df(6,:) = qbr(3,:)*qbr(1,:)

        do k = 1,nlayers

            ! Baroclinic variables at the dofs (nodal points)
            q_df(1:3,:,k) = qr(1:3,:,k)
            q_df(4,:,k) = q_df(2,:,k)*q_df(1,:,k)
            q_df(5,:,k) = q_df(3,:,k)*q_df(1,:,k)

            ! Prime variables at the dofs (nodal points)
            qprime_df(:,:,k) = qprimer(:,:,k)
            

            do Iq = 1,npoin_q 
                do ip = 1,npts 

                    I = indexq(ip,Iq)

                    q(:,Iq,k) = q(:,Iq,k) + psih(ip,Iq)*q_df(:,I,k)
                    qprime(:,Iq,k) = qprime(:,Iq,k) + psih(ip,Iq)*qprime_df(:,I,k)

                end do
            end do 
        end do 

        do Iq = 1,npoin_q 
            do ip = 1,npts 

                I = indexq(ip,Iq) 

                qb(:,Iq) = qb(:,Iq) + psih(ip,Iq)*qb_df(:,I)
            end do 
        end do 

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

                q_face(:,1,iquad,iface,:) = q(:,I,:)
                qprime_face(:,1,iquad,iface,:) = qprime(:,I,:)
                qb_face(:,1,iquad,iface) = qb(:,I)

                if(er > 0) then

                    ir = imapr_q(1,iquad,1,iface)
                    jr = imapr_q(2,iquad,1,iface)
                    kr = imapr_q(3,iquad,1,iface)

                    I = intma_dg_quad(ir,jr,kr,er)

                    q_face(:,2,iquad,iface,:) = q(:,I,:)
                    qprime_face(:,2,iquad,iface,:) = qprime(:,I,:)
                    qb_face(:,2,iquad,iface) = qb(:,I)

                else

                    q_face(:,2,iquad,iface,:) = q_face(:,1,iquad,iface,:)
                    qprime_face(:,2,iquad,iface,:) = qprime_face(:,1,iquad,iface,:)
                    qb_face(:,2,iquad,iface) = qb_face(:,1,iquad,iface)

                end if
            end do 
        end do 
    
    end subroutine restart_mlswe_variales

    subroutine restart_mlswe(q_df,qb_df,q,qb,qprime,qprime_df,q_face,qprime_face,qb_face, qp_df_out, q_df_read, qb_df_read)

        use mod_grid, only : npoin_q, npoin, intma_dg_quad, intma,nface,face
        use mod_basis, only: npts
        use mod_input, only: nlayers
        use mod_initial, only: psih, indexq, pbprime_df, alpha_mlswe, pbprime, zbot_df
        use mod_face, only: imapl_q, imapr_q
        use mod_constants, only: gravity
        use mod_layer_terms, only: evaluate_mom, evaluate_mom_face, evaluate_dp, evaluate_dp_face

        real, dimension(3,npoin), intent(in) :: qb_df_read
        real, dimension(3,npoin,nlayers), intent(in) :: q_df_read

        real, dimension(4,npoin), intent(out) :: qb_df
        real, dimension(3,npoin,nlayers), intent(out) :: q_df
        real, dimension(4,npoin_q), intent(out) :: qb
        real, dimension(3,npoin_q,nlayers), intent(out) :: q

        real, dimension(3,npoin_q,nlayers), intent(out) :: qprime
        real, dimension(3,npoin,nlayers), intent(out) :: qprime_df
        real, dimension(3,2,nq,nface,nlayers), intent(out) :: q_face
        real, dimension(3,2,nq,nface,nlayers), intent(out) :: qprime_face
        real, dimension(4,2,nq,nface), intent(out) :: qb_face
        real, dimension(5,npoin,nlayers) :: qp_df_out

        real, dimension(npoin,nlayers+1) :: mslwe_elevation

        integer :: I, Iq, ip, iface, el, er, iquad, ir, jr, kr, il, jl, kl, n,k, ilocl, ilocr
        real :: hi
        real, dimension(npoin) :: one_plus_eta_temp

        q = 0.0
        qb = 0.0
        qprime = 0.0
        q_face = 0.0
        qprime_face = 0.0
        qb_face = 0.0

        ! Barotropic variables at the dofs (nodal points)

        qb_df(1,:) = qb_df_read(1,:)
        qb_df(3:4,:) = qb_df_read(2:3,:)
        qb_df(2,:) = qb_df(1,:) - pbprime_df(:)

        ! Interpolate to the quadrature points

        call btp_evaluate_mom_dp(qb,qb_df)
        call btp_evaluate_mom_dp_face(qb_face, qb)

        ! Baroclinic variables at the dofs (nodal points)

        !q_df(:,:,:) = q_df_read(:,:,:)

        do k = 1,nlayers

            q_df(1,:,k) = (gravity/alpha_mlswe(k))*q_df_read(1,:,k)
            q_df(2,:,k) = q_df_read(2,:,k)*q_df(1,:,k)
            q_df(3,:,k) = q_df_read(3,:,k)*q_df(1,:,k)

        end do

        ! Interpolate to the quadrature points

        call evaluate_dp(q,qprime,q_df, pbprime)
        call evaluate_dp_face(q_face, qprime_face,q, qprime)

        call evaluate_mom(q,q_df)
        call evaluate_mom_face(q_face, q)

        ! Prime variables at the dofs (nodal points) and quadrature points
        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)

        do k = 1,nlayers

            qprime_df(1,:,k) = q_df(1,:,k) / one_plus_eta_temp(:)

            qprime(2,:,k) = q(2,:,k)/q(1,:,k) - qb(3,:)/qb(1,:)
            qprime(3,:,k) = q(3,:,k)/q(1,:,k) - qb(4,:)/qb(1,:)
            qprime_df(2,:,k) = q_df(2,:,k)/q_df(1,:,k) - qb_df(3,:)/qb_df(1,:)
            qprime_df(3,:,k) = q_df(3,:,k)/q_df(1,:,k) - qb_df(4,:)/qb_df(1,:)

            qprime_face(2,:,:,:,k) = q_face(2,:,:,:,k)/q_face(1,:,:,:,k) - qb_face(3,:,:,:)/qb_face(1,:,:,:)
            qprime_face(3,:,:,:,k) = q_face(3,:,:,:,k)/q_face(1,:,:,:,k) - qb_face(4,:,:,:)/qb_face(1,:,:,:)
        end do

        ! Prepare output variables

        do k = 1,nlayers
            qp_df_out(1,:,k) = (alpha_mlswe(k)/gravity)*q_df(1,:,k)
            qp_df_out(2,:,k) = q_df(2,:,k) / q_df(1,:,k)
            qp_df_out(3,:,k) = q_df(3,:,k) / q_df(1,:,k)
        end do

        mslwe_elevation(:,nlayers+1) = zbot_df

        do k = nlayers,1,-1
            mslwe_elevation(:,k) = mslwe_elevation(:,k+1) + qp_df_out(1,:,k)
        end do

        qp_df_out(4,:,1) = qb_df(3,:)
        qp_df_out(4,:,2) = qb_df(3,:)

        qp_df_out(5,:,1) = mslwe_elevation(:,1)
        qp_df_out(5,:,2) = mslwe_elevation(:,2)
    
    end subroutine restart_mlswe

    subroutine compute_gradient_uv(grad_uv,uv)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad
        use mod_basis, only: nq, psiq, psiqx, psiqy, psiqz, dpsiq, dpsiqx, dpsiqy, dpsiqz
        use mod_metrics, only: ksiq_x,ksiq_y,ksiq_z, etaq_x,etaq_y,etaq_z, zetaq_x,zetaq_y,zetaq_z

        use mod_initial, only: dpsidx,dpsidy, indexq

        implicit none

        real, dimension(2,npoin), intent(in) :: uv
        real, dimension(2,2,npoin_q), intent(out) :: grad_uv
        integer :: Iq, I, ip
        real :: dhdx, dhdy
        
        grad_uv = 0.0
        
        !Construct Volume Integral Contributions

        do Iq = 1, npoin_q

            do ip = 1,npts

                I = indexq(ip,Iq)
                dhdx = dpsidx(ip,Iq)
                dhdy = dpsidy(ip,Iq)

                grad_uv(1,1,Iq) = grad_uv(1,1,Iq) + dhdx*uv(1,I)
                grad_uv(1,2,Iq) = grad_uv(1,2,Iq) + dhdy*uv(1,I)

                grad_uv(2,1,Iq) = grad_uv(2,1,Iq) + dhdx*uv(2,I)
                grad_uv(2,2,Iq) = grad_uv(2,2,Iq) + dhdy*uv(2,I)

            end do
        end do

    end subroutine compute_gradient_uv

    subroutine compute_gradient_uv_v1(grad_uv,uv)

        use mod_basis, only:  npts
        use mod_grid, only:  npoin

        use mod_initial, only: dpsidx_df,dpsidy_df, index_df

        implicit none

        real, dimension(2,npoin), intent(in) :: uv
        real, dimension(4,npoin), intent(out) :: grad_uv
        integer :: Iq, I, ip
        real :: dhdx, dhdy
        
        grad_uv = 0.0
        
        !Construct Volume Integral Contributions

        do Iq = 1, npoin

            do ip = 1,npts

                I = index_df(ip,Iq)
                dhdx = dpsidx_df(ip,Iq)
                dhdy = dpsidy_df(ip,Iq)

                grad_uv(1,Iq) = grad_uv(1,Iq) + dhdx*uv(1,I)
                grad_uv(2,Iq) = grad_uv(2,Iq) + dhdy*uv(1,I)

                grad_uv(3,Iq) = grad_uv(3,Iq) + dhdx*uv(2,I)
                grad_uv(4,Iq) = grad_uv(4,Iq) + dhdy*uv(2,I)

            end do
        end do

    end subroutine compute_gradient_uv_v1

end module mod_barotropic_terms
            