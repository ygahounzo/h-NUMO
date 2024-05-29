! ===========================================================================================================================
! This module contains the routines for the barotropic-baroclinic splitting
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: March 27, 2023
! ==========================================================================================================================

module mod_barotropic_terms

    
    use mod_grid, only: npoin_q, nface, intma_dg_quad
    use mod_basis, only: nqx, nqy, nqz, nq
    use mod_input, only: nlayers, mass_exact
        
    implicit none

    public :: btp_evaluate_mom, &
                btp_evaluate_mom_face, btp_evaluate_pb, btp_evaluate_pb_face, btp_mass_advection_terms, &
                btp_bcl_coeffs, btp_mom_boundary_df, &
                limiter_Positivity_Preserving_mom, limiter_Positivity_Preserving_mass, evaluate_quprime, &
                compute_btp_terms, btp_evaluate_mom_dp_face, btp_evaluate_mom_dp, btp_mom_boundary, &
                evaluate_quprime2, evaluate_visc_terms, restart_mlswe_variales, restart_mlswe, massinv_rhs, evaluate_quprime1, &
                compute_gradient_uv, compute_btp_mom_terms, compute_gradient_uv_v2

    contains

    subroutine compute_btp_terms(qb,qb_face, qprime, qb_df)

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

        btp_mass_flux = qb(3:4,:)

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

            ! Compute H_face at each element face.

            H_face(:,iface) = (one_plus_eta_edge_2(:,iface)**2) * 0.5*(H_bcl_edge(1,:,iface) + H_bcl_edge(2,:,iface))

            if(er < 0) then 

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

            Qu_face(1,:,iface) = 0.5*(ul*upl + ur*upr) + 0.5*one_plus_eta_edge(:,iface) * (Q_uu_dp_edge(1,:, iface) + Q_uu_dp_edge(2,:, iface))
            Qu_face(2,:,iface) = 0.5*(vl*upl + vr*upr) + 0.5*one_plus_eta_edge(:,iface) * (Q_uv_dp_edge(1,:, iface) + Q_uv_dp_edge(2,:, iface))

            Qv_face(1,:,iface) = 0.5*(ul*vpl + ur*vpr) + 0.5*one_plus_eta_edge(:,iface) * (Q_uv_dp_edge(1,:, iface) + Q_uv_dp_edge(2,:, iface))
            Qv_face(2,:,iface) = 0.5*(vl*vpl + vr*vpr) + 0.5*one_plus_eta_edge(:,iface) * (Q_vv_dp_edge(1,:, iface) + Q_vv_dp_edge(2,:, iface))

        end do

    end subroutine compute_btp_terms

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

            H_face(:,iface) = (one_plus_eta_edge_2(:,iface)**2) * 0.5*(H_bcl_edge(1,:,iface) + H_bcl_edge(2,:,iface))

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

            Qu_face(1,:,iface) = 0.5*(ul*upl + ur*upr) + 0.5*one_plus_eta_edge(:,iface) * (Q_uu_dp_edge(1,:, iface) + Q_uu_dp_edge(2,:, iface))
            Qu_face(2,:,iface) = 0.5*(vl*upl + vr*upr) + 0.5*one_plus_eta_edge(:,iface) * (Q_uv_dp_edge(1,:, iface) + Q_uv_dp_edge(2,:, iface))

            Qv_face(1,:,iface) = 0.5*(ul*vpl + ur*vpr) + 0.5*one_plus_eta_edge(:,iface) * (Q_uv_dp_edge(1,:, iface) + Q_uv_dp_edge(2,:, iface))
            Qv_face(2,:,iface) = 0.5*(vl*vpl + vr*vpr) + 0.5*one_plus_eta_edge(:,iface) * (Q_vv_dp_edge(1,:, iface) + Q_vv_dp_edge(2,:, iface))

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

                I = indexq(Iq,ip)
                hn = psih(Iq,ip)

                qb(2,Iq) = qb(2,Iq) + hn*qb_df(2,I)

            end do

            qb(1,Iq) = qb(2,Iq) + pbprime(Iq)

        end do
        
    end subroutine btp_evaluate_pb


    subroutine btp_evaluate_mom(qb,qb_df)


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
                
                I = indexq(Iq,ip)
                hi = psih(Iq,ip)
                
                qb(3,Iq) = qb(3,Iq) + hi*qb_df(3,I)
                qb(4,Iq) = qb(4,Iq) + hi*qb_df(4,I)
                
            end do

            if (qb(1,Iq) <= 0.0) then
                write(*,*) 'Nonpositive depth in barotropic.m'
                stop
            end if
        end do

         call btp_mom_boundary(qb)
    
    end subroutine btp_evaluate_mom

    subroutine btp_evaluate_mom_dp(qb,qb_df)


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
                
                I = indexq(Iq,ip)
                hi = psih(Iq,ip)

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

    subroutine btp_mom_boundary(qb)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl,nq
        use mod_grid, only:  npoin_q, intma_dg_quad, nface, face,mod_grid_get_face_nq
        use mod_initial, only: pbprime_face
        use mod_face, only: imapl_q, imapr_q, normal_vector_q

        implicit none

        real, intent(inout) :: qb(4,npoin_q)

        integer :: iface, ilr, iquad,jquad, m, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr,plane_ij,nq_i,nq_j
        integer :: bcflag
        real :: nx, ny, unl, upnl

      
        do iface = 1, nface                  !i specifies the grid cell

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            if(er == -4) then

                do iquad = 1, nq

                    il=imapl_q(1,iquad,1,iface)
                    jl=imapl_q(2,iquad,1,iface)
                    kl=imapl_q(3,iquad,1,iface)
                    I=intma_dg_quad(il,jl,kl,el)

                    nx = normal_vector_q(1,iquad,1,iface)
                    ny = normal_vector_q(2,iquad,1,iface)

                    unl = qb(3,I)*nx + qb(4,I)*ny

                    qb(3,I) = qb(3,I) - unl*nx
                    qb(4,I) = qb(4,I) - unl*ny
                    
                end do
            elseif(er == -2) then

                do iquad = 1, nq

                    il=imapl_q(1,iquad,1,iface)
                    jl=imapl_q(2,iquad,1,iface)
                    kl=imapl_q(3,iquad,1,iface)
                    I=intma_dg_quad(il,jl,kl,el)

                    qb(3,I) = 0.0
                    qb(4,I) = 0.0
                    
                end do

            end if

      
        end do
      
    end subroutine btp_mom_boundary

    subroutine limiter_Positivity_Preserving_mass(qb_df)

        use mod_basis, only : npts, nglx, ngly
        use mod_grid, only : npoin, nelem, intma
        
        implicit none

        !global arrays
        real, intent(inout) :: qb_df(npoin)

        !local arrays
        real, dimension(npts) :: qlocal, qlim
        real :: theta, qmean, qmin
        real :: eps, rho
        integer :: e, i, j, k, l, m, inodes(npts), ii, Ip

        !Define Constants
        eps=1.0e-13

        !loop through elements
        do e=1,nelem

            !Store local DOFs and compute Mean
            ii=0
            do j=1,ngly
                do i=1,nglx
                    Ip = intma(i,j,1,e)
                    ii=ii+1
                    inodes(ii) = Ip
                    !Store Primitive Variables
                    qlocal(ii) = qb_df(Ip)
                end do !i
            end do !j

            qmean = sum(qlocal(:)) / npts !area of canonical element in 3D (4 in 2D)

            !Compute Weight
            qmin = minval(qlocal(:))
            theta = min(1.0, qmean/(abs(qmean - qmin)  + eps ))

            !Apply Weighted combination of Mean Value and High-Order Solution
            do ii=1,npts
                qlim(ii) = qmean + abs(theta)*(qlocal(ii) - qmean)
                if (qmin == 0.0 .and. qmean == 0.0) qlim(ii) = 0.0
            end do !ii

            !Check for Positivity and Store Values
            do ii=1,npts
                Ip = inodes(ii)
                ! qlim(ii) = max(0.0, qlim(ii))
                qb_df(Ip) = qlim(ii) 
            end do !ii

        end do !e

    end subroutine limiter_Positivity_Preserving_mass

    subroutine limiter_Positivity_Preserving_mom(qb_df)

        use mod_basis, only : npts, nglx, ngly
        use mod_grid, only : npoin, nelem, intma

        implicit none

        !global arrays
        real, intent(inout) :: qb_df(3,npoin)

        !local arrays
        real, dimension(3,npts) :: qlocal, qlim
        real, dimension(3) :: qmean, qmin
        real :: eps, rho, theta
        integer :: e, i, j, k, l, m, inodes(npts), ii, Ip

        !Define Constants
        eps=1.0e-13

        !loop through elements
        do e = 1,nelem

            !Store local DOFs and compute Mean
            ii=0
            do j=1,ngly
                do i=1,nglx
                    Ip = intma(i,j,1,e)
                    ii=ii+1
                    inodes(ii) = Ip
                    !Store Primitive Variables
                    do m = 1,3
                        qlocal(m,ii) = qb_df(m,Ip)
                    end do
                end do !i
            end do !j

            qmean(:) = sum(qlocal(:,:),dim=2) / npts !area of canonical element in 3D (4 in 2D)

            !Compute Weight
            do m = 1,3
                qmin(m) = minval(qlocal(m,:))
            end do

            theta = min(1.0, qmean(1)/(abs(qmean(1) - qmin(1)) + eps))

            !Apply Weighted combination of Mean Value and High-Order Solution
            do ii=1,npts
                do m = 2,3
                    qlim(m,ii) = qmean(m) + abs(theta)*(qlocal(m,ii) - qmean(m))
                    if (qmin(m) == 0.0 .and. qmean(m) == 0.0) qlim(m,ii) = 0.0d0
                end do
            end do !ii

            !Check for Positivity and Store Values
            do ii=1,npts
                Ip = inodes(ii)
                do m = 2,3
                    qb_df(m,Ip) = qlim(m,ii) 
                end do
            end do !ii

        end do !e

    end subroutine limiter_Positivity_Preserving_mom       

    subroutine btp_bcl_coeffs(qprime,qprime_face)
        
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

            Q_uu_dp_edge(1,:,iface) = left_uudp
            Q_uu_dp_edge(2,:,iface) = right_uudp
            Q_uv_dp_edge(1,:,iface) = left_uvdp 
            Q_uv_dp_edge(2,:,iface) = right_uvdp
            Q_vv_dp_edge(1,:,iface) = left_vvdp 
            Q_vv_dp_edge(2,:,iface) = right_vvdp

            H_bcl_edge(1,:,iface) = left_dp
            H_bcl_edge(2,:,iface) = right_dp

        end do

    end subroutine btp_bcl_coeffs

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

                I = indexq(Iq,ip)

                !Xi derivatives
                dhdx = dpsidx(Iq,ip)
                !Eta derivatives
                dhdy = dpsidy(Iq,ip)

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

        ! if(mass_exact) then 
            ! call massinv_rhs(qu_df,2)
            ! call massinv_rhs(qv_df,2)
        ! else
            qu_df(1,:) = massinv(:)*qu_df(1,:)
            qu_df(2,:) = massinv(:)*qu_df(2,:)
            qv_df(1,:) = massinv(:)*qv_df(1,:)
            qv_df(2,:) = massinv(:)*qv_df(2,:)
        ! end if

        do Iq = 1,npoin_q
            do ip = 1,npts
                
                I = indexq(Iq,ip)
                hi = psih(Iq,ip)
                
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

                I = indexq(Iq,ip)

                !Xi derivatives
                dhdx = dpsidx(Iq,ip)
                !Eta derivatives
                dhdy = dpsidy(Iq,ip)

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
                
                I = indexq(Iq,ip)
                hi = psih(Iq,ip)
                
                quvprime(1,1,Iq) = quvprime(1,1,Iq) + hi*qu_df(1,I)
                quvprime(1,2,Iq) = quvprime(1,2,Iq) + hi*qu_df(2,I)
                
                quvprime(2,1,Iq) = quvprime(2,1,Iq) + hi*qv_df(1,I)
                quvprime(2,2,Iq) = quvprime(2,2,Iq) + hi*qv_df(2,I)
                
            end do
        end do


    end subroutine evaluate_quprime2

    subroutine evaluate_quprime1(qprime1, qp, qp_face)

        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq, nface,face, intma, npoin
        use mod_basis, only: npts, nq, ngl, psiq
        use mod_metrics, only: massinv
        use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac
        use mod_face, only: imapl, imapr, normal_vector_q, jac_faceq

        implicit none

        real, dimension(2, npoin_q), intent(out) :: qprime1

        real, dimension(npoin_q), intent(in) :: qp
        real, dimension(2, nq, nface), intent(in) :: qp_face

        real, dimension(2, npoin) :: qu_df

        real :: nxl, nyl, dhdx, dhdy, wq, var_u, var_v, hi, um, vm
        integer :: Iq, ip, I, iface, el, er, iquad, ir, jr, kr, il, jl, kl, n

        qprime1 = 0.0
        qu_df = 0.0

        do Iq = 1,npoin_q
            
            wq = wjac(Iq)
            var_u = qp(Iq)

            do ip = 1, npts

                I = indexq(Iq,ip)

                !Xi derivatives
                dhdx = dpsidx(Iq,ip)
                !Eta derivatives
                dhdy = dpsidy(Iq,ip)

                qu_df(1,I) = qu_df(1,I) - wq*var_u*dhdx
                qu_df(2,I) = qu_df(2,I) - wq*var_u*dhdy

            end do
        end do

        do iface = 1,nface 

            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                wq = jac_faceq(iquad,1,iface)

                nxl = normal_vector_q(1,iquad,1,iface)
                nyl = normal_vector_q(2,iquad,1,iface)

                um = 0.5*(qp_face(1,iquad,iface) + qp_face(2,iquad,iface))

                do n = 1, ngl

                    hi = psiq(n,iquad)
                    
                    il = imapl(1,n,1,iface)
                    jl = imapl(2,n,1,iface)
                    kl = imapl(3,n,1,iface)

                    I = intma(il,jl,kl,el)

                    qu_df(1,I) = qu_df(1,I) + wq*hi*nxl*um
                    qu_df(2,I) = qu_df(2,I) + wq*hi*nyl*um

                    if(er > 0) then
                            
                        ir = imapr(1,n,1,iface)
                        jr = imapr(2,n,1,iface)
                        kr = imapr(3,n,1,iface)

                        I = intma(ir,jr,kr,er)

                        qu_df(1,I) = qu_df(1,I) - wq*hi*nxl*um
                        qu_df(2,I) = qu_df(2,I) - wq*hi*nyl*um

                    end if

                end do
            end do 
        end do

        call massinv_rhs(qu_df,2)


        do Iq = 1,npoin_q
            do ip = 1,npts
                
                I = indexq(Iq,ip)
                hi = psih(Iq,ip)
                
                qprime1(1,Iq) = qprime1(1,Iq) + hi*qu_df(1,I)
                qprime1(2,Iq) = qprime1(2,Iq) + hi*qu_df(2,I)
                
            end do
        end do


    end subroutine evaluate_quprime1


    ! subroutine btp_wet_dry(qb,qb_df,neg_pb_pos,neg_pb_pos_q)

    !     use mod_grid, only: npoin, npoin_q
    !     use mod_input, only: depth_dry
    !     use mod_initial, only: pbprime, pbprime_df, alpha_mlswe
    !     use mod_constants, only: gravity
    
    !     implicit none
        
    !     real, intent(inout) :: qb(6,npoin_q)
    !     real, intent(inout) :: qb_df(6,npoin)
    !     integer, dimension(npoin),intent(out) :: neg_pb_pos
    !     integer, dimension(npoin_q),intent(out) :: neg_pb_pos_q

    !     integer :: I, Iq
    !     real :: h_limit, h, pb_limit
    !     real :: theta, qmin 
    !     real, dimension(3,npts) :: qlocal
    !     real, dimension(3) :: qmean

    !     !Define Constants
    !     eps=1.0e-3

    !     !loop through elements
    !     do e=1,nelem

    !         !Store local DOFs and compute Mean
    !         ii=0
    !         do j=1,ngly
    !             do i=1,nglx
    !                 Ip = intma(i,j,1,e)
    !                 ii=ii+1
    !                 inodes(ii) = Ip
    !                 !Store Primitive Variables
    !                 qlocal(1,ii) = qb_df(1,Ip)
    !                 qlocal(2,ii) = qb_df(5,Ip)
    !                 qlocal(3,ii) = qb_df(6,Ip)
    !             end do !i
    !         end do !j

    !         !Compute Weight
    !         qmin = minval(qlocal(1,:))

    !         if(qmin < depth_dry) then 

    !             qmean = sum(qlocal(:)) / npts !area of canonical element in 3D (4 in 2D)
    !             if(qmean < depth_dry) then
    !                 do ii=1,npts
    !                     Ip = inodes(ii)
    !                     qb_df(1,Ip) = (gravity/alpha_mlswe(nlayers)) * depth_dry
    !                     qb_df(5:6,Ip) = 0.0
    !                 end do !ii
    !             else 

    !                 theta = min(1.0, qmean/(qmean - qmin))

    !                 !Apply Weighted combination of Mean Value and High-Order Solution
    !                 do ii=1,npts
    !                     Ip = inodes(ii)
    !                     qb_df(1,Ip) = qmean(1) + theta*(qlocal(1,ii) - qmean(1))
    !                     qb_df(5,Ip) = qmean(2) + theta*(qlocal(5,ii) - qmean(2))
    !                     qb_df(6,Ip) = qmean(3) + theta*(qlocal(6,ii) - qmean(3))
    !                 end do !ii
    !             end if 
    !         end if 

    !     end do !e

    !     qb_df(2,:) = qb_df(1,:) - pbprime_df(:)

    !     neg_pb_pos = 0
    !     neg_pb_pos_q = 0

    !     h_limit = 50.0
    !     pb_limit = (gravity/sum(alpha_mlswe(:)))*h_limit

    !     ! Compute the wetting and drying treatment for the barotropic

    !     do I = 1,npoin
    !         if(qb_df(1,I) < pb_limit) then 
    !             neg_pb_pos(I) = 1
    !             qb_df(1,I) = pb_limit
    !             qb_df(2,I) = qb_df(1,I) - pbprime_df(I)
    !             qb_df(3:6,I) = 0.0
    !         end if
    !     end do

    !     do Iq = 1,npoin_q
    !         if (qb(1,Iq) < pb_limit) then
    !             neg_pb_pos_q(I) = 1
    !             qb(1,Iq) = pb_limit
    !             qb(2,Iq) = qb(1,Iq) - pbprime(Iq)
    !             qb(3:6,Iq) = 0.0
    !         end if
    !     end do
    
    ! end subroutine btp_wet_dry

    subroutine btp_wet_dry_mom(qb,qb_df,neg_pb_pos,neg_pb_pos_q)

        use mod_grid, only: npoin, npoin_q

        use mod_initial, only: pbprime, pbprime_df, alpha_mlswe
        use mod_constants, only: gravity
    
        implicit none
        
        real, intent(inout) :: qb(6,npoin_q)
        real, intent(inout) :: qb_df(6,npoin)
        integer, dimension(npoin),intent(in) :: neg_pb_pos
        integer, dimension(npoin_q),intent(in) :: neg_pb_pos_q

        integer :: I, Iq
        real :: h_limit

        h_limit = 0.1

        ! Compute the wetting and drying treatment for the barotropic

        do I = 1,npoin
            if(neg_pb_pos(I) == 1) then 
                qb_df(3:6,I) = 0.0
            end if
        end do

        do Iq = 1,npoin_q
            
            if (neg_pb_pos_q(Iq) == 1) then
                qb(3:6,Iq) = 0.0
            end if
        end do
    
    end subroutine btp_wet_dry_mom

    subroutine btp_wet_dry_mass(qb,qb_df)

        use mod_grid, only: npoin, npoin_q

        use mod_initial, only: pbprime, pbprime_df
    
        implicit none
        
        real, intent(inout) :: qb(6,npoin_q)
        real, intent(inout) :: qb_df(6,npoin)

        integer :: I, Iq

        ! Compute the wetting and drying treatment for the barotropic

        do I = 1,npoin
            if(qb_df(1,I) < 10250.0) then 
                qb_df(1,I) = 10250.0
                qb_df(2,I) = qb_df(1,I) - pbprime_df(I)
            end if
        end do

        do Iq = 1,npoin_q
            
            if (qb(1,Iq) < 10250.0) then
                qb(1,Iq) = 10250.0
                qb(2,Iq) = qb(1,Iq) - pbprime(Iq)
            end if
        end do
    
    end subroutine btp_wet_dry_mass


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

                    I = indexq(Iq,ip) 

                    q(:,Iq,k) = q(:,Iq,k) + psih(Iq,ip)*q_df(:,I,k)
                    qprime(:,Iq,k) = qprime(:,Iq,k) + psih(Iq,ip)*qprime_df(:,I,k)

                end do
            end do 
        end do 

        do Iq = 1,npoin_q 
            do ip = 1,npts 

                I = indexq(Iq,ip) 

                qb(:,Iq) = qb(:,Iq) + psih(Iq,ip)*qb_df(:,I)
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

    subroutine massinv_rhs(rhs,nvarb)

        use mod_basis, only : npts, nglx, ngly
        use mod_grid, only : npoin, nelem, intma
        use mod_metrics, only : massinv_e

        implicit none

        !global arrays
        integer, intent(in) :: nvarb
        real, intent(inout) :: rhs(nvarb,npoin)

        !local arrays
        real, dimension(nvarb,npts) :: qlocal, rhs_e
        integer :: e, i, j, k, l, m, inodes(npts), ii, Ip

        !loop through elements
        do e = 1,nelem

            !Store local DOFs
            ii = 0
            do j=1,ngly
                do i=1,nglx
                    Ip = intma(i,j,1,e)
                    ii=ii+1
                    inodes(ii) = Ip
                    !Store Primitive Variables

                    qlocal(:,ii) = rhs(:,Ip)

                end do !i
            end do !j

            ! Multiply the rhs for each element by element inverse mass matrix
            do m = 1, nvarb
                rhs_e(m,:) = matmul(massinv_e,qlocal(m,:))
            end do

            ! Store the element rhs into the global rhs
            do l = 1,npts
                Ip = inodes(l)
                rhs(:,Ip) = rhs_e(:,l) 
            end do !ii
        end do
    end subroutine massinv_rhs

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

                I = indexq(Iq,ip)
                dhdx = dpsidx(Iq,ip)
                dhdy = dpsidy(Iq,ip)

                grad_uv(1,1,Iq) = grad_uv(1,1,Iq) + dhdx*uv(1,I)
                grad_uv(1,2,Iq) = grad_uv(1,2,Iq) + dhdy*uv(1,I)

                grad_uv(2,1,Iq) = grad_uv(2,1,Iq) + dhdx*uv(2,I)
                grad_uv(2,2,Iq) = grad_uv(2,2,Iq) + dhdy*uv(2,I)

            end do
        end do

    end subroutine compute_gradient_uv

    subroutine compute_gradient_uv_v1(grad_uv,uv)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad
        use mod_basis, only: nq, psi, psix, psiy, psiz, dpsi, dpsix, dpsiy, dpsiz
        use mod_metrics, only: ksi_x,ksi_y,ksi_z, eta_x,eta_y,eta_z, zeta_x,zeta_y,zeta_z, jac

        use mod_initial, only: dpsidx_df,dpsidy_df, index_df

        implicit none

        real, dimension(2,npoin), intent(in) :: uv
        real, dimension(2,2,npoin), intent(out) :: grad_uv
        integer :: Iq, I, ip, e, iquad, jquad, n, m
        real :: dhdx, dhdy, hi
        real :: h_e, h_n
        real :: e_x, e_y, n_x, n_y
        
        grad_uv = 0.0
        
        !Construct Volume Integral Contributions

        do e = 1,nelem

            do jquad = 1,ngly
                do iquad = 1,nglx
                    
                    Iq = intma(iquad,jquad,1,e)

                    e_x = ksi_x(iquad,jquad,1,e); e_y = ksi_y(iquad,jquad,1,e); 
                    n_x = eta_x(iquad,jquad,1,e); n_y = eta_y(iquad,jquad,1,e);

                    ip = 0
                    
                    do m = 1, ngly 
                        do n = 1, nglx 

                            I = intma(n,m,1,e)
                
                            ! Xi derivatives
                            h_e = dpsix(n, iquad) * psiy(m, jquad)

                            ! Eta derivatives
                            h_n = psix(n, iquad) * dpsiy(m, jquad)
                
                            ! Pressure terms
                            dhdx = h_e * e_x + h_n * n_x
                            dhdy = h_e * e_y + h_n * n_y

                            grad_uv(1,1,Iq) = grad_uv(1,1,Iq) + dhdx*uv(1,I)
                            grad_uv(1,2,Iq) = grad_uv(1,2,Iq) + dhdy*uv(1,I)

                            grad_uv(2,1,Iq) = grad_uv(2,1,Iq) + dhdx*uv(2,I)
                            grad_uv(2,2,Iq) = grad_uv(2,2,Iq) + dhdy*uv(2,I)

                        end do !n
                    end do !m

                end do !iquad
            end do !jquad

        end do !e

    end subroutine compute_gradient_uv_v1

    subroutine compute_gradient_uv_v2(grad_uv,uv)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad
        use mod_basis, only: nq, psiq, psiqx, psiqy, psiqz, dpsiq, dpsiqx, dpsiqy, dpsiqz
        use mod_metrics, only: ksiq_x,ksiq_y,ksiq_z, etaq_x,etaq_y,etaq_z, zetaq_x,zetaq_y,zetaq_z

        use mod_initial, only: dpsidx_df,dpsidy_df, index_df

        implicit none

        real, dimension(2,npoin), intent(in) :: uv
        real, dimension(2,2,npoin), intent(out) :: grad_uv
        integer :: Iq, I, ip
        real :: dhdx, dhdy
        
        grad_uv = 0.0
        
        !Construct Volume Integral Contributions

        do Iq = 1, npoin

            do ip = 1,npts

                I = index_df(Iq,ip)
                dhdx = dpsidx_df(Iq,ip)
                dhdy = dpsidy_df(Iq,ip)

                grad_uv(1,1,Iq) = grad_uv(1,1,Iq) + dhdx*uv(1,I)
                grad_uv(1,2,Iq) = grad_uv(1,2,Iq) + dhdy*uv(1,I)

                grad_uv(2,1,Iq) = grad_uv(2,1,Iq) + dhdx*uv(2,I)
                grad_uv(2,2,Iq) = grad_uv(2,2,Iq) + dhdy*uv(2,I)

            end do
        end do

    end subroutine compute_gradient_uv_v2

end module mod_barotropic_terms
            