
module mod_splitting

    ! ===========================================================================================================================
    ! This module contains the routines for the barotropic-baroclinic splitting
    !   Author: Yao Gahounzo 
    !   Date: March 27, 2023
    ! It contains the following routines:
    ! - ti_barotropic: barotropic substem for splitting system using two-level time integration (see Higdon et al. 2005)
    ! - thickness: baroclinic substem for splitting system using two-level time integration 
    ! - momentum: baroclinic substem for splitting system using two-level time integration
    !
    ! ===========================================================================================================================

    use mod_initial, only: coriolis_df, coriolis_quad, tau_wind, pbprime, pbprime_face
    use mod_initial, only: N_btp, coriolis_quad, pbprime_df, fdt_btp, fdt2_btp, a_btp, b_btp, tau_wind, alpha_mlswe
    use mod_grid, only: npoin, npoin_q, nface, face, intma_dg_quad, intma
    use mod_basis, only: nqx, nqy, nqz, nq, ngl
    use mod_input, only: nlayers, dt_btp, ifilter, nlayers, dt, ad_mlswe
    use mod_create_rhs_mlswe, only: create_rhs_btp_momentum, btp_mass_advection_rhs, create_rhs_btp_momentum_new, &
                                    layer_mass_advection_rhs, layer_mass_advection_rhs1
    use mod_barotropic_terms, only: btp_evaluate_mom, btp_evaluate_mom_face, btp_evaluate_pb, btp_evaluate_pb_face, &
                                    btp_mass_advection_terms, btp_bcl_coeffs, btp_evaluate_mom_dp, btp_evaluate_mom_dp_face

    use mod_input, only: method_visc
    use mod_barotropic_terms, only: btp_mom_boundary_df, compute_btp_terms, compute_btp_mom_terms
    use mod_layer_terms, only: filter_mlswe, layer_mass_advection_terms, evaluate_dp, evaluate_dp_face, consistency_mass_terms1
    use mod_laplacian_quad, only: create_laplacian_mlswe_v3, create_laplacian_mlswe_v4

    use mod_initial, only: zbot_face, fdt_bcl, fdt2_bcl, a_bcl, b_bcl,  coriolis_quad, a_bclp, b_bclp
    use mod_create_rhs_mlswe, only: layer_momentum_rhs, interpolate_layer_from_quad_to_node, rhs_layer_shear_stress
    use mod_layer_terms, only: compute_momentum_edge_values, layer_momentum_advec_terms, layer_pressure_terms, layer_windbot_stress_terms, &
                                shear_stress_system, layer_mom_boundary_df, &
                                filter_mlswe, evaluate_mom, velocity_df, evaluate_mom_face

    use mod_layer_terms, only: bcl_wet_dry_mom_df,bcl_wet_dry_mom
    use mod_face, only: imapl_q, imapr_q, normal_vector_q, imapl
    use mod_laplacian_quad, only: create_laplacian_mlswe_layer_v3
        
    implicit none

    public :: ti_barotropic, thickness, momentum, momentum_mass
    
    contains

    subroutine ti_barotropic(one_plus_eta_out,one_plus_eta_edge_2_ave,uvb_ave, ope_ave, ope2_ave, &
        H_ave, Qu_ave, Qv_ave, Quv_ave, btp_mass_flux_ave, ope_ave_df, uvb_face_ave, &
        ope_face_ave, btp_mass_flux_face_ave, H_face_ave, &
        Qu_face_ave, Qv_face_ave, Quv_face_ave, tau_wind_ave, tau_bot_ave,qb,qb_face,&
        qb_df, qprime,qprime_face, qprime_df, flag_pred, uvb_df_ave)

        ! ===========================================================================================================================
        ! This subroutine predicts and corrects the barotropic quantities for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the barotropic variables (pb_df(or pbpert), pbub_df, pbvb_df) are stored in qb_df
        ! The quadrature points of the barotropic variables (pb, pbpert, ub, vb, pbub, pbvb) are stored in qb. pb = pbprime + pb_df
        ! The face values of the barotropic variables are stored in qb_face
        ! The quadrature points variables of(dpprime, uprime, vprime) are stored in qprime
        !
        ! During the prediction step, use the baroclinic quantities from baroclinic time level  n. 
        ! During the correction step, use averages of the predicted values and values from baroclinic time level n.
        ! 
        ! Return the next baroclinic time step values of the barotropic variables and the barotropic substep averages.
        ! ===========================================================================================================================

        implicit none

        real, dimension(4,npoin_q), intent(inout) :: qb
        real, dimension(4,2,nq,nface), intent(inout) :: qb_face
        real, dimension(4,npoin), intent(inout) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface, nlayers), intent(in) :: qprime_face
        real, dimension(nq,nface), intent(out) :: one_plus_eta_edge_2_ave
        real, dimension(npoin_q), intent(out) :: ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave 
        real, dimension(2,npoin_q), intent(out) :: btp_mass_flux_ave, uvb_ave
        real, dimension(npoin), intent(out) :: ope_ave_df
        real, dimension(2,2,nq,nface), intent(out) :: uvb_face_ave
        real, dimension(2,nq,nface), intent(out) :: btp_mass_flux_face_ave, ope_face_ave
        real, dimension(nq,nface), intent(out) :: H_face_ave
        real, dimension(2,nq,nface), intent(out) :: Qu_face_ave, Qv_face_ave, Quv_face_ave
        real, dimension(npoin), intent(out) :: one_plus_eta_out
        real, dimension(2,npoin_q), intent(out) :: tau_wind_ave, tau_bot_ave
        integer, intent(in) :: flag_pred
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(2,npoin), intent(out) :: uvb_df_ave

        real, dimension(npoin_q) :: one_plus_eta
        real, dimension(4,npoin) :: qb_df_pred
        real, dimension(4,npoin_q) :: qb_init
        real, dimension(4,2,nq,nface) :: qb_face_init
        real, dimension(2,npoin) :: rhs_mom, rhs_visc_bcl, rhs_visc_bcl1
        real, dimension(2,nq,nface) :: Qu_face, Qv_face, one_plus_eta_face, flux_edge, H_bcl_edge, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge
        real, dimension(nq,nface) :: H_face, one_plus_eta_edge_2, one_plus_eta_edge
        real, dimension(npoin_q) :: Quu, Qvv, Quv, H, Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, pbprime_visc
        real, dimension(npoin) :: pb_advec, tempu, tempv, one_plus_eta_df
        real, dimension(2,npoin_q) :: btp_mass_flux, tau_bot

        real, dimension(2,nq,nface) :: qb_com1
        real, dimension(2,2,nq,nface) :: qb_com2
        real, dimension(3,npoin) :: rhs_mom1
        real, dimension(2,npoin) :: qb_df_temp, uvb_df

        integer :: mstep, I, Iq, iquad, iface, ilr, k
        real :: N_inv

        one_plus_eta_edge_2_ave = 0.0
        uvb_ave  = 0.0
        ope_ave = 0.0
        btp_mass_flux_ave = 0.0
        H_ave = 0.0
        Qu_ave = 0.0
        Qv_ave = 0.0
        Quv_ave = 0.0
        ope_ave_df = 0.0
        uvb_face_ave  = 0.0
        ope_face_ave = 0.0
        btp_mass_flux_face_ave = 0.0
        H_face_ave = 0.0
        Qu_face_ave = 0.0
        Qv_face_ave = 0.0
        Quv_face_ave = 0.0
        tau_wind_ave = 0.0
        tau_bot_ave = 0.0
        rhs_visc_bcl = 0.0
        ope2_ave = 0.0
        uvb_df_ave = 0.0

        ! Compute baroclinic coefficients in the barotropic momentum fluxes, barotropic pressure forcing, and barotropic
        ! horizontal viscosity terms.  These are needed for the barotropic momentum equation.

        call btp_bcl_coeffs(Q_uu_dp,Q_uv_dp,H_bcl,Q_vv_dp,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge,H_bcl_edge, qprime,qprime_face)

        qb_init = qb
        qb_face_init = qb_face

        ! ========== The time loop begins here ==================

        do mstep = 1, N_btp

            ! ============================================ Prediction steps ==========================================================

            ! Communicate the grid cell face values of the barotropic variables to the neighboring processors.

            call create_communicator_quad(qb_face_init,4)

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta,one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, flux_edge, btp_mass_flux, qb_init,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face_init,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb_df)

            ! Compute RHS viscosity terms

            if(method_visc == 1) then
                call create_laplacian_mlswe_v3(rhs_visc_bcl,qprime,qprime_face,qb_init,qb_face_init)
            elseif(method_visc == 2) then
                uvb_df(1,:) = qb_df(3,:)/qb_df(1,:)
                uvb_df(2,:) = qb_df(4,:)/qb_df(1,:)
                call create_laplacian_mlswe_v4(rhs_visc_bcl,qprime_df,uvb_df,qprime(1,:,:), qprime_face, qb_face)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_momentum_new(rhs_mom1,Quu,Qvv,Quv,H,qb_init,H_face,Qu_face,Qv_face,tau_bot,rhs_visc_bcl,btp_mass_flux, flux_edge,qb_face_init)

            ! Predict barotropic 

            qb_df_pred(3,:) = qb_df(3,:) + fdt_btp(:)*qb_df(4,:) + dt_btp*(rhs_mom1(1,:) + rhs_visc_bcl(1,:))
            qb_df_pred(4,:) = qb_df(4,:) - fdt_btp(:)*qb_df(3,:) + dt_btp*(rhs_mom1(2,:) + rhs_visc_bcl(2,:))

            call btp_mom_boundary_df(qb_df_pred(3:4,:))

            ! ********** Barotropic mass equation ***********

            qb_df_pred(2,:) = qb_df(2,:) + dt_btp * rhs_mom1(3,:)
            qb_df_pred(1,:) = qb_df_pred(2,:) + pbprime_df(:)

            ! Use basis functions and the predicted degrees of freedom to compute predicted values of
            ! pb*ubar, pb*vbar, pbpert, and pb at endpoints and quadrature points in the grid cells.
            ! Then use division to compute ubar and vbar at those points.

            call btp_evaluate_mom_dp(qb,qb_df_pred)
            call btp_evaluate_mom_dp_face(qb_face, qb)

            ! Communication of the faces values within the processors
            call create_communicator_quad(qb_face,4)

            ! Accumulate sums for time averaging
            uvb_ave(1,:) = uvb_ave(1,:) + qb_init(3,:)/qb_init(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb_init(4,:)/qb_init(1,:)
            uvb_face_ave(1,:,:,:) = uvb_face_ave(1,:,:,:) + qb_face_init(3,:,:,:)/qb_face_init(1,:,:,:)
            uvb_face_ave(2,:,:,:) = uvb_face_ave(2,:,:,:) + qb_face_init(4,:,:,:)/qb_face_init(1,:,:,:)

            ope_ave = ope_ave + one_plus_eta
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2
            ope2_ave = ope2_ave + one_plus_eta**2 
            ! ============================================ Correction steps ==========================================================

            ! ******** Barotropic mass equation ***********

            ! Compute the mass advection term, using data from barotropic previous time and predicted data for barotropic level.

            btp_mass_flux = qb(3:4,:)

            call btp_mass_advection_terms(flux_edge,qb_face)

            ! Compute RHS for the barotropic mass equation

            call btp_mass_advection_rhs(pb_advec, btp_mass_flux, flux_edge)

            ! Correct barotropic mass

            qb_df(2,:) = qb_df(2,:) + 0.5*dt_btp * (pb_advec(:) + rhs_mom1(3,:))
            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            ! Evaluate the corrected barotropic mass at the quadrature points and face 

            call btp_evaluate_pb(qb, qb_df)
            call btp_evaluate_pb_face(qb_face, qb)

            ! Communication of the dp values within the processors

            qb_com1(:,:,:) = qb_face(2,:,:,:)

            call create_communicator_quad_1var(qb_face(2,:,:,:))

            ! ******** Barotropic momentum equation ***********

            call compute_btp_mom_terms(Quu,Qvv,Quv,Qu_face,Qv_face, H, H_face,tau_bot, one_plus_eta, one_plus_eta_edge_2, one_plus_eta_df, &
                one_plus_eta_face, qb,Q_uu_dp,Q_uv_dp,Q_vv_dp,qb_face,Q_uu_dp_edge,Q_uv_dp_edge,Q_vv_dp_edge, qprime, &
                H_bcl, H_bcl_edge, qb_df)

            ! Compute RHS viscosity terms

            if(method_visc == 1) then
                call create_laplacian_mlswe_v3(rhs_visc_bcl1,qprime,qprime_face,qb,qb_face)
            elseif(method_visc == 2) then
                uvb_df(1,:) = qb_df_pred(3,:)/qb_df_pred(1,:)
                uvb_df(2,:) = qb_df_pred(4,:)/qb_df_pred(1,:)
                call create_laplacian_mlswe_v4(rhs_visc_bcl1,qprime_df,uvb_df,qprime(1,:,:),qprime_face, qb_face)
            end if

            ! Compute RHS for the barotropic momentum equation

            call create_rhs_btp_momentum(rhs_mom,Quu,Qvv,Quv,H,qb, H_face,Qu_face,Qv_face,tau_bot,rhs_visc_bcl, qb_face)

            rhs_mom = 0.5*(rhs_mom1(1:2,:) + rhs_mom)
            rhs_visc_bcl = 0.5*(rhs_visc_bcl1 + rhs_visc_bcl)

            ! Correct barotropic momentum

            qb_df_temp(1,:) = qb_df(3,:) + dt_btp*(rhs_mom(1,:) + rhs_visc_bcl(1,:))
            qb_df_temp(2,:) = qb_df(4,:) + dt_btp*(rhs_mom(2,:) + rhs_visc_bcl(2,:))

            call btp_mom_boundary_df(qb_df_temp(:,:))

            tempu = qb_df_temp(1,:) + fdt2_btp(:)*qb_df(4,:)
            tempv = qb_df_temp(2,:) - fdt2_btp(:)*qb_df(3,:)

            qb_df(3,:) = a_btp(:) * tempu + b_btp(:) * tempv
            qb_df(4,:) = -b_btp(:) * tempu + a_btp(:) * tempv

            call btp_mom_boundary_df(qb_df(3:4,:))

            if(mstep == N_btp .and. ifilter > 0 .and. flag_pred == 0) then 
                call filter_mlswe(qb_df(3:4,:),2)
                call btp_mom_boundary_df(qb_df(3:4,:))
            end if

            ! Evaluate the corrected barotropic momentum at the quadrature points and face

            call btp_evaluate_mom(qb,qb_df)
            call btp_evaluate_mom_face(qb_face, qb)

            ! Communication of the variables within the processors
            qb_com2(1,:,:,:) = qb_face(3,:,:,:)/qb_face(1,:,:,:)
            qb_com2(2,:,:,:) = qb_face(4,:,:,:)/qb_face(1,:,:,:)
            call create_communicator_quad(qb_com2,2)

            qb_face(2,:,:,:)= qb_com1(:,:,:)

            ! Accumulate sums for time averaging

            uvb_ave(1,:) = uvb_ave(1,:) + qb(3,:)/qb(1,:)
            uvb_ave(2,:) = uvb_ave(2,:) + qb(4,:)/qb(1,:)

            uvb_face_ave = uvb_face_ave + qb_com2(1:2,:,:,:)

            ope_ave = ope_ave + one_plus_eta; tau_bot_ave = tau_bot_ave + tau_bot
            btp_mass_flux_ave = btp_mass_flux_ave + btp_mass_flux
            H_ave = H_ave + H;  Qu_ave = Qu_ave + Quu
            Qv_ave = Qv_ave + Qvv;  Quv_ave = Quv_ave + Quv; ope_ave_df = ope_ave_df + one_plus_eta_df

            ope_face_ave = ope_face_ave + one_plus_eta_face
            btp_mass_flux_face_ave = btp_mass_flux_face_ave + flux_edge
            H_face_ave = H_face_ave + H_face; Qu_face_ave = Qu_face_ave + Qu_face; Qv_face_ave = Qv_face_ave + Qv_face

            one_plus_eta_edge_2_ave = one_plus_eta_edge_2_ave + one_plus_eta_edge_2

            tau_wind_ave = tau_wind_ave + tau_wind
            ope2_ave = ope2_ave + one_plus_eta**2

            uvb_df_ave(1,:) = uvb_df_ave(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_df_ave(2,:) = uvb_df_ave(2,:) + qb_df(4,:)/qb_df(1,:)

            qb_init = qb
            qb_face_init = qb_face

        end do

        ! Update the barotropic variables at the baroclinic time level n+1

        one_plus_eta_out = one_plus_eta_df

        N_inv = 1.0 / real(2*N_btp)

        uvb_ave = N_inv*uvb_ave

        ope_ave = N_inv*ope_ave
        H_ave = N_inv*H_ave 
        Qu_ave = N_inv*Qu_ave
        Qv_ave = N_inv*Qv_ave 
        Quv_ave = N_inv*Quv_ave 
        btp_mass_flux_ave = N_inv*btp_mass_flux_ave

        uvb_face_ave = N_inv*uvb_face_ave 

        ope_face_ave = N_inv*ope_face_ave 
        H_face_ave = N_inv*H_face_ave 
        Qu_face_ave = N_inv*Qu_face_ave 
        Qv_face_ave = N_inv*Qv_face_ave 
        Quv_face_ave = N_inv*Quv_face_ave 
        btp_mass_flux_face_ave = N_inv*btp_mass_flux_face_ave 

        one_plus_eta_edge_2_ave = N_inv*one_plus_eta_edge_2_ave 
        ope2_ave = ope2_ave*N_inv

        ope_ave_df = N_inv*ope_ave_df
        uvb_df_ave = N_inv*uvb_df_ave

        tau_wind_ave = tau_wind_ave / real(N_btp)
        tau_bot_ave = tau_bot_ave / real(N_btp)

    end subroutine ti_barotropic

    subroutine thickness(q,qprime, q_df, q_face, qprime_face, u_edge, v_edge, uvb_ave, btp_mass_flux_ave, ope_ave, uvb_face_ave, ope_face_ave, &
        btp_mass_flux_face_ave,dpprime_df, flag_pred, qb_df)

        ! ===========================================================================================================================
        ! This subroutine is used to predict or correct the layer thickness for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the layer thickness dpprime_df is stored in q_df(1,:,:)
        ! The quadrature points of the layer thickness dpprime is stored in qprime(1,:,:)
        ! The face values of the layer thickness dpprime_face is stored in qprime_face(1,:,:,:,:)
        ! The quadrature points of the layer thickness dp is stored in qb(1,:,:)
        ! The face values of the layer thickness dp_face is stored in qb_face(1,:,:,:,:)
        !
        ! Enforce consistency between the layer masses and the barotropic mass.
        ! ===========================================================================================================================

        implicit none
    
        ! Input variables
        real, intent(inout) :: qprime(3,npoin_q,nlayers), q(3,npoin_q,nlayers)
        real, intent(inout) :: q_df(3,npoin,nlayers), qprime_face(3,2,nq,nface,nlayers), q_face(3,2,nq,nface,nlayers)
        real, intent(in)    :: uvb_ave(2,npoin_q), ope_ave(npoin_q), btp_mass_flux_ave(2,npoin_q)

        real, intent(in)    :: uvb_face_ave(2,2,nq,nface), ope_face_ave(2,nq,nface)
        real, intent(in)    :: btp_mass_flux_face_ave(2,nq,nface)
        integer, intent(in) :: flag_pred
        real, intent(in)    :: qb_df(4,npoin)
    
        ! Output variables
        real, dimension(2,nq, nface, nlayers), intent(out)   :: u_edge, v_edge
        real, intent(out) :: dpprime_df(npoin,nlayers)
    
        ! Other variables
        real :: sum_layer_mass_flux(2,npoin_q), uvdp_temp(2,npoin_q,nlayers), flux_edge(2,nq,nface,nlayers), flux_adjustment(2,npoin_q,nlayers)
        real :: flux_adjust_edge(2,nq,nface,nlayers), one_plus_eta_temp(npoin), sum_layer_mass_flux_face(2,nq,nface), dp_advec(npoin,nlayers)
        real :: one_plus_eta_temp1
        integer :: k, I, iquad, iface, ilr,Iq
        real :: adjust_mass(npoin)
        real :: disp(nq,nface,nlayers)

        ! =========================================== layer mass ===============================================================
    
        ! Compute the mass advection term.   

        call layer_mass_advection_terms(sum_layer_mass_flux,sum_layer_mass_flux_face,u_edge,v_edge,uvdp_temp,flux_edge, &
                qprime,uvb_ave,ope_ave,qprime_face,uvb_face_ave,ope_face_ave, disp)

        ! Compute RHS for layer thickness, return the result in array  dp_advec .

        call layer_mass_advection_rhs1(dp_advec, uvdp_temp, flux_edge, q_face, disp)

        ! Compute the tentative values of the predicted or corrected 
        ! degrees of freedom for  dp.  These would be the values at
        ! baroclinic time level  n+1,  except for the need to enforce
        ! consistency between the layer masses and the barotropic mass 
        ! and for the possible need to use a variation limiter to 
        ! prevent pointwise negative layer thicknesses. 
        
        do k = 1, nlayers

            q_df(1,:,k) = q_df(1,:,k) + dt*dp_advec(:,k)

            ! Check for negative layer thicknesses.  If any are found, stop the program.
            if(any(q_df(1, :, k) < 0.0)) then
                write(*,*) 'Negative mass in thickness at some points'
                stop
            endif
        end do

 
        adjust_mass(:) = (1.0/real(nlayers))*(qb_df(1,:) - sum(q_df(1,:,:),dim=2))

        ! Apply filter to the thickness
        do k = 1,nlayers
            q_df(1,:,k) = q_df(1,:,k) + adjust_mass(:)
            
            if(ifilter > 0 .and. flag_pred == 0) then 
                call filter_mlswe(q_df(1,:,k),1)
            end if
        end do

        ! Use the adjusted degrees of freedom q_df(1,:,:) to compute revised values of  dp  and  dp' at cell edges and quadrature points.

        call evaluate_dp(q,qprime,q_df, pbprime)
        call evaluate_dp_face(q_face, qprime_face,q, qprime)

        ! Store the degree of freedom (nodal points) values of dpprime_df ( or q_df(1,:,:)) in the array dpprime_df for use in layer_pressure_terms.

        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)

        do k = 1,nlayers
            dpprime_df(:,k) = q_df(1,:,k) / one_plus_eta_temp(:)
        end do
        
    end subroutine thickness


    subroutine momentum(q,qprime,q_df,q_face,qprime_face,qb,qb_face,ope_ave,one_plus_eta_edge_2_ave,uvb_ave,u_edge,v_edge,Qu_ave,&
        Qv_ave,Quv_ave,H_ave,uvb_face_ave,ope_face_ave,Qu_face_ave,Qv_face_ave,Quv_face_ave,H_face_ave,&
        qprime_df,ope_ave_df,tau_bot_ave,tau_wind_ave,qprime_face2,flag_pred,&
        q2,q_face2,qb_df,qprime2, qprime_face3, ope2_ave, uvb_df_ave)

        ! ===========================================================================================================================
        ! This subroutine is used to predict or correct the layer momentum for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the layer momentum is stored in q_df(2:3,:,:)
        ! The quadrature points of the layer momentum are stored in qprime(2:3,:,:)
        ! The face values of the layer momentum are stored in qprime_face(2:3,:,:,:,:)
        ! The quadrature points of the layer momentum are stored in qb(2:5,:,:)
        ! The face values of the layer momentum are stored in qb_face(2:5,:,:,:)
        ! ===========================================================================================================================

        implicit none

        ! Input variables
        real, dimension(3,npoin_q,nlayers), intent(inout) :: q
        real, dimension(3,npoin_q,nlayers), intent(inout) :: qprime
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(3,2,nq,nface,nlayers), intent(inout) :: q_face
        real, dimension(3,2,nq,nface,nlayers), intent(inout) :: qprime_face, qprime_face2, qprime_face3
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(4,2,nq,nface), intent(in) :: qb_face
        
        real, dimension(npoin_q), intent(in) :: ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave
        real, dimension(2,2,nq,nface), intent(in) :: uvb_face_ave
        real, dimension(nq,nface), intent(in) :: H_face_ave, one_plus_eta_edge_2_ave
        real, dimension(2, npoin_q), intent(in) :: tau_wind_ave, tau_bot_ave, uvb_ave
        real, dimension(2, nq, nface), intent(in) :: Qu_face_ave, Qv_face_ave, Quv_face_ave, ope_face_ave
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(npoin), intent(in) :: ope_ave_df
        integer, intent(in) :: flag_pred
        real, dimension(2,nq, nface, nlayers), intent(in)   :: u_edge, v_edge
        real, dimension(3,npoin_q,nlayers), intent(in) :: q2
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: q_face2
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime2
        real, dimension(2,npoin), intent(in) :: uvb_df_ave

        ! Local variables

        real, dimension(2,npoin_q,nlayers+1)   :: grad_z
        real, dimension(2,npoin_q,nlayers)     :: coriolis, uvdp,uv
        real, dimension(2,npoin,nlayers)       :: q_df_temp
        real, dimension(npoin)                 :: fdt, fdt2, a, b, zk
        real, dimension(nq,nface,nlayers)      :: udp_left, vdp_left, udp_right, vdp_right, disp
        real, dimension(npoin_q, nlayers)      :: u_udp_temp, v_vdp_temp
        real, dimension(2, npoin_q, nlayers)   :: u_vdp_temp
        real, dimension(2, nq, nface, nlayers) :: udp_flux_edge, vdp_flux_edge
        real, dimension(npoin_q,nlayers)       :: H_r
        real, dimension(2,nq,nface,nlayers)    :: H_r_face
        real, dimension(2,2,nq,nface,nlayers)    :: uv_face, uvdp_face
        real, dimension(npoin_q,nlayers+1)     :: p
        real, dimension(npoin,nlayers+1)       :: z_elev
        real, dimension(2,npoin_q,nlayers+1)     :: tau_wind_int, tau_bot_int
        real, dimension(2,npoin,nlayers) :: rhs_mom, rhs_stress, rhs_visc_bcl, rhs_c
        real, dimension(3,npoin,nlayers) :: q_df3
        
        real, dimension(2,npoin_q) :: z_xy, qpvisc
        real, dimension(npoin_q,nlayers) :: dp
        integer :: k, Iq, iquad, I, iface, ilr, n
        real :: tempu, tempv
        real, dimension(npoin) :: tempu1, tempv1, ub_df_ave, vb_df_ave
        real :: nx, ny, un, E(2,npoin), uvb_ave_face(2,2,nq,nface)
        integer :: el, er, il, jl, kl, kr, jr, ir
        real, dimension(2,2,npoin_q) :: grad_qpvisc
        real, dimension(2,2,nq,nface) :: qpvisc_face
        real, dimension(4,2,nq,nface) :: grad_qpvisc_face
        real, dimension(2,npoin_q) :: grad_ub, grad_vb
        real, dimension(2,npoin) :: rhs_visc_bcl_temp, q_vic
        real, dimension(3,npoin,nlayers) :: qprime_df2

        grad_z = 0.0

        q_df_temp = 0.0
        coriolis = 0.0

        u_udp_temp = 0.0
        v_vdp_temp = 0.0
        u_vdp_temp = 0.0
        udp_flux_edge = 0.0
        vdp_flux_edge = 0.0

        rhs_visc_bcl = 0.0

        ! ==== layer momentum ====

        ! Compute the pressure terms

        call layer_pressure_terms(H_r, H_r_face, p, z_elev, qprime, qprime_face, ope_ave, H_ave, ope_face_ave, zbot_face, H_face_ave, qprime_df, one_plus_eta_edge_2_ave, ope_ave_df, grad_z, ope2_ave)

        call compute_momentum_edge_values(udp_left, vdp_left, udp_right, vdp_right, qprime_face3, uvb_face_ave, ope_face_ave, disp)

        call layer_momentum_advec_terms(u_udp_temp, u_vdp_temp, v_vdp_temp, udp_flux_edge, vdp_flux_edge, &
            q2, qprime, uvb_ave, ope_ave, u_edge, v_edge, Qu_ave, Qv_ave, Quv_ave, &
            udp_left, vdp_left, udp_right, vdp_right, Qu_face_ave, Qv_face_ave, Quv_face_ave)

        ! Compute the wind stress terms

        call layer_windbot_stress_terms(tau_wind_int, tau_bot_int,qprime, tau_bot_ave, tau_wind_ave)

        ! Compute the RHS viscosity terms
        if(method_visc > 0) then 
            call create_laplacian_mlswe_layer_v3(rhs_visc_bcl,qprime2,qprime_face2,uvb_ave,uvb_face_ave)
        end if

        ! Compute the RHS of the layer momentum equation

        call layer_momentum_rhs(rhs_mom, H_r,p,grad_z,u_udp_temp,u_vdp_temp,v_vdp_temp,udp_flux_edge,vdp_flux_edge,H_r_face,tau_wind_int,tau_bot_int,rhs_visc_bcl, q_face, disp)

        ! Compute the momentum equation variables for the next time step

        do k = 1,nlayers
            q_df_temp(1,:,k) = q_df(2,:,k) + dt*(rhs_mom(1,:,k) + rhs_visc_bcl(1,:,k))
            q_df_temp(2,:,k) = q_df(3,:,k) + dt*(rhs_mom(2,:,k) + rhs_visc_bcl(2,:,k))
        end do

        ! Compute the shear stress terms
        ! if(ad_mlswe > 0.0) then 
            
        !     do k = 1,nlayers
        !         do I = 1,npoin
        !             tempu1(I) = q_df_temp(1,I,k) + fdt2_bcl(I)*q_df(5,I,k)
        !             tempv1(I) = q_df_temp(2,I,k) - fdt2_bcl(I)*q_df(4,I,k)

        !             q_df3(4,I,k) = a_bcl(I)*tempu1(I) + b_bcl(I)*tempv1(I)
        !             q_df3(5,I,k) = - b_bcl(I)*tempu1(I) + a_bcl(I)*tempv1(I)
        !         end do
        !         q_df3(1,:,k) = q_df(1,:,k)
        !     end do

        !     call layer_mom_boundary_df(q_df3(4:5,:,:))

        !     ! Extract velocity from the momentum
        !     call velocity_df(q_df3, qb_df)
        !     ! Evaluate velocity and momentum at the quad points
        !     call evaluate_mom(q,q_df3)

        !     ! Compute the vertical stress terms

        !     call shear_stress_system(q)

        !     call rhs_layer_shear_stress(rhs_stress,q)

        !     do k = 1,nlayers
        !         q_df_temp(1,:,k) = q_df_temp(1,:,k) + dt*rhs_stress(1,:,k)
        !         q_df_temp(2,:,k) = q_df_temp(2,:,k) + dt*rhs_stress(2,:,k)
        !     end do

        ! end if ! ad_mlswe > 0.0
        
        ! Add the Coriolis term

        do k = 1,nlayers

            tempu1(:) = q_df_temp(1,:,k) + fdt2_bcl(:)*q_df(3,:,k)
            tempv1(:) = q_df_temp(2,:,k) - fdt2_bcl(:)*q_df(2,:,k)

            q_df(2,:,k) = a_bcl(:)*tempu1(:) + b_bcl(:)*tempv1(:)
            q_df(3,:,k) = - b_bcl(:)*tempu1(:) + a_bcl(:)*tempv1(:)
        end do

        call layer_mom_boundary_df(q_df(2:3,:,:))

        ! Extract velocity from the momentum
        call velocity_df(q_df, qb_df, flag_pred)
        ! Evaluate velocity and momentum at the quad points
        call evaluate_mom(q,q_df)
        ! Extract faces values
        call evaluate_mom_face(q_face, q)

        ! Compute uprime and vprime at the quad and nodal points

        do k = 1,nlayers

            qprime(2,:,k) = q(2,:,k)/q(1,:,k) - qb(3,:)/qb(1,:)
            qprime(3,:,k) = q(3,:,k)/q(1,:,k) - qb(4,:)/qb(1,:)
            qprime_df(2,:,k) = q_df(2,:,k)/q_df(1,:,k) - qb_df(3,:)/qb_df(1,:)
            qprime_df(3,:,k) = q_df(3,:,k)/q_df(1,:,k) - qb_df(4,:)/qb_df(1,:)

            qprime_face(2,:,:,:,k) = q_face(2,:,:,:,k)/q_face(1,:,:,:,k) - qb_face(3,:,:,:)/qb_face(1,:,:,:)
            qprime_face(3,:,:,:,k) = q_face(3,:,:,:,k)/q_face(1,:,:,:,k) - qb_face(4,:,:,:)/qb_face(1,:,:,:)
        end do 

    end subroutine momentum


    subroutine momentum_mass(q,qprime,q_df,q_face,qprime_face,qb,qb_face,ope_ave,one_plus_eta_edge_2_ave,uvb_ave,Qu_ave,&
        Qv_ave,Quv_ave,H_ave,uvb_face_ave,ope_face_ave,Qu_face_ave,Qv_face_ave,Quv_face_ave,H_face_ave,&
        qprime_df,ope_ave_df,tau_bot_ave,tau_wind_ave,qprime_face2,flag_pred,&
        q2,q_face2,qb_df,qprime2, qprime_face3, ope2_ave, btp_mass_flux_ave, btp_mass_flux_face_ave, uvb_df_ave)

        ! ===========================================================================================================================
        ! This subroutine is used to predict or correct the layer momentum for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the layer momentum is stored in q_df(2:3,:,:)
        ! The quadrature points of the layer momentum are stored in qprime(2:3,:,:)
        ! The face values of the layer momentum are stored in qprime_face(2:3,:,:,:,:)
        ! The quadrature points of the layer momentum are stored in qb(2:5,:,:)
        ! The face values of the layer momentum are stored in qb_face(2:5,:,:,:)
        ! ===========================================================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma_dg_quad, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, dt, ad_mlswe, explt_coriolis, method_visc, ifilter, is_mlswe_linear
        use mod_initial, only: coriolis_df, zbot_face, fdt_bcl, fdt2_bcl, a_bcl, b_bcl,  coriolis_quad, a_bclp, b_bclp, pbprime_df, pbprime
        use mod_create_rhs_mlswe, only: layer_momentum_rhs, interpolate_layer_from_quad_to_node, rhs_layer_shear_stress
        use mod_layer_terms, only: compute_momentum_edge_values, layer_momentum_advec_terms, layer_pressure_terms, layer_windbot_stress_terms, &
                                    shear_stress_system, layer_mom_boundary_df, &
                                    filter_mlswe, evaluate_mom, velocity_df, &
                                    evaluate_mom_face, evaluate_dp, evaluate_dp_face

        use mod_layer_terms, only: bcl_wet_dry_mom_df,bcl_wet_dry_mom
        use mod_face, only: imapl_q, imapr_q, normal_vector_q, imapl
        use mod_laplacian_quad, only: create_laplacian_mlswe_layer_v3

        implicit none

        ! Input variables
        real, dimension(3,npoin_q,nlayers), intent(inout) :: q
        real, dimension(3,npoin_q,nlayers), intent(inout) :: qprime
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(3,2,nq,nface,nlayers), intent(inout) :: q_face
        real, dimension(3,2,nq,nface,nlayers), intent(inout) :: qprime_face
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face2, qprime_face3
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(4,2,nq,nface), intent(in) :: qb_face
        
        real, dimension(npoin_q), intent(in) :: ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave
        real, dimension(2,2,nq,nface), intent(in) :: uvb_face_ave
        real, dimension(nq,nface), intent(in) :: H_face_ave, one_plus_eta_edge_2_ave
        real, dimension(2, npoin_q), intent(in) :: tau_wind_ave, tau_bot_ave, uvb_ave
        real, dimension(2, nq, nface), intent(in) :: Qu_face_ave, Qv_face_ave, Quv_face_ave, ope_face_ave
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(npoin), intent(in) :: ope_ave_df
        integer, intent(in) :: flag_pred
        real, dimension(2,npoin), intent(in) :: uvb_df_ave
        
        real, dimension(3,npoin_q,nlayers), intent(in) :: q2
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: q_face2
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime2
        real, intent(in)    :: btp_mass_flux_ave(2,npoin_q)
        real, intent(in)    :: btp_mass_flux_face_ave(2,nq,nface)

        ! Local variables

        real, dimension(2,nq, nface, nlayers):: u_edge, v_edge
        real, dimension(npoin,nlayers) :: dp_advec, dpprime_df
        real, dimension(2,npoin,nlayers)       :: q_df_temp
        real, dimension(2,npoin,nlayers) :: rhs_mom, rhs_visc_bcl, rhs_stress
        real, dimension(npoin) :: one_plus_eta_temp, adjust_mass

        real, dimension(3,npoin_q,nlayers) :: qprime_temp 
        real, dimension(3,2,nq,nface,nlayers) :: qprime_face_temp 
        real, dimension(3,npoin,nlayers) :: q_df3
        real, dimension(2,npoin_q,nlayers) :: uv
        
        integer :: k, Iq, iquad, I, iface
        real :: tempu, tempv
        real, dimension(npoin) :: tempu1, tempv1
        real, dimension(2,nq,nface,nlayers) :: dprime_face_corr
        real :: sum_layer_mass_flux_face(2,nq,nface), sum_layer_mass_flux(2,npoin_q)

        q_df_temp = 0.0

        ! =========================================== layer mass ===============================================================
    
        call rhs_thickness(dp_advec, sum_layer_mass_flux, sum_layer_mass_flux_face, u_edge, v_edge, qprime, qprime_face, uvb_ave, ope_ave, uvb_face_ave, ope_face_ave, q_face)

        do k = 1, nlayers

            q_df(1,:,k) = q_df(1,:,k) + dt*dp_advec(:,k)

            ! Check for negative layer thicknesses.  If any are found, stop the program.
            if(any(q_df(1, :, k) < 0.0)) then
                write(*,*) 'Negative mass in thickness at some points'
                stop
            endif

        end do

        adjust_mass(:) = (1.0/real(nlayers))*(qb_df(1,:) - sum(q_df(1,:,:),dim=2))

        !call apply_consistency(dp_advec, q_df, btp_mass_flux_ave, ope_face_ave, &
        !    btp_mass_flux_face_ave, sum_layer_mass_flux, sum_layer_mass_flux_face)

        ! Apply filter to the thickness
        do k = 1,nlayers
            q_df(1,:,k) = q_df(1,:,k) + adjust_mass(:)
            
            if(flag_pred == 0 .and. ifilter > 0) then 
                call filter_mlswe(q_df(1,:,k),1)
            end if
        end do
        
        ! Use the adjusted degrees of freedom q_df(1,:,:) to compute revised values of  dp  and  dp' at cell edges and quadrature points.

        call evaluate_dp(q,qprime_temp,q_df, pbprime)
        call evaluate_dp_face(q_face, qprime_face_temp,q, qprime_temp)

        ! Store the degree of freedom (nodal points) values of dpprime_df ( or q_df(1,:,:)) in the array dpprime_df for use in layer_pressure_terms.

        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)

        do k = 1,nlayers
            dpprime_df(:,k) = q_df(1,:,k) / one_plus_eta_temp(:)
        end do

        ! If correction step, compute the average values of the layer thickness with values from time step n and and the prediction step
        if(flag_pred == 0) then 

            dprime_face_corr = qprime_face_temp(1,:,:,:,:)
            ! Communication of qprime_face values within the processor boundary
            call bcl_create_communicator(qprime_face_temp(1,:,:,:,:),1,nlayers,nq)

            ! Compute average values 
            qprime_df(1,:,:) = 0.5*(qprime_df(1,:,:) + dpprime_df(:,:))
            qprime(1,:,:) = 0.5*(qprime(1,:,:) + qprime_temp(1,:,:))
            qprime_face(1,:,:,:,:) = 0.5*(qprime_face(1,:,:,:,:) + qprime_face_temp(1,:,:,:,:))
        end if

        ! ==================================== layer momentum ==========================

        call rhs_momentum(rhs_mom, rhs_visc_bcl, qprime,q_face,qprime_face,ope_ave,one_plus_eta_edge_2_ave,uvb_ave,u_edge,v_edge,Qu_ave,&
            Qv_ave,Quv_ave,H_ave,uvb_face_ave,ope_face_ave,Qu_face_ave,Qv_face_ave,Quv_face_ave,H_face_ave,&
            qprime_df,ope_ave_df,tau_bot_ave,tau_wind_ave,qprime_face2,&
            q,qprime2, qprime_face, ope2_ave, uvb_df_ave)

        ! Compute the momentum equation variables for the next time step

        do k = 1,nlayers
            q_df_temp(1,:,k) = q_df(2,:,k) + dt*(rhs_mom(1,:,k) + rhs_visc_bcl(1,:,k))
            q_df_temp(2,:,k) = q_df(3,:,k) + dt*(rhs_mom(2,:,k) + rhs_visc_bcl(2,:,k))
        end do

        ! Compute the shear stress terms
        if(ad_mlswe > 0.0) then 
            
            do k = 1,nlayers
                do I = 1,npoin
                    tempu1(I) = q_df_temp(1,I,k) + fdt2_bcl(I)*q_df(3,I,k)
                    tempv1(I) = q_df_temp(2,I,k) - fdt2_bcl(I)*q_df(2,I,k)

                    q_df3(2,I,k) = a_bcl(I)*tempu1(I) + b_bcl(I)*tempv1(I)
                    q_df3(3,I,k) = - b_bcl(I)*tempu1(I) + a_bcl(I)*tempv1(I)
                end do
                q_df3(1,:,k) = q_df(1,:,k)
            end do

            !call layer_mom_boundary_df(q_df3(2:3,:,:))

            ! Extract velocity from the momentum
            call velocity_df(q_df3, qb_df, flag_pred)
            ! Evaluate velocity and momentum at the quad points
            call evaluate_mom(q,q_df3)

            ! Compute the vertical stress terms

            call shear_stress_system(uv,q)

            call rhs_layer_shear_stress(rhs_stress,uv)

            do k = 1,nlayers
                q_df_temp(1,:,k) = q_df_temp(1,:,k) + dt*rhs_stress(1,:,k)
                q_df_temp(2,:,k) = q_df_temp(2,:,k) + dt*rhs_stress(2,:,k)
            end do

        end if ! ad_mlswe >
        
        ! Add the Coriolis term

        do k = 1,nlayers

            tempu1(:) = q_df_temp(1,:,k) + fdt2_bcl(:)*q_df(3,:,k)
            tempv1(:) = q_df_temp(2,:,k) - fdt2_bcl(:)*q_df(2,:,k)

            q_df(2,:,k) = a_bcl(:)*tempu1(:) + b_bcl(:)*tempv1(:)
            q_df(3,:,k) = - b_bcl(:)*tempu1(:) + a_bcl(:)*tempv1(:)
        end do

        call layer_mom_boundary_df(q_df(2:3,:,:))

        ! Extract velocity from the momentum
        call velocity_df(q_df, qb_df, flag_pred)
        ! Evaluate velocity and momentum at the quad points
        call evaluate_mom(q,q_df)
        ! Extract faces values
        call evaluate_mom_face(q_face, q)

        ! Compute uprime and vprime at the quad and nodal points

        do k = 1,nlayers
            qprime(2,:,k) = q(2,:,k)/q(1,:,k) - qb(3,:)/qb(1,:)
            qprime(3,:,k) = q(3,:,k)/q(1,:,k) - qb(4,:)/qb(1,:)
            qprime_df(2,:,k) = q_df(2,:,k)/q_df(1,:,k) - qb_df(3,:)/qb_df(1,:)
            qprime_df(3,:,k) = q_df(3,:,k)/q_df(1,:,k) - qb_df(4,:)/qb_df(1,:)

            qprime_face(2,:,:,:,k) = q_face(2,:,:,:,k)/q_face(1,:,:,:,k) - qb_face(3,:,:,:)/qb_face(1,:,:,:)
            qprime_face(3,:,:,:,k) = q_face(3,:,:,:,k)/q_face(1,:,:,:,k) - qb_face(4,:,:,:)/qb_face(1,:,:,:)

            qprime(1,:,k) = qprime_temp(1,:,k)
            qprime_df(1,:,k) = dpprime_df(:,k)

            if(flag_pred == 1) then 
                qprime_face(1,:,:,:,k) = qprime_face_temp(1,:,:,:,k)
            else 
                qprime_face(1,:,:,:,k) = dprime_face_corr(:,:,:,k)
            end if
        end do 

    end subroutine momentum_mass


    subroutine rhs_thickness(dp_advec, sum_layer_mass_flux, sum_layer_mass_flux_face, u_edge, v_edge, qprime, qprime_face, uvb_ave, ope_ave, uvb_face_ave, ope_face_ave, q_face)

        ! ===========================================================================================================================
        ! This subroutine is used to predict or correct the layer thickness for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the layer thickness dpprime_df is stored in q_df(1,:,:)
        ! The quadrature points of the layer thickness dpprime is stored in qprime(1,:,:)
        ! The face values of the layer thickness dpprime_face is stored in qprime_face(1,:,:,:,:)
        ! The quadrature points of the layer thickness dp is stored in qb(1,:,:)
        ! The face values of the layer thickness dp_face is stored in qb_face(1,:,:,:,:)
        !
        ! Enforce consistency between the layer masses and the barotropic mass.
        ! ===========================================================================================================================

        use mod_input, only: nlayers, dt, icase, ifilter
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq
        use mod_initial, only: pbprime_df, pbprime, pbprime_face, q_df_mlswe_init, alpha_mlswe
        use mod_create_rhs_mlswe, only: layer_mass_advection_rhs, layer_mass_advection_rhs1
        use mod_layer_terms, only: layer_mass_advection_terms

        use mod_mpi_utilities, only : irank
        use mod_constants, only: gravity


        implicit none
    
        ! Input variables
        real, intent(inout) :: qprime(3,npoin_q,nlayers)
        real, intent(inout) :: qprime_face(3,2,nq,nface,nlayers)
        real, intent(in)    :: uvb_ave(2,npoin_q), ope_ave(npoin_q)
        real, intent(in)    :: uvb_face_ave(2,2,nq,nface), ope_face_ave(2,nq,nface)
        real, intent(in) :: q_face(3,2,nq,nface,nlayers)
    
        ! Output variables
        real, dimension(2,nq, nface, nlayers), intent(out)   :: u_edge, v_edge
        real, intent(out) :: dp_advec(npoin,nlayers), sum_layer_mass_flux_face(2,nq,nface), sum_layer_mass_flux(2,npoin_q)
    
        ! Other variables
        real :: uvdp_temp(2,npoin_q,nlayers), flux_edge(2,nq,nface,nlayers)
        real :: flux_adjust_edge(2,nq,nface,nlayers)
        real :: disp(nq,nface,nlayers)

        ! =========================================== layer mass ===============================================================
    
        ! Compute the mass advection term.   

        call layer_mass_advection_terms(sum_layer_mass_flux,sum_layer_mass_flux_face,u_edge,v_edge,uvdp_temp,flux_edge, &
                qprime,uvb_ave,ope_ave,qprime_face,uvb_face_ave,ope_face_ave, disp)

        ! Compute RHS for layer thickness, return the result in array  dp_advec .

        ! call layer_mass_advection_rhs(dp_advec, uvdp_temp, flux_edge)

        call layer_mass_advection_rhs1(dp_advec, uvdp_temp, flux_edge, q_face, disp)
        
    end subroutine rhs_thickness


    subroutine rhs_momentum(rhs_mom, rhs_visc_bcl, qprime,q_face,qprime_face,ope_ave,one_plus_eta_edge_2_ave,uvb_ave,u_edge,v_edge,Qu_ave,&
        Qv_ave,Quv_ave,H_ave,uvb_face_ave,ope_face_ave,Qu_face_ave,Qv_face_ave,Quv_face_ave,H_face_ave,&
        qprime_df,ope_ave_df,tau_bot_ave,tau_wind_ave,qprime_face2,&
        q2,qprime2, qprime_face3, ope2_ave, uvb_df_ave)

        ! ===========================================================================================================================
        ! This subroutine is used to predict or correct the layer momentum for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the layer momentum is stored in q_df(2:3,:,:)
        ! The quadrature points of the layer momentum are stored in qprime(2:3,:,:)
        ! The face values of the layer momentum are stored in qprime_face(2:3,:,:,:,:)
        ! The quadrature points of the layer momentum are stored in qb(2:5,:,:)
        ! The face values of the layer momentum are stored in qb_face(2:5,:,:,:)
        ! ===========================================================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma_dg_quad, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, dt, ad_mlswe, explt_coriolis, method_visc, ifilter, is_mlswe_linear
        use mod_initial, only: coriolis_df, zbot_face, fdt_bcl, fdt2_bcl, a_bcl, b_bcl,  coriolis_quad, a_bclp, b_bclp
        use mod_create_rhs_mlswe, only: layer_momentum_rhs, interpolate_layer_from_quad_to_node, rhs_layer_shear_stress
        use mod_layer_terms, only: compute_momentum_edge_values, layer_momentum_advec_terms, layer_pressure_terms, layer_windbot_stress_terms

        use mod_layer_terms, only: bcl_wet_dry_mom_df,bcl_wet_dry_mom
        use mod_face, only: imapl_q, imapr_q, normal_vector_q, imapl
        use mod_laplacian_quad, only: create_laplacian_mlswe_layer_v3

        implicit none

        ! Input variables
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: q_face
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face, qprime_face2, qprime_face3
        
        real, dimension(npoin_q), intent(in) :: ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave
        real, dimension(2,2,nq,nface), intent(in) :: uvb_face_ave
        real, dimension(nq,nface), intent(in) :: H_face_ave, one_plus_eta_edge_2_ave
        real, dimension(2, npoin_q), intent(in) :: tau_wind_ave, tau_bot_ave, uvb_ave
        real, dimension(2, nq, nface), intent(in) :: Qu_face_ave, Qv_face_ave, Quv_face_ave, ope_face_ave
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(npoin), intent(in) :: ope_ave_df
        real, dimension(2,nq, nface, nlayers), intent(in)   :: u_edge, v_edge
        real, dimension(3,npoin_q,nlayers), intent(in) :: q2
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime2
        real, dimension(2,npoin), intent(in) :: uvb_df_ave

        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom, rhs_visc_bcl

        ! Local variables

        real, dimension(2,npoin_q,nlayers+1)   :: grad_z
        real, dimension(2,npoin_q,nlayers)     :: coriolis
        real, dimension(2,npoin,nlayers)       :: q_df_temp
        real, dimension(nq,nface,nlayers)      :: udp_left, vdp_left, udp_right, vdp_right, disp
        real, dimension(npoin_q, nlayers)      :: u_udp_temp, v_vdp_temp
        real, dimension(2, npoin_q, nlayers)   :: u_vdp_temp
        real, dimension(2, nq, nface, nlayers) :: udp_flux_edge, vdp_flux_edge
        real, dimension(npoin_q,nlayers)       :: H_r
        real, dimension(2,nq,nface,nlayers)    :: H_r_face
        real, dimension(npoin_q,nlayers+1)     :: p
        real, dimension(npoin,nlayers+1)       :: z_elev
        real, dimension(2,npoin_q,nlayers+1)     :: tau_wind_int, tau_bot_int

        u_udp_temp = 0.0
        v_vdp_temp = 0.0
        u_vdp_temp = 0.0
        udp_flux_edge = 0.0
        vdp_flux_edge = 0.0

        rhs_visc_bcl = 0.0

        ! ==== layer momentum ====

        ! Compute the pressure terms

        call layer_pressure_terms(H_r, H_r_face, p, z_elev, qprime, qprime_face, ope_ave, H_ave, ope_face_ave, &
            zbot_face, H_face_ave, qprime_df, one_plus_eta_edge_2_ave, ope_ave_df, grad_z, ope2_ave)

        call compute_momentum_edge_values(udp_left, vdp_left, udp_right, vdp_right, qprime_face3, uvb_face_ave, ope_face_ave, disp)

        call layer_momentum_advec_terms(u_udp_temp, u_vdp_temp, v_vdp_temp, udp_flux_edge, vdp_flux_edge, &
            q2, qprime, uvb_ave, ope_ave, u_edge, v_edge, Qu_ave, Qv_ave, Quv_ave, &
            udp_left, vdp_left, udp_right, vdp_right, Qu_face_ave, Qv_face_ave, Quv_face_ave)

        ! Compute the wind stress terms

        call layer_windbot_stress_terms(tau_wind_int, tau_bot_int,qprime, tau_bot_ave, tau_wind_ave)

        ! Compute the RHS viscosity terms
        if(method_visc > 0) then 
            call create_laplacian_mlswe_layer_v3(rhs_visc_bcl,qprime2,qprime_face2,uvb_ave,uvb_face_ave)
        end if

        ! Compute the RHS of the layer momentum equation

        call layer_momentum_rhs(rhs_mom, H_r,p,grad_z,u_udp_temp,u_vdp_temp,v_vdp_temp,udp_flux_edge,vdp_flux_edge,H_r_face,tau_wind_int,tau_bot_int,rhs_visc_bcl, q_face, disp)

    end subroutine rhs_momentum


    subroutine apply_consistency(dp_advec, q_df, btp_mass_flux_ave, ope_face_ave, &
        btp_mass_flux_face_ave, sum_layer_mass_flux, sum_layer_mass_flux_face)

        ! ===========================================================================================================================
        ! This subroutine is used to predict or correct the layer thickness for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the layer thickness dpprime_df is stored in q_df(1,:,:)
        ! The quadrature points of the layer thickness dpprime is stored in qprime(1,:,:)
        ! The face values of the layer thickness dpprime_face is stored in qprime_face(1,:,:,:,:)
        ! The quadrature points of the layer thickness dp is stored in qb(1,:,:)
        ! The face values of the layer thickness dp_face is stored in qb_face(1,:,:,:,:)
        !
        ! Enforce consistency between the layer masses and the barotropic mass.
        ! ===========================================================================================================================

        use mod_input, only: nlayers
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq
        use mod_initial, only: pbprime
        use mod_create_rhs_mlswe, only: layer_mass_advection_rhs
        use mod_layer_terms, only: evaluate_dp, evaluate_dp_face, consistency_mass_terms1

        implicit none
    
        ! Input variables
        
        real, intent(in) :: q_df(3,npoin,nlayers)
        real, intent(in)    :: btp_mass_flux_ave(2,npoin_q)
        real, intent(in)    :: sum_layer_mass_flux(2,npoin_q), ope_face_ave(2,nq,nface)
        real, intent(in)    :: btp_mass_flux_face_ave(2,nq,nface), sum_layer_mass_flux_face(2,nq,nface)
    
        ! Output variables
        real, intent(out) :: dp_advec(npoin,nlayers)
    
        ! Other variables
        real :: flux_adjustment(2,npoin_q,nlayers)
        real :: flux_adjust_edge(2,nq,nface,nlayers)
        real :: qprime(3,npoin_q,nlayers), q(3,npoin_q,nlayers), qprime_face(3,2,nq,nface,nlayers), q_face(3,2,nq,nface,nlayers)
        
        call evaluate_dp(q,qprime,q_df, pbprime)
        call evaluate_dp_face(q_face, qprime_face,q, qprime)

        call bcl_create_communicator(qprime_face(1,:,:,:,:),1,nlayers,nq)

        call consistency_mass_terms1(flux_adjustment, flux_adjust_edge, q_df, qprime, &
            sum_layer_mass_flux, btp_mass_flux_ave, ope_face_ave, sum_layer_mass_flux_face, &
            btp_mass_flux_face_ave, qprime_face)

        call layer_mass_advection_rhs(dp_advec, flux_adjustment, flux_adjust_edge)
        
    end subroutine apply_consistency


end module mod_splitting
