
module mod_splitting

    ! ===========================================================================================================================
    ! This module contains the routines for the barotropic-baroclinic splitting
    !   Author: Yao Gahounzo 
    !   Computing PhD 
    !   Boise State University
    !   Date: March 27, 2023
    ! It contains the following routines:
    ! - ti_barotropic: barotropic substem for splitting system using two-level time integration (see Higdon et al. 2005)
    ! - thickness: baroclinic substem for splitting system using two-level time integration 
    ! - momentum: baroclinic substem for splitting system using two-level time integration
    ! These routines are based on Prof. Higdon 1D MLSWE code
    !
    ! ===========================================================================================================================

    use mod_initial, only: coriolis_df, coriolis_quad, tau_wind
        
    implicit none

    public :: ti_barotropic, thickness, momentum, momentum_mass
    
    contains

    subroutine ti_barotropic(qb,qb_face,qb_df, qprime,qprime_face, qprime_df)

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

        use mod_initial, only: N_btp, coriolis_quad, pbprime_df, fdt_btp, fdt2_btp, a_btp, b_btp, tau_wind
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nqx, nqy, nqz, nq
        use mod_input, only: nlayers, dt_btp, ifilter, nlayers
        use mod_rhs_btp, only: create_rhs_btp_momentum, btp_mass_advection_rhs, create_rhs_btp_mom_mass
        use mod_barotropic_terms, only: btp_evaluate_mom, btp_evaluate_mom_face, btp_evaluate_pb, btp_evaluate_pb_face, &
                                        btp_mass_advection_terms, btp_bcl_coeffs, btp_evaluate_mom_dp, btp_evaluate_mom_dp_face, &
                                        btp_bcl_grad_coeffs

        use mod_variables, only: one_plus_eta_edge_2_ave, ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave, &
                                ope_ave_df, uvb_face_ave, btp_mass_flux_face_ave, ope_face_ave, H_face_ave, &
                                Qu_face_ave, Qv_face_ave, Quv_face_ave, one_plus_eta_out, tau_wind_ave, tau_bot_ave, &
                                btp_mass_flux_ave, uvb_ave, one_plus_eta, &
                                Qu_face, Qv_face, one_plus_eta_face, flux_edge, H_bcl_edge, Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, &
                                H_face, one_plus_eta_edge_2, one_plus_eta_edge, &
                                Quu, Qvv, Quv, H, Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, one_plus_eta_df, &
                                btp_mass_flux, tau_bot, uvb_ave_df

        use mod_variables, only: graduvb_face_ave, graduvb_ave

        use mod_input, only: method_visc
        use mod_barotropic_terms, only: btp_mom_boundary_df, compute_btp_terms, compute_btp_mom_terms
        use mod_layer_terms, only: filter_mlswe
        use mod_laplacian_quad, only: btp_create_laplacian_v1, btp_create_laplacian, btp_create_laplacian_v3


        implicit none

        real, dimension(4,npoin_q), intent(inout) :: qb
        real, dimension(4,2,nq,nface), intent(inout) :: qb_face
        real, dimension(4,npoin), intent(inout) :: qb_df
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface, nlayers), intent(in) :: qprime_face
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        real, dimension(4,npoin) :: qb_df_pred
        real, dimension(4,npoin_q) :: qb_init
        real, dimension(4,2,nq,nface) :: qb_face_init
        real, dimension(2,npoin) :: rhs_mom, rhs_visc_btp, rhs_visc_btp1
        real, dimension(npoin) :: pb_advec, tempu, tempv

        real, dimension(2,nq,nface) :: qb_com1
        real, dimension(2,2,nq,nface) :: qb_com2
        real, dimension(3,npoin) :: rhs_mom1
        real, dimension(2,npoin) :: qb_df_temp

        integer :: mstep, I, Iq, iquad, iface, ilr, k
        real :: N_inv

        one_plus_eta_edge_2_ave = 0.0
        uvb_ave  = 0.0
        uvb_ave_df  = 0.0
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
        rhs_visc_btp = 0.0
        ope2_ave = 0.0
        graduvb_face_ave = 0.0
        graduvb_ave = 0.0

        ! Compute baroclinic coefficients in the barotropic momentum fluxes, barotropic pressure forcing, and barotropic
        ! horizontal viscosity terms.  These are needed for the barotropic momentum equation.

        call btp_bcl_coeffs(qprime,qprime_face)
        call btp_bcl_grad_coeffs(qprime_df)

        qb_init = qb
        qb_face_init = qb_face

        ! ========== The time loop begins here ==================

        do mstep = 1, N_btp

            ! ============================================ Prediction steps ==========================================================

            ! Communicate the grid cell face values of the barotropic variables to the neighboring processors.

            call create_communicator_quad(qb_face_init,4)

            ! ******** Barotropic mass & momentum equation ***********

            call compute_btp_terms(qb_init,qb_face_init, qprime, qb_df)

            ! Compute RHS viscosity terms

            if(method_visc == 1) then
                call btp_create_laplacian_v1(rhs_visc_btp,qprime,qprime_face,qb_init,qb_face_init)
            elseif(method_visc == 2) then
                call btp_create_laplacian(rhs_visc_btp,qprime_df,qb_df)
            elseif(method_visc == 3) then 
                call btp_create_laplacian_v3(rhs_visc_btp,qb_df)
            end if

            ! Compute RHS for the barotropic

            call create_rhs_btp_mom_mass(rhs_mom1,qb_init,qb_face_init)

            ! Predict barotropic 

            qb_df_pred(2,:) = qb_df(2,:) + dt_btp * rhs_mom1(1,:)
            qb_df_pred(3,:) = qb_df(3,:) + fdt_btp(:)*qb_df(4,:) + dt_btp*(rhs_mom1(2,:) + rhs_visc_btp(1,:))
            qb_df_pred(4,:) = qb_df(4,:) - fdt_btp(:)*qb_df(3,:) + dt_btp*(rhs_mom1(3,:) + rhs_visc_btp(2,:))

            qb_df_pred(1,:) = qb_df_pred(2,:) + pbprime_df(:)

            call btp_mom_boundary_df(qb_df_pred(3:4,:))

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

            uvb_ave_df(1,:) = uvb_ave_df(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_ave_df(2,:) = uvb_ave_df(2,:) + qb_df(4,:)/qb_df(1,:)

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

            call btp_mass_advection_terms(qb_face)

            ! Compute RHS for the barotropic mass equation

            call btp_mass_advection_rhs(pb_advec)

            ! Correct barotropic mass

            qb_df(2,:) = qb_df(2,:) + 0.5*dt_btp * (pb_advec(:) + rhs_mom1(1,:))
            qb_df(1,:) = qb_df(2,:) + pbprime_df(:)

            ! Evaluate the corrected barotropic mass at the quadrature points and face 

            call btp_evaluate_pb(qb, qb_df)
            call btp_evaluate_pb_face(qb_face, qb)

            ! Communication of the dp values within the processors

            qb_com1(:,:,:) = qb_face(2,:,:,:)

            call create_communicator_quad_1var(qb_face(2,:,:,:))

            ! ******** Barotropic momentum equation ***********

            call compute_btp_mom_terms(qb,qb_face,qprime,qb_df)

            ! Compute RHS viscosity terms

            if(method_visc == 1) then
                call btp_create_laplacian_v1(rhs_visc_btp1,qprime,qprime_face,qb,qb_face)
            elseif(method_visc == 2) then
                call btp_create_laplacian(rhs_visc_btp1,qprime_df,qb_df)
            elseif(method_visc == 3) then 
                call btp_create_laplacian_v3(rhs_visc_btp1,qb_df)
            end if

            ! Compute RHS for the barotropic momentum equation

            call create_rhs_btp_momentum(rhs_mom,qb,qb_face)

            rhs_mom = 0.5*(rhs_mom1(2:3,:) + rhs_mom)
            rhs_visc_btp = 0.5*(rhs_visc_btp1 + rhs_visc_btp)

            ! Correct barotropic momentum

            qb_df_temp(1,:) = qb_df(3,:) + dt_btp*(rhs_mom(1,:) + rhs_visc_btp(1,:))
            qb_df_temp(2,:) = qb_df(4,:) + dt_btp*(rhs_mom(2,:) + rhs_visc_btp(2,:))

            call btp_mom_boundary_df(qb_df_temp(:,:))

            tempu = qb_df_temp(1,:) + fdt2_btp(:)*qb_df(4,:)
            tempv = qb_df_temp(2,:) - fdt2_btp(:)*qb_df(3,:)

            qb_df(3,:) = a_btp(:) * tempu + b_btp(:) * tempv
            qb_df(4,:) = -b_btp(:) * tempu + a_btp(:) * tempv

            call btp_mom_boundary_df(qb_df(3:4,:))

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

            uvb_ave_df(1,:) = uvb_ave_df(1,:) + qb_df(3,:)/qb_df(1,:)
            uvb_ave_df(2,:) = uvb_ave_df(2,:) + qb_df(4,:)/qb_df(1,:)

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

        tau_wind_ave = tau_wind_ave / real(N_btp)
        tau_bot_ave = tau_bot_ave / real(N_btp)
        graduvb_face_ave = N_inv*graduvb_face_ave 
        graduvb_ave = N_inv*graduvb_ave

    end subroutine ti_barotropic

    subroutine thickness(qprime, q_df, qprime_face, dpprime_df, qb_df)

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
        use mod_initial, only: pbprime_df, pbprime
        use mod_layer_terms, only: evaluate_dp, evaluate_dp_face, evaluate_dpp, evaluate_dpp_face
        use mod_variables, only: sum_layer_mass_flux_face, sum_layer_mass_flux
        use mod_create_rhs_mlswe, only: layer_mass_rhs

        implicit none
    
        ! Input variables
        real, dimension(3,npoin_q,nlayers), intent(inout) :: qprime
        real, dimension(3,2,nq,nface,nlayers), intent(inout) :: qprime_face
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, intent(in)    :: qb_df(4,npoin)
    
        ! Output variables
        real, intent(out) :: dpprime_df(npoin,nlayers)
    
        ! Other variables
        real :: one_plus_eta_temp(npoin), dp_advec(npoin,nlayers)
        integer :: k

        ! =========================================== layer mass ===============================================================
    
        ! Compute the mass advection term and RHS for mass equation, return the result in array  dp_advec. 

        call layer_mass_rhs(dp_advec, qprime, qprime_face)

        ! Compute the tentative values of the predicted or corrected 
        ! degrees of freedom for  dp (q_df(1,:,:)).  These would be the values at
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

        ! consistency through flux adjustment 
        call apply_consistency(q_df,sum_layer_mass_flux, sum_layer_mass_flux_face) 

        ! Store the degree of freedom (nodal points) values of dpprime_df ( or q_df(1,:,:)) in the array dpprime_df for use in layer_pressure_terms.

        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)

        do k = 1,nlayers
            dpprime_df(:,k) = q_df(1,:,k) / one_plus_eta_temp(:)
        end do

        ! Use the adjusted degrees of freedom q_df(1,:,:) to compute revised values of dp' at cell edges and quadrature points.
        call evaluate_dpp(qprime,dpprime_df)
        call evaluate_dpp_face(qprime_face,qprime)
        
    end subroutine thickness


    subroutine momentum(qprime,q_df,qprime_face,qprime_df,qb_df, qprime_face3)

        ! ===========================================================================================================================
        ! This subroutine is used to correct the layer momentum for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the layer momentum is stored in q_df(2:3,:,:)
        ! The quadrature points of the layer momentum are stored in qprime(2:3,:,:)
        ! The face values of the layer momentum are stored in qprime_face(2:3,:,:,:,:)
        ! ===========================================================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma_dg_quad, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, dt, ad_mlswe, ifilter
        use mod_initial, only: fdt_bcl, fdt2_bcl, a_bcl, b_bcl
        use mod_create_rhs_mlswe, only: rhs_layer_shear_stress
        use mod_layer_terms, only: layer_mom_boundary_df, evaluate_mom, evaluate_mom_face_v1, &
                                    interpolate_mom, velocity, velocity_face, evaluate_mom_face_all, &
                                    shear_stress_system, velocity_df, evaluate_bcl, evaluate_bcl_v1, &
                                    evaluate_mom_face, evaluate_bcl_v3

        implicit none

        ! Input variables
        
        real, dimension(3,npoin_q,nlayers), intent(inout) :: qprime
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(3,2,nq,nface,nlayers), intent(inout) :: qprime_face, qprime_face3
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(4,npoin), intent(in) :: qb_df

        ! Local variables

        real, dimension(2,npoin_q,nlayers)     :: uv
        real, dimension(2,npoin,nlayers)       :: q_df_temp
        real, dimension(2,npoin,nlayers) :: rhs_mom, rhs_stress
        real, dimension(3,npoin,nlayers) :: q_df3
        real, dimension(3,npoin_q,nlayers) :: q
        
        integer :: k, I
        real, dimension(npoin) :: tempu, tempv

        q_df_temp = 0.0

        ! ==== layer momentum ====

        call rhs_momentum(rhs_mom, qprime,qprime_face,qprime_df,qprime_face3, q_df)

        ! Compute the momentum equation variables for the next time step

        do k = 1,nlayers
            q_df_temp(1,:,k) = q_df(2,:,k) + dt*rhs_mom(1,:,k)
            q_df_temp(2,:,k) = q_df(3,:,k) + dt*rhs_mom(2,:,k)
        end do

        ! Compute the shear stress terms
        if(ad_mlswe > 0.0) then 
      
            do k = 1,nlayers
                do I = 1,npoin
                    tempu(I) = q_df_temp(1,I,k) + fdt2_bcl(I)*q_df(3,I,k)
                    tempv(I) = q_df_temp(2,I,k) - fdt2_bcl(I)*q_df(2,I,k)

                    q_df3(2,I,k) = a_bcl(I)*tempu(I) + b_bcl(I)*tempv(I)
                    q_df3(3,I,k) = - b_bcl(I)*tempu(I) + a_bcl(I)*tempv(I)
                end do
                q_df3(1,:,k) = q_df(1,:,k)
            end do

            !call layer_mom_boundary_df(q_df3(2:3,:,:))

            !Extract velocity from the momentum at nodal pts
            call velocity_df(q_df3, qb_df)

            ! Evaluate velocity and momentum at the quad points
            call evaluate_mom(q,q_df3)

            ! Compute the vertical stress terms

            call shear_stress_system(uv,q)

            call rhs_layer_shear_stress(rhs_stress,uv)

            do k = 1,nlayers
                q_df_temp(1,:,k) = q_df_temp(1,:,k) + dt*rhs_stress(1,:,k)
                q_df_temp(2,:,k) = q_df_temp(2,:,k) + dt*rhs_stress(2,:,k)
            end do

        end if ! ad_mlswe > 0.0
        
        ! Add the Coriolis term

        do k = 1,nlayers

            tempu(:) = q_df_temp(1,:,k) + fdt2_bcl(:)*q_df(3,:,k)
            tempv(:) = q_df_temp(2,:,k) - fdt2_bcl(:)*q_df(2,:,k)

            q_df(2,:,k) = a_bcl(:)*tempu(:) + b_bcl(:)*tempv(:)
            q_df(3,:,k) = - b_bcl(:)*tempu(:) + a_bcl(:)*tempv(:)
        end do

        call layer_mom_boundary_df(q_df(2:3,:,:))

        ! Compute dpprime, uprime and vprime at the quad and nodal points

        call evaluate_bcl_v3(q_df, qprime_df, qb_df)

    end subroutine momentum


    subroutine momentum_mass(qprime,q_df,qprime_face,qprime_df,qb_df,qprime_face3)

        ! ===========================================================================================================================
        ! This subroutine is used to predict the layer mass and momentum for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the layer mass is stored in q_df(1,:,:), momentum is stored in q_df(2:3,:,:)
        ! The quadrature points of the layer mass is stored in qprime(1,:,:), momentum are stored in qprime(2:3,:,:)
        ! The face values of the layer mass is stored in qprime_face(1,:,:,:,:), momentum are stored in qprime_face(2:3,:,:,:,:)
        ! ===========================================================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma_dg_quad, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, dt, ad_mlswe, ifilter
        use mod_initial, only: fdt_bcl, fdt2_bcl, a_bcl, b_bcl, pbprime_df, pbprime
        use mod_create_rhs_mlswe, only: rhs_layer_shear_stress, layer_mass_rhs
        use mod_layer_terms, only: shear_stress_system, layer_mom_boundary_df, filter_mlswe, evaluate_mom, &
                                    evaluate_mom_face_v1, evaluate_dp, evaluate_dp_face, &
                                    velocity_df, evaluate_bcl, evaluate_bcl_v1, evaluate_dpp, evaluate_dpp_face, evaluate_bcl_v2

        use mod_layer_terms, only: interpolate_mom, velocity, velocity_face, evaluate_mom_face_all, evaluate_mom_face

        use mod_variables, only: sum_layer_mass_flux_face, sum_layer_mass_flux

        implicit none

        ! Input variables
        real, dimension(3,npoin_q,nlayers), intent(inout) :: qprime
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(3,2,nq,nface,nlayers), intent(inout) :: qprime_face
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face3
        
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(4,npoin), intent(in) :: qb_df

        ! Local variables

        real, dimension(npoin,nlayers) :: dp_advec, dpprime_df
        real, dimension(2,npoin,nlayers) :: q_df_temp
        real, dimension(2,npoin,nlayers) :: rhs_mom, rhs_stress
        real, dimension(npoin) :: one_plus_eta_temp

        real, dimension(3,npoin_q,nlayers) :: qprime_temp, q
        real, dimension(3,2,nq,nface,nlayers) :: qprime_face_temp 
        real, dimension(3,npoin,nlayers) :: q_df3
        real, dimension(2,npoin_q,nlayers) :: uv
        
        integer :: k, I
        real, dimension(npoin) :: tempu, tempv

        q_df_temp = 0.0

        ! =========================================== layer mass ===============================================================

        call layer_mass_rhs(dp_advec, qprime, qprime_face)

        do k = 1, nlayers

            q_df(1,:,k) = q_df(1,:,k) + dt*dp_advec(:,k)

            ! Check for negative layer thicknesses.  If any are found, stop the program.
            if(any(q_df(1, :, k) < 0.0)) then
                write(*,*) 'Negative mass in thickness at some points'
                stop
            endif

        end do

        ! consistency through flux adjustment 
        call apply_consistency(q_df,sum_layer_mass_flux, sum_layer_mass_flux_face)

        ! ==================================== layer momentum ==========================

        call rhs_momentum(rhs_mom, qprime, qprime_face, qprime_df,qprime_face, q_df)

        ! Compute the momentum equation variables for the next time step

        do k = 1,nlayers
            q_df_temp(1,:,k) = q_df(2,:,k) + dt*rhs_mom(1,:,k)
            q_df_temp(2,:,k) = q_df(3,:,k) + dt*rhs_mom(2,:,k)
        end do

        ! Compute the shear stress terms
        if(ad_mlswe > 0.0) then 
            
            do k = 1,nlayers
                do I = 1,npoin
                    tempu(I) = q_df_temp(1,I,k) + fdt2_bcl(I)*q_df(3,I,k)
                    tempv(I) = q_df_temp(2,I,k) - fdt2_bcl(I)*q_df(2,I,k)

                    q_df3(2,I,k) = a_bcl(I)*tempu(I) + b_bcl(I)*tempv(I)
                    q_df3(3,I,k) = - b_bcl(I)*tempu(I) + a_bcl(I)*tempv(I)
                end do
                q_df3(1,:,k) = q_df(1,:,k)
            end do

            !Extract velocity from the momentum at nodal pts
            call velocity_df(q_df3, qb_df)

            ! Evaluate velocity and momentum at the quad points
            call evaluate_mom(q,q_df3)

            ! Compute the vertical stress terms

            call shear_stress_system(uv,q)

            call rhs_layer_shear_stress(rhs_stress,uv)

            do k = 1,nlayers
                q_df_temp(1,:,k) = q_df_temp(1,:,k) + dt*rhs_stress(1,:,k)
                q_df_temp(2,:,k) = q_df_temp(2,:,k) + dt*rhs_stress(2,:,k)
            end do

        end if ! ad_mlswe > 0
        
        ! Add the Coriolis term

        do k = 1,nlayers

            tempu(:) = q_df_temp(1,:,k) + fdt2_bcl(:)*q_df(3,:,k)
            tempv(:) = q_df_temp(2,:,k) - fdt2_bcl(:)*q_df(2,:,k)

            q_df(2,:,k) = a_bcl(:)*tempu(:) + b_bcl(:)*tempv(:)
            q_df(3,:,k) = - b_bcl(:)*tempu(:) + a_bcl(:)*tempv(:)
        end do

        call layer_mom_boundary_df(q_df(2:3,:,:))

        ! Compute dpprime, uprime and vprime at the quad and nodal points

        call evaluate_bcl_v2(qprime, qprime_face, q_df, qprime_df, qb_df)

    end subroutine momentum_mass


    subroutine rhs_thickness(dp_advec, qprime, qprime_face)

        ! ===========================================================================================================================
        ! This subroutine compute the RHS of the layer mass equation and stored it in dp_advec

        ! ===========================================================================================================================

        use mod_input, only: nlayers, dt, icase, ifilter
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq
        use mod_create_rhs_mlswe, only: layer_mass_advection_rhs
        use mod_layer_terms, only: layer_mass_advection_terms
        use mod_variables, only: u_edge, v_edge, sum_layer_mass_flux_face, sum_layer_mass_flux, uvdp_temp, flux_edge_bcl

        implicit none
    
        ! Input variables
        real, intent(in) :: qprime(3,npoin_q,nlayers)
        real, intent(in) :: qprime_face(3,2,nq,nface,nlayers)
    
        ! Output variables
        real, intent(out) :: dp_advec(npoin,nlayers)

        ! =========================================== layer mass ===============================================================
    
        ! Compute the mass advection term.   

        call layer_mass_advection_terms(qprime,qprime_face)

        ! Compute RHS for layer thickness, return the result in array  dp_advec .

        call layer_mass_advection_rhs(dp_advec, uvdp_temp, flux_edge_bcl)
        
    end subroutine rhs_thickness


    subroutine rhs_momentum(rhs_mom, qprime,qprime_face,qprime_df, qprime_face3, q_df)

        ! ===========================================================================================================================
        ! This subroutine computes the RHS of the layer momentum equation and viscosity terms and stores them in 
        ! rhs_mom and rhs_visc_bcl
        ! ===========================================================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma_dg_quad, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, method_visc
        use mod_create_rhs_mlswe, only: layer_momentum_rhs, rhs_layer_shear_stress
        use mod_layer_terms, only: compute_momentum_edge_values, layer_momentum_advec_terms, layer_pressure_terms, &
                                    layer_windbot_stress_terms, layer_momentum_advec_terms_upwind, layer_momentum_advec_terms_upwind_v1
        use mod_laplacian_quad, only: bcl_create_laplacian, bcl_create_laplacian_v1, bcl_create_laplacian_v2

        use mod_variables, only: H_r,u_udp_temp, v_vdp_temp, p, z_elev, udp_left, &
                    vdp_left, udp_right, vdp_right, u_vdp_temp, grad_z, u_edge, v_edge, udp_flux_edge, vdp_flux_edge, &
                    H_r_face, tau_wind_int, tau_bot_int

        implicit none

        ! Input variables
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face, qprime_face3
        
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df, q_df

        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom

        real, dimension(2,npoin,nlayers) :: rhs_visc_bcl

        rhs_visc_bcl = 0.0

        ! Compute the pressure terms

        call layer_pressure_terms(qprime, qprime_face, qprime_df)

        call compute_momentum_edge_values(qprime_face3)

        call layer_momentum_advec_terms_upwind_v1(q_df, qprime)

        ! Compute the wind stress terms

        call layer_windbot_stress_terms(qprime)

        if(method_visc == 3) then
            call bcl_create_laplacian_v2(rhs_visc_bcl)
        end if

        ! Compute the RHS of the layer momentum equation

        call layer_momentum_rhs(rhs_mom, rhs_visc_bcl)

    end subroutine rhs_momentum


    subroutine apply_consistency(q_df,sum_layer_mass_flux, sum_layer_mass_flux_face)

        ! ===========================================================================================================================
        ! This subroutine enforce consistency between the layer masses and the barotropic mass through flux adjustment (see Higdon (2015)).
        ! ===========================================================================================================================

        use mod_input, only: nlayers, dt
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq
        use mod_initial, only: pbprime_df
        use mod_create_rhs_mlswe, only: layer_mass_advection_rhs, consistency_mass_rhs
        use mod_layer_terms, only: consistency_mass_terms, evaluate_dpp, evaluate_dpp_face
        use mod_variables, only: btp_mass_flux_face_ave, flux_adjustment, flux_adjust_edge

        implicit none
    
        ! Input variables
        
        real, intent(inout)    :: q_df(3,npoin,nlayers)
        real, dimension(2,npoin_q), intent(in)  :: sum_layer_mass_flux
        real, dimension(2,nq,nface), intent(in) :: sum_layer_mass_flux_face
    
        ! Other variables

        real :: qprime(3,npoin_q,nlayers), qprime_face(3,2,nq,nface,nlayers)
        real, dimension(2,2,nq,nface)  :: flux_deficit_mass_face
        integer :: k
        real :: one_plus_eta_temp(npoin), dpprime_df(npoin,nlayers), dp_advec(npoin,nlayers)

        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)

        do k = 1,nlayers
            dpprime_df(:,k) = q_df(1,:,k) / one_plus_eta_temp(:)
        end do

        call evaluate_dpp(qprime,dpprime_df)
        call evaluate_dpp_face(qprime_face,qprime)

        ! Communicate the interface values within the neighboring processors
        call bcl_create_communicator(qprime_face(1,:,:,:,:),1,nlayers,nq)

        flux_deficit_mass_face(1,1,:,:) = btp_mass_flux_face_ave(1,:,:) - sum_layer_mass_flux_face(1,:,:)
        flux_deficit_mass_face(2,1,:,:) = btp_mass_flux_face_ave(2,:,:) - sum_layer_mass_flux_face(2,:,:)

        flux_deficit_mass_face(1,2,:,:) = flux_deficit_mass_face(1,1,:,:)
        flux_deficit_mass_face(2,2,:,:) = flux_deficit_mass_face(2,1,:,:)

        call create_communicator_quad(flux_deficit_mass_face,2)

        ! Consistency flux terms 
        !call consistency_mass_terms(qprime, qprime_face, flux_deficit_mass_face)

        ! RHS of the consistency terms
        !call layer_mass_advection_rhs(dp_advec, flux_adjustment, flux_adjust_edge)

        call consistency_mass_rhs(dp_advec, qprime, qprime_face, flux_deficit_mass_face)

        ! Apply consistency to the thickness

        do k = 1,nlayers
            q_df(1,:,k) = q_df(1,:,k) + dt*dp_advec(:,k)
        end do
        
    end subroutine apply_consistency

end module mod_splitting