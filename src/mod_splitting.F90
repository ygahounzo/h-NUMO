
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

    public :: thickness, momentum, momentum_mass
    
    contains

    subroutine thickness(qprime_df, q_df, qb_df, qprime_df_face)

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

        use mod_input, only: nlayers, dt
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq, ngl
        use mod_initial, only: pbprime_df, pbprime
        use mod_layer_terms, only: evaluate_dpp, evaluate_dpp_face, extract_dprime_df_face
        use mod_create_rhs_mlswe, only: layer_mass_rhs

        implicit none
    
        ! Input variables
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(3,2,ngl,nface,nlayers), intent(inout) :: qprime_df_face
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, intent(in)    :: qb_df(4,npoin)
    
        ! Other variables
        real :: one_plus_eta_temp(npoin), dp_advec(npoin,nlayers)
        integer :: k

        ! =========================================== layer mass ===============================================================
    
        ! Compute the mass advection term and RHS for mass equation, return the result in array  dp_advec. 

        call layer_mass_rhs(dp_advec, qprime_df, qprime_df_face)

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
        call apply_consistency(q_df) 

        ! Store the degree of freedom (nodal points) values of dpprime_df ( or q_df(1,:,:)) in the array dpprime_df for use in layer_pressure_terms.

        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)
        do k = 1,nlayers
            qprime_df(1,:,k) = q_df(1,:,k) / one_plus_eta_temp(:)
        end do

        call extract_dprime_df_face(qprime_df_face(1,:,:,:,:), qprime_df(1,:,:))
        
    end subroutine thickness


    subroutine momentum(q_df,qprime_df,qb_df,qprime_df_face)

        ! ===========================================================================================================================
        ! This subroutine is used to correct the layer momentum for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the layer momentum is stored in q_df(2:3,:,:)
        ! The quadrature points of the layer momentum are stored in qprime(2:3,:,:)
        ! The face values of the layer momentum are stored in qprime_face(2:3,:,:,:,:)
        ! ===========================================================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma_dg_quad, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, dt, ad_mlswe
        use mod_initial, only: fdt_bcl, fdt2_bcl, a_bcl, b_bcl
        use mod_create_rhs_mlswe, only: rhs_layer_shear_stress
        use mod_layer_terms, only: layer_mom_boundary_df, evaluate_mom, &
                                    shear_stress_system, velocity_df, evaluate_bcl_v1

        implicit none

        ! Input variables
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(3,2,ngl,nface,nlayers), intent(in) :: qprime_df_face
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
        call rhs_momentum(rhs_mom, qprime_df,q_df, qprime_df_face)

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
        call evaluate_bcl_v1(q_df, qprime_df, qb_df)

    end subroutine momentum


    subroutine momentum_mass(q_df,qprime_df_face,qprime_df,qb_df)

        ! ===========================================================================================================================
        ! This subroutine is used to predict the layer mass and momentum for the splitting system using two-level time integration
        ! The nodal points or degree of freedom of the layer mass is stored in q_df(1,:,:), momentum is stored in q_df(2:3,:,:)
        ! The quadrature points of the layer mass is stored in qprime(1,:,:), momentum are stored in qprime(2:3,:,:)
        ! The face values of the layer mass is stored in qprime_face(1,:,:,:,:), momentum are stored in qprime_face(2:3,:,:,:,:)
        ! ===========================================================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma_dg_quad, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, dt, ad_mlswe
        use mod_initial, only: fdt_bcl, fdt2_bcl, a_bcl, b_bcl, pbprime_df, pbprime
        use mod_create_rhs_mlswe, only: rhs_layer_shear_stress, layer_mass_rhs
        use mod_layer_terms, only: shear_stress_system, layer_mom_boundary_df, filter_mlswe, evaluate_mom, &
                                    velocity_df, evaluate_bcl

        implicit none

        ! Input variables
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(3,2,ngl,nface,nlayers), intent(inout) :: qprime_df_face
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

        call layer_mass_rhs(dp_advec, qprime_df, qprime_df_face)

        do k = 1, nlayers
            q_df(1,:,k) = q_df(1,:,k) + dt*dp_advec(:,k)

            ! Check for negative layer thicknesses.  If any are found, stop the program.
            if(any(q_df(1, :, k) < 0.0)) then
                write(*,*) 'Negative mass in thickness at some points'
                stop
            endif

        end do

        ! consistency through flux adjustment 
        call apply_consistency(q_df)

        ! ==================================== layer momentum ==========================

        call rhs_momentum(rhs_mom, qprime_df, q_df, qprime_df_face)

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
        call evaluate_bcl(qprime_df_face, q_df, qprime_df, qb_df)

    end subroutine momentum_mass

    subroutine rhs_momentum(rhs_mom, qprime_df, q_df, qprime_df_face)

        ! ===========================================================================================================================
        ! This subroutine computes the RHS of the layer momentum equation and viscosity terms and stores them in 
        ! rhs_mom and rhs_visc_bcl
        ! ===========================================================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma_dg_quad, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, method_visc
        use mod_create_rhs_mlswe, only: layer_momentum_rhs
        use mod_layer_terms, only: layer_pressure_terms, layer_momentum_advec_terms_upwind, compute_momentum_edge_values
        use mod_laplacian_quad, only: bcl_create_laplacian

        use mod_variables, only: H_r,u_udp_temp, v_vdp_temp, p, z_elev, udp_left, &
                    vdp_left, udp_right, vdp_right, u_vdp_temp, grad_z, u_edge, v_edge, udp_flux_edge, vdp_flux_edge, &
                    H_r_face, tau_wind_int, tau_bot_int

        implicit none

        ! Input variables
        real, dimension(3,2,ngl,nface,nlayers), intent(in) :: qprime_df_face
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df, q_df
        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom

        real, dimension(2,npoin,nlayers) :: rhs_visc_bcl

        rhs_visc_bcl = 0.0

        if(method_visc > 0) call bcl_create_laplacian(rhs_visc_bcl)

        ! Compute the RHS of the layer momentum equation
        call layer_momentum_rhs(rhs_mom, rhs_visc_bcl, qprime_df, q_df, qprime_df_face)

    end subroutine rhs_momentum

    subroutine apply_consistency(q_df)

        ! ===========================================================================================================================
        ! This subroutine enforce consistency between the layer masses and the barotropic mass through flux adjustment (see Higdon (2015)).
        ! ===========================================================================================================================

        use mod_input, only: nlayers, dt
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq
        use mod_initial, only: pbprime_df
        use mod_create_rhs_mlswe, only: consistency_mass_rhs
        use mod_layer_terms, only: evaluate_dpp, evaluate_dpp_face, evaluate_consistency_face
        use mod_variables, only: btp_mass_flux_face_ave, flux_adjustment, flux_adjust_edge
        use mod_metrics, only: massinv
        
        implicit none
    
        ! Input variables
        real, intent(inout)    :: q_df(3,npoin,nlayers)
    
        ! Local variables
        real :: qprime(3,npoin_q,nlayers)
        integer :: k
        real :: one_plus_eta_temp(npoin), dpprime_df(npoin,nlayers), dp_advec(npoin,nlayers)
        real :: mass_deficit_mass_face(2,2,nq,nface,nlayers)

        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)

        do k = 1,nlayers
            dpprime_df(:,k) = q_df(1,:,k) / one_plus_eta_temp(:)
        end do

        call evaluate_dpp(qprime,dpprime_df)

        ! Consistency flux terms 
        call evaluate_consistency_face(mass_deficit_mass_face,qprime)

        ! RHS of the consistency terms
        call consistency_mass_rhs(dp_advec, qprime, mass_deficit_mass_face)

        ! Apply consistency to the thickness
        do k = 1,nlayers
            q_df(1,:,k) = q_df(1,:,k) + dt*massinv(:)*dp_advec(:,k)
        end do
        
    end subroutine apply_consistency

end module mod_splitting
