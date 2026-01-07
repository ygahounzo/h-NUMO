
module mod_splitting

    ! ============================================================================================
    ! This module contains the routines for the barotropic-baroclinic splitting
    !   Author: Yao Gahounzo 
    !   Computing PhD 
    !   Boise State University
    !   Date: March 27, 2023
    ! It contains the following routines:
    ! - thickness: baroclinic substem for splitting system using two-level time integration 
    ! - momentum: baroclinic substem for splitting system using two-level time integration
    ! These routines are based on Prof. Higdon 1D MLSWE code
    !
    ! ============================================================================================

    use mod_initial, only: coriolis_df, coriolis_quad, tau_wind
        
    implicit none

    public :: thickness, momentum, momentum_mass
    
    contains

    subroutine thickness(qprime_df, q_df, qb_df)

        ! ========================================================================================
        ! This subroutine is used to predict or correct the layer thickness for the splitting 
        ! system using two-level time integration
        ! The nodal points or degree of freedom of the layer thickness is stored in qprime_df(1,:,:)
        ! The quadrature points of the layer thickness dp is stored in q_df(1,:,:)
        !
        ! Enforce consistency between the layer masses and the barotropic mass.
        ! ========================================================================================

        use mod_input, only: nlayers, dt
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq, ngl
        use mod_initial, only: pbprime_df, alpha_mlswe
        use mod_layer_terms, only: extract_dprime_df_face
        use mod_create_rhs_mlswe, only: layer_mass_rhs

        implicit none
    
        ! Input variables
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, intent(in)    :: qb_df(4,npoin)
    
        ! Other variables
        real :: one_plus_eta_temp(npoin), dp_advec(npoin,nlayers)
        integer :: k

        ! =========================================== layer mass =================================
    
        ! Compute the mass advection term and RHS for mass equation, 
        ! return the result in array  dp_advec. 

        call layer_mass_rhs(dp_advec, qprime_df)

        ! Compute the tentative values of the predicted or corrected 
        ! degrees of freedom for  dp (q_df(1,:,:)).  These would be the values at
        ! baroclinic time level  n+1,  except for the need to enforce
        ! consistency between the layer masses and the barotropic mass 
        ! and for the possible need to use a variation limiter to 
        ! prevent pointwise negative layer thicknesses. 
        
        do k = 1, nlayers

            q_df(1,:,k) = q_df(1,:,k) + dt*dp_advec(:,k)

            ! Check for negative layer thicknesses.  If any are found, stop the program.
            ! if(any(q_df(1, :, k) < 0.0)) then
            !     write(*,*) 'Negative mass in thickness at some points'
            !     stop
            ! endif
        end do

        ! Store the degree of freedom (nodal points) values of dpprime_df
        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)
        do k = 1,nlayers
            qprime_df(1,:,k) = q_df(1,:,k) / one_plus_eta_temp(:)
        end do
        
    end subroutine thickness


    subroutine momentum(q_df,qprime_df,qb_df)

        ! ========================================================================================
        ! This subroutine is used to correct the layer momentum for the splitting system 
        ! The nodal points or degree of freedom of the layer momentum is stored in q_df(2:3,:,:)
        ! ========================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma_dg_quad, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, dt, ad_mlswe
        use mod_initial, only: fdt_bcl, fdt2_bcl, a_bcl, b_bcl, alpha_mlswe
        use mod_create_rhs_mlswe, only: rhs_layer_shear_stress
        use mod_layer_terms, only: layer_mom_boundary_df, &
                                    velocity_df, extract_velocity
        use mod_metrics, only: massinv

        implicit none

        ! Input variables
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(4,npoin), intent(in) :: qb_df

        ! Local variables
        real, dimension(2,npoin_q,nlayers)     :: uv
        real, dimension(2,npoin,nlayers)       :: q_df_temp, uv_df
        real, dimension(2,npoin,nlayers) :: rhs_mom, rhs_stress
        real, dimension(3,npoin,nlayers) :: q_df3
        real, dimension(3,npoin_q,nlayers) :: q
        
        integer :: k, I
        real, dimension(npoin) :: tempu, tempv

        q_df_temp = 0.0

        ! ==== layer momentum ====
        call rhs_momentum(rhs_mom, qprime_df,q_df)

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

            !call layer_mom_boundary_df(q_df3)

            ! Velocity smoothing from the momentum
            call velocity_df(q_df3, qb_df)

            ! Compute the vertical stress terms
            call rhs_layer_shear_stress(rhs_stress,uv)

            do k = 1,nlayers
                q_df_temp(1,:,k) = q_df_temp(1,:,k) + dt*(massinv(:)*rhs_stress(1,:,k))
                q_df_temp(2,:,k) = q_df_temp(2,:,k) + dt*(massinv(:)*rhs_stress(2,:,k))
            end do
        end if ! ad_mlswe > 0.0
        
        ! Add the Coriolis term
        do k = 1,nlayers

            tempu(:) = q_df_temp(1,:,k) + fdt2_bcl(:)*q_df(3,:,k)
            tempv(:) = q_df_temp(2,:,k) - fdt2_bcl(:)*q_df(2,:,k)
            q_df(2,:,k) = a_bcl(:)*tempu(:) + b_bcl(:)*tempv(:)
            q_df(3,:,k) = - b_bcl(:)*tempu(:) + a_bcl(:)*tempv(:)
        end do

        call layer_mom_boundary_df(q_df)

        ! Compute dpprime, uprime and vprime at the nodal points
        call extract_velocity(uv_df, q_df, qb_df)

        do k = 1,nlayers
            q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
            q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
        end do

    end subroutine momentum

    subroutine momentum_mass(q_df, qprime_df, qb_df)

        ! ========================================================================================
        ! This subroutine is used to predict the layer mass and momentum for the splitting system
        ! The nodal points or degree of freedom of the layer mass is stored in q_df(1,:,:),
        ! momentum is stored in q_df(2:3,:,:)
        ! ========================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma_dg_quad, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, dt, ad_mlswe
        use mod_initial, only: fdt_bcl, fdt2_bcl, a_bcl, b_bcl, pbprime_df, alpha_mlswe
        use mod_create_rhs_mlswe, only: rhs_layer_shear_stress
        use mod_layer_terms, only: layer_mom_boundary_df, &
                                    velocity_df, extract_qprime_df_face,extract_velocity
        use mod_metrics, only: massinv

        implicit none

        ! Input variables
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(4,npoin), intent(in) :: qb_df

        ! Local variables
        real, dimension(npoin,nlayers) :: dp_advec, dpprime_df
        real, dimension(2,npoin,nlayers) :: uv_df
        real, dimension(2,npoin,nlayers) :: rhs_stress
        real, dimension(npoin) :: one_plus_eta_temp
        real, dimension(3,npoin_q,nlayers) :: qprime_temp, q
        real, dimension(3,npoin,nlayers) :: q_df3
        real, dimension(2,npoin_q,nlayers) :: uv
        integer :: k, I
        real, dimension(npoin) :: tempu, tempv
        real, dimension(3,npoin,nlayers) :: rhs, q_df_temp

        call create_rhs_bcl(rhs, qprime_df, q_df)

        ! Compute the momentum equation variables for the next time step
        do k = 1,nlayers
            q_df_temp(1,:,k) = q_df(1,:,k) + dt*rhs(1,:,k)
            q_df_temp(2,:,k) = q_df(2,:,k) + dt*rhs(2,:,k)
            q_df_temp(3,:,k) = q_df(3,:,k) + dt*rhs(3,:,k)

            ! Check for negative layer thicknesses.  If any are found, stop the program.
            if(any(q_df(1, :, k) < 0.0)) then
                write(*,*) 'Negative mass in thickness at some points'
                stop
            endif
        end do

        ! Compute the shear stress terms
        if(ad_mlswe > 0.0) then 
            do k = 1,nlayers
                do I = 1,npoin
                    tempu(I) = q_df_temp(2,I,k) + fdt2_bcl(I)*q_df(3,I,k)
                    tempv(I) = q_df_temp(3,I,k) - fdt2_bcl(I)*q_df(2,I,k)
                    q_df3(2,I,k) = a_bcl(I)*tempu(I) + b_bcl(I)*tempv(I)
                    q_df3(3,I,k) = - b_bcl(I)*tempu(I) + a_bcl(I)*tempv(I)
                end do
                q_df3(1,:,k) = q_df(1,:,k)
            end do

            !call layer_mom_boundary_df(q_df3)

            ! Velocity smoothing from the momentum
            call velocity_df(q_df3, qb_df)

            ! Compute the vertical stress terms
            call rhs_layer_shear_stress(rhs_stress,q_df3)

            do k = 1,nlayers
                q_df_temp(2,:,k) = q_df_temp(2,:,k) + dt*(massinv(:)*rhs_stress(1,:,k))
                q_df_temp(3,:,k) = q_df_temp(3,:,k) + dt*(massinv(:)*rhs_stress(2,:,k))
            end do
        end if ! ad_mlswe > 0
        
        ! Add the Coriolis term
        do k = 1,nlayers

            tempu(:) = q_df_temp(2,:,k) + fdt2_bcl(:)*q_df(3,:,k)
            tempv(:) = q_df_temp(3,:,k) - fdt2_bcl(:)*q_df(2,:,k)
            q_df(2,:,k) = a_bcl(:)*tempu(:) + b_bcl(:)*tempv(:)
            q_df(3,:,k) = - b_bcl(:)*tempu(:) + a_bcl(:)*tempv(:)

            q_df(1,:,k) = q_df_temp(1,:,k)
        end do

        call layer_mom_boundary_df(q_df)

        ! Compute dpprime, uprime and vprime at the quad and nodal points
        call extract_velocity(uv_df, q_df, qb_df)

        do k = 1,nlayers
            q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
            q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
        end do

        call extract_qprime_df_face(qprime_df, q_df, qb_df)

    end subroutine momentum_mass

    subroutine rhs_momentum(rhs_mom, qprime_df, q_df)

        ! =======================================================================================
        ! This subroutine computes the RHS of the layer momentum equation and viscosity terms
        ! and stores them in rhs_mom and rhs_visc_bcl
        ! =======================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, method_visc
        use mod_create_rhs_mlswe, only: layer_momentum_rhs
        use mod_laplacian_quad, only: bcl_create_laplacian, bcl_create_laplacian_v2

        implicit none

        ! Input variables
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df, q_df
        real, dimension(2,npoin,nlayers), intent(out) :: rhs_mom

        real, dimension(2,npoin,nlayers) :: rhs_visc_bcl

        rhs_visc_bcl = 0.0

        if (method_visc > 1) then
            call bcl_create_laplacian(rhs_visc_bcl)
        endif 

        ! Compute the RHS of the layer momentum equation
        call layer_momentum_rhs(rhs_mom, rhs_visc_bcl, qprime_df, q_df)

    end subroutine rhs_momentum

    subroutine create_rhs_bcl(rhs, qprime_df, q_df)

        ! =======================================================================================
        ! This subroutine computes the RHS of the layer momentum equation and viscosity terms
        ! and stores them in rhs_mom and rhs_visc_bcl
        ! =======================================================================================

        use mod_grid, only: npoin, npoin_q, nface, face, intma
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers, method_visc
        use mod_create_rhs_mlswe, only: bcl_rhs
        use mod_laplacian_quad, only: bcl_create_laplacian, bcl_create_laplacian_v2

        implicit none

        ! Input variables
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df, q_df
        real, dimension(3,npoin,nlayers), intent(out) :: rhs

        real, dimension(2,npoin,nlayers) :: rhs_visc_bcl

        rhs_visc_bcl = 0.0

        if (method_visc > 1) then
            call bcl_create_laplacian(rhs_visc_bcl)
        endif 

        ! Compute the RHS of the layer momentum equation
        call bcl_rhs(rhs, rhs_visc_bcl, qprime_df, q_df)

    end subroutine create_rhs_bcl

end module mod_splitting
