
module mod_splitting_v2

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

    use mod_initial, only: coriolis_df, coriolis_quad, tau_wind
        
    implicit none

    public :: thickness_v2, momentum_v2, momentum_mass_v2
    
    contains


    subroutine thickness_v2(q,qprime, q_df, q_face, qprime_face, u_edge, v_edge, uvb_ave, btp_mass_flux_ave, ope_ave, uvb_face_ave, ope_face_ave, &
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

        use mod_input, only: nlayers, dt, icase, ifilter
        use mod_grid, only: npoin, npoin_q, nface
        use mod_basis, only: nq
        use mod_initial, only: pbprime_df, pbprime, pbprime_face, q_df_mlswe_init, alpha_mlswe
        use mod_create_rhs_mlswe, only: layer_mass_advection_rhs, layer_mass_advection_rhs1
        use mod_layer_terms, only: layer_mass_advection_terms, evaluate_dp, evaluate_dp_face, consistency_mass_terms1, &
                                    filter_mlswe

        use mod_mpi_utilities, only : irank
        use mod_constants, only: gravity


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

        call layer_mass_advection_rhs(dp_advec, uvdp_temp, flux_edge)

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
        
    end subroutine thickness_v2


    subroutine momentum_v2(q,qprime,q_df,q_face,qprime_face,qb,ope_ave,one_plus_eta_edge_2_ave,uvb_ave,u_edge,v_edge,Qu_ave,&
        Qv_ave,Quv_ave,H_ave,uvb_face_ave,ope_face_ave,Qu_face_ave,Qv_face_ave,Quv_face_ave,H_face_ave,&
        qprime_df,ope_ave_df,tau_bot_ave,tau_wind_ave,qprime_face2,flag_pred,&
        q2,q_face2,qb_df,qprime2, qprime_face3, ope2_ave, uvb_df_ave, dpprime_df)

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
        use mod_layer_terms, only: compute_momentum_edge_values, layer_momentum_advec_terms, layer_pressure_terms, layer_windbot_stress_terms, &
                                    shear_stress_system, layer_mom_boundary_df, &
                                    filter_mlswe, evaluate_mom, velocity_df, evaluate_mom_face

        use mod_layer_terms, only: bcl_wet_dry_mom_df,bcl_wet_dry_mom
        use mod_face, only: imapl_q, imapr_q, normal_vector_q, imapl
        use mod_laplacian_quad, only: create_laplacian_mlswe_layer_v3, create_laplacian_mlswe_layer_v5, bcl_create_laplacian

        implicit none

        ! Input variables
        real, dimension(3,npoin_q,nlayers), intent(inout) :: q
        real, dimension(3,npoin_q,nlayers), intent(inout) :: qprime
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(3,2,nq,nface,nlayers), intent(inout) :: q_face
        real, dimension(3,2,nq,nface,nlayers), intent(inout) :: qprime_face, qprime_face2, qprime_face3
        real, dimension(4,npoin_q), intent(in) :: qb
        
        real, dimension(npoin_q), intent(in) :: ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave
        real, dimension(2,2,nq,nface), intent(in) :: uvb_face_ave
        real, dimension(nq,nface), intent(in) :: H_face_ave, one_plus_eta_edge_2_ave
        real, dimension(2, npoin_q), intent(in) :: tau_wind_ave, tau_bot_ave, uvb_ave
        real, dimension(2, nq, nface), intent(in) :: Qu_face_ave, Qv_face_ave, Quv_face_ave, ope_face_ave
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(npoin,nlayers), intent(in) :: dpprime_df
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
        
        real, dimension(npoin_q,nlayers) :: dp
        integer :: k, Iq, iquad, I, iface, ilr, n
        real, dimension(npoin) :: tempu1, tempv1

        q_df_temp = 0.0

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
            !call create_laplacian_mlswe_layer_v5(rhs_visc_bcl,qprime2,qprime_face2,uvb_ave,uvb_face_ave)
            !call bcl_create_laplacian(rhs_visc_bcl,qprime2,qprime_face2,uvb_ave,uvb_face_ave, dpprime_df(:,:))
        end if

        ! Compute the RHS of the layer momentum equation

        call layer_momentum_rhs(rhs_mom, H_r,p,grad_z,u_udp_temp,u_vdp_temp,v_vdp_temp,udp_flux_edge,vdp_flux_edge,H_r_face,tau_wind_int,tau_bot_int,rhs_visc_bcl, q_face, disp)

        ! Compute the momentum equation variables for the next time step

        do k = 1,nlayers
            q_df_temp(1,:,k) = q_df(2,:,k) + dt*(rhs_mom(1,:,k) + rhs_visc_bcl(1,:,k))
            q_df_temp(2,:,k) = q_df(3,:,k) + dt*(rhs_mom(2,:,k) + rhs_visc_bcl(2,:,k))
        end do
        
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

    end subroutine momentum_v2


    subroutine momentum_mass_v2(q,qprime,q_df,q_face,qprime_face,qb,qb_face,ope_ave,one_plus_eta_edge_2_ave,uvb_ave,Qu_ave,&
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

        ! Apply filter to the thickness
        do k = 1,nlayers
            q_df(1,:,k) = q_df(1,:,k) + adjust_mass(:)
            
            if(flag_pred == 0 .and. ifilter > 0) then 
                call filter_mlswe(q_df(1,:,k),1)
            end if
        end do
        
        ! Use the adjusted degrees of freedom q_df(1,:,:) to compute revised values of  dp  and  dp' at cell edges and quadrature points.

        call evaluate_dp(q,qprime_temp,q_df, pbprime)
        !call evaluate_dp_face(q_face, qprime_face_temp,q, qprime_temp)

        ! Store the degree of freedom (nodal points) values of dpprime_df ( or q_df(1,:,:)) in the array dpprime_df for use in layer_pressure_terms.

        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)

        do k = 1,nlayers
            dpprime_df(:,k) = q_df(1,:,k) / one_plus_eta_temp(:)
        end do

        ! ==================================== layer momentum ==========================

        call rhs_momentum(rhs_mom, rhs_visc_bcl, qprime,q_face,qprime_face,ope_ave,one_plus_eta_edge_2_ave,uvb_ave,u_edge,v_edge,Qu_ave,&
            Qv_ave,Quv_ave,H_ave,uvb_face_ave,ope_face_ave,Qu_face_ave,Qv_face_ave,Quv_face_ave,H_face_ave,&
            qprime_df,ope_ave_df,tau_bot_ave,tau_wind_ave,qprime_face2,&
            q,qprime2, qprime_face, ope2_ave, uvb_df_ave, qprime_df)

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
        !call evaluate_mom_face(q_face, q)

        ! Compute uprime and vprime at the quad and nodal points

        do k = 1,nlayers
            qprime(1,:,k) = qprime_temp(1,:,k)
            qprime(2,:,k) = q(2,:,k)/q(1,:,k) - qb(3,:)/qb(1,:)
            qprime(3,:,k) = q(3,:,k)/q(1,:,k) - qb(4,:)/qb(1,:)

            qprime_df(1,:,k) = dpprime_df(:,k)
            qprime_df(2,:,k) = q_df(2,:,k)/q_df(1,:,k) - qb_df(3,:)/qb_df(1,:)
            qprime_df(3,:,k) = q_df(3,:,k)/q_df(1,:,k) - qb_df(4,:)/qb_df(1,:)

            !qprime_face(2,:,:,:,k) = q_face(2,:,:,:,k)/q_face(1,:,:,:,k) - qb_face(3,:,:,:)/qb_face(1,:,:,:)
            !qprime_face(3,:,:,:,k) = q_face(3,:,:,:,k)/q_face(1,:,:,:,k) - qb_face(4,:,:,:)/qb_face(1,:,:,:)

            !qprime_face(1,:,:,:,k) = qprime_face_temp(1,:,:,:,k)
        end do 

    end subroutine momentum_mass_v2


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
        q2,qprime2, qprime_face3, ope2_ave, uvb_df_ave, qprime_df2)

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
        use mod_laplacian_quad, only: create_laplacian_mlswe_layer_v3, create_laplacian_mlswe_layer_v5, bcl_create_laplacian

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
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df, qprime_df2
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
            !call create_laplacian_mlswe_layer_v5(rhs_visc_bcl,qprime2,qprime_face2,uvb_ave,uvb_face_ave)

            !call bcl_create_laplacian(rhs_visc_bcl,qprime2,qprime_face2,uvb_ave,uvb_face_ave, qprime_df2(1,:,:))
        end if

        ! Compute the RHS of the layer momentum equation

        call layer_momentum_rhs(rhs_mom, H_r,p,grad_z,u_udp_temp,u_vdp_temp,v_vdp_temp,udp_flux_edge,vdp_flux_edge,H_r_face,tau_wind_int,tau_bot_int,rhs_visc_bcl, q_face, disp)

    end subroutine rhs_momentum


end module mod_splitting_v2
