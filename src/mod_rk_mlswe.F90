module mod_rk_mlswe

    ! ============================================================================================
    ! This module contains the routine for the barotropic solver
    !   Author: Yao Gahounzo
    !   Date: March 27, 2023
    ! It contains the following routine:
    ! - ti_barotropic_ssprk_mlswe: based on SSPRK time integration
    !   (see Higdon et al. 2005)
    !
    ! =============================================================================================

    implicit none

    public :: ti_barotropic_ssprk_mlswe

    contains

    subroutine ti_barotropic_ssprk_mlswe_v0(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, &
                                             qb_df, qprime_df)

        use mod_grid,              only: grid
        use mod_input,             only: input
        use mod_initial,           only: initial
        use mod_basis,             only: basis
        use mod_face,              only: face_CS
        use mod_variables,         only: btp_CS
        use mod_parallel,          only: parallel_CS
        use mod_ref,               only: mref
        use mod_mpi_communicator,  only: mpi_communicator
        use mod_metrics,           only: metrics
        use mod_tensor,            only: tensor_CS
        use mod_rhs_btp,           only: create_rhs_btp
        use mod_barotropic_terms,  only: btp_mom_boundary_df

        implicit none

        type(grid),             intent(in)    :: G
        type(input),            intent(in)    :: inp
        type(basis),            intent(in)    :: b
        type(face_CS),          intent(in)    :: mf
        type(parallel_CS),      intent(in)    :: par
        type(initial),          intent(in)    :: init
        type(mref),             intent(inout) :: ref
        type(mpi_communicator), intent(inout)    :: mpic
        type(metrics),          intent(in)    :: mt
        type(tensor_CS),        intent(in)    :: tsp
        type(btp_CS),           intent(inout) :: btp

        real, dimension(4,G%npoin),             intent(inout) :: qb_df
        real, dimension(3,G%npoin,inp%nlayers), intent(in)    :: qprime_df

        real, dimension(3,G%npoin) :: rhs
        real, dimension(4,G%npoin) :: qb0_df, qb1_df, qb2_df
        integer :: mstep, I, Iq, iquad, iface, ilr, k, ik
        real :: N_inv, a0, a1, a2, beta, dtt

        btp%one_plus_eta_edge_2_ave = 0.0
        btp%uvb_ave                 = 0.0
        btp%uvb_ave_df              = 0.0
        btp%ope_ave                 = 0.0
        btp%btp_mass_flux_ave       = 0.0

        btp%H_ave    = 0.0
        btp%Qu_ave   = 0.0
        btp%Qv_ave   = 0.0
        btp%Quv_ave  = 0.0

        btp%ope2_ave_df            = 0.0
        btp%uvb_face_ave           = 0.0
        btp%ope_face_ave           = 0.0
        btp%ope2_face_ave          = 0.0
        btp%btp_mass_flux_face_ave = 0.0

        btp%H_face_ave    = 0.0
        btp%Qu_face_ave   = 0.0
        btp%Qv_face_ave   = 0.0
        btp%Quv_face_ave  = 0.0

        btp%tau_wind_ave     = 0.0
        btp%tau_bot_ave      = 0.0
        btp%ope2_ave         = 0.0
        btp%graduvb_face_ave = 0.0
        btp%graduvb_ave      = 0.0

        qb2_df = 0.0

        ! ========== The time loop begins here ==================

        do mstep = 1, init%N_btp

            qb0_df = qb_df
            qb1_df = qb_df

            do ik=1,inp%kstages

                dtt = inp%dt_btp*init%ssprk_beta(ik)
                btp%ope2_ave_df     = btp%ope2_ave_df     + (1.0 + qb1_df(2,:)/init%pbprime_df(:))**2
                btp%uvb_ave_df(1,:) = btp%uvb_ave_df(1,:) + qb1_df(3,:)/qb1_df(1,:)
                btp%uvb_ave_df(2,:) = btp%uvb_ave_df(2,:) + qb1_df(4,:)/qb1_df(1,:)

                ! Compute RHS for the barotropic
                call create_rhs_btp(G, inp, b, mf, par, btp, init, ref, mpic, mt, tsp, &
                                    rhs, qb1_df, qprime_df)

                ! Update barotropic variables
                qb_df(2,:) = init%ssprk_a(ik,1)*qb0_df(2,:) + init%ssprk_a(ik,2)*qb1_df(2,:) + &
                             init%ssprk_a(ik,3)*qb2_df(2,:) + dtt*rhs(1,:)
                qb_df(3,:) = init%ssprk_a(ik,1)*qb0_df(3,:) + init%ssprk_a(ik,2)*qb1_df(3,:) + &
                             init%ssprk_a(ik,3)*qb2_df(3,:) + dtt*rhs(2,:)
                qb_df(4,:) = init%ssprk_a(ik,1)*qb0_df(4,:) + init%ssprk_a(ik,2)*qb1_df(4,:) + &
                             init%ssprk_a(ik,3)*qb2_df(4,:) + dtt*rhs(3,:)

                qb_df(1,:) = qb_df(2,:) + init%pbprime_df(:)

                call btp_mom_boundary_df(G, b, mf, qb_df)

                qb1_df = qb_df

                if (inp%kstages == 5 .and. ik == 2) qb2_df = qb_df

            end do
            btp%tau_wind_ave = btp%tau_wind_ave + init%tau_wind

        end do

        ! Compute time averages over all barotropic substeps of the baroclinic time interval.

        N_inv = 1.0 / real(inp%kstages*init%N_btp)

        btp%uvb_ave_df           = N_inv * btp%uvb_ave_df
        btp%graduvb_face_ave     = N_inv * btp%graduvb_face_ave
        btp%graduvb_ave          = N_inv * btp%graduvb_ave
        btp%tau_wind_ave         = btp%tau_wind_ave / real(init%N_btp)
        btp%ope2_ave_df          = N_inv * btp%ope2_ave_df
        btp%ope2_ave             = N_inv * btp%ope2_ave
        btp%ope_ave              = N_inv * btp%ope_ave
        btp%H_ave                = N_inv * btp%H_ave
        btp%Qu_ave               = N_inv * btp%Qu_ave
        btp%Qv_ave               = N_inv * btp%Qv_ave
        btp%Quv_ave              = N_inv * btp%Quv_ave
        btp%btp_mass_flux_ave    = N_inv * btp%btp_mass_flux_ave
        btp%tau_bot_ave          = N_inv * btp%tau_bot_ave
        btp%ope_face_ave             = N_inv * btp%ope_face_ave
        btp%ope2_face_ave            = N_inv * btp%ope2_face_ave
        btp%H_face_ave               = N_inv * btp%H_face_ave
        btp%Qu_face_ave              = N_inv * btp%Qu_face_ave
        btp%Qv_face_ave              = N_inv * btp%Qv_face_ave
        btp%btp_mass_flux_face_ave   = N_inv * btp%btp_mass_flux_face_ave
        btp%one_plus_eta_edge_2_ave  = N_inv * btp%one_plus_eta_edge_2_ave
        btp%uvb_ave                  = N_inv * btp%uvb_ave
        btp%uvb_face_ave             = N_inv * btp%uvb_face_ave

    end subroutine ti_barotropic_ssprk_mlswe_v0

    subroutine ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, &
                                          qb_df, qprime_df)

        use mod_grid,              only: grid
        use mod_input,             only: input
        use mod_initial,           only: initial
        use mod_basis,             only: basis
        use mod_face,              only: face_CS
        use mod_variables,         only: btp_CS
        use mod_parallel,          only: parallel_CS
        use mod_ref,               only: mref
        use mod_mpi_communicator,  only: mpi_communicator
        use mod_metrics,           only: metrics
        use mod_tensor,            only: tensor_CS
        use mod_rhs_btp,           only: create_rhs_btp
        use mod_barotropic_terms,  only: btp_mom_boundary_df

        implicit none

        type(grid),             intent(in)    :: G
        type(input),            intent(in)    :: inp
        type(basis),            intent(in)    :: b
        type(face_CS),          intent(in)    :: mf
        type(parallel_CS),      intent(in)    :: par
        type(initial),          intent(in)    :: init
        type(mref),             intent(inout) :: ref
        type(mpi_communicator), intent(inout)    :: mpic
        type(metrics),          intent(in)    :: mt
        type(tensor_CS),        intent(in)    :: tsp
        type(btp_CS),           intent(inout) :: btp

        real, dimension(4,G%npoin),             intent(inout) :: qb_df
        real, dimension(3,G%npoin,inp%nlayers), intent(in)    :: qprime_df

        real, dimension(4,G%npoin) :: qb0_df, qb1_df, qb2_df
        real, dimension(G%npoin)   :: inv_pb1
        integer :: mstep, ik
        real    :: N_inv, a0, a1, a2, dtt

        ! Zero-initialise all accumulation buffers
        btp%one_plus_eta_edge_2_ave = 0.0;  btp%uvb_ave             = 0.0
        btp%uvb_ave_df              = 0.0;  btp%ope_ave             = 0.0
        btp%btp_mass_flux_ave       = 0.0;  btp%H_ave               = 0.0
        btp%Qu_ave                  = 0.0;  btp%Qv_ave              = 0.0
        btp%Quv_ave                 = 0.0;  btp%ope2_ave_df         = 0.0
        btp%uvb_face_ave            = 0.0;  btp%ope_face_ave        = 0.0
        btp%ope2_face_ave           = 0.0;  btp%btp_mass_flux_face_ave = 0.0
        btp%H_face_ave              = 0.0;  btp%Qu_face_ave         = 0.0
        btp%Qv_face_ave             = 0.0;  btp%Quv_face_ave        = 0.0
        btp%tau_wind_ave            = 0.0;  btp%tau_bot_ave         = 0.0
        btp%ope2_ave                = 0.0;  btp%graduvb_face_ave    = 0.0
        btp%graduvb_ave             = 0.0

        ! Push zeroed accumulators to the device so GPU accumulation starts from 0.
        !$acc update device(btp%H_ave, btp%Qu_ave, btp%Qv_ave, btp%Quv_ave, btp%tau_bot_ave, &
        !$acc               btp%ope_ave, btp%ope2_ave, btp%btp_mass_flux_ave, btp%uvb_ave)

        qb2_df = 0.0

        !$acc update device(qb_df)

        ! Time loop for the barotropic solver, with SSPRK time integration.
        do mstep = 1, init%N_btp

            qb0_df = qb_df

            do ik = 1, inp%kstages

                a0  = init%ssprk_a(ik,1)
                a1  = init%ssprk_a(ik,2)
                a2  = init%ssprk_a(ik,3)
                dtt = inp%dt_btp * init%ssprk_beta(ik)

                ! Compute reciprocal of qb_df(1,:)
                inv_pb1 = 1.0 / qb_df(1,:)

                btp%ope2_ave_df     = btp%ope2_ave_df     + (1.0 + qb_df(2,:)/init%pbprime_df)**2
                btp%uvb_ave_df(1,:) = btp%uvb_ave_df(1,:) + qb_df(3,:) * inv_pb1
                btp%uvb_ave_df(2,:) = btp%uvb_ave_df(2,:) + qb_df(4,:) * inv_pb1

                call create_rhs_btp(G, inp, b, mf, par, btp, init, ref, mpic, mt, tsp, &
                                    btp%rhs_btp, qb_df, qprime_df)
                ! Download GPU volume-integral result so CPU face-flux routine
                ! and the RK update below both see the correct rhs_btp.
                !$acc update host(btp%rhs_btp)

                ! Update barotropic variables using SSPRK formula
                qb1_df(2:4,:) = a0*qb0_df(2:4,:) + a1*qb_df(2:4,:) + a2*qb2_df(2:4,:) + dtt*btp%rhs_btp

                qb1_df(1,:) = qb1_df(2,:) + init%pbprime_df

                call btp_mom_boundary_df(G, b, mf, qb1_df)

                qb_df = qb1_df
                !$acc update device(qb_df)

                if (inp%kstages == 5 .and. ik == 2) qb2_df = qb_df

            end do

            btp%tau_wind_ave = btp%tau_wind_ave + init%tau_wind

        end do

        ! Download GPU-accumulated averaging variables to host before normalisation.
        !$acc update host(btp%H_ave, btp%Qu_ave, btp%Qv_ave, btp%Quv_ave, btp%tau_bot_ave, &
        !$acc             btp%ope_ave, btp%ope2_ave, btp%btp_mass_flux_ave, btp%uvb_ave)

        ! Normalise accumulators to get time averages
        N_inv = 1.0 / real(inp%kstages * init%N_btp)

        btp%uvb_ave_df               = N_inv * btp%uvb_ave_df
        btp%graduvb_face_ave         = N_inv * btp%graduvb_face_ave
        btp%graduvb_ave              = N_inv * btp%graduvb_ave
        btp%ope2_ave_df              = N_inv * btp%ope2_ave_df
        btp%ope2_ave                 = N_inv * btp%ope2_ave
        btp%ope_ave                  = N_inv * btp%ope_ave
        btp%H_ave                    = N_inv * btp%H_ave
        btp%Qu_ave                   = N_inv * btp%Qu_ave
        btp%Qv_ave                   = N_inv * btp%Qv_ave
        btp%Quv_ave                  = N_inv * btp%Quv_ave
        btp%btp_mass_flux_ave        = N_inv * btp%btp_mass_flux_ave
        btp%tau_bot_ave              = N_inv * btp%tau_bot_ave
        btp%ope_face_ave             = N_inv * btp%ope_face_ave
        btp%ope2_face_ave            = N_inv * btp%ope2_face_ave
        btp%H_face_ave               = N_inv * btp%H_face_ave
        btp%Qu_face_ave              = N_inv * btp%Qu_face_ave
        btp%Qv_face_ave              = N_inv * btp%Qv_face_ave
        btp%btp_mass_flux_face_ave   = N_inv * btp%btp_mass_flux_face_ave
        btp%one_plus_eta_edge_2_ave  = N_inv * btp%one_plus_eta_edge_2_ave
        btp%uvb_ave                  = N_inv * btp%uvb_ave
        btp%uvb_face_ave             = N_inv * btp%uvb_face_ave

        btp%tau_wind_ave = btp%tau_wind_ave / real(init%N_btp)

    end subroutine ti_barotropic_ssprk_mlswe


end module mod_rk_mlswe
