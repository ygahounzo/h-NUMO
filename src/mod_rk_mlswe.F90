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
      integer :: mstep, ik, I
      real    :: N_inv, a0, a1, a2, dtt

      !$acc declare create(qb0_df, qb1_df, qb2_df)

      ! Zero-initialise all accumulation buffers on device.
      !$acc kernels present(btp%one_plus_eta_edge_2_ave, btp%uvb_ave, btp%uvb_ave_df,      &
      !$acc                  btp%ope_ave, btp%btp_mass_flux_ave, btp%H_ave,                &
      !$acc                  btp%Qu_ave, btp%Qv_ave, btp%Quv_ave, btp%ope2_ave_df,        &
      !$acc                  btp%uvb_face_ave, btp%ope_face_ave, btp%ope2_face_ave,        &
      !$acc                  btp%btp_mass_flux_face_ave, btp%H_face_ave, btp%Qu_face_ave,  &
      !$acc                  btp%Qv_face_ave, btp%Quv_face_ave, btp%tau_wind_ave,          &
      !$acc                  btp%tau_bot_ave, btp%ope2_ave, btp%graduvb_face_ave,          &
      !$acc                  btp%graduvb_ave)
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
      !$acc end kernels

      !$acc kernels present(qb2_df)
      qb2_df = 0.0
      !$acc end kernels

      !$acc update device(qb_df)

      ! Time loop for the barotropic solver, with SSPRK time integration.
      do mstep = 1, init%N_btp

         !$acc parallel loop present(qb0_df, qb_df)
         do I = 1, G%npoin
            qb0_df(1,I) = qb_df(1,I); qb0_df(2,I) = qb_df(2,I)
            qb0_df(3,I) = qb_df(3,I); qb0_df(4,I) = qb_df(4,I)
         end do
         !$acc end parallel loop

         do ik = 1, inp%kstages

            a0  = init%ssprk_a(ik,1)
            a1  = init%ssprk_a(ik,2)
            a2  = init%ssprk_a(ik,3)
            dtt = inp%dt_btp * init%ssprk_beta(ik)

            !$acc parallel loop present(btp, qb_df, init)
            do I = 1, G%npoin
               btp%ope2_ave_df(I)  = btp%ope2_ave_df(I)  + (1.0 + qb_df(2,I)/init%pbprime_df(I))**2
               btp%uvb_ave_df(1,I) = btp%uvb_ave_df(1,I) + qb_df(3,I) / qb_df(1,I)
               btp%uvb_ave_df(2,I) = btp%uvb_ave_df(2,I) + qb_df(4,I) / qb_df(1,I)
            end do
            !$acc end parallel loop

            call create_rhs_btp(G, inp, b, mf, par, btp, init, ref, mpic, mt, tsp, &
               btp%rhs_btp, qb_df, qprime_df)

            !$acc parallel loop present(qb0_df, qb1_df, qb_df, qb2_df, btp, init, mt)
            do I = 1, G%npoin
               qb1_df(2,I) = a0*qb0_df(2,I) + a1*qb_df(2,I) + a2*qb2_df(2,I) &
                             + dtt * (mt%massinv(I) * btp%rhs_btp(1,I))
               qb1_df(3,I) = a0*qb0_df(3,I) + a1*qb_df(3,I) + a2*qb2_df(3,I) &
                             + dtt * (mt%massinv(I) * (btp%rhs_btp(2,I) + inp%visc_mlswe*btp%rhs_btp_visc(1,I)))
               qb1_df(4,I) = a0*qb0_df(4,I) + a1*qb_df(4,I) + a2*qb2_df(4,I) &
                             + dtt * (mt%massinv(I) * (btp%rhs_btp(3,I) + inp%visc_mlswe*btp%rhs_btp_visc(2,I)))
               qb1_df(1,I) = qb1_df(2,I) + init%pbprime_df(I)
            end do
            !$acc end parallel loop

            call btp_mom_boundary_df(G, b, mf, qb1_df)

            !$acc parallel loop present(qb_df, qb1_df)
            do I = 1, G%npoin
               qb_df(1,I) = qb1_df(1,I); qb_df(2,I) = qb1_df(2,I)
               qb_df(3,I) = qb1_df(3,I); qb_df(4,I) = qb1_df(4,I)
            end do
            !$acc end parallel loop

            if (inp%kstages == 5 .and. ik == 2) then
               !$acc parallel loop present(qb2_df, qb_df)
               do I = 1, G%npoin
                  qb2_df(1,I) = qb_df(1,I); qb2_df(2,I) = qb_df(2,I)
                  qb2_df(3,I) = qb_df(3,I); qb2_df(4,I) = qb_df(4,I)
               end do
               !$acc end parallel loop
            end if

         end do

         btp%tau_wind_ave = btp%tau_wind_ave + init%tau_wind

      end do

      ! Normalise accumulators on GPU — no host round-trip needed.
      N_inv = 1.0 / real(inp%kstages * init%N_btp)

      !$acc kernels present(btp%uvb_ave_df, btp%graduvb_face_ave, btp%graduvb_ave,         &
      !$acc                  btp%ope2_ave_df, btp%ope2_ave, btp%ope_ave,                   &
      !$acc                  btp%H_ave, btp%Qu_ave, btp%Qv_ave, btp%Quv_ave,               &
      !$acc                  btp%btp_mass_flux_ave, btp%tau_bot_ave,                        &
      !$acc                  btp%ope_face_ave, btp%ope2_face_ave, btp%H_face_ave,           &
      !$acc                  btp%Qu_face_ave, btp%Qv_face_ave, btp%btp_mass_flux_face_ave, &
      !$acc                  btp%one_plus_eta_edge_2_ave, btp%uvb_ave, btp%uvb_face_ave)
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
      !$acc end kernels

      !$acc kernels present(btp%tau_wind_ave)
      btp%tau_wind_ave = btp%tau_wind_ave / real(init%N_btp)
      !$acc end kernels

      ! qb_df was updated on GPU; callers use it on CPU after this returns.
      !!$acc update host(qb_df)

   end subroutine ti_barotropic_ssprk_mlswe


end module mod_rk_mlswe
