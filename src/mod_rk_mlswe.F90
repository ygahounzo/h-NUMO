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
      qb_df, qprime_df, dt_btp_in)

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
      use mod_laplacian_quad,    only: btp_create_laplacian

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

      real, dimension(inp%nvar_btp,G%npoin),             intent(inout) :: qb_df
      real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers), intent(in)    :: qprime_df
      real,                                               intent(in)    :: dt_btp_in

      integer :: mstep, ik, I, iv, nvarb_f
      real    :: N_inv, a0, a1, a2, dtt, visc_term
      logical :: has_w

      has_w   = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex
      nvarb_f = inp%nvar_btp

      ! Zero-initialise all accumulation buffers on device.
      !$acc kernels present(btp%one_plus_eta_edge_2_ave, btp%uvb_ave, btp%uvb_ave_df,      &
      !$acc                  btp%ope_ave, btp%btp_mass_flux_ave, btp%H_ave,                &
      !$acc                  btp%Qu_ave, btp%Qv_ave, btp%Qw_ave, btp%ope2_ave_df,          &
      !$acc                  btp%uvb_face_ave, btp%ope_face_ave, btp%ope2_face_ave,        &
      !$acc                  btp%btp_mass_flux_face_ave, btp%H_face_ave, btp%Qu_face_ave,  &
      !$acc                  btp%Qv_face_ave, btp%Qw_face_ave, btp%tau_wind_ave,           &
      !$acc                  btp%tau_bot_ave, btp%ope2_ave, btp%graduvb_face_ave,          &
      !$acc                  btp%graduvb_ave, btp%qb2_df, btp%rhs_btp_visc)
      btp%one_plus_eta_edge_2_ave = 0.0;  btp%uvb_ave             = 0.0
      btp%uvb_ave_df              = 0.0;  btp%ope_ave             = 0.0
      btp%btp_mass_flux_ave       = 0.0;  btp%H_ave               = 0.0
      btp%Qu_ave                  = 0.0;  btp%Qv_ave              = 0.0
      btp%Qw_ave                  = 0.0;  btp%ope2_ave_df         = 0.0
      btp%uvb_face_ave            = 0.0;  btp%ope_face_ave        = 0.0
      btp%ope2_face_ave           = 0.0;  btp%btp_mass_flux_face_ave = 0.0
      btp%H_face_ave              = 0.0;  btp%Qu_face_ave         = 0.0
      btp%Qv_face_ave             = 0.0;  btp%Qw_face_ave         = 0.0
      btp%tau_wind_ave            = 0.0;  btp%tau_bot_ave         = 0.0
      btp%ope2_ave                = 0.0;  btp%graduvb_face_ave    = 0.0
      btp%graduvb_ave             = 0.0;  btp%qb2_df              = 0.0
      btp%rhs_btp_visc            = 0.0
      !$acc end kernels

      ! Frozen viscosity: compute rhs_btp_visc once per BCL stage using the
      ! initial qb_df, then hold it fixed across all BTP substeps.
      ! Valid when the viscous diffusion number nu*dt_btp/dx^2 << 1, which is
      ! satisfied whenever dt_btp is set by the gravity-wave CFL (dt_btp ~ dx/c):
      !   nu*dt_btp/dx^2 ~ nu/(c*dx)
      ! For c~200 m/s, dx>=1 km, nu<=500: this ratio is ~2.5e-3 — negligible.
      ! This eliminates N_btp*kstages-1 Laplacian solves and MPI halo exchanges.
      if (inp%method_visc > 0) &
         call btp_create_laplacian(G, inp, b, mf, par, btp, init, ref, mpic, tsp, &
                                   btp%rhs_btp_visc, qb_df)

      ! Time loop for the barotropic solver, with SSPRK time integration.
      do mstep = 1, init%N_btp

         !$acc parallel loop present(btp%qb0_df, qb_df) private(iv) firstprivate(nvarb_f)
         do I = 1, G%npoin
            !$acc loop seq
            do iv = 1, nvarb_f
               btp%qb0_df(iv,I) = qb_df(iv,I)
            end do
         end do
         !$acc end parallel loop

         do ik = 1, inp%kstages

            a0  = init%ssprk_a(ik,1)
            a1  = init%ssprk_a(ik,2)
            a2  = init%ssprk_a(ik,3)
            dtt = dt_btp_in * init%ssprk_beta(ik)

            call create_rhs_btp(G, inp, b, mf, par, btp, init, ref, mpic, mt, tsp, &
               btp%rhs_btp, qb_df, qprime_df)

            ! Fused: accumulate averages from old qb_df, then SSPRK update in-place.
            ! Reading qb_df(I) for averages and a1-coefficient happens before the
            ! write to qb_df(I) within the same thread — no race across nodes.
            !$acc parallel loop present(btp%qb0_df, qb_df, btp%qb2_df, btp, init, mt) &
            !$acc    private(iv, visc_term) firstprivate(a0, a1, a2, dtt, nvarb_f)
            do I = 1, G%npoin
               btp%ope2_ave_df(I) = btp%ope2_ave_df(I) + (1.0 + qb_df(2,I)/init%pbprime_df(I))**2

               ! Momentum components: u,v[,w] -> uvb_ave_df(1..nvarb_f-2), rhs_btp(2..nvarb_f-1).
               ! Viscosity (rhs_btp_visc) is only defined for u,v (indices 3,4); w has none yet.
               !$acc loop seq
               do iv = 3, nvarb_f
                  btp%uvb_ave_df(iv-2,I) = btp%uvb_ave_df(iv-2,I) + qb_df(iv,I) / qb_df(1,I)
                  visc_term = 0.0
                  if (iv <= 4) visc_term = inp%visc_mlswe * btp%rhs_btp_visc(iv-2,I)
                  qb_df(iv,I) = a0*btp%qb0_df(iv,I) + a1*qb_df(iv,I) + a2*btp%qb2_df(iv,I) &
                               + dtt * (mt%massinv(I) * (btp%rhs_btp(iv-1,I) + visc_term))
               end do

               qb_df(2,I) = a0*btp%qb0_df(2,I) + a1*qb_df(2,I) + a2*btp%qb2_df(2,I) &
                            + dtt * (mt%massinv(I) * btp%rhs_btp(1,I))
               qb_df(1,I) = qb_df(2,I) + init%pbprime_df(I)
            end do
            !$acc end parallel loop

            call btp_mom_boundary_df(G, b, mf, init, inp, qb_df)

            if (inp%kstages == 5 .and. ik == 2) then
               !$acc kernels present(btp%qb2_df, qb_df)
               btp%qb2_df = qb_df
               !$acc end kernels
            end if

         end do

      end do

      ! tau_wind is constant over sub-steps; skip the broken CPU accumulation loop.
      !$acc kernels present(btp%tau_wind_ave, init%tau_wind)
      btp%tau_wind_ave = init%tau_wind
      !$acc end kernels

      ! Normalise accumulators on GPU — no host round-trip needed.
      N_inv = 1.0 / real(inp%kstages * init%N_btp)

      !$acc kernels present(btp%uvb_ave_df, btp%ope2_ave_df, btp%ope2_ave, btp%ope_ave,    &
      !$acc                  btp%H_ave, btp%Qu_ave, btp%Qv_ave, btp%Qw_ave,                &
      !$acc                  btp%btp_mass_flux_ave, btp%tau_bot_ave,                        &
      !$acc                  btp%ope_face_ave, btp%ope2_face_ave, btp%H_face_ave,           &
      !$acc                  btp%Qu_face_ave, btp%Qv_face_ave, btp%Qw_face_ave,             &
      !$acc                  btp%btp_mass_flux_face_ave,                                    &
      !$acc                  btp%one_plus_eta_edge_2_ave, btp%uvb_ave, btp%uvb_face_ave)
      btp%uvb_ave_df               = N_inv * btp%uvb_ave_df
      btp%ope2_ave_df              = N_inv * btp%ope2_ave_df
      btp%ope2_ave                 = N_inv * btp%ope2_ave
      btp%ope_ave                  = N_inv * btp%ope_ave
      btp%H_ave                    = N_inv * btp%H_ave
      btp%Qu_ave                   = N_inv * btp%Qu_ave
      btp%Qv_ave                   = N_inv * btp%Qv_ave
      btp%Qw_ave                   = N_inv * btp%Qw_ave
      btp%btp_mass_flux_ave        = N_inv * btp%btp_mass_flux_ave
      btp%tau_bot_ave              = N_inv * btp%tau_bot_ave
      btp%ope_face_ave             = N_inv * btp%ope_face_ave
      btp%ope2_face_ave            = N_inv * btp%ope2_face_ave
      btp%H_face_ave               = N_inv * btp%H_face_ave
      btp%Qu_face_ave              = N_inv * btp%Qu_face_ave
      btp%Qv_face_ave              = N_inv * btp%Qv_face_ave
      btp%Qw_face_ave              = N_inv * btp%Qw_face_ave
      btp%btp_mass_flux_face_ave   = N_inv * btp%btp_mass_flux_face_ave
      btp%one_plus_eta_edge_2_ave  = N_inv * btp%one_plus_eta_edge_2_ave
      btp%uvb_ave                  = N_inv * btp%uvb_ave
      btp%uvb_face_ave             = N_inv * btp%uvb_face_ave
      !$acc end kernels

   end subroutine ti_barotropic_ssprk_mlswe


end module mod_rk_mlswe
