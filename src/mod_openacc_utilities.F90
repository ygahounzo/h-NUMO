module mod_openacc_utilities

   use openacc
   use mod_mpi_utilities
#ifdef _OPENACC
   use mod_parallel, only: nproc
#endif
   use mod_grid,      only: grid
   use mod_basis,     only: basis
   use mod_initial,   only: initial
   use mod_metrics,   only: metrics
   use mod_face,      only: face_CS
   use mod_tensor,    only: tensor_CS
   use mod_variables, only: btp_CS, bcl_CS

   implicit none

   integer :: numgpu, igpu
   logical :: is_openacc_initialized = .false.

   public :: numgpu, igpu
   public :: openacc_initialize, openacc_enter_data

contains

   ! =========================================================================
   subroutine openacc_initialize()
      ! =========================================================================
      !  Assigns each MPI rank a GPU (round-robin) and prints a diagnostic line.
      !  Must be called after MPI_Init and before any OpenACC data regions.
      ! =========================================================================
#ifdef _OPENACC
      numgpu = acc_get_num_devices(acc_device_nvidia)
      if (numgpu <= 0) then
         print *, 'Error in openacc_initialize: no NVIDIA GPUs found, numgpu = ', numgpu
         call exit(-99)
      end if

      ! mod(irank, numgpu) is always in [0, numgpu-1] because irank >= 0
      ! and numgpu > 0 — no need to check igpu < 0 afterwards.
      igpu = mod(irank, numgpu)
      call acc_set_device_num(igpu, acc_device_nvidia)

      is_openacc_initialized = .true.

      print *, '------------------------------------------'
      write(*,'("irank igpu nproc numgpu = ",4(i7,1x))') irank, igpu, nproc, numgpu
      print *, '------------------------------------------'
#endif
   end subroutine openacc_initialize

   ! =========================================================================
   subroutine openacc_enter_data(G, b, mf, btp, bcl, init, tsp, mt)
      ! =========================================================================
      !  Uploads all static/read-only type-member arrays to the GPU once at startup.
      !  Read-write accumulation arrays (btp, bcl) are also entered here so
      !  they persist across barotropic substep calls.
      !
      !  All directives use async to overlap transfers; a blocking !$acc wait
      !  at the end ensures everything is resident before the first kernel launch.
      ! =========================================================================
#ifdef _OPENACC
      use mod_ref,      only: recv_data
      use mod_parallel, only: nbh_send_recv
#endif
      use mod_bc,       only: bc_count, bc_list

      implicit none

      type(grid),      intent(in) :: G
      type(basis),     intent(in) :: b
      type(face_CS),   intent(in) :: mf
      type(btp_CS),    intent(in) :: btp
      type(bcl_CS),    intent(in) :: bcl
      type(initial),   intent(in) :: init
      type(tensor_CS), intent(in) :: tsp
      type(metrics),   intent(in) :: mt

      ! mod_basis
      !$acc enter data copyin(b%nglx, b%ngly, b%nglz, b%ngl, b%npts,          &
      !$acc                   b%dpsix, b%dpsiy, b%dpsiz,                       &
      !$acc                   b%dpsix_tr, b%dpsiy_tr, b%dpsiz_tr) async
      !$acc enter data copyin(b%nqx, b%nqy, b%nqz, b%nq,                      &
      !$acc                   b%dpsiqx, b%dpsiqy, b%dpsiqz, b%psiq) async

      ! mod_grid
      !$acc enter data copyin(G%intma, G%face, G%intma_table, G%face_type) async

      ! mod_initial
      !$acc enter data copyin(init%kvector, init%coriolis_constant) async
      !$acc enter data copyin(tsp%wjac, tsp%psih, tsp%dpsidx, tsp%dpsidy,     &
      !$acc                   tsp%indexq,                                       &
      !$acc                   tsp%wjac_df, tsp%psih_df, tsp%dpsidx_df,         &
      !$acc                   tsp%dpsidy_df, tsp%index_df) async
      !$acc enter data copyin(init%q_df_mlswe_init, init%pbprime_df,           &
      !$acc                   init%qb_df_mlswe_init) async
      !$acc enter data copyin(init%alpha_mlswe, init%pbprime_df_face,          &
      !$acc                   init%zbot_df, init%tau_wind_df, init%zbot_face,  &
      !$acc                   init%grad_zbot_quad) async
      !$acc enter data copyin(init%tau_wind, init%coriolis_df,                 &
      !$acc                   init%coriolis_quad, init%fdt_bcl, init%fdt2_bcl, &
      !$acc                   init%a_bcl, init%b_bcl) async
      !$acc enter data copyin(btp%qb_df, bcl%qprime_df, btp%rhs_btp) async

      ! mod_metrics
      !$acc enter data copyin(mt%ksi_x, mt%ksi_y, mt%ksi_z,                   &
      !$acc                   mt%eta_x, mt%eta_y, mt%eta_z,                    &
      !$acc                   mt%zeta_x, mt%zeta_y, mt%zeta_z,                 &
      !$acc                   mt%jac, mt%xjac, mt%massinv) async
      !$acc enter data copyin(mt%ksiq_x, mt%ksiq_y, mt%ksiq_z,                &
      !$acc                   mt%etaq_x, mt%etaq_y, mt%etaq_z,                 &
      !$acc                   mt%zetaq_x, mt%zetaq_y, mt%zetaq_z, mt%jacq) async

      ! mod_ref
      ! recv_data is receive buffer: allocate on device, no initial copy needed
      !$acc enter data create(recv_data) async

      ! mod_parallel
      !$acc enter data copyin(nbh_send_recv) async

      ! mod_face
      !$acc enter data copyin(mf%imapl, mf%imapr, mf%normal_vector,           &
      !$acc                   mf%normal_vector_q, mf%jac_faceq) async

      ! mod_bc
      !$acc enter data copyin(bc_count, bc_list) async

      ! mod_variables (btp_CS)
      !$acc enter data copyin(btp%ope_ave, btp%H_ave, btp%Qu_ave, btp%Qv_ave, &
      !$acc                   btp%Quv_ave, btp%ope2_ave, btp%ope2_ave_df,      &
      !$acc                   btp%btp_mass_flux_ave, btp%uvb_ave,              &
      !$acc                   btp%ope2_face_ave) async
      !$acc enter data copyin(btp%one_plus_eta_edge_2_ave, btp%H_face_ave,     &
      !$acc                   btp%tau_wind_ave, btp%tau_bot_ave,               &
      !$acc                   btp%uvb_face_ave) async
      !$acc enter data copyin(btp%btp_mass_flux_face_ave, btp%ope_face_ave,    &
      !$acc                   btp%Qu_face_ave, btp%Qv_face_ave,                &
      !$acc                   btp%Quv_face_ave) async
      !$acc enter data copyin(btp%uvb_ave_df, btp%dpprime_visc, btp%dpp_graduv, &
      !$acc                   btp%btp_dpp_graduv, btp%pbprime_visc,            &
      !$acc                   btp%dpp_uvp) async
      !$acc enter data copyin(btp%graduvb_ave, btp%graduvb_face_ave) async

      ! Block until all transfers complete before first kernel launch
      !$acc wait

   end subroutine openacc_enter_data

end module mod_openacc_utilities
