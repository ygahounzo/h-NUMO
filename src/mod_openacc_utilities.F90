module mod_openacc_utilities

   use openacc
   use mod_mpi_utilities
   use mod_grid,      only: grid
   use mod_basis,     only: basis
   use mod_initial,   only: initial
   use mod_metrics,   only: metrics
   use mod_face,      only: face_CS
   use mod_tensor,    only: tensor_CS
   use mod_variables, only: btp_CS, bcl_CS
   use mod_parallel,  only: parallel_CS
   use mod_ref,       only: mref

   implicit none

   integer :: numgpu, igpu
   logical :: is_openacc_initialized = .false.

   public :: numgpu, igpu
   public :: openacc_initialize, openacc_enter_data

contains

   ! =========================================================================
   subroutine openacc_initialize(par)
      ! =========================================================================
      !  Assigns each MPI rank a GPU (round-robin) and prints a diagnostic line.
      !  Must be called after MPI_Init and before any OpenACC data regions.
      ! =========================================================================
      type(parallel_CS), intent(in) :: par

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
      write(*,'("irank igpu nproc numgpu = ",4(i7,1x))') irank, igpu, par%nproc, numgpu
      print *, '------------------------------------------'
#endif
   end subroutine openacc_initialize

   ! =========================================================================
   subroutine openacc_enter_data(G, b, mf, btp, bcl, init, tsp, mt, par, ref)
      ! =========================================================================
      !  Uploads all static/read-only type-member arrays to the GPU once at startup.
      !  Read-write accumulation arrays (btp, bcl) are also entered here so
      !  they persist across barotropic substep calls.
      !
      !  All directives use async to overlap transfers; a blocking !$acc wait
      !  at the end ensures everything is resident before the first kernel launch.
      ! =========================================================================
      use mod_bc, only: bc_count, bc_list

      implicit none

      type(grid),        intent(in) :: G
      type(basis),       intent(in) :: b
      type(face_CS),     intent(in) :: mf
      type(btp_CS),      intent(in) :: btp
      type(bcl_CS),      intent(in) :: bcl
      type(initial),     intent(in) :: init
      type(tensor_CS),   intent(in) :: tsp
      type(metrics),     intent(in) :: mt
      type(parallel_CS), intent(in) :: par
      type(mref),        intent(in) :: ref

      ! NVHPC two-step deep copy: first copyin the parent struct (puts scalars
      ! and array descriptors on the device), then copyin each allocatable
      ! member array (creates device memory and patches device-side descriptors).
      ! Without the first step, kernel access to struct%array dereferences an
      ! invalid device pointer and causes CUDA error 700.

      ! mod_basis (scalar members are firstprivate in each kernel — NOT entered)
      !$acc enter data copyin(b)
      !$acc enter data copyin(b%dpsix, b%dpsiy, b%dpsiz,                       &
      !$acc                   b%dpsix_tr, b%dpsiy_tr, b%dpsiz_tr)
      !$acc enter data copyin(b%dpsiqx, b%dpsiqy, b%dpsiqz, b%psiq)
      !$acc enter data copyin(b%psi)
      !$acc enter data copyin(b%wglx, b%wgly)

      ! mod_grid
      !$acc enter data copyin(G)
      !$acc enter data copyin(G%intma, G%face, G%intma_table, G%face_type)

      ! mod_tensor
      !$acc enter data copyin(tsp)
      !$acc enter data copyin(tsp%wjac, tsp%psih, tsp%dpsidx, tsp%dpsidy,     &
      !$acc                   tsp%indexq,  tsp%indexq_e,                       &
      !$acc                   tsp%wjac_df, tsp%psih_df, tsp%dpsidx_df,         &
      !$acc                   tsp%dpsidy_df, tsp%index_df, tsp%index_df_elt)

      ! mod_initial
      !$acc enter data copyin(init)
      !$acc enter data copyin(init%kvector, init%coriolis_constant)
      !$acc enter data copyin(init%q_df, init%pbprime_df,                      &
      !$acc                   init%qb_df)
      !$acc enter data copyin(init%alpha_mlswe, init%pbprime_df_face,          &
      !$acc                   init%zbot_df, init%tau_wind_df, init%zbot_face,  &
      !$acc                   init%grad_zbot_quad)
      !$acc enter data copyin(init%tau_wind, init%coriolis_df,                 &
      !$acc                   init%coriolis_quad, init%fdt_bcl, init%fdt2_bcl, &
      !$acc                   init%a_bcl, init%b_bcl)

      ! mod_metrics
      !$acc enter data copyin(mt)
      !$acc enter data copyin(mt%ksi_x, mt%ksi_y, mt%ksi_z,                   &
      !$acc                   mt%eta_x, mt%eta_y, mt%eta_z,                    &
      !$acc                   mt%zeta_x, mt%zeta_y, mt%zeta_z,                 &
      !$acc                   mt%jac, mt%xjac, mt%massinv)
      !$acc enter data copyin(mt%ksiq_x, mt%ksiq_y, mt%ksiq_z,                &
      !$acc                   mt%etaq_x, mt%etaq_y, mt%etaq_z,                 &
      !$acc                   mt%zetaq_x, mt%zetaq_y, mt%zetaq_z, mt%jacq)

      ! mod_face
      !$acc enter data copyin(mf)
      !$acc enter data copyin(mf%imapl, mf%imapr, mf%normal_vector,           &
      !$acc                   mf%normal_vector_q, mf%jac_faceq, mf%jac_face)

      ! mod_variables (btp_CS)
      !$acc enter data copyin(btp)
      !$acc enter data copyin(btp%qb_df, btp%rhs_btp, btp%rhs_btp_visc)
      !$acc enter data create(btp%qb0_df, btp%qb2_df)
      !$acc enter data copyin(btp%ope_ave, btp%H_ave, btp%Qu_ave, btp%Qv_ave, &
      !$acc                   btp%Qw_ave,                                      &
      !$acc                   btp%ope2_ave, btp%ope2_ave_df,      &
      !$acc                   btp%btp_mass_flux_ave, btp%uvb_ave,              &
      !$acc                   btp%ope2_face_ave)
      !$acc enter data copyin(btp%one_plus_eta_edge_2_ave, btp%H_face_ave,     &
      !$acc                   btp%tau_wind_ave, btp%tau_bot_ave,               &
      !$acc                   btp%uvb_face_ave)
      !$acc enter data copyin(btp%btp_mass_flux_face_ave, btp%ope_face_ave,    &
      !$acc                   btp%Qu_face_ave, btp%Qv_face_ave,                &
      !$acc                   btp%Qw_face_ave)
      !$acc enter data copyin(btp%uvb_ave_df, btp%btp_dpp_graduvw,             &
      !$acc                   btp%pbprime_visc)
      !$acc enter data copyin(btp%graduvb_ave, btp%graduvb_face_ave)
      !$acc enter data create(btp%bcl_H, btp%bcl_flux, btp%bcl_btp_flux, btp%pbq)

      ! mod_variables (bcl_CS)
      !$acc enter data copyin(bcl)
      !$acc enter data create(bcl%rhs_bcl)
      !$acc enter data create(bcl%rhs_visc_bcl)
      !$acc enter data copyin(bcl%q_df)
      !$acc enter data copyin(bcl%qprime_df, bcl%dpprime_visc, bcl%dpp_graduvw, &
      !$acc                   bcl%dpp_uvp)
      !$acc enter data create(bcl%q0_df, bcl%q1_df, bcl%uv_df, bcl%qbp_df)

      ! mod_ref
      ! recv_data is a receive buffer: allocate on device, no initial copy needed.
      ! q_send/q_recv hold packed MPI boundary data; updated each sub-step via
      ! !$acc update device after unpack so create_nbhs_face_df can run on GPU.
      !$acc enter data copyin(ref)
      !$acc enter data create(ref%recv_data)
      !$acc enter data copyin(ref%q_send, ref%q_recv)
      !$acc enter data create(ref%send_data_dg)
      !$acc enter data create(ref%recv_data_dg)
      !$acc enter data create(ref%send_data_dg_lap)
      !$acc enter data create(ref%recv_data_dg_lap)
      !$acc enter data create(ref%q_send_lap, ref%q_recv_lap)
      !$acc enter data copyin(ref%face_pack_list)
      !$acc enter data create(ref%send_data_bcl, ref%recv_data_bcl)
      !$acc enter data create(ref%q_send_bcl, ref%q_recv_bcl)
      !$acc enter data create(ref%send_data_lap_bcl, ref%recv_data_lap_bcl)
      !$acc enter data create(ref%q_send_lap_bcl, ref%q_recv_lap_bcl)

      ! mod_parallel
      !$acc enter data copyin(par)
      !$acc enter data copyin(par%nbh_send_recv, par%num_send_recv, par%nbh_send_recv_multi)

      ! mod_bc
      !$acc enter data copyin(bc_count, bc_list)

      ! Block until all transfers complete before first kernel launch
      !$acc wait

   end subroutine openacc_enter_data

end module mod_openacc_utilities
