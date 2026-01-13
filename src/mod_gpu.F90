module mod_gpu

  use openacc
  use mpi_utilities

  implicit none

  integer :: numgpu, igpu
  logical :: is_openacc_initialized=.false.

  public :: numgpu, igpu
  public :: openacc_initialize, openacc_enter_data

contains

  subroutine openacc_initialize()

#ifdef _OPENACC
    !Get number of GPUs
    numgpu = acc_get_num_devices(acc_device_nvidia)
    if(numgpu <= 0) then
      print *,'Error in Initialize_OpenACC_Init: numgpu = ',numgpu
      call exit(-99)
    end if
    is_openacc_initialized = .true.

    !Set the GPU device mapped to this MPI rank
    igpu = mod(irank, numgpu)
    call acc_set_device_num(igpu, acc_device_nvidia)

    print*,'------------------------------------------'
    write(*,'("irank igpu nproc numgpu = ",4(i7,1x))') irank,igpu,nproc,numgpu
    print*,'------------------------------------------'
    if(igpu < 0) then
      print *,'Error in Initialize_OpenACC_Init: igpu = ',igpu
      call exit(-99)
    end if
#endif

  end subroutine openacc_initialize

  subroutine openacc_enter_data(sim)
    use mod_type_numa_simulation
    type(numa_simulation), intent(inout) :: sim(nsim_rank)

    !$acc enter data copyin(sim(1), sim(1)%q0, sim(1)%q, sim(1)%press, sim(1)%sd )
    
    !$acc enter data copyin( sim(1)%ti, &
    !$acc           sim(1)%ti%itime, sim(1)%ti%time, sim(1)%ti%Nstages, sim(1)%ti%dt, &
    !$acc           sim(1)%ti%delta, sim(1)%ti%lambda, &
    !$acc           sim(1)%ti%q1, sim(1)%ti%rhs, sim(1)%ti%RK_a, sim(1)%ti%RK_beta, &
    !$acc           sim(1)%ti%A, sim(1)%ti%At, sim(1)%ti%b, sim(1)%ti%q_tt, sim(1)%ti%b_tt, &
    !$acc           sim(1)%ti%q_ARK_hat, sim(1)%ti%q_ARK_hat_correction, sim(1)%ti%q_Stage, sim(1)%ti%q_Stage_R)

    !$acc enter data copyin(sim(1)%ti%sol, &
    !$acc           sim(1)%ti%sol%nkrylov, sim(1)%ti%sol%maxiter, sim(1)%ti%sol%etol, &
    !$acc           sim(1)%ti%sol%krylov, sim(1)%ti%sol%krylow, sim(1)%ti%sol%weight, &
    !$acc           sim(1)%ti%sol%gamma0, sim(1)%ti%sol%gamma0_norm, sim(1)%ti%sol%gamma1, &
    !$acc           sim(1)%ti%sol%solver_tol, sim(1)%ti%sol%ndof, &
    !$acc           sim(1)%ti%sol%tvec1, sim(1)%ti%sol%tvec2, sim(1)%ti%sol%tscal, sim(1)%ti%sol%tscal_sum, &
    !$acc           sim(1)%ti%sol%h, sim(1)%ti%sol%c, sim(1)%ti%sol%s )

    !$acc enter data copyin( sim(1)%inp, sim(1)%inp%nvar, &
    !$acc       sim(1)%inp%x_boundary, sim(1)%inp%y_boundary, sim(1)%inp%z_boundary, &
    !$acc       sim(1)%inp%nDirichletBC, sim(1)%inp%lmountain, sim(1)%inp%irestart )

    !$acc enter data copyin( sim(1)%sd%lc, &
    !$acc       sim(1)%sd%lc%qq, sim(1)%sd%lc%divrhou, &
    !$acc       sim(1)%sd%lc%rhoue, sim(1)%sd%lc%rhoun, sim(1)%sd%lc%rhouc, &
    !$acc       sim(1)%sd%lc%qq_e, sim(1)%sd%lc%qq_n, sim(1)%sd%lc%qq_c, &
    !$acc       sim(1)%sd%lc%qq_1, sim(1)%sd%lc%qq_2, sim(1)%sd%lc%qq_3, &
    !$acc       sim(1)%sd%lc%bq, sim(1)%sd%lc%couple, sim(1)%sd%lc%rhs_continuous, sim(1)%sd%lc%rhs_temp, &
    !$acc       sim(1)%sd%lc%q_visc, sim(1)%sd%lc%rhs_visc, sim(1)%sd%lc%temp )

    !$acc enter data copyin( sim(1)%sd%g, &
    !$acc       sim(1)%sd%g%npoin, sim(1)%sd%g%intma_cg )

    !$acc enter data copyin( sim(1)%sd%b, &
    !$acc       sim(1)%sd%b%dpsix, sim(1)%sd%b%dpsiy, sim(1)%sd%b%dpsiz, &
    !$acc       sim(1)%sd%b%dpsix_tr, sim(1)%sd%b%dpsiy_tr, sim(1)%sd%b%dpsiz_tr )

    !$acc enter data copyin( sim(1)%sd%bc, &
    !$acc       sim(1)%sd%bc%damping_parameter )

    !$acc enter data copyin( sim(1)%sd%init, &
    !$acc       sim(1)%sd%init%q_ref, sim(1)%sd%init%q_ref_initial, &
    !$acc       sim(1)%sd%init%press_ref, sim(1)%sd%init%grad_press_ref, &
    !$acc       sim(1)%sd%init%grad_rho_ref, sim(1)%sd%init%grad_theta_ref, &
    !$acc       sim(1)%sd%init%kvector, sim(1)%sd%init%coriolis_constant, &
    !$acc       sim(1)%sd%init%gravity_vector )

    !$acc enter data copyin( sim(1)%sd%v, &
    !$acc       sim(1)%sd%v%visc_on, sim(1)%sd%v%metrics_visc )

    !$acc enter data copyin( sim(1)%sd%m, &
    !$acc       sim(1)%sd%m%ksi_x, sim(1)%sd%m%ksi_y, sim(1)%sd%m%ksi_z, &
    !$acc       sim(1)%sd%m%eta_x, sim(1)%sd%m%eta_y, sim(1)%sd%m%eta_z, &
    !$acc       sim(1)%sd%m%zeta_x, sim(1)%sd%m%zeta_y, sim(1)%sd%m%zeta_z, &
    !$acc       sim(1)%sd%m%jac, sim(1)%sd%m%massinv )

    !$acc enter data copyin( sim(1)%sd%f, &
    !$acc       sim(1)%sd%f%npoin_face, sim(1)%sd%f%nid_face_1, sim(1)%sd%f%nid_face_2, sim(1)%sd%f%nid_face_3, &
    !$acc       sim(1)%sd%f%nid_face_4, sim(1)%sd%f%nid_face_5, sim(1)%sd%f%nid_face_6, & 
    !$acc       sim(1)%sd%f%normal_face, sim(1)%sd%f%normal_face_bottom )

    !$acc enter data copyin( sim(1)%sd%nbh, &
    !$acc       sim(1)%sd%nbh%num_send_recv, sim(1)%sd%nbh%num_send_recv_total, &
    !$acc       sim(1)%sd%nbh%nbh_send_recv, sim(1)%sd%nbh%nbh_proc )

    !$acc enter data copyin( sim(1)%sd%h, &
    !$acc       sim(1)%sd%h%heating, sim(1)%sd%h%nt )

  end subroutine openacc_enter_data

end module mod_gpu