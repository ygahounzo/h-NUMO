module mod_openacc_utilities

    use openacc
    use mpi_utilities

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
    subroutine openacc_enter_data()
      ! =========================================================================
      !  Uploads all static/read-only module arrays to the GPU once at startup.
      !  Read-write accumulation arrays (mod_variables) are also entered here so
      !  they persist across barotropic substep calls.
      !
      !  Subroutine arguments (rhs, qb, qprime_df, …) are NOT entered here;
      !  the caller is responsible for managing those per-call.
      !
      !  All directives use async to overlap transfers; a blocking !$acc wait
      !  at the end ensures everything is resident before the first kernel launch.
      ! =========================================================================
        use mod_basis, only: nglx, ngly, nglz, ngl, npts,            &
                              dpsix, dpsiy, dpsiz,                    &
                              dpsix_tr, dpsiy_tr, dpsiz_tr
        use mod_basis, only: nqx, nqy, nqz, nq,                      &
                              dpsiqx, dpsiqy, dpsiqz,                 &
                              psiq                      
        use mod_grid, only: intma, intma_table, face,                 &
                             face_type                
        use mod_initial, only: kvector, coriolis_constant
        use mod_initial, only: wjac, psih, dpsidx, dpsidy, indexq,   &
                                wjac_df, psih_df, dpsidx_df, dpsidy_df, index_df
        use mod_initial, only: q_df_mlswe_init, pbprime_df, qb_df_mlswe_init
        use mod_initial, only: alpha_mlswe, pbprime_df_face,          &
                                zbot_df, tau_wind_df, zbot_face, grad_zbot_quad
        use mod_initial, only: tau_wind, coriolis_df, coriolis_quad,  &
                                fdt_bcl, fdt2_bcl, a_bcl, b_bcl

        use mod_metrics, only: ksi_x, ksi_y, ksi_z,                  &
                                eta_x, eta_y, eta_z,                  &
                                zeta_x, zeta_y, zeta_z,               &
                                jac, xjac, massinv
        use mod_metrics, only: ksiq_x, ksiq_y, ksiq_z,               &
                                etaq_x, etaq_y, etaq_z,               &
                                zetaq_x, zetaq_y, zetaq_z, jacq

        use mod_ref,      only: recv_data
        use mod_parallel, only: nbh_send_recv
        use mod_face, only: imapl, imapr, normal_vector, normal_vector_q, &
                             jac_faceq                 
        use mod_bc, only: bc_count, bc_list
        use mod_variables, only: ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave,    &
                                  ope2_ave, ope2_ave_df, btp_mass_flux_ave,   &
                                  uvb_ave, ope2_face_ave
        use mod_variables, only: one_plus_eta_edge_2_ave, H_face_ave,         &
                                  tau_wind_ave, tau_bot_ave, uvb_face_ave
        use mod_variables, only: btp_mass_flux_face_ave, ope_face_ave,        &
                                  Qu_face_ave, Qv_face_ave, Quv_face_ave
        use mod_variables, only: uvb_ave_df, dpprime_visc, dpp_graduv,        &
                                  btp_dpp_graduv, pbprime_visc, dpp_uvp
        use mod_variables, only: graduvb_ave, graduvb_face_ave

        implicit none

        ! mod_basis
        !$acc enter data copyin(nglx, ngly, nglz, ngl, npts,                  &
        !$acc                   dpsix, dpsiy, dpsiz,                           &
        !$acc                   dpsix_tr, dpsiy_tr, dpsiz_tr) async
        !$acc enter data copyin(nqx, nqy, nqz, nq,                            &
        !$acc                   dpsiqx, dpsiqy, dpsiqz, psiq) async

        ! mod_grid
        !$acc enter data copyin(intma, face, intma_table, face_type) async

        ! mod_initial
        !$acc enter data copyin(kvector, coriolis_constant) async
        !$acc enter data copyin(wjac, psih, dpsidx, dpsidy, indexq,           &
        !$acc                   wjac_df, psih_df, dpsidx_df, dpsidy_df,       &
        !$acc                   index_df) async
        !$acc enter data copyin(q_df_mlswe_init, pbprime_df,                  &
        !$acc                   qb_df_mlswe_init) async
        !$acc enter data copyin(alpha_mlswe, pbprime_df_face,                 &
        !$acc                   zbot_df, tau_wind_df, zbot_face,              &
        !$acc                   grad_zbot_quad) async
        !$acc enter data copyin(tau_wind, coriolis_df, coriolis_quad,         &
        !$acc                   fdt_bcl, fdt2_bcl, a_bcl, b_bcl) async

        ! mod_metrics
        !$acc enter data copyin(ksi_x, ksi_y, ksi_z,                          &
        !$acc                   eta_x, eta_y, eta_z,                           &
        !$acc                   zeta_x, zeta_y, zeta_z, jac, xjac, massinv) async
        !$acc enter data copyin(ksiq_x, ksiq_y, ksiq_z,                       &
        !$acc                   etaq_x, etaq_y, etaq_z,                        &
        !$acc                   zetaq_x, zetaq_y, zetaq_z, jacq) async

        ! mod_ref
        ! recv_data is receive buffer: allocate on device, no initial copy needed
        !$acc enter data create(recv_data) async

        ! mod_parallel
        !$acc enter data copyin(nbh_send_recv) async

        ! mod_face
        !$acc enter data copyin(imapl, imapr, normal_vector,                  &
        !$acc                   normal_vector_q, jac_faceq) async

        ! mod_bc
        !$acc enter data copyin(bc_count, bc_list) async

        ! mod_variables
        ! These are read-write: entered with copyin so existing values
        ! are preserved; updated in place on device.
        !$acc enter data copyin(ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave,      &
        !$acc                   ope2_ave, ope2_ave_df, btp_mass_flux_ave,      &
        !$acc                   uvb_ave, ope2_face_ave) async
        !$acc enter data copyin(one_plus_eta_edge_2_ave, H_face_ave,           &
        !$acc                   tau_wind_ave, tau_bot_ave, uvb_face_ave) async
        !$acc enter data copyin(btp_mass_flux_face_ave, ope_face_ave,          &
        !$acc                   Qu_face_ave, Qv_face_ave, Quv_face_ave) async
        !$acc enter data copyin(uvb_ave_df, dpprime_visc, dpp_graduv,          &
        !$acc                   btp_dpp_graduv, pbprime_visc, dpp_uvp) async
        !$acc enter data copyin(graduvb_ave, graduvb_face_ave) async

        ! Block until all transfers complete before first kernel launch
        !$acc wait

    end subroutine openacc_enter_data

end module mod_openacc_utilities