
module mod_openacc_utilities
    
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



subroutine openacc_enter_data()
  use mod_basis, only: nglx,ngly,nglz,dpsix,dpsiy,dpsiz
  use mod_basis, only: nqx,nqy,nq,dpsiqx,dpsiqy,dpsiqz
  use mod_grid, only: intma,intma_table,face
  use mod_initial, only: kvector, coriolis_constant
  use mod_metrics, only: ksi_x,ksi_y,ksi_z,eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,jac,xjac, massinv
  use mod_metrics, only: ksiq_x,ksiq_y,ksiq_z,etaq_x,etaq_y,etaq_z,zetaq_x,zetaq_y,zetaq_z,jacq
  use mod_ref, only: recv_data
  use mod_parallel, only: nbh_send_recv
  use mod_face, only: imapl, normal_vector, imapr, normal_vector_q
  use mod_bc, only: bc_count, bc_list
  use mod_initial, only: wjac, psih, dpsidx, dpsidy, indexq, wjac_df, psih_df, dpsidx_df, dpsidy_df, index_df
  use mod_initial, only: q_df_mlswe_init, pbprime_df, qb_df_mlswe_init
  use mod_initial, only: alpha_mlswe, pbprime_df_face, zbot_df, tau_wind_df, zbot_face, grad_zbot_quad
  use mod_initial, only: tau_wind, coriolis_df, coriolis_quad, fdt_bcl, fdt2_bcl, a_bcl, b_bcl
  use mod_variables, only: ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave, ope2_ave_df, btp_mass_flux_ave, uvb_ave, ope2_face_ave
  use mod_variables, only: one_plus_eta_edge_2_ave, H_face_ave, tau_wind_ave, tau_bot_ave, uvb_face_ave
  use mod_variables, only: btp_mass_flux_face_ave, ope_face_ave, Qu_face_ave, Qv_face_ave, Quv_face_ave
  use mod_variables, only: dpprime_visc, dpp_graduv, btp_dpp_graduv, pbprime_visc, dpp_uvp
  use mod_variables, only: graduvb_ave, graduvb_face_ave

  implicit none


  !$acc enter data copyin(nglx,ngly,nglz,ngl,npts,nqx,nqy,nqz,nq)
  !$acc enter data copyin(dpsix,dpsiy,dpsiz,dpsiqx,dpsiqy,dpsiqz)

  ! mod_grid
  !$acc enter data copyin(intma,face,intma_table)
  ! mod_initial

  !$acc enter data copyin(kvector,coriolis_constant)
  !$acc enter data copyin(wjac, psih, dpsidx, dpsidy, indexq, wjac_df, psih_df, dpsidx_df, dpsidy_df, index_df)
  !$acc enter data copyin(q_df_mlswe_init, pbprime_df, qb_df_mlswe_init)
  !$acc enter data copyin(alpha_mlswe, pbprime_df_face, zbot_df, tau_wind_df, zbot_face, grad_zbot_quad)
  !$acc enter data copyin(tau_wind, coriolis_df, coriolis_quad, fdt_bcl, fdt2_bcl, a_bcl, b_bcl)
  ! mod_metrics

  !$acc enter data copyin(ksi_x,ksi_y,ksi_z,eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,jac,xjac,massinv)
  !$acc enter data copyin(ksiq_x,ksiq_y,ksiq_z,etaq_x,etaq_y,etaq_z,zetaq_x,zetaq_y,zetaq_z,jacq)

  ! mod_ref
  !$acc enter data create(recv_data)
  ! mod_parallel
  !$acc enter data copyin(nbh_send_recv)

  ! mod_face
  !$acc enter data copyin(imapl, imapr, normal_vector, normal_vector_q)

  ! mod_bc
  !$acc enter data copyin(bc_count, bc_list)

  ! mod_variables
  !$acc enter data copyin(ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave, ope2_ave_df, btp_mass_flux_ave, uvb_ave, ope2_face_ave)
  !$acc enter data copyin(one_plus_eta_edge_2_ave, H_face_ave, tau_wind_ave, tau_bot_ave, uvb_face_ave) 
  !$acc enter data copyin(btp_mass_flux_face_ave, ope_face_ave, Qu_face_ave, Qv_face_ave, Quv_face_ave)
  !$acc enter data copyin(uvb_ave_df, dpprime_visc, dpp_graduv, btp_dpp_graduv, pbprime_visc, dpp_uvp)
  !$acc enter data copyin(graduvb_ave, graduvb_face_ave)

end subroutine openacc_enter_data

end module mod_openacc_utilities
