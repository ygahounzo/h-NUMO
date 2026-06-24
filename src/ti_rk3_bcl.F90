! ===========================================================================================================================
! This module contains the routines for the predictor-corrector (for baroclinic) and the RK35 time integration (for barotropic) methods
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: October 27, 2023
! ==========================================================================================================================

subroutine ti_rk3_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, q_df, qb_df)

  ! q: layer variable dp, u*dp, v*dp at quad points and their face values: q_face
  ! q_df : layer variable dp, u*dp, v*dp at nodal (dof) point
  ! qprime: value dp', u' and v' at quad points and their face values: qprime_face
  ! qb : barotopic variable pb, pb_pert = pb'*eta, ub*pb, vb*pb at quad points and their face values: qb_face
  ! qb_df: : barotopic variable pb, pb_pert = pb'*eta, ub*pb, vb*pb at nodal points
  ! qprime_df: value dp', u' and v' at nodal points
  ! qp_df_out: output variable, thickness h_k, velocity u_k,v_k, free surface ssh

  use mod_grid,              only: grid
  use mod_input,             only: input
  use mod_basis,             only: basis
  use mod_face,              only: face_CS
  use mod_parallel,          only: parallel_CS
  use mod_initial,           only: initial
  use mod_ref,               only: mref
  use mod_mpi_communicator,  only: mpi_communicator
  use mod_metrics,           only: metrics
  use mod_tensor,         only: tensor_CS
  use mod_variables,         only: bcl_CS, btp_CS
  use mod_splitting,         only: create_rhs_bcl
  use mod_rk_mlswe,          only: ti_barotropic_ssprk_mlswe
  use mod_barotropic_terms,  only: btp_bcl_coeffs_qdf
  use mod_layer_terms,       only: extract_qprime_df_face, layer_mom_boundary_df, extract_velocity
  use mod_initial_mlswe,    only: poslimiter, find_dry_elements
  use mod_constants,        only: gravity
  use ieee_arithmetic,      only: ieee_is_nan

  implicit none

  type(grid),             intent(in)    :: G
  type(input),            intent(in)    :: inp
  type(basis),            intent(in)    :: b
  type(face_CS),          intent(in)    :: mf
  type(parallel_CS),      intent(in)    :: par
  type(initial),          intent(inout) :: init
  type(mref),             intent(inout) :: ref
  type(mpi_communicator), intent(inout)    :: mpic
  type(metrics),          intent(in)    :: mt
  type(tensor_CS),        intent(in)    :: tsp
  type(bcl_CS),           intent(inout) :: bcl
  type(btp_CS),           intent(inout) :: btp

  real, dimension(4,G%npoin),             intent(inout) :: qb_df
  real, dimension(3,G%npoin,inp%nlayers), intent(inout) :: q_df

  real, dimension(4,G%npoin)             :: qbp_df
  integer :: k, ik, I
  real, dimension(3,G%npoin,inp%nlayers) :: q0_df, q1_df, rhs
  real, dimension(2,G%npoin,inp%nlayers) :: uv_df
  real :: dtt

  q0_df = q_df
  q1_df = q_df

  qbp_df = qb_df

  do ik = 1, inp%kstages_bcl

    dtt = inp%dt * init%ssprk_beta_bcl(ik)

    call find_dry_elements(G, inp, b, tsp, btp%ope2_ave_df, bcl%q_df, &
                           init%alpha_mlswe, bcl%dry_flg)

    call extract_qprime_df_face(G, inp, b, mt, init, bcl, bcl%qprime_df, q1_df, qb_df)

    call btp_bcl_coeffs_qdf(G, inp, b, tsp, bcl, btp, bcl%qprime_df)

    ! Always sync — create_rhs_bcl reads bcl%qprime_df from device every stage.
    !$acc update device(bcl%qprime_df)

    if (ik == 1 .and. inp%rk_bcl_FS) then
      ! FS: BTP at stage 1 only, full N_btp
      call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, &
                                      qb_df, bcl%qprime_df)

    elseif (.not. inp%rk_bcl_FS) then
      ! ES: BTP at every stage with N_btp scaled to the stage dt
      qb_df      = qbp_df
      init%N_btp = max(1, ceiling(dtt / inp%dt_btp))
      call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, &
                                      qb_df, bcl%qprime_df)
    endif

    call create_rhs_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs, bcl%qprime_df, q1_df)

    do k = 1, inp%nlayers
      q_df(1,:,k) = init%ssprk_a_bcl(ik,1)*q0_df(1,:,k) + init%ssprk_a_bcl(ik,2)*q1_df(1,:,k) + dtt*rhs(1,:,k)
      q_df(2,:,k) = init%ssprk_a_bcl(ik,1)*q0_df(2,:,k) + init%ssprk_a_bcl(ik,2)*q1_df(2,:,k) + dtt*rhs(2,:,k)
      q_df(3,:,k) = init%ssprk_a_bcl(ik,1)*q0_df(3,:,k) + init%ssprk_a_bcl(ik,2)*q1_df(3,:,k) + dtt*rhs(3,:,k)
    end do

    call layer_mom_boundary_df(G, inp, b, mf, q_df)

    ! Limit dp before velocity extraction so dry-cell nodes with nonzero momentum
    ! from the DG integral do not produce large velocities via division.
    if (inp%lposlimiter) call poslimiter(b, G, inp, mt, q_df, init%alpha_mlswe)

    ! Zero momentum wherever the velocity after the BCL update would be unphysical.
    ! Use the same combined cap as the qprime clamp above:
    !   smaller of (10x local wave speed) or (200 m/s absolute).
    ! For q_df: u = q(2)/q(1), c^2 = alpha*q(1), so u^2 > (10c)^2 means
    !   q(2)^2 + q(3)^2 > (10c)^2 * q(1)^2 = 100 * alpha * q(1)^3.
    ! do k = 1, inp%nlayers
    !   do ie = 1, G%nelem
    !     if (bcl%dry_flg(ie, k) == 2) cycle
    !     Iq = tsp%indexq_e(1, ie)
    !     do ip = 1, b%npts
    !       I = tsp%indexq(ip, Iq)
    !       if (q_df(1,I,k) > 0.0) then
    !         if (q_df(2,I,k)**2 + q_df(3,I,k)**2 > &
    !             min((10.0)**2 * init%alpha_mlswe(k) * q_df(1,I,k)**3, &
    !                 (200.0 * q_df(1,I,k))**2)) then
    !           q_df(2,I,k) = 0.0
    !           q_df(3,I,k) = 0.0
    !         end if
    !       end if
    !     end do
    !   end do
    ! end do

    ! Compute dpprime, uprime and vprime at the quad and nodal points
    call extract_velocity(G, inp, b, mt, init, bcl, uv_df, q_df, qb_df)

    do k = 1,inp%nlayers
      q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
      q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
    end do

    q1_df = q_df

  end do
end subroutine ti_rk3_bcl
