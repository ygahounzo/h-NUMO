! ===========================================================================================================================
! LSRK3 (Wicker & Skamarock, 2002) time integration for the baroclinic (BCL) mode:
!   phi^(n+1/3) = phi^n + (dt/3) * F(phi^n)
!   phi^(n+1/2) = phi^n + (dt/2) * F(phi^(n+1/3))
!   phi^(n+1)   = phi^n + dt     * F(phi^(n+1/2))
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
! ==========================================================================================================================

subroutine ti_lsrk3_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, q_df, qb_df)

  use mod_grid,              only: grid
  use mod_input,             only: input
  use mod_basis,             only: basis
  use mod_face,              only: face_CS
  use mod_parallel,          only: parallel_CS
  use mod_initial,           only: initial
  use mod_ref,               only: mref
  use mod_mpi_communicator,  only: mpi_communicator
  use mod_metrics,           only: metrics
  use mod_tensor,            only: tensor_CS
  use mod_variables,         only: bcl_CS, btp_CS
  use mod_splitting,         only: create_rhs_bcl
  use mod_rk_mlswe,          only: ti_barotropic_ssprk_mlswe
  use mod_barotropic_terms,  only: btp_bcl_coeffs_qdf
  use mod_layer_terms,       only: extract_qprime_df_face, layer_mom_boundary_df, extract_velocity
  use mod_initial_mlswe,    only: poslimiter, find_dry_elements
  use mod_constants,        only: gravity

  implicit none

  type(grid),             intent(in)    :: G
  type(input),            intent(in)    :: inp
  type(basis),            intent(in)    :: b
  type(face_CS),          intent(in)    :: mf
  type(parallel_CS),      intent(in)    :: par
  type(initial),          intent(inout) :: init
  type(mref),             intent(inout) :: ref
  type(mpi_communicator), intent(inout) :: mpic
  type(metrics),          intent(in)    :: mt
  type(tensor_CS),        intent(in)    :: tsp
  type(bcl_CS),           intent(inout) :: bcl
  type(btp_CS),           intent(inout) :: btp

  real, dimension(4,G%npoin),             intent(inout) :: qb_df
  real, dimension(3,G%npoin,inp%nlayers), intent(inout) :: q_df

  real, dimension(4,G%npoin)             :: qbp_df
  integer :: k, ik, I
  real, dimension(3,G%npoin,inp%nlayers) :: q0_df, rhs
  real, dimension(2,G%npoin,inp%nlayers) :: uv_df
  real :: dtt
  real, parameter :: lsrk3_beta(3) = (/ 1.0/3.0, 1.0/2.0, 1.0 /)

  q0_df  = q_df
  qbp_df = qb_df

  do ik = 1, 3

    dtt = inp%dt * lsrk3_beta(ik)

    call find_dry_elements(G, inp, b, tsp, btp%ope2_ave_df, q_df, &
                           init%alpha_mlswe, bcl%dry_flg)

    call extract_qprime_df_face(G, inp, b, mt, init, bcl, bcl%qprime_df, q_df, qb_df)
    call btp_bcl_coeffs_qdf(G, inp, b, tsp, bcl, btp, bcl%qprime_df)

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

    call create_rhs_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs, bcl%qprime_df, q_df)

    do k = 1, inp%nlayers
      q_df(1,:,k) = q0_df(1,:,k) + dtt * rhs(1,:,k)
      q_df(2,:,k) = q0_df(2,:,k) + dtt * rhs(2,:,k)
      q_df(3,:,k) = q0_df(3,:,k) + dtt * rhs(3,:,k)
    end do

    call layer_mom_boundary_df(G, inp, b, mf, q_df)

    if (inp%lposlimiter) call poslimiter(b, G, inp, mt, q_df, init%alpha_mlswe)

    ! Zero momentum at dry-threshold nodes to prevent large velocities
    ! from pressure-gradient forces on nearly-dry cells in mixed elements.
    do k = 1, inp%nlayers
      do I = 1, G%npoin
        if (q_df(1,I,k) <= (gravity/init%alpha_mlswe(k)) * inp%dry_cutoff) then
          q_df(2,I,k) = 0.0
          q_df(3,I,k) = 0.0
        end if
      end do
    end do

    call extract_velocity(G, inp, b, mt, init, bcl, uv_df, q_df, qb_df)

    do k = 1, inp%nlayers
      q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
      q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
    end do

  end do
end subroutine ti_lsrk3_bcl
