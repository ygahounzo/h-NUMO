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
  use mod_initial_mlswe,    only: poslimiter, check_layer_thickness
  use mod_create_rhs_mlswe, only: bcl_apply_implicit_vertical_viscosity

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

  real, dimension(inp%nvar_btp,G%npoin),             intent(inout) :: qb_df
  real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers), intent(inout) :: q_df

  integer :: k, ik, I
  real :: dtt, dt_btp_in, beta_ik
  real :: rl, rm, rn, omg, Px, Py, Pz, d_inv
  real :: a, bb, tempu, tempv
  logical :: has_w
  real, parameter :: lsrk3_beta(3) = (/ 1.0/3.0, 1.0/2.0, 1.0 /)

  has_w = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

  ! Save stage-0 qb_df on GPU for the ES branch reset each stage.
  !$acc kernels present(bcl%qbp_df, qb_df)
  bcl%qbp_df = qb_df
  !$acc end kernels

  ! Initialise stage array on device.
  !$acc kernels present(bcl%q0_df, q_df)
  bcl%q0_df = q_df
  !$acc end kernels

  do ik = 1, 3

    dtt = inp%dt * lsrk3_beta(ik)

    ! GPU: compute baroclinic primed variables from the current stage state.
    call extract_qprime_df_face(G, inp, b, mt, tsp, init, bcl, bcl%qprime_df, q_df, qb_df)
    ! GPU: compute bcl/btp coefficient arrays (fused single kernel).
    call btp_bcl_coeffs_qdf(G, inp, b, tsp, bcl, btp, bcl%qprime_df, init%alpha_mlswe, init%pbprime_df)

    if (ik == 1 .and. inp%rk_bcl_FS) then
      ! FS: BTP at stage 1 only, full N_btp.
      call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, &
                                      qb_df, bcl%qprime_df, inp%dt_btp)

    elseif (.not. inp%rk_bcl_FS) then
      ! ES: restore qb_df to the stage-0 value on GPU.
      !$acc kernels present(qb_df, bcl%qbp_df)
      qb_df      = bcl%qbp_df
      !$acc end kernels
      init%N_btp = max(1, ceiling(dtt / inp%dt_btp))
      dt_btp_in = dtt / real(init%N_btp)
      call ti_barotropic_ssprk_mlswe(G, inp, b, mf, par, init, ref, mpic, mt, tsp, btp, &
                                      qb_df, bcl%qprime_df, dt_btp_in)
    endif

    ! GPU: build baroclinic RHS; rhs_bcl downloaded to host inside but device copy valid.
    call create_rhs_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, bcl%rhs_bcl, bcl%qprime_df, q_df)

    ! GPU LSRK3 update: q = q0 + dt * rhs.
    !$acc kernels present(q_df, bcl%q0_df, bcl%rhs_bcl)
    do k = 1, inp%nlayers
      q_df(1,:,k) = bcl%q0_df(1,:,k) + dtt * bcl%rhs_bcl(1,:,k)
      q_df(2,:,k) = bcl%q0_df(2,:,k) + dtt * bcl%rhs_bcl(2,:,k)
      q_df(3,:,k) = bcl%q0_df(3,:,k) + dtt * bcl%rhs_bcl(3,:,k)
      if (has_w) q_df(4,:,k) = bcl%q0_df(4,:,k) + dtt * bcl%rhs_bcl(4,:,k)
    end do
    !$acc end kernels

    ! Semi-implicit Coriolis rotation (GPU), always applied — Coriolis is
    ! never added explicitly in the RHS (mod_create_rhs_mlswe.F90). Applies
    ! (I + ω K_r)^{-1} acting on the LSRK3-updated momentum, where the
    ! "P-vector" blends the updated momentum with the stage-0 cross-product
    ! term (bcl%q0_df, fixed for the whole stage loop, present on device).
    beta_ik = lsrk3_beta(ik)
    if (has_w) then
      ! 3D spherical rotation.
      !$acc parallel loop gang collapse(2) &
      !$acc    present(init%kvector, init%fdt2_bcl, q_df, bcl%q0_df) &
      !$acc    private(rl, rm, rn, omg, Px, Py, Pz, d_inv) &
      !$acc    firstprivate(beta_ik)
      do k = 1, inp%nlayers
        do I = 1, G%npoin
          rl    = init%kvector(1,I)
          rm    = init%kvector(2,I)
          rn    = init%kvector(3,I)
          omg   = init%fdt2_bcl(I) * beta_ik
          Px    = q_df(2,I,k) + omg*(rn*bcl%q0_df(3,I,k) - rm*bcl%q0_df(4,I,k))
          Py    = q_df(3,I,k) + omg*(rl*bcl%q0_df(4,I,k) - rn*bcl%q0_df(2,I,k))
          Pz    = q_df(4,I,k) + omg*(rm*bcl%q0_df(2,I,k) - rl*bcl%q0_df(3,I,k))
          d_inv = 1.0 / (1.0 + omg**2)
          q_df(2,I,k) = (Px - omg*(rm*Pz - rn*Py)) * d_inv
          q_df(3,I,k) = (Py - omg*(rn*Px - rl*Pz)) * d_inv
          q_df(4,I,k) = (Pz - omg*(rl*Py - rm*Px)) * d_inv
        end do
      end do
      !$acc end parallel loop
    else
      ! 2D Cartesian beta-plane rotation — has_w=(0,0,1) specialization of
      ! the 3D formula above.
      !$acc parallel loop gang collapse(2) &
      !$acc    present(init%fdt2_bcl, q_df, bcl%q0_df) &
      !$acc    private(omg, a, bb, tempu, tempv) &
      !$acc    firstprivate(beta_ik)
      do k = 1, inp%nlayers
        do I = 1, G%npoin
          omg   = init%fdt2_bcl(I) * beta_ik
          a     = 1.0 / (1.0 + omg**2)
          bb    = omg / (1.0 + omg**2)
          tempu = q_df(2,I,k) + omg*bcl%q0_df(3,I,k)
          tempv = q_df(3,I,k) - omg*bcl%q0_df(2,I,k)
          q_df(2,I,k) =  a*tempu + bb*tempv
          q_df(3,I,k) = -bb*tempu + a*tempv
        end do
      end do
      !$acc end parallel loop
    end if

    ! GPU: wall BC — gang over faces, atomic updates for corner nodes.
    call layer_mom_boundary_df(G, inp, b, mf, init, q_df)

    ! Floor thickness before extract_velocity: if q_df(1) went negative from the
    ! RK update, the division by thickness inside extract_velocity would produce
    ! astronomical velocities (denominator = thickness + 1e-20 ≈ 1e-20 when
    ! thickness < 0), corrupting momentum and triggering a NaN cascade.
    call poslimiter(b, G, inp, mt, q_df, init%alpha_mlswe)

    ! GPU: extract baroclinic velocity (removes barotropic component).
    call extract_velocity(G, inp, b, mt, tsp, init, bcl, bcl%uv_df, q_df, qb_df)

    ! GPU: reconstruct momentum from corrected velocity.
    !$acc kernels present(q_df, bcl%uv_df)
    do k = 1, inp%nlayers
      q_df(2,:,k) = bcl%uv_df(1,:,k) * q_df(1,:,k)
      q_df(3,:,k) = bcl%uv_df(2,:,k) * q_df(1,:,k)
      if (has_w) q_df(4,:,k) = bcl%uv_df(3,:,k) * q_df(1,:,k)
    end do
    !$acc end kernels

    ! GPU: Zhang-Shu positivity limiter (element-local, no cross-element races).
    call poslimiter(b, G, inp, mt, q_df, init%alpha_mlswe)

    call check_layer_thickness(b, G, inp, q_df, 'ti_lsrk3_bcl', ik)

  end do

  ! Inter-layer (ad_mlswe) vertical viscosity, applied implicitly once per
  ! full step rather than embedded in the RHS above.
  ! No-op internally when inp%ad_mlswe<=0.
  ! call bcl_apply_implicit_vertical_viscosity(G, inp, b, mf, init, tsp, mt, bcl, q_df, qb_df)

end subroutine ti_lsrk3_bcl
