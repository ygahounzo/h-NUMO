! =================================================================================================
! This module contains the routines for the barotropic and baroclinic horizontal viscosity terms
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: July 24, 2023
! =================================================================================================

module mod_laplacian_quad

   implicit none

   public :: bcl_create_laplacian, btp_create_laplacian, compute_bcl_lap_z, &
             bcl_compute_nu_smag, btp_compute_nu_smag, diagnose_min_grid_scale

contains

   !----------------------------------------------------------------------!
   !>@brief One-time startup diagnostic: scans every element in the mesh
   !> (across all ranks) for the smallest element-edge-based filter width
   !> Delta = min(dx_len, dy_len)/(N+1) -- same convention as the paper's
   !> Delta_bar before Eq. 6 -- and reports the resulting explicit
   !> diffusive-CFL ceiling A_H <= Delta^2/(16*dt) for both the BCL and BTP
   !> viscous operators.  This is the actual visc_mlswe ceiling for THIS
   !> mesh: the domain-average element size is not representative when the
   !> mesh has non-uniform element quality (e.g. sphere_ico), and exceeding
   !> the ceiling at even a single worst-case element causes a fast, violent
   !> linear blow-up localized at that element -- as opposed to running out
   !> of (nonlinear, slower) dissipation for a wave-breaking event.
   !----------------------------------------------------------------------!
   subroutine diagnose_min_grid_scale(G, b, tsp, inp)

      use mod_grid,          only: grid
      use mod_basis,         only: basis
      use mod_tensor,        only: tensor_CS
      use mod_input,         only: input
      use mod_mpi_utilities, only: irank, irank0, MPI_PRECISION
      use mpi

      implicit none

      type(grid),      intent(in) :: G
      type(basis),     intent(in) :: b
      type(tensor_CS), intent(in) :: tsp
      type(input),     intent(in) :: inp

      integer :: ie, nglx_l, ngly_l, npts_l
      integer :: I00, IN0, I0N
      real    :: dx_len, dy_len, delta, delta2
      real    :: delta2_min_local, delta2_min_global
      real    :: mu_cfl_bcl, mu_cfl_btp
      integer :: ierr_mpi

      nglx_l = b%nglx
      ngly_l = b%ngly
      npts_l = b%npts

      delta2_min_local = huge(1.0)

      do ie = 1, G%nelem
         I00 = tsp%index_df_elt(1, ie)
         IN0 = tsp%index_df_elt(nglx_l, ie)
         I0N = tsp%index_df_elt((ngly_l-1)*nglx_l + 1, ie)

         dx_len = sqrt((G%coord(1,I00)-G%coord(1,IN0))**2 + &
                       (G%coord(2,I00)-G%coord(2,IN0))**2 + &
                       (G%coord(3,I00)-G%coord(3,IN0))**2)
         dy_len = sqrt((G%coord(1,I00)-G%coord(1,I0N))**2 + &
                       (G%coord(2,I00)-G%coord(2,I0N))**2 + &
                       (G%coord(3,I00)-G%coord(3,I0N))**2)

         delta  = min(dx_len/real(nglx_l), dy_len/real(ngly_l))
         delta2 = delta*delta

         delta2_min_local = min(delta2_min_local, delta2)
      end do

      call mpi_allreduce(delta2_min_local, delta2_min_global, 1, MPI_PRECISION, &
                          mpi_min, mpi_comm_world, ierr_mpi)

      mu_cfl_bcl = delta2_min_global / (16.0 * inp%dt)
      mu_cfl_btp = delta2_min_global / (16.0 * inp%dt_btp)

      if (irank == irank0) then
         write(*,'(A)') ' ========================================================================'
         write(*,'(A)') ' Grid-scale diagnostic (worst-case element across whole mesh, all ranks):'
         write(*,'(A,ES14.6,A)') '   min Delta^2  = ', delta2_min_global, ' m^2'
         write(*,'(A,ES14.6,A)') '   min Delta    = ', sqrt(delta2_min_global), ' m'
         write(*,'(A,ES14.6,A,ES14.6,A)') '   explicit visc_mlswe ceiling (BCL, dt=', inp%dt, 's): ', mu_cfl_bcl, ' m^2/s'
         write(*,'(A,ES14.6,A,ES14.6,A)') '   explicit visc_mlswe ceiling (BTP, dt_btp=', inp%dt_btp, 's): ', mu_cfl_btp, ' m^2/s'
         write(*,'(A)') ' ========================================================================'
      end if

   end subroutine diagnose_min_grid_scale

   subroutine btp_create_laplacian(G, inp, b, mf, par, btp, init, ref, mpic, tsp, rhs_btp_visc, qb_df)

      use mod_grid,             only: grid
      use mod_input,            only: input
      use mod_basis,            only: basis
      use mod_face,             only: face_CS
      use mod_parallel,         only: parallel_CS
      use mod_variables,        only: btp_CS
      use mod_initial,          only: initial
      use mod_ref,              only: mref
      use mod_mpi_communicator, only: mpi_communicator
      use mod_tensor,           only: tensor_CS

      implicit none

      type(grid),             intent(in)    :: G
      type(input),            intent(in)    :: inp
      type(basis),            intent(in)    :: b
      type(face_CS),          intent(in)    :: mf
      type(parallel_CS),      intent(in)    :: par
      type(btp_CS),           intent(inout) :: btp
      type(initial),          intent(in)    :: init
      type(mref),             intent(inout) :: ref
      type(mpi_communicator), intent(inout) :: mpic
      type(tensor_CS),        intent(in)    :: tsp

      real, intent(out) :: rhs_btp_visc(inp%nvar_btp-2,G%npoin)
      real, dimension(inp%nvar_btp,G%npoin), intent(in) :: qb_df

      real, dimension(inp%nvar_btp-2,G%npoin) :: Uk
      real, dimension(inp%ngraduvw_var,G%npoin) :: graduv
      integer :: I, Iq, ip, iw, nw_btp_l
      real    :: dhdx, dhdy, dhdz, A_H_btp
      real    :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz
      logical :: has_w

      has_w = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex

      !$acc data create(Uk, graduv)

      ! Compute barotropic velocity — qb_df lives on device during the BTP loop.
      !$acc parallel loop present(Uk, qb_df) firstprivate(has_w)
      do I = 1, G%npoin
         Uk(1,I) = qb_df(3,I) / qb_df(1,I)
         Uk(2,I) = qb_df(4,I) / qb_df(1,I)
         if (has_w) Uk(3,I) = qb_df(5,I) / qb_df(1,I)
      end do
      !$acc end parallel loop

      ! Compute velocity gradient (inlined from compute_gradient_uv), using the
      ! same shell-embedding correction (dpsidz_df_x/y/z, from mt%zeta_x/y/z)
      ! as btp_bcl_coeffs_qdf so graduv/dpp_graduvw stay index-consistent.
      !$acc kernels present(graduv)
      graduv = 0.0
      !$acc end kernels

      !$acc parallel loop present(graduv, Uk, tsp, b, G) &
      !$acc    private(dhdx, dhdy, dhdz, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, &
      !$acc            dw_dx, dw_dy, dw_dz, I) firstprivate(has_w)
      do Iq = 1, G%npoin
         du_dx = 0.0; du_dy = 0.0; du_dz = 0.0
         dv_dx = 0.0; dv_dy = 0.0; dv_dz = 0.0
         dw_dx = 0.0; dw_dy = 0.0; dw_dz = 0.0
         !$acc loop seq
         do ip = 1, b%npts
            I    = tsp%index_df(ip,Iq)
            dhdx = tsp%dpsidx_df(ip,Iq)
            dhdy = tsp%dpsidy_df(ip,Iq)
            dhdz = 0.0
            if (has_w) then
               dhdx = dhdx + tsp%dpsidz_df_x(ip,Iq)
               dhdy = dhdy + tsp%dpsidz_df_y(ip,Iq)
               dhdz = tsp%dpsidz_df(ip,Iq) + tsp%dpsidz_df_z(ip,Iq)
            end if
            du_dx = du_dx + dhdx*Uk(1,I)
            du_dy = du_dy + dhdy*Uk(1,I)
            dv_dx = dv_dx + dhdx*Uk(2,I)
            dv_dy = dv_dy + dhdy*Uk(2,I)
            if (has_w) then
               du_dz = du_dz + dhdz*Uk(1,I)
               dv_dz = dv_dz + dhdz*Uk(2,I)
               dw_dx = dw_dx + dhdx*Uk(3,I)
               dw_dy = dw_dy + dhdy*Uk(3,I)
               dw_dz = dw_dz + dhdz*Uk(3,I)
            end if
         end do
         graduv(1,Iq) = du_dx; graduv(2,Iq) = du_dy
         graduv(3,Iq) = dv_dx; graduv(4,Iq) = dv_dy
         if (has_w) then
            graduv(5,Iq) = du_dz; graduv(6,Iq) = dv_dz
            graduv(7,Iq) = dw_dx; graduv(8,Iq) = dw_dy; graduv(9,Iq) = dw_dz
         end if
      end do
      !$acc end parallel loop

      ! Smagorinsky eddy viscosity
      if (inp%C_smag > 0.0) &
         call btp_compute_nu_smag(G, inp, btp, tsp, graduv)

      ! Accumulate the graduv time-average, and write pbprime_visc_scaled/
      ! btp_dpp_graduvw_scaled = A_H×(pbprime_visc/btp_dpp_graduvw), A_H =
      ! (visc + nu_smag), for correct ∇·(A_H×Δp×∇u) with {A_H×Δp} face
      ! averaging.  Writing into separate _scaled fields (rather than
      ! mutating pbprime_visc/btp_dpp_graduvw in place) avoids needing to
      ! restore them afterward -- those are shared fields, fixed for the
      ! whole BTP sub-step loop by btp_bcl_coeffs_qdf, so corrupting them
      ! here would corrupt every subsequent BTP sub-stage this BCL step.
      ! When C_smag=0, A_H≡visc_mlswe is a true spatial constant and this is
      ! a no-op relative to scaling rhs_btp_visc afterward instead; it
      ! matters once nu_smag varies in space.
      nw_btp_l = inp%ngraduvw_var
      !$acc parallel loop gang present(btp, graduv, btp%pbprime_visc, btp%pbprime_visc_scaled, &
      !$acc                            btp%btp_dpp_graduvw, btp%btp_dpp_graduvw_scaled, btp%nu_smag) &
      !$acc    private(A_H_btp, iw) firstprivate(has_w, nw_btp_l)
      do I = 1, G%npoin
         btp%graduvb_ave(1,I) = btp%graduvb_ave(1,I) + graduv(1,I)
         btp%graduvb_ave(2,I) = btp%graduvb_ave(2,I) + graduv(2,I)
         btp%graduvb_ave(3,I) = btp%graduvb_ave(3,I) + graduv(3,I)
         btp%graduvb_ave(4,I) = btp%graduvb_ave(4,I) + graduv(4,I)
         if (has_w) then
            btp%graduvb_ave(5,I) = btp%graduvb_ave(5,I) + graduv(5,I)
            btp%graduvb_ave(6,I) = btp%graduvb_ave(6,I) + graduv(6,I)
            btp%graduvb_ave(7,I) = btp%graduvb_ave(7,I) + graduv(7,I)
            btp%graduvb_ave(8,I) = btp%graduvb_ave(8,I) + graduv(8,I)
            btp%graduvb_ave(9,I) = btp%graduvb_ave(9,I) + graduv(9,I)
         end if

         A_H_btp = inp%visc_mlswe + btp%nu_smag(I)
         btp%pbprime_visc_scaled(I) = A_H_btp * btp%pbprime_visc(I)
         !$acc loop seq
         do iw = 1, nw_btp_l
             btp%btp_dpp_graduvw_scaled(iw,I) = A_H_btp * btp%btp_dpp_graduvw(iw,I)
         end do
      end do
      !$acc end parallel loop

      ! MPI precommunicator reads graduv from host — download once before packing.
      !$acc update host(graduv)

      call btp_lap_create_precommunicator(G, b, mf, init, par, btp, ref, mpic, graduv, inp%ngraduvw_var, Uk, inp%nvar_btp-2)

      ! rhs_btp_visc already zeroed on device by caller (create_rhs_btp).
      call btp_compute_laplacian_qp(G, inp, b, btp, tsp, rhs_btp_visc, graduv)
      call create_rhs_laplacian_flux(G, inp, b, mf, btp, tsp, Uk, rhs_btp_visc, graduv)

      call create_rhs_lap_postcommunicator_df(G, inp, b, mf, par, btp, ref, mpic, tsp, rhs_btp_visc, inp%nvar_btp-2)

      ! rhs_btp_visc stays on device; read by the SSPRK GPU kernel in the caller.
      !$acc end data

   end subroutine btp_create_laplacian

   subroutine bcl_create_laplacian(G, inp, b, mf, par, btp, bcl, ref, mpic, mt, tsp, rhs_lap)

      use mod_grid,             only: grid
      use mod_input,            only: input
      use mod_basis,            only: basis
      use mod_face,             only: face_CS
      use mod_parallel,         only: parallel_CS
      use mod_variables,        only: btp_CS, bcl_CS
      use mod_ref,              only: mref
      use mod_mpi_communicator, only: mpi_communicator
      use mod_metrics,          only: metrics
      use mod_tensor,           only: tensor_CS

      implicit none

      type(grid),             intent(in)    :: G
      type(input),            intent(in)    :: inp
      type(basis),            intent(in)    :: b
      type(face_CS),          intent(in)    :: mf
      type(parallel_CS),      intent(in)    :: par
      type(btp_CS),           intent(in)    :: btp
      type(bcl_CS),           intent(inout) :: bcl
      type(mref),             intent(inout) :: ref
      type(mpi_communicator), intent(inout) :: mpic
      type(metrics),          intent(in)    :: mt
      type(tensor_CS),        intent(in)    :: tsp

      real, intent(out) :: rhs_lap(inp%nvar_bcl-1,G%npoin,inp%nlayers)

      integer :: k, I, iw, nw_bcl_l, nlayers_l, npoin_l
      real    :: A_H_bcl
      logical :: has_w

      has_w = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex
      nw_bcl_l  = inp%ngraduvw_var
      nlayers_l = inp%nlayers
      npoin_l   = G%npoin

      ! rhs_lap is created on device; all kernels run GPU-to-GPU.
      !$acc data create(rhs_lap)

      ! Smagorinsky eddy viscosity -- must run before the scaling below, which
      ! reads bcl%nu_smag; only depends on the already-stable btp%graduvb_ave
      ! (BTP sub-stepping for this stage has already completed) and the raw
      ! (unscaled) dpprime_visc/dpp_graduvw.
      if (inp%C_smag > 0.0) call bcl_compute_nu_smag(G, inp, b, btp, bcl, tsp)

      ! Write dpprime_visc_scaled/dpp_graduvw_scaled = A_H*(dpprime_visc/dpp_graduvw),
      ! A_H = (visc + nu_smag), for correct ∇·(A_H×Δp×∇u) with {A_H×Δp} face
      ! averaging (DG-correct form, same as btp_create_laplacian).  Writing into
      ! separate _scaled fields avoids mutating the shared dpprime_visc/dpp_graduvw
      ! in place.  When C_smag=0, A_H≡visc_mlswe is a true spatial constant and
      ! this is a no-op relative to scaling rhs_lap afterward instead.
      !$acc parallel loop gang collapse(2) &
      !$acc    present(bcl%dpprime_visc, bcl%dpprime_visc_scaled, bcl%dpp_graduvw, &
      !$acc             bcl%dpp_graduvw_scaled, bcl%nu_smag) &
      !$acc    private(A_H_bcl, iw) firstprivate(nlayers_l, npoin_l, nw_bcl_l)
      do k = 1, nlayers_l
         do I = 1, npoin_l
            A_H_bcl = inp%visc_mlswe + bcl%nu_smag(I,k)
            bcl%dpprime_visc_scaled(I,k) = A_H_bcl * bcl%dpprime_visc(I,k)
            !$acc loop seq
            do iw = 1, nw_bcl_l
                bcl%dpp_graduvw_scaled(iw,I,k) = A_H_bcl * bcl%dpp_graduvw(iw,I,k)
            end do
         end do
      end do
      !$acc end parallel loop

      ! Pre-comm GPU pack reads the scaled dpp_graduvw/dpprime_visc from device,
      ! so the DG-correct {A_H×Δp} face averaging is consistent across MPI
      ! boundaries too.
      call bcl_lap_create_precommunicator(G, inp, b, mf, par, ref, mpic, bcl%dpp_graduvw_scaled, bcl%dpprime_visc_scaled, bcl%qprime_df)
      call bcl_compute_laplacian(G, inp, b, btp, bcl, tsp, rhs_lap)
      call bcl_create_rhs_laplacian_flux(G, inp, b, mf, btp, bcl, tsp, rhs_lap)
      ! CPU mpi_waitall inside post-comm, followed by GPU unpack + GPU face scatter.
      call bcl_create_rhs_lap_postcommunicator_df(G, inp, b, mf, par, btp, ref, mpic, tsp, rhs_lap)

      ! GPU: apply mass-inverse scaling while rhs_lap is still on device.
      ! Viscosity is already baked into rhs_lap via the pre-scaled dpprime_visc/
      ! dpp_graduvw above.
      !$acc kernels present(rhs_lap, mt%massinv) firstprivate(has_w)
      do k = 1, inp%nlayers
         rhs_lap(1,:,k) = mt%massinv(:)*rhs_lap(1,:,k)
         rhs_lap(2,:,k) = mt%massinv(:)*rhs_lap(2,:,k)
         if (has_w) rhs_lap(3,:,k) = mt%massinv(:)*rhs_lap(3,:,k)
      end do
      !$acc end kernels

      !$acc update host(rhs_lap)
      !$acc end data

   end subroutine bcl_create_laplacian

   subroutine btp_compute_laplacian(G, b, btp, tsp, rhs_btp_visc, grad_dpuvp)

      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_variables, only: btp_CS
      use mod_tensor,    only: tensor_CS

      implicit none

      type(grid),      intent(in)  :: G
      type(basis),     intent(in)  :: b
      type(btp_CS),    intent(in)  :: btp
      type(tensor_CS), intent(in)  :: tsp

      real, intent(out) :: rhs_btp_visc(2,G%npoin)
      real, dimension(4,G%npoin), intent(in) :: grad_dpuvp

      integer :: Iq, I, ip, ie, ijq, npts_f, nelem_f
      real :: wq, qq(4)
      real :: lap_q_loc(2, b%npts)

      npts_f  = b%npts
      nelem_f = G%nelem

      ! One gang per element.  Accumulate into gang-private lap_q_loc, then
      ! scatter once to global lap_q.  DG nodes belong to exactly one element
      ! so the scatter is a direct assignment — no atomics, no races.
      !$acc parallel loop gang private(lap_q_loc, qq, wq, Iq, I, ip, ijq) &
      !$acc    present(tsp%wjac_df, tsp%index_df, tsp%index_df_elt,        &
      !$acc            tsp%dpsidx_df, tsp%dpsidy_df,                        &
      !$acc            btp%pbprime_visc, btp%btp_dpp_graduvw,               &
      !$acc            grad_dpuvp, rhs_btp_visc)                                   &
      !$acc    firstprivate(npts_f, nelem_f)
      do ie = 1, nelem_f

         !$acc loop seq
         do ip = 1, npts_f
            lap_q_loc(1, ip) = 0.0
            lap_q_loc(2, ip) = 0.0
         end do

         !$acc loop seq
         do ijq = 1, npts_f
            Iq    = tsp%index_df_elt(ijq, ie)
            wq    = tsp%wjac_df(Iq)
            qq(1) = btp%pbprime_visc(Iq)*grad_dpuvp(1,Iq) + btp%btp_dpp_graduvw(1,Iq)
            qq(2) = btp%pbprime_visc(Iq)*grad_dpuvp(2,Iq) + btp%btp_dpp_graduvw(2,Iq)
            qq(3) = btp%pbprime_visc(Iq)*grad_dpuvp(3,Iq) + btp%btp_dpp_graduvw(3,Iq)
            qq(4) = btp%pbprime_visc(Iq)*grad_dpuvp(4,Iq) + btp%btp_dpp_graduvw(4,Iq)
            !$acc loop seq
            do ip = 1, npts_f
               lap_q_loc(1,ip) = lap_q_loc(1,ip) - &
                  wq*(tsp%dpsidx_df(ip,Iq)*qq(1) + tsp%dpsidy_df(ip,Iq)*qq(2))
               lap_q_loc(2,ip) = lap_q_loc(2,ip) - &
                  wq*(tsp%dpsidx_df(ip,Iq)*qq(3) + tsp%dpsidy_df(ip,Iq)*qq(4))
            end do
         end do

         ! index_df(ip, Iq) is element-local: same for any Iq in element ie.
         Iq = tsp%index_df_elt(1, ie)
         !$acc loop seq
         do ip = 1, npts_f
            I = tsp%index_df(ip, Iq)
            rhs_btp_visc(1, I) = lap_q_loc(1, ip)
            rhs_btp_visc(2, I) = lap_q_loc(2, ip)
         end do

      end do
      !$acc end parallel loop

   end subroutine btp_compute_laplacian

   subroutine btp_compute_laplacian_qp(G, inp, b, btp, tsp, rhs_btp_visc, grad_dpuvp)

      use mod_grid,      only: grid
      use mod_input,     only: input
      use mod_basis,     only: basis
      use mod_variables, only: btp_CS
      use mod_tensor,    only: tensor_CS

      implicit none

      type(grid),      intent(in)  :: G
      type(input),     intent(in)  :: inp
      type(basis),     intent(in)  :: b
      type(btp_CS),    intent(in)  :: btp
      type(tensor_CS), intent(in)  :: tsp

      real, intent(out) :: rhs_btp_visc(inp%nvar_btp-2,G%npoin)
      real, dimension(inp%ngraduvw_var,G%npoin), intent(in) :: grad_dpuvp

      integer :: Iq, I, ip
      real :: wq, qq(9), dhdx, dhdy, dhdz
      integer :: npoin_l, npts_l
      logical :: has_w

      npoin_l = G%npoin
      npts_l  = b%npts
      has_w   = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex

      !$acc data present(tsp%wjac_df, tsp%dpsidx_df, tsp%dpsidy_df, tsp%index_df,  &
      !$acc              tsp%dpsidz_df, tsp%dpsidz_df_x, tsp%dpsidz_df_y, tsp%dpsidz_df_z, &
      !$acc              btp%pbprime_visc_scaled, btp%btp_dpp_graduvw_scaled, grad_dpuvp, &
      !$acc              rhs_btp_visc)

      !$acc kernels present(rhs_btp_visc)
      rhs_btp_visc = 0.0
      !$acc end kernels

      !$acc parallel loop gang                                     &
      !$acc   private(qq, wq, I, ip, dhdx, dhdy, dhdz)              &
      !$acc   firstprivate(npoin_l, npts_l, has_w)
      do Iq = 1, npoin_l

         wq    = tsp%wjac_df(Iq)
         qq(1) = btp%pbprime_visc_scaled(Iq)*grad_dpuvp(1,Iq) + btp%btp_dpp_graduvw_scaled(1,Iq)
         qq(2) = btp%pbprime_visc_scaled(Iq)*grad_dpuvp(2,Iq) + btp%btp_dpp_graduvw_scaled(2,Iq)
         qq(3) = btp%pbprime_visc_scaled(Iq)*grad_dpuvp(3,Iq) + btp%btp_dpp_graduvw_scaled(3,Iq)
         qq(4) = btp%pbprime_visc_scaled(Iq)*grad_dpuvp(4,Iq) + btp%btp_dpp_graduvw_scaled(4,Iq)
         if (has_w) then
            qq(5) = btp%pbprime_visc_scaled(Iq)*grad_dpuvp(5,Iq) + btp%btp_dpp_graduvw_scaled(5,Iq)
            qq(6) = btp%pbprime_visc_scaled(Iq)*grad_dpuvp(6,Iq) + btp%btp_dpp_graduvw_scaled(6,Iq)
            qq(7) = btp%pbprime_visc_scaled(Iq)*grad_dpuvp(7,Iq) + btp%btp_dpp_graduvw_scaled(7,Iq)
            qq(8) = btp%pbprime_visc_scaled(Iq)*grad_dpuvp(8,Iq) + btp%btp_dpp_graduvw_scaled(8,Iq)
            qq(9) = btp%pbprime_visc_scaled(Iq)*grad_dpuvp(9,Iq) + btp%btp_dpp_graduvw_scaled(9,Iq)
         end if

         !$acc loop seq
         do ip = 1, npts_l
            I = tsp%index_df(ip,Iq)
            dhdx = tsp%dpsidx_df(ip,Iq)
            dhdy = tsp%dpsidy_df(ip,Iq)
            dhdz = 0.0
            if (has_w) then
               dhdx = dhdx + tsp%dpsidz_df_x(ip,Iq)
               dhdy = dhdy + tsp%dpsidz_df_y(ip,Iq)
               dhdz = tsp%dpsidz_df(ip,Iq) + tsp%dpsidz_df_z(ip,Iq)
            end if
            !$acc atomic update
            rhs_btp_visc(1,I) = rhs_btp_visc(1,I) - wq*(dhdx*qq(1) + dhdy*qq(2))
            !$acc atomic update
            rhs_btp_visc(2,I) = rhs_btp_visc(2,I) - wq*(dhdx*qq(3) + dhdy*qq(4))
            if (has_w) then
               !$acc atomic update
               rhs_btp_visc(1,I) = rhs_btp_visc(1,I) - wq*dhdz*qq(5)
               !$acc atomic update
               rhs_btp_visc(2,I) = rhs_btp_visc(2,I) - wq*dhdz*qq(6)
               !$acc atomic update
               rhs_btp_visc(3,I) = rhs_btp_visc(3,I) - wq*(dhdx*qq(7) + dhdy*qq(8) + dhdz*qq(9))
            end if
         end do
      end do
      !$acc end parallel loop

      !$acc end data

   end subroutine btp_compute_laplacian_qp

   ! Compute Smagorinsky eddy viscosity ν_s = (C_s·Δ)²·|S| at every BCL DOF node.
   ! Uses the velocity gradient already assembled in bcl%dpp_graduvw by the
   ! pre-communicator phase.  Result written to bcl%nu_smag (device-resident).
   subroutine bcl_compute_nu_smag(G, inp, b, btp, bcl, tsp)

      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_input,     only: input
      use mod_variables, only: btp_CS, bcl_CS
      use mod_tensor,    only: tensor_CS

      implicit none

      type(grid),      intent(in)    :: G
      type(basis),     intent(in)    :: b
      type(input),     intent(in)    :: inp
      type(btp_CS),    intent(in)    :: btp
      type(bcl_CS),    intent(inout) :: bcl
      type(tensor_CS), intent(in)    :: tsp

      integer :: ie, iq_local, Iq, Iq0, I, ip, k
      integer :: npts_l, nelem_l, nlayers_l
      real    :: dp_inv, g11, g12, g21, g22, g13, g23, g31, g32, g33
      real    :: two_S2, wjac_sum, Csm2
      real    :: nu_loc(b%npts, inp%nlayers)
      logical :: has_w
      real, parameter :: eps1 = 1.0e-20

      npts_l    = b%npts
      nelem_l   = G%nelem
      nlayers_l = inp%nlayers
      has_w     = (inp%nvar_bcl == 4)
      Csm2      = inp%C_smag**2

      !$acc data present(tsp%wjac_df, tsp%index_df_elt, tsp%index_df,               &
      !$acc              bcl%dpprime_visc, bcl%dpp_graduvw,                          &
      !$acc              btp%graduvb_ave, bcl%nu_smag)

      !$acc parallel loop gang                                                        &
      !$acc   private(nu_loc, wjac_sum, dp_inv,                                      &
      !$acc           g11, g12, g21, g22, g13, g23, g31, g32, g33,                  &
      !$acc           two_S2, Iq, Iq0, I, ip, k, iq_local)                          &
      !$acc   firstprivate(npts_l, nelem_l, nlayers_l, has_w, Csm2)
      do ie = 1, nelem_l

         ! Element scale: (C_s·Δ)² where Δ² = mean Jacobian (area per DOF).
         wjac_sum = 0.0
         !$acc loop seq
         do iq_local = 1, npts_l
            Iq = tsp%index_df_elt(iq_local, ie)
            wjac_sum = wjac_sum + tsp%wjac_df(Iq)
         end do
         wjac_sum = Csm2 * wjac_sum / real(npts_l)

         ! ν_s at each DOF: ν_s = (C_s·Δ)² · √(2 S_ij S_ij)
         !$acc loop seq
         do k = 1, nlayers_l
            !$acc loop seq
            do iq_local = 1, npts_l
               Iq     = tsp%index_df_elt(iq_local, ie)
               dp_inv = 1.0 / max(bcl%dpprime_visc(Iq,k), eps1)
               g11 = (bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(1,Iq) + bcl%dpp_graduvw(1,Iq,k)) * dp_inv
               g12 = (bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(2,Iq) + bcl%dpp_graduvw(2,Iq,k)) * dp_inv
               g21 = (bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(3,Iq) + bcl%dpp_graduvw(3,Iq,k)) * dp_inv
               g22 = (bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(4,Iq) + bcl%dpp_graduvw(4,Iq,k)) * dp_inv
               ! 2 S_ij S_ij (2D terms)
               two_S2 = 2.0*(g11**2 + g22**2) + (g12+g21)**2
               if (has_w) then
                  ! Additional 3D terms for sphere (indices 5-9: du_dz,dv_dz,dw_dx,dw_dy,dw_dz)
                  g13 = (bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(5,Iq) + bcl%dpp_graduvw(5,Iq,k)) * dp_inv
                  g23 = (bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(6,Iq) + bcl%dpp_graduvw(6,Iq,k)) * dp_inv
                  g31 = (bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(7,Iq) + bcl%dpp_graduvw(7,Iq,k)) * dp_inv
                  g32 = (bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(8,Iq) + bcl%dpp_graduvw(8,Iq,k)) * dp_inv
                  g33 = (bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(9,Iq) + bcl%dpp_graduvw(9,Iq,k)) * dp_inv
                  two_S2 = two_S2 + 2.0*g33**2 + (g13+g31)**2 + (g23+g32)**2
               end if
               nu_loc(iq_local,k) = wjac_sum * sqrt(two_S2)
            end do
         end do

         ! Scatter to global DOF array (elements are disjoint in the DG layout).
         Iq0 = tsp%index_df_elt(1, ie)
         !$acc loop seq
         do k = 1, nlayers_l
            !$acc loop seq
            do ip = 1, npts_l
               I = tsp%index_df(ip, Iq0)
               bcl%nu_smag(I,k) = nu_loc(ip,k)
            end do
         end do

      end do
      !$acc end parallel loop

      !$acc end data

   end subroutine bcl_compute_nu_smag

   subroutine btp_compute_nu_smag(G, inp, btp, tsp, graduv)

      use mod_grid,      only: grid
      use mod_input,     only: input
      use mod_variables, only: btp_CS
      use mod_tensor,    only: tensor_CS

      implicit none

      type(grid),      intent(in)    :: G
      type(input),     intent(in)    :: inp
      type(btp_CS),    intent(inout) :: btp
      type(tensor_CS), intent(in)    :: tsp

      real, intent(in) :: graduv(inp%ngraduvw_var,G%npoin)

      integer :: Iq
      real    :: g11, g12, g21, g22, g13, g23, g31, g32, g33
      real    :: two_S2, Csm2
      logical :: has_w

      has_w = (inp%nvar_btp == 5)
      Csm2  = inp%C_smag**2

      !$acc parallel loop present(graduv, btp%nu_smag, tsp%wjac_df)           &
      !$acc   private(g11, g12, g21, g22, g13, g23, g31, g32, g33, two_S2)   &
      !$acc   firstprivate(has_w, Csm2)
      do Iq = 1, G%npoin
         g11 = graduv(1,Iq); g12 = graduv(2,Iq)
         g21 = graduv(3,Iq); g22 = graduv(4,Iq)
         two_S2 = 2.0*(g11**2 + g22**2) + (g12+g21)**2
         if (has_w) then
            g13 = graduv(5,Iq); g23 = graduv(6,Iq)
            g31 = graduv(7,Iq); g32 = graduv(8,Iq); g33 = graduv(9,Iq)
            two_S2 = two_S2 + 2.0*g33**2 + (g13+g31)**2 + (g23+g32)**2
         end if
         btp%nu_smag(Iq) = Csm2 * tsp%wjac_df(Iq) * sqrt(two_S2)
      end do
      !$acc end parallel loop

   end subroutine btp_compute_nu_smag

   subroutine bcl_compute_laplacian(G, inp, b, btp, bcl, tsp, lap_q)

      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_input,     only: input
      use mod_variables, only: btp_CS, bcl_CS
      use mod_tensor,    only: tensor_CS

      implicit none

      type(grid),      intent(in)  :: G
      type(basis),     intent(in)  :: b
      type(input),     intent(in)  :: inp
      type(btp_CS),    intent(in)  :: btp
      type(bcl_CS),    intent(in)  :: bcl
      type(tensor_CS), intent(in)  :: tsp

      real, intent(out) :: lap_q(inp%nvar_bcl-1,G%npoin,inp%nlayers)

      integer :: ie, iq_local, Iq, Iq0, I, ip, k
      integer :: npts_l, nelem_l, nlayers_l
      real :: wq, qq(9), dhdx, dhdy, dhdz
      real :: lap_loc(3, b%npts, inp%nlayers)
      logical :: has_w

      npts_l    = b%npts
      nelem_l   = G%nelem
      nlayers_l = inp%nlayers
      has_w     = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

      !$acc data present(tsp%wjac_df, tsp%dpsidx_df, tsp%dpsidy_df,        &
      !$acc              tsp%dpsidz_df, tsp%dpsidz_df_x, tsp%dpsidz_df_y, tsp%dpsidz_df_z, &
      !$acc              tsp%index_df_elt, tsp%index_df,                     &
      !$acc              btp%graduvb_ave, bcl%dpprime_visc_scaled, bcl%dpp_graduvw_scaled, lap_q)

      !$acc kernels present(lap_q)
      lap_q = 0.0
      !$acc end kernels

      !$acc parallel loop gang                                                    &
      !$acc   private(lap_loc, qq, wq, dhdx, dhdy, dhdz, Iq, Iq0, I, ip, k, iq_local) &
      !$acc   firstprivate(npts_l, nelem_l, nlayers_l, has_w)
      do ie = 1, nelem_l

         !$acc loop seq
         do k = 1, nlayers_l
            !$acc loop seq
            do ip = 1, npts_l
               lap_loc(1,ip,k) = 0.0
               lap_loc(2,ip,k) = 0.0
               if (has_w) lap_loc(3,ip,k) = 0.0
            end do
         end do

         !$acc loop seq
         do k = 1, nlayers_l
            !$acc loop seq
            do iq_local = 1, npts_l
               Iq = tsp%index_df_elt(iq_local, ie)
               wq = tsp%wjac_df(Iq)
               qq(1) = bcl%dpprime_visc_scaled(Iq,k)*btp%graduvb_ave(1,Iq) + bcl%dpp_graduvw_scaled(1,Iq,k)
               qq(2) = bcl%dpprime_visc_scaled(Iq,k)*btp%graduvb_ave(2,Iq) + bcl%dpp_graduvw_scaled(2,Iq,k)
               qq(3) = bcl%dpprime_visc_scaled(Iq,k)*btp%graduvb_ave(3,Iq) + bcl%dpp_graduvw_scaled(3,Iq,k)
               qq(4) = bcl%dpprime_visc_scaled(Iq,k)*btp%graduvb_ave(4,Iq) + bcl%dpp_graduvw_scaled(4,Iq,k)
               if (has_w) then
                  qq(5) = bcl%dpprime_visc_scaled(Iq,k)*btp%graduvb_ave(5,Iq) + bcl%dpp_graduvw_scaled(5,Iq,k)
                  qq(6) = bcl%dpprime_visc_scaled(Iq,k)*btp%graduvb_ave(6,Iq) + bcl%dpp_graduvw_scaled(6,Iq,k)
                  qq(7) = bcl%dpprime_visc_scaled(Iq,k)*btp%graduvb_ave(7,Iq) + bcl%dpp_graduvw_scaled(7,Iq,k)
                  qq(8) = bcl%dpprime_visc_scaled(Iq,k)*btp%graduvb_ave(8,Iq) + bcl%dpp_graduvw_scaled(8,Iq,k)
                  qq(9) = bcl%dpprime_visc_scaled(Iq,k)*btp%graduvb_ave(9,Iq) + bcl%dpp_graduvw_scaled(9,Iq,k)
               end if
               !$acc loop seq
               do ip = 1, npts_l
                  dhdx = tsp%dpsidx_df(ip,Iq)
                  dhdy = tsp%dpsidy_df(ip,Iq)
                  dhdz = 0.0
                  if (has_w) then
                     dhdx = dhdx + tsp%dpsidz_df_x(ip,Iq)
                     dhdy = dhdy + tsp%dpsidz_df_y(ip,Iq)
                     dhdz = tsp%dpsidz_df(ip,Iq) + tsp%dpsidz_df_z(ip,Iq)
                  end if
                  lap_loc(1,ip,k) = lap_loc(1,ip,k) - wq*(dhdx*qq(1) + dhdy*qq(2))
                  lap_loc(2,ip,k) = lap_loc(2,ip,k) - wq*(dhdx*qq(3) + dhdy*qq(4))
                  if (has_w) then
                     lap_loc(1,ip,k) = lap_loc(1,ip,k) - wq*dhdz*qq(5)
                     lap_loc(2,ip,k) = lap_loc(2,ip,k) - wq*dhdz*qq(6)
                     lap_loc(3,ip,k) = lap_loc(3,ip,k) - wq*(dhdx*qq(7) + dhdy*qq(8) + dhdz*qq(9))
                  end if
               end do
            end do
         end do

         ! Elements are disjoint — scatter without atomics.
         Iq0 = tsp%index_df_elt(1, ie)
         !$acc loop seq
         do k = 1, nlayers_l
            !$acc loop seq
            do ip = 1, npts_l
               I = tsp%index_df(ip, Iq0)
               lap_q(1,I,k) = lap_loc(1,ip,k)
               lap_q(2,I,k) = lap_loc(2,ip,k)
               if (has_w) lap_q(3,I,k) = lap_loc(3,ip,k)
            end do
         end do

      end do
      !$acc end parallel loop

      !$acc end data

   end subroutine bcl_compute_laplacian

   subroutine bcl_compute_laplacian_qp(G, inp, b, btp, bcl, tsp, lap_q)

      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_input,     only: input
      use mod_variables, only: btp_CS, bcl_CS
      use mod_tensor,    only: tensor_CS

      implicit none

      type(grid),      intent(in)  :: G
      type(basis),     intent(in)  :: b
      type(input),     intent(in)  :: inp
      type(btp_CS),    intent(in)  :: btp
      type(bcl_CS),    intent(in)  :: bcl
      type(tensor_CS), intent(in)  :: tsp

      real, intent(out) :: lap_q(2,G%npoin,inp%nlayers)

      integer :: Iq, I, ip, k
      real :: wq, qq(4)

      lap_q = 0.0

      do k = 1, inp%nlayers
         do Iq = 1, G%npoin

            wq = tsp%wjac_df(Iq)

            qq(1) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(1,Iq) + bcl%dpp_graduvw(1,Iq,k)
            qq(2) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(2,Iq) + bcl%dpp_graduvw(2,Iq,k)
            qq(3) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(3,Iq) + bcl%dpp_graduvw(3,Iq,k)
            qq(4) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(4,Iq) + bcl%dpp_graduvw(4,Iq,k)

            do ip = 1, b%npts
               I = tsp%index_df(ip,Iq)
               lap_q(1,I,k) = lap_q(1,I,k) - wq*(tsp%dpsidx_df(ip,Iq)*qq(1) + tsp%dpsidy_df(ip,Iq)*qq(2))
               lap_q(2,I,k) = lap_q(2,I,k) - wq*(tsp%dpsidx_df(ip,Iq)*qq(3) + tsp%dpsidy_df(ip,Iq)*qq(4))
            end do
         end do
      end do

   end subroutine bcl_compute_laplacian_qp

   subroutine create_rhs_laplacian_flux(G, inp, b, mf, btp, tsp, Uk, rhs, gradq)
      !=========================================================================
      !  Face flux contribution to the barotropic viscous RHS.
      !
      !  GPU parallelism: one gang per face.  The quadrature loop (iquad) and
      !  the basis-node scatter loops (i) run sequentially within each gang,
      !  which avoids intra-gang atomics on nodes shared by adjacent iquad pts.
      !
      !  Race conditions:
      !    rhs(m,ip): multiple faces share boundary nodes — !$acc atomic update.
      !    btp%graduvb_face_ave(ivar,side,iquad,iface): each (iquad,iface) pair
      !      is touched by exactly one gang×seq-step — no atomics needed.
      !
      !  Sphere (has_w): gradq/btp_dpp_graduvw carry 9 components
      !  (du_dx,du_dy,dv_dx,dv_dy,du_dz,dv_dz,dw_dx,dw_dy,dw_dz); the wall
      !  reflection (ier==-4) mirrors each velocity's full (dx,dy,dz) gradient
      !  vector across the 3-component face normal, and w gets its own flux
      !  row (rhs component 3) built from (dw_dx,dw_dy,dw_dz).
      !=========================================================================
      use mod_grid,      only: grid
      use mod_input,     only: input
      use mod_basis,     only: basis
      use mod_face,      only: face_CS
      use mod_variables, only: btp_CS
      use mod_tensor,    only: tensor_CS

      implicit none

      type(grid),      intent(in)    :: G
      type(input),     intent(in)    :: inp
      type(basis),     intent(in)    :: b
      type(face_CS),   intent(in)    :: mf
      type(btp_CS),    intent(inout) :: btp
      type(tensor_CS), intent(in)    :: tsp

      real, intent(in)    :: Uk(inp%nvar_btp-2, G%npoin)
      real, intent(inout) :: rhs(inp%nvar_btp-2,G%npoin)
      real, intent(in)    :: gradq(inp%ngraduvw_var,G%npoin)

      real :: qu_mean(3), qv_mean(3), qw_mean(3)
      real :: flux_uv_visc_face(9,2)
      real :: qul(3), qur(3), qvl(3), qvr(3), qwl(3), qwr(3)
      real :: nx, ny, nz, wq, un
      integer :: iface, i, il, jl, kl, ir, jr, kr
      integer :: iel, ier, ip, iquad, ivar, nw
      real :: flux_qu, flux_qv, flux_qw, hi, alpha, beta, iflux
      real :: pconst, sigma, pb_avg, du_pen, dv_pen, dw_pen
      real    :: ql_s(9), qr_s(9), btp_ql_s(10), btp_qr_s(10)
      integer, dimension(b%ngl) :: I_l, I_r
      integer :: ngl_f, nface_f
      logical :: has_w

      beta  = 0.5
      alpha = 1.0 - beta
      iflux = 0.0

      ngl_f   = b%ngl
      nface_f = G%nface
      nw      = inp%ngraduvw_var
      has_w   = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex
      pconst  = real((ngl_f+1)*(ngl_f+2)) / 2.0 * real(inp%SIPG_constant)

      !$acc data present(G%face, G%face_type, G%intma,                       &
      !$acc               mf%imapl, mf%imapr, mf%normal_vector, mf%jac_face, &
      !$acc               b%psi,                                               &
      !$acc               btp%btp_dpp_graduvw_scaled, btp%pbprime_visc_scaled, &
      !$acc               btp%graduvb_face_ave,                                &
      !$acc               tsp%wjac_df, Uk,                                     &
      !$acc               rhs, gradq)

      !$acc parallel loop gang                                                &
      !$acc   private(iel, ier, il, jl, kl, ir, jr, kr, ip,                  &
      !$acc           I_l, I_r,                                               &
      !$acc           nx, ny, nz, un, wq,                                     &
      !$acc           ql_s, qr_s, btp_ql_s, btp_qr_s,                       &
      !$acc           flux_uv_visc_face,                                      &
      !$acc           qul, qur, qvl, qvr, qwl, qwr,                          &
      !$acc           qu_mean, qv_mean, qw_mean,                              &
      !$acc           sigma, pb_avg, du_pen, dv_pen, dw_pen,                  &
      !$acc           flux_qu, flux_qv, flux_qw, hi, ivar, i, iquad)          &
      !$acc   firstprivate(alpha, beta, iflux, ngl_f, nface_f, nw, has_w, pconst)
      do iface = 1, nface_f

         if (G%face_type(iface) == 2) cycle

         iel = G%face(7, iface)
         ier = G%face(8, iface)

         ! Precompute left global node indices once per face.
         !$acc loop seq
         do i = 1, ngl_f
            il = mf%imapl(1,i,1,iface)
            jl = mf%imapl(2,i,1,iface)
            kl = mf%imapl(3,i,1,iface)
            I_l(i) = G%intma(il,jl,kl,iel)
         end do

         if (ier > 0) then
            !$acc loop seq
            do i = 1, ngl_f
               ir = mf%imapr(1,i,1,iface)
               jr = mf%imapr(2,i,1,iface)
               kr = mf%imapr(3,i,1,iface)
               I_r(i) = G%intma(ir,jr,kr,ier)
            end do
         end if

         !$acc loop seq
         do iquad = 1, ngl_f

            nx = mf%normal_vector(1,iquad,1,iface)
            ny = mf%normal_vector(2,iquad,1,iface)
            nz = mf%normal_vector(3,iquad,1,iface)

            ip = I_l(iquad)

            !$acc loop seq
            do ivar = 1, nw
               ql_s(ivar)     = gradq(ivar,ip)
               btp_ql_s(ivar) = btp%btp_dpp_graduvw_scaled(ivar,ip)
            end do
            btp_ql_s(nw+1) = btp%pbprime_visc_scaled(ip)

            if (ier > 0) then

               ip = I_r(iquad)

               !$acc loop seq
               do ivar = 1, nw
                  qr_s(ivar)     = gradq(ivar,ip)
                  btp_qr_s(ivar) = btp%btp_dpp_graduvw_scaled(ivar,ip)
               end do
               btp_qr_s(nw+1) = btp%pbprime_visc_scaled(ip)

            else
               !$acc loop seq
               do ivar = 1, nw
                  qr_s(ivar)     = ql_s(ivar)
                  btp_qr_s(ivar) = btp_ql_s(ivar)
               end do
               btp_qr_s(nw+1) = btp_ql_s(nw+1)

               if (ier == -4) then

                  ! Mirror each velocity's full gradient vector across the
                  ! face normal. u: (dx,dy[,dz]) = indices (1,2[,5]);
                  ! v: (3,4[,6]); w (sphere only): (7,8,9).
                  if (has_w) then
                     un = ql_s(1)*nx + ql_s(2)*ny + ql_s(5)*nz
                     qr_s(1) = ql_s(1) - 2.0*un*nx
                     qr_s(2) = ql_s(2) - 2.0*un*ny
                     qr_s(5) = ql_s(5) - 2.0*un*nz

                     un = ql_s(3)*nx + ql_s(4)*ny + ql_s(6)*nz
                     qr_s(3) = ql_s(3) - 2.0*un*nx
                     qr_s(4) = ql_s(4) - 2.0*un*ny
                     qr_s(6) = ql_s(6) - 2.0*un*nz

                     un = ql_s(7)*nx + ql_s(8)*ny + ql_s(9)*nz
                     qr_s(7) = ql_s(7) - 2.0*un*nx
                     qr_s(8) = ql_s(8) - 2.0*un*ny
                     qr_s(9) = ql_s(9) - 2.0*un*nz

                     un = btp_ql_s(1)*nx + btp_ql_s(2)*ny + btp_ql_s(5)*nz
                     btp_qr_s(1) = btp_ql_s(1) - 2.0*un*nx
                     btp_qr_s(2) = btp_ql_s(2) - 2.0*un*ny
                     btp_qr_s(5) = btp_ql_s(5) - 2.0*un*nz

                     un = btp_ql_s(3)*nx + btp_ql_s(4)*ny + btp_ql_s(6)*nz
                     btp_qr_s(3) = btp_ql_s(3) - 2.0*un*nx
                     btp_qr_s(4) = btp_ql_s(4) - 2.0*un*ny
                     btp_qr_s(6) = btp_ql_s(6) - 2.0*un*nz

                     un = btp_ql_s(7)*nx + btp_ql_s(8)*ny + btp_ql_s(9)*nz
                     btp_qr_s(7) = btp_ql_s(7) - 2.0*un*nx
                     btp_qr_s(8) = btp_ql_s(8) - 2.0*un*ny
                     btp_qr_s(9) = btp_ql_s(9) - 2.0*un*nz
                  else
                     un = ql_s(1)*nx + ql_s(2)*ny
                     qr_s(1) = ql_s(1) - 2.0*un*nx
                     qr_s(2) = ql_s(2) - 2.0*un*ny

                     un = ql_s(3)*nx + ql_s(4)*ny
                     qr_s(3) = ql_s(3) - 2.0*un*nx
                     qr_s(4) = ql_s(4) - 2.0*un*ny

                     un = btp_ql_s(1)*nx + btp_ql_s(2)*ny
                     btp_qr_s(1) = btp_ql_s(1) - 2.0*un*nx
                     btp_qr_s(2) = btp_ql_s(2) - 2.0*un*ny

                     un = btp_ql_s(3)*nx + btp_ql_s(4)*ny
                     btp_qr_s(3) = btp_ql_s(3) - 2.0*un*nx
                     btp_qr_s(4) = btp_ql_s(4) - 2.0*un*ny
                  end if

               end if
            end if

            ! (iquad,iface) unique per gang×seq-step — no atomics needed.
            !$acc loop seq
            do ivar = 1, nw
               btp%graduvb_face_ave(ivar,1,iquad,iface) = btp%graduvb_face_ave(ivar,1,iquad,iface) + ql_s(ivar)
               btp%graduvb_face_ave(ivar,2,iquad,iface) = btp%graduvb_face_ave(ivar,2,iquad,iface) + qr_s(ivar)
            end do

            !$acc loop seq
            do ivar = 1, nw
               flux_uv_visc_face(ivar,1) = btp_ql_s(nw+1)*ql_s(ivar) + btp_ql_s(ivar)
               flux_uv_visc_face(ivar,2) = btp_qr_s(nw+1)*qr_s(ivar) + btp_qr_s(ivar)
            end do

            qul(1) = flux_uv_visc_face(1,1)
            qul(2) = flux_uv_visc_face(2,1)
            qvl(1) = flux_uv_visc_face(3,1)
            qvl(2) = flux_uv_visc_face(4,1)

            qur(1) = flux_uv_visc_face(1,2)
            qur(2) = flux_uv_visc_face(2,2)
            qvr(1) = flux_uv_visc_face(3,2)
            qvr(2) = flux_uv_visc_face(4,2)

            qul(3) = 0.0; qvl(3) = 0.0; qur(3) = 0.0; qvr(3) = 0.0
            qwl = 0.0; qwr = 0.0
            if (has_w) then
               qul(3) = flux_uv_visc_face(5,1); qvl(3) = flux_uv_visc_face(6,1)
               qur(3) = flux_uv_visc_face(5,2); qvr(3) = flux_uv_visc_face(6,2)
               qwl(1) = flux_uv_visc_face(7,1); qwl(2) = flux_uv_visc_face(8,1); qwl(3) = flux_uv_visc_face(9,1)
               qwr(1) = flux_uv_visc_face(7,2); qwr(2) = flux_uv_visc_face(8,2); qwr(3) = flux_uv_visc_face(9,2)
            end if

            qu_mean(1) = alpha*qul(1) + beta*qur(1)
            qu_mean(2) = alpha*qul(2) + beta*qur(2)
            qu_mean(3) = alpha*qul(3) + beta*qur(3)
            qv_mean(1) = alpha*qvl(1) + beta*qvr(1)
            qv_mean(2) = alpha*qvl(2) + beta*qvr(2)
            qv_mean(3) = alpha*qvl(3) + beta*qvr(3)
            qw_mean(1) = alpha*qwl(1) + beta*qwr(1)
            qw_mean(2) = alpha*qwl(2) + beta*qwr(2)
            qw_mean(3) = alpha*qwl(3) + beta*qwr(3)

            wq = mf%jac_face(iquad,1,iface)

            flux_qu = (qu_mean(1) - iflux*qul(1))*nx + (qu_mean(2) - iflux*qul(2))*ny &
                    + (qu_mean(3) - iflux*qul(3))*nz
            flux_qv = (qv_mean(1) - iflux*qvl(1))*nx + (qv_mean(2) - iflux*qvl(2))*ny &
                    + (qv_mean(3) - iflux*qvl(3))*nz
            flux_qw = (qw_mean(1) - iflux*qwl(1))*nx + (qw_mean(2) - iflux*qwl(2))*ny &
                    + (qw_mean(3) - iflux*qwl(3))*nz

            ! SIP interior penalty on interior faces only (ier > 0).
            ! Jump [ū] = Uk_L - Uk_R; average p̄_b_avg = 0.5*(pbprime_L + pbprime_R).
            if (ier > 0) then
               sigma  = pconst * wq / min(tsp%wjac_df(I_l(iquad)), tsp%wjac_df(I_r(iquad)))
               pb_avg = 0.5 * (btp%pbprime_visc_scaled(I_l(iquad)) + btp%pbprime_visc_scaled(I_r(iquad)))
               du_pen = Uk(1, I_l(iquad)) - Uk(1, I_r(iquad))
               dv_pen = Uk(2, I_l(iquad)) - Uk(2, I_r(iquad))
               flux_qu = flux_qu - sigma * pb_avg * du_pen
               flux_qv = flux_qv - sigma * pb_avg * dv_pen
               if (has_w) then
                  dw_pen  = Uk(3, I_l(iquad)) - Uk(3, I_r(iquad))
                  flux_qw = flux_qw - sigma * pb_avg * dw_pen
               end if
            end if

            ! rhs scatter — atomics: multiple faces share boundary nodes.
            !$acc loop seq
            do i = 1, ngl_f
               hi = b%psi(i,iquad)
               ip = I_l(i)
               !$acc atomic update
               rhs(1,ip) = rhs(1,ip) + wq*hi*flux_qu
               !$acc atomic update
               rhs(2,ip) = rhs(2,ip) + wq*hi*flux_qv
               if (has_w) then
                  !$acc atomic update
                  rhs(3,ip) = rhs(3,ip) + wq*hi*flux_qw
               end if
            end do

            if (ier > 0) then
               !$acc loop seq
               do i = 1, ngl_f
                  hi = b%psi(i,iquad)
                  ip = I_r(i)
                  !$acc atomic update
                  rhs(1,ip) = rhs(1,ip) - wq*hi*flux_qu
                  !$acc atomic update
                  rhs(2,ip) = rhs(2,ip) - wq*hi*flux_qv
                  if (has_w) then
                     !$acc atomic update
                     rhs(3,ip) = rhs(3,ip) - wq*hi*flux_qw
                  end if
               end do
            end if

         end do  ! iquad

      end do  ! iface
      !$acc end data

   end subroutine create_rhs_laplacian_flux

   subroutine bcl_create_rhs_laplacian_flux(G, inp, b, mf, btp, bcl, tsp, rhs)

      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_input,     only: input
      use mod_face,      only: face_CS
      use mod_variables, only: btp_CS, bcl_CS
      use mod_tensor,    only: tensor_CS

      implicit none

      type(grid),      intent(in)    :: G
      type(basis),     intent(in)    :: b
      type(input),     intent(in)    :: inp
      type(face_CS),   intent(in)    :: mf
      type(btp_CS),    intent(in)    :: btp
      type(bcl_CS),    intent(in)    :: bcl
      type(tensor_CS), intent(in)    :: tsp

      real, intent(inout) :: rhs(inp%nvar_bcl-1,G%npoin,inp%nlayers)

      real, dimension(3) :: qu_mean, qv_mean, qw_mean
      real, dimension(9,2) :: flux_uv_visc_face
      real, dimension(3) :: qul, qur, qvl, qvr, qwl, qwr

      real :: nx, ny, nz, un
      real :: wq
      integer :: iface, i, il, jl, kl, ir, jr, kr
      integer :: iel, ier, ip
      integer :: iquad, ivar, k, nw
      real :: flux_qu, flux_qv, flux_qw, hi, alpha, beta, iflux
      real :: ql_cur(10), qr_cur(10)
      integer :: ngl_f, nlayers_f, nface_f
      logical :: has_w
      integer :: ip_face_L, ip_face_R

      beta  = 0.5
      alpha = 1.0 - beta
      iflux = 0.0

      ngl_f     = b%ngl
      nlayers_f = inp%nlayers
      nface_f   = G%nface
      nw        = inp%ngraduvw_var
      has_w     = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

      !$acc data present(G%face, G%face_type, G%intma,                          &
      !$acc              mf%imapl, mf%imapr, mf%normal_vector, mf%jac_face,     &
      !$acc              b%psi,                                                   &
      !$acc              btp%graduvb_face_ave,                                   &
      !$acc              bcl%dpp_graduvw_scaled, bcl%dpprime_visc_scaled, bcl%qprime_df, &
      !$acc              tsp%wjac_df,                                             &
      !$acc              rhs)

      !$acc parallel loop gang                                                    &
      !$acc   private(ql_cur, qr_cur, flux_uv_visc_face, qul, qur, qvl, qvr, qwl, qwr, &
      !$acc           qu_mean, qv_mean, qw_mean,                                     &
      !$acc           nx, ny, nz, un, wq, flux_qu, flux_qv, flux_qw, hi,            &
      !$acc           iel, ier, ip, il, jl, kl, ir, jr, kr, iquad, k, i, ivar,  &
      !$acc           ip_face_L, ip_face_R) &
      !$acc   firstprivate(alpha, beta, iflux, ngl_f, nlayers_f, nface_f, nw, has_w)
      do iface = 1, nface_f

         if (G%face_type(iface) == 2) cycle

         iel = G%face(7,iface)
         ier = G%face(8,iface)

         !$acc loop seq
         do k = 1, nlayers_f
            !$acc loop seq
            do iquad = 1, ngl_f

               nx = mf%normal_vector(1,iquad,1,iface)
               ny = mf%normal_vector(2,iquad,1,iface)
               nz = mf%normal_vector(3,iquad,1,iface)

               il = mf%imapl(1,iquad,1,iface)
               jl = mf%imapl(2,iquad,1,iface)
               kl = mf%imapl(3,iquad,1,iface)
               ip = G%intma(il,jl,kl,iel)
               ip_face_L = ip
               ip_face_R = ip   ! default = same node (zero jump for boundary faces)

               !$acc loop seq
               do ivar = 1, nw
                  ql_cur(ivar) = bcl%dpp_graduvw_scaled(ivar,ip,k)
               end do
               ql_cur(nw+1) = bcl%dpprime_visc_scaled(ip,k)

               if (ier > 0) then
                  ir = mf%imapr(1,iquad,1,iface)
                  jr = mf%imapr(2,iquad,1,iface)
                  kr = mf%imapr(3,iquad,1,iface)
                  ip = G%intma(ir,jr,kr,ier)
                  ip_face_R = ip
                  !$acc loop seq
                  do ivar = 1, nw
                     qr_cur(ivar) = bcl%dpp_graduvw_scaled(ivar,ip,k)
                  end do
                  qr_cur(nw+1) = bcl%dpprime_visc_scaled(ip,k)
               else
                  !$acc loop seq
                  do ivar = 1, nw+1
                     qr_cur(ivar) = ql_cur(ivar)
                  end do
                  if (ier == -4) then
                     if (has_w) then
                        un = ql_cur(1)*nx + ql_cur(2)*ny + ql_cur(5)*nz
                        qr_cur(1) = ql_cur(1) - 2.0*un*nx
                        qr_cur(2) = ql_cur(2) - 2.0*un*ny
                        qr_cur(5) = ql_cur(5) - 2.0*un*nz

                        un = ql_cur(3)*nx + ql_cur(4)*ny + ql_cur(6)*nz
                        qr_cur(3) = ql_cur(3) - 2.0*un*nx
                        qr_cur(4) = ql_cur(4) - 2.0*un*ny
                        qr_cur(6) = ql_cur(6) - 2.0*un*nz

                        un = ql_cur(7)*nx + ql_cur(8)*ny + ql_cur(9)*nz
                        qr_cur(7) = ql_cur(7) - 2.0*un*nx
                        qr_cur(8) = ql_cur(8) - 2.0*un*ny
                        qr_cur(9) = ql_cur(9) - 2.0*un*nz
                     else
                        un = ql_cur(1)*nx + ql_cur(2)*ny
                        qr_cur(1) = ql_cur(1) - 2.0*un*nx
                        qr_cur(2) = ql_cur(2) - 2.0*un*ny
                        un = ql_cur(3)*nx + ql_cur(4)*ny
                        qr_cur(3) = ql_cur(3) - 2.0*un*nx
                        qr_cur(4) = ql_cur(4) - 2.0*un*ny
                     end if
                  end if
               end if

               !$acc loop seq
               do ivar = 1, nw
                  flux_uv_visc_face(ivar,1) = ql_cur(nw+1)*btp%graduvb_face_ave(ivar,1,iquad,iface) + ql_cur(ivar)
                  flux_uv_visc_face(ivar,2) = qr_cur(nw+1)*btp%graduvb_face_ave(ivar,2,iquad,iface) + qr_cur(ivar)
               end do

               qul(1) = flux_uv_visc_face(1,1);  qul(2) = flux_uv_visc_face(2,1)
               qvl(1) = flux_uv_visc_face(3,1);  qvl(2) = flux_uv_visc_face(4,1)
               qur(1) = flux_uv_visc_face(1,2);  qur(2) = flux_uv_visc_face(2,2)
               qvr(1) = flux_uv_visc_face(3,2);  qvr(2) = flux_uv_visc_face(4,2)

               qul(3) = 0.0; qvl(3) = 0.0; qur(3) = 0.0; qvr(3) = 0.0
               qwl = 0.0; qwr = 0.0
               if (has_w) then
                  qul(3) = flux_uv_visc_face(5,1); qvl(3) = flux_uv_visc_face(6,1)
                  qur(3) = flux_uv_visc_face(5,2); qvr(3) = flux_uv_visc_face(6,2)
                  qwl(1) = flux_uv_visc_face(7,1); qwl(2) = flux_uv_visc_face(8,1); qwl(3) = flux_uv_visc_face(9,1)
                  qwr(1) = flux_uv_visc_face(7,2); qwr(2) = flux_uv_visc_face(8,2); qwr(3) = flux_uv_visc_face(9,2)
               end if

               qu_mean(1) = alpha*qul(1) + beta*qur(1)
               qu_mean(2) = alpha*qul(2) + beta*qur(2)
               qu_mean(3) = alpha*qul(3) + beta*qur(3)
               qv_mean(1) = alpha*qvl(1) + beta*qvr(1)
               qv_mean(2) = alpha*qvl(2) + beta*qvr(2)
               qv_mean(3) = alpha*qvl(3) + beta*qvr(3)
               qw_mean(1) = alpha*qwl(1) + beta*qwr(1)
               qw_mean(2) = alpha*qwl(2) + beta*qwr(2)
               qw_mean(3) = alpha*qwl(3) + beta*qwr(3)

               wq = mf%jac_face(iquad,1,iface)
               flux_qu = (qu_mean(1) - iflux*qul(1))*nx + (qu_mean(2) - iflux*qul(2))*ny &
                       + (qu_mean(3) - iflux*qul(3))*nz
               flux_qv = (qv_mean(1) - iflux*qvl(1))*nx + (qv_mean(2) - iflux*qvl(2))*ny &
                       + (qv_mean(3) - iflux*qvl(3))*nz
               flux_qw = (qw_mean(1) - iflux*qwl(1))*nx + (qw_mean(2) - iflux*qwl(2))*ny &
                       + (qw_mean(3) - iflux*qwl(3))*nz

               !$acc loop seq
               do i = 1, ngl_f
                  hi = b%psi(i,iquad)
                  il = mf%imapl(1,i,1,iface)
                  jl = mf%imapl(2,i,1,iface)
                  kl = mf%imapl(3,i,1,iface)
                  ip = G%intma(il,jl,kl,iel)
                  !$acc atomic update
                  rhs(1,ip,k) = rhs(1,ip,k) + wq*hi*flux_qu
                  !$acc atomic update
                  rhs(2,ip,k) = rhs(2,ip,k) + wq*hi*flux_qv
                  if (has_w) then
                     !$acc atomic update
                     rhs(3,ip,k) = rhs(3,ip,k) + wq*hi*flux_qw
                  end if
               end do

               if (ier > 0) then
                  !$acc loop seq
                  do i = 1, ngl_f
                     hi = b%psi(i,iquad)
                     ir = mf%imapr(1,i,1,iface)
                     jr = mf%imapr(2,i,1,iface)
                     kr = mf%imapr(3,i,1,iface)
                     ip = G%intma(ir,jr,kr,ier)
                     !$acc atomic update
                     rhs(1,ip,k) = rhs(1,ip,k) - wq*hi*flux_qu
                     !$acc atomic update
                     rhs(2,ip,k) = rhs(2,ip,k) - wq*hi*flux_qv
                     if (has_w) then
                        !$acc atomic update
                        rhs(3,ip,k) = rhs(3,ip,k) - wq*hi*flux_qw
                     end if
                  end do
               end if

            end do !iquad
         end do !k
      end do !iface
      !$acc end parallel loop

      !$acc end data

   end subroutine bcl_create_rhs_laplacian_flux

   ! ===========================================================================
   !  compute_bcl_lap_z
   !
   !  LDG (Local Discontinuous Galerkin) Laplacian of the layer-interface
   !  height z_k for internal interfaces k = 2..nlayers.
   !
   !  Three-phase CPU implementation:
   !
   !  Phase 1 — Gradient: element-local ∇z_k stored globally in grad_z.
   !    grad_z(d, J, k) = Σ_ip dpsid_df(ip, J) × z_k(J_ip)   d=1,2[,3]
   !
   !  [Pre-communicate grad_z at MPI faces — TODO: add MPI exchange here
   !   using the same pack/send/recv pattern as bcl_lap_create_precommunicator
   !   so that MPI boundary faces also get the correct central-flux gradient.]
   !
   !  Phase 2 — Central-flux correction: for every interior face, replace
   !    grad_z at the shared boundary DOF nodes with the element average
   !    {∇z_k} = ½(∇z_k^L + ∇z_k^R), giving the correct LDG flux.
   !    MPI faces (face_type == 2) are skipped — they await the pre-comm TODO.
   !
   !  Phase 3 — Volume integral: element-local divergence of the corrected
   !    gradient gives the LDG Laplacian:
   !    lap_z(J_ip0) = Σ_ip dpsid_df(ip, J_ip0) × grad_z(d, J_ip, k)
   !
   !  CPU-only (no OpenACC): call !$acc update device(bcl%lap_z_df) afterwards.
   ! ===========================================================================
   ! ===========================================================================
   !  bcl_create_rhs_gradz_flux
   !  Interior-face LDG central-flux correction for the lap_z_df accumulator.
   !  Adds {grad_z}*n face term to lap_z_df at interior face DOFs (both sides).
   ! ===========================================================================
   subroutine bcl_create_rhs_gradz_flux(G, inp, b, mf, bcl, lap_z_df)

      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_input,     only: input
      use mod_face,      only: face_CS
      use mod_variables, only: bcl_CS

      implicit none

      type(grid),      intent(in)    :: G
      type(basis),     intent(in)    :: b
      type(input),     intent(in)    :: inp
      type(face_CS),   intent(in)    :: mf
      type(bcl_CS),    intent(in)    :: bcl

      real, intent(inout) :: lap_z_df(G%npoin, inp%nlayers+1)

      integer :: iface, iel, ier, k, iquad, i, ip, il, jl, kl, ir, jr, kr
      integer :: ip_face_L, ip_face_R
      real    :: wq, nx, ny, nz, flux_gradz, hi
      real    :: gzl_x, gzl_y, gzl_z, gzr_x, gzr_y, gzr_z
      logical :: has_w

      has_w = (inp%nvar_bcl == 4)

      do iface = 1, G%nface
         if (G%face_type(iface) == 2) cycle
         iel = G%face(7, iface)
         ier = G%face(8, iface)

         do k = 2, inp%nlayers
            do iquad = 1, b%ngl

               il = mf%imapl(1,iquad,1,iface);  jl = mf%imapl(2,iquad,1,iface);  kl = mf%imapl(3,iquad,1,iface)
               ip_face_L = G%intma(il, jl, kl, iel)

               nx = mf%normal_vector(1, iquad, 1, iface)
               ny = mf%normal_vector(2, iquad, 1, iface)
               wq = mf%jac_face(iquad, 1, iface)

               gzl_x = bcl%grad_z_df(1, ip_face_L, k)
               gzl_y = bcl%grad_z_df(2, ip_face_L, k)
               if (has_w) then
                  gzl_z = bcl%grad_z_df(3, ip_face_L, k)
                  nz    = mf%normal_vector(3, iquad, 1, iface)
               end if

               if (ier > 0) then
                  ir = mf%imapr(1,iquad,1,iface);  jr = mf%imapr(2,iquad,1,iface);  kr = mf%imapr(3,iquad,1,iface)
                  ip_face_R = G%intma(ir, jr, kr, ier)
                  gzr_x = bcl%grad_z_df(1, ip_face_R, k)
                  gzr_y = bcl%grad_z_df(2, ip_face_R, k)
                  if (has_w) gzr_z = bcl%grad_z_df(3, ip_face_R, k)
               else
                  ! Neumann: mirror gradient so {grad_z}*n = gzl*n at boundary
                  gzr_x = gzl_x
                  gzr_y = gzl_y
                  if (has_w) then
                     gzr_z = gzl_z
                     if (ier == -4) then
                        ! Slip wall: zero normal gradient
                        nz = mf%normal_vector(3, iquad, 1, iface)
                        gzr_x = gzl_x - 2.0*(gzl_x*nx + gzl_y*ny + gzl_z*nz)*nx
                        gzr_y = gzl_y - 2.0*(gzl_x*nx + gzl_y*ny + gzl_z*nz)*ny
                        gzr_z = gzl_z - 2.0*(gzl_x*nx + gzl_y*ny + gzl_z*nz)*nz
                     end if
                  else
                     if (ier == -4) then
                        ! Slip wall: zero normal gradient
                        gzr_x = gzl_x - 2.0*(gzl_x*nx + gzl_y*ny)*nx
                        gzr_y = gzl_y - 2.0*(gzl_x*nx + gzl_y*ny)*ny
                     end if
                  end if
               end if

               flux_gradz = 0.5*(gzl_x + gzr_x)*nx + 0.5*(gzl_y + gzr_y)*ny
               if (has_w) flux_gradz = flux_gradz + 0.5*(gzl_z + gzr_z)*nz

               do i = 1, b%ngl
                  hi = b%psi(i, iquad)
                  il = mf%imapl(1,i,1,iface);  jl = mf%imapl(2,i,1,iface);  kl = mf%imapl(3,i,1,iface)
                  ip = G%intma(il, jl, kl, iel)
                  lap_z_df(ip, k) = lap_z_df(ip, k) + wq*hi*flux_gradz
               end do

               if (ier > 0) then
                  do i = 1, b%ngl
                     hi = b%psi(i, iquad)
                     ir = mf%imapr(1,i,1,iface);  jr = mf%imapr(2,i,1,iface);  kr = mf%imapr(3,i,1,iface)
                     ip = G%intma(ir, jr, kr, ier)
                     lap_z_df(ip, k) = lap_z_df(ip, k) - wq*hi*flux_gradz
                  end do
               end if

            end do   ! iquad
         end do   ! k
      end do   ! iface

   end subroutine bcl_create_rhs_gradz_flux

   ! ===========================================================================
   !  compute_bcl_lap_z  —  full LDG Laplacian of interface height z_k
   !
   !  Mirrors the bcl_create_laplacian 4-call pattern:
   !    Phase 1 (gradient)   → bcl%grad_z_df
   !    pre-comm             → bcl_gradz_create_precommunicator
   !    Phase 3 (weak vol.)  → lap_z_df (mass-weighted integral)
   !    interior flux        → bcl_create_rhs_gradz_flux
   !    post-comm            → bcl_create_rhs_gradz_postcommunicator_df
   !    mass inverse         → lap_z_df = massinv * lap_z_df
   !
   !  CPU-only: caller must !$acc update host before and
   !            !$acc update device(bcl%lap_z_df) after.
   ! ===========================================================================
   subroutine compute_bcl_lap_z(G, inp, b, mf, par, ref, mpic, init, btp, bcl, mt, tsp, qprime_df, lap_z_df)

      use mod_grid,             only: grid
      use mod_input,            only: input
      use mod_basis,            only: basis
      use mod_face,             only: face_CS
      use mod_parallel,         only: parallel_CS
      use mod_ref,              only: mref
      use mod_mpi_communicator, only: mpi_communicator
      use mod_initial,          only: initial
      use mod_variables,        only: btp_CS, bcl_CS
      use mod_metrics,          only: metrics
      use mod_tensor,           only: tensor_CS
      use mod_constants,        only: gravity

      implicit none

      type(grid),             intent(in)    :: G
      type(input),            intent(in)    :: inp
      type(basis),            intent(in)    :: b
      type(face_CS),          intent(in)    :: mf
      type(parallel_CS),      intent(in)    :: par
      type(mref),             intent(inout) :: ref
      type(mpi_communicator), intent(inout) :: mpic
      type(initial),          intent(in)    :: init
      type(btp_CS),           intent(in)    :: btp
      type(bcl_CS),           intent(inout) :: bcl
      type(metrics),          intent(in)    :: mt
      type(tensor_CS),        intent(in)    :: tsp
      real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: qprime_df
      real, dimension(G%npoin, inp%nlayers+1),              intent(out) :: lap_z_df

      integer :: ie, iq_local, ip, k, I
      integer :: Iq, Iq0
      integer :: npts_l, nelem_l, nlayers_l
      real    :: z_loc(b%npts, inp%nlayers+1)
      real    :: lap_loc(b%npts, inp%nlayers)
      real    :: dhdx, dhdy, dhdz, wq, gzx, gzy, gzz
      logical :: has_w

      has_w     = (inp%nvar_bcl == 4)
      npts_l    = b%npts
      nelem_l   = G%nelem
      nlayers_l = inp%nlayers

      lap_z_df(:,:)        = 0.0
      bcl%grad_z_df(:,:,:) = 0.0

      ! -----------------------------------------------------------------------
      ! Phase 1: element-local gradient of z_k -> bcl%grad_z_df.
      ! Outer index ip = test function / evaluation point (nodal collocation).
      ! Inner index iq_local = basis function summed over.
      ! -----------------------------------------------------------------------
      do ie = 1, nelem_l

         do ip = 1, npts_l
            Iq = tsp%index_df_elt(ip, ie)
            z_loc(ip, nlayers_l+1) = init%zbot_df(Iq)
            do k = nlayers_l, 1, -1
               z_loc(ip, k) = z_loc(ip, k+1) + &
                  (init%alpha_mlswe(k)/gravity) * sqrt(btp%ope2_ave_df(Iq)) * qprime_df(1, Iq, k)
            end do
         end do

         do k = 2, nlayers_l
            do ip = 1, npts_l
               Iq = tsp%index_df_elt(ip, ie)
               bcl%grad_z_df(1, Iq, k) = 0.0
               bcl%grad_z_df(2, Iq, k) = 0.0
               if (has_w) bcl%grad_z_df(3, Iq, k) = 0.0
               do iq_local = 1, npts_l
                  dhdx = tsp%dpsidx_df(iq_local, Iq)
                  dhdy = tsp%dpsidy_df(iq_local, Iq)
                  if (has_w) then
                     dhdx = dhdx + tsp%dpsidz_df_x(iq_local, Iq)
                     dhdy = dhdy + tsp%dpsidz_df_y(iq_local, Iq)
                     dhdz = tsp%dpsidz_df(iq_local, Iq) + tsp%dpsidz_df_z(iq_local, Iq)
                  end if
                  bcl%grad_z_df(1, Iq, k) = bcl%grad_z_df(1, Iq, k) + dhdx * z_loc(iq_local, k)
                  bcl%grad_z_df(2, Iq, k) = bcl%grad_z_df(2, Iq, k) + dhdy * z_loc(iq_local, k)
                  if (has_w) bcl%grad_z_df(3, Iq, k) = bcl%grad_z_df(3, Iq, k) + dhdz * z_loc(iq_local, k)
               end do
            end do
         end do

      end do  ! ie (Phase 1)

      ! Send grad_z_df at MPI boundary face DOFs while Phase 3 runs.
      call bcl_gradz_create_precommunicator(G, inp, b, mf, par, ref, mpic, bcl%grad_z_df)

      ! -----------------------------------------------------------------------
      ! Phase 3: weak-form volume integral  -wjac * dpsi/dx * grad_z_df.
      ! Mirrors bcl_compute_laplacian: outer loop = quadrature points (Iq),
      ! inner loop = test functions (ip).
      ! -----------------------------------------------------------------------
      do ie = 1, nelem_l

         do k = 2, nlayers_l
            do ip = 1, npts_l
               lap_loc(ip, k) = 0.0
            end do
         end do

         do k = 2, nlayers_l
            do iq_local = 1, npts_l
               Iq  = tsp%index_df_elt(iq_local, ie)
               wq  = tsp%wjac_df(Iq)
               gzx = bcl%grad_z_df(1, Iq, k)
               gzy = bcl%grad_z_df(2, Iq, k)
               if (has_w) gzz = bcl%grad_z_df(3, Iq, k)
               do ip = 1, npts_l
                  dhdx = tsp%dpsidx_df(ip, Iq)
                  dhdy = tsp%dpsidy_df(ip, Iq)
                  lap_loc(ip, k) = lap_loc(ip, k) - wq*(dhdx*gzx + dhdy*gzy)
                  if (has_w) then
                     dhdz = tsp%dpsidz_df(ip, Iq) + tsp%dpsidz_df_z(ip, Iq)
                     lap_loc(ip, k) = lap_loc(ip, k) - wq*dhdz*gzz
                  end if
               end do
            end do
         end do

         Iq0 = tsp%index_df_elt(1, ie)
         do k = 2, nlayers_l
            do ip = 1, npts_l
               I = tsp%index_df(ip, Iq0)
               lap_z_df(I, k) = lap_loc(ip, k)
            end do
         end do

      end do  ! ie (Phase 3)

      ! Interior face flux correction: {grad_z}*n face term.
      call bcl_create_rhs_gradz_flux(G, inp, b, mf, bcl, lap_z_df)

      ! Wait for MPI recv; apply same correction at MPI boundary faces.
      call bcl_create_rhs_gradz_postcommunicator_df(G, inp, b, mf, par, btp, ref, mpic, tsp, lap_z_df)

      ! Mass inverse: convert weak-form integral to pointwise Laplacian value.
      do k = 2, nlayers_l
         lap_z_df(:, k) = mt%massinv(:) * lap_z_df(:, k)
      end do

   end subroutine compute_bcl_lap_z

end module mod_laplacian_quad
