! ===============================================================================================
! This module contains the barotropic right-hand side by GLL collocation
! (btp_collocation = .true., default .false.).
!
!   The barotropic fluxes and sources are evaluated at the N+1 GLL nodes of each element
!   and face instead of the 2N+1 over-integration points. The volume, face and MPI-face
!   kernels below follow create_rhs_btp_volume_qdf_qp, create_btp_fluxes_qdf and
!   create_nbhs_face_df with the nodal tensor arrays, normals and fields; the volume
!   derivative terms are projected onto the test functions by sum factorisation.
!
!   The time averages the baroclinic step needs are accumulated at the nodes with the
!   effective SSPRK stage weights (init%ssprk_w) and interpolated once per stage to the
!   quadrature-point arrays (btp_colloc_to_quad), so the baroclinic code is unchanged.
!   The layer sum is the reference for the column mass: the barotropic column is driven
!   toward it with a mass source, as MOM6's eta_cor (btp_colloc_set_etacor).
!
!   The nodal and face-node arrays live in btp_CS (mod_variables).
!   Date: November 20, 2026
! ===============================================================================================

module mod_btp_colloc

   implicit none

   private

   public :: btp_colloc_check, btp_colloc_begin, create_rhs_btp_colloc, &
             btp_colloc_to_quad, btp_colloc_set_etacor, btp_colloc_measure, &
             btp_colloc_report

   ! Largest |p_b - sum_k dp_k| / p_b since the last report (btp_colloc_diag)
   real :: fix_max = 0.0

   ! 1 / (number of barotropic sub-stages in the current cycle): sums -> averages
   real :: colloc_ninv = 1.0

   ! Stage weight kstages*b_i of the current stage (init%ssprk_w), set by the driver
   real :: wst = 1.0

contains

   subroutine btp_colloc_check(inp)

      ! Stop if btp_collocation is combined with options it does not support yet

      use mod_input,         only: input
      use mod_mpi_utilities, only: irank, irank0, exit_all
      use mpi

      implicit none

      type(input), intent(in) :: inp

      character(len=100) :: msg
      integer            :: ierr

      if (.not. inp%btp_collocation) return

      msg = ''
      if (inp%lgpu_kernels .and. inp%c_APE > 0.0) msg = 'APE term not available with lgpu_kernels'
      if (trim(inp%bcl_time_method) /= 'lsrk3' .and. trim(inp%bcl_time_method) /= 'rk3' .and. &
          trim(inp%bcl_time_method) /= '2levels') &
         msg = "needs bcl_time_method = 'lsrk3', 'rk3' or '2levels'"

      ! MPI_ABORT stops every rank (error stop on its own can leave ranks waiting).
      ! All ranks see the same input, so they all reach the barrier; it lets rank 0
      ! print the message before the abort.
      if (len_trim(msg) > 0) then
         if (irank == irank0) then
            print *, 'ERROR btp_collocation: ', trim(msg)
            flush(6)
         end if
         call mpi_barrier(mpi_comm_world, ierr)
         call exit_all(1)
      end if

   end subroutine btp_colloc_check

   subroutine btp_colloc_begin(G, inp, b, init, btp, qprime_df)

      ! Start of a barotropic sub-cycle: zero the nodal accumulators and build the
      ! nodal layer coefficients from qprime_df (fixed during the cycle). Same as the
      ! quad-point loop of btp_bcl_coeffs_qdf, at the nodes.

      use mod_grid,      only: grid
      use mod_input,     only: input
      use mod_basis,     only: basis
      use mod_initial,   only: initial
      use mod_constants, only: gravity
      use mod_variables, only: btp_CS

      implicit none

      type(grid),    intent(in)    :: G
      type(input),   intent(in)    :: inp
      type(basis),   intent(in)    :: b
      type(initial), intent(in)    :: init
      type(btp_CS),  intent(inout) :: btp

      real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers), intent(in) :: qprime_df

      integer :: I, k, nv, npoin_l, nlayers_l
      real    :: pp_k, up_k, vp_k, wp_k, pprime_k, pprime_k1, H_b, flux(3,3), g_l, dry_l
      logical :: has_w

      nv        = inp%nvar_btp - 2
      has_w     = (inp%nvar_bcl == 4)
      npoin_l   = G%npoin
      nlayers_l = inp%nlayers
      g_l       = gravity
      dry_l     = inp%dry_cutoff

      ! Accumulators of this sub-cycle (ope2_ave_df and uvb_ave_df: mod_rk_mlswe)
      !$acc kernels async(1) present(btp%ope_ave_df, btp%Qu_ave_df, btp%Qv_ave_df, btp%Qw_ave_df,    &
      !$acc                  btp%tau_bot_ave_df, btp%btp_mass_flux_ave_df,                         &
      !$acc                  btp%btp_mass_flux_face_ave_df, btp%one_plus_eta_edge_2_ave_df,        &
      !$acc                  btp%Qu_face_ave_df, btp%Qv_face_ave_df, btp%Qw_face_ave_df,           &
      !$acc                  btp%uvb_face_ave_df, btp%ope_face_ave_df, btp%ope2_face_ave_df)
      btp%ope_ave_df             = 0.0
      btp%Qu_ave_df              = 0.0;  btp%Qv_ave_df                  = 0.0
      btp%Qw_ave_df              = 0.0;  btp%tau_bot_ave_df             = 0.0
      btp%btp_mass_flux_ave_df   = 0.0;  btp%btp_mass_flux_face_ave_df  = 0.0
      btp%one_plus_eta_edge_2_ave_df = 0.0
      btp%Qu_face_ave_df         = 0.0;  btp%Qv_face_ave_df             = 0.0
      btp%Qw_face_ave_df         = 0.0;  btp%uvb_face_ave_df            = 0.0
      btp%ope_face_ave_df        = 0.0;  btp%ope2_face_ave_df           = 0.0
      !$acc end kernels

      colloc_ninv = 1.0 / real(inp%kstages * init%N_btp)

      ! Split-form pressure weights beta_m/2 (interface m = k+1):
      ! beta_m = alpha_{m-1} - alpha_m, beta_{L+1} = alpha_L, so H_b = sum_m (beta_m/2) A_m^2
      do k = 1, inp%nlayers - 1
         btp%bcl_dalpha(k) = 0.5 * (init%alpha_mlswe(k) - init%alpha_mlswe(k+1))
      end do
      btp%bcl_dalpha(inp%nlayers) = 0.5 * init%alpha_mlswe(inp%nlayers)
      !$acc update device(btp%bcl_dalpha) async(1)

      ! Layer coefficients at the nodes
      !$acc parallel loop gang vector async(1)                                              &
      !$acc    present(qprime_df, init%alpha_mlswe, btp%bcl_H_df, btp%bcl_flux_df,          &
      !$acc            btp%bcl_btp_flux_df, btp%bcl_pprime_df)                              &
      !$acc    private(k, pp_k, up_k, vp_k, wp_k, pprime_k, pprime_k1, H_b, flux)           &
      !$acc    firstprivate(nlayers_l, has_w, nv, g_l, dry_l)
      do I = 1, npoin_l

         H_b      = 0.0
         pprime_k = 0.0
         flux     = 0.0
         pp_k     = 0.0
         up_k     = 0.0
         vp_k     = 0.0
         wp_k     = 0.0

         !$acc loop seq
         do k = 1, nlayers_l

            pp_k = qprime_df(1,I,k)
            up_k = qprime_df(2,I,k)
            vp_k = qprime_df(3,I,k)
            wp_k = 0.0
            if (has_w) wp_k = qprime_df(4,I,k)

            ! Dry layers: floor the thickness, no velocity (as btp_bcl_coeffs_qdf)
            if (pp_k < (g_l/init%alpha_mlswe(k)) * dry_l) then
               pp_k = (g_l/init%alpha_mlswe(k)) * dry_l
               up_k = 0.0
               vp_k = 0.0
               wp_k = 0.0
            end if

            pprime_k1 = pprime_k + pp_k
            H_b = H_b + 0.5 * init%alpha_mlswe(k) * (pprime_k1**2 - pprime_k**2)
            btp%bcl_pprime_df(k,I) = pprime_k1

            flux(1,1) = flux(1,1) + up_k*(pp_k*up_k)
            flux(2,1) = flux(2,1) + vp_k*(pp_k*up_k)
            flux(1,2) = flux(1,2) + up_k*(pp_k*vp_k)
            flux(2,2) = flux(2,2) + vp_k*(pp_k*vp_k)
            if (has_w) then
               flux(3,1) = flux(3,1) + wp_k*(pp_k*up_k)
               flux(3,2) = flux(3,2) + wp_k*(pp_k*vp_k)
               flux(1,3) = flux(1,3) + up_k*(pp_k*wp_k)
               flux(2,3) = flux(2,3) + vp_k*(pp_k*wp_k)
               flux(3,3) = flux(3,3) + wp_k*(pp_k*wp_k)
            end if

            pprime_k = pprime_k1

         end do

         ! After the layer loop pp_k, up_k, vp_k, wp_k hold the bottom-layer values
         btp%bcl_H_df(I)          = H_b
         btp%bcl_flux_df(:,:,I)   = flux(1:nv,1:nv)
         btp%bcl_btp_flux_df(1,I) = pp_k
         btp%bcl_btp_flux_df(2,I) = up_k
         btp%bcl_btp_flux_df(3,I) = vp_k
         if (has_w) btp%bcl_btp_flux_df(4,I) = wp_k

      end do
      !$acc end parallel loop

   end subroutine btp_colloc_begin

   subroutine create_rhs_btp_colloc(G, inp, b, mf, par, btp, init, ref, mpic, tsp, mt, ik, &
                                    rhs_btp, qb_df, qprime_df)

      ! Barotropic RHS by collocation (replaces create_rhs_btp for stage ik)

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
      use mod_metrics,          only: metrics

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
      type(metrics),          intent(in)    :: mt
      integer,                intent(in)    :: ik

      real, dimension(inp%nvar_btp-1,G%npoin),           intent(out) :: rhs_btp
      real, dimension(inp%nvar_btp,G%npoin),             intent(in)  :: qb_df
      real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers), intent(in)  :: qprime_df

      integer :: I, npoin_l

      ! Weight of this stage in the time averages
      wst = real(inp%kstages) * init%ssprk_w(ik)

      ! 1. MPI halo exchange
      call btp_create_precommunicator(G, inp, b, mf, init, par, ref, mpic, qb_df, qprime_df, inp%nvar_btp)

      ! 2. Volume terms (zeros rhs_btp), then 3. local face fluxes
      if (inp%lgpu_kernels) then
         call btp_colloc_volume_gpu(G, b, inp, init, btp, tsp, mt, rhs_btp, qb_df)
      else
         call btp_colloc_volume_cpu(G, b, inp, init, btp, tsp, mt, rhs_btp, qb_df)
      end if
      call btp_colloc_faces(G, b, inp, mf, init, btp, rhs_btp, qb_df, qprime_df)

      ! 4. Receive the neighbour states
      call mpi_waitall(mpic%nreq, mpic%ireq, mpic%status, mpic%ierr)
      call unpack_data_dg_general_df(G, b, par, ref%q_send, ref%q_recv, ref%send_data_dg, ref%recv_data_dg, &
                                     ref%nbtp_var, ref%nboun_valid)

      ! 5. MPI face fluxes
      call btp_colloc_nbhs_faces(G, inp, b, mf, par, btp, init, ref, rhs_btp)

      ! 6. eta_cor: mass source toward the layer sum (btp_colloc_set_etacor)
      npoin_l = G%npoin
      !$acc parallel loop gang vector async(1) present(rhs_btp, btp%etacor_src_df)
      do I = 1, npoin_l
         rhs_btp(1,I) = rhs_btp(1,I) + btp%etacor_src_df(I)
      end do
      !$acc end parallel loop

   end subroutine create_rhs_btp_colloc

   subroutine btp_colloc_to_quad(G, inp, b, tsp, btp)

      ! Interpolate the nodal accumulators to the quadrature-point arrays read by the
      ! baroclinic step (before ti_barotropic_ssprk_mlswe normalises them).

      use mod_grid,      only: grid
      use mod_input,     only: input
      use mod_basis,     only: basis
      use mod_tensor,    only: tensor_CS
      use mod_variables, only: btp_CS

      implicit none

      type(grid),      intent(in)    :: G
      type(input),     intent(in)    :: inp
      type(basis),     intent(in)    :: b
      type(tensor_CS), intent(in)    :: tsp
      type(btp_CS),    intent(inout) :: btp

      if (inp%lgpu_kernels) then
         call btp_colloc_to_quad_gpu(G, inp, b, btp)
      else
         call btp_colloc_to_quad_cpu(G, inp, b, tsp, btp)
      end if

   end subroutine btp_colloc_to_quad

   subroutine btp_colloc_to_quad_cpu(G, inp, b, tsp, btp)

      ! Interpolation to the quadrature points (CPU).
      ! Volume: sum factorisation per element, x then y. Faces: along each face.

      use mod_grid,      only: grid
      use mod_input,     only: input
      use mod_basis,     only: basis
      use mod_tensor,    only: tensor_CS
      use mod_variables, only: btp_CS

      implicit none

      type(grid),      intent(in)    :: G
      type(input),     intent(in)    :: inp
      type(basis),     intent(in)    :: b
      type(tensor_CS), intent(in)    :: tsp
      type(btp_CS),    intent(inout) :: btp

      ! Element work arrays. Scalar fields: 1 (1+eta), 2 (1+eta)^2.
      ! Vector fields (nvar_btp-2 components): 1 Qu, 2 Qv, 3 Qw, 4 tau_bot,
      ! 5 mass flux, 6 velocity.
      real :: Ws(2, b%nglx, b%ngly), Wv(inp%nvar_btp-2, 6, b%nglx, b%ngly)
      real :: Ts(2, b%nqx,  b%ngly), Tv(inp%nvar_btp-2, 6, b%nqx,  b%ngly)
      real :: Qs(2, b%nqx,  b%nqy),  Qv(inp%nvar_btp-2, 6, b%nqx,  b%nqy)

      ! Interface pressures A_m: nodes, interpolated in x, at the quadrature points
      real :: Vp(inp%nlayers, b%nglx, b%ngly), Tp(inp%nlayers, b%nqx, b%ngly), Qp(inp%nlayers, b%nqx, b%nqy)
      real :: pprime_ql(inp%nlayers), pprime_qr(inp%nlayers)

      integer :: Iq, I, iface, iquad, n, m, ie, kq, lq
      real    :: hi

      do ie = 1, G%nelem

         ! Nodal values of the element
         do m = 1, b%ngly
            do n = 1, b%nglx
               I = G%intma(n, m, 1, ie)
               Ws(1, n, m)   = btp%ope_ave_df(I)
               Ws(2, n, m)   = btp%ope2_ave_df(I)
               Wv(:, 1, n, m) = btp%Qu_ave_df(:,I)
               Wv(:, 2, n, m) = btp%Qv_ave_df(:,I)
               Wv(:, 3, n, m) = btp%Qw_ave_df(:,I)
               Wv(:, 4, n, m) = btp%tau_bot_ave_df(:,I)
               Wv(:, 5, n, m) = btp%btp_mass_flux_ave_df(:,I)
               Wv(:, 6, n, m) = btp%uvb_ave_df(:,I)
               Vp(:, n, m)    = btp%bcl_pprime_df(:,I)
            end do
         end do

         ! Interpolate in x
         Ts = 0.0
         Tv = 0.0
         Tp = 0.0
         do m = 1, b%ngly
            do kq = 1, b%nqx
               do n = 1, b%nglx
                  Ts(:, kq, m)    = Ts(:, kq, m)    + b%psiqx(n, kq) * Ws(:, n, m)
                  Tv(:, :, kq, m) = Tv(:, :, kq, m) + b%psiqx(n, kq) * Wv(:, :, n, m)
                  Tp(:, kq, m)    = Tp(:, kq, m)    + b%psiqx(n, kq) * Vp(:, n, m)
               end do
            end do
         end do

         ! Interpolate in y
         Qs = 0.0
         Qv = 0.0
         Qp = 0.0
         do lq = 1, b%nqy
            do m = 1, b%ngly
               Qs(:, :, lq)    = Qs(:, :, lq)    + b%psiqy(m, lq) * Ts(:, :, m)
               Qv(:, :, :, lq) = Qv(:, :, :, lq) + b%psiqy(m, lq) * Tv(:, :, :, m)
               Qp(:, :, lq)    = Qp(:, :, lq)    + b%psiqy(m, lq) * Tp(:, :, m)
            end do
         end do

         ! Store at the quadrature points
         do lq = 1, b%nqy
            do kq = 1, b%nqx
               Iq = G%intma_dg_quad(kq, lq, 1, ie)
               ! H_b is fixed in the sub-cycle: H_ave = ope2_ave * H_b, with H_b from the
               ! interface pressures interpolated to this point (as the layers build theirs)
               btp%H_ave(Iq)               = Qs(2, kq, lq) * sum(btp%bcl_dalpha(:) * Qp(:, kq, lq)**2)
               btp%ope_ave(Iq)             = Qs(1, kq, lq)
               btp%ope2_ave(Iq)            = Qs(2, kq, lq)
               btp%Qu_ave(:,Iq)            = Qv(:, 1, kq, lq)
               btp%Qv_ave(:,Iq)            = Qv(:, 2, kq, lq)
               btp%Qw_ave(:,Iq)            = Qv(:, 3, kq, lq)
               btp%tau_bot_ave(:,Iq)       = Qv(:, 4, kq, lq)
               btp%btp_mass_flux_ave(:,Iq) = Qv(:, 5, kq, lq)
               btp%uvb_ave(:,Iq)           = Qv(:, 6, kq, lq)
            end do
         end do

      end do

      ! Faces: face nodes -> face quadrature points
      btp%btp_mass_flux_face_ave  = 0.0
      btp%H_face_ave              = 0.0
      btp%Qu_face_ave             = 0.0
      btp%Qv_face_ave             = 0.0
      btp%Qw_face_ave             = 0.0
      btp%ope_face_ave            = 0.0
      btp%ope2_face_ave           = 0.0
      btp%one_plus_eta_edge_2_ave = 0.0
      btp%uvb_face_ave            = 0.0

      do iface = 1, G%nface
         do iquad = 1, b%nq
            do n = 1, b%ngl

               hi = b%psiq(n, iquad)

               btp%btp_mass_flux_face_ave(:,iquad,iface) = btp%btp_mass_flux_face_ave(:,iquad,iface) &
                                                         + hi * btp%btp_mass_flux_face_ave_df(:,n,iface)
               btp%Qu_face_ave(:,iquad,iface)   = btp%Qu_face_ave(:,iquad,iface)   + hi * btp%Qu_face_ave_df(:,n,iface)
               btp%Qv_face_ave(:,iquad,iface)   = btp%Qv_face_ave(:,iquad,iface)   + hi * btp%Qv_face_ave_df(:,n,iface)
               btp%Qw_face_ave(:,iquad,iface)   = btp%Qw_face_ave(:,iquad,iface)   + hi * btp%Qw_face_ave_df(:,n,iface)
               btp%ope_face_ave(:,iquad,iface)  = btp%ope_face_ave(:,iquad,iface)  + hi * btp%ope_face_ave_df(:,n,iface)
               btp%ope2_face_ave(:,iquad,iface) = btp%ope2_face_ave(:,iquad,iface) + hi * btp%ope2_face_ave_df(:,n,iface)
               btp%one_plus_eta_edge_2_ave(iquad,iface) = btp%one_plus_eta_edge_2_ave(iquad,iface) &
                                                        + hi * btp%one_plus_eta_edge_2_ave_df(n,iface)
               btp%uvb_face_ave(:,:,iquad,iface) = btp%uvb_face_ave(:,:,iquad,iface) &
                                                 + hi * btp%uvb_face_ave_df(:,:,n,iface)

            end do
         end do
      end do

      ! Face H_face_ave from the interface pressures of each side, interpolated along
      ! the face, times the averaged (1+eta)^2 (as for the volume H_ave above)
      do iface = 1, G%nface
         do iquad = 1, b%nq

            pprime_ql = 0.0
            pprime_qr = 0.0
            do n = 1, b%ngl
               hi = b%psiq(n, iquad)
               pprime_ql = pprime_ql + hi * btp%bcl_pprime_face_df(:,1,n,iface)
               pprime_qr = pprime_qr + hi * btp%bcl_pprime_face_df(:,2,n,iface)
            end do

            btp%H_face_ave(iquad,iface) = btp%one_plus_eta_edge_2_ave(iquad,iface) &
                                        * 0.5 * sum(btp%bcl_dalpha(:) * (pprime_ql**2 + pprime_qr**2))

         end do
      end do

   end subroutine btp_colloc_to_quad_cpu

   subroutine btp_colloc_to_quad_gpu(G, inp, b, btp)

      ! Interpolate the nodal accumulators to the quadrature-point arrays read by the
      ! baroclinic step (GPU): one thread per volume quadrature point and per face
      ! quadrature point. H_ave and H_face_ave are built from the interface pressures
      ! interpolated to the point, as in the CPU btp_colloc_to_quad.

      use mod_grid,      only: grid
      use mod_input,     only: input
      use mod_basis,     only: basis
      use mod_variables, only: btp_CS

      implicit none

      type(grid),   intent(in)    :: G
      type(input),  intent(in)    :: inp
      type(basis),  intent(in)    :: b
      type(btp_CS), intent(inout) :: btp

      integer :: ie, kq, lq, n, m, iv, k, I, Iq, iface, iquad, nv
      integer :: nelem_l, nglx_l, ngly_l, nqx_l, nqy_l, nlayers_l, nface_l, nq_l, ngl_l
      real    :: hi, s_ope, s_ope2, s_H, p_k, pl, pr
      real    :: s_qu(3), s_qv(3), s_qw(3), s_tb(3), s_mf(3), s_uv(3)
      real    :: s_mff(3), s_quf(3), s_qvf(3), s_qwf(3), s_opef(2), s_ope2f(2), s_oe2, s_uvf(3,2)

      nv        = inp%nvar_btp - 2
      nelem_l   = G%nelem
      nglx_l    = b%nglx
      ngly_l    = b%ngly
      nqx_l     = b%nqx
      nqy_l     = b%nqy
      nlayers_l = inp%nlayers
      nface_l   = G%nface
      nq_l      = b%nq
      ngl_l     = b%ngl

      ! Volume quadrature points
      !$acc parallel loop gang vector collapse(3) async(1)                                     &
      !$acc    present(G%intma, G%intma_dg_quad, b%psiqx, b%psiqy, btp%ope_ave_df,             &
      !$acc            btp%ope2_ave_df, btp%Qu_ave_df, btp%Qv_ave_df, btp%Qw_ave_df,           &
      !$acc            btp%tau_bot_ave_df, btp%btp_mass_flux_ave_df, btp%uvb_ave_df,           &
      !$acc            btp%bcl_pprime_df, btp%bcl_dalpha,                                      &
      !$acc            btp%H_ave, btp%ope_ave, btp%ope2_ave, btp%Qu_ave, btp%Qv_ave,           &
      !$acc            btp%Qw_ave, btp%tau_bot_ave, btp%btp_mass_flux_ave, btp%uvb_ave)        &
      !$acc    private(Iq, I, n, m, iv, k, hi, s_ope, s_ope2, s_H, p_k,                        &
      !$acc            s_qu, s_qv, s_qw, s_tb, s_mf, s_uv)                                     &
      !$acc    firstprivate(nv, nglx_l, ngly_l, nlayers_l)
      do ie = 1, nelem_l
         do lq = 1, nqy_l
            do kq = 1, nqx_l

               s_ope  = 0.0
               s_ope2 = 0.0
               s_qu = 0.0;  s_qv = 0.0;  s_qw = 0.0;  s_tb = 0.0;  s_mf = 0.0;  s_uv = 0.0
               !$acc loop seq
               do m = 1, ngly_l
                  !$acc loop seq
                  do n = 1, nglx_l
                     I  = G%intma(n, m, 1, ie)
                     hi = b%psiqx(n, kq) * b%psiqy(m, lq)
                     s_ope  = s_ope  + hi * btp%ope_ave_df(I)
                     s_ope2 = s_ope2 + hi * btp%ope2_ave_df(I)
                     !$acc loop seq
                     do iv = 1, nv
                        s_qu(iv) = s_qu(iv) + hi * btp%Qu_ave_df(iv,I)
                        s_qv(iv) = s_qv(iv) + hi * btp%Qv_ave_df(iv,I)
                        s_qw(iv) = s_qw(iv) + hi * btp%Qw_ave_df(iv,I)
                        s_tb(iv) = s_tb(iv) + hi * btp%tau_bot_ave_df(iv,I)
                        s_mf(iv) = s_mf(iv) + hi * btp%btp_mass_flux_ave_df(iv,I)
                        s_uv(iv) = s_uv(iv) + hi * btp%uvb_ave_df(iv,I)
                     end do
                  end do
               end do

               ! H_b from the interface pressures interpolated to this point
               s_H = 0.0
               !$acc loop seq
               do k = 1, nlayers_l
                  p_k = 0.0
                  !$acc loop seq
                  do m = 1, ngly_l
                     !$acc loop seq
                     do n = 1, nglx_l
                        p_k = p_k + b%psiqx(n, kq) * b%psiqy(m, lq) * btp%bcl_pprime_df(k, G%intma(n, m, 1, ie))
                     end do
                  end do
                  s_H = s_H + btp%bcl_dalpha(k) * p_k * p_k
               end do

               Iq = G%intma_dg_quad(kq, lq, 1, ie)
               btp%H_ave(Iq)    = s_ope2 * s_H
               btp%ope_ave(Iq)  = s_ope
               btp%ope2_ave(Iq) = s_ope2
               !$acc loop seq
               do iv = 1, nv
                  btp%Qu_ave(iv,Iq)            = s_qu(iv)
                  btp%Qv_ave(iv,Iq)            = s_qv(iv)
                  btp%Qw_ave(iv,Iq)            = s_qw(iv)
                  btp%tau_bot_ave(iv,Iq)       = s_tb(iv)
                  btp%btp_mass_flux_ave(iv,Iq) = s_mf(iv)
                  btp%uvb_ave(iv,Iq)           = s_uv(iv)
               end do

            end do
         end do
      end do
      !$acc end parallel loop

      ! Face quadrature points
      !$acc parallel loop gang vector collapse(2) async(1)                                     &
      !$acc    present(b%psiq, btp%btp_mass_flux_face_ave_df, btp%Qu_face_ave_df,              &
      !$acc            btp%Qv_face_ave_df, btp%Qw_face_ave_df, btp%ope_face_ave_df,            &
      !$acc            btp%ope2_face_ave_df, btp%one_plus_eta_edge_2_ave_df,                   &
      !$acc            btp%uvb_face_ave_df, btp%bcl_pprime_face_df, btp%bcl_dalpha,            &
      !$acc            btp%btp_mass_flux_face_ave, btp%H_face_ave, btp%Qu_face_ave,            &
      !$acc            btp%Qv_face_ave, btp%Qw_face_ave, btp%ope_face_ave, btp%ope2_face_ave,  &
      !$acc            btp%one_plus_eta_edge_2_ave, btp%uvb_face_ave)                          &
      !$acc    private(n, iv, k, hi, s_mff, s_quf, s_qvf, s_qwf, s_opef, s_ope2f, s_oe2,       &
      !$acc            s_uvf, s_H, pl, pr)                                                     &
      !$acc    firstprivate(nv, ngl_l, nlayers_l)
      do iface = 1, nface_l
         do iquad = 1, nq_l

            s_mff = 0.0;  s_quf = 0.0;  s_qvf = 0.0;  s_qwf = 0.0
            s_opef = 0.0;  s_ope2f = 0.0;  s_oe2 = 0.0;  s_uvf = 0.0
            !$acc loop seq
            do n = 1, ngl_l
               hi = b%psiq(n, iquad)
               !$acc loop seq
               do iv = 1, nv
                  s_mff(iv)    = s_mff(iv)    + hi * btp%btp_mass_flux_face_ave_df(iv,n,iface)
                  s_quf(iv)    = s_quf(iv)    + hi * btp%Qu_face_ave_df(iv,n,iface)
                  s_qvf(iv)    = s_qvf(iv)    + hi * btp%Qv_face_ave_df(iv,n,iface)
                  s_qwf(iv)    = s_qwf(iv)    + hi * btp%Qw_face_ave_df(iv,n,iface)
                  s_uvf(iv,1)  = s_uvf(iv,1)  + hi * btp%uvb_face_ave_df(iv,1,n,iface)
                  s_uvf(iv,2)  = s_uvf(iv,2)  + hi * btp%uvb_face_ave_df(iv,2,n,iface)
               end do
               s_opef(1)  = s_opef(1)  + hi * btp%ope_face_ave_df(1,n,iface)
               s_opef(2)  = s_opef(2)  + hi * btp%ope_face_ave_df(2,n,iface)
               s_ope2f(1) = s_ope2f(1) + hi * btp%ope2_face_ave_df(1,n,iface)
               s_ope2f(2) = s_ope2f(2) + hi * btp%ope2_face_ave_df(2,n,iface)
               s_oe2      = s_oe2      + hi * btp%one_plus_eta_edge_2_ave_df(n,iface)
            end do

            ! H_face_ave from the interface pressures of each side at this point
            s_H = 0.0
            !$acc loop seq
            do k = 1, nlayers_l
               pl = 0.0
               pr = 0.0
               !$acc loop seq
               do n = 1, ngl_l
                  pl = pl + b%psiq(n, iquad) * btp%bcl_pprime_face_df(k,1,n,iface)
                  pr = pr + b%psiq(n, iquad) * btp%bcl_pprime_face_df(k,2,n,iface)
               end do
               s_H = s_H + btp%bcl_dalpha(k) * (pl*pl + pr*pr)
            end do

            !$acc loop seq
            do iv = 1, nv
               btp%btp_mass_flux_face_ave(iv,iquad,iface) = s_mff(iv)
               btp%Qu_face_ave(iv,iquad,iface)            = s_quf(iv)
               btp%Qv_face_ave(iv,iquad,iface)            = s_qvf(iv)
               btp%Qw_face_ave(iv,iquad,iface)            = s_qwf(iv)
               btp%uvb_face_ave(iv,1,iquad,iface)         = s_uvf(iv,1)
               btp%uvb_face_ave(iv,2,iquad,iface)         = s_uvf(iv,2)
            end do
            btp%ope_face_ave(1,iquad,iface)          = s_opef(1)
            btp%ope_face_ave(2,iquad,iface)          = s_opef(2)
            btp%ope2_face_ave(1,iquad,iface)         = s_ope2f(1)
            btp%ope2_face_ave(2,iquad,iface)         = s_ope2f(2)
            btp%one_plus_eta_edge_2_ave(iquad,iface) = s_oe2
            btp%H_face_ave(iquad,iface)              = s_oe2 * 0.5 * s_H

         end do
      end do
      !$acc end parallel loop

   end subroutine btp_colloc_to_quad_gpu

   subroutine btp_colloc_set_etacor(G, inp, mt, btp, q_df, qb_df)

      ! MOM6-style eta_cor (bt_mass_source): at the start of a long step the mismatch
      ! d = sum_k dp_k - p_b becomes a constant barotropic mass source d/dt for this
      ! step, so that p_b is driven toward the layer sum (the layers are the reference).
      ! Stored before M^-1, as the barotropic update multiplies rhs_btp by massinv.

      use mod_grid,      only: grid
      use mod_input,     only: input
      use mod_metrics,   only: metrics
      use mod_variables, only: btp_CS

      implicit none

      type(grid),    intent(in)    :: G
      type(input),   intent(in)    :: inp
      type(metrics), intent(in)    :: mt
      type(btp_CS),  intent(inout) :: btp

      real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers), intent(in) :: q_df
      real, dimension(inp%nvar_btp,G%npoin),             intent(in) :: qb_df

      integer :: I, npoin_l
      real    :: dt_l

      npoin_l = G%npoin
      dt_l    = inp%dt

      !$acc parallel loop gang vector present(q_df, qb_df, mt%massinv, btp%etacor_src_df) &
      !$acc    firstprivate(dt_l)
      do I = 1, npoin_l
         btp%etacor_src_df(I) = (sum(q_df(1,I,:)) - qb_df(1,I)) / (dt_l * mt%massinv(I))
      end do
      !$acc end parallel loop

   end subroutine btp_colloc_set_etacor

   subroutine btp_colloc_measure(G, inp, q_df, qb_df)

      ! Diagnostic (btp_colloc_diag): largest |p_b - sum_k dp_k| / p_b at the nodes

      use mod_grid,  only: grid
      use mod_input, only: input

      implicit none

      type(grid),  intent(in) :: G
      type(input), intent(in) :: inp

      real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers), intent(in) :: q_df
      real, dimension(inp%nvar_btp,G%npoin),             intent(in) :: qb_df

      integer :: I
      real    :: s

      do I = 1, G%npoin
         s = sum(q_df(1,I,:))
         if (qb_df(1,I) > 0.0) fix_max = max(fix_max, abs(qb_df(1,I) - s) / qb_df(1,I))
      end do

   end subroutine btp_colloc_measure

   subroutine btp_colloc_report(inp, itime)

      ! Print (rank 0) the largest |p_b - sum_k dp_k| / p_b since the last call

      use mpi
      use mod_input,         only: input
      use mod_mpi_utilities, only: irank, irank0, MPI_PRECISION

      implicit none

      type(input), intent(in) :: inp
      integer,     intent(in) :: itime

      real    :: gmax
      integer :: ierr

      if (.not. inp%btp_colloc_diag) return

      call mpi_reduce(fix_max, gmax, 1, MPI_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      if (irank == irank0) write(*,'(a,i8,a,es10.3)') ' btp_collocation: step ', itime, &
         '  max |p_b - sum dp_k|/p_b since last report = ', gmax

      fix_max = 0.0

   end subroutine btp_colloc_report

   subroutine btp_colloc_volume_cpu(G, b, inp, init, btp, tsp, mt, rhs_btp, qb_df)
      !===========================================================================
      !  Collocated barotropic volume term (CPU prototype). Same physics as
      !  create_rhs_btp_volume_qdf_qp, evaluated at the GLL nodes of each element;
      !  the weak-form derivative terms are projected onto the test functions by sum
      !  factorisation (row and column of the 1D derivative matrix) with the nodal
      !  metrics.
      !===========================================================================
      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_constants, only: gravity
      use mod_initial,   only: initial
      use mod_variables, only: btp_CS
      use mod_input,     only: input
      use mod_tensor,    only: tensor_CS
      use mod_initial_mlswe, only: check_btp_thickness
      use mod_ape,       only: ape_on, ape_force_btp_df
      use mod_metrics,   only: metrics

      implicit none

      type(grid),      intent(in)    :: G
      type(basis),     intent(in)    :: b
      type(input),     intent(in)    :: inp
      type(initial),   intent(in)    :: init
      type(btp_CS),    intent(inout) :: btp
      type(tensor_CS), intent(in)    :: tsp
      type(metrics),   intent(in)    :: mt

      real, dimension(inp%nvar_btp,G%npoin),   intent(in)  :: qb_df
      real, dimension(inp%nvar_btp-1,G%npoin), intent(out) :: rhs_btp

      real :: sc_x, sc_y, sc_z, Hq, qu, quv, qvu, qv, quz, qvz, qwu, qwv, qw, H_b
      real :: wq, hi, dhdx, dhdy, dhdz, tb_u, tb_v, tb_w, ope, ope2
      real :: dp, dpp, udp, vdp, wdp, ub, vb, wb, ubot, vbot, wbot, spd, pbq
      real :: pp_k, up_k, vp_k, wp_k
      real :: sum_up2, sum_uv, sum_vu, sum_vp2, sum_uw, sum_vw, sum_ww, sum_wu, sum_wv
      real :: cfx, cfy, cfz
      integer :: I, Iq, ie, kx, ky, l, neq
      real    :: e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z
      real    :: F(3,4), S(4), Fe(4), Fn(4), Fc(4)
      ! Split-form pressure: element node values (ope_pp_e = ope * A_m), metrics and the
      ! antisymmetric pair weights h_e(i,j) = w_j D_ji - w_i D_ij (h_n likewise)
      real    :: ope_pp_e(inp%nlayers, b%nglx, b%ngly), jac_e(3, b%nglx, b%ngly), jac_n(3, b%nglx, b%ngly)
      real    :: h_e(b%nglx, b%nglx), h_n(b%ngly, b%ngly)
      real    :: H_bcl_qij, f_ij, t_ij, jac
      integer :: nv, m, ic, Il
      integer :: npts_l, botfr_l
      real    :: cd_l, alpha_bot_l
      logical :: has_w, ape_l

      npts_l      = b%npts
      botfr_l     = inp%botfr
      cd_l        = real(inp%cd_mlswe)
      alpha_bot_l = real(init%alpha_mlswe(inp%nlayers))
      has_w       = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex
      ape_l       = ape_on        ! Chen APE: column sum of the layer forces (mod_ape)


      rhs_btp = 0.0

      ! Abort with element/rank diagnostics rather than SIGFPE-ing on ub=udp/dp
      ! below if the barotropic depth has gone negative or NaN (e.g. interface
      ! outcropping near steep topography feeding through into qb_df).
      call check_btp_thickness(b, G, inp, qb_df, 'btp_colloc_volume_cpu')

      neq = 3
      if (has_w) neq = 4
      nv  = neq - 1

      ! Pair weights of the split-form pressure (basis only): h(j,i) = -h(i,j)
      do kx = 1, b%nglx
         do l = 1, b%nglx
            h_e(kx,l) = b%wglx(l)*b%dpsix(kx,l) - b%wglx(kx)*b%dpsix(l,kx)
         end do
      end do
      do ky = 1, b%ngly
         do l = 1, b%ngly
            h_n(ky,l) = b%wgly(l)*b%dpsiy(ky,l) - b%wgly(ky)*b%dpsiy(l,ky)
         end do
      end do

      do ie = 1, G%nelem
         do ky = 1, b%ngly
            do kx = 1, b%nglx
               Iq = G%intma(kx, ky, 1, ie)

               dp = 0.0;  dpp = 0.0;  udp = 0.0;  vdp = 0.0;  wdp = 0.0

               ! Collocation: the barotropic state at the node itself.
               dp  = qb_df(1, Iq);  dpp = qb_df(2, Iq)
               udp = qb_df(3, Iq);  vdp = qb_df(4, Iq)
               if (has_w) wdp = qb_df(5, Iq)

               ! Load BCL layer integrals precomputed once per RK stage.
               H_b     = btp%bcl_H_df(Iq)
               sum_up2 = btp%bcl_flux_df(1,1,Iq)
               sum_uv  = btp%bcl_flux_df(1,2,Iq)
               sum_vu  = btp%bcl_flux_df(2,1,Iq)
               sum_vp2 = btp%bcl_flux_df(2,2,Iq)
               pp_k    = btp%bcl_btp_flux_df(1,Iq)
               up_k    = btp%bcl_btp_flux_df(2,Iq)
               vp_k    = btp%bcl_btp_flux_df(3,Iq)
               pbq     = init%pbprime_df(Iq)

               sum_uw = 0.0;  sum_vw = 0.0;  sum_ww = 0.0;  wp_k = 0.0
               sum_wu = 0.0;  sum_wv = 0.0
               if (has_w) then
                  sum_uw = btp%bcl_flux_df(3,1,Iq)
                  sum_vw = btp%bcl_flux_df(3,2,Iq)
                  sum_wu = btp%bcl_flux_df(1,3,Iq)
                  sum_wv = btp%bcl_flux_df(2,3,Iq)
                  sum_ww = btp%bcl_flux_df(3,3,Iq)
                  wp_k   = btp%bcl_btp_flux_df(4,Iq)
               end if

               wq = tsp%wjac_df(Iq)
               ub = udp / dp;  vb = vdp / dp;  wb = wdp / dp

               ! 3D Coriolis vector (f * r̂); reduces to the f-plane pair (0,0,f) on
               ! cartesian since kvector=(0,0,1) there — see wind_stress_coriolis.
               cfx = init%coriolis_df(Iq)*init%kvector(1, Iq)
               cfy = init%coriolis_df(Iq)*init%kvector(2, Iq)
               cfz = init%coriolis_df(Iq)*init%kvector(3, Iq)

               tb_u = 0.0;  tb_v = 0.0;  tb_w = 0.0
               if (botfr_l == 1) then
                  ubot = up_k + ub;  vbot = vp_k + vb;  wbot = wp_k + wb
                  spd  = (cd_l / gravity) * pp_k
                  tb_u = spd * ubot;  tb_v = spd * vbot;  tb_w = spd * wbot
               elseif (botfr_l == 2) then
                  ubot = up_k + ub;  vbot = vp_k + vb;  wbot = wp_k + wb
                  spd  = (cd_l / alpha_bot_l) * sqrt(ubot**2 + vbot**2 + wbot**2)
                  tb_u = spd * ubot;  tb_v = spd * vbot;  tb_w = spd * wbot
               end if

               ! Source terms (3D Coriolis cross-product -f×p*u; reduces exactly to
               ! the old f-plane form on cartesian since cfx=cfy=0 and wdp=0 there).
               sc_x = -(cfy*wdp - cfz*vdp)                                  &
                  + gravity * (init%tau_wind_df(1, Iq) - tb_u)                 &
                  - gravity * dp * init%grad_zbot_df(1, Iq)
               sc_y = -(cfz*udp - cfx*wdp)                                  &
                  + gravity * (init%tau_wind_df(2, Iq) - tb_v)                 &
                  - gravity * dp * init%grad_zbot_df(2, Iq)

               ope  = 1.0 + (dpp / pbq)
               ope2 = ope * ope
               Hq   = ope2 * H_b

               qu  = ub * udp + ope * sum_up2
               quv = vb * udp + ope * sum_uv
               qvu = ub * vdp + ope * sum_vu
               qv  = vb * vdp + ope * sum_vp2

               ! Time-average accumulators
               btp%Qu_ave_df(1,Iq) = btp%Qu_ave_df(1,Iq) + wst*(qu)
               btp%Qu_ave_df(2,Iq) = btp%Qu_ave_df(2,Iq) + wst*(quv)
               btp%Qv_ave_df(1,Iq) = btp%Qv_ave_df(1,Iq) + wst*(qvu)
               btp%Qv_ave_df(2,Iq) = btp%Qv_ave_df(2,Iq) + wst*(qv)
               btp%tau_bot_ave_df(1,Iq) = btp%tau_bot_ave_df(1,Iq) + wst*(tb_u)
               btp%tau_bot_ave_df(2,Iq) = btp%tau_bot_ave_df(2,Iq) + wst*(tb_v)
               btp%ope_ave_df(Iq) = btp%ope_ave_df(Iq) + wst*(ope)
               btp%btp_mass_flux_ave_df(1,Iq) = btp%btp_mass_flux_ave_df(1,Iq) + wst*(udp)
               btp%btp_mass_flux_ave_df(2,Iq) = btp%btp_mass_flux_ave_df(2,Iq) + wst*(vdp)
               ! Sphere-only: z-momentum flux/source terms and their time averages.
               sc_z = 0.0;  quz = 0.0;  qvz = 0.0;  qwu = 0.0;  qwv = 0.0;  qw = 0.0
               if (has_w) then
                  sc_z = -(cfx*vdp - cfy*udp)                                &
                       - gravity * tb_w                                        &
                       - gravity * dp * init%grad_zbot_df(3, Iq)

                  quz = wb * udp + ope * sum_uw
                  qvz = wb * vdp + ope * sum_vw
                  qwu = ub * wdp + ope * sum_wu
                  qwv = vb * wdp + ope * sum_wv
                  qw  = wb * wdp + ope * sum_ww

                  btp%Qu_ave_df(3,Iq) = btp%Qu_ave_df(3,Iq) + wst*(quz)
                  btp%Qv_ave_df(3,Iq) = btp%Qv_ave_df(3,Iq) + wst*(qvz)
                  btp%Qw_ave_df(1,Iq) = btp%Qw_ave_df(1,Iq) + wst*(qwu)
                  btp%Qw_ave_df(2,Iq) = btp%Qw_ave_df(2,Iq) + wst*(qwv)
                  btp%Qw_ave_df(3,Iq) = btp%Qw_ave_df(3,Iq) + wst*(qw)
                  btp%tau_bot_ave_df(3,Iq) = btp%tau_bot_ave_df(3,Iq) + wst*(tb_w)
                  btp%btp_mass_flux_ave_df(3,Iq) = btp%btp_mass_flux_ave_df(3,Iq) + wst*(wdp)
               end if

               if (ape_l) then
                  sc_x = sc_x + ape_force_btp_df(1, Iq)
                  sc_y = sc_y + ape_force_btp_df(2, Iq)
                  sc_z = sc_z + ape_force_btp_df(3, Iq)
               end if

               ! Sum-factorised projection. The weak-form weight of node J at this node is
               ! h_e (grad xi . F) + h_n (grad eta . F) + hi (zeta . F + S): h_e is nonzero
               ! only on this node's row (dpsix(l,kx)), h_n only on its column (dpsiy(l,ky)),
               ! hi only at the node itself.

               ! Nodal metrics
               e_x = mt%ksi_x(kx,ky,1,ie);  e_y = mt%ksi_y(kx,ky,1,ie);  e_z = mt%ksi_z(kx,ky,1,ie)
               n_x = mt%eta_x(kx,ky,1,ie);  n_y = mt%eta_y(kx,ky,1,ie);  n_z = mt%eta_z(kx,ky,1,ie)

               ! Node values for the split-form pressure (below): ope * interface pressures,
               ! J grad xi and J grad eta
               jac = wq / (b%wglx(kx) * b%wgly(ky))
               do m = 1, inp%nlayers
                  ope_pp_e(m,kx,ky) = ope * btp%bcl_pprime_df(m,Iq)
               end do
               jac_e(1,kx,ky) = jac * e_x;  jac_e(2,kx,ky) = jac * e_y;  jac_e(3,kx,ky) = jac * e_z
               jac_n(1,kx,ky) = jac * n_x;  jac_n(2,kx,ky) = jac * n_y;  jac_n(3,kx,ky) = jac * n_z

               ! Fluxes F(direction, equation) and sources S(equation): mass, u, v, w.
               ! The pressure Hq is not in F: it enters in split form below (xi, eta)
               ! and pointwise in the zeta term.
               F(1,1) = udp;       F(2,1) = vdp;       F(3,1) = wdp
               F(1,2) = qu;        F(2,2) = quv;       F(3,2) = quz
               F(1,3) = qvu;       F(2,3) = qv;        F(3,3) = qvz
               F(1,4) = qwu;       F(2,4) = qwv;       F(3,4) = qw

               S(1) = 0.0;  S(2) = sc_x;  S(3) = sc_y;  S(4) = sc_z

               ! Weighted contravariant fluxes
               if (has_w) then
                  c_x = mt%zeta_x(kx,ky,1,ie);  c_y = mt%zeta_y(kx,ky,1,ie);  c_z = mt%zeta_z(kx,ky,1,ie)
                  Fe(1:neq) = wq * (e_x*F(1,1:neq) + e_y*F(2,1:neq) + e_z*F(3,1:neq))
                  Fn(1:neq) = wq * (n_x*F(1,1:neq) + n_y*F(2,1:neq) + n_z*F(3,1:neq))
                  Fc(1:neq) = wq * (c_x*F(1,1:neq) + c_y*F(2,1:neq) + c_z*F(3,1:neq) + S(1:neq))
                  Fc(2) = Fc(2) + wq * c_x * Hq
                  Fc(3) = Fc(3) + wq * c_y * Hq
                  Fc(4) = Fc(4) + wq * c_z * Hq
               else
                  Fe(1:neq) = wq * (e_x*F(1,1:neq) + e_y*F(2,1:neq))
                  Fn(1:neq) = wq * (n_x*F(1,1:neq) + n_y*F(2,1:neq))
                  Fc(1:neq) = wq * S(1:neq)
               end if

               ! dHde Terms: nodes on the xi line through this node
               do l = 1, b%nglx
                  I = G%intma(l, ky, 1, ie)
                  rhs_btp(1:neq, I) = rhs_btp(1:neq, I) + b%dpsix(l, kx) * Fe(1:neq)
               end do

               ! dHdn Terms: nodes on the eta line through this node
               do l = 1, b%ngly
                  I = G%intma(kx, l, 1, ie)
                  rhs_btp(1:neq, I) = rhs_btp(1:neq, I) + b%dpsiy(l, ky) * Fn(1:neq)
               end do

               ! dHdc Terms and sources: this node only
               rhs_btp(1:neq, Iq) = rhs_btp(1:neq, Iq) + Fc(1:neq)

            end do
         end do

         ! Pressure in split form (well-balanced at rest). For each node pair (i,j) on a
         ! line: (w_j D_ji - w_i D_ij) {J grad xi}_ij H_bcl_qij, with
         ! H_bcl_qij = sum_m (beta_m/2) (ope A_m)_i (ope A_m)_j (docs/redesign/btp_collocation_wellbalanced).
         ! The weight is antisymmetric and H_bcl_qij and the metric average are symmetric, so
         ! each pair is formed once: added to node i, subtracted from node j.

         ! dHde Terms: pairs on each xi line
         do ky = 1, b%ngly
            do kx = 1, b%nglx - 1
               I = G%intma(kx, ky, 1, ie)
               do l = kx + 1, b%nglx
                  Il = G%intma(l, ky, 1, ie)

                  H_bcl_qij = 0.0
                  do m = 1, inp%nlayers
                     H_bcl_qij = H_bcl_qij + btp%bcl_dalpha(m) * ope_pp_e(m,kx,ky) * ope_pp_e(m,l,ky)
                  end do
                  f_ij = 0.5 * b%wgly(ky) * h_e(kx,l) * H_bcl_qij

                  do ic = 1, nv
                     t_ij = f_ij * (jac_e(ic,kx,ky) + jac_e(ic,l,ky))
                     rhs_btp(ic+1, I)  = rhs_btp(ic+1, I)  + t_ij
                     rhs_btp(ic+1, Il) = rhs_btp(ic+1, Il) - t_ij
                  end do

               end do
            end do
         end do

         ! dHdn Terms: pairs on each eta line
         do kx = 1, b%nglx
            do ky = 1, b%ngly - 1
               I = G%intma(kx, ky, 1, ie)
               do l = ky + 1, b%ngly
                  Il = G%intma(kx, l, 1, ie)

                  H_bcl_qij = 0.0
                  do m = 1, inp%nlayers
                     H_bcl_qij = H_bcl_qij + btp%bcl_dalpha(m) * ope_pp_e(m,kx,ky) * ope_pp_e(m,kx,l)
                  end do
                  f_ij = 0.5 * b%wglx(kx) * h_n(ky,l) * H_bcl_qij

                  do ic = 1, nv
                     t_ij = f_ij * (jac_n(ic,kx,ky) + jac_n(ic,kx,l))
                     rhs_btp(ic+1, I)  = rhs_btp(ic+1, I)  + t_ij
                     rhs_btp(ic+1, Il) = rhs_btp(ic+1, Il) - t_ij
                  end do

               end do
            end do
         end do

      end do

   end subroutine btp_colloc_volume_cpu

   subroutine btp_colloc_volume_gpu(G, b, inp, init, btp, tsp, mt, rhs_btp, qb_df)

      ! Collocated barotropic volume term (GPU). Same terms as btp_colloc_volume_cpu, in two passes:
      ! 1. at each node: physics, time averages, weighted contravariant fluxes Fe, Fn and the
      !    node term Fc (sources and zeta pressure), and ope*A_m for the split-form pressure;
      ! 2. at each node: Fc + its xi row (dHde) + its eta column (dHdn) + its split-form
      !    pressure pairs on both lines. Each thread writes only its own node.

      use mod_grid,          only: grid
      use mod_basis,         only: basis
      use mod_constants,     only: gravity
      use mod_initial,       only: initial
      use mod_variables,     only: btp_CS
      use mod_input,         only: input
      use mod_tensor,        only: tensor_CS
      use mod_initial_mlswe, only: check_btp_thickness
      use mod_metrics,       only: metrics

      implicit none

      type(grid),      intent(in)    :: G
      type(basis),     intent(in)    :: b
      type(input),     intent(in)    :: inp
      type(initial),   intent(in)    :: init
      type(btp_CS),    intent(inout) :: btp
      type(tensor_CS), intent(in)    :: tsp
      type(metrics),   intent(in)    :: mt

      real, dimension(inp%nvar_btp,G%npoin),   intent(in)  :: qb_df
      real, dimension(inp%nvar_btp-1,G%npoin), intent(out) :: rhs_btp

      real    :: sc_x, sc_y, sc_z, Hq, qu, quv, qvu, qv, quz, qvz, qwu, qwv, qw, H_b
      real    :: wq, tb_u, tb_v, tb_w, ope, ope2, g_l, wst_l
      real    :: dp, dpp, udp, vdp, wdp, ub, vb, wb, ubot, vbot, wbot, spd, pbq
      real    :: pp_k, up_k, vp_k, wp_k
      real    :: sum_up2, sum_uv, sum_vu, sum_vp2, sum_uw, sum_vw, sum_ww, sum_wu, sum_wv
      real    :: cfx, cfy, cfz
      real    :: e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z
      real    :: F(3,4), S(4)
      real    :: r(4), h_e, h_n, H_bcl_qij, f_ij, jac_i(3), jac_j(3)
      integer :: I, Iq, ie, kx, ky, l, m, ic, ieq, neq, nv, nelem_l, nglx_l, ngly_l, nlayers_l
      integer :: botfr_l
      real    :: cd_l, alpha_bot_l
      logical :: has_w

      botfr_l     = inp%botfr
      cd_l        = real(inp%cd_mlswe)
      alpha_bot_l = real(init%alpha_mlswe(inp%nlayers))
      has_w       = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex
      g_l         = gravity
      wst_l       = wst
      nelem_l     = G%nelem
      nglx_l      = b%nglx
      ngly_l      = b%ngly
      nlayers_l   = inp%nlayers
      neq = 3
      if (has_w) neq = 4
      nv  = neq - 1

      ! Abort with element/rank diagnostics if the barotropic depth went negative or NaN
      call check_btp_thickness(b, G, inp, qb_df, 'btp_colloc_volume_gpu')

      ! 1. Physics, time averages and weighted fluxes at each node
      !$acc parallel loop gang vector collapse(3) async(1)                                       &
      !$acc    present(G%intma, qb_df, init%pbprime_df, init%coriolis_df, init%kvector,          &
      !$acc            init%tau_wind_df, init%grad_zbot_df, tsp%wjac_df,                         &
      !$acc            mt%ksi_x, mt%ksi_y, mt%ksi_z, mt%eta_x, mt%eta_y, mt%eta_z,               &
      !$acc            mt%zeta_x, mt%zeta_y, mt%zeta_z,                                          &
      !$acc            btp%bcl_H_df, btp%bcl_flux_df, btp%bcl_btp_flux_df, btp%bcl_pprime_df,    &
      !$acc            btp%ope_ave_df, btp%Qu_ave_df, btp%Qv_ave_df, btp%Qw_ave_df,              &
      !$acc            btp%tau_bot_ave_df, btp%btp_mass_flux_ave_df,                             &
      !$acc            btp%colloc_flux_df, btp%ope_pp_df)                                        &
      !$acc    private(Iq, m, dp, dpp, udp, vdp, wdp, ub, vb, wb, ubot, vbot, wbot, spd, pbq,    &
      !$acc            pp_k, up_k, vp_k, wp_k, H_b, Hq, ope, ope2, wq,                           &
      !$acc            sum_up2, sum_uv, sum_vu, sum_vp2, sum_uw, sum_vw, sum_ww, sum_wu, sum_wv, &
      !$acc            cfx, cfy, cfz, tb_u, tb_v, tb_w, sc_x, sc_y, sc_z,                        &
      !$acc            qu, quv, qvu, qv, quz, qvz, qwu, qwv, qw,                                 &
      !$acc            e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, F, S, ieq)                    &
      !$acc    firstprivate(has_w, neq, botfr_l, cd_l, alpha_bot_l, g_l, wst_l, nlayers_l)
      do ie = 1, nelem_l
         do ky = 1, ngly_l
            do kx = 1, nglx_l

               Iq = G%intma(kx, ky, 1, ie)

               ! Collocation: the barotropic state at the node itself
               dp  = qb_df(1, Iq);  dpp = qb_df(2, Iq)
               udp = qb_df(3, Iq);  vdp = qb_df(4, Iq)
               wdp = 0.0
               if (has_w) wdp = qb_df(5, Iq)

               ! Layer coefficients of this barotropic sub-cycle
               H_b     = btp%bcl_H_df(Iq)
               sum_up2 = btp%bcl_flux_df(1,1,Iq)
               sum_uv  = btp%bcl_flux_df(1,2,Iq)
               sum_vu  = btp%bcl_flux_df(2,1,Iq)
               sum_vp2 = btp%bcl_flux_df(2,2,Iq)
               pp_k    = btp%bcl_btp_flux_df(1,Iq)
               up_k    = btp%bcl_btp_flux_df(2,Iq)
               vp_k    = btp%bcl_btp_flux_df(3,Iq)
               pbq     = init%pbprime_df(Iq)

               sum_uw = 0.0;  sum_vw = 0.0;  sum_ww = 0.0;  wp_k = 0.0
               sum_wu = 0.0;  sum_wv = 0.0
               if (has_w) then
                  sum_uw = btp%bcl_flux_df(3,1,Iq)
                  sum_vw = btp%bcl_flux_df(3,2,Iq)
                  sum_wu = btp%bcl_flux_df(1,3,Iq)
                  sum_wv = btp%bcl_flux_df(2,3,Iq)
                  sum_ww = btp%bcl_flux_df(3,3,Iq)
                  wp_k   = btp%bcl_btp_flux_df(4,Iq)
               end if

               wq = tsp%wjac_df(Iq)
               ub = udp / dp;  vb = vdp / dp;  wb = wdp / dp

               ! 3D Coriolis vector (f * r̂); the f-plane pair (0,0,f) on cartesian
               cfx = init%coriolis_df(Iq)*init%kvector(1, Iq)
               cfy = init%coriolis_df(Iq)*init%kvector(2, Iq)
               cfz = init%coriolis_df(Iq)*init%kvector(3, Iq)

               tb_u = 0.0;  tb_v = 0.0;  tb_w = 0.0
               if (botfr_l == 1) then
                  ubot = up_k + ub;  vbot = vp_k + vb;  wbot = wp_k + wb
                  spd  = (cd_l / g_l) * pp_k
                  tb_u = spd * ubot;  tb_v = spd * vbot;  tb_w = spd * wbot
               elseif (botfr_l == 2) then
                  ubot = up_k + ub;  vbot = vp_k + vb;  wbot = wp_k + wb
                  spd  = (cd_l / alpha_bot_l) * sqrt(ubot**2 + vbot**2 + wbot**2)
                  tb_u = spd * ubot;  tb_v = spd * vbot;  tb_w = spd * wbot
               end if

               ! Sources (3D Coriolis cross-product, wind, bottom stress, bottom slope)
               sc_x = -(cfy*wdp - cfz*vdp)                                  &
                  + g_l * (init%tau_wind_df(1, Iq) - tb_u)                  &
                  - g_l * dp * init%grad_zbot_df(1, Iq)
               sc_y = -(cfz*udp - cfx*wdp)                                  &
                  + g_l * (init%tau_wind_df(2, Iq) - tb_v)                  &
                  - g_l * dp * init%grad_zbot_df(2, Iq)

               ope  = 1.0 + (dpp / pbq)
               ope2 = ope * ope
               Hq   = ope2 * H_b

               qu  = ub * udp + ope * sum_up2
               quv = vb * udp + ope * sum_uv
               qvu = ub * vdp + ope * sum_vu
               qv  = vb * vdp + ope * sum_vp2

               ! Time-average accumulators (one thread per node)
               btp%Qu_ave_df(1,Iq) = btp%Qu_ave_df(1,Iq) + wst_l*(qu)
               btp%Qu_ave_df(2,Iq) = btp%Qu_ave_df(2,Iq) + wst_l*(quv)
               btp%Qv_ave_df(1,Iq) = btp%Qv_ave_df(1,Iq) + wst_l*(qvu)
               btp%Qv_ave_df(2,Iq) = btp%Qv_ave_df(2,Iq) + wst_l*(qv)
               btp%tau_bot_ave_df(1,Iq) = btp%tau_bot_ave_df(1,Iq) + wst_l*(tb_u)
               btp%tau_bot_ave_df(2,Iq) = btp%tau_bot_ave_df(2,Iq) + wst_l*(tb_v)
               btp%ope_ave_df(Iq) = btp%ope_ave_df(Iq) + wst_l*(ope)
               btp%btp_mass_flux_ave_df(1,Iq) = btp%btp_mass_flux_ave_df(1,Iq) + wst_l*(udp)
               btp%btp_mass_flux_ave_df(2,Iq) = btp%btp_mass_flux_ave_df(2,Iq) + wst_l*(vdp)

               ! Sphere-only: z-momentum flux/source terms and their time averages
               sc_z = 0.0;  quz = 0.0;  qvz = 0.0;  qwu = 0.0;  qwv = 0.0;  qw = 0.0
               if (has_w) then
                  sc_z = -(cfx*vdp - cfy*udp)                                &
                       - g_l * tb_w                                          &
                       - g_l * dp * init%grad_zbot_df(3, Iq)

                  quz = wb * udp + ope * sum_uw
                  qvz = wb * vdp + ope * sum_vw
                  qwu = ub * wdp + ope * sum_wu
                  qwv = vb * wdp + ope * sum_wv
                  qw  = wb * wdp + ope * sum_ww

                  btp%Qu_ave_df(3,Iq) = btp%Qu_ave_df(3,Iq) + wst_l*(quz)
                  btp%Qv_ave_df(3,Iq) = btp%Qv_ave_df(3,Iq) + wst_l*(qvz)
                  btp%Qw_ave_df(1,Iq) = btp%Qw_ave_df(1,Iq) + wst_l*(qwu)
                  btp%Qw_ave_df(2,Iq) = btp%Qw_ave_df(2,Iq) + wst_l*(qwv)
                  btp%Qw_ave_df(3,Iq) = btp%Qw_ave_df(3,Iq) + wst_l*(qw)
                  btp%tau_bot_ave_df(3,Iq) = btp%tau_bot_ave_df(3,Iq) + wst_l*(tb_w)
                  btp%btp_mass_flux_ave_df(3,Iq) = btp%btp_mass_flux_ave_df(3,Iq) + wst_l*(wdp)
               end if

               ! ope * A_m for the split-form pressure (pass 2)
               !$acc loop seq
               do m = 1, nlayers_l
                  btp%ope_pp_df(m,Iq) = ope * btp%bcl_pprime_df(m,Iq)
               end do

               ! Nodal metrics
               e_x = mt%ksi_x(kx,ky,1,ie);  e_y = mt%ksi_y(kx,ky,1,ie);  e_z = mt%ksi_z(kx,ky,1,ie)
               n_x = mt%eta_x(kx,ky,1,ie);  n_y = mt%eta_y(kx,ky,1,ie);  n_z = mt%eta_z(kx,ky,1,ie)

               ! Fluxes F(direction, equation) and sources S(equation): mass, u, v, w.
               ! The pressure Hq is not in F: split form in pass 2, pointwise in zeta.
               F(1,1) = udp;       F(2,1) = vdp;       F(3,1) = wdp
               F(1,2) = qu;        F(2,2) = quv;       F(3,2) = quz
               F(1,3) = qvu;       F(2,3) = qv;        F(3,3) = qvz
               F(1,4) = qwu;       F(2,4) = qwv;       F(3,4) = qw

               S(1) = 0.0;  S(2) = sc_x;  S(3) = sc_y;  S(4) = sc_z

               ! Weighted contravariant fluxes colloc_flux_df(direction, equation, node)
               if (has_w) then
                  c_x = mt%zeta_x(kx,ky,1,ie);  c_y = mt%zeta_y(kx,ky,1,ie);  c_z = mt%zeta_z(kx,ky,1,ie)
                  !$acc loop seq
                  do ieq = 1, neq
                     btp%colloc_flux_df(1,ieq,Iq) = wq * (e_x*F(1,ieq) + e_y*F(2,ieq) + e_z*F(3,ieq))
                     btp%colloc_flux_df(2,ieq,Iq) = wq * (n_x*F(1,ieq) + n_y*F(2,ieq) + n_z*F(3,ieq))
                     btp%colloc_flux_df(3,ieq,Iq) = wq * (c_x*F(1,ieq) + c_y*F(2,ieq) + c_z*F(3,ieq) + S(ieq))
                  end do
                  btp%colloc_flux_df(3,2,Iq) = btp%colloc_flux_df(3,2,Iq) + wq * c_x * Hq
                  btp%colloc_flux_df(3,3,Iq) = btp%colloc_flux_df(3,3,Iq) + wq * c_y * Hq
                  btp%colloc_flux_df(3,4,Iq) = btp%colloc_flux_df(3,4,Iq) + wq * c_z * Hq
               else
                  !$acc loop seq
                  do ieq = 1, neq
                     btp%colloc_flux_df(1,ieq,Iq) = wq * (e_x*F(1,ieq) + e_y*F(2,ieq))
                     btp%colloc_flux_df(2,ieq,Iq) = wq * (n_x*F(1,ieq) + n_y*F(2,ieq))
                     btp%colloc_flux_df(3,ieq,Iq) = wq * S(ieq)
                  end do
               end if

            end do
         end do
      end do
      !$acc end parallel loop

      ! 2. Each node collects its own term, row, column and split-form pressure pairs
      !$acc parallel loop gang vector collapse(3) async(1)                                       &
      !$acc    present(G%intma, b%dpsix, b%dpsiy, b%wglx, b%wgly, tsp%wjac_df,                   &
      !$acc            mt%ksi_x, mt%ksi_y, mt%ksi_z, mt%eta_x, mt%eta_y, mt%eta_z,               &
      !$acc            btp%colloc_flux_df, btp%ope_pp_df, btp%bcl_dalpha, rhs_btp)               &
      !$acc    private(Iq, I, l, m, ic, ieq, r, h_e, h_n, H_bcl_qij, f_ij, jac_i, jac_j)          &
      !$acc    firstprivate(neq, nv, nglx_l, ngly_l, nlayers_l)
      do ie = 1, nelem_l
         do ky = 1, ngly_l
            do kx = 1, nglx_l

               Iq = G%intma(kx, ky, 1, ie)

               ! dHdc Terms and sources: this node
               !$acc loop seq
               do ieq = 1, neq
                  r(ieq) = btp%colloc_flux_df(3,ieq,Iq)
               end do

               ! dHde Terms: the xi line through this node
               !$acc loop seq
               do l = 1, nglx_l
                  I = G%intma(l, ky, 1, ie)
                  !$acc loop seq
                  do ieq = 1, neq
                     r(ieq) = r(ieq) + b%dpsix(kx, l) * btp%colloc_flux_df(1,ieq,I)
                  end do
               end do

               ! dHdn Terms: the eta line through this node
               !$acc loop seq
               do l = 1, ngly_l
                  I = G%intma(kx, l, 1, ie)
                  !$acc loop seq
                  do ieq = 1, neq
                     r(ieq) = r(ieq) + b%dpsiy(ky, l) * btp%colloc_flux_df(2,ieq,I)
                  end do
               end do

               ! Split-form pressure: pairs (this node, l) on the xi line,
               ! 0.5 w_eta (w_l D_il - w_i D_li) {J grad xi} H_bcl_qij
               jac_i(1) = tsp%wjac_df(Iq) / (b%wglx(kx) * b%wgly(ky)) * mt%ksi_x(kx,ky,1,ie)
               jac_i(2) = tsp%wjac_df(Iq) / (b%wglx(kx) * b%wgly(ky)) * mt%ksi_y(kx,ky,1,ie)
               jac_i(3) = tsp%wjac_df(Iq) / (b%wglx(kx) * b%wgly(ky)) * mt%ksi_z(kx,ky,1,ie)
               !$acc loop seq
               do l = 1, nglx_l
                  if (l /= kx) then
                     I = G%intma(l, ky, 1, ie)
                     H_bcl_qij = 0.0
                     !$acc loop seq
                     do m = 1, nlayers_l
                        H_bcl_qij = H_bcl_qij + btp%bcl_dalpha(m) * btp%ope_pp_df(m,Iq) * btp%ope_pp_df(m,I)
                     end do
                     h_e  = b%wglx(l)*b%dpsix(kx,l) - b%wglx(kx)*b%dpsix(l,kx)
                     f_ij = 0.5 * b%wgly(ky) * h_e * H_bcl_qij
                     jac_j(1) = tsp%wjac_df(I) / (b%wglx(l) * b%wgly(ky)) * mt%ksi_x(l,ky,1,ie)
                     jac_j(2) = tsp%wjac_df(I) / (b%wglx(l) * b%wgly(ky)) * mt%ksi_y(l,ky,1,ie)
                     jac_j(3) = tsp%wjac_df(I) / (b%wglx(l) * b%wgly(ky)) * mt%ksi_z(l,ky,1,ie)
                     !$acc loop seq
                     do ic = 1, nv
                        r(ic+1) = r(ic+1) + f_ij * (jac_i(ic) + jac_j(ic))
                     end do
                  end if
               end do

               ! Split-form pressure: pairs (this node, l) on the eta line
               jac_i(1) = tsp%wjac_df(Iq) / (b%wglx(kx) * b%wgly(ky)) * mt%eta_x(kx,ky,1,ie)
               jac_i(2) = tsp%wjac_df(Iq) / (b%wglx(kx) * b%wgly(ky)) * mt%eta_y(kx,ky,1,ie)
               jac_i(3) = tsp%wjac_df(Iq) / (b%wglx(kx) * b%wgly(ky)) * mt%eta_z(kx,ky,1,ie)
               !$acc loop seq
               do l = 1, ngly_l
                  if (l /= ky) then
                     I = G%intma(kx, l, 1, ie)
                     H_bcl_qij = 0.0
                     !$acc loop seq
                     do m = 1, nlayers_l
                        H_bcl_qij = H_bcl_qij + btp%bcl_dalpha(m) * btp%ope_pp_df(m,Iq) * btp%ope_pp_df(m,I)
                     end do
                     h_n  = b%wgly(l)*b%dpsiy(ky,l) - b%wgly(ky)*b%dpsiy(l,ky)
                     f_ij = 0.5 * b%wglx(kx) * h_n * H_bcl_qij
                     jac_j(1) = tsp%wjac_df(I) / (b%wglx(kx) * b%wgly(l)) * mt%eta_x(kx,l,1,ie)
                     jac_j(2) = tsp%wjac_df(I) / (b%wglx(kx) * b%wgly(l)) * mt%eta_y(kx,l,1,ie)
                     jac_j(3) = tsp%wjac_df(I) / (b%wglx(kx) * b%wgly(l)) * mt%eta_z(kx,l,1,ie)
                     !$acc loop seq
                     do ic = 1, nv
                        r(ic+1) = r(ic+1) + f_ij * (jac_i(ic) + jac_j(ic))
                     end do
                  end if
               end do

               ! Store RHS
               !$acc loop seq
               do ieq = 1, neq
                  rhs_btp(ieq, Iq) = r(ieq)
               end do

            end do
         end do
      end do
      !$acc end parallel loop

   end subroutine btp_colloc_volume_gpu

   subroutine btp_colloc_faces(G, b, inp, mf, init, btp, rhs_btp, qb, qprime_df)
      !===========================================================================
      !  Collocated barotropic face term (Riemann solver). Same physics as
      !  create_btp_fluxes_qdf, evaluated at the GLL nodes of each face.
      !  CPU and GPU: on the GPU one thread per face node; the RHS update is
      !  atomic (corner nodes are shared by two faces of an element).
      !===========================================================================
      use mod_basis,     only: basis
      use mod_grid,      only: grid
      use mod_face,      only: face_CS
      use mod_variables, only: btp_CS
      use mod_initial,   only: initial
      use mod_constants, only: gravity
      use mod_input,     only: input

      implicit none

      type(grid),    intent(in)    :: G
      type(basis),   intent(in)    :: b
      type(input),   intent(in)    :: inp
      type(face_CS), intent(in)    :: mf
      type(initial), intent(in)    :: init
      type(btp_CS),  intent(inout) :: btp

      real, dimension(inp%nvar_btp-1,G%npoin),             intent(inout) :: rhs_btp
      real, dimension(inp%nvar_btp,G%npoin),               intent(in)    :: qb
      real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers),   intent(in)    :: qprime_df

      integer :: iface, iquad, el, er, I_l, I_r, k
      real    :: wq, nxl, nyl, nzl, un
      real    :: pbl, pbr, clam, one_eta, one_eta2, half_clam, quarter_clam
      real    :: pU_L, pU_R, half_clam_dqb2, pbpert_edge
      real    :: qbl(5), qbr(5)
      real    :: c_minus, c_plus
      real    :: ul, ur, vl, vr, wl, wr, opl, opr
      real    :: oe2_Hql, oe2_Hqr
      real    :: H_bcl_qlr   ! split-form pressure: sum_m (beta_m/2) A_m,L A_m,R
      real    :: Qu_ql1, Qu_ql2, Qu_ql3, Qu_qr1, Qu_qr2, Qu_qr3
      real    :: Qv_ql1, Qv_ql2, Qv_ql3, Qv_qr1, Qv_qr2, Qv_qr3
      real    :: Qw_ql1, Qw_ql2, Qw_ql3, Qw_qr1, Qw_qr2, Qw_qr3
      real    :: flux_edge_x, flux_edge_y, flux_edge_z, flux_pb, flux_u, flux_v, flux_w
      real    :: pprime_lk, pprime_rk, pprime_lk1, pprime_rk1
      real    :: pkl, pkr, ukl, ukr, vkl, vkr, wkl, wkr
      real    :: ope_ppl_k, ope_ppr_k, uv_cross_l, uv_cross_r
      real    :: wq_flux_pb, wq_flux_u, wq_flux_v, wq_flux_w
      real    :: g_l, dry_l, wst_l
      integer :: ngl_f, nlayers_f, nface_f
      logical :: has_w

      ! Scalar copies (firstprivate on the GPU)
      ngl_f     = b%ngl
      nlayers_f = inp%nlayers
      nface_f   = G%nface
      has_w     = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex
      g_l       = gravity
      dry_l     = inp%dry_cutoff
      wst_l     = wst

      !$acc parallel loop gang vector collapse(2) async(1)                                     &
      !$acc    present(G%face, G%face_type, G%intma, mf%imapl, mf%imapr, mf%normal_vector,     &
      !$acc            mf%jac_face, qb, qprime_df, init%pbprime_df, init%alpha_mlswe,          &
      !$acc            btp%bcl_dalpha, btp%bcl_pprime_face_df, rhs_btp,                        &
      !$acc            btp%btp_mass_flux_face_ave_df, btp%Qu_face_ave_df, btp%Qv_face_ave_df,  &
      !$acc            btp%Qw_face_ave_df, btp%ope_face_ave_df, btp%ope2_face_ave_df,          &
      !$acc            btp%one_plus_eta_edge_2_ave_df, btp%uvb_face_ave_df)                    &
      !$acc    private(el, er, I_l, I_r, k, wq, nxl, nyl, nzl, un, pbl, pbr, clam,            &
      !$acc            one_eta, one_eta2, half_clam, quarter_clam, pU_L, pU_R, half_clam_dqb2, &
      !$acc            pbpert_edge, qbl, qbr, c_minus, c_plus, ul, ur, vl, vr, wl, wr,         &
      !$acc            opl, opr, oe2_Hql, oe2_Hqr, H_bcl_qlr,                                  &
      !$acc            Qu_ql1, Qu_ql2, Qu_ql3, Qu_qr1, Qu_qr2, Qu_qr3,                         &
      !$acc            Qv_ql1, Qv_ql2, Qv_ql3, Qv_qr1, Qv_qr2, Qv_qr3,                         &
      !$acc            Qw_ql1, Qw_ql2, Qw_ql3, Qw_qr1, Qw_qr2, Qw_qr3,                         &
      !$acc            flux_edge_x, flux_edge_y, flux_edge_z, flux_pb, flux_u, flux_v, flux_w, &
      !$acc            pprime_lk, pprime_rk, pprime_lk1, pprime_rk1,                           &
      !$acc            pkl, pkr, ukl, ukr, vkl, vkr, wkl, wkr,                                 &
      !$acc            ope_ppl_k, ope_ppr_k, uv_cross_l, uv_cross_r,                           &
      !$acc            wq_flux_pb, wq_flux_u, wq_flux_v, wq_flux_w)                            &
      !$acc    firstprivate(nlayers_f, has_w, g_l, dry_l, wst_l)
      do iface = 1, nface_f
         do iquad = 1, ngl_f

            if (G%face_type(iface) /= 2) then

               el = G%face(7, iface)
               er = G%face(8, iface)

               nxl = mf%normal_vector(1, iquad, 1, iface)
               nyl = mf%normal_vector(2, iquad, 1, iface)
               nzl = mf%normal_vector(3, iquad, 1, iface)

               ! LEFT state: collocation, the face node itself
               I_l = G%intma(mf%imapl(1,iquad,1,iface), mf%imapl(2,iquad,1,iface), mf%imapl(3,iquad,1,iface), el)
               qbl(1) = qb(1,I_l);  qbl(2) = qb(2,I_l)
               qbl(3) = qb(3,I_l);  qbl(4) = qb(4,I_l)
               qbl(5) = 0.0
               if (has_w) qbl(5) = qb(5,I_l)
               pbl = init%pbprime_df(I_l)

               ! RIGHT state or boundary condition
               I_r = I_l
               if (er > 0) then
                  I_r = G%intma(mf%imapr(1,iquad,1,iface), mf%imapr(2,iquad,1,iface), mf%imapr(3,iquad,1,iface), er)
                  qbr(1) = qb(1,I_r);  qbr(2) = qb(2,I_r)
                  qbr(3) = qb(3,I_r);  qbr(4) = qb(4,I_r)
                  qbr(5) = 0.0
                  if (has_w) qbr(5) = qb(5,I_r)
                  pbr = init%pbprime_df(I_r)
               else
                  qbr(1) = qbl(1);  qbr(2) = qbl(2)
                  qbr(3) = qbl(3);  qbr(4) = qbl(4);  qbr(5) = qbl(5)
                  pbr = pbl
                  if (er == -4) then          ! slip wall
                     un     = nxl*qbl(3) + nyl*qbl(4) + nzl*qbl(5)
                     qbr(3) = qbl(3) - 2.0*un*nxl
                     qbr(4) = qbl(4) - 2.0*un*nyl
                     if (has_w) qbr(5) = qbl(5) - 2.0*un*nzl
                  else if (er == -2) then     ! no-slip wall
                     qbr(3) = -qbl(3);  qbr(4) = -qbl(4)
                     if (has_w) qbr(5) = -qbl(5)
                  end if
               end if

               ! Wave speeds
               pU_L = nxl*qbl(3) + nyl*qbl(4) + nzl*qbl(5)
               pU_R = -(nxl*qbr(3) + nyl*qbr(4) + nzl*qbr(5))

               c_minus      = sqrt(init%alpha_mlswe(nlayers_f) * pbr)
               c_plus       = sqrt(init%alpha_mlswe(nlayers_f) * pbl)
               clam         = max(c_minus, c_plus)
               half_clam    = 0.5  * clam
               quarter_clam = 0.25 * clam

               pbpert_edge = 0.5*(qbl(2) + qbr(2)) + (0.5/clam)*(pU_L + pU_R)
               one_eta     = 1.0 + pbpert_edge / pbl
               one_eta2    = one_eta * one_eta

               ul = qbl(3) / qbl(1);  ur = qbr(3) / qbr(1)
               vl = qbl(4) / qbl(1);  vr = qbr(4) / qbr(1)
               wl = 0.0;  wr = 0.0
               if (has_w) then
                  wl = qbl(5) / qbl(1);  wr = qbr(5) / qbr(1)
               end if

               ! Initialise flux tensors from barotropic state
               Qu_ql1 = ul * qbl(3);  Qu_ql2 = vl * qbl(3);  Qu_ql3 = 0.0
               Qu_qr1 = ur * qbr(3);  Qu_qr2 = vr * qbr(3);  Qu_qr3 = 0.0
               Qv_ql1 = ul * qbl(4);  Qv_ql2 = vl * qbl(4);  Qv_ql3 = 0.0
               Qv_qr1 = ur * qbr(4);  Qv_qr2 = vr * qbr(4);  Qv_qr3 = 0.0
               Qw_ql1 = 0.0;  Qw_ql2 = 0.0;  Qw_ql3 = 0.0
               Qw_qr1 = 0.0;  Qw_qr2 = 0.0;  Qw_qr3 = 0.0
               if (has_w) then
                  Qu_ql3 = wl * qbl(3);  Qu_qr3 = wr * qbr(3)
                  Qv_ql3 = wl * qbl(4);  Qv_qr3 = wr * qbr(4)
                  Qw_ql1 = ul * qbl(5);  Qw_ql2 = vl * qbl(5);  Qw_ql3 = wl * qbl(5)
                  Qw_qr1 = ur * qbr(5);  Qw_qr2 = vr * qbr(5);  Qw_qr3 = wr * qbr(5)
               end if
               H_bcl_qlr = 0.0
               pprime_lk = 0.0;  pprime_rk = 0.0

               ! Layer contributions at the face node
               !$acc loop seq
               do k = 1, nlayers_f

                  pkl = qprime_df(1, I_l, k)
                  ukl = qprime_df(2, I_l, k)
                  vkl = qprime_df(3, I_l, k)
                  wkl = 0.0
                  if (has_w) wkl = qprime_df(4, I_l, k)

                  if (er > 0) then
                     pkr = qprime_df(1, I_r, k)
                     ukr = qprime_df(2, I_r, k)
                     vkr = qprime_df(3, I_r, k)
                     wkr = 0.0
                     if (has_w) wkr = qprime_df(4, I_r, k)
                  else
                     pkr = pkl
                     wkr = 0.0
                     if (er == -4) then
                        un  = nxl * ukl + nyl * vkl + nzl * wkl
                        ukr = ukl - 2.0 * un * nxl
                        vkr = vkl - 2.0 * un * nyl
                        if (has_w) wkr = wkl - 2.0 * un * nzl
                     else if (er == -2) then
                        ukr = -ukl;  vkr = -vkl
                        if (has_w) wkr = -wkl
                     else
                        ukr = ukl;  vkr = vkl;  wkr = wkl
                     end if
                  end if

                  ! Dry-cell protection, matching the volume/face BCL flux routines: the
                  ! thickness is clamped to a floor (the interface pressures stay continuous
                  ! across dry_cutoff) and the velocity is zeroed, which removes this layer's
                  ! Qu/Qv/Qw flux terms.
                  if (pkl < (g_l/init%alpha_mlswe(k)) * dry_l) then
                     pkl = (g_l/init%alpha_mlswe(k)) * dry_l
                     ukl = 0.0;  vkl = 0.0;  wkl = 0.0
                  end if
                  if (pkr < (g_l/init%alpha_mlswe(k)) * dry_l) then
                     pkr = (g_l/init%alpha_mlswe(k)) * dry_l
                     ukr = 0.0;  vkr = 0.0;  wkr = 0.0
                  end if

                  ope_ppl_k  = one_eta * pkl
                  ope_ppr_k  = one_eta * pkr
                  uv_cross_l = ukl * vkl * ope_ppl_k
                  uv_cross_r = ukr * vkr * ope_ppr_k

                  Qu_ql1 = Qu_ql1 + ukl * ukl * ope_ppl_k
                  Qu_ql2 = Qu_ql2 + uv_cross_l
                  Qv_ql1 = Qv_ql1 + uv_cross_l
                  Qv_ql2 = Qv_ql2 + vkl * vkl * ope_ppl_k

                  Qu_qr1 = Qu_qr1 + ukr * ukr * ope_ppr_k
                  Qu_qr2 = Qu_qr2 + uv_cross_r
                  Qv_qr1 = Qv_qr1 + uv_cross_r
                  Qv_qr2 = Qv_qr2 + vkr * vkr * ope_ppr_k

                  if (has_w) then
                     Qu_ql3 = Qu_ql3 + ukl * wkl * ope_ppl_k;  Qu_qr3 = Qu_qr3 + ukr * wkr * ope_ppr_k
                     Qv_ql3 = Qv_ql3 + vkl * wkl * ope_ppl_k;  Qv_qr3 = Qv_qr3 + vkr * wkr * ope_ppr_k
                     Qw_ql1 = Qw_ql1 + ukl * wkl * ope_ppl_k;  Qw_qr1 = Qw_qr1 + ukr * wkr * ope_ppr_k
                     Qw_ql2 = Qw_ql2 + vkl * wkl * ope_ppl_k;  Qw_qr2 = Qw_qr2 + vkr * wkr * ope_ppr_k
                     Qw_ql3 = Qw_ql3 + wkl * wkl * ope_ppl_k;  Qw_qr3 = Qw_qr3 + wkr * wkr * ope_ppr_k
                  end if

                  pprime_lk1 = pprime_lk + pkl
                  pprime_rk1 = pprime_rk + pkr
                  H_bcl_qlr  = H_bcl_qlr + btp%bcl_dalpha(k) * pprime_lk1 * pprime_rk1
                  btp%bcl_pprime_face_df(k,1,iquad,iface) = pprime_lk1
                  btp%bcl_pprime_face_df(k,2,iquad,iface) = pprime_rk1
                  pprime_lk  = pprime_lk1
                  pprime_rk  = pprime_rk1

               end do  ! k

               half_clam_dqb2 = half_clam * (qbl(2) - qbr(2))
               flux_edge_x = 0.5*(qbl(3) + qbr(3)) + nxl * half_clam_dqb2
               flux_edge_y = 0.5*(qbl(4) + qbr(4)) + nyl * half_clam_dqb2
               flux_edge_z = 0.0
               if (has_w) flux_edge_z = 0.5*(qbl(5) + qbr(5)) + nzl * half_clam_dqb2
               flux_pb     = nxl*flux_edge_x + nyl*flux_edge_y + nzl*flux_edge_z

               ! Time-average accumulators (one per face node)
               btp%btp_mass_flux_face_ave_df(1,iquad,iface) = btp%btp_mass_flux_face_ave_df(1,iquad,iface) + wst_l*(flux_edge_x)
               btp%btp_mass_flux_face_ave_df(2,iquad,iface) = btp%btp_mass_flux_face_ave_df(2,iquad,iface) + wst_l*(flux_edge_y)
               btp%Qu_face_ave_df(1,iquad,iface) = btp%Qu_face_ave_df(1,iquad,iface) + wst_l*(0.5*(Qu_ql1 + Qu_qr1))
               btp%Qu_face_ave_df(2,iquad,iface) = btp%Qu_face_ave_df(2,iquad,iface) + wst_l*(0.5*(Qu_ql2 + Qu_qr2))
               btp%Qv_face_ave_df(1,iquad,iface) = btp%Qv_face_ave_df(1,iquad,iface) + wst_l*(0.5*(Qv_ql1 + Qv_qr1))
               btp%Qv_face_ave_df(2,iquad,iface) = btp%Qv_face_ave_df(2,iquad,iface) + wst_l*(0.5*(Qv_ql2 + Qv_qr2))
               opl = 1.0 + qbl(2)/pbl
               opr = 1.0 + qbr(2)/pbr
               btp%ope_face_ave_df(1,iquad,iface) = btp%ope_face_ave_df(1,iquad,iface) + wst_l*(opl)
               btp%ope_face_ave_df(2,iquad,iface) = btp%ope_face_ave_df(2,iquad,iface) + wst_l*(opr)
               btp%ope2_face_ave_df(1,iquad,iface) = btp%ope2_face_ave_df(1,iquad,iface) + wst_l*(opl*opl)
               btp%ope2_face_ave_df(2,iquad,iface) = btp%ope2_face_ave_df(2,iquad,iface) + wst_l*(opr*opr)
               btp%one_plus_eta_edge_2_ave_df(iquad,iface) = btp%one_plus_eta_edge_2_ave_df(iquad,iface) + wst_l*(one_eta2)
               btp%uvb_face_ave_df(1,1,iquad,iface) = btp%uvb_face_ave_df(1,1,iquad,iface) + wst_l*(ul)
               btp%uvb_face_ave_df(1,2,iquad,iface) = btp%uvb_face_ave_df(1,2,iquad,iface) + wst_l*(ur)
               btp%uvb_face_ave_df(2,1,iquad,iface) = btp%uvb_face_ave_df(2,1,iquad,iface) + wst_l*(vl)
               btp%uvb_face_ave_df(2,2,iquad,iface) = btp%uvb_face_ave_df(2,2,iquad,iface) + wst_l*(vr)
               if (has_w) then
                  btp%btp_mass_flux_face_ave_df(3,iquad,iface) = btp%btp_mass_flux_face_ave_df(3,iquad,iface) + wst_l*(flux_edge_z)
                  btp%Qu_face_ave_df(3,iquad,iface) = btp%Qu_face_ave_df(3,iquad,iface) + wst_l*(0.5*(Qu_ql3 + Qu_qr3))
                  btp%Qv_face_ave_df(3,iquad,iface) = btp%Qv_face_ave_df(3,iquad,iface) + wst_l*(0.5*(Qv_ql3 + Qv_qr3))
                  btp%Qw_face_ave_df(1,iquad,iface) = btp%Qw_face_ave_df(1,iquad,iface) + wst_l*(0.5*(Qw_ql1 + Qw_qr1))
                  btp%Qw_face_ave_df(2,iquad,iface) = btp%Qw_face_ave_df(2,iquad,iface) + wst_l*(0.5*(Qw_ql2 + Qw_qr2))
                  btp%Qw_face_ave_df(3,iquad,iface) = btp%Qw_face_ave_df(3,iquad,iface) + wst_l*(0.5*(Qw_ql3 + Qw_qr3))
                  btp%uvb_face_ave_df(3,1,iquad,iface) = btp%uvb_face_ave_df(3,1,iquad,iface) + wst_l*(wl)
                  btp%uvb_face_ave_df(3,2,iquad,iface) = btp%uvb_face_ave_df(3,2,iquad,iface) + wst_l*(wr)
               end if

               ! Momentum fluxes
               ! Pressure part of the face flux in split form: 0.5*(oe2_Hql + oe2_Hqr) = one_eta2*H_bcl_qlr
               oe2_Hql = one_eta2 * H_bcl_qlr
               oe2_Hqr = one_eta2 * H_bcl_qlr

               Qu_ql1 = Qu_ql1 + oe2_Hql
               Qu_qr1 = Qu_qr1 + oe2_Hqr
               Qv_ql2 = Qv_ql2 + oe2_Hql
               Qv_qr2 = Qv_qr2 + oe2_Hqr
               if (has_w) then
                  Qw_ql3 = Qw_ql3 + oe2_Hql
                  Qw_qr3 = Qw_qr3 + oe2_Hqr
               end if

               flux_u = 0.5 * (nxl*(Qu_ql1 + Qu_qr1) + nyl*(Qu_ql2 + Qu_qr2) + nzl*(Qu_ql3 + Qu_qr3)) &
                  - quarter_clam * (qbr(3) - qbl(3))
               flux_v = 0.5 * (nxl*(Qv_ql1 + Qv_qr1) + nyl*(Qv_ql2 + Qv_qr2) + nzl*(Qv_ql3 + Qv_qr3)) &
                  - quarter_clam * (qbr(4) - qbl(4))
               flux_w = 0.0
               if (has_w) flux_w = 0.5 * (nxl*(Qw_ql1 + Qw_qr1) + nyl*(Qw_ql2 + Qw_qr2) + nzl*(Qw_ql3 + Qw_qr3)) &
                  - quarter_clam * (qbr(5) - qbl(5))

               wq         = mf%jac_face(iquad, 1, iface)
               wq_flux_pb = wq * flux_pb
               wq_flux_u  = wq * flux_u
               wq_flux_v  = wq * flux_v
               wq_flux_w  = wq * flux_w

               ! Store RHS: left and right sides
               !$acc atomic update
               rhs_btp(1,I_l) = rhs_btp(1,I_l) - wq_flux_pb
               !$acc atomic update
               rhs_btp(2,I_l) = rhs_btp(2,I_l) - wq_flux_u
               !$acc atomic update
               rhs_btp(3,I_l) = rhs_btp(3,I_l) - wq_flux_v
               if (has_w) then
                  !$acc atomic update
                  rhs_btp(4,I_l) = rhs_btp(4,I_l) - wq_flux_w
               end if

               if (er > 0) then
                  !$acc atomic update
                  rhs_btp(1,I_r) = rhs_btp(1,I_r) + wq_flux_pb
                  !$acc atomic update
                  rhs_btp(2,I_r) = rhs_btp(2,I_r) + wq_flux_u
                  !$acc atomic update
                  rhs_btp(3,I_r) = rhs_btp(3,I_r) + wq_flux_v
                  if (has_w) then
                     !$acc atomic update
                     rhs_btp(4,I_r) = rhs_btp(4,I_r) + wq_flux_w
                  end if
               end if

            end if

         end do  ! iquad
      end do  ! iface
      !$acc end parallel loop

   end subroutine btp_colloc_faces

subroutine btp_colloc_nbhs_faces(G, inp, b, mf, par, btp, init, ref, rhs)

   ! Collocated barotropic MPI-face term, left side only; the right state comes from the
   ! neighbouring rank (ref%q_recv). CPU and GPU: on the GPU one thread per face node;
   ! the RHS update is atomic (corner nodes are shared by two faces of an element).

   use mod_grid,      only: grid
   use mod_input,     only: input
   use mod_basis,     only: basis
   use mod_face,      only: face_CS
   use mod_parallel,  only: parallel_CS
   use mod_variables, only: btp_CS
   use mod_initial,   only: initial
   use mod_ref,       only: mref
   use mod_constants, only: gravity

   implicit none

   type(grid),        intent(in)    :: G
   type(input),       intent(in)    :: inp
   type(basis),       intent(in)    :: b
   type(face_CS),     intent(in)    :: mf
   type(parallel_CS), intent(in)    :: par
   type(btp_CS),      intent(inout) :: btp
   type(initial),     intent(in)    :: init
   type(mref),        intent(in)    :: ref

   real,           intent(inout) :: rhs(inp%nvar_btp-1, G%npoin)

   integer :: kk, iface, iquad, el, I, k, ii, iv, nboun_valid
   real :: wq, nxl, nyl, nzl, nxr, nyr, nzr
   real :: ul, ur, vl, vr, wl, wr, pbl, pbr, clam, one_eta
   real :: pU_L, pU_R, pbpert_edge, c_minus, c_plus
   real :: qbl(5), qbr(5)
   real :: ppl, ppr, upl, upr, vpl, vpr, wpl, wpr
   real :: pprime_lk, pprime_rk
   real :: flux_pb, flux_u, flux_v, flux_w, fxl, fxr
   real :: oe2_Hql, oe2_Hqr
   real :: H_bcl_qlr   ! split-form pressure: sum_m (beta_m/2) A_m,L A_m,R
   real :: flux_edge_x, flux_edge_y, flux_edge_z
   real, dimension(3) :: Qu_ql, Qu_qr, Qv_ql, Qv_qr, Qw_ql, Qw_qr
   integer :: ngl_f, nlayers_f, nvarb_f, bcl_stride_f, pb_idx_f
   real    :: dry_cutoff_f, g_l, wst_l
   logical :: has_w

   ngl_f        = b%ngl
   nlayers_f    = inp%nlayers
   nboun_valid  = ref%nboun_valid
   nvarb_f      = inp%nvar_btp
   bcl_stride_f = inp%nvar_bcl
   dry_cutoff_f = inp%dry_cutoff   ! scalar copy: inp itself is not on the device
   pb_idx_f     = inp%nvar_btp + 1  ! pbprime is packed right after the momentum block
   has_w        = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex
   g_l          = gravity
   wst_l        = wst

   if (nboun_valid < 1) return

   ! Store RHS: left side only; the right element belongs to the neighbouring rank.
   !$acc parallel loop gang vector collapse(2) async(1)                                     &
   !$acc    present(ref%face_pack_list, ref%q_send, ref%q_recv, G%face, G%intma,            &
   !$acc            mf%imapl, mf%normal_vector, mf%jac_face, init%alpha_mlswe, rhs,         &
   !$acc            btp%bcl_dalpha, btp%bcl_pprime_face_df,                                 &
   !$acc            btp%btp_mass_flux_face_ave_df, btp%Qu_face_ave_df, btp%Qv_face_ave_df,  &
   !$acc            btp%Qw_face_ave_df, btp%ope_face_ave_df, btp%ope2_face_ave_df,          &
   !$acc            btp%one_plus_eta_edge_2_ave_df, btp%uvb_face_ave_df)                    &
   !$acc    private(iface, el, I, k, ii, iv, wq, nxl, nyl, nzl, nxr, nyr, nzr,             &
   !$acc            ul, ur, vl, vr, wl, wr, pbl, pbr, clam, one_eta, pU_L, pU_R,            &
   !$acc            pbpert_edge, c_minus, c_plus, qbl, qbr,                                 &
   !$acc            ppl, ppr, upl, upr, vpl, vpr, wpl, wpr, pprime_lk, pprime_rk,           &
   !$acc            flux_pb, flux_u, flux_v, flux_w, fxl, fxr, oe2_Hql, oe2_Hqr, H_bcl_qlr, &
   !$acc            flux_edge_x, flux_edge_y, flux_edge_z,                                  &
   !$acc            Qu_ql, Qu_qr, Qv_ql, Qv_qr, Qw_ql, Qw_qr)                               &
   !$acc    firstprivate(nlayers_f, nvarb_f, bcl_stride_f, pb_idx_f, has_w, g_l,            &
   !$acc                 dry_cutoff_f, wst_l)
   do kk = 1, nboun_valid
      do iquad = 1, ngl_f

         iface = ref%face_pack_list(kk)
         el    = G%face(7, iface)

         nxl = mf%normal_vector(1,iquad,1,iface)
         nyl = mf%normal_vector(2,iquad,1,iface)
         nzl = 0.0
         if (has_w) nzl = mf%normal_vector(3,iquad,1,iface)
         nxr = -nxl; nyr = -nyl; nzr = -nzl

         ! Collocation: the face node itself, this rank (left) and the neighbour (right)
         qbl = 0.0; qbr = 0.0
         !$acc loop seq
         do iv = 1, nvarb_f
            qbl(iv) = ref%q_send(iv,iquad,kk)
            qbr(iv) = ref%q_recv(iv,iquad,kk)
         end do
         pbl = ref%q_send(pb_idx_f,iquad,kk)
         pbr = ref%q_recv(pb_idx_f,iquad,kk)

         pU_L = nxl*qbl(3) + nyl*qbl(4) + nzl*qbl(5)
         pU_R = nxr*qbr(3) + nyr*qbr(4) + nzr*qbr(5)

         c_minus = sqrt(init%alpha_mlswe(nlayers_f) * pbr)
         c_plus  = sqrt(init%alpha_mlswe(nlayers_f) * pbl)
         clam    = max(c_minus, c_plus)

         pbpert_edge = 0.5*(qbl(2) + qbr(2)) + (0.5/clam)*(pU_L + pU_R)
         one_eta     = 1.0 + (pbpert_edge/pbl)

         ul = qbl(3)/qbl(1); ur = qbr(3)/qbr(1)
         vl = qbl(4)/qbl(1); vr = qbr(4)/qbr(1)
         wl = 0.0; wr = 0.0
         if (has_w) then
            wl = qbl(5)/qbl(1); wr = qbr(5)/qbr(1)
         end if

         Qu_ql(1) = ul*qbl(3); Qu_ql(2) = vl*qbl(3); Qu_ql(3) = 0.0
         Qu_qr(1) = ur*qbr(3); Qu_qr(2) = vr*qbr(3); Qu_qr(3) = 0.0
         Qv_ql(1) = ul*qbl(4); Qv_ql(2) = vl*qbl(4); Qv_ql(3) = 0.0
         Qv_qr(1) = ur*qbr(4); Qv_qr(2) = vr*qbr(4); Qv_qr(3) = 0.0
         Qw_ql = 0.0; Qw_qr = 0.0
         if (has_w) then
            Qu_ql(3) = wl*qbl(3); Qu_qr(3) = wr*qbr(3)
            Qv_ql(3) = wl*qbl(4); Qv_qr(3) = wr*qbr(4)
            Qw_ql(1) = ul*qbl(5); Qw_ql(2) = vl*qbl(5); Qw_ql(3) = wl*qbl(5)
            Qw_qr(1) = ur*qbr(5); Qw_qr(2) = vr*qbr(5); Qw_qr(3) = wr*qbr(5)
         end if

         pprime_lk = 0.0; pprime_rk = 0.0
         H_bcl_qlr = 0.0

         !$acc loop seq
         do k = 1, nlayers_f
            ii = pb_idx_f + (k-1)*bcl_stride_f
            ppl = ref%q_send(ii+1,iquad,kk)
            ppr = ref%q_recv(ii+1,iquad,kk)
            upl = ref%q_send(ii+2,iquad,kk)
            upr = ref%q_recv(ii+2,iquad,kk)
            vpl = ref%q_send(ii+3,iquad,kk)
            vpr = ref%q_recv(ii+3,iquad,kk)
            wpl = 0.0; wpr = 0.0
            if (has_w) then
               wpl = ref%q_send(ii+4,iquad,kk)
               wpr = ref%q_recv(ii+4,iquad,kk)
            end if

            ! Dry-cell protection, matching the volume/face BCL flux routines: the
            ! thickness is clamped to a floor (the interface pressures stay continuous
            ! across dry_cutoff) and the velocity is zeroed, which removes this layer's
            ! Qu/Qv/Qw flux terms.
            if (ppl < (g_l/init%alpha_mlswe(k)) * dry_cutoff_f) then
               ppl = (g_l/init%alpha_mlswe(k)) * dry_cutoff_f
               upl = 0.0;  vpl = 0.0;  wpl = 0.0
            end if
            if (ppr < (g_l/init%alpha_mlswe(k)) * dry_cutoff_f) then
               ppr = (g_l/init%alpha_mlswe(k)) * dry_cutoff_f
               upr = 0.0;  vpr = 0.0;  wpr = 0.0
            end if

            Qu_ql(1) = Qu_ql(1) + upl*(upl*(one_eta*ppl))
            Qu_qr(1) = Qu_qr(1) + upr*(upr*(one_eta*ppr))
            Qu_ql(2) = Qu_ql(2) + vpl*(upl*(one_eta*ppl))
            Qu_qr(2) = Qu_qr(2) + vpr*(upr*(one_eta*ppr))
            Qv_ql(1) = Qv_ql(1) + upl*(vpl*(one_eta*ppl))
            Qv_qr(1) = Qv_qr(1) + upr*(vpr*(one_eta*ppr))
            Qv_ql(2) = Qv_ql(2) + vpl*(vpl*(one_eta*ppl))
            Qv_qr(2) = Qv_qr(2) + vpr*(vpr*(one_eta*ppr))

            if (has_w) then
               Qu_ql(3) = Qu_ql(3) + wpl*(upl*(one_eta*ppl));  Qu_qr(3) = Qu_qr(3) + wpr*(upr*(one_eta*ppr))
               Qv_ql(3) = Qv_ql(3) + wpl*(vpl*(one_eta*ppl));  Qv_qr(3) = Qv_qr(3) + wpr*(vpr*(one_eta*ppr))
               Qw_ql(1) = Qw_ql(1) + upl*(wpl*(one_eta*ppl));  Qw_qr(1) = Qw_qr(1) + upr*(wpr*(one_eta*ppr))
               Qw_ql(2) = Qw_ql(2) + vpl*(wpl*(one_eta*ppl));  Qw_qr(2) = Qw_qr(2) + vpr*(wpr*(one_eta*ppr))
               Qw_ql(3) = Qw_ql(3) + wpl*(wpl*(one_eta*ppl));  Qw_qr(3) = Qw_qr(3) + wpr*(wpr*(one_eta*ppr))
            end if

            pprime_lk = pprime_lk + ppl
            pprime_rk = pprime_rk + ppr

            H_bcl_qlr = H_bcl_qlr + btp%bcl_dalpha(k) * pprime_lk * pprime_rk
            btp%bcl_pprime_face_df(k,1,iquad,iface) = pprime_lk
            btp%bcl_pprime_face_df(k,2,iquad,iface) = pprime_rk
         end do

         flux_edge_x  = 0.5*(qbl(3)+qbr(3)) + (0.5*clam)*(nxl*qbl(2)+nxr*qbr(2))
         flux_edge_y  = 0.5*(qbl(4)+qbr(4)) + (0.5*clam)*(nyl*qbl(2)+nyr*qbr(2))
         flux_edge_z  = 0.0
         if (has_w) flux_edge_z = 0.5*(qbl(5)+qbr(5)) + (0.5*clam)*(nzl*qbl(2)+nzr*qbr(2))
         flux_pb = nxl*flux_edge_x + nyl*flux_edge_y + nzl*flux_edge_z

         ! Pressure part of the face flux in split form
         oe2_Hql = (one_eta**2) * H_bcl_qlr
         oe2_Hqr = (one_eta**2) * H_bcl_qlr

         ! Time-average accumulators (one per face node)
         btp%btp_mass_flux_face_ave_df(1,iquad,iface) = btp%btp_mass_flux_face_ave_df(1,iquad,iface) + wst_l*(flux_edge_x)
         btp%btp_mass_flux_face_ave_df(2,iquad,iface) = btp%btp_mass_flux_face_ave_df(2,iquad,iface) + wst_l*(flux_edge_y)
         btp%Qu_face_ave_df(1,iquad,iface)            = btp%Qu_face_ave_df(1,iquad,iface) + wst_l*(0.5*(Qu_ql(1)+Qu_qr(1)))
         btp%Qu_face_ave_df(2,iquad,iface)            = btp%Qu_face_ave_df(2,iquad,iface) + wst_l*(0.5*(Qu_ql(2)+Qu_qr(2)))
         btp%Qv_face_ave_df(1,iquad,iface)            = btp%Qv_face_ave_df(1,iquad,iface) + wst_l*(0.5*(Qv_ql(1)+Qv_qr(1)))
         btp%Qv_face_ave_df(2,iquad,iface)            = btp%Qv_face_ave_df(2,iquad,iface) + wst_l*(0.5*(Qv_ql(2)+Qv_qr(2)))
         btp%ope_face_ave_df(1,iquad,iface) = btp%ope_face_ave_df(1,iquad,iface) + wst_l*((1.0 + qbl(2)/pbl))
         btp%ope_face_ave_df(2,iquad,iface) = btp%ope_face_ave_df(2,iquad,iface) + wst_l*((1.0 + qbr(2)/pbr))
         btp%ope2_face_ave_df(1,iquad,iface) = btp%ope2_face_ave_df(1,iquad,iface) + wst_l*((1.0+qbl(2)/pbl)**2)
         btp%ope2_face_ave_df(2,iquad,iface) = btp%ope2_face_ave_df(2,iquad,iface) + wst_l*((1.0+qbr(2)/pbr)**2)
         btp%one_plus_eta_edge_2_ave_df(iquad,iface) = btp%one_plus_eta_edge_2_ave_df(iquad,iface) + wst_l*(one_eta**2)
         btp%uvb_face_ave_df(1,1,iquad,iface) = btp%uvb_face_ave_df(1,1,iquad,iface) + wst_l*(ul)
         btp%uvb_face_ave_df(1,2,iquad,iface) = btp%uvb_face_ave_df(1,2,iquad,iface) + wst_l*(ur)
         btp%uvb_face_ave_df(2,1,iquad,iface) = btp%uvb_face_ave_df(2,1,iquad,iface) + wst_l*(vl)
         btp%uvb_face_ave_df(2,2,iquad,iface) = btp%uvb_face_ave_df(2,2,iquad,iface) + wst_l*(vr)
         if (has_w) then
            btp%btp_mass_flux_face_ave_df(3,iquad,iface) = btp%btp_mass_flux_face_ave_df(3,iquad,iface) + wst_l*(flux_edge_z)
            btp%Qu_face_ave_df(3,iquad,iface) = btp%Qu_face_ave_df(3,iquad,iface) + wst_l*(0.5*(Qu_ql(3)+Qu_qr(3)))
            btp%Qv_face_ave_df(3,iquad,iface) = btp%Qv_face_ave_df(3,iquad,iface) + wst_l*(0.5*(Qv_ql(3)+Qv_qr(3)))
            btp%Qw_face_ave_df(1,iquad,iface) = btp%Qw_face_ave_df(1,iquad,iface) + wst_l*(0.5*(Qw_ql(1)+Qw_qr(1)))
            btp%Qw_face_ave_df(2,iquad,iface) = btp%Qw_face_ave_df(2,iquad,iface) + wst_l*(0.5*(Qw_ql(2)+Qw_qr(2)))
            btp%Qw_face_ave_df(3,iquad,iface) = btp%Qw_face_ave_df(3,iquad,iface) + wst_l*(0.5*(Qw_ql(3)+Qw_qr(3)))
            btp%uvb_face_ave_df(3,1,iquad,iface) = btp%uvb_face_ave_df(3,1,iquad,iface) + wst_l*(wl)
            btp%uvb_face_ave_df(3,2,iquad,iface) = btp%uvb_face_ave_df(3,2,iquad,iface) + wst_l*(wr)
         end if

         Qu_ql(1) = Qu_ql(1) + oe2_Hql
         Qu_qr(1) = Qu_qr(1) + oe2_Hqr
         fxl = nxl*Qu_ql(1) + nyl*Qu_ql(2) + nzl*Qu_ql(3)
         fxr = nxr*Qu_qr(1) + nyr*Qu_qr(2) + nzr*Qu_qr(3)
         flux_u = 0.5*(fxl-fxr) - (0.25*clam)*(qbr(3)-qbl(3))

         Qv_ql(2) = Qv_ql(2) + oe2_Hql
         Qv_qr(2) = Qv_qr(2) + oe2_Hqr
         fxl = nxl*Qv_ql(1) + nyl*Qv_ql(2) + nzl*Qv_ql(3)
         fxr = nxr*Qv_qr(1) + nyr*Qv_qr(2) + nzr*Qv_qr(3)
         flux_v = 0.5*(fxl-fxr) - (0.25*clam)*(qbr(4)-qbl(4))

         flux_w = 0.0
         if (has_w) then
            Qw_ql(3) = Qw_ql(3) + oe2_Hql
            Qw_qr(3) = Qw_qr(3) + oe2_Hqr
            fxl = nxl*Qw_ql(1) + nyl*Qw_ql(2) + nzl*Qw_ql(3)
            fxr = nxr*Qw_qr(1) + nyr*Qw_qr(2) + nzr*Qw_qr(3)
            flux_w = 0.5*(fxl-fxr) - (0.25*clam)*(qbr(5)-qbl(5))
         end if

         ! Store RHS: left side
         wq = mf%jac_face(iquad,1,iface)
         I  = G%intma(mf%imapl(1,iquad,1,iface), mf%imapl(2,iquad,1,iface), mf%imapl(3,iquad,1,iface), el)
         !$acc atomic update
         rhs(1,I) = rhs(1,I) - wq*flux_pb
         !$acc atomic update
         rhs(2,I) = rhs(2,I) - wq*flux_u
         !$acc atomic update
         rhs(3,I) = rhs(3,I) - wq*flux_v
         if (has_w) then
            !$acc atomic update
            rhs(4,I) = rhs(4,I) - wq*flux_w
         end if

      end do
   end do
   !$acc end parallel loop

end subroutine btp_colloc_nbhs_faces

end module mod_btp_colloc
