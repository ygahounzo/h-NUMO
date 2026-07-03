! =================================================================================================
! This module contains the routines for the barotropic and baroclinic horizontal viscosity terms
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: July 24, 2023
! =================================================================================================

module mod_laplacian_quad

   implicit none

   public :: bcl_create_laplacian, btp_create_laplacian

contains

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

      real, intent(out) :: rhs_btp_visc(2,G%npoin)
      real, dimension(4,G%npoin), intent(in) :: qb_df

      real, dimension(2,G%npoin) :: Uk
      real, dimension(4,G%npoin) :: graduv
      integer :: I, Iq, ip
      real    :: dhdx, dhdy

      !$acc data create(Uk, graduv)

      ! Compute barotropic velocity — qb_df lives on device during the BTP loop.
      !$acc parallel loop present(Uk, qb_df)
      do I = 1, G%npoin
         Uk(1,I) = qb_df(3,I) / qb_df(1,I)
         Uk(2,I) = qb_df(4,I) / qb_df(1,I)
      end do
      !$acc end parallel loop

      ! Compute velocity gradient (inlined from compute_gradient_uv).
      !$acc kernels present(graduv)
      graduv = 0.0
      !$acc end kernels

      !$acc parallel loop present(graduv, Uk, tsp, b, G)
      do Iq = 1, G%npoin
         !$acc loop seq
         do ip = 1, b%npts
            I    = tsp%index_df(ip,Iq)
            dhdx = tsp%dpsidx_df(ip,Iq)
            dhdy = tsp%dpsidy_df(ip,Iq)
            graduv(1,Iq) = graduv(1,Iq) + dhdx*Uk(1,I)
            graduv(2,Iq) = graduv(2,Iq) + dhdy*Uk(1,I)
            graduv(3,Iq) = graduv(3,Iq) + dhdx*Uk(2,I)
            graduv(4,Iq) = graduv(4,Iq) + dhdy*Uk(2,I)
         end do
      end do
      !$acc end parallel loop

      ! Accumulate time-average on device.
      !$acc parallel loop present(btp, graduv)
      do I = 1, G%npoin
         btp%graduvb_ave(1,I) = btp%graduvb_ave(1,I) + graduv(1,I)
         btp%graduvb_ave(2,I) = btp%graduvb_ave(2,I) + graduv(2,I)
         btp%graduvb_ave(3,I) = btp%graduvb_ave(3,I) + graduv(3,I)
         btp%graduvb_ave(4,I) = btp%graduvb_ave(4,I) + graduv(4,I)
      end do
      !$acc end parallel loop

      ! MPI precommunicator reads graduv from host — download once before packing.
      !$acc update host(graduv)

      call btp_lap_create_precommunicator(G, b, mf, init, par, btp, ref, mpic, graduv, 4)

      ! rhs_btp_visc already zeroed on device by caller (create_rhs_btp).
      call btp_compute_laplacian_qp(G, b, btp, tsp, rhs_btp_visc, graduv)
      call create_rhs_laplacian_flux(G, b, mf, btp, rhs_btp_visc, graduv)

      call create_rhs_lap_postcommunicator_df(G, b, mf, par, btp, ref, mpic, rhs_btp_visc, 4)

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
      type(bcl_CS),           intent(in)    :: bcl
      type(mref),             intent(inout) :: ref
      type(mpi_communicator), intent(inout) :: mpic
      type(metrics),          intent(in)    :: mt
      type(tensor_CS),        intent(in)    :: tsp

      real, intent(out) :: rhs_lap(2,G%npoin,inp%nlayers)

      integer :: k

      ! Pre-comm GPU pack reads dpp_graduvw/dpprime_visc directly from device.
      ! rhs_lap is created on device; all kernels run GPU-to-GPU.
      ! Only the cartesian-equivalent 4 components cross MPI boundaries today —
      ! the sphere-only slots 5:9 of dpp_graduvw aren't yet consumed by the
      ! viscosity Laplacian (same rationale as btp_dpp_graduvw).
      !$acc data create(rhs_lap)
      call bcl_lap_create_precommunicator(G, inp, b, mf, par, ref, mpic, bcl%dpp_graduvw(1:4,:,:), bcl%dpprime_visc)
      call bcl_compute_laplacian(G, inp, b, btp, bcl, tsp, rhs_lap)
      call bcl_create_rhs_laplacian_flux(G, inp, b, mf, btp, bcl, rhs_lap)
      ! CPU mpi_waitall inside post-comm, followed by GPU unpack + GPU face scatter.
      call bcl_create_rhs_lap_postcommunicator_df(G, inp, b, mf, par, btp, ref, mpic, rhs_lap)

      ! GPU: apply viscous mass-inverse scaling while rhs_lap is still on device.
      !$acc kernels present(rhs_lap, mt%massinv)
      do k = 1, inp%nlayers
         rhs_lap(1,:,k) = inp%visc_mlswe*mt%massinv(:)*rhs_lap(1,:,k)
         rhs_lap(2,:,k) = inp%visc_mlswe*mt%massinv(:)*rhs_lap(2,:,k)
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

   subroutine btp_compute_laplacian_qp(G, b, btp, tsp, rhs_btp_visc, grad_dpuvp)

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

      integer :: Iq, I, ip
      real :: wq, qq(4)
      integer :: npoin_l, npts_l

      npoin_l = G%npoin
      npts_l  = b%npts

      !$acc data present(tsp%wjac_df, tsp%dpsidx_df, tsp%dpsidy_df, tsp%index_df,  &
      !$acc              btp%pbprime_visc, btp%btp_dpp_graduvw, grad_dpuvp,           &
      !$acc              rhs_btp_visc)

      !$acc kernels present(rhs_btp_visc)
      rhs_btp_visc = 0.0
      !$acc end kernels

      !$acc parallel loop gang                                     &
      !$acc   private(qq, wq, I, ip)                               &
      !$acc   firstprivate(npoin_l, npts_l)
      do Iq = 1, npoin_l

         wq    = tsp%wjac_df(Iq)
         qq(1) = btp%pbprime_visc(Iq)*grad_dpuvp(1,Iq) + btp%btp_dpp_graduvw(1,Iq)
         qq(2) = btp%pbprime_visc(Iq)*grad_dpuvp(2,Iq) + btp%btp_dpp_graduvw(2,Iq)
         qq(3) = btp%pbprime_visc(Iq)*grad_dpuvp(3,Iq) + btp%btp_dpp_graduvw(3,Iq)
         qq(4) = btp%pbprime_visc(Iq)*grad_dpuvp(4,Iq) + btp%btp_dpp_graduvw(4,Iq)

         !$acc loop seq
         do ip = 1, npts_l
            I = tsp%index_df(ip,Iq)
            !$acc atomic update
            rhs_btp_visc(1,I) = rhs_btp_visc(1,I) - wq*(tsp%dpsidx_df(ip,Iq)*qq(1) + tsp%dpsidy_df(ip,Iq)*qq(2))
            !$acc atomic update
            rhs_btp_visc(2,I) = rhs_btp_visc(2,I) - wq*(tsp%dpsidx_df(ip,Iq)*qq(3) + tsp%dpsidy_df(ip,Iq)*qq(4))
         end do
      end do
      !$acc end parallel loop

      !$acc end data

   end subroutine btp_compute_laplacian_qp

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

      real, intent(out) :: lap_q(2,G%npoin,inp%nlayers)

      integer :: ie, iq_local, Iq, Iq0, I, ip, k
      integer :: npts_l, nelem_l, nlayers_l
      real :: wq, qq(4)
      real :: lap_loc(2, b%npts, inp%nlayers)

      npts_l    = b%npts
      nelem_l   = G%nelem
      nlayers_l = inp%nlayers

      !$acc data present(tsp%wjac_df, tsp%dpsidx_df, tsp%dpsidy_df,        &
      !$acc              tsp%index_df_elt, tsp%index_df,                     &
      !$acc              btp%graduvb_ave, bcl%dpprime_visc, bcl%dpp_graduvw, lap_q)

      !$acc kernels present(lap_q)
      lap_q = 0.0
      !$acc end kernels

      !$acc parallel loop gang                                                    &
      !$acc   private(lap_loc, qq, wq, Iq, Iq0, I, ip, k, iq_local)            &
      !$acc   firstprivate(npts_l, nelem_l, nlayers_l)
      do ie = 1, nelem_l

         !$acc loop seq
         do k = 1, nlayers_l
            !$acc loop seq
            do ip = 1, npts_l
               lap_loc(1,ip,k) = 0.0
               lap_loc(2,ip,k) = 0.0
            end do
         end do

         !$acc loop seq
         do k = 1, nlayers_l
            !$acc loop seq
            do iq_local = 1, npts_l
               Iq = tsp%index_df_elt(iq_local, ie)
               wq = tsp%wjac_df(Iq)
               qq(1) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(1,Iq) + bcl%dpp_graduvw(1,Iq,k)
               qq(2) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(2,Iq) + bcl%dpp_graduvw(2,Iq,k)
               qq(3) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(3,Iq) + bcl%dpp_graduvw(3,Iq,k)
               qq(4) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(4,Iq) + bcl%dpp_graduvw(4,Iq,k)
               !$acc loop seq
               do ip = 1, npts_l
                  lap_loc(1,ip,k) = lap_loc(1,ip,k) - wq*(tsp%dpsidx_df(ip,Iq)*qq(1) + tsp%dpsidy_df(ip,Iq)*qq(2))
                  lap_loc(2,ip,k) = lap_loc(2,ip,k) - wq*(tsp%dpsidx_df(ip,Iq)*qq(3) + tsp%dpsidy_df(ip,Iq)*qq(4))
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

   subroutine create_rhs_laplacian_flux(G, b, mf, btp, rhs, gradq)
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
      !=========================================================================
      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_face,      only: face_CS
      use mod_variables, only: btp_CS

      implicit none

      type(grid),    intent(in)    :: G
      type(basis),   intent(in)    :: b
      type(face_CS), intent(in)    :: mf
      type(btp_CS),  intent(inout) :: btp

      real, intent(inout) :: rhs(2,G%npoin)
      real, intent(in)    :: gradq(4,G%npoin)

      real :: qu_mean(2), qv_mean(2)
      real :: flux_uv_visc_face(4,2)
      real :: qul(2), qur(2), qvl(2), qvr(2)
      real :: nx, ny, wq, un
      integer :: iface, i, il, jl, kl, ir, jr, kr
      integer :: iel, ier, ip, iquad, ivar
      real :: flux_qu, flux_qv, hi, alpha, beta, iflux
      real    :: ql_s(4), qr_s(4), btp_ql_s(5), btp_qr_s(5)
      integer, dimension(b%ngl) :: I_l, I_r
      integer :: ngl_f, nface_f

      beta  = 0.5
      alpha = 1.0 - beta
      iflux = 0.0

      ngl_f   = b%ngl
      nface_f = G%nface

      !$acc data present(G%face, G%face_type, G%intma,                       &
      !$acc               mf%imapl, mf%imapr, mf%normal_vector, mf%jac_face, &
      !$acc               b%psi,                                               &
      !$acc               btp%btp_dpp_graduvw, btp%pbprime_visc,               &
      !$acc               btp%graduvb_face_ave,                                &
      !$acc               rhs, gradq)

      !$acc parallel loop gang                                                &
      !$acc   private(iel, ier, il, jl, kl, ir, jr, kr, ip,                  &
      !$acc           I_l, I_r,                                               &
      !$acc           nx, ny, un, wq,                                         &
      !$acc           ql_s, qr_s, btp_ql_s, btp_qr_s,                       &
      !$acc           flux_uv_visc_face,                                      &
      !$acc           qul, qur, qvl, qvr, qu_mean, qv_mean,                  &
      !$acc           flux_qu, flux_qv, hi, ivar, i, iquad)                  &
      !$acc   firstprivate(alpha, beta, iflux, ngl_f, nface_f)
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

            ip = I_l(iquad)

            ql_s(1) = gradq(1,ip)
            ql_s(2) = gradq(2,ip)
            ql_s(3) = gradq(3,ip)
            ql_s(4) = gradq(4,ip)

            btp_ql_s(1) = btp%btp_dpp_graduvw(1,ip)
            btp_ql_s(2) = btp%btp_dpp_graduvw(2,ip)
            btp_ql_s(3) = btp%btp_dpp_graduvw(3,ip)
            btp_ql_s(4) = btp%btp_dpp_graduvw(4,ip)
            btp_ql_s(5) = btp%pbprime_visc(ip)

            if (ier > 0) then

               ip = I_r(iquad)

               qr_s(1) = gradq(1,ip)
               qr_s(2) = gradq(2,ip)
               qr_s(3) = gradq(3,ip)
               qr_s(4) = gradq(4,ip)

               btp_qr_s(1) = btp%btp_dpp_graduvw(1,ip)
               btp_qr_s(2) = btp%btp_dpp_graduvw(2,ip)
               btp_qr_s(3) = btp%btp_dpp_graduvw(3,ip)
               btp_qr_s(4) = btp%btp_dpp_graduvw(4,ip)
               btp_qr_s(5) = btp%pbprime_visc(ip)

            else
               qr_s(1) = ql_s(1)
               qr_s(2) = ql_s(2)
               qr_s(3) = ql_s(3)
               qr_s(4) = ql_s(4)

               btp_qr_s(1) = btp_ql_s(1)
               btp_qr_s(2) = btp_ql_s(2)
               btp_qr_s(3) = btp_ql_s(3)
               btp_qr_s(4) = btp_ql_s(4)
               btp_qr_s(5) = btp_ql_s(5)

               if (ier == -4) then

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

            ! (iquad,iface) unique per gang×seq-step — no atomics needed.
            btp%graduvb_face_ave(1,1,iquad,iface) = btp%graduvb_face_ave(1,1,iquad,iface) + ql_s(1)
            btp%graduvb_face_ave(2,1,iquad,iface) = btp%graduvb_face_ave(2,1,iquad,iface) + ql_s(2)
            btp%graduvb_face_ave(3,1,iquad,iface) = btp%graduvb_face_ave(3,1,iquad,iface) + ql_s(3)
            btp%graduvb_face_ave(4,1,iquad,iface) = btp%graduvb_face_ave(4,1,iquad,iface) + ql_s(4)

            btp%graduvb_face_ave(1,2,iquad,iface) = btp%graduvb_face_ave(1,2,iquad,iface) + qr_s(1)
            btp%graduvb_face_ave(2,2,iquad,iface) = btp%graduvb_face_ave(2,2,iquad,iface) + qr_s(2)
            btp%graduvb_face_ave(3,2,iquad,iface) = btp%graduvb_face_ave(3,2,iquad,iface) + qr_s(3)
            btp%graduvb_face_ave(4,2,iquad,iface) = btp%graduvb_face_ave(4,2,iquad,iface) + qr_s(4)

            !$acc loop seq
            do ivar = 1, 4
               flux_uv_visc_face(ivar,1) = btp_ql_s(5)*ql_s(ivar) + btp_ql_s(ivar)
               flux_uv_visc_face(ivar,2) = btp_qr_s(5)*qr_s(ivar) + btp_qr_s(ivar)
            end do

            qul(1) = flux_uv_visc_face(1,1)
            qul(2) = flux_uv_visc_face(2,1)
            qvl(1) = flux_uv_visc_face(3,1)
            qvl(2) = flux_uv_visc_face(4,1)

            qur(1) = flux_uv_visc_face(1,2)
            qur(2) = flux_uv_visc_face(2,2)
            qvr(1) = flux_uv_visc_face(3,2)
            qvr(2) = flux_uv_visc_face(4,2)

            qu_mean(1) = alpha*qul(1) + beta*qur(1)
            qu_mean(2) = alpha*qul(2) + beta*qur(2)
            qv_mean(1) = alpha*qvl(1) + beta*qvr(1)
            qv_mean(2) = alpha*qvl(2) + beta*qvr(2)

            wq = mf%jac_face(iquad,1,iface)

            flux_qu = (qu_mean(1) - iflux*qul(1))*nx + (qu_mean(2) - iflux*qul(2))*ny
            flux_qv = (qv_mean(1) - iflux*qvl(1))*nx + (qv_mean(2) - iflux*qvl(2))*ny

            ! rhs scatter — atomics: multiple faces share boundary nodes.
            !$acc loop seq
            do i = 1, ngl_f
               hi = b%psi(i,iquad)
               ip = I_l(i)
               !$acc atomic update
               rhs(1,ip) = rhs(1,ip) + wq*hi*flux_qu
               !$acc atomic update
               rhs(2,ip) = rhs(2,ip) + wq*hi*flux_qv
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
               end do
            end if

         end do  ! iquad

      end do  ! iface
      !$acc end data

   end subroutine create_rhs_laplacian_flux

   subroutine bcl_create_rhs_laplacian_flux(G, inp, b, mf, btp, bcl, rhs)

      use mod_grid,      only: grid
      use mod_basis,     only: basis
      use mod_input,     only: input
      use mod_face,      only: face_CS
      use mod_variables, only: btp_CS, bcl_CS

      implicit none

      type(grid),    intent(in)    :: G
      type(basis),   intent(in)    :: b
      type(input),   intent(in)    :: inp
      type(face_CS), intent(in)    :: mf
      type(btp_CS),  intent(in)    :: btp
      type(bcl_CS),  intent(in)    :: bcl

      real, intent(inout) :: rhs(2,G%npoin,inp%nlayers)

      real, dimension(2) :: qu_mean, qv_mean
      real, dimension(4,2) :: flux_uv_visc_face
      real, dimension(2) :: qul, qur
      real, dimension(2) :: qvl, qvr

      real :: nx, ny, un
      real :: wq
      integer :: iface, i, il, jl, kl, ir, jr, kr
      integer :: iel, ier, ip
      integer :: iquad, ivar, k
      real :: flux_qu, flux_qv, hi, alpha, beta, iflux
      real :: ql_cur(5), qr_cur(5)
      integer :: ngl_f, nlayers_f, nface_f

      beta  = 0.5
      alpha = 1.0 - beta
      iflux = 0.0

      ngl_f     = b%ngl
      nlayers_f = inp%nlayers
      nface_f   = G%nface

      !$acc data present(G%face, G%face_type, G%intma,                          &
      !$acc              mf%imapl, mf%imapr, mf%normal_vector, mf%jac_face,     &
      !$acc              b%psi,                                                   &
      !$acc              btp%graduvb_face_ave,                                   &
      !$acc              bcl%dpp_graduvw, bcl%dpprime_visc,                       &
      !$acc              rhs)

      !$acc parallel loop gang                                                    &
      !$acc   private(ql_cur, qr_cur, flux_uv_visc_face, qul, qur, qvl, qvr,   &
      !$acc           qu_mean, qv_mean,                                           &
      !$acc           nx, ny, un, wq, flux_qu, flux_qv, hi,                      &
      !$acc           iel, ier, ip, il, jl, kl, ir, jr, kr, iquad, k, i, ivar)  &
      !$acc   firstprivate(alpha, beta, iflux, ngl_f, nlayers_f, nface_f)
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

               il = mf%imapl(1,iquad,1,iface)
               jl = mf%imapl(2,iquad,1,iface)
               kl = mf%imapl(3,iquad,1,iface)
               ip = G%intma(il,jl,kl,iel)

               ql_cur(1) = bcl%dpp_graduvw(1,ip,k)
               ql_cur(2) = bcl%dpp_graduvw(2,ip,k)
               ql_cur(3) = bcl%dpp_graduvw(3,ip,k)
               ql_cur(4) = bcl%dpp_graduvw(4,ip,k)
               ql_cur(5) = bcl%dpprime_visc(ip,k)

               if (ier > 0) then
                  ir = mf%imapr(1,iquad,1,iface)
                  jr = mf%imapr(2,iquad,1,iface)
                  kr = mf%imapr(3,iquad,1,iface)
                  ip = G%intma(ir,jr,kr,ier)
                  qr_cur(1) = bcl%dpp_graduvw(1,ip,k)
                  qr_cur(2) = bcl%dpp_graduvw(2,ip,k)
                  qr_cur(3) = bcl%dpp_graduvw(3,ip,k)
                  qr_cur(4) = bcl%dpp_graduvw(4,ip,k)
                  qr_cur(5) = bcl%dpprime_visc(ip,k)
               else
                  qr_cur = ql_cur
                  if (ier == -4) then
                     un = ql_cur(1)*nx + ql_cur(2)*ny
                     qr_cur(1) = ql_cur(1) - 2.0*un*nx
                     qr_cur(2) = ql_cur(2) - 2.0*un*ny
                     un = ql_cur(3)*nx + ql_cur(4)*ny
                     qr_cur(3) = ql_cur(3) - 2.0*un*nx
                     qr_cur(4) = ql_cur(4) - 2.0*un*ny
                  end if
               end if

               !$acc loop seq
               do ivar = 1, 4
                  flux_uv_visc_face(ivar,1) = ql_cur(5)*btp%graduvb_face_ave(ivar,1,iquad,iface) + ql_cur(ivar)
                  flux_uv_visc_face(ivar,2) = qr_cur(5)*btp%graduvb_face_ave(ivar,2,iquad,iface) + qr_cur(ivar)
               end do

               qul(1) = flux_uv_visc_face(1,1);  qul(2) = flux_uv_visc_face(2,1)
               qvl(1) = flux_uv_visc_face(3,1);  qvl(2) = flux_uv_visc_face(4,1)
               qur(1) = flux_uv_visc_face(1,2);  qur(2) = flux_uv_visc_face(2,2)
               qvr(1) = flux_uv_visc_face(3,2);  qvr(2) = flux_uv_visc_face(4,2)

               qu_mean(1) = alpha*qul(1) + beta*qur(1)
               qu_mean(2) = alpha*qul(2) + beta*qur(2)
               qv_mean(1) = alpha*qvl(1) + beta*qvr(1)
               qv_mean(2) = alpha*qvl(2) + beta*qvr(2)

               wq = mf%jac_face(iquad,1,iface)
               flux_qu = (qu_mean(1) - iflux*qul(1))*nx + (qu_mean(2) - iflux*qul(2))*ny
               flux_qv = (qv_mean(1) - iflux*qvl(1))*nx + (qv_mean(2) - iflux*qvl(2))*ny

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
                  end do
               end if

            end do !iquad
         end do !k
      end do !iface
      !$acc end parallel loop

      !$acc end data

   end subroutine bcl_create_rhs_laplacian_flux

end module mod_laplacian_quad
