! =================================================================================================
! This module contains the routines for the barotropic and baroclinic horizontal viscosity terms
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: July 24, 2023
! =================================================================================================

module mod_laplacian_quad

   implicit none

   public :: bcl_create_laplacian, btp_create_laplacian, compute_bcl_lap_z

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

      real, intent(out) :: rhs_btp_visc(inp%nvar_btp-2,G%npoin)
      real, dimension(inp%nvar_btp,G%npoin), intent(in) :: qb_df

      real, dimension(inp%nvar_btp-2,G%npoin) :: Uk
      real, dimension(inp%ngraduvw_var,G%npoin) :: graduv
      integer :: I, Iq, ip
      real    :: dhdx, dhdy, dhdz
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

      ! Accumulate time-average on device.
      !$acc parallel loop present(btp, graduv) firstprivate(has_w)
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
      type(bcl_CS),           intent(in)    :: bcl
      type(mref),             intent(inout) :: ref
      type(mpi_communicator), intent(inout) :: mpic
      type(metrics),          intent(in)    :: mt
      type(tensor_CS),        intent(in)    :: tsp

      real, intent(out) :: rhs_lap(inp%nvar_bcl-1,G%npoin,inp%nlayers)

      integer :: k
      logical :: has_w

      has_w = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

      ! Pre-comm GPU pack reads dpp_graduvw/dpprime_visc directly from device.
      ! rhs_lap is created on device; all kernels run GPU-to-GPU.
      !$acc data create(rhs_lap)
      call bcl_lap_create_precommunicator(G, inp, b, mf, par, ref, mpic, bcl%dpp_graduvw, bcl%dpprime_visc, bcl%qprime_df)
      call bcl_compute_laplacian(G, inp, b, btp, bcl, tsp, rhs_lap)
      call bcl_create_rhs_laplacian_flux(G, inp, b, mf, btp, bcl, tsp, rhs_lap)
      ! CPU mpi_waitall inside post-comm, followed by GPU unpack + GPU face scatter.
      call bcl_create_rhs_lap_postcommunicator_df(G, inp, b, mf, par, btp, ref, mpic, tsp, rhs_lap)

      ! GPU: apply viscous mass-inverse scaling while rhs_lap is still on device.
      !$acc kernels present(rhs_lap, mt%massinv) firstprivate(has_w)
      do k = 1, inp%nlayers
         rhs_lap(1,:,k) = inp%visc_mlswe*mt%massinv(:)*rhs_lap(1,:,k)
         rhs_lap(2,:,k) = inp%visc_mlswe*mt%massinv(:)*rhs_lap(2,:,k)
         if (has_w) rhs_lap(3,:,k) = inp%visc_mlswe*mt%massinv(:)*rhs_lap(3,:,k)
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
      !$acc              btp%pbprime_visc, btp%btp_dpp_graduvw, grad_dpuvp,           &
      !$acc              rhs_btp_visc)

      !$acc kernels present(rhs_btp_visc)
      rhs_btp_visc = 0.0
      !$acc end kernels

      !$acc parallel loop gang                                     &
      !$acc   private(qq, wq, I, ip, dhdx, dhdy, dhdz)              &
      !$acc   firstprivate(npoin_l, npts_l, has_w)
      do Iq = 1, npoin_l

         wq    = tsp%wjac_df(Iq)
         qq(1) = btp%pbprime_visc(Iq)*grad_dpuvp(1,Iq) + btp%btp_dpp_graduvw(1,Iq)
         qq(2) = btp%pbprime_visc(Iq)*grad_dpuvp(2,Iq) + btp%btp_dpp_graduvw(2,Iq)
         qq(3) = btp%pbprime_visc(Iq)*grad_dpuvp(3,Iq) + btp%btp_dpp_graduvw(3,Iq)
         qq(4) = btp%pbprime_visc(Iq)*grad_dpuvp(4,Iq) + btp%btp_dpp_graduvw(4,Iq)
         if (has_w) then
            qq(5) = btp%pbprime_visc(Iq)*grad_dpuvp(5,Iq) + btp%btp_dpp_graduvw(5,Iq)
            qq(6) = btp%pbprime_visc(Iq)*grad_dpuvp(6,Iq) + btp%btp_dpp_graduvw(6,Iq)
            qq(7) = btp%pbprime_visc(Iq)*grad_dpuvp(7,Iq) + btp%btp_dpp_graduvw(7,Iq)
            qq(8) = btp%pbprime_visc(Iq)*grad_dpuvp(8,Iq) + btp%btp_dpp_graduvw(8,Iq)
            qq(9) = btp%pbprime_visc(Iq)*grad_dpuvp(9,Iq) + btp%btp_dpp_graduvw(9,Iq)
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
      !$acc              btp%graduvb_ave, bcl%dpprime_visc, bcl%dpp_graduvw, lap_q)

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
               qq(1) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(1,Iq) + bcl%dpp_graduvw(1,Iq,k)
               qq(2) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(2,Iq) + bcl%dpp_graduvw(2,Iq,k)
               qq(3) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(3,Iq) + bcl%dpp_graduvw(3,Iq,k)
               qq(4) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(4,Iq) + bcl%dpp_graduvw(4,Iq,k)
               if (has_w) then
                  qq(5) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(5,Iq) + bcl%dpp_graduvw(5,Iq,k)
                  qq(6) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(6,Iq) + bcl%dpp_graduvw(6,Iq,k)
                  qq(7) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(7,Iq) + bcl%dpp_graduvw(7,Iq,k)
                  qq(8) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(8,Iq) + bcl%dpp_graduvw(8,Iq,k)
                  qq(9) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(9,Iq) + bcl%dpp_graduvw(9,Iq,k)
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
      !$acc               btp%btp_dpp_graduvw, btp%pbprime_visc,               &
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
               btp_ql_s(ivar) = btp%btp_dpp_graduvw(ivar,ip)
            end do
            btp_ql_s(nw+1) = btp%pbprime_visc(ip)

            if (ier > 0) then

               ip = I_r(iquad)

               !$acc loop seq
               do ivar = 1, nw
                  qr_s(ivar)     = gradq(ivar,ip)
                  btp_qr_s(ivar) = btp%btp_dpp_graduvw(ivar,ip)
               end do
               btp_qr_s(nw+1) = btp%pbprime_visc(ip)

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
               pb_avg = 0.5 * (btp%pbprime_visc(I_l(iquad)) + btp%pbprime_visc(I_r(iquad)))
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
      real :: pconst, sigma, dp_avg_pen, du_pen, dv_pen, dw_pen

      beta  = 0.5
      alpha = 1.0 - beta
      iflux = 0.0

      ngl_f     = b%ngl
      nlayers_f = inp%nlayers
      nface_f   = G%nface
      nw        = inp%ngraduvw_var
      has_w     = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

      ! Shahbazi 2D SIPG constant: pconst = (p+1)(p+2)/2, p = nopx = ngl_f - 1.
      ! For ngl_f=5 (nopx=4): pconst = 15.  sigma = pconst * jac_face / jac_vol
      ! gives correct 1/h scaling for the SIP penalty.
      ! inp%SIPG_constant scales the penalty: 0 = LDG (no penalty), 1 = standard Shahbazi SIP.
      pconst = real((ngl_f+1)*(ngl_f+2)) / 2.0 * real(inp%SIPG_constant)

      !$acc data present(G%face, G%face_type, G%intma,                          &
      !$acc              mf%imapl, mf%imapr, mf%normal_vector, mf%jac_face,     &
      !$acc              b%psi,                                                   &
      !$acc              btp%graduvb_face_ave,                                   &
      !$acc              bcl%dpp_graduvw, bcl%dpprime_visc, bcl%qprime_df,       &
      !$acc              tsp%wjac_df,                                             &
      !$acc              rhs)

      !$acc parallel loop gang                                                    &
      !$acc   private(ql_cur, qr_cur, flux_uv_visc_face, qul, qur, qvl, qvr, qwl, qwr, &
      !$acc           qu_mean, qv_mean, qw_mean,                                     &
      !$acc           nx, ny, nz, un, wq, flux_qu, flux_qv, flux_qw, hi,            &
      !$acc           iel, ier, ip, il, jl, kl, ir, jr, kr, iquad, k, i, ivar,  &
      !$acc           ip_face_L, ip_face_R, sigma, dp_avg_pen, du_pen, dv_pen, dw_pen) &
      !$acc   firstprivate(alpha, beta, iflux, ngl_f, nlayers_f, nface_f, nw, has_w, pconst)
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
                  ql_cur(ivar) = bcl%dpp_graduvw(ivar,ip,k)
               end do
               ql_cur(nw+1) = bcl%dpprime_visc(ip,k)

               if (ier > 0) then
                  ir = mf%imapr(1,iquad,1,iface)
                  jr = mf%imapr(2,iquad,1,iface)
                  kr = mf%imapr(3,iquad,1,iface)
                  ip = G%intma(ir,jr,kr,ier)
                  ip_face_R = ip
                  !$acc loop seq
                  do ivar = 1, nw
                     qr_cur(ivar) = bcl%dpp_graduvw(ivar,ip,k)
                  end do
                  qr_cur(nw+1) = bcl%dpprime_visc(ip,k)
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

               ! SIP interior penalty: stabilises the variable-coefficient DG Laplacian.
               ! For constant dp_k the central-flux term has a zero eigenvalue for the
               ! 2-Delta-x checkerboard mode; for variable dp_k it becomes POSITIVE
               ! (anti-diffusive).  The penalty -eta*dp_avg/h*[u'_k] restores coercivity.
               ! Applied only for interior faces (ip_face_R /= ip_face_L after the above).
               ! Boundary faces retain ip_face_R = ip_face_L, giving zero jump and no effect.
               sigma = pconst * wq / min(tsp%wjac_df(ip_face_L), tsp%wjac_df(ip_face_R))
               dp_avg_pen = 0.5 * (ql_cur(nw+1) + qr_cur(nw+1))
               du_pen = bcl%qprime_df(2,ip_face_L,k) - bcl%qprime_df(2,ip_face_R,k)
               dv_pen = bcl%qprime_df(3,ip_face_L,k) - bcl%qprime_df(3,ip_face_R,k)
               flux_qu = flux_qu - sigma * dp_avg_pen * du_pen
               flux_qv = flux_qv - sigma * dp_avg_pen * dv_pen
               if (has_w) then
                  dw_pen = bcl%qprime_df(4,ip_face_L,k) - bcl%qprime_df(4,ip_face_R,k)
                  flux_qw = flux_qw - sigma * dp_avg_pen * dw_pen
               end if

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
   !  Computes the element-local Laplacian of the layer-interface height z_k
   !  at every DOF node for internal interfaces k = 2..nlayers.  The result
   !  is stored in bcl%lap_z_df(I, k) and later used by the Chen (2025) APE
   !  stabilization in create_rhs_dynamics_volume_bcl_sphere.
   !
   !  Method: two sequential element-local gradient passes with the
   !  DG derivative matrices dpsidx_df / dpsidy_df / dpsidz_df:
   !    1. gz(d, ip0) = Σ_ip  dpsid_df(ip, J_ip0) × z_k(J_ip)   = ∇z_k|_{J_ip0}
   !    2. lap_z(ip0) = Σ_ip  dpsid_df(ip, J_ip0) × gz(d, ip)   = ∇²z_k|_{J_ip0}
   !  Both sums run over the npts canonical nodes of the element.
   !
   !  Because DG nodes are element-local (no shared DOFs), the scatter to the
   !  global lap_z_df array needs no atomics.
   !
   !  CPU-only (no OpenACC): call !$acc update device(bcl%lap_z_df) afterwards.
   ! ===========================================================================
   subroutine compute_bcl_lap_z(G, inp, b, init, btp, tsp, qprime_df, lap_z_df)

      use mod_grid,      only: grid
      use mod_input,     only: input
      use mod_basis,     only: basis
      use mod_initial,   only: initial
      use mod_variables, only: btp_CS
      use mod_tensor,    only: tensor_CS
      use mod_constants, only: gravity

      implicit none

      type(grid),      intent(in)  :: G
      type(input),     intent(in)  :: inp
      type(basis),     intent(in)  :: b
      type(initial),   intent(in)  :: init
      type(btp_CS),    intent(in)  :: btp
      type(tensor_CS), intent(in)  :: tsp
      real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: qprime_df
      real, dimension(G%npoin, inp%nlayers+1),              intent(out) :: lap_z_df

      integer :: ie, ip0, ip, k
      integer :: J_ip0, J_ip
      integer :: npts_l, nelem_l, nlayers_l
      real    :: z_loc(b%npts, inp%nlayers+1)
      real    :: gz_x(b%npts), gz_y(b%npts), gz_z(b%npts)
      real    :: lap_z_loc(b%npts)
      real    :: dhdx, dhdy, dhdz
      logical :: has_w

      has_w     = (inp%nvar_bcl == 4)
      npts_l    = b%npts
      nelem_l   = G%nelem
      nlayers_l = inp%nlayers

      lap_z_df(:,:) = 0.0

      do ie = 1, nelem_l

         ! Step 1: gather z_k at all local nodes for k = 2..nlayers
         do ip = 1, npts_l
            J_ip = tsp%index_df_elt(ip, ie)
            z_loc(ip, nlayers_l+1) = init%zbot_df(J_ip)
            do k = nlayers_l, 1, -1
               z_loc(ip, k) = z_loc(ip, k+1) + &
                  (init%alpha_mlswe(k)/gravity) * sqrt(btp%ope2_ave_df(J_ip)) * qprime_df(1, J_ip, k)
            end do
         end do

         ! Steps 2–3: for each internal interface k, compute the element-local
         ! Laplacian at each node ip0 via two gradient applications.
         do k = 2, nlayers_l

            ! Step 2: gradient of z_k at each local node ip0
            do ip0 = 1, npts_l
               J_ip0  = tsp%index_df_elt(ip0, ie)
               gz_x(ip0) = 0.0
               gz_y(ip0) = 0.0
               if (has_w) gz_z(ip0) = 0.0
               do ip = 1, npts_l
                  dhdx = tsp%dpsidx_df(ip, J_ip0)
                  dhdy = tsp%dpsidy_df(ip, J_ip0)
                  if (has_w) then
                     dhdx = dhdx + tsp%dpsidz_df_x(ip, J_ip0)
                     dhdy = dhdy + tsp%dpsidz_df_y(ip, J_ip0)
                     dhdz = tsp%dpsidz_df(ip, J_ip0) + tsp%dpsidz_df_z(ip, J_ip0)
                  end if
                  gz_x(ip0) = gz_x(ip0) + dhdx * z_loc(ip, k)
                  gz_y(ip0) = gz_y(ip0) + dhdy * z_loc(ip, k)
                  if (has_w) gz_z(ip0) = gz_z(ip0) + dhdz * z_loc(ip, k)
               end do
            end do

            ! Step 3: divergence of gradient at each local node ip0 = ∇²z_k
            do ip0 = 1, npts_l
               J_ip0      = tsp%index_df_elt(ip0, ie)
               lap_z_loc(ip0) = 0.0
               do ip = 1, npts_l
                  dhdx = tsp%dpsidx_df(ip, J_ip0)
                  dhdy = tsp%dpsidy_df(ip, J_ip0)
                  if (has_w) then
                     dhdx = dhdx + tsp%dpsidz_df_x(ip, J_ip0)
                     dhdy = dhdy + tsp%dpsidz_df_y(ip, J_ip0)
                     dhdz = tsp%dpsidz_df(ip, J_ip0) + tsp%dpsidz_df_z(ip, J_ip0)
                  end if
                  lap_z_loc(ip0) = lap_z_loc(ip0) + dhdx*gz_x(ip) + dhdy*gz_y(ip)
                  if (has_w) lap_z_loc(ip0) = lap_z_loc(ip0) + dhdz*gz_z(ip)
               end do
            end do

            ! Step 4: scatter (no atomics — each DOF belongs to exactly one element)
            do ip0 = 1, npts_l
               J_ip0 = tsp%index_df_elt(ip0, ie)
               lap_z_df(J_ip0, k) = lap_z_loc(ip0)
            end do

         end do ! k

      end do ! ie

   end subroutine compute_bcl_lap_z

end module mod_laplacian_quad
