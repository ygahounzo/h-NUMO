subroutine create_nbhs_face_df(G, inp, b, mf, par, btp, init, ref, rhs)

   use mod_grid,      only: grid
   use mod_input,     only: input
   use mod_basis,     only: basis
   use mod_face,      only: face_CS
   use mod_parallel,  only: parallel_CS
   use mod_variables, only: btp_CS
   use mod_initial,   only: initial
   use mod_ref,       only: mref

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

   integer :: kk, iface, iquad, el, il, jl, I, kl, n, k, ii, nboun_valid
   real :: wq, hi, nxl, nyl, nzl, nxr, nyr, nzr
   real :: ul, ur, vl, vr, wl, wr, pbl, pbr, clam, one_eta
   real :: pU_L, pU_R, pbpert_edge, c_minus, c_plus
   real :: qbl(5), qbr(5), flux(4,b%nq)
   real :: ppl, ppr, upl, upr, vpl, vpr, wpl, wpr
   real, dimension(inp%nlayers+1) :: pprime_l, pprime_r
   real :: H_bcl_ql, H_bcl_qr, H_bcl_q, flux_pb, flux_u, flux_v, flux_w, fxl, fxr
   real :: oe2_Hql, oe2_Hqr
   real :: flux_edge_x, flux_edge_y, flux_edge_z
   real, dimension(3) :: Qu_ql, Qu_qr, Qv_ql, Qv_qr, Qw_ql, Qw_qr
   integer :: nq_f, ngl_f, nlayers_f, nvarb_f, bcl_stride_f, pb_idx_f
   logical :: has_w

   nq_f       = b%nq
   ngl_f      = b%ngl
   nlayers_f  = inp%nlayers
   nboun_valid = ref%nboun_valid
   nvarb_f      = inp%nvar_btp
   bcl_stride_f = inp%nvar_bcl
   pb_idx_f     = inp%nvar_btp + 1  ! pbprime is packed right after the momentum block
   has_w        = (inp%nvar_btp == 5) ! w (vertical momentum) is only carried on sphere_hex

   ! One gang per valid MPI face; iquad is vectorised within each face.
   ! rhs scatter uses atomic: different iquad lanes and different faces can
   ! write to the same volume node through imapl.
   !$acc parallel loop gang private(iface, el, flux) &
   !$acc    present(G%face, G%intma, mf%normal_vector_q, mf%jac_faceq, mf%imapl, &
   !$acc            b%psiq, init%alpha_mlswe, ref%q_send, ref%q_recv, rhs, &
   !$acc            btp%btp_mass_flux_face_ave, btp%H_face_ave, &
   !$acc            btp%Qu_face_ave, btp%Qv_face_ave, btp%Qw_face_ave, btp%ope_face_ave, &
   !$acc            btp%ope2_face_ave, btp%one_plus_eta_edge_2_ave, btp%uvb_face_ave, &
   !$acc            ref%face_pack_list) &
   !$acc    firstprivate(nq_f, ngl_f, nlayers_f, nboun_valid, nvarb_f, bcl_stride_f, pb_idx_f, has_w)
   do kk = 1, nboun_valid
      iface = ref%face_pack_list(kk)
      el    = G%face(7, iface)

      ! Compute fluxes and accumulate face-averaged quantities
      !$acc loop vector &
      !$acc    private(qbl, qbr, pprime_l, pprime_r, Qu_ql, Qu_qr, Qv_ql, Qv_qr, Qw_ql, Qw_qr)
      do iquad = 1, nq_f

         nxl = mf%normal_vector_q(1,iquad,1,iface)
         nyl = mf%normal_vector_q(2,iquad,1,iface)
         nzl = 0.0
         if (has_w) nzl = mf%normal_vector_q(3,iquad,1,iface)
         nxr = -nxl; nyr = -nyl; nzr = -nzl

         qbl = 0.0; qbr = 0.0; pbl = 0.0; pbr = 0.0
         !$acc loop seq
         do n = 1, ngl_f
            hi           = b%psiq(n,iquad)
            qbl(1:nvarb_f) = qbl(1:nvarb_f) + hi*ref%q_send(1:nvarb_f,n,kk)
            pbl          = pbl      + hi*ref%q_send(pb_idx_f,n,kk)
            qbr(1:nvarb_f) = qbr(1:nvarb_f) + hi*ref%q_recv(1:nvarb_f,n,kk)
            pbr          = pbr      + hi*ref%q_recv(pb_idx_f,n,kk)
         end do

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

         pprime_l(:) = 0.0; pprime_r(:) = 0.0
         H_bcl_ql    = 0.0; H_bcl_qr   = 0.0; H_bcl_q = 0.0

         !$acc loop seq
         do k = 1, nlayers_f
            ppl = 0.0; ppr = 0.0; upl = 0.0; upr = 0.0; vpl = 0.0; vpr = 0.0
            wpl = 0.0; wpr = 0.0
            ii = pb_idx_f + (k-1)*bcl_stride_f
            !$acc loop seq
            do n = 1, ngl_f
               hi  = b%psiq(n,iquad)
               ppl = ppl + hi*ref%q_send(ii+1,n,kk)
               ppr = ppr + hi*ref%q_recv(ii+1,n,kk)
               upl = upl + hi*ref%q_send(ii+2,n,kk)
               upr = upr + hi*ref%q_recv(ii+2,n,kk)
               vpl = vpl + hi*ref%q_send(ii+3,n,kk)
               vpr = vpr + hi*ref%q_recv(ii+3,n,kk)
               if (has_w) then
                  wpl = wpl + hi*ref%q_send(ii+4,n,kk)
                  wpr = wpr + hi*ref%q_recv(ii+4,n,kk)
               end if
            end do

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

            pprime_l(k+1) = pprime_l(k) + ppl
            pprime_r(k+1) = pprime_r(k) + ppr

            H_bcl_ql = H_bcl_ql + 0.5*init%alpha_mlswe(k)*(pprime_l(k+1)**2 - pprime_l(k)**2)
            H_bcl_qr = H_bcl_qr + 0.5*init%alpha_mlswe(k)*(pprime_r(k+1)**2 - pprime_r(k)**2)
         end do

         flux_edge_x  = 0.5*(qbl(3)+qbr(3)) + (0.5*clam)*(nxl*qbl(2)+nxr*qbr(2))
         flux_edge_y  = 0.5*(qbl(4)+qbr(4)) + (0.5*clam)*(nyl*qbl(2)+nyr*qbr(2))
         flux_edge_z  = 0.0
         if (has_w) flux_edge_z = 0.5*(qbl(5)+qbr(5)) + (0.5*clam)*(nzl*qbl(2)+nzr*qbr(2))
         flux(1,iquad) = nxl*flux_edge_x + nyl*flux_edge_y + nzl*flux_edge_z

         H_bcl_q = (one_eta**2) * (0.5*(H_bcl_ql + H_bcl_qr))
         oe2_Hql = (one_eta**2) * H_bcl_ql
         oe2_Hqr = (one_eta**2) * H_bcl_qr

         btp%btp_mass_flux_face_ave(1,iquad,iface) = btp%btp_mass_flux_face_ave(1,iquad,iface) + flux_edge_x
         btp%btp_mass_flux_face_ave(2,iquad,iface) = btp%btp_mass_flux_face_ave(2,iquad,iface) + flux_edge_y
         btp%H_face_ave(iquad,iface)               = btp%H_face_ave(iquad,iface) + H_bcl_q
         btp%Qu_face_ave(1,iquad,iface)            = btp%Qu_face_ave(1,iquad,iface) + 0.5*(Qu_ql(1)+Qu_qr(1))
         btp%Qu_face_ave(2,iquad,iface)            = btp%Qu_face_ave(2,iquad,iface) + 0.5*(Qu_ql(2)+Qu_qr(2))
         btp%Qv_face_ave(1,iquad,iface)            = btp%Qv_face_ave(1,iquad,iface) + 0.5*(Qv_ql(1)+Qv_qr(1))
         btp%Qv_face_ave(2,iquad,iface)            = btp%Qv_face_ave(2,iquad,iface) + 0.5*(Qv_ql(2)+Qv_qr(2))
         btp%ope_face_ave(1,iquad,iface)           = btp%ope_face_ave(1,iquad,iface) + (1.0 + qbl(2)/pbl)
         btp%ope_face_ave(2,iquad,iface)           = btp%ope_face_ave(2,iquad,iface) + (1.0 + qbr(2)/pbr)
         btp%ope2_face_ave(1,iquad,iface)          = btp%ope2_face_ave(1,iquad,iface) + (1.0+qbl(2)/pbl)**2
         btp%ope2_face_ave(2,iquad,iface)          = btp%ope2_face_ave(2,iquad,iface) + (1.0+qbr(2)/pbr)**2
         btp%one_plus_eta_edge_2_ave(iquad,iface)  = btp%one_plus_eta_edge_2_ave(iquad,iface) + one_eta**2
         btp%uvb_face_ave(1,1,iquad,iface)         = btp%uvb_face_ave(1,1,iquad,iface) + ul
         btp%uvb_face_ave(1,2,iquad,iface)         = btp%uvb_face_ave(1,2,iquad,iface) + ur
         btp%uvb_face_ave(2,1,iquad,iface)         = btp%uvb_face_ave(2,1,iquad,iface) + vl
         btp%uvb_face_ave(2,2,iquad,iface)         = btp%uvb_face_ave(2,2,iquad,iface) + vr

         if (has_w) then
            btp%btp_mass_flux_face_ave(3,iquad,iface) = btp%btp_mass_flux_face_ave(3,iquad,iface) + flux_edge_z
            btp%Qu_face_ave(3,iquad,iface) = btp%Qu_face_ave(3,iquad,iface) + 0.5*(Qu_ql(3)+Qu_qr(3))
            btp%Qv_face_ave(3,iquad,iface) = btp%Qv_face_ave(3,iquad,iface) + 0.5*(Qv_ql(3)+Qv_qr(3))
            btp%Qw_face_ave(1,iquad,iface) = btp%Qw_face_ave(1,iquad,iface) + 0.5*(Qw_ql(1)+Qw_qr(1))
            btp%Qw_face_ave(2,iquad,iface) = btp%Qw_face_ave(2,iquad,iface) + 0.5*(Qw_ql(2)+Qw_qr(2))
            btp%Qw_face_ave(3,iquad,iface) = btp%Qw_face_ave(3,iquad,iface) + 0.5*(Qw_ql(3)+Qw_qr(3))
            btp%uvb_face_ave(3,1,iquad,iface) = btp%uvb_face_ave(3,1,iquad,iface) + wl
            btp%uvb_face_ave(3,2,iquad,iface) = btp%uvb_face_ave(3,2,iquad,iface) + wr
         end if

         Qu_ql(1) = Qu_ql(1) + oe2_Hql
         Qu_qr(1) = Qu_qr(1) + oe2_Hqr
         fxl = nxl*Qu_ql(1) + nyl*Qu_ql(2) + nzl*Qu_ql(3)
         fxr = nxr*Qu_qr(1) + nyr*Qu_qr(2) + nzr*Qu_qr(3)
         flux(2,iquad) = 0.5*(fxl-fxr) - (0.25*clam)*(qbr(3)-qbl(3))

         Qv_ql(2) = Qv_ql(2) + oe2_Hql
         Qv_qr(2) = Qv_qr(2) + oe2_Hqr
         fxl = nxl*Qv_ql(1) + nyl*Qv_ql(2) + nzl*Qv_ql(3)
         fxr = nxr*Qv_qr(1) + nyr*Qv_qr(2) + nzr*Qv_qr(3)
         flux(3,iquad) = 0.5*(fxl-fxr) - (0.25*clam)*(qbr(4)-qbl(4))

         flux(4,iquad) = 0.0
         if (has_w) then
            Qw_ql(3) = Qw_ql(3) + oe2_Hql
            Qw_qr(3) = Qw_qr(3) + oe2_Hqr
            fxl = nxl*Qw_ql(1) + nyl*Qw_ql(2) + nzl*Qw_ql(3)
            fxr = nxr*Qw_qr(1) + nyr*Qw_qr(2) + nzr*Qw_qr(3)
            flux(4,iquad) = 0.5*(fxl-fxr) - (0.25*clam)*(qbr(5)-qbl(5))
         end if

      end do ! iquad pass 1

      ! Scatter flux into rhs
      ! Implicit barrier after pass 1 ensures flux is fully written before read.
      !$acc loop vector
      do iquad = 1, nq_f
         wq      = mf%jac_faceq(iquad,1,iface)
         flux_pb = flux(1,iquad)
         flux_u  = flux(2,iquad)
         flux_v  = flux(3,iquad)
         flux_w  = flux(4,iquad)
         !$acc loop seq
         do n = 1, ngl_f
            hi = b%psiq(n,iquad)
            il = mf%imapl(1,n,1,iface)
            jl = mf%imapl(2,n,1,iface)
            kl = mf%imapl(3,n,1,iface)
            I  = G%intma(il,jl,kl,el)
            !$acc atomic update
            rhs(1,I) = rhs(1,I) - wq*hi*flux_pb
            !$acc atomic update
            rhs(2,I) = rhs(2,I) - wq*hi*flux_u
            !$acc atomic update
            rhs(3,I) = rhs(3,I) - wq*hi*flux_v
            if (has_w) then
               !$acc atomic update
               rhs(4,I) = rhs(4,I) - wq*hi*flux_w
            end if
         end do
      end do ! iquad pass 2

   end do
   !$acc end parallel loop

end subroutine create_nbhs_face_df

!----------------------------------------------------------------------!
!>@brief Constructs the Inter-Processor Fluxes that lie on the boundary of each processor element
!>@details Requires calling send_receive_boundary, which gets the data lying off-processor and
!> populates q_send and q_recv
!>@author James F. Kelly
!>@date 16 November 2010
!>@date September 7, 2015 by F.X. Giraldo to reuse new ELEMENTAL_FLUX and RUSANOV
!>routines that use long vectors
!>@date January 2016 by M.A. Kopera - accounts for non-conforming faces
!>@date April 2024 modified by Yao Gahounzo
!>@date May 2026 GPU port by Yao Gahounzo
!>  One gang per MPI boundary face.  q_send/q_recv must be present on device
!>  before calling (upload via !$acc update device after mpi_waitall + unpack).
!>  rhs and btp%*_face_ave are already on device.  rhs is scattered with
!>  atomics; face_ave entries are (iquad,iface)-unique per gang — no atomics.
!----------------------------------------------------------------------!

subroutine create_nbhs_face_df_v0(G, inp, b, mf, par, btp, init, rhs, q_send, q_recv, nvarb)

   use mod_grid,      only: grid
   use mod_input,     only: input
   use mod_basis,     only: basis
   use mod_face,      only: face_CS
   use mod_parallel,  only: parallel_CS
   use mod_variables, only: btp_CS
   use mod_initial,   only: initial
   use mod_constants, only: gravity

   implicit none

   type(grid),        intent(in)    :: G
   type(input),       intent(in)    :: inp
   type(basis),       intent(in)    :: b
   type(face_CS),     intent(in)    :: mf
   type(parallel_CS), intent(in)    :: par
   type(btp_CS),      intent(inout) :: btp
   type(initial),     intent(in)    :: init

   integer, intent(in) :: nvarb
   real, intent(in)    :: q_send(nvarb, b%ngl, G%nboun)
   real, intent(in)    :: q_recv(nvarb, b%ngl, G%nboun)
   real, intent(inout) :: rhs(3, G%npoin)

   integer :: jj, iface, iquad, el, il, jl, kl, I, n, k
   real    :: wq, hi, nxl, nyl
   real    :: ul, ur, vl, vr, pbl, pbr, clam, one_eta, one_eta2
   real    :: pU_L, pU_R, pbpert_edge, half_clam, quarter_clam, half_clam_dqb2
   real    :: qbl(4), qbr(4)
   real    :: c_minus, c_plus
   real    :: H_bcl_ql, H_bcl_qr, H_bcl_q, oe2_Hql, oe2_Hqr
   real    :: Qu_ql1, Qu_ql2, Qu_qr1, Qu_qr2
   real    :: Qv_ql1, Qv_ql2, Qv_qr1, Qv_qr2
   real    :: flux_edge_x, flux_edge_y, flux_pb, flux_u, flux_v
   real    :: pprime_lk, pprime_rk, pprime_lk1, pprime_rk1
   real    :: pkl, pkr, ukl, ukr, vkl, vkr
   real    :: ope_ppl_k, ope_ppr_k, uv_cross_l, uv_cross_r
   real    :: opl, opr, wq_flux_pb, wq_flux_u, wq_flux_v
   real,    dimension(b%ngl) :: hi_c
   integer, dimension(b%ngl) :: I_l

   ! Capture derived-type scalars for GPU firstprivate.
   integer :: ngl_f, nq_f, nlayers_f, nboun_f
   ngl_f     = b%ngl
   nq_f      = b%nq
   nlayers_f = inp%nlayers
   nboun_f   = G%nboun

   if (nboun_f == 0) return

   !$acc data present(btp, rhs, G%face, G%intma, par%nbh_send_recv,         &
   !$acc              mf%imapl, mf%normal_vector_q, mf%jac_faceq,           &
   !$acc              b%psiq, init%alpha_mlswe,                              &
   !$acc              btp%H_face_ave, btp%ope_face_ave, btp%ope2_face_ave,   &
   !$acc              btp%btp_mass_flux_face_ave, btp%Qu_face_ave,           &
   !$acc              btp%Qv_face_ave, btp%one_plus_eta_edge_2_ave,          &
   !$acc              btp%uvb_face_ave, q_send, q_recv)

   !$acc parallel loop gang                                                   &
   !$acc   private(iface, el, il, jl, kl, I, n, k,                           &
   !$acc           I_l, hi_c,                                                 &
   !$acc           qbl, qbr, pbl, pbr,                                        &
   !$acc           nxl, nyl, wq, hi,                                          &
   !$acc           clam, one_eta, one_eta2, half_clam, quarter_clam,          &
   !$acc           pU_L, pU_R, half_clam_dqb2, pbpert_edge,                   &
   !$acc           c_minus, c_plus, ul, ur, vl, vr, opl, opr,                 &
   !$acc           H_bcl_ql, H_bcl_qr, H_bcl_q, oe2_Hql, oe2_Hqr,           &
   !$acc           Qu_ql1, Qu_ql2, Qu_qr1, Qu_qr2,                           &
   !$acc           Qv_ql1, Qv_ql2, Qv_qr1, Qv_qr2,                           &
   !$acc           flux_edge_x, flux_edge_y, flux_pb, flux_u, flux_v,         &
   !$acc           pprime_lk, pprime_rk, pprime_lk1, pprime_rk1,              &
   !$acc           pkl, pkr, ukl, ukr, vkl, vkr,                              &
   !$acc           ope_ppl_k, ope_ppr_k, uv_cross_l, uv_cross_r,             &
   !$acc           wq_flux_pb, wq_flux_u, wq_flux_v)
   do jj = 1, nboun_f

      iface = par%nbh_send_recv(jj)
      el    = G%face(7, iface)

      ! Precompute left element node indices once per face.
      !$acc loop seq
      do n = 1, ngl_f
         il     = mf%imapl(1, n, 1, iface)
         jl     = mf%imapl(2, n, 1, iface)
         kl     = mf%imapl(3, n, 1, iface)
         I_l(n) = G%intma(il, jl, kl, el)
      end do

      !$acc loop seq
      do iquad = 1, nq_f

         nxl = mf%normal_vector_q(1, iquad, 1, iface)
         nyl = mf%normal_vector_q(2, iquad, 1, iface)

         !$acc loop seq
         do n = 1, ngl_f
            hi_c(n) = b%psiq(n, iquad)
         end do

         ! Project LEFT state from packed send buffer (this processor's boundary data).
         qbl(1) = 0.0;  qbl(2) = 0.0;  qbl(3) = 0.0;  qbl(4) = 0.0
         pbl = 0.0
         !$acc loop seq
         do n = 1, ngl_f
            hi     = hi_c(n)
            qbl(1) = qbl(1) + hi * q_send(1, n, jj)
            qbl(2) = qbl(2) + hi * q_send(2, n, jj)
            qbl(3) = qbl(3) + hi * q_send(3, n, jj)
            qbl(4) = qbl(4) + hi * q_send(4, n, jj)
            pbl    = pbl    + hi * q_send(5, n, jj)
         end do

         ! Project RIGHT state from receive buffer (neighbour's boundary data).
         qbr(1) = 0.0;  qbr(2) = 0.0;  qbr(3) = 0.0;  qbr(4) = 0.0
         pbr = 0.0
         !$acc loop seq
         do n = 1, ngl_f
            hi     = hi_c(n)
            qbr(1) = qbr(1) + hi * q_recv(1, n, jj)
            qbr(2) = qbr(2) + hi * q_recv(2, n, jj)
            qbr(3) = qbr(3) + hi * q_recv(3, n, jj)
            qbr(4) = qbr(4) + hi * q_recv(4, n, jj)
            pbr    = pbr    + hi * q_recv(5, n, jj)
         end do

         ! Wave speeds.
         pU_L = nxl*qbl(3) + nyl*qbl(4)
         pU_R = -(nxl*qbr(3) + nyl*qbr(4))

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

         ! Initialise flux tensors from barotropic state.
         Qu_ql1 = ul * qbl(3);  Qu_ql2 = vl * qbl(3)
         Qu_qr1 = ur * qbr(3);  Qu_qr2 = vr * qbr(3)
         Qv_ql1 = ul * qbl(4);  Qv_ql2 = vl * qbl(4)
         Qv_qr1 = ur * qbr(4);  Qv_qr2 = vr * qbr(4)
         H_bcl_ql  = 0.0;  H_bcl_qr  = 0.0
         pprime_lk = 0.0;  pprime_rk = 0.0

         ! Layer projection + flux accumulation.
         !$acc loop seq
         do k = 1, nlayers_f
            pkl = 0.0;  pkr = 0.0
            ukl = 0.0;  ukr = 0.0
            vkl = 0.0;  vkr = 0.0
            !$acc loop seq
            do n = 1, ngl_f
               hi  = hi_c(n)
               pkl = pkl + hi * q_send(5 + (k-1)*3 + 1, n, jj)
               pkr = pkr + hi * q_recv(5 + (k-1)*3 + 1, n, jj)
               ukl = ukl + hi * q_send(5 + (k-1)*3 + 2, n, jj)
               ukr = ukr + hi * q_recv(5 + (k-1)*3 + 2, n, jj)
               vkl = vkl + hi * q_send(5 + (k-1)*3 + 3, n, jj)
               vkr = vkr + hi * q_recv(5 + (k-1)*3 + 3, n, jj)
            end do

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

            pprime_lk1 = pprime_lk + pkl
            pprime_rk1 = pprime_rk + pkr
            H_bcl_ql   = H_bcl_ql + 0.5*init%alpha_mlswe(k) * (pprime_lk1 + pprime_lk) * pkl
            H_bcl_qr   = H_bcl_qr + 0.5*init%alpha_mlswe(k) * (pprime_rk1 + pprime_rk) * pkr
            pprime_lk  = pprime_lk1
            pprime_rk  = pprime_rk1
         end do  ! k

         H_bcl_q = 0.5 * (H_bcl_ql + H_bcl_qr)

         half_clam_dqb2 = half_clam * (qbl(2) - qbr(2))
         flux_edge_x    = 0.5*(qbl(3) + qbr(3)) + nxl * half_clam_dqb2
         flux_edge_y    = 0.5*(qbl(4) + qbr(4)) + nyl * half_clam_dqb2
         flux_pb        = nxl*flux_edge_x + nyl*flux_edge_y

         H_bcl_q = one_eta2 * H_bcl_q

         ! Time-average accumulators — (iquad,iface) unique per gang, no atomics.
         btp%btp_mass_flux_face_ave(1,iquad,iface) = btp%btp_mass_flux_face_ave(1,iquad,iface) + flux_edge_x
         btp%btp_mass_flux_face_ave(2,iquad,iface) = btp%btp_mass_flux_face_ave(2,iquad,iface) + flux_edge_y
         btp%H_face_ave(iquad,iface)    = btp%H_face_ave(iquad,iface)    + H_bcl_q
         btp%Qu_face_ave(1,iquad,iface) = btp%Qu_face_ave(1,iquad,iface) + 0.5*(Qu_ql1 + Qu_qr1)
         btp%Qu_face_ave(2,iquad,iface) = btp%Qu_face_ave(2,iquad,iface) + 0.5*(Qu_ql2 + Qu_qr2)
         btp%Qv_face_ave(1,iquad,iface) = btp%Qv_face_ave(1,iquad,iface) + 0.5*(Qv_ql1 + Qv_qr1)
         btp%Qv_face_ave(2,iquad,iface) = btp%Qv_face_ave(2,iquad,iface) + 0.5*(Qv_ql2 + Qv_qr2)
         opl = 1.0 + qbl(2)/pbl
         opr = 1.0 + qbr(2)/pbr
         btp%ope_face_ave(1,iquad,iface)  = btp%ope_face_ave(1,iquad,iface)  + opl
         btp%ope_face_ave(2,iquad,iface)  = btp%ope_face_ave(2,iquad,iface)  + opr
         btp%ope2_face_ave(1,iquad,iface) = btp%ope2_face_ave(1,iquad,iface) + opl*opl
         btp%ope2_face_ave(2,iquad,iface) = btp%ope2_face_ave(2,iquad,iface) + opr*opr
         btp%one_plus_eta_edge_2_ave(iquad,iface) = btp%one_plus_eta_edge_2_ave(iquad,iface) + one_eta2
         btp%uvb_face_ave(1,1,iquad,iface) = btp%uvb_face_ave(1,1,iquad,iface) + ul
         btp%uvb_face_ave(1,2,iquad,iface) = btp%uvb_face_ave(1,2,iquad,iface) + ur
         btp%uvb_face_ave(2,1,iquad,iface) = btp%uvb_face_ave(2,1,iquad,iface) + vl
         btp%uvb_face_ave(2,2,iquad,iface) = btp%uvb_face_ave(2,2,iquad,iface) + vr

         ! Momentum fluxes.
         oe2_Hql = one_eta2 * H_bcl_ql
         oe2_Hqr = one_eta2 * H_bcl_qr

         Qu_ql1 = Qu_ql1 + oe2_Hql
         Qu_qr1 = Qu_qr1 + oe2_Hqr
         flux_u = 0.5*(nxl*(Qu_ql1 + Qu_qr1) + nyl*(Qu_ql2 + Qu_qr2)) &
            - quarter_clam * (qbr(3) - qbl(3))

         Qv_ql2 = Qv_ql2 + oe2_Hql
         Qv_qr2 = Qv_qr2 + oe2_Hqr
         flux_v = 0.5*(nxl*(Qv_ql1 + Qv_qr1) + nyl*(Qv_ql2 + Qv_qr2)) &
            - quarter_clam * (qbr(4) - qbl(4))

         wq         = mf%jac_faceq(iquad, 1, iface)
         wq_flux_pb = wq * flux_pb
         wq_flux_u  = wq * flux_u
         wq_flux_v  = wq * flux_v

         ! RHS scatter — left element only (right element is on the MPI neighbour).
         ! Atomics: different MPI boundary faces of the same element share nodes.
         !$acc loop seq
         do n = 1, ngl_f
            hi = hi_c(n)
            I  = I_l(n)
            !$acc atomic update
            rhs(1,I) = rhs(1,I) - hi * wq_flux_pb
            !$acc atomic update
            rhs(2,I) = rhs(2,I) - hi * wq_flux_u
            !$acc atomic update
            rhs(3,I) = rhs(3,I) - hi * wq_flux_v
         end do

      end do  ! iquad

   end do  ! jj (MPI boundary faces)
   !$acc end data

end subroutine create_nbhs_face_df_v0

subroutine create_nbhs_face_df_lap(G, b, mf, par, btp, rhs, ref, nvarb)

   use mod_grid,      only: grid
   use mod_basis,     only: basis
   use mod_face,      only: face_CS
   use mod_parallel,  only: parallel_CS
   use mod_variables, only: btp_CS
   use mod_ref,       only: mref

   implicit none

   type(grid),        intent(in)    :: G
   type(basis),       intent(in)    :: b
   type(face_CS),     intent(in)    :: mf
   type(parallel_CS), intent(in)    :: par
   type(btp_CS),      intent(inout) :: btp
   type(mref),        intent(in)    :: ref

   !global arrays
   integer, intent(in) :: nvarb
   real, intent(inout) :: rhs(2, G%npoin)

   !local variables

   integer :: kk, ivar, iquad, i, ip, iface, el, il, jl, kl
   integer :: ngl_f, nboun_f
   real :: wq, iflux, alpha, beta, nxl, nyl, hi, flux_qu, flux_qv

   real, dimension(2) :: qu_mean, qv_mean
   real, dimension(4,2) :: flux_uv_visc_face
   real, dimension(2) :: qul, qur, qvl, qvr
   real, dimension(5) :: grad_uvb_pb_l, grad_uvb_pb_r

   ngl_f   = b%ngl
   nboun_f = G%nboun
   beta    = 0.5
   alpha   = 1.0 - beta
   iflux   = 0.0

   !$acc parallel loop gang                                                  &
   !$acc    private(iface, el, iquad, i, ip, ivar,                           &
   !$acc            nxl, nyl, wq, hi, flux_qu, flux_qv,                      &
   !$acc            qu_mean, qv_mean, flux_uv_visc_face,                     &
   !$acc            qul, qur, qvl, qvr, grad_uvb_pb_l, grad_uvb_pb_r,       &
   !$acc            il, jl, kl)                                              &
   !$acc    present(G%face, G%intma, mf%imapl, mf%normal_vector,             &
   !$acc            mf%jac_face, b%psi,                                      &
   !$acc            par%nbh_send_recv,                                       &
   !$acc            ref%q_send_lap, ref%q_recv_lap,                          &
   !$acc            btp%graduvb_face_ave, rhs)                               &
   !$acc    firstprivate(ngl_f, nboun_f, alpha, beta, iflux)
   do kk = 1, nboun_f
      iface = par%nbh_send_recv(kk)
      el    = G%face(7, iface)

      !$acc loop seq
      do iquad = 1, ngl_f

         do ivar = 6, 10
            grad_uvb_pb_l(ivar-5) = ref%q_send_lap(ivar, iquad, kk)
            grad_uvb_pb_r(ivar-5) = ref%q_recv_lap(ivar, iquad, kk)
         end do

         do ivar = 1, 4
            flux_uv_visc_face(ivar,1) = grad_uvb_pb_l(5)*ref%q_send_lap(ivar,iquad,kk) + grad_uvb_pb_l(ivar)
            flux_uv_visc_face(ivar,2) = grad_uvb_pb_r(5)*ref%q_recv_lap(ivar,iquad,kk) + grad_uvb_pb_r(ivar)
            btp%graduvb_face_ave(ivar,1,iquad,iface) = btp%graduvb_face_ave(ivar,1,iquad,iface) + ref%q_send_lap(ivar,iquad,kk)
            btp%graduvb_face_ave(ivar,2,iquad,iface) = btp%graduvb_face_ave(ivar,2,iquad,iface) + ref%q_recv_lap(ivar,iquad,kk)
         end do

         nxl = mf%normal_vector(1, iquad, 1, iface)
         nyl = mf%normal_vector(2, iquad, 1, iface)

         qul(1) = flux_uv_visc_face(1,1);  qul(2) = flux_uv_visc_face(2,1)
         qvl(1) = flux_uv_visc_face(3,1);  qvl(2) = flux_uv_visc_face(4,1)
         qur(1) = flux_uv_visc_face(1,2);  qur(2) = flux_uv_visc_face(2,2)
         qvr(1) = flux_uv_visc_face(3,2);  qvr(2) = flux_uv_visc_face(4,2)

         ! Central flux (beta=0.5)
         qu_mean(1) = alpha*qul(1) + beta*qur(1);  qu_mean(2) = alpha*qul(2) + beta*qur(2)
         qv_mean(1) = alpha*qvl(1) + beta*qvr(1);  qv_mean(2) = alpha*qvl(2) + beta*qvr(2)

         wq = mf%jac_face(iquad, 1, iface)

         flux_qu = (qu_mean(1) - iflux*qul(1))*nxl + (qu_mean(2) - iflux*qul(2))*nyl
         flux_qv = (qv_mean(1) - iflux*qvl(1))*nxl + (qv_mean(2) - iflux*qvl(2))*nyl

         !$acc loop seq
         do i = 1, ngl_f
            hi = b%psi(i, iquad)
            il = mf%imapl(1, i, 1, iface)
            jl = mf%imapl(2, i, 1, iface)
            kl = mf%imapl(3, i, 1, iface)
            ip = G%intma(il, jl, kl, el)
            !$acc atomic update
            rhs(1, ip) = rhs(1, ip) + wq*hi*flux_qu
            !$acc atomic update
            rhs(2, ip) = rhs(2, ip) + wq*hi*flux_qv
         end do

      end do  ! iquad
   end do  ! kk
   !$acc end parallel loop

end subroutine create_nbhs_face_df_lap

subroutine create_nbhs_face_bcl(G, inp, b, mf, par, btp, init, ref, rhs, q_send_bcl, q_recv_bcl)

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

   real, intent(inout) :: rhs(3, G%npoin, inp%nlayers)
   real, intent(in)    :: q_send_bcl(inp%nvar_bcl*inp%nlayers, b%ngl, par%num_send_recv_total)
   real, intent(in)    :: q_recv_bcl(inp%nvar_bcl*inp%nlayers, b%ngl, par%num_send_recv_total)

   integer :: kk, iquad, k, ktemp, n, ivar, iface, el, il, jl, kl, I
   integer :: ngl_f, nq_f, nlayers_f, nvarb_f, nboun_valid_f

   real, dimension(inp%nlayers)   :: alpha_over_g, g_over_alpha
   real, dimension(2,inp%nlayers+1) :: p_face, z_face
   real, dimension(inp%nlayers+1) :: p_edge_plus, p_edge_minus, z_edge_plus, z_edge_minus
   real, dimension(inp%nvar_bcl,inp%nlayers) :: ql, qr
   real, dimension(2,inp%nlayers) :: dp_flux, udp_flux, vdp_flux
   real, dimension(inp%nlayers)   :: H_face_q, udpl, udpr, vdpl, vdpr, dp_lr_l, dp_lr_r
   real :: dp_deficit(2), uu_dp_flux_deficit(2), vv_dp_flux_deficit(2)
   real :: ope_l, ope_r, one_plus_eta_edge
   real :: ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv
   real :: wq, hi, weight, acceleration
   real :: z_intersect_top, z_intersect_bot, dz_intersect, H_r_plus, H_r_minus
   real :: p_intersect_bot, p_intersect_top, H_corr1, p_inc1
   real :: flux, flux_x, flux_y
   real, parameter :: eps1 = 1.0e-20

   ngl_f         = b%ngl
   nq_f          = b%nq
   nlayers_f     = inp%nlayers
   nvarb_f       = inp%nvar_bcl
   nboun_valid_f = ref%nboun_valid

   !$acc data present(ref%face_pack_list, G%face, G%intma,                        &
   !$acc              mf%imapl, mf%normal_vector_q, mf%jac_faceq,                 &
   !$acc              b%psiq,                                                       &
   !$acc              btp%ope_face_ave, btp%uvb_face_ave, btp%ope2_face_ave,       &
   !$acc              btp%one_plus_eta_edge_2_ave, btp%H_face_ave,                 &
   !$acc              btp%btp_mass_flux_face_ave, btp%Qu_face_ave, btp%Qv_face_ave,&
   !$acc              init%alpha_mlswe, init%zbot_face,                            &
   !$acc              q_send_bcl, q_recv_bcl, rhs)

   !$acc parallel loop gang                                                         &
   !$acc   private(ql, qr, alpha_over_g, g_over_alpha,                             &
   !$acc           p_face, z_face, p_edge_plus, p_edge_minus, z_edge_plus, z_edge_minus, &
   !$acc           dp_flux, udp_flux, vdp_flux, H_face_q,                          &
   !$acc           dp_lr_l, dp_lr_r, udpl, udpr, vdpl, vdpr,                       &
   !$acc           dp_deficit, uu_dp_flux_deficit, vv_dp_flux_deficit,             &
   !$acc           ope_l, ope_r, one_plus_eta_edge,                                &
   !$acc           ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv,                    &
   !$acc           wq, hi, weight, acceleration,                                    &
   !$acc           z_intersect_top, z_intersect_bot, dz_intersect,                 &
   !$acc           H_r_plus, H_r_minus, p_intersect_bot, p_intersect_top,          &
   !$acc           H_corr1, p_inc1, flux, flux_x, flux_y,                          &
   !$acc           kk, iquad, k, ktemp, n, ivar, iface, el, il, jl, kl, I)        &
   !$acc   firstprivate(ngl_f, nq_f, nlayers_f, nvarb_f, nboun_valid_f)
   do kk = 1, nboun_valid_f

      iface = ref%face_pack_list(kk)
      el    = G%face(7, iface)

      !$acc loop seq
      do k = 1, nlayers_f
         alpha_over_g(k) = init%alpha_mlswe(k) / gravity
         g_over_alpha(k) = gravity / init%alpha_mlswe(k)
      end do

      !$acc loop seq
      do iquad = 1, nq_f

         nxl = mf%normal_vector_q(1,iquad,1,iface)
         nyl = mf%normal_vector_q(2,iquad,1,iface)

         ql = 0.0; qr = 0.0
         !$acc loop seq
         do k = 1, nlayers_f
            !$acc loop seq
            do n = 1, ngl_f
               hi = b%psiq(n,iquad)
               !$acc loop seq
               do ivar = 1, nvarb_f
                  ql(ivar,k) = ql(ivar,k) + hi*q_send_bcl(nvarb_f*(k-1)+ivar, n, kk)
                  qr(ivar,k) = qr(ivar,k) + hi*q_recv_bcl(nvarb_f*(k-1)+ivar, n, kk)
               end do
            end do

            dpl = btp%ope_face_ave(1,iquad,iface) * ql(1,k)
            dpr = btp%ope_face_ave(2,iquad,iface) * qr(1,k)
            dp_lr_l(k) = dpl
            dp_lr_r(k) = dpr

            ul = ql(2,k) + btp%uvb_face_ave(1,1,iquad,iface)
            ur = qr(2,k) + btp%uvb_face_ave(1,2,iquad,iface)
            vl = ql(3,k) + btp%uvb_face_ave(2,1,iquad,iface)
            vr = qr(3,k) + btp%uvb_face_ave(2,2,iquad,iface)

            uu = 0.5*(ul+ur)
            vv = 0.5*(vl+vr)
            udpl(k) = ul*dpl;  udpr(k) = ur*dpr
            vdpl(k) = vl*dpl;  vdpr(k) = vr*dpr

            if (uu*nxl > 0.0) then
               dp_flux(1,k)  = uu * dpl;  udp_flux(1,k) = uu*(ul*dpl);  vdp_flux(1,k) = uu*(vl*dpl)
            else
               dp_flux(1,k)  = uu * dpr;  udp_flux(1,k) = uu*(ur*dpr);  vdp_flux(1,k) = uu*(vr*dpr)
            end if
            if (vv*nyl > 0.0) then
               dp_flux(2,k)  = vv * dpl;  udp_flux(2,k) = vv*(ul*dpl);  vdp_flux(2,k) = vv*(vl*dpl)
            else
               dp_flux(2,k)  = vv * dpr;  udp_flux(2,k) = vv*(ur*dpr);  vdp_flux(2,k) = vv*(vr*dpr)
            end if
         end do  ! k

         uu_dp_flux_deficit(1) = btp%Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,:))
         uu_dp_flux_deficit(2) = btp%Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,:))
         vv_dp_flux_deficit(1) = btp%Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,:))
         vv_dp_flux_deficit(2) = btp%Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,:))
         dp_deficit(1) = btp%btp_mass_flux_face_ave(1,iquad,iface) - sum(dp_flux(1,:))
         dp_deficit(2) = btp%btp_mass_flux_face_ave(2,iquad,iface) - sum(dp_flux(2,:))

         !$acc loop seq
         do k = 1, nlayers_f
            weight = dp_lr_l(k) / (sum(abs(dp_lr_l(:))+eps1))
            if (dp_deficit(1)*nxl < 0.0) weight = dp_lr_r(k) / (sum(abs(dp_lr_r(:))+eps1))
            dp_flux(1,k) = dp_flux(1,k) + weight * dp_deficit(1)

            weight = dp_lr_l(k) / (sum(abs(dp_lr_l(:))+eps1))
            if (dp_deficit(2)*nyl < 0.0) weight = dp_lr_r(k) / (sum(abs(dp_lr_r(:))+eps1))
            dp_flux(2,k) = dp_flux(2,k) + weight * dp_deficit(2)

            weight = abs(udpl(k)) / (sum(abs(udpl(:))+eps1))
            if (uu_dp_flux_deficit(1)*nxl < 0.0) weight = abs(udpr(k)) / (sum(abs(udpr(:))+eps1))
            udp_flux(1,k) = udp_flux(1,k) + weight * uu_dp_flux_deficit(1)

            weight = abs(udpl(k)) / (sum(abs(udpl(:))+eps1))
            if (uu_dp_flux_deficit(2)*nyl < 0.0) weight = abs(udpr(k)) / (sum(abs(udpr(:))+eps1))
            udp_flux(2,k) = udp_flux(2,k) + weight * uu_dp_flux_deficit(2)

            weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
            if (vv_dp_flux_deficit(1)*nxl < 0.0) weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
            vdp_flux(1,k) = vdp_flux(1,k) + weight * vv_dp_flux_deficit(1)

            weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
            if (vv_dp_flux_deficit(2)*nyl < 0.0) weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
            vdp_flux(2,k) = vdp_flux(2,k) + weight * vv_dp_flux_deficit(2)
         end do

         z_face = 0.0;  p_face = 0.0
         z_edge_plus = 0.0;  z_edge_minus = 0.0
         p_edge_plus = 0.0;  p_edge_minus = 0.0

         ope_l = sqrt(btp%ope2_face_ave(1,iquad,iface))
         ope_r = sqrt(btp%ope2_face_ave(2,iquad,iface))
         p_face(1,1) = 0.0;  p_face(2,1) = 0.0
         !$acc loop seq
         do k = 1, nlayers_f
            p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,k)
            p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,k)
         end do

         one_plus_eta_edge = sqrt(btp%one_plus_eta_edge_2_ave(iquad,iface))
         z_face(1,nlayers_f+1)    = init%zbot_face(1,iquad,iface)
         z_face(2,nlayers_f+1)    = init%zbot_face(2,iquad,iface)
         z_edge_plus(nlayers_f+1) = init%zbot_face(1,iquad,iface)
         z_edge_minus(nlayers_f+1)= init%zbot_face(2,iquad,iface)
         !$acc loop seq
         do k = nlayers_f, 1, -1
            z_face(1,k)    = z_face(1,k+1)    + alpha_over_g(k)*(ope_l * ql(1,k))
            z_face(2,k)    = z_face(2,k+1)    + alpha_over_g(k)*(ope_r * qr(1,k))
            z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k)*(one_plus_eta_edge * ql(1,k))
            z_edge_minus(k)= z_edge_minus(k+1)+ alpha_over_g(k)*(one_plus_eta_edge * qr(1,k))
         end do

         p_edge_plus(2)  = one_plus_eta_edge * ql(1,1)
         p_edge_minus(2) = one_plus_eta_edge * qr(1,1)
         !$acc loop seq
         do k = 2, nlayers_f
            p_edge_plus(k+1)  = p_edge_plus(k)  + one_plus_eta_edge * ql(1,k)
            p_edge_minus(k+1) = p_edge_minus(k) + one_plus_eta_edge * qr(1,k)
         end do

         !$acc loop seq
         do k = 1, nlayers_f
            H_r_plus  = 0.5*init%alpha_mlswe(k)*(p_edge_plus(k+1)**2 - p_edge_plus(k)**2)
            H_r_minus = 0.0
            !$acc loop seq
            do ktemp = 1, nlayers_f
               z_intersect_top = min(z_edge_minus(ktemp),   z_edge_plus(k))
               z_intersect_bot = max(z_edge_minus(ktemp+1), z_edge_plus(k+1))
               dz_intersect    = z_intersect_top - z_intersect_bot
               if (dz_intersect > 0.0) then
                  p_intersect_bot = p_edge_minus(ktemp+1) &
                     - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_minus(ktemp+1))
                  p_intersect_top = p_edge_minus(ktemp+1) &
                     - g_over_alpha(ktemp)*(z_intersect_top - z_edge_minus(ktemp+1))
                  H_r_minus = H_r_minus + &
                     0.5*init%alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
               end if
            end do
            H_face_q(k) = 0.5*(H_r_plus + H_r_minus)
         end do

         !$acc loop seq
         do k = 1, nlayers_f-1
            p_inc1  = g_over_alpha(k)*(z_face(1,k+1) - z_edge_plus(k+1))
            H_corr1 = 0.5*init%alpha_mlswe(k)*((p_face(1,k+1) + p_inc1)**2 - p_face(1,k+1)**2)
            H_face_q(k)   = H_face_q(k)   - H_corr1
            H_face_q(k+1) = H_face_q(k+1) + H_corr1
         end do

         weight = 1.0
         acceleration = sum(H_face_q(:))
         if (acceleration > 0.0) weight = btp%H_face_ave(iquad,iface) / acceleration
         H_face_q(:) = H_face_q(:) * weight

         wq = mf%jac_faceq(iquad,1,iface)
         !$acc loop seq
         do k = 1, nlayers_f
            flux   = nxl*dp_flux(1,k)                      + nyl*dp_flux(2,k)
            flux_x = nxl*(udp_flux(1,k) + H_face_q(k)) + nyl*udp_flux(2,k)
            flux_y = nxl*vdp_flux(1,k)  + nyl*(vdp_flux(2,k) + H_face_q(k))
            !$acc loop seq
            do n = 1, ngl_f
               hi = b%psiq(n,iquad)
               il = mf%imapl(1,n,1,iface)
               jl = mf%imapl(2,n,1,iface)
               kl = mf%imapl(3,n,1,iface)
               I  = G%intma(il,jl,kl,el)
               !$acc atomic update
               rhs(1,I,k) = rhs(1,I,k) - wq*hi*flux
               !$acc atomic update
               rhs(2,I,k) = rhs(2,I,k) - wq*hi*flux_x
               !$acc atomic update
               rhs(3,I,k) = rhs(3,I,k) - wq*hi*flux_y
            end do
         end do

      end do  ! iquad

   end do  ! kk
   !$acc end parallel loop
   !$acc end data

end subroutine create_nbhs_face_bcl

! sphere_hex counterpart of create_nbhs_face_bcl: adds the w-momentum row
! (mass,u,v,w -> 4-wide ql/qr, 3-component flux tensors, nzl normal, rhs(4,...)).
! Forked rather than has_w-gated per the reference (hnumo-sphere) design, since
! the baroclinic RHS pipeline there is itself split cartesian/sphere.
subroutine create_nbhs_face_bcl_sphere(G, inp, b, mf, par, btp, init, ref, rhs, q_send_bcl, q_recv_bcl)

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

   real, intent(inout) :: rhs(inp%nvar_bcl, G%npoin, inp%nlayers)
   real, intent(in)    :: q_send_bcl(inp%nvar_bcl*inp%nlayers, b%ngl, par%num_send_recv_total)
   real, intent(in)    :: q_recv_bcl(inp%nvar_bcl*inp%nlayers, b%ngl, par%num_send_recv_total)

   integer :: kk, iquad, k, ktemp, n, ivar, iface, el, il, jl, kl, I
   integer :: ngl_f, nq_f, nlayers_f, nvarb_f, nboun_valid_f

   real, dimension(inp%nlayers)     :: alpha_over_g, g_over_alpha
   real, dimension(2,inp%nlayers+1) :: p_face, z_face
   real, dimension(inp%nlayers+1)   :: p_edge_plus, p_edge_minus, z_edge_plus, z_edge_minus
   real, dimension(inp%nvar_bcl,inp%nlayers) :: ql, qr
   real, dimension(3,inp%nlayers)   :: dp_flux, udp_flux, vdp_flux, wdp_flux
   real, dimension(inp%nlayers)     :: H_face_q, udpl, udpr, vdpl, vdpr, wdpl, wdpr, dp_lr_l, dp_lr_r
   real :: dp_deficit(3), uu_dp_flux_deficit(3), vv_dp_flux_deficit(3), ww_dp_flux_deficit(3)
   real :: ope_l, ope_r, one_plus_eta_edge
   real :: ul, ur, vl, vr, wl, wr, dpl, dpr, nxl, nyl, nzl, uu, vv, ww, un
   real :: wq, hi, weight, acceleration
   real :: z_intersect_top, z_intersect_bot, dz_intersect, H_r_plus, H_r_minus
   real :: p_intersect_bot, p_intersect_top, H_corr1, p_inc1
   real :: flux, flux_x, flux_y, flux_z
   real, parameter :: eps1 = 1.0e-20

   ngl_f         = b%ngl
   nq_f          = b%nq
   nlayers_f     = inp%nlayers
   nvarb_f       = inp%nvar_bcl
   nboun_valid_f = ref%nboun_valid

   !$acc data present(ref%face_pack_list, G%face, G%intma,                        &
   !$acc              mf%imapl, mf%normal_vector_q, mf%jac_faceq,                 &
   !$acc              b%psiq,                                                       &
   !$acc              btp%ope_face_ave, btp%uvb_face_ave, btp%ope2_face_ave,       &
   !$acc              btp%one_plus_eta_edge_2_ave, btp%H_face_ave,                 &
   !$acc              btp%btp_mass_flux_face_ave, btp%Qu_face_ave, btp%Qv_face_ave,&
   !$acc              btp%Qw_face_ave,                                             &
   !$acc              init%alpha_mlswe, init%zbot_face,                            &
   !$acc              q_send_bcl, q_recv_bcl, rhs)

   !$acc parallel loop gang                                                         &
   !$acc   private(ql, qr, alpha_over_g, g_over_alpha,                             &
   !$acc           p_face, z_face, p_edge_plus, p_edge_minus, z_edge_plus, z_edge_minus, &
   !$acc           dp_flux, udp_flux, vdp_flux, wdp_flux, H_face_q,                &
   !$acc           dp_lr_l, dp_lr_r, udpl, udpr, vdpl, vdpr, wdpl, wdpr,           &
   !$acc           dp_deficit, uu_dp_flux_deficit, vv_dp_flux_deficit, ww_dp_flux_deficit, &
   !$acc           ope_l, ope_r, one_plus_eta_edge,                                &
   !$acc           ul, ur, vl, vr, wl, wr, dpl, dpr, nxl, nyl, nzl, uu, vv, ww, un, &
   !$acc           wq, hi, weight, acceleration,                                    &
   !$acc           z_intersect_top, z_intersect_bot, dz_intersect,                 &
   !$acc           H_r_plus, H_r_minus, p_intersect_bot, p_intersect_top,          &
   !$acc           H_corr1, p_inc1, flux, flux_x, flux_y, flux_z,                  &
   !$acc           kk, iquad, k, ktemp, n, ivar, iface, el, il, jl, kl, I)        &
   !$acc   firstprivate(ngl_f, nq_f, nlayers_f, nvarb_f, nboun_valid_f)
   do kk = 1, nboun_valid_f

      iface = ref%face_pack_list(kk)
      el    = G%face(7, iface)

      !$acc loop seq
      do k = 1, nlayers_f
         alpha_over_g(k) = init%alpha_mlswe(k) / gravity
         g_over_alpha(k) = gravity / init%alpha_mlswe(k)
      end do

      !$acc loop seq
      do iquad = 1, nq_f

         nxl = mf%normal_vector_q(1,iquad,1,iface)
         nyl = mf%normal_vector_q(2,iquad,1,iface)
         nzl = mf%normal_vector_q(3,iquad,1,iface)

         ql = 0.0; qr = 0.0
         !$acc loop seq
         do k = 1, nlayers_f
            !$acc loop seq
            do n = 1, ngl_f
               hi = b%psiq(n,iquad)
               !$acc loop seq
               do ivar = 1, nvarb_f
                  ql(ivar,k) = ql(ivar,k) + hi*q_send_bcl(nvarb_f*(k-1)+ivar, n, kk)
                  qr(ivar,k) = qr(ivar,k) + hi*q_recv_bcl(nvarb_f*(k-1)+ivar, n, kk)
               end do
            end do

            dpl = btp%ope_face_ave(1,iquad,iface) * ql(1,k)
            dpr = btp%ope_face_ave(2,iquad,iface) * qr(1,k)
            dp_lr_l(k) = dpl
            dp_lr_r(k) = dpr

            ul = ql(2,k) + btp%uvb_face_ave(1,1,iquad,iface)
            ur = qr(2,k) + btp%uvb_face_ave(1,2,iquad,iface)
            vl = ql(3,k) + btp%uvb_face_ave(2,1,iquad,iface)
            vr = qr(3,k) + btp%uvb_face_ave(2,2,iquad,iface)
            wl = ql(4,k) + btp%uvb_face_ave(3,1,iquad,iface)
            wr = qr(4,k) + btp%uvb_face_ave(3,2,iquad,iface)

            uu = 0.5*(ul+ur)
            vv = 0.5*(vl+vr)
            ww = 0.5*(wl+wr)
            udpl(k) = ul*dpl;  udpr(k) = ur*dpr
            vdpl(k) = vl*dpl;  vdpr(k) = vr*dpr
            wdpl(k) = wl*dpl;  wdpr(k) = wr*dpr

            un = uu*nxl + vv*nyl + ww*nzl
            if (un > 0.0) then
               dp_flux(1,k)  = uu * dpl;  udp_flux(1,k) = uu*(ul*dpl);  vdp_flux(1,k) = uu*(vl*dpl);  wdp_flux(1,k) = uu*(wl*dpl)
               dp_flux(2,k)  = vv * dpl;  udp_flux(2,k) = vv*(ul*dpl);  vdp_flux(2,k) = vv*(vl*dpl);  wdp_flux(2,k) = vv*(wl*dpl)
               dp_flux(3,k)  = ww * dpl;  udp_flux(3,k) = ww*(ul*dpl);  vdp_flux(3,k) = ww*(vl*dpl);  wdp_flux(3,k) = ww*(wl*dpl)
            else
               dp_flux(1,k)  = uu * dpr;  udp_flux(1,k) = uu*(ur*dpr);  vdp_flux(1,k) = uu*(vr*dpr);  wdp_flux(1,k) = uu*(wr*dpr)
               dp_flux(2,k)  = vv * dpr;  udp_flux(2,k) = vv*(ur*dpr);  vdp_flux(2,k) = vv*(vr*dpr);  wdp_flux(2,k) = vv*(wr*dpr)
               dp_flux(3,k)  = ww * dpr;  udp_flux(3,k) = ww*(ur*dpr);  vdp_flux(3,k) = ww*(vr*dpr);  wdp_flux(3,k) = ww*(wr*dpr)
            end if
         end do  ! k

         uu_dp_flux_deficit(1) = btp%Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,:))
         uu_dp_flux_deficit(2) = btp%Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,:))
         uu_dp_flux_deficit(3) = btp%Qu_face_ave(3,iquad,iface) - sum(udp_flux(3,:))
         vv_dp_flux_deficit(1) = btp%Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,:))
         vv_dp_flux_deficit(2) = btp%Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,:))
         vv_dp_flux_deficit(3) = btp%Qv_face_ave(3,iquad,iface) - sum(vdp_flux(3,:))
         ww_dp_flux_deficit(1) = btp%Qw_face_ave(1,iquad,iface) - sum(wdp_flux(1,:))
         ww_dp_flux_deficit(2) = btp%Qw_face_ave(2,iquad,iface) - sum(wdp_flux(2,:))
         ww_dp_flux_deficit(3) = btp%Qw_face_ave(3,iquad,iface) - sum(wdp_flux(3,:))
         dp_deficit(1) = btp%btp_mass_flux_face_ave(1,iquad,iface) - sum(dp_flux(1,:))
         dp_deficit(2) = btp%btp_mass_flux_face_ave(2,iquad,iface) - sum(dp_flux(2,:))
         dp_deficit(3) = btp%btp_mass_flux_face_ave(3,iquad,iface) - sum(dp_flux(3,:))

         !$acc loop seq
         do k = 1, nlayers_f
            weight = dp_lr_l(k) / (sum(abs(dp_lr_l(:))+eps1))
            if ((dp_deficit(1)*nxl + dp_deficit(2)*nyl + dp_deficit(3)*nzl) < 0.0) &
               weight = dp_lr_r(k) / (sum(abs(dp_lr_r(:))+eps1))
            dp_flux(1,k) = dp_flux(1,k) + weight * dp_deficit(1)
            dp_flux(2,k) = dp_flux(2,k) + weight * dp_deficit(2)
            dp_flux(3,k) = dp_flux(3,k) + weight * dp_deficit(3)

            weight = abs(udpl(k)) / (sum(abs(udpl(:))+eps1))
            if ((uu_dp_flux_deficit(1)*nxl + uu_dp_flux_deficit(2)*nyl + uu_dp_flux_deficit(3)*nzl) < 0.0) &
               weight = abs(udpr(k)) / (sum(abs(udpr(:))+eps1))
            udp_flux(1,k) = udp_flux(1,k) + weight * uu_dp_flux_deficit(1)
            udp_flux(2,k) = udp_flux(2,k) + weight * uu_dp_flux_deficit(2)
            udp_flux(3,k) = udp_flux(3,k) + weight * uu_dp_flux_deficit(3)

            weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
            if ((vv_dp_flux_deficit(1)*nxl + vv_dp_flux_deficit(2)*nyl + vv_dp_flux_deficit(3)*nzl) < 0.0) &
               weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
            vdp_flux(1,k) = vdp_flux(1,k) + weight * vv_dp_flux_deficit(1)
            vdp_flux(2,k) = vdp_flux(2,k) + weight * vv_dp_flux_deficit(2)
            vdp_flux(3,k) = vdp_flux(3,k) + weight * vv_dp_flux_deficit(3)

            weight = abs(wdpl(k)) / (sum(abs(wdpl(:))+eps1))
            if ((ww_dp_flux_deficit(1)*nxl + ww_dp_flux_deficit(2)*nyl + ww_dp_flux_deficit(3)*nzl) < 0.0) &
               weight = abs(wdpr(k)) / (sum(abs(wdpr(:))+eps1))
            wdp_flux(1,k) = wdp_flux(1,k) + weight * ww_dp_flux_deficit(1)
            wdp_flux(2,k) = wdp_flux(2,k) + weight * ww_dp_flux_deficit(2)
            wdp_flux(3,k) = wdp_flux(3,k) + weight * ww_dp_flux_deficit(3)
         end do

         z_face = 0.0;  p_face = 0.0
         z_edge_plus = 0.0;  z_edge_minus = 0.0
         p_edge_plus = 0.0;  p_edge_minus = 0.0

         ope_l = sqrt(btp%ope2_face_ave(1,iquad,iface))
         ope_r = sqrt(btp%ope2_face_ave(2,iquad,iface))
         p_face(1,1) = 0.0;  p_face(2,1) = 0.0
         !$acc loop seq
         do k = 1, nlayers_f
            p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,k)
            p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,k)
         end do

         one_plus_eta_edge = sqrt(btp%one_plus_eta_edge_2_ave(iquad,iface))
         z_face(1,nlayers_f+1)    = init%zbot_face(1,iquad,iface)
         z_face(2,nlayers_f+1)    = init%zbot_face(2,iquad,iface)
         z_edge_plus(nlayers_f+1) = init%zbot_face(1,iquad,iface)
         z_edge_minus(nlayers_f+1)= init%zbot_face(2,iquad,iface)
         !$acc loop seq
         do k = nlayers_f, 1, -1
            z_face(1,k)    = z_face(1,k+1)    + alpha_over_g(k)*(ope_l * ql(1,k))
            z_face(2,k)    = z_face(2,k+1)    + alpha_over_g(k)*(ope_r * qr(1,k))
            z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k)*(one_plus_eta_edge * ql(1,k))
            z_edge_minus(k)= z_edge_minus(k+1)+ alpha_over_g(k)*(one_plus_eta_edge * qr(1,k))
         end do

         p_edge_plus(2)  = one_plus_eta_edge * ql(1,1)
         p_edge_minus(2) = one_plus_eta_edge * qr(1,1)
         !$acc loop seq
         do k = 2, nlayers_f
            p_edge_plus(k+1)  = p_edge_plus(k)  + one_plus_eta_edge * ql(1,k)
            p_edge_minus(k+1) = p_edge_minus(k) + one_plus_eta_edge * qr(1,k)
         end do

         !$acc loop seq
         do k = 1, nlayers_f
            H_r_plus  = 0.5*init%alpha_mlswe(k)*(p_edge_plus(k+1)**2 - p_edge_plus(k)**2)
            H_r_minus = 0.0
            !$acc loop seq
            do ktemp = 1, nlayers_f
               z_intersect_top = min(z_edge_minus(ktemp),   z_edge_plus(k))
               z_intersect_bot = max(z_edge_minus(ktemp+1), z_edge_plus(k+1))
               dz_intersect    = z_intersect_top - z_intersect_bot
               if (dz_intersect > 0.0) then
                  p_intersect_bot = p_edge_minus(ktemp+1) &
                     - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_minus(ktemp+1))
                  p_intersect_top = p_edge_minus(ktemp+1) &
                     - g_over_alpha(ktemp)*(z_intersect_top - z_edge_minus(ktemp+1))
                  H_r_minus = H_r_minus + &
                     0.5*init%alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
               end if
            end do
            H_face_q(k) = 0.5*(H_r_plus + H_r_minus)
         end do

         !$acc loop seq
         do k = 1, nlayers_f-1
            p_inc1  = g_over_alpha(k)*(z_face(1,k+1) - z_edge_plus(k+1))
            H_corr1 = 0.5*init%alpha_mlswe(k)*((p_face(1,k+1) + p_inc1)**2 - p_face(1,k+1)**2)
            H_face_q(k)   = H_face_q(k)   - H_corr1
            H_face_q(k+1) = H_face_q(k+1) + H_corr1
         end do

         weight = 1.0
         acceleration = sum(H_face_q(:))
         if (acceleration > 0.0) weight = btp%H_face_ave(iquad,iface) / acceleration
         H_face_q(:) = H_face_q(:) * weight

         wq = mf%jac_faceq(iquad,1,iface)
         !$acc loop seq
         do k = 1, nlayers_f
            flux   = nxl*dp_flux(1,k)  + nyl*dp_flux(2,k)  + nzl*dp_flux(3,k)
            flux_x = nxl*(udp_flux(1,k) + H_face_q(k)) + nyl*udp_flux(2,k) + nzl*udp_flux(3,k)
            flux_y = nxl*vdp_flux(1,k) + nyl*(vdp_flux(2,k) + H_face_q(k)) + nzl*vdp_flux(3,k)
            flux_z = nxl*wdp_flux(1,k) + nyl*wdp_flux(2,k) + nzl*(wdp_flux(3,k) + H_face_q(k))
            !$acc loop seq
            do n = 1, ngl_f
               hi = b%psiq(n,iquad)
               il = mf%imapl(1,n,1,iface)
               jl = mf%imapl(2,n,1,iface)
               kl = mf%imapl(3,n,1,iface)
               I  = G%intma(il,jl,kl,el)
               !$acc atomic update
               rhs(1,I,k) = rhs(1,I,k) - wq*hi*flux
               !$acc atomic update
               rhs(2,I,k) = rhs(2,I,k) - wq*hi*flux_x
               !$acc atomic update
               rhs(3,I,k) = rhs(3,I,k) - wq*hi*flux_y
               !$acc atomic update
               rhs(4,I,k) = rhs(4,I,k) - wq*hi*flux_z
            end do
         end do

      end do  ! iquad

   end do  ! kk
   !$acc end parallel loop
   !$acc end data

end subroutine create_nbhs_face_bcl_sphere

subroutine create_nbhs_face_bcl_continuity(G, inp, b, mf, par, btp, rhs, q_send, q_recv, multirate)

   use mod_grid,      only: grid
   use mod_input,     only: input
   use mod_basis,     only: basis
   use mod_face,      only: face_CS
   use mod_parallel,  only: parallel_CS
   use mod_variables, only: btp_CS

   implicit none

   type(grid),        intent(in)    :: G
   type(input),       intent(in)    :: inp
   type(basis),       intent(in)    :: b
   type(face_CS),     intent(in)    :: mf
   type(parallel_CS), intent(in)    :: par
   type(btp_CS),      intent(inout) :: btp
   integer,           intent(in)    :: multirate

   !global arrays
   real, intent(inout) :: rhs(G%npoin, inp%nlayers)
   real, intent(in)    :: q_send(inp%nvar_bcl*inp%nlayers, b%ngl, par%num_send_recv_total)
   real, intent(in)    :: q_recv(inp%nvar_bcl*inp%nlayers, b%ngl, par%num_send_recv_total)

   !local variables

   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar, imulti

   integer :: k, iface, iquad, el, er, il, jl, ir, jr, I
   integer :: kl, kr, jquad, n, m, index
   real :: wq, nxl, nyl, hi
   real :: dpl, dpr, uu, vv, flux, ul, ur, vl, vr, weight
   real :: dp_lr(2,inp%nlayers), dp_deficit(2)
   real, dimension(b%nq,inp%nlayers) :: flux_edge_u, flux_edge_v
   real, dimension(inp%nvar_bcl) :: ql, qr
   real, dimension(3,b%nq) :: qbl, qbr
   real, parameter :: eps = 1.0e-10

   jj=1
   kk=1
   imm = 0

   do inbh = 1, par%num_nbh
      do ib=1,par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         el=G%face(7,iface)
         er=G%face(8,iface)

         qbl(1,:) = btp%ope_face_ave(1,:,iface)
         qbl(2,:) = btp%uvb_face_ave(1,1,:,iface)
         qbl(3,:) = btp%uvb_face_ave(2,1,:,iface)
         qbr(1,:) = btp%ope_face_ave(2,:,iface)
         qbr(2,:) = btp%uvb_face_ave(1,2,:,iface)
         qbr(3,:) = btp%uvb_face_ave(2,2,:,iface)

         do iquad = 1,b%nq

            nxl = mf%normal_vector_q(1,iquad,1,iface)
            nyl = mf%normal_vector_q(2,iquad,1,iface)

            do k = 1,inp%nlayers

               index = inp%nvar_bcl*(k-1)

               ql = 0.0; qr = 0.0
               do n = 1, b%ngl
                  hi = b%psiq(n,iquad)
                  do ivar = 1,inp%nvar_bcl
                     !Left Element
                     ql(ivar) = ql(ivar) + hi*q_send(index+ivar,n,kk)
                     !Right Element
                     qr(ivar) = qr(ivar) + hi*q_recv(index+ivar,n,kk)
                  end do
               enddo

               dpl = qbl(1,iquad) * ql(1)
               dpr = qbr(1,iquad) * qr(1)

               dp_lr(1,k) = dpl
               dp_lr(2,k) = dpr

               ul = ql(2) + qbl(2,iquad)
               ur = qr(2) + qbr(2,iquad)
               vl = ql(3) + qbl(3,iquad)
               vr = qr(3) + qbr(3,iquad)

               uu = 0.5*(ul + ur)
               vv = 0.5*(vl + vr)

               if(uu*nxl > 0.0) then
                  flux_edge_u(iquad,k) = uu * dpl
               else
                  flux_edge_u(iquad,k) = uu * dpr
               endif

               if(vv*nyl > 0.0) then
                  flux_edge_v(iquad,k) = vv * dpl
               else
                  flux_edge_v(iquad,k) = vv * dpr
               endif

            end do ! k

            dp_deficit(1) = btp%btp_mass_flux_face_ave(1,iquad,iface) - sum(flux_edge_u(iquad,:))
            dp_deficit(2) = btp%btp_mass_flux_face_ave(2,iquad,iface) - sum(flux_edge_v(iquad,:))

            do k = 1, inp%nlayers

               weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:))+eps))
               if (dp_deficit(1)*nxl < 0.0) &
                  weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:))+eps))
               flux_edge_u(iquad,k) = flux_edge_u(iquad,k) + weight*dp_deficit(1)

               weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:))+eps))
               if (dp_deficit(1)*nxl < 0.0) &
                  weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:))+eps))
               flux_edge_v(iquad,k) = flux_edge_v(iquad,k) + weight*dp_deficit(2)

            end do ! k
         end do ! iquad

         ! Do Gauss-Lobatto Integration
         do k = 1,inp%nlayers

            do iquad = 1, b%nq

               wq = mf%jac_faceq(iquad,1,iface)
               nxl = mf%normal_vector_q(1,iquad,1,iface)
               nyl = mf%normal_vector_q(2,iquad,1,iface)

               flux = nxl*flux_edge_u(iquad,k) + nyl*flux_edge_v(iquad,k)

               do n = 1, b%ngl

                  hi = b%psiq(n,iquad)
                  il = mf%imapl(1,n,1,iface)
                  jl = mf%imapl(2,n,1,iface)
                  kl = mf%imapl(3,n,1,iface)
                  I = G%intma(il,jl,kl,el)

                  rhs(I,k) = rhs(I,k) - wq*hi*flux
               end do
            end do
         end do

         kk=kk+1
         jj=jj+1
      end do

   end do !iface

end subroutine create_nbhs_face_bcl_continuity

subroutine create_nbhs_face_bcl_momentum(G, inp, b, mf, par, btp, init, rhs, q_send, q_recv, multirate)

   use mod_grid,      only: grid
   use mod_input,     only: input
   use mod_basis,     only: basis
   use mod_face,      only: face_CS
   use mod_parallel,  only: parallel_CS
   use mod_variables, only: btp_CS
   use mod_initial,   only: initial
   use mod_constants, only: gravity

   implicit none

   type(grid),        intent(in)    :: G
   type(input),       intent(in)    :: inp
   type(basis),       intent(in)    :: b
   type(face_CS),     intent(in)    :: mf
   type(parallel_CS), intent(in)    :: par
   type(btp_CS),      intent(inout) :: btp
   type(initial),     intent(in)    :: init
   integer,           intent(in)    :: multirate

   !global arrays
   real, intent(inout) :: rhs(2, G%npoin, inp%nlayers)
   real, intent(in)    :: q_send(inp%nvar_bcl*inp%nlayers, b%ngl, par%num_send_recv_total)
   real, intent(in)    :: q_recv(inp%nvar_bcl*inp%nlayers, b%ngl, par%num_send_recv_total)

   !local variables

   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar, imulti

   real, dimension(inp%nlayers) :: alpha_over_g, g_over_alpha
   real, dimension(2,inp%nlayers+1) :: p_face, z_face
   real, dimension(inp%nlayers+1) :: p_edge_plus, p_edge_minus, p2l, p2r, z_edge_plus, z_edge_minus
   real, dimension(inp%nvar_bcl,inp%nlayers) :: ql, qr
   real, dimension(3,b%nq) :: qbl, qbr
   integer :: iface, ilr, k, iquad, ktemp, I
   real :: z_intersect_top,z_intersect_bot, dz_intersect, H_r_plus, H_r_minus, acceleration
   real :: p_intersect_bot, p_intersect_top, one_plus_eta_edge
   real :: H_corr,p_inc, weight, H_corr1,p_inc1, H_corr2,p_inc2, temp, ope_l, ope_r
   integer :: Iq, el, er
   real ::  ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv
   real, dimension(inp%nlayers) :: udpl, udpr, vdpl, vdpr
   real :: uu_dp_flux_deficit(2), vv_dp_flux_deficit(2)
   real, parameter :: eps1 = 1.0e-20
   integer :: il, jl, ir, jr, kl, kr, jquad, n, m, index
   real :: wq, hi, hlx_k, hly_k, hrx_k, hry_k, flux_x, flux_y, hx_k, hy_k, flux
   real, dimension(2,b%nq,inp%nlayers) :: udp_flux, vdp_flux
   real, dimension(b%nq,inp%nlayers) :: H_face, flux_ul, flux_vl

   jj=1
   kk=1
   imm = 0

   do k=1,inp%nlayers
      alpha_over_g(k) = init%alpha_mlswe(k)/gravity
      g_over_alpha(k) = gravity/init%alpha_mlswe(k)
   enddo

   do inbh = 1, par%num_nbh
      do ib=1,par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         el=G%face(7,iface)
         er=G%face(8,iface)

         qbl(1,:) = btp%ope_face_ave(1,:,iface)
         qbl(2,:) = btp%uvb_face_ave(1,1,:,iface)
         qbl(3,:) = btp%uvb_face_ave(2,1,:,iface)
         qbr(1,:) = btp%ope_face_ave(2,:,iface)
         qbr(2,:) = btp%uvb_face_ave(1,2,:,iface)
         qbr(3,:) = btp%uvb_face_ave(2,2,:,iface)

         do iquad = 1,b%nq

            nxl = mf%normal_vector_q(1,iquad,1,iface)
            nyl = mf%normal_vector_q(2,iquad,1,iface)

            ql = 0.0; qr = 0.0

            do k = 1,inp%nlayers

               index = inp%nvar_bcl*(k-1)

               do n = 1, b%ngl
                  hi = b%psiq(n,iquad)
                  do ivar = 1,inp%nvar_bcl
                     !Left Element
                     ql(ivar,k) = ql(ivar,k) + hi*q_send(index+ivar,n,kk)
                     !Right Element
                     qr(ivar,k) = qr(ivar,k) + hi*q_recv(index+ivar,n,kk)
                  end do
               enddo

               ! Compute the fluxes
               dpl = qbl(1,iquad) * ql(1,k)
               dpr = qbr(1,iquad) * qr(1,k)
               ul = ql(2,k) + qbl(2,iquad)
               ur = qr(2,k) + qbr(2,iquad)
               vl = ql(3,k) + qbl(3,iquad)
               vr = qr(3,k) + qbr(3,iquad)

               uu = 0.5*(ul+ur)
               vv = 0.5*(vl+vr)
               udpl(k) = ul*dpl
               udpr(k) = ur*dpr
               vdpl(k) = vl*dpl
               vdpr(k) = vr*dpr

               if(uu*nxl > 0.0) then
                  udp_flux(1,iquad,k) = uu * (ul*dpl)
                  vdp_flux(1,iquad,k) = uu * (vl*dpl)
               else
                  udp_flux(1,iquad,k) = uu * (ur*dpr)
                  vdp_flux(1,iquad,k) = uu * (vr*dpr)
               endif
               if(vv*nyl > 0.0) then
                  udp_flux(2,iquad,k) = vv * (ul*dpl)
                  vdp_flux(2,iquad,k) = vv * (vl*dpl)
               else
                  udp_flux(2,iquad,k) = vv * (ur*dpr)
                  vdp_flux(2,iquad,k) = vv * (vr*dpr)
               endif

            enddo !k

            uu_dp_flux_deficit(1) = btp%Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,iquad,:))
            uu_dp_flux_deficit(2) = btp%Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,iquad,:))
            vv_dp_flux_deficit(1) = btp%Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,iquad,:))
            vv_dp_flux_deficit(2) = btp%Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,iquad,:))

            do k = 1,inp%nlayers

               ! Adjust the fluxes for the u-momentum equation
               !x-direction
               weight = abs(udpl(k)) / (sum(abs(udpl(:))+eps1))
               if(uu_dp_flux_deficit(1)*nxl < 0.0) &
                  weight = abs(udpr(k)) / (sum(abs(udpr(:))+eps1))
               udp_flux(1,iquad,k) = udp_flux(1,iquad,k) + weight * uu_dp_flux_deficit(1)

               !y-direction
               weight = abs(udpl(k)) / (sum(abs(udpl(:))+eps1))
               if(uu_dp_flux_deficit(2)*nyl < 0.0) &
                  weight = abs(udpr(k)) / (sum(abs(udpr(:))+eps1))
               udp_flux(2,iquad,k) = udp_flux(2,iquad,k) + weight * uu_dp_flux_deficit(2)

               ! Adjust the fluxes for the v-momentum equation
               !x-direction
               weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
               if(vv_dp_flux_deficit(1)*nxl < 0.0) &
                  weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
               vdp_flux(1,iquad,k) = vdp_flux(1,iquad,k) + weight * vv_dp_flux_deficit(1)

               !y-direction
               weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
               if(vv_dp_flux_deficit(2)*nyl < 0.0) &
                  weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
               vdp_flux(2,iquad,k) = vdp_flux(2,iquad,k) + weight * vv_dp_flux_deficit(2)
            end do

            z_face = 0.0 ; p_face = 0.0
            z_edge_plus = 0.0 ; z_edge_minus = 0.0
            p_edge_plus = 0.0; p_edge_minus = 0.0

            !Store Left Side Variables
            ope_l = sqrt(btp%ope2_face_ave(1,iquad,iface))
            ope_r = sqrt(btp%ope2_face_ave(2,iquad,iface))
            p_face(1,1) = 0.0
            p_face(2,1) = 0.0
            do k=1,inp%nlayers
               p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,k)
               p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,k)
            end do

            one_plus_eta_edge = sqrt(btp%one_plus_eta_edge_2_ave(iquad,iface))
            z_face(1,inp%nlayers+1) = init%zbot_face(1,iquad,iface)
            z_face(2,inp%nlayers+1) = init%zbot_face(2,iquad,iface)
            z_edge_plus(inp%nlayers+1) = init%zbot_face(1,iquad,iface)
            z_edge_minus(inp%nlayers+1) = init%zbot_face(2,iquad,iface)
            do k=inp%nlayers,1,-1
               z_face(1,k) = z_face(1,k+1) + alpha_over_g(k) * (ope_l * ql(1,k))
               z_face(2,k) = z_face(2,k+1) + alpha_over_g(k) * (ope_r * qr(1,k))
               z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k) * &
                  (one_plus_eta_edge * ql(1,k))
               z_edge_minus(k) = z_edge_minus(k+1) + alpha_over_g(k) * &
                  (one_plus_eta_edge * qr(1,k))
            end do

            p_edge_plus(2) = one_plus_eta_edge * ql(1,1)
            p_edge_minus(2) = one_plus_eta_edge * qr(1,1)
            do k = 2,inp%nlayers
               p_edge_plus(k+1) = p_edge_plus(k) + one_plus_eta_edge * ql(1,k)
               p_edge_minus(k+1) = p_edge_minus(k) + one_plus_eta_edge * qr(1,k)
            end do

            do k = 1, inp%nlayers

               ! Computation from + side for layer k
               H_r_plus = 0.5*init%alpha_mlswe(k)*(p_edge_plus(k+1)**2 - p_edge_plus(k)**2)

               ! Computation from - side for layer k
               H_r_minus = 0.0
               do ktemp = 1, inp%nlayers

                  z_intersect_top = min(z_edge_minus(ktemp), z_edge_plus(k))
                  z_intersect_bot = max(z_edge_minus(ktemp+1), z_edge_plus(k+1))
                  dz_intersect = z_intersect_top - z_intersect_bot

                  if (dz_intersect > 0.0) then
                     p_intersect_bot = p_edge_minus(ktemp+1) &
                        - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_minus(ktemp+1))
                     p_intersect_top = p_edge_minus(ktemp+1) &
                        - g_over_alpha(ktemp)*(z_intersect_top - z_edge_minus(ktemp+1))
                     H_r_minus = H_r_minus + &
                        0.5*init%alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)

                  end if
               end do
               H_face(iquad,k) = 0.5*(H_r_plus + H_r_minus) !computation of H_r for the left side

            end do !k

            do k = 1, inp%nlayers-1          ! interface at the bottom of layer k
               ! Corrections at the left side of a face.
               p_inc1 = g_over_alpha(k)*(z_face(1,k+1) - z_edge_plus(k+1))
               H_corr1 = 0.5 * init%alpha_mlswe(k) * ((p_face(1,k+1) &
                  + p_inc1)**2 - p_face(1,k+1)**2)
               H_face(iquad,k) = H_face(iquad,k) - H_corr1
               H_face(iquad,k+1) = H_face(iquad,k+1) + H_corr1

            end do

            ! Left side of face
            weight = 1.0
            acceleration = sum(H_face(iquad,:))
            if(acceleration > 0.0) then
               weight = btp%H_face_ave(iquad,iface) / acceleration
            end if
            H_face(iquad,:) = H_face(iquad,:) * weight

            do k = 1,inp%nlayers
               flux_ul(iquad,k) = nxl*(udp_flux(1,iquad,k) + H_face(iquad,k)) + nyl*udp_flux(2,iquad,k)
               flux_vl(iquad,k) = nxl*vdp_flux(1,iquad,k) + nyl*(vdp_flux(2,iquad,k) + H_face(iquad,k))
            enddo

            uu_dp_flux_deficit(1) = btp%Qu_face_ave(1,iquad,iface) - sum(flux_ul(iquad,:))
            vv_dp_flux_deficit(1) = btp%Qv_face_ave(1,iquad,iface) - sum(flux_vl(iquad,:))

         end do ! iquad

         ! Do Gauss-Lobatto Integration
         do k = 1,inp%nlayers

            do iquad = 1, b%nq

               wq = mf%jac_faceq(iquad,1,iface)

               flux_x = flux_ul(iquad,k)
               flux_y = flux_vl(iquad,k)

               do n = 1, b%ngl

                  hi = b%psiq(n,iquad)
                  il = mf%imapl(1,n,1,iface)
                  jl = mf%imapl(2,n,1,iface)
                  kl = mf%imapl(3,n,1,iface)
                  I = G%intma(il,jl,kl,el)

                  rhs(1,I,k) = rhs(1,I,k) - wq*hi*flux_x
                  rhs(2,I,k) = rhs(2,I,k) - wq*hi*flux_y
               end do
            end do
         end do

         kk=kk+1
         jj=jj+1
      end do

   end do !iface

end subroutine create_nbhs_face_bcl_momentum

subroutine create_nbhs_face_df_lap_bcl(G, b, mf, par, btp, ref, rhs, q_send, q_recv, nlayers, multirate)

   use mod_grid,      only: grid
   use mod_basis,     only: basis
   use mod_face,      only: face_CS
   use mod_parallel,  only: parallel_CS
   use mod_variables, only: btp_CS
   use mod_ref,       only: mref

   implicit none

   type(grid),        intent(in)    :: G
   type(basis),       intent(in)    :: b
   type(face_CS),     intent(in)    :: mf
   type(parallel_CS), intent(in)    :: par
   type(btp_CS),      intent(inout) :: btp
   type(mref),        intent(in)    :: ref
   integer,           intent(in)    :: nlayers
   integer,           intent(in)    :: multirate

   !global arrays
   real, intent(inout) :: rhs(2, G%npoin, nlayers)
   real, intent(in)    :: q_send(5*nlayers, b%ngl, par%num_send_recv_total)
   real, intent(in)    :: q_recv(5*nlayers, b%ngl, par%num_send_recv_total)

   integer :: kk, k, iquad, i, ivar, iface, el, il, jl, kl, ip, index
   integer :: ngl_f, nlayers_f, nboun_valid_f
   real    :: wq, nxl, nyl, hi, flux_qu, flux_qv, alpha, beta, iflux
   real, dimension(4,2) :: flux_uv_visc_face
   real, dimension(2)   :: qu_mean, qv_mean, qul, qur, qvl, qvr

   ngl_f         = b%ngl
   nlayers_f     = nlayers
   nboun_valid_f = ref%nboun_valid
   beta  = 0.5
   alpha = 1.0 - beta
   iflux = 0.0

   !$acc parallel loop gang &
   !$acc    present(ref%face_pack_list, G%face, G%intma, mf%imapl, mf%jac_face, &
   !$acc            mf%normal_vector, b%psi, btp%graduvb_face_ave, q_send, q_recv, rhs) &
   !$acc    firstprivate(ngl_f, nlayers_f, nboun_valid_f, alpha, beta, iflux) &
   !$acc    private(flux_uv_visc_face, qu_mean, qv_mean, qul, qur, qvl, qvr)
   do kk = 1, nboun_valid_f
      iface = ref%face_pack_list(kk)
      el    = G%face(7, iface)

      !$acc loop seq
      do k = 1, nlayers_f
         index = 5*(k-1)

         !$acc loop seq
         do iquad = 1, ngl_f
            do ivar = 1, 4
               flux_uv_visc_face(ivar,1) = q_send(index+5,iquad,kk)* &
                  btp%graduvb_face_ave(ivar,1,iquad,iface) + q_send(index+ivar,iquad,kk)
               flux_uv_visc_face(ivar,2) = q_recv(index+5,iquad,kk)* &
                  btp%graduvb_face_ave(ivar,2,iquad,iface) + q_recv(index+ivar,iquad,kk)
            end do

            nxl = mf%normal_vector(1,iquad,1,iface)
            nyl = mf%normal_vector(2,iquad,1,iface)

            qul(1) = flux_uv_visc_face(1,1); qul(2) = flux_uv_visc_face(2,1)
            qvl(1) = flux_uv_visc_face(3,1); qvl(2) = flux_uv_visc_face(4,1)
            qur(1) = flux_uv_visc_face(1,2); qur(2) = flux_uv_visc_face(2,2)
            qvr(1) = flux_uv_visc_face(3,2); qvr(2) = flux_uv_visc_face(4,2)

            qu_mean(1) = alpha*qul(1) + beta*qur(1)
            qu_mean(2) = alpha*qul(2) + beta*qur(2)
            qv_mean(1) = alpha*qvl(1) + beta*qvr(1)
            qv_mean(2) = alpha*qvl(2) + beta*qvr(2)

            wq = mf%jac_face(iquad,1,iface)

            flux_qu = (qu_mean(1) - iflux*qul(1))*nxl + (qu_mean(2) - iflux*qul(2))*nyl
            flux_qv = (qv_mean(1) - iflux*qvl(1))*nxl + (qv_mean(2) - iflux*qvl(2))*nyl

            do i = 1, ngl_f
               hi = b%psi(i, iquad)
               il = mf%imapl(1,i,1,iface)
               jl = mf%imapl(2,i,1,iface)
               kl = mf%imapl(3,i,1,iface)
               ip = G%intma(il,jl,kl,el)
               !$acc atomic update
               rhs(1,ip,k) = rhs(1,ip,k) + wq*hi*flux_qu
               !$acc atomic update
               rhs(2,ip,k) = rhs(2,ip,k) + wq*hi*flux_qv
            end do
         end do
      end do
   end do
   !$acc end parallel loop

end subroutine create_nbhs_face_df_lap_bcl

subroutine create_nbhs_face_quad_layer(G, b, par, q_face, q_send, q_recv, nvarb, nlayers, nq)

   use mod_grid,     only: grid, mod_grid_get_face_nq
   use mod_basis,    only: basis
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in)    :: G
   type(basis),       intent(in)    :: b
   type(parallel_CS), intent(in)    :: par

   !global arrays
   integer, intent(in) :: nvarb, nlayers, nq
   real, intent(inout) :: q_face(nvarb, 2, nq, G%nface, nlayers)
   real, intent(in)    :: q_send(nvarb, nq, par%num_send_recv_total, nlayers)
   real, intent(in)    :: q_recv(nvarb, nq, par%num_send_recv_total, nlayers)

   !local variables

   integer ::iface,k
   integer :: ilocl, nq_i, nq_j, plane_ij

   integer :: jj, imm, inbh, ib, kk
   integer :: ftype, imulti, ll, inode, iel

   jj=1
   kk=1
   imm = 0

   do inbh = 1, par%num_nbh
      do ib=1,par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         ftype = G%face_type(iface)

         if (ftype==2) then
            ilocl=G%face(5,iface)
            iel=G%face(7,iface)

            !-------------------------------------
            !Store Left Side Variables
            !-------------------------------------
            call mod_grid_get_face_nq(b, ilocl, nq_i, nq_j, plane_ij)

            do ll = 1,nlayers
               do inode = 1,nq
                  do k=1,nvarb
                     !Left Element
                     q_face(k,1,inode,iface,ll)=q_send(k,inode,kk,ll)

                     !Right Element
                     q_face(k,2,inode,iface,ll)=q_recv(k,inode,jj,ll)
                  end do
               end do
            end do

            kk=kk+1
         end if
         jj=jj+1
      end do

   end do !iface

end subroutine create_nbhs_face_quad_layer
