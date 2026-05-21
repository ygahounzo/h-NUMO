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
!----------------------------------------------------------------------!

 subroutine create_nbhs_face_df(G, inp, b, mf, par, btp, init, rhs, q_send, q_recv, nvarb)

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

   !global arrays
   integer, intent(in) :: nvarb
   real, intent(in)    :: q_send(nvarb, b%ngl, G%nboun)
   real, intent(in)    :: q_recv(nvarb, b%ngl, G%nboun)
   real, intent(inout) :: rhs(3, G%npoin)
 
   !local variables
 
   integer :: jj, imm, inbh, ib, kk, ii
   integer :: imulti
   integer :: iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n
   real :: wq, hi, nxl, nyl, nxr, nyr, un
   real :: ul, ur, vl, vr, pbl, pbr, clam, one_eta
   real :: pU_L, pU_R, pbpert_edge
   real :: qbl(4), qbr(4), flux(3,b%nq)
   real :: c_minus, c_plus
   real :: ppl, ppr, upl, upr, vpl, vpr
   real, dimension(inp%nlayers+1) :: pprime_l, pprime_r
   integer :: itype, k
   real :: H_bcl_ql, H_bcl_qr, H_bcl_q, flux_pb, flux_u, flux_v, fxl, fxr, flux_edge_x, flux_edge_y
   real, dimension(2) :: Qu_ql, Qu_qr, Qv_ql, Qv_qr

   jj=1
   kk=1
   imm = 0

   do inbh = 1, par%num_nbh
      do ib=1,par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         el = G%face(7,iface)
         er = G%face(8,iface)
 
         !-------------------------------------
         !Store Left Side Variables
         !-------------------------------------

         do iquad = 1, b%nq

            nxl = mf%normal_vector_q(1,iquad,1,iface)
            nyl = mf%normal_vector_q(2,iquad,1,iface)

            nxr = -nxl; nyr = -nyl

            qbl = 0.0; qbr = 0.0 ; pbl = 0.0; pbr = 0.0
            do n = 1,b%ngl
               hi = b%psiq(n,iquad)
               !Left Element
               qbl(1:4) = qbl(1:4) + hi*q_send(1:4,n,kk)
               pbl = pbl + hi*q_send(5,n,kk)
               !Right Element
               qbr(1:4) = qbr(1:4) + hi*q_recv(1:4,n,kk)
               pbr = pbr + hi*q_recv(5,n,kk)
            end do

            pU_L = nxl * qbl(3) + nyl * qbl(4)
            pU_R = nxr * qbr(3) + nyr * qbr(4)

            c_minus = sqrt(init%alpha_mlswe(inp%nlayers) * pbr)
            c_plus  = sqrt(init%alpha_mlswe(inp%nlayers) * pbl)
            clam = max(c_minus, c_plus)

            pbpert_edge = 0.5*(qbl(2) + qbr(2)) + (0.5/clam) * (pU_L + pU_R)
            one_eta = 1.0 + (pbpert_edge/pbl)

            ul = qbl(3)/qbl(1); ur = qbr(3)/qbr(1)
            vl = qbl(4)/qbl(1); vr = qbr(4)/qbr(1)

            Qu_ql(1) = ul*qbl(3) ; Qu_ql(2) = vl*qbl(3)
            Qu_qr(1) = ur*qbr(3) ; Qu_qr(2) = vr*qbr(3)

            Qv_ql(1) = ul*qbl(4) ; Qv_ql(2) = vl*qbl(4)
            Qv_qr(1) = ur*qbr(4) ; Qv_qr(2) = vr*qbr(4)

            pprime_l(:) = 0.0; pprime_r(:) = 0.0 ; H_bcl_q = 0.0
            H_bcl_ql = 0.0 ; H_bcl_qr = 0.0
            do k = 1, inp%nlayers
               ppl = 0.0; ppr = 0.0 ; upl = 0.0; upr = 0.0 ; vpl = 0.0; vpr = 0.0

               ii = 5 + (k-1)*3
               do n = 1, b%ngl
                  il = mf%imapl(1,n,1,iface)
                  jl = mf%imapl(2,n,1,iface)
                  kl = mf%imapl(3,n,1,iface)
                  I = G%intma(il,jl,kl,el)

                  hi = b%psiq(n,iquad)
                  ppl = ppl + hi*q_send(ii+1,n,kk)
                  ppr = ppr + hi*q_recv(ii+1,n,kk)

                  upl = upl + hi*q_send(ii+2,n,kk)
                  upr = upr + hi*q_recv(ii+2,n,kk)

                  vpl = vpl + hi*q_send(ii+3,n,kk)
                  vpr = vpr + hi*q_recv(ii+3,n,kk)
               enddo

               Qu_ql(1) = Qu_ql(1) + (upl*(upl*(one_eta*ppl)))
               Qu_qr(1) = Qu_qr(1) + (upr*(upr*(one_eta*ppr)))
               Qu_ql(2) = Qu_ql(2) + (vpl*(upl*(one_eta*ppl)))
               Qu_qr(2) = Qu_qr(2) + (vpr*(upr*(one_eta*ppr)))

               Qv_ql(1) = Qv_ql(1) + (upl*(vpl*(one_eta*ppl)))
               Qv_qr(1) = Qv_qr(1) + (upr*(vpr*(one_eta*ppr)))
               Qv_ql(2) = Qv_ql(2) + (vpl*(vpl*(one_eta*ppl)))
               Qv_qr(2) = Qv_qr(2) + (vpr*(vpr*(one_eta*ppr)))

               pprime_l(k+1) = pprime_l(k) + ppl
               pprime_r(k+1) = pprime_r(k) + ppr

               H_bcl_ql  = H_bcl_ql + 0.5*init%alpha_mlswe(k)*(pprime_l(k+1)**2 - pprime_l(k)**2)
               H_bcl_qr  = H_bcl_qr + 0.5*init%alpha_mlswe(k)*(pprime_r(k+1)**2 - pprime_r(k)**2)
            end do

            ! Compute mass fluxes at each element face.

            flux_edge_x = 0.5*(qbl(3) + qbr(3)) + (0.5*clam)*(nxl * qbl(2) + nxr * qbr(2))
            flux_edge_y = 0.5*(qbl(4) + qbr(4)) + (0.5*clam)*(nyl * qbl(2) + nyr * qbr(2))
            flux(1,iquad) = nxl*flux_edge_x + nyl*flux_edge_y

            ! Compute pressure forcing H_face at each element face.
            H_bcl_q = (one_eta**2) * (0.5 * (H_bcl_ql + H_bcl_qr))

            ! Accumulate sums for time averaging

            btp%btp_mass_flux_face_ave(1,iquad,iface) = btp%btp_mass_flux_face_ave(1,iquad,iface) + flux_edge_x
            btp%btp_mass_flux_face_ave(2,iquad,iface) = btp%btp_mass_flux_face_ave(2,iquad,iface) + flux_edge_y

            btp%H_face_ave(iquad,iface) = btp%H_face_ave(iquad,iface) + H_bcl_q
            btp%Qu_face_ave(1,iquad,iface) = btp%Qu_face_ave(1,iquad,iface) + 0.5*(Qu_ql(1) + Qu_qr(1))
            btp%Qu_face_ave(2,iquad,iface) = btp%Qu_face_ave(2,iquad,iface) + 0.5*(Qu_ql(2) + Qu_qr(2))
            btp%Qv_face_ave(1,iquad,iface) = btp%Qv_face_ave(1,iquad,iface) + 0.5*(Qv_ql(1) + Qv_qr(1))
            btp%Qv_face_ave(2,iquad,iface) = btp%Qv_face_ave(2,iquad,iface) + 0.5*(Qv_ql(2) + Qv_qr(2))
            btp%ope_face_ave(1,iquad,iface) = btp%ope_face_ave(1,iquad,iface) + (1.0 + (qbl(2)/pbl))
            btp%ope_face_ave(2,iquad,iface) = btp%ope_face_ave(2,iquad,iface) + (1.0 + (qbr(2)/pbr))
            btp%ope2_face_ave(1,iquad,iface) = btp%ope2_face_ave(1,iquad,iface) + (1.0 + (qbl(2)/pbl))**2
            btp%ope2_face_ave(2,iquad,iface) = btp%ope2_face_ave(2,iquad,iface) + (1.0 + (qbr(2)/pbr))**2
            btp%one_plus_eta_edge_2_ave(iquad,iface) = btp%one_plus_eta_edge_2_ave(iquad,iface) &
                                                + one_eta**2
            btp%uvb_face_ave(1,1,iquad,iface) = btp%uvb_face_ave(1,1,iquad,iface) + ul
            btp%uvb_face_ave(1,2,iquad,iface) = btp%uvb_face_ave(1,2,iquad,iface) + ur
            btp%uvb_face_ave(2,1,iquad,iface) = btp%uvb_face_ave(2,1,iquad,iface) + vl
            btp%uvb_face_ave(2,2,iquad,iface) = btp%uvb_face_ave(2,2,iquad,iface) + vr

            Qu_ql(1) = Qu_ql(1) + (one_eta**2)*H_bcl_ql
            Qu_qr(1) = Qu_qr(1) + (one_eta**2)*H_bcl_qr

            fxl = nxl*Qu_ql(1) + nyl*Qu_ql(2)
            fxr = nxr*Qu_qr(1) + nyr*Qu_qr(2)

            flux(2,iquad) = 0.5*(fxl - fxr) - (0.25*clam)*(qbr(3) - qbl(3))

            Qv_ql(2) = Qv_ql(2) + (one_eta**2)*H_bcl_ql
            Qv_qr(2) = Qv_qr(2) + (one_eta**2)*H_bcl_qr

            fxl = nxl*Qv_ql(1) + nyl*Qv_ql(2)
            fxr = nxr*Qv_qr(1) + nyr*Qv_qr(2)

            flux(3,iquad) = 0.5*(fxl - fxr) - (0.25*clam)*(qbr(4) - qbl(4))

         end do !iquad

         do iquad = 1, b%nq

            wq = mf%jac_faceq(iquad,1,iface)

            flux_pb = flux(1,iquad)
            flux_u = flux(2,iquad)
            flux_v = flux(3,iquad)

            do n = 1, b%ngl

               hi = b%psiq(n,iquad)
               il = mf%imapl(1,n,1,iface)
               jl = mf%imapl(2,n,1,iface)
               kl = mf%imapl(3,n,1,iface)
               I = G%intma(il,jl,kl,el)

               rhs(1,I) = rhs(1,I) - wq*hi*flux_pb
               rhs(2,I) = rhs(2,I) - wq*hi*flux_u
               rhs(3,I) = rhs(3,I) - wq*hi*flux_v
            end do
         end do !iquad

         kk=kk+1
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_df

 subroutine create_nbhs_face_df_lap(G, b, mf, par, btp, rhs, q_send, q_recv, nvarb)

   use mod_grid,      only: grid
   use mod_basis,     only: basis
   use mod_face,      only: face_CS
   use mod_parallel,  only: parallel_CS
   use mod_variables, only: btp_CS

   implicit none

   type(grid),        intent(in)    :: G
   type(basis),       intent(in)    :: b
   type(face_CS),     intent(in)    :: mf
   type(parallel_CS), intent(in)    :: par
   type(btp_CS),      intent(inout) :: btp

   !global arrays
   integer, intent(in) :: nvarb
   real, intent(inout) :: rhs(2, G%npoin)
   real, intent(in)    :: q_send(10, b%ngl, G%nboun)
   real, intent(in)    :: q_recv(10, b%ngl, G%nboun)
 
   !local variables
 
   real :: wq, iflux
   integer ::iface,k
   integer :: iel, ier, ilocl, ilocr
   integer :: ifaceb, nq_i, nq_j, plane_ij
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
   integer :: ftype, pface, subface, imulti
   real :: nxl, nyl, nxr, nyr
   integer :: el, iquad, n, I, il, jl, kl, er, itype, ip

   real, dimension(2) :: qu_mean, qv_mean
   real, dimension(4,2) :: flux_uv_visc_face
   real, dimension(2) :: qul,qur
   real, dimension(2) :: qvl,qvr
   real, dimension(5) :: grad_uvb_pb_l, grad_uvb_pb_r
   real :: flux_qu, flux_qv, hi, mul, mur,c_jump, alpha, beta

   jj=1
   kk=1
   imm = 0

   beta = 0.5
   alpha = 1.0 - beta
   iflux = 0.0

   do inbh = 1, par%num_nbh
      do ib=1,par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         el=G%face(7,iface)
         er=G%face(8,iface)
         
         do iquad = 1,b%ngl

            do ivar = 6, 10
               grad_uvb_pb_l(ivar-5) =  q_send(ivar,iquad,kk)
               grad_uvb_pb_r(ivar-5) =  q_recv(ivar,iquad,kk)
            end do

            do ivar = 1,4
               flux_uv_visc_face(ivar,1) = grad_uvb_pb_l(5)* q_send(ivar,iquad,kk) + &
                                          grad_uvb_pb_l(ivar)

               flux_uv_visc_face(ivar,2) = grad_uvb_pb_r(5)* q_recv(ivar,iquad,kk) + &
                                          grad_uvb_pb_r(ivar)

               btp%graduvb_face_ave(ivar,1,iquad,iface) = btp%graduvb_face_ave(ivar,1,iquad,iface) + q_send(ivar,iquad,kk)
               btp%graduvb_face_ave(ivar,2,iquad,iface) = btp%graduvb_face_ave(ivar,2,iquad,iface) + q_recv(ivar,iquad,kk)

            end do

            nxl = mf%normal_vector(1,iquad,1,iface)
            nyl = mf%normal_vector(2,iquad,1,iface)

            qul(1) = flux_uv_visc_face(1,1)
            qul(2) = flux_uv_visc_face(2,1)
            qvl(1) = flux_uv_visc_face(3,1)
            qvl(2) = flux_uv_visc_face(4,1)

            qur(1) = flux_uv_visc_face(1,2)
            qur(2) = flux_uv_visc_face(2,2)
            qvr(1) = flux_uv_visc_face(3,2)
            qvr(2) = flux_uv_visc_face(4,2)

            ! The Flip-Flop flux of Cockburn & Shu 
            ! NOTE: beta=0.5 is the central flux
            qu_mean(1) = alpha*qul(1) + beta*qur(1)
            qu_mean(2) = alpha*qul(2) + beta*qur(2)
            qv_mean(1) = alpha*qvl(1) + beta*qvr(1)
            qv_mean(2) = alpha*qvl(2) + beta*qvr(2)

            wq=mf%jac_face(iquad,1,iface)

            flux_qu = (qu_mean(1) - iflux*qul(1))*nxl + (qu_mean(2) - iflux*qul(2))*nyl
            flux_qv = (qv_mean(1) - iflux*qvl(1))*nxl + (qv_mean(2) - iflux*qvl(2))*nyl

            !  Do Gauss-Lobatto Integration
            do i=1,b%ngl

               hi = b%psi(i,iquad)

               il=mf%imapl(1,i,1,iface)
               jl=mf%imapl(2,i,1,iface)
               kl=mf%imapl(3,i,1,iface)
               ip=G%intma(il,jl,kl,el)
               
               !Update Flux
               rhs(1,ip) = rhs(1,ip) + wq*hi*flux_qu
               rhs(2,ip) = rhs(2,ip) + wq*hi*flux_qv

            end do !i
         end do !iquad

         kk=kk+1
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_df_lap

 subroutine create_nbhs_face_bcl(G, inp, b, mf, par, btp, init, rhs, q_send, q_recv)

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

   !global arrays
   real, intent(inout) :: rhs(3, G%npoin, inp%nlayers)
   real, intent(in)    :: q_send(3*inp%nlayers, b%ngl, G%nboun)
   real, intent(in)    :: q_recv(3*inp%nlayers, b%ngl, G%nboun)
 
   !local variables
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar, imulti
   
   real, dimension(inp%nlayers) :: alpha_over_g, g_over_alpha
   real, dimension(2,inp%nlayers+1) :: p_face, z_face
   real, dimension(inp%nlayers+1) :: p_edge_plus, p_edge_minus, p2l, p2r, z_edge_plus, z_edge_minus
   real, dimension(3,inp%nlayers) :: ql, qr
   real, dimension(3,b%nq) :: qbl, qbr
   integer :: iface, ilr, k, iquad, ktemp, I
   real :: z_intersect_top,z_intersect_bot, dz_intersect, H_r_plus, H_r_minus, acceleration
   real :: p_intersect_bot, p_intersect_top, one_plus_eta_edge
   real :: H_corr,p_inc, weight, H_corr1,p_inc1, H_corr2,p_inc2, temp, ope_l, ope_r
   integer :: Iq, el, er
   real ::  ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv
   real, dimension(inp%nlayers) :: udpl, udpr, vdpl, vdpr
   real :: uu_dp_flux_deficit(2), vv_dp_flux_deficit(2)
   real, parameter :: eps1 = 1.0e-20 !  Parameter used to prevent division by zero.
   integer :: il, jl, ir, jr, kl, kr, jquad, n, m, index
   real :: wq, hi, hlx_k, hly_k, hrx_k, hry_k, flux_x, flux_y, hx_k, hy_k, flux
   real, dimension(2,b%nq,inp%nlayers) :: udp_flux, vdp_flux, dp_flux
   real :: dp_lr(2,inp%nlayers), dp_deficit(2)
   real, dimension(b%nq,inp%nlayers) :: H_face, flux_ul, flux_vl, flux_dp

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

               index = 3*(k-1)

               do n = 1, b%ngl
                  hi = b%psiq(n,iquad)
                  do ivar = 1,3
                     !Left Element
                     ql(ivar,k) = ql(ivar,k) + hi*q_send(index+ivar,n,kk)
                     !Right Element
                     qr(ivar,k) = qr(ivar,k) + hi*q_recv(index+ivar,n,kk)
                  end do
               enddo

               ! Left side of the edge
               dpl = qbl(1,iquad) * ql(1,k)
               dpr = qbr(1,iquad) * qr(1,k)
               dp_lr(1,k) = dpl
               dp_lr(2,k) = dpr

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
                  dp_flux(1,iquad,k) = uu * dpl
                  udp_flux(1,iquad,k) = uu * (ul*dpl)
                  vdp_flux(1,iquad,k) = uu * (vl*dpl)
               else
                  dp_flux(1,iquad,k) = uu * dpr
                  udp_flux(1,iquad,k) = uu * (ur*dpr)
                  vdp_flux(1,iquad,k) = uu * (vr*dpr)
               endif
               if(vv*nyl > 0.0) then
                  dp_flux(2,iquad,k) = vv * dpl
                  udp_flux(2,iquad,k) = vv * (ul*dpl)
                  vdp_flux(2,iquad,k) = vv * (vl*dpl)
               else
                  dp_flux(2,iquad,k) = vv * dpr
                  udp_flux(2,iquad,k) = vv * (ur*dpr)
                  vdp_flux(2,iquad,k) = vv * (vr*dpr)
               endif

            enddo !k

            uu_dp_flux_deficit(1) = btp%Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,iquad,:))
            uu_dp_flux_deficit(2) = btp%Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,iquad,:))
            vv_dp_flux_deficit(1) = btp%Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,iquad,:))
            vv_dp_flux_deficit(2) = btp%Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,iquad,:))

            dp_deficit(1) = btp%btp_mass_flux_face_ave(1,iquad,iface) - sum(dp_flux(1,iquad,:))
            dp_deficit(2) = btp%btp_mass_flux_face_ave(2,iquad,iface) - sum(dp_flux(2,iquad,:))

            do k = 1,inp%nlayers

               weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:))+eps1))
               if (dp_deficit(1)*nxl < 0.0) &
                  weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:))+eps1))
               dp_flux(1,iquad,k) = dp_flux(1,iquad,k) + weight * dp_deficit(1)

               weight = dp_lr(1,k) / (sum(abs(dp_lr(1,:))+eps1))
               if (dp_deficit(2)*nyl < 0.0) &
                  weight = dp_lr(2,k) / (sum(abs(dp_lr(2,:))+eps1))
               dp_flux(2,iquad,k) = dp_flux(2,iquad,k) + weight * dp_deficit(2)

               ! Adjust the fluxes for the u-momentum equation
               ! x-direction
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
               flux_dp(iquad,k) = nxl*dp_flux(1,iquad,k) + nyl*dp_flux(2,iquad,k)
               flux_ul(iquad,k) = nxl*(udp_flux(1,iquad,k) + H_face(iquad,k)) + nyl*udp_flux(2,iquad,k)
               flux_vl(iquad,k) = nxl*vdp_flux(1,iquad,k) + nyl*(vdp_flux(2,iquad,k) + H_face(iquad,k))
            enddo

         end do ! iquad
         
         ! Do Gauss-Lobatto Integration
         do k = 1,inp%nlayers

            do iquad = 1, b%nq

               wq = mf%jac_faceq(iquad,1,iface)

               flux = flux_dp(iquad,k)
               flux_x = flux_ul(iquad,k)
               flux_y = flux_vl(iquad,k)

               do n = 1, b%ngl

                  hi = b%psiq(n,iquad)
                  il = mf%imapl(1,n,1,iface)
                  jl = mf%imapl(2,n,1,iface)
                  kl = mf%imapl(3,n,1,iface)
                  I = G%intma(il,jl,kl,el)

                  rhs(1,I,k) = rhs(1,I,k) - wq*hi*flux
                  rhs(2,I,k) = rhs(2,I,k) - wq*hi*flux_x
                  rhs(3,I,k) = rhs(3,I,k) - wq*hi*flux_y
               end do
            end do
         end do

         kk=kk+1
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_bcl

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
   real, intent(in)    :: q_send(3*inp%nlayers, b%ngl, G%nboun)
   real, intent(in)    :: q_recv(3*inp%nlayers, b%ngl, G%nboun)

   !local variables

   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar, imulti

   integer :: k, iface, iquad, el, er, il, jl, ir, jr, I
   integer :: kl, kr, jquad, n, m, index
   real :: wq, nxl, nyl, hi
   real :: dpl, dpr, uu, vv, flux, ul, ur, vl, vr, weight
   real :: dp_lr(2,inp%nlayers), dp_deficit(2)
   real, dimension(b%nq,inp%nlayers) :: flux_edge_u, flux_edge_v
   real, dimension(3) :: ql, qr
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

               index = 3*(k-1)

               ql = 0.0; qr = 0.0
               do n = 1, b%ngl
                  hi = b%psiq(n,iquad)
                  do ivar = 1,3
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
   real, intent(in)    :: q_send(3*inp%nlayers, b%ngl, G%nboun)
   real, intent(in)    :: q_recv(3*inp%nlayers, b%ngl, G%nboun)

   !local variables

   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar, imulti

   real, dimension(inp%nlayers) :: alpha_over_g, g_over_alpha
   real, dimension(2,inp%nlayers+1) :: p_face, z_face
   real, dimension(inp%nlayers+1) :: p_edge_plus, p_edge_minus, p2l, p2r, z_edge_plus, z_edge_minus
   real, dimension(3,inp%nlayers) :: ql, qr
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

               index = 3*(k-1)

               do n = 1, b%ngl
                  hi = b%psiq(n,iquad)
                  do ivar = 1,3
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

  subroutine create_nbhs_face_df_lap_bcl(G, b, mf, par, btp, rhs, q_send, q_recv, nlayers, multirate)

   use mod_grid,      only: grid
   use mod_basis,     only: basis
   use mod_face,      only: face_CS
   use mod_parallel,  only: parallel_CS
   use mod_variables, only: btp_CS

   implicit none

   type(grid),        intent(in)    :: G
   type(basis),       intent(in)    :: b
   type(face_CS),     intent(in)    :: mf
   type(parallel_CS), intent(in)    :: par
   type(btp_CS),      intent(inout) :: btp
   integer,           intent(in)    :: nlayers
   integer,           intent(in)    :: multirate

   !global arrays
   real, intent(inout) :: rhs(2, G%npoin, nlayers)
   real, intent(in)    :: q_send(5*nlayers, b%ngl, G%nboun)
   real, intent(in)    :: q_recv(5*nlayers, b%ngl, G%nboun)

   !local variables

   real :: wq, a_constant, iflux
   integer ::iface,k
   integer :: iel, ier, ilocl, ilocr
   integer :: ifaceb, nq_i, nq_j, plane_ij

   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
   integer :: ftype, pface, subface, imulti
   real :: nxl, nyl, nxr, nyr
   integer :: el, iquad, n, I, il, jl, kl, er, itype, ip, index

   real, dimension(2) :: qu_mean, qv_mean
   real, dimension(4,2) :: flux_uv_visc_face
   real, dimension(2) :: qul,qur
   real, dimension(2) :: qvl,qvr
   real, dimension(5) :: grad_uvb_pb_l, grad_uvb_pb_r
   real :: flux_qu, flux_qv, hi, mul, mur,c_jump, alpha, beta

   jj=1
   kk=1
   imm = 0

   beta = 0.5
   alpha = 1.0 - beta
   iflux = 0.0

   do inbh = 1, par%num_nbh
      do ib=1,par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         el=G%face(7,iface)
         er=G%face(8,iface)

         do k = 1, nlayers

            index = 5*(k-1)

            do iquad = 1,b%ngl

               do ivar = 1,4
                  flux_uv_visc_face(ivar,1) = q_send(index+5,iquad,kk)* &
                                                btp%graduvb_face_ave(ivar,1,iquad,iface) + &
                                                q_send(index+ivar,iquad,kk)

                  flux_uv_visc_face(ivar,2) = q_recv(index+5,iquad,kk)* &
                                                btp%graduvb_face_ave(ivar,2,iquad,iface) + &
                                                q_recv(index+ivar,iquad,kk)

               end do

               nxl = mf%normal_vector(1,iquad,1,iface)
               nyl = mf%normal_vector(2,iquad,1,iface)

               qul(1) = flux_uv_visc_face(1,1)
               qul(2) = flux_uv_visc_face(2,1)
               qvl(1) = flux_uv_visc_face(3,1)
               qvl(2) = flux_uv_visc_face(4,1)

               qur(1) = flux_uv_visc_face(1,2)
               qur(2) = flux_uv_visc_face(2,2)
               qvr(1) = flux_uv_visc_face(3,2)
               qvr(2) = flux_uv_visc_face(4,2)

               ! The Flip-Flop flux of Cockburn & Shu
               ! NOTE: beta=0.5 is the central flux
               qu_mean(1) = alpha*qul(1) + beta*qur(1)
               qu_mean(2) = alpha*qul(2) + beta*qur(2)
               qv_mean(1) = alpha*qvl(1) + beta*qvr(1)
               qv_mean(2) = alpha*qvl(2) + beta*qvr(2)

               wq=mf%jac_face(iquad,1,iface)

               flux_qu = (qu_mean(1) - iflux*qul(1))*nxl + (qu_mean(2) - iflux*qul(2))*nyl
               flux_qv = (qv_mean(1) - iflux*qvl(1))*nxl + (qv_mean(2) - iflux*qvl(2))*nyl

               !  Do Gauss-Lobatto Integration
               do i=1,b%ngl

                  hi = b%psi(i,iquad)

                  il=mf%imapl(1,i,1,iface)
                  jl=mf%imapl(2,i,1,iface)
                  kl=mf%imapl(3,i,1,iface)
                  ip=G%intma(il,jl,kl,el)

                  !Update Flux
                  rhs(1,ip,k) = rhs(1,ip,k) + wq*hi*flux_qu
                  rhs(2,ip,k) = rhs(2,ip,k) + wq*hi*flux_qv

               end do !i
            end do !iquad
         enddo !klayers

         kk=kk+1
         jj=jj+1
      end do

   end do !iface

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
   real, intent(in)    :: q_send(nvarb, nq, G%nboun, nlayers)
   real, intent(in)    :: q_recv(nvarb, nq, G%nboun, nlayers)

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
