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

subroutine create_nbhs_face_quad(q_face,q_send,q_recv,nvarb,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
 
   use mod_initial, only: nvar
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, &
                        gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                     num_nbh, num_send_recv
 
   use mod_ref, only: nmessage
 
   implicit none
 
   !global arrays
   real, intent(inout) :: q_face(nvarb,2,nq,nface)
   real,intent(in):: q_send(nvarb,nq,nboun)
   real,intent(in):: q_recv(nvarb,nq,nboun)
   integer, intent(in) :: nvarb
 
   !local variables
 
   real :: wq, a_constant, iflux
   integer ::iface,k
   integer :: iel, ier, ilocl, ilocr
   integer :: ifaceb, nq_i, nq_j, plane_ij
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
   integer :: ftype, pface, subface, imulti
   integer :: multirate
 
   !Constants
 
   jj=1
   kk=1
   imm = 0

   ! do ifaceb =1,nboun
   !    iface = face_send(ifaceb)  ! Get Local Face
 
   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)
 
         ftype = face_type(iface)
 
         do im = 1,imulti
 
            ! if (ftype==2) then
               ilocl=face(5,iface)
               iel=face(7,iface)
 
               ! if (multirate==1 .and. lev_list(iel)==0) then
               !    kk=kk+1
               !    cycle
               ! end if
            ! end if
 
            !-------------------------------------
            !Store Left Side Variables
            !-------------------------------------
            call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)
            
            do k=1,nvarb
               !Left Element
               q_face(k,1,:,iface)=q_send(k,:,kk)

               !Right Element
               q_face(k,2,:,iface)=q_recv(k,:,kk)
            end do
 
            kk=kk+1

            ! end if
         end do
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_quad

 subroutine create_nbhs_face_df(rhs,q_send,q_recv,nvarb,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq, psiq
   use mod_face, only: normal_vector_q, jac_faceq, imapl, imapr, face_send
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
   use mod_initial, only: nvar, alpha_mlswe
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv
 
   use mod_ref, only: nbtp_var
   use mod_input, only: nlayers, dry_cutoff
   use mod_constants, only: gravity
   use mod_variables, only: H_face_ave,ope_face_ave,btp_mass_flux_face_ave, &
                                Qu_face_ave, Qv_face_ave, one_plus_eta_edge_2_ave, &
                                uvb_face_ave, ope2_face_ave

   implicit none
 
   !global arrays
   real, intent(inout) :: rhs(3,npoin)
   real,intent(in):: q_send(nbtp_var,ngl,nboun)
   real,intent(in):: q_recv(nbtp_var,ngl,nboun)
   integer, intent(in) :: nvarb
 
   !local variables
 
   integer :: jj, imm, inbh, ib, kk, ivar, ii
   integer :: imulti
   integer :: multirate
   integer :: iface, iquad, el, er, il, jl, ir, jr, I, kl, kr, n
   real :: wq, hi, nxl, nyl, nxr, nyr, un
   real :: ul, ur, vl, vr, pbl, pbr, clam, one_eta
   real :: pU_L, pU_R, pbpert_edge
   real :: qbl(4), qbr(4), flux(3,nq)
   real :: c_minus, c_plus
   real :: ppl, ppr, upl, upr, vpl, vpr
   real, dimension(nlayers+1) :: pprime_l, pprime_r
   integer :: itype, k
   real :: H_bcl_ql, H_bcl_qr, H_bcl_q, flux_pb, flux_u, flux_v, fxl, fxr, flux_edge_x, flux_edge_y
   real, dimension(2) :: Qu_ql, Qu_qr, Qv_ql, Qv_qr

   jj=1
   kk=1
   imm = 0

   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)

         el=face(7,iface)
         er=face(8,iface)
 
         !-------------------------------------
         !Store Left Side Variables
         !-------------------------------------

         do iquad = 1, nq

            nxl = normal_vector_q(1,iquad,1,iface)
            nyl = normal_vector_q(2,iquad,1,iface)

            nxr = -nxl; nyr = -nyl

            qbl = 0.0; qbr = 0.0 ; pbl = 0.0; pbr = 0.0
            do n = 1,ngl
               hi = psiq(n,iquad)
               !Left Element
               qbl(1:4) = qbl(1:4) + hi*q_send(1:4,n,kk)
               pbl = pbl + hi*q_send(5,n,kk)
               !Right Element
               qbr(1:4) = qbr(1:4) + hi*q_recv(1:4,n,kk)
               pbr = pbr + hi*q_recv(5,n,kk)
            end do

            pU_L = nxl * qbl(3) + nyl * qbl(4)
            pU_R = nxr * qbr(3) + nyr * qbr(4)

            c_minus = sqrt(alpha_mlswe(nlayers) * pbr)
            c_plus  = sqrt(alpha_mlswe(nlayers) * pbl)
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
            do k = 1, nlayers
               ppl = 0.0; ppr = 0.0 ; upl = 0.0; upr = 0.0 ; vpl = 0.0; vpr = 0.0

               ii = 5 + (k-1)*3
               do n = 1, ngl
                  il = imapl(1,n,1,iface)
                  jl = imapl(2,n,1,iface)
                  kl = imapl(3,n,1,iface)
                  I = intma(il,jl,kl,el)

                  hi = psiq(n,iquad)
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

               H_bcl_ql  = H_bcl_ql + 0.5*alpha_mlswe(k)*(pprime_l(k+1)**2 - pprime_l(k)**2)
               H_bcl_qr  = H_bcl_qr + 0.5*alpha_mlswe(k)*(pprime_r(k+1)**2 - pprime_r(k)**2)
            end do

            ! Compute mass fluxes at each element face.

            flux_edge_x = 0.5*(qbl(3) + qbr(3)) + (0.5*clam)*(nxl * qbl(2) + nxr * qbr(2))
            flux_edge_y = 0.5*(qbl(4) + qbr(4)) + (0.5*clam)*(nyl * qbl(2) + nyr * qbr(2))
            ! flux(1,iquad) = nxl*flux_edge_x(iquad) + nyl*flux_edge_y(iquad)

            fxl = nxl*qbl(3) + nyl*qbl(4)
            fxr = nxr*qbr(3) + nyr*qbr(4)
            flux(1,iquad) = 0.5*(fxl - fxr) - 0.5*clam*(qbr(2) - qbl(2))

            ! Compute pressure forcing H_face at each element face.
            H_bcl_q = (one_eta**2) * (0.5 * (H_bcl_ql + H_bcl_qr))

            ! Accumulate sums for time averaging

            btp_mass_flux_face_ave(1,iquad,iface) = btp_mass_flux_face_ave(1,iquad,iface) + flux_edge_x
            btp_mass_flux_face_ave(2,iquad,iface) = btp_mass_flux_face_ave(2,iquad,iface) + flux_edge_y

            ! btp_mass_flux_face_ave(1,iquad,iface) = btp_mass_flux_face_ave(1,iquad,iface) + flux(1,iquad)
            ! btp_mass_flux_face_ave(2,iquad,iface) = btp_mass_flux_face_ave(2,iquad,iface) + flux(1,iquad)

            H_face_ave(iquad,iface) = H_face_ave(iquad,iface) + H_bcl_q
            Qu_face_ave(1,iquad,iface) = Qu_face_ave(1,iquad,iface) + 0.5*(Qu_ql(1) + Qu_qr(1))
            Qu_face_ave(2,iquad,iface) = Qu_face_ave(2,iquad,iface) + 0.5*(Qu_ql(2) + Qu_qr(2))
            Qv_face_ave(1,iquad,iface) = Qv_face_ave(1,iquad,iface) + 0.5*(Qv_ql(1) + Qv_qr(1))
            Qv_face_ave(2,iquad,iface) = Qv_face_ave(2,iquad,iface) + 0.5*(Qv_ql(2) + Qv_qr(2))
            ope_face_ave(1,iquad,iface) = ope_face_ave(1,iquad,iface) + (1.0 + (qbl(2)/pbl))
            ope_face_ave(2,iquad,iface) = ope_face_ave(2,iquad,iface) + (1.0 + (qbr(2)/pbr))
            ope2_face_ave(1,iquad,iface) = ope2_face_ave(1,iquad,iface) + (1.0 + (qbl(2)/pbl))**2
            ope2_face_ave(2,iquad,iface) = ope2_face_ave(2,iquad,iface) + (1.0 + (qbr(2)/pbr))**2
            one_plus_eta_edge_2_ave(iquad,iface) = one_plus_eta_edge_2_ave(iquad,iface) &
                                                + one_eta**2
            uvb_face_ave(1,1,iquad,iface) = uvb_face_ave(1,1,iquad,iface) + ul
            uvb_face_ave(1,2,iquad,iface) = uvb_face_ave(1,2,iquad,iface) + ur
            uvb_face_ave(2,1,iquad,iface) = uvb_face_ave(2,1,iquad,iface) + vl
            uvb_face_ave(2,2,iquad,iface) = uvb_face_ave(2,2,iquad,iface) + vr

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

         do iquad = 1, nq

            wq = jac_faceq(iquad,1,iface)

            flux_pb = flux(1,iquad)
            flux_u = flux(2,iquad)
            flux_v = flux(3,iquad)

            do n = 1, ngl

               hi = psiq(n,iquad)
               il = imapl(1,n,1,iface)
               jl = imapl(2,n,1,iface)
               kl = imapl(3,n,1,iface)
               I = intma(il,jl,kl,el)

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

 subroutine create_nbhs_face_df_lap(rhs,q_send,q_recv,nvarb,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq, psi
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, face_type,nface
 
   use mod_metrics, only: jac, massinv
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv

   use mod_variables, only: btp_graduv_dpp_face, graduvb_face_ave
 
   implicit none
 
   !global arrays
   real, intent(inout) :: rhs(2,npoin)
   real,intent(in):: q_send(10,ngl,nboun)
   real,intent(in):: q_recv(10,ngl,nboun)
   integer, intent(in) :: nvarb
 
   !local variables
 
   real :: wq, a_constant, iflux
   integer ::iface,k
   integer :: iel, ier, ilocl, ilocr
   integer :: ifaceb, nq_i, nq_j, plane_ij
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
   integer :: ftype, pface, subface, imulti
   integer :: multirate
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

   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)

         el=face(7,iface)
         er=face(8,iface)
         
         do iquad = 1,ngl

            do ivar = 6, 10
               grad_uvb_pb_l(ivar-5) =  q_send(ivar,iquad,kk)
               grad_uvb_pb_r(ivar-5) =  q_recv(ivar,iquad,kk)
            end do

            do ivar = 1,4
               flux_uv_visc_face(ivar,1) = grad_uvb_pb_l(5)* q_send(ivar,iquad,kk) + &
                                          grad_uvb_pb_l(ivar)

               flux_uv_visc_face(ivar,2) = grad_uvb_pb_r(5)* q_recv(ivar,iquad,kk) + &
                                          grad_uvb_pb_r(ivar)

               graduvb_face_ave(ivar,1,iquad,iface) = graduvb_face_ave(ivar,1,iquad,iface) + q_send(ivar,iquad,kk)
               graduvb_face_ave(ivar,2,iquad,iface) = graduvb_face_ave(ivar,2,iquad,iface) + q_recv(ivar,iquad,kk)

            end do

            nxl = normal_vector(1,iquad,1,iface)
            nyl = normal_vector(2,iquad,1,iface)

            qul(:) = flux_uv_visc_face(1:2,1)
            qvl(:) = flux_uv_visc_face(3:4,1)
            qur(:) = flux_uv_visc_face(1:2,2)
            qvr(:) = flux_uv_visc_face(3:4,2)

            ! The Flip-Flop flux of Cockburn & Shu 
            ! NOTE: beta=0.5 is the central flux
            qu_mean(:) = alpha*qul(:) + beta*qur(:)
            qv_mean(:) = alpha*qvl(:) + beta*qvr(:)

            wq=jac_face(iquad,1,iface)

            flux_qu = (qu_mean(1) - qul(1)*nxl) + (qu_mean(2) - qul(2)*nyl)
            flux_qv = (qv_mean(1) - qvl(1)*nxl) + (qv_mean(2) - qvl(2)*nyl)

            !  Do Gauss-Lobatto Integration
            do i=1,ngl

               hi = psi(i,iquad)

               il=imapl(1,i,1,iface)
               jl=imapl(2,i,1,iface)
               kl=imapl(3,i,1,iface)
               ip=intma(il,jl,kl,el)
               
               !Update Flux
               rhs(1,ip) = rhs(1,ip) + wq*hi*flux_qu
               rhs(2,ip) = rhs(2,ip) + wq*hi*flux_qv

            end do !i
         end do !iquad

            kk=kk+1
         ! end do
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_df_lap

 subroutine create_nbhs_face_bcl(rhs,q_send,q_recv,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq, psi, psiq
   use mod_face, only: normal_vector_q, jac_faceq, imapl, imapr, face_send
   use mod_grid, only: nelem, npoin, intma, face, nboun, face_type,nface
   use mod_metrics, only: jac, massinv
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv
   use mod_input, only: nlayers, dry_cutoff
   use mod_constants, only : gravity
   use mod_initial, only : alpha_mlswe, zbot_face
   use mod_variables, only: ope_face_ave, H_face_ave, one_plus_eta_edge_2_ave, &
                                uvb_face_ave, Quv_face_ave, Qu_face_ave, Qv_face_ave, ope2_face_ave
   use mod_variables, only: sum_layer_mass_flux_face, btp_mass_flux_face_ave
 
   implicit none
 
   !global arrays
   real, intent(inout) :: rhs(3,npoin,nlayers)
   real,intent(in):: q_send(3*nlayers,ngl,nboun)
   real,intent(in):: q_recv(3*nlayers,ngl,nboun)
 
   !local variables
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar, imulti
   integer :: multirate
   
   real, dimension(nlayers) :: alpha_over_g, g_over_alpha
   real, dimension(2,nlayers+1) :: p_face, z_face
   real, dimension(nlayers+1) :: p_edge_plus, p_edge_minus, p2l, p2r, z_edge_plus, z_edge_minus
   real, dimension(3,nlayers) :: ql, qr
   real, dimension(3,nq) :: qbl, qbr
   integer :: iface, ilr, k, iquad, ktemp, I
   real :: z_intersect_top,z_intersect_bot, dz_intersect, H_r_plus, H_r_minus, acceleration
   real :: p_intersect_bot, p_intersect_top, one_plus_eta_edge
   real :: H_corr,p_inc, weight, H_corr1,p_inc1, H_corr2,p_inc2, temp, ope_l, ope_r
   integer :: Iq, el, er
   real ::  ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv
   real, dimension(nlayers) :: udpl, udpr, vdpl, vdpr
   real :: uu_dp_flux_deficit(2), vv_dp_flux_deficit(2)
   real, parameter :: eps1 = 1.0e-20 !  Parameter used to prevent division by zero.
   integer :: il, jl, ir, jr, kl, kr, jquad, n, m, index
   real :: wq, hi, hlx_k, hly_k, hrx_k, hry_k, flux_x, flux_y, hx_k, hy_k, flux
   real, dimension(2,nq,nlayers) :: H_face, udp_flux, vdp_flux, dp_flux
   real :: dp_lr(2,nlayers), dp_deficit(2)

   jj=1
   kk=1
   imm = 0

   do k=1,nlayers
      alpha_over_g(k) = alpha_mlswe(k)/gravity
      g_over_alpha(k) = gravity/alpha_mlswe(k)
   enddo

   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)

         el=face(7,iface)
         er=face(8,iface)

         qbl(1,:) = ope_face_ave(1,:,iface)
         qbl(2,:) = uvb_face_ave(1,1,:,iface)
         qbl(3,:) = uvb_face_ave(2,1,:,iface)
         qbr(1,:) = ope_face_ave(2,:,iface)
         qbr(2,:) = uvb_face_ave(1,2,:,iface)
         qbr(3,:) = uvb_face_ave(2,2,:,iface)

         do iquad = 1,nq

            nxl = normal_vector_q(1,iquad,1,iface)
            nyl = normal_vector_q(2,iquad,1,iface)

            ql = 0.0; qr = 0.0

            do k = 1,nlayers

               index = 3*(k-1)

               do n = 1, ngl
                  hi = psiq(n,iquad)
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

            ! sum_layer_mass_flux_face(1,iquad,iface) = sum_layer_mass_flux_face(1,iquad,iface)  + &
            !                                           sum(udp_flux(1,iquad,:))
            ! sum_layer_mass_flux_face(2,iquad,iface) = sum_layer_mass_flux_face(2,iquad,iface)  + &
            !                                           sum(vdp_flux(1,iquad,:))

            uu_dp_flux_deficit(1) = Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,iquad,:))
            uu_dp_flux_deficit(2) = Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,iquad,:))
            vv_dp_flux_deficit(1) = Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,iquad,:))
            vv_dp_flux_deficit(2) = Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,iquad,:))

            dp_deficit(1) = btp_mass_flux_face_ave(1,iquad,iface) - sum(dp_flux(1,iquad,:))
            dp_deficit(2) = btp_mass_flux_face_ave(2,iquad,iface) - sum(dp_flux(2,iquad,:))

            do k = 1,nlayers

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
               weight = abs(vdpl(k)) / (sum(abs(vdpl(:))+eps1))
               if(uu_dp_flux_deficit(2)*nyl < 0.0) &
                  weight = abs(vdpr(k)) / (sum(abs(vdpr(:))+eps1))
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
            ope_l = sqrt(ope2_face_ave(1,iquad,iface))
            ope_r = sqrt(ope2_face_ave(2,iquad,iface))
            p_face(1,1) = 0.0
            p_face(2,1) = 0.0
            do k=1,nlayers
               p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,k)
               p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,k)
            end do

            one_plus_eta_edge = sqrt(one_plus_eta_edge_2_ave(iquad,iface))
            z_face(1,nlayers+1) = zbot_face(1,iquad,iface)
            z_face(2,nlayers+1) = zbot_face(2,iquad,iface)
            z_edge_plus(nlayers+1) = zbot_face(1,iquad,iface)
            z_edge_minus(nlayers+1) = zbot_face(2,iquad,iface)
            do k=nlayers,1,-1
               z_face(1,k) = z_face(1,k+1) + alpha_over_g(k) * (ope_l * ql(1,k))
               z_face(2,k) = z_face(2,k+1) + alpha_over_g(k) * (ope_r * qr(1,k))
               z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k) * &
                                                   (one_plus_eta_edge * ql(1,k))
               z_edge_minus(k) = z_edge_minus(k+1) + alpha_over_g(k) * &
                                                   (one_plus_eta_edge * qr(1,k))
            end do

            p_edge_plus(2) = one_plus_eta_edge * ql(1,1)
            p_edge_minus(2) = one_plus_eta_edge * qr(1,1)
            do k = 2,nlayers
               p_edge_plus(k+1) = p_edge_plus(k) + one_plus_eta_edge * ql(1,k)
               p_edge_minus(k+1) = p_edge_minus(k) + one_plus_eta_edge * qr(1,k)
            end do

            do k = 1, nlayers

               ! Computation from + side for layer k
               H_r_plus = 0.5*alpha_mlswe(k)*(p_edge_plus(k+1)**2 - p_edge_plus(k)**2)

               ! Computation from - side for layer k
               H_r_minus = 0.0
               do ktemp = 1, nlayers

                  z_intersect_top = min(z_edge_minus(ktemp), z_edge_plus(k))
                  z_intersect_bot = max(z_edge_minus(ktemp+1), z_edge_plus(k+1))
                  dz_intersect = z_intersect_top - z_intersect_bot

                  if (dz_intersect > 0.0) then
                     p_intersect_bot = p_edge_minus(ktemp+1) &
                           - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_minus(ktemp+1))
                     p_intersect_top = p_edge_minus(ktemp+1) &
                           - g_over_alpha(ktemp)*(z_intersect_top - z_edge_minus(ktemp+1))
                     H_r_minus = H_r_minus + &
                           0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)

                  end if
               end do
               H_face(1,iquad,k) = 0.5*(H_r_plus + H_r_minus) !computation of H_r for the left side
               ! Computation from - side for layer k
               H_r_minus = 0.5*alpha_mlswe(k)*(p_edge_minus(k+1)**2 - p_edge_minus(k)**2)

               ! Computation from + side for layer k
               H_r_plus = 0.0
               do ktemp = 1, nlayers

                  z_intersect_top = min(z_edge_plus(ktemp), z_edge_minus(k))
                  z_intersect_bot = max(z_edge_plus(ktemp+1), z_edge_minus(k+1))
                  dz_intersect = z_intersect_top - z_intersect_bot

                  if (dz_intersect > 0.0) then
                     p_intersect_bot = p_edge_plus(ktemp+1) &
                        - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_plus(ktemp+1))
                     p_intersect_top = p_edge_plus(ktemp+1) &
                        - g_over_alpha(ktemp)*(z_intersect_top - z_edge_plus(ktemp+1))
                     H_r_plus = H_r_plus &
                        + 0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
                  end if
               end do
               H_face(2,iquad,k) = 0.5*(H_r_plus + H_r_minus) ! H_r for the right side
            end do !k

            do k = 1, nlayers-1          ! interface at the bottom of layer k
               ! Corrections at the left side of a face.
               p_inc1 = g_over_alpha(k)*(z_face(1,k+1) - z_edge_plus(k+1))
               H_corr1 = 0.5 * alpha_mlswe(k) * ((p_face(1,k+1) &
                           + p_inc1)**2 - p_face(1,k+1)**2)
               H_face(1,iquad,k) = H_face(1,iquad,k) - H_corr1
               H_face(1,iquad,k+1) = H_face(1,iquad,k+1) + H_corr1

               ! Corrections at the right side of a face.
               p_inc2 = g_over_alpha(k)*(z_face(2,k+1) - z_edge_minus(k+1))
               H_corr2 = 0.5 * alpha_mlswe(k) * ((p_face(2,k+1) &
                           + p_inc2)**2 - p_face(2,k+1)**2)
               H_face(2,iquad,k) = H_face(2,iquad,k) - H_corr2
               H_face(2,iquad,k+1) = H_face(2,iquad,k+1) + H_corr2

            end do

            ! Left side of face
            weight = 1.0
            acceleration = sum(H_face(1,iquad,:))
            if(acceleration > 0.0) then
               weight = H_face_ave(iquad,iface) / acceleration
            end if
            H_face(1,iquad,:) = H_face(1,iquad,:) * weight

            ! Right side of face
            weight = 1.0
            acceleration = sum(H_face(2,iquad,:))
            if(acceleration > 0.0) then
               weight = H_face_ave(iquad,iface) / acceleration
            end if
            H_face(2,iquad,:) = H_face(2,iquad,:) * weight

         end do ! iquad
         
         ! Do Gauss-Lobatto Integration
         do k = 1,nlayers

            do iquad = 1, nq

               wq = jac_faceq(iquad,1,iface)
               nxl = normal_vector_q(1,iquad,1,iface)
               nyl = normal_vector_q(2,iquad,1,iface)

               hlx_k = nxl*H_face(1,iquad,k)
               hrx_k = nxl*H_face(2,iquad,k)
               hly_k = nyl*H_face(1,iquad,k)
               hry_k = nyl*H_face(2,iquad,k)
               flux = nxl*dp_flux(1,iquad,k) + nyl*dp_flux(2,iquad,k)
               flux_x = nxl*udp_flux(1,iquad,k) + nyl*udp_flux(2,iquad,k)
               flux_y = nxl*vdp_flux(1,iquad,k) + nyl*vdp_flux(2,iquad,k)

               do n = 1, ngl

                  hi = psiq(n,iquad)
                  il = imapl(1,n,1,iface)
                  jl = imapl(2,n,1,iface)
                  kl = imapl(3,n,1,iface)
                  I = intma(il,jl,kl,el)

                  rhs(1,I,k) = rhs(1,I,k) - wq*hi*flux
                  rhs(2,I,k) = rhs(2,I,k) - wq*hi*(hlx_k + flux_x)
                  rhs(3,I,k) = rhs(3,I,k) - wq*hi*(hly_k + flux_y)
               end do
            end do
         end do

         kk=kk+1
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_bcl

  subroutine create_nbhs_face_bcl_continuity(rhs,q_send,q_recv,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq, psi, psiq
   use mod_face, only: normal_vector_q, jac_faceq, imapl, imapr, face_send
   use mod_grid, only: nelem, npoin, intma, face, nboun, face_type,nface
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv
   use mod_input, only: nlayers, dry_cutoff
   use mod_variables, only: ope_face_ave, uvb_face_ave, sum_layer_mass_flux_face, btp_mass_flux_face_ave
   use mod_constants, only : gravity
   use mod_initial, only : alpha_mlswe
 
   implicit none
 
   !global arrays
   real, intent(inout) :: rhs(npoin,nlayers)
   real,intent(in):: q_send(3*nlayers,ngl,nboun)
   real,intent(in):: q_recv(3*nlayers,ngl,nboun)
 
   !local variables
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar, imulti
   integer :: multirate
   
   integer :: k, iface, iquad, el, er, il, jl, ir, jr, I
   integer :: kl, kr, jquad, n, m, index
   real :: wq, nxl, nyl, hi
   real :: dpl, dpr, uu, vv, flux, ul, ur, vl, vr, weight, dp_lr(2,nlayers), dp_deficit(2)
   real, dimension(nq,nlayers) :: flux_edge_u, flux_edge_v
   real, dimension(3) :: ql, qr
   real, dimension(3,nq) :: qbl, qbr
   real, parameter :: eps = 1.0e-10

   jj=1
   kk=1
   imm = 0

   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)

         el=face(7,iface)
         er=face(8,iface)

         qbl(1,:) = ope_face_ave(1,:,iface)
         qbl(2,:) = uvb_face_ave(1,1,:,iface)
         qbl(3,:) = uvb_face_ave(2,1,:,iface)
         qbr(1,:) = ope_face_ave(2,:,iface)
         qbr(2,:) = uvb_face_ave(1,2,:,iface)
         qbr(3,:) = uvb_face_ave(2,2,:,iface)

         do iquad = 1,nq

            nxl = normal_vector_q(1,iquad,1,iface)
            nyl = normal_vector_q(2,iquad,1,iface)

            do k = 1,nlayers

               index = 3*(k-1)

               ql = 0.0; qr = 0.0
               do n = 1, ngl
                  hi = psiq(n,iquad)
                  do ivar = 1,3
                     !Left Element
                     ql(ivar) = ql(ivar) + hi*q_send(index+ivar,n,kk)
                     !Right Element
                     qr(ivar) = qr(ivar) + hi*q_recv(index+ivar,n,kk)
                  end do
               enddo

               dpl = qbl(1,iquad) * ql(1)
               dpr = qbr(1,iquad) * qr(1)

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

            dp_deficit(1) = btp_mass_flux_face_ave(1,iquad,iface) - sum(flux_edge_u(iquad,:))
            dp_deficit(2) = btp_mass_flux_face_ave(2,iquad,iface) - sum(flux_edge_v(iquad,:))

            do k = 1, nlayers

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
         do k = 1,nlayers

            do iquad = 1, nq

               wq = jac_faceq(iquad,1,iface)
               nxl = normal_vector_q(1,iquad,1,iface)
               nyl = normal_vector_q(2,iquad,1,iface)

               flux = nxl*flux_edge_u(iquad,k) + nyl*flux_edge_v(iquad,k)

               do n = 1, ngl

                  hi = psiq(n,iquad)
                  il = imapl(1,n,1,iface)
                  jl = imapl(2,n,1,iface)
                  kl = imapl(3,n,1,iface)
                  I = intma(il,jl,kl,el)

                  rhs(I,k) = rhs(I,k) - wq*hi*flux
               end do
            end do
         end do

         kk=kk+1
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_bcl_continuity

 subroutine create_nbhs_face_bcl_momentum(rhs,q_send,q_recv,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq, psi, psiq
   use mod_face, only: normal_vector_q, jac_faceq, imapl, imapr, face_send
   use mod_grid, only: nelem, npoin, intma, face, nboun, face_type,nface
   use mod_metrics, only: jac, massinv
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv
   use mod_input, only: nlayers, dry_cutoff
   use mod_constants, only : gravity
   use mod_initial, only : alpha_mlswe, zbot_face
   use mod_variables, only: ope_face_ave, H_face_ave, one_plus_eta_edge_2_ave, &
                                uvb_face_ave, Quv_face_ave, Qu_face_ave, Qv_face_ave, ope2_face_ave
 
   implicit none
 
   !global arrays
   real, intent(inout) :: rhs(2,npoin,nlayers)
   real,intent(in):: q_send(3*nlayers,ngl,nboun)
   real,intent(in):: q_recv(3*nlayers,ngl,nboun)
 
   !local variables
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar, imulti
   integer :: multirate
   
   real, dimension(nlayers) :: alpha_over_g, g_over_alpha
   real, dimension(2,nlayers+1) :: p_face, z_face
   real, dimension(nlayers+1) :: p_edge_plus, p_edge_minus, p2l, p2r, z_edge_plus, z_edge_minus
   real, dimension(3,nq,nlayers) :: ql, qr
   real, dimension(3,nq) :: qbl, qbr
   integer :: iface, ilr, k, iquad, ktemp, I
   real :: z_intersect_top,z_intersect_bot, dz_intersect, H_r_plus, H_r_minus, acceleration
   real :: p_intersect_bot, p_intersect_top, one_plus_eta_edge
   real :: H_corr,p_inc, weight, H_corr1,p_inc1, H_corr2,p_inc2, temp, ope_l, ope_r
   integer :: Iq, el, er
   real ::  ul, ur, vl, vr, dpl, dpr, nxl, nyl, uu, vv
   real, dimension(nq,nlayers) :: udpl, udpr, vdpl, vdpr
   real :: one_over_sum_l, one_over_sum_r, uu_dp_flux_deficit, uv_dp_flux_deficit
   real :: vu_dp_flux_deficit, vv_dp_flux_deficit
   real, parameter :: eps1 = 1.0e-20 !  Parameter used to prevent division by zero.
   integer :: il, jl, ir, jr, kl, kr, jquad, n, m, index
   real :: wq, hi, hlx_k, hly_k, hrx_k, hry_k, flux_x, flux_y, hx_k, hy_k, flux
   real, dimension(2,nq,nlayers) :: H_face, udp_flux, vdp_flux

   jj=1
   kk=1
   imm = 0

   do k=1,nlayers
      alpha_over_g(k) = alpha_mlswe(k)/gravity
      g_over_alpha(k) = gravity/alpha_mlswe(k)
   enddo

   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)

         el=face(7,iface)
         er=face(8,iface)

         qbl(1,:) = ope_face_ave(1,:,iface)
         qbl(2,:) = uvb_face_ave(1,1,:,iface)
         qbl(3,:) = uvb_face_ave(2,1,:,iface)
         qbr(1,:) = ope_face_ave(2,:,iface)
         qbr(2,:) = uvb_face_ave(1,2,:,iface)
         qbr(3,:) = uvb_face_ave(2,2,:,iface)

         ql = 0.0; qr = 0.0
         do iquad = 1,nq

            nxl = normal_vector_q(1,iquad,1,iface)
            nyl = normal_vector_q(2,iquad,1,iface)

            do k = 1,nlayers

               index = 3*(k-1)

               do n = 1, ngl
                  hi = psiq(n,iquad)
                  do ivar = 1,3
                     !Left Element
                     ql(ivar,iquad,k) = ql(ivar,iquad,k) + hi*q_send(index+ivar,n,kk)
                     !Right Element
                     qr(ivar,iquad,k) = qr(ivar,iquad,k) + hi*q_recv(index+ivar,n,kk)
                  end do
               enddo

               if (ql(1,iquad,k) < g_over_alpha(k)*dry_cutoff) then
                  ql(1,iquad,k) = g_over_alpha(k)*dry_cutoff
                  ql(2:3,iquad,k) = 0.0
               end if
               if (qr(1,iquad,k) < g_over_alpha(k)*dry_cutoff) then
                  qr(1,iquad,k) = g_over_alpha(k)*dry_cutoff
                  qr(2:3,iquad,k) = 0.0
               end if

               ! Left side of the edge
               dpl = qbl(1,iquad) * ql(1,iquad,k)
               dpr = qbr(1,iquad) * qr(1,iquad,k)
               ul = ql(2,iquad,k)+qbl(2,iquad)
               ur = qr(2,iquad,k)+qbr(2,iquad)
               vl = ql(3,iquad,k)+qbl(3,iquad)
               vr = qr(3,iquad,k)+qbr(3,iquad)

               uu = 0.5*(ul+ur)
               vv = 0.5*(vl+vr)
               udpl(iquad,k) = ul*dpl
               udpr(iquad,k) = ur*dpr
               vdpl(iquad,k) = vl*dpl
               vdpr(iquad,k) = vr*dpr

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

            uu_dp_flux_deficit = Qu_face_ave(1,iquad,iface) - sum(udp_flux(1,iquad,:))
            uv_dp_flux_deficit = Qu_face_ave(2,iquad,iface) - sum(udp_flux(2,iquad,:))
            vu_dp_flux_deficit = Qv_face_ave(1,iquad,iface) - sum(vdp_flux(1,iquad,:))
            vv_dp_flux_deficit = Qv_face_ave(2,iquad,iface) - sum(vdp_flux(2,iquad,:))

            ! Adjust the fluxes for the u-momentum equation
            one_over_sum_l = 1.0 / sum(abs(udpl(iquad,:))+eps1)
            one_over_sum_r = 1.0 / sum(abs(udpr(iquad,:))+eps1)
            !x-direction
            if(uu_dp_flux_deficit*nxl > 0.0) then
               do k = 1,nlayers
                  weight = abs(udpl(iquad,k)) * one_over_sum_l
                  udp_flux(1,iquad,k) = udp_flux(1,iquad,k) + weight * uu_dp_flux_deficit
               end do
            else
               do k = 1,nlayers
                  weight = abs(udpr(iquad,k)) * one_over_sum_r
                  udp_flux(1,iquad,k) = udp_flux(1,iquad,k) + weight * uu_dp_flux_deficit
               end do
            end if
            !y-direction
            if(uv_dp_flux_deficit*nyl > 0.0) then
               do k = 1,nlayers
                  weight = abs(udpl(iquad,k)) * one_over_sum_l
                  udp_flux(2,iquad,k) = udp_flux(2,iquad,k) + weight * uv_dp_flux_deficit
               end do
            else
               do k = 1,nlayers
                  weight = abs(udpr(iquad,k)) * one_over_sum_r
                  udp_flux(2,iquad,k) = udp_flux(2,iquad,k) + weight * uv_dp_flux_deficit
               end do
            end if
            ! Adjust the fluxes for the v-momentum equation
            one_over_sum_l = 1.0 / sum(abs(vdpl(iquad,:))+eps1)
            one_over_sum_r = 1.0 / sum(abs(vdpr(iquad,:))+eps1)
            !x-direction
            if(vu_dp_flux_deficit*nxl > 0.0) then
               do k = 1,nlayers
                  weight = abs(vdpl(iquad,k)) * one_over_sum_l
                  vdp_flux(2,iquad,k) = vdp_flux(2,iquad,k) + weight * vu_dp_flux_deficit
               end do
            else
               do k = 1,nlayers
                  weight = abs(vdpr(iquad,k)) * one_over_sum_r
                  vdp_flux(2,iquad,k) = vdp_flux(2,iquad,k) + weight * vu_dp_flux_deficit
               end do
            end if
            !y-direction
            if(vv_dp_flux_deficit*nyl > 0.0) then
               do k = 1,nlayers
                  weight = abs(vdpl(iquad,k)) * one_over_sum_l
                  vdp_flux(2,iquad,k) = vdp_flux(2,iquad,k) + weight * vv_dp_flux_deficit
               end do
            else
               do k = 1,nlayers
                  weight = abs(vdpr(iquad,k)) * one_over_sum_r
                  vdp_flux(2,iquad,k) = vdp_flux(2,iquad,k) + weight * vv_dp_flux_deficit
               end do
            end if

            z_face = 0.0 ; p_face = 0.0
            z_edge_plus = 0.0 ; z_edge_minus = 0.0
            p_edge_plus = 0.0; p_edge_minus = 0.0

            !Store Left Side Variables
            ope_l = sqrt(ope2_face_ave(1,iquad,iface))
            ope_r = sqrt(ope2_face_ave(2,iquad,iface))
            p_face(1,1) = 0.0
            p_face(2,1) = 0.0
            do k=1,nlayers
               p_face(1,k+1) = p_face(1,k) + ope_l * ql(1,iquad,k)
               p_face(2,k+1) = p_face(2,k) + ope_r * qr(1,iquad,k)
            end do

            one_plus_eta_edge = sqrt(one_plus_eta_edge_2_ave(iquad,iface))
            z_face(1,nlayers+1) = zbot_face(1,iquad,iface)
            z_face(2,nlayers+1) = zbot_face(2,iquad,iface)
            z_edge_plus(nlayers+1) = zbot_face(1,iquad,iface)
            z_edge_minus(nlayers+1) = zbot_face(2,iquad,iface)
            do k=nlayers,1,-1
               z_face(1,k) = z_face(1,k+1) + alpha_over_g(k) * (ope_l * ql(1,iquad,k))
               z_face(2,k) = z_face(2,k+1) + alpha_over_g(k) * (ope_r * qr(1,iquad,k))
               z_edge_plus(k) = z_edge_plus(k+1) + alpha_over_g(k) * &
                                                   (one_plus_eta_edge * ql(1,iquad,k))
               z_edge_minus(k) = z_edge_minus(k+1) + alpha_over_g(k) * &
                                                   (one_plus_eta_edge * qr(1,iquad,k))
            end do

            p_edge_plus(2) = one_plus_eta_edge * ql(1,iquad,1)
            p_edge_minus(2) = one_plus_eta_edge * qr(1,iquad,1)
            do k = 2,nlayers
               p_edge_plus(k+1) = p_edge_plus(k) + one_plus_eta_edge * ql(1,iquad,k)
               p_edge_minus(k+1) = p_edge_minus(k) + one_plus_eta_edge * qr(1,iquad,k)
            end do

            do k = 1, nlayers

               ! Computation from + side for layer k
               H_r_plus = 0.5*alpha_mlswe(k)*(p_edge_plus(k+1)**2 - p_edge_plus(k)**2)

               ! Computation from - side for layer k
               H_r_minus = 0.0
               do ktemp = 1, nlayers

                  z_intersect_top = min(z_edge_minus(ktemp), z_edge_plus(k))
                  z_intersect_bot = max(z_edge_minus(ktemp+1), z_edge_plus(k+1))
                  dz_intersect = z_intersect_top - z_intersect_bot

                  if (dz_intersect > 0.0) then
                     p_intersect_bot = p_edge_minus(ktemp+1) &
                           - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_minus(ktemp+1))
                     p_intersect_top = p_edge_minus(ktemp+1) &
                           - g_over_alpha(ktemp)*(z_intersect_top - z_edge_minus(ktemp+1))
                     H_r_minus = H_r_minus + &
                           0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)

                  end if
               end do
               H_face(1,iquad,k) = 0.5*(H_r_plus + H_r_minus) !computation of H_r for the left side
               ! Computation from - side for layer k
               H_r_minus = 0.5*alpha_mlswe(k)*(p_edge_minus(k+1)**2 - p_edge_minus(k)**2)

               ! Computation from + side for layer k
               H_r_plus = 0.0
               do ktemp = 1, nlayers

                  z_intersect_top = min(z_edge_plus(ktemp), z_edge_minus(k))
                  z_intersect_bot = max(z_edge_plus(ktemp+1), z_edge_minus(k+1))
                  dz_intersect = z_intersect_top - z_intersect_bot

                  if (dz_intersect > 0.0) then
                     p_intersect_bot = p_edge_plus(ktemp+1) &
                        - g_over_alpha(ktemp)*(z_intersect_bot - z_edge_plus(ktemp+1))
                     p_intersect_top = p_edge_plus(ktemp+1) &
                        - g_over_alpha(ktemp)*(z_intersect_top - z_edge_plus(ktemp+1))
                     H_r_plus = H_r_plus &
                        + 0.5*alpha_mlswe(ktemp)*(p_intersect_bot**2 - p_intersect_top**2)
                  end if
               end do
               H_face(2,iquad,k) = 0.5*(H_r_plus + H_r_minus) ! H_r for the right side
            end do !k

            do k = 1, nlayers-1          ! interface at the bottom of layer k
               ! Corrections at the left side of a face.
               p_inc1 = g_over_alpha(k)*(z_face(1,k+1) - z_edge_plus(k+1))
               H_corr1 = 0.5 * alpha_mlswe(k) * ((p_face(1,k+1) &
                           + p_inc1)**2 - p_face(1,k+1)**2)
               H_face(1,iquad,k) = H_face(1,iquad,k) - H_corr1
               H_face(1,iquad,k+1) = H_face(1,iquad,k+1) + H_corr1

               ! Corrections at the right side of a face.
               p_inc2 = g_over_alpha(k)*(z_face(2,k+1) - z_edge_minus(k+1))
               H_corr2 = 0.5 * alpha_mlswe(k) * ((p_face(2,k+1) &
                           + p_inc2)**2 - p_face(2,k+1)**2)
               H_face(2,iquad,k) = H_face(2,iquad,k) - H_corr2
               H_face(2,iquad,k+1) = H_face(2,iquad,k+1) + H_corr2

            end do

            ! Left side of face
            weight = 1.0
            acceleration = sum(H_face(1,iquad,:))
            if(acceleration > 0.0) then
               weight = H_face_ave(iquad,iface) / acceleration
            end if
            H_face(1,iquad,:) = H_face(1,iquad,:) * weight

            ! Right side of face
            weight = 1.0
            acceleration = sum(H_face(2,iquad,:))
            if(acceleration > 0.0) then
               weight = H_face_ave(iquad,iface) / acceleration
            end if
            H_face(2,iquad,:) = H_face(2,iquad,:) * weight

         end do ! iquad
         
         ! Do Gauss-Lobatto Integration
         do k = 1,nlayers

            do iquad = 1, nq

               wq = jac_faceq(iquad,1,iface)
               nxl = normal_vector_q(1,iquad,1,iface)
               nyl = normal_vector_q(2,iquad,1,iface)

               hlx_k = nxl*H_face(1,iquad,k)
               hrx_k = nxl*H_face(2,iquad,k)
               hly_k = nyl*H_face(1,iquad,k)
               hry_k = nyl*H_face(2,iquad,k)
               flux_x = nxl*udp_flux(1,iquad,k) + nyl*udp_flux(2,iquad,k)
               flux_y = nxl*vdp_flux(1,iquad,k) + nyl*vdp_flux(2,iquad,k)

               do n = 1, ngl

                  hi = psiq(n,iquad)
                  il = imapl(1,n,1,iface)
                  jl = imapl(2,n,1,iface)
                  kl = imapl(3,n,1,iface)
                  I = intma(il,jl,kl,el)

                  rhs(1,I,k) = rhs(1,I,k) - wq*hi*(hlx_k + flux_x)
                  rhs(2,I,k) = rhs(2,I,k) - wq*hi*(hly_k + flux_y)
               end do
            end do
         end do

         kk=kk+1
         jj=jj+1
      end do
 
   end do !iface

 end subroutine create_nbhs_face_bcl_momentum

  subroutine create_nbhs_face_df_lap_bcl(rhs,q_send,q_recv,nlayers,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq, psi
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, face_type,nface
 
   use mod_metrics, only: jac, massinv
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv

   use mod_variables, only: graduvb_face_ave
 
   implicit none
 
   !global arrays
   real, intent(inout) :: rhs(2,npoin,nlayers)
   real,intent(in):: q_send(5*nlayers,ngl,nboun)
   real,intent(in):: q_recv(5*nlayers,ngl,nboun)
   integer, intent(in) :: nlayers
 
   !local variables
 
   real :: wq, a_constant, iflux
   integer ::iface,k
   integer :: iel, ier, ilocl, ilocr
   integer :: ifaceb, nq_i, nq_j, plane_ij
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
   integer :: ftype, pface, subface, imulti
   integer :: multirate
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

   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)

         el=face(7,iface)
         er=face(8,iface)

         do k = 1, nlayers

            index = 5*(k-1)
         
            do iquad = 1,ngl

               do ivar = 1,4
                  flux_uv_visc_face(ivar,1) = q_send(index+5,iquad,kk)* &
                                                graduvb_face_ave(ivar,1,iquad,iface) + &
                                                q_send(index+ivar,iquad,kk)

                  flux_uv_visc_face(ivar,2) = q_recv(index+5,iquad,kk)* &
                                                graduvb_face_ave(ivar,2,iquad,iface) + &
                                                q_recv(index+ivar,iquad,kk)

               end do 

               nxl = normal_vector(1,iquad,1,iface)
               nyl = normal_vector(2,iquad,1,iface)

               qul(:) = flux_uv_visc_face(1:2,1)
               qvl(:) = flux_uv_visc_face(3:4,1)
               qur(:) = flux_uv_visc_face(1:2,2)
               qvr(:) = flux_uv_visc_face(3:4,2)

               ! The Flip-Flop flux of Cockburn & Shu 
               ! NOTE: beta=0.5 is the central flux
               qu_mean(:) = alpha*qul(:) + beta*qur(:)
               qv_mean(:) = alpha*qvl(:) + beta*qvr(:)

               wq=jac_face(iquad,1,iface)

               flux_qu = (qu_mean(1) - qul(1)*nxl) + (qu_mean(2) - qul(2)*nyl)
               flux_qv = (qv_mean(1) - qvl(1)*nxl) + (qv_mean(2) - qvl(2)*nyl)

               !  Do Gauss-Lobatto Integration
               do i=1,ngl

                  hi = psi(i,iquad)

                  il=imapl(1,i,1,iface)
                  jl=imapl(2,i,1,iface)
                  kl=imapl(3,i,1,iface)
                  ip=intma(il,jl,kl,el)
                  
                  !Update Flux
                  rhs(1,ip,k) = rhs(1,ip,k) + wq*hi*flux_qu
                  rhs(2,ip,k) = rhs(2,ip,k) + wq*hi*flux_qv

               end do !i
            end do !iquad
         enddo !klayers

            kk=kk+1
         ! end do
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_df_lap_bcl

 subroutine create_nbhs_face_consistency(rhs,q_send,q_recv,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq, psiq
   use mod_face, only: normal_vector_q, jac_faceq, imapl, imapr, face_send
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
   use mod_initial, only: nvar, alpha_mlswe
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv
 
   use mod_input, only: nlayers
   use mod_constants, only: gravity
   use mod_variables, only: btp_mass_flux_face_ave, sum_layer_mass_flux_face
 
   implicit none
 
   !global arrays
   real, intent(inout) :: rhs(npoin,nlayers)
   real,intent(in):: q_send(nlayers+1,ngl,nboun)
   real,intent(in):: q_recv(nlayers+1,ngl,nboun)
 
   !local variables
 
   real :: wq
   integer :: iface, k
 
   integer :: jj, imm, inbh, ib, kk, ivar, ii
   integer :: imulti
   integer :: multirate
   real :: nxl, nyl, nxr, nyr
   integer :: el, iquad, n, I, il, jl, kl, er, itype
   real, dimension(nlayers) :: qpl, qpr
   real, dimension(nq,nlayers) :: flux_edge_u, flux_edge_v
   real :: mass_deficit_u_l, mass_deficit_v_l
   real :: mass_deficit_u_r, mass_deficit_v_r
   real :: weight, flux, pbl, pbr, hi

   jj=1
   kk=1
   imm = 0

   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)

         el=face(7,iface)
         er=face(8,iface)
 
         !-------------------------------------
         !Store Left Side Variables
         !-------------------------------------

         do iquad = 1, nq

            nxl = normal_vector_q(1,iquad,1,iface)
            nyl = normal_vector_q(2,iquad,1,iface)

            qpl = 0.0 ; qpr = 0.0 ; pbl = 0.0 ; pbr = 0.0
            do n = 1, ngl
               il = imapl(1,n,1,iface)
               jl = imapl(2,n,1,iface)
               kl = imapl(3,n,1,iface)
               I = intma(il,jl,kl,el)

               hi = psiq(n,iquad)

               do k = 1,nlayers
                  qpl(k) = qpl(k) + hi*q_send(k,n,kk)
                  qpr(k) = qpr(k) + hi*q_recv(k,n,kk)
               end do
               pbl = pbl + hi*q_send(nlayers+1,n,kk)
               pbr = pbr + hi*q_recv(nlayers+1,n,kk)
            enddo

            do k = 1, nlayers

               weight = qpl(k) / pbl
               mass_deficit_u_l = weight * (btp_mass_flux_face_ave(1,iquad,iface) &
                                                        - sum_layer_mass_flux_face(1,iquad,iface))

               mass_deficit_v_l = weight * (btp_mass_flux_face_ave(2,iquad,iface) &
                                                        - sum_layer_mass_flux_face(2,iquad,iface))

               weight = qpr(k) / pbr
               mass_deficit_u_r = weight * (btp_mass_flux_face_ave(1,iquad,iface) &
                                                        - sum_layer_mass_flux_face(1,iquad,iface))
               mass_deficit_v_r = weight * (btp_mass_flux_face_ave(2,iquad,iface) &
                                                        - sum_layer_mass_flux_face(2,iquad,iface))

               if(mass_deficit_u_l*nxl > 0.0) then 

                  flux_edge_u(iquad,k) = mass_deficit_u_l
               else 
                  flux_edge_u(iquad,k) = mass_deficit_u_r
               end if 

               if(mass_deficit_v_l*nyl > 0.0) then 
                  flux_edge_v(iquad,k) = mass_deficit_v_l
               else 
                  flux_edge_v(iquad,k) = mass_deficit_v_r
               end if 
            enddo

         end do !iquad

         do k = 1, nlayers

            do iquad = 1, nq

               wq = jac_faceq(iquad,1,iface)

               nxl = normal_vector_q(1,iquad,1,iface)
               nyl = normal_vector_q(2,iquad,1,iface)

               flux = nxl*flux_edge_u(iquad,k) + nyl*flux_edge_v(iquad,k)

               do n = 1, ngl

                  hi = psiq(n,iquad)
                  il = imapl(1,n,1,iface)
                  jl = imapl(2,n,1,iface)
                  kl = imapl(3,n,1,iface)
                  I = intma(il,jl,kl,el)

                  rhs(I,k) = rhs(I,k) - wq*hi*flux
               end do
            end do !iquad
         end do !klayers

         kk=kk+1
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_consistency

 subroutine create_nbhs_face_quad_all(q_face,grad_uvdp_face,q_send,q_recv,nvarb,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
 
   use mod_initial, only: nvar
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, &
                        gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv
 
   use mod_ref, only: nmessage
 
   implicit none
 
   !global arrays
   real, intent(inout) :: q_face(nvarb,2,nq,nface), grad_uvdp_face(nvarb,2,nq,nface)
   real,intent(in):: q_send(2*nvarb,nq,nboun)
   real,intent(in):: q_recv(2*nvarb,nq,nboun)
   integer, intent(in) :: nvarb
 
   !local variables
 
   real :: wq, a_constant, iflux
   integer ::iface,k
   integer :: iel, ier, ilocl, ilocr
   integer :: ifaceb, nq_i, nq_j, plane_ij
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
   integer :: ftype, pface, subface, imulti
   integer :: multirate
 
   !Constants
 
   jj=1
   kk=1
   imm = 0

   ! do ifaceb =1,nboun
   !    iface = face_send(ifaceb)  ! Get Local Face
 
   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)
 
         ftype = face_type(iface)
 
         do im = 1,imulti
 
            ! if (ftype==2) then
               ilocl=face(5,iface)
               iel=face(7,iface)
 
               ! if (multirate==1 .and. lev_list(iel)==0) then
               !    kk=kk+1
               !    cycle
               ! end if
            ! end if
 
            !-------------------------------------
            !Store Left Side Variables
            !-------------------------------------
            call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)
            
            do k=1,nvarb
               !Left Element
               q_face(k,1,:,iface)=q_send(k,:,kk)

               !Right Element
               q_face(k,2,:,iface)=q_recv(k,:,kk)
            end do

            do k=nvarb+1,2*nvarb
               !Left Element
               grad_uvdp_face(k-nvarb,1,:,iface)=q_send(k,:,kk)

               !Right Element
               grad_uvdp_face(k-nvarb,2,:,iface)=q_recv(k,:,kk)
            end do
 
            kk=kk+1

            ! end if
         end do
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_quad_all

 subroutine create_nbhs_face_lap_quad_ip(rhs,q_send,q_recv,nvarb,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq, psiq
 
   use mod_face, only: normal_vector_q, jac_face, imapl, imapr, face_send, jac_faceq
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
 
   use mod_initial, only: nvar
 
   !use mod_metrics, only: jac_faceq
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, &
                        gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv
 
   use mod_ref, only: nmessage
 
   implicit none
 
   !global arrays
   real, intent(inout) :: rhs(2,npoin)
   real,intent(in):: q_send(nvarb,nq,nboun)
   real,intent(in):: q_recv(nvarb,nq,nboun)
   integer, intent(in) :: nvarb
 
   !local variables
 
   real :: wq, a_constant, iflux, nx, ny, hi 
   integer ::iface,k, il, jl, kl, ip, i
   integer :: iel, ier, ilocl, ilocr
   integer :: ifaceb, nq_i, nq_j, plane_ij
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
   integer :: ftype, pface, subface, imulti
   integer :: multirate, iquad
   real, dimension(2) :: qu_mean, qv_mean
   real :: flux_qu, flux_qv
 
   !Constants
 
   jj=1
   kk=1
   imm = 0

   ! do ifaceb =1,nboun
   !    iface = face_send(ifaceb)  ! Get Local Face
 
   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)
 
         ftype = face_type(iface)
 
         do im = 1,imulti
 
            ilocl=face(5,iface)
            iel=face(7,iface)
 
            !-------------------------------------
            !Store Left Side Variables
            !-------------------------------------
            ! call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)
            
            ! do k=1,nvarb
            !    !Left Element
            !    q_face(k,1,:,iface)=q_send(k,:,kk)

            !    !Right Element
            !    q_face(k,2,:,iface)=q_recv(k,:,kk)
            ! end do

            !----------------------------Left Element
            do iquad = 1,nq

                qu_mean(:) = 0.5*(q_send(1:2,iquad,kk) + q_recv(1:2,iquad,kk))
                qv_mean(:) = 0.5*(q_send(3:4,iquad,kk) + q_recv(3:4,iquad,kk))

                wq=jac_faceq(iquad,1,iface)

                !Store normals
                nx = normal_vector_q(1,iquad,1,iface)
                ny = normal_vector_q(2,iquad,1,iface)

                flux_qu = nx*qu_mean(1) + ny*qu_mean(2)
                flux_qv = nx*qv_mean(1) + ny*qv_mean(2)

                !---------------------------------
                !  Do Gauss-Lobatto Integration
                !---------------------------------
                !--------------Left Side------------------!

                do i=1,ngl

                  hi = psiq(i,iquad)

                  !Pointers
                  if(ftype==21) then 
                     il=imapl(1,i,1,iface)
                     jl=imapl(2,i,1,iface)
                     kl=imapl(3,i,1,iface)
                  else
                     il=imapl(1,i,1,iface)
                     jl=imapl(2,i,1,iface)
                     kl=imapl(3,i,1,iface)
                  end if
                  ip=intma(il,jl,kl,iel)
                  
                  !Update Flux
                  rhs(1,ip) = rhs(1,ip) + wq*hi*flux_qu
                  rhs(2,ip) = rhs(2,ip) + wq*hi*flux_qv

               end do !i
            
            end do !iquad
 
            kk=kk+1

            ! end if
         end do
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_lap_quad_ip

 subroutine create_nbhs_face_quad1(q_face,face_sendrecv,q_send,q_recv,nvarb,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
 
   use mod_initial, only: nvar
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, &
                        gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv
 
   use mod_ref, only: nmessage
 
   implicit none
 
   !global arrays
   real, intent(inout) :: q_face(nvarb,2,nq,nface)
   real,intent(in):: q_send(nvarb,nq+1,nboun)
   real,intent(in):: q_recv(nvarb,nq+1,nboun)
   integer, intent(in) :: nvarb
   integer, intent(out) :: face_sendrecv(nvarb,2,nboun)
 
   !local variables
 
   real :: wq, a_constant, iflux
   integer ::iface,k
   integer :: iel, ier, ilocl, ilocr
   integer :: ifaceb, nq_i, nq_j, plane_ij
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
   integer :: ftype, pface, subface, imulti
   integer :: multirate,inode
 
   !Constants
 
   jj=1
   kk=1
   imm = 0

   ! do ifaceb =1,nboun
   !    iface = face_send(ifaceb)  ! Get Local Face
 
   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)
 
         ftype = face_type(iface)
 
         do im = 1,imulti
 
            ! if (ftype==2) then
               ilocl=face(5,iface)
               iel=face(7,iface)
 
               ! if (multirate==1 .and. lev_list(iel)==0) then
               !    kk=kk+1
               !    cycle
               ! end if
            ! end if
 
            !-------------------------------------
            !Store Left Side Variables
            !-------------------------------------
            call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)
            
            do inode = 1,nq
               do k=1,nvarb
                  !Left Element
                  !q_face(k,1,:,iface)=q_send(k,:,kk)

                  !Right Element
                  q_face(k,2,inode,iface)=q_recv(k,inode,jj)
                  face_sendrecv(nvarb,1,jj) = q_send(nvarb,nq+1,jj)
                  face_sendrecv(nvarb,2,jj) = q_recv(nvarb,nq+1,jj)
               end do
            end do

            kk=kk+1

            ! end if
         end do
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_quad1

 subroutine create_nbhs_face_quad_layer(q_face,q_send,q_recv,nvarb,nlayers,nq)

   use mod_basis, only: FACE_CHILDREN
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
 
   use mod_initial, only: nvar
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, &
                        gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                     num_nbh, num_send_recv
 
   use mod_ref, only: nmessage
 
   implicit none
 
   !global arrays
   real, intent(inout) :: q_face(nvarb,2,nq,nface,nlayers)
   real,intent(in):: q_send(nvarb,nq,nboun,nlayers)
   real,intent(in):: q_recv(nvarb,nq,nboun,nlayers)
   integer, intent(in) :: nvarb,nlayers,nq
 
   !local variables
 
   real :: wq, a_constant, iflux
   integer ::iface,k
   integer :: iel, ier, ilocl, ilocr
   integer :: ifaceb, nq_i, nq_j, plane_ij
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
   integer :: ftype, pface, subface, imulti, ll, inode
 
   !Constants
 
   jj=1
   kk=1
   imm = 0
 
   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)
 
         ftype = face_type(iface)
 
         ! do im = 1,imulti
 
            if (ftype==2) then
               ilocl=face(5,iface)
               iel=face(7,iface)
 
            !-------------------------------------
            !Store Left Side Variables
            !-------------------------------------
            call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)
            
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
         ! end do
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_quad_layer

  subroutine create_nbhs_face_quad_layer_all(q_face,qprime_face,q_send,q_recv,nvarb, &
               nlayers,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
 
   use mod_initial, only: nvar
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, &
                        gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv
 
   use mod_ref, only: nmessage
 
   implicit none
 
   !global arrays
   real, dimension(nvarb,2,nq,nface,nlayers), intent(inout) :: q_face, qprime_face
   real,intent(in):: q_send(2*nvarb,nq,nboun,nlayers)
   real,intent(in):: q_recv(2*nvarb,nq,nboun,nlayers)
   integer, intent(in) :: nvarb,nlayers
 
   !local variables
 
   real :: wq, a_constant, iflux
   integer ::iface,k
   integer :: iel, ier, ilocl, ilocr
   integer :: ifaceb, nq_i, nq_j, plane_ij
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
   integer :: ftype, pface, subface, imulti, ll
   integer :: multirate
 
   !Constants
 
   jj=1
   kk=1
   imm = 0
 
   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)
 
         ftype = face_type(iface)

         !print*, 'iface', iface, 'ftype', ftype, 'imulti', imulti
 
         ! do im = 1,imulti
 
            if (ftype==2) then
               ilocl=face(5,iface)
               iel=face(7,iface)
 
            !    if (multirate==1 .and. lev_list(iel)==0) then
            !       kk=kk+1
            !       cycle
            !    end if
            ! end if
 
            !-------------------------------------
            !Store Left Side Variables
            !-------------------------------------
            call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)
            
            do ll = 1,nlayers
               do k=1,nvarb
                  !Left Element
                  q_face(k,1,:,iface,ll)=q_send(k,:,kk,ll)

                  !Right Element
                  q_face(k,2,:,iface,ll)=q_recv(k,:,jj,ll)
               end do

               do k=nvarb+1,2*nvarb
                  !Left Element
                  qprime_face(k-nvarb,1,:,iface,ll)=q_send(k,:,kk,ll)

                  !Right Element
                  qprime_face(k-nvarb,2,:,iface,ll)=q_recv(k,:,jj,ll)
               end do
            end do
 
            kk=kk+1
            end if
         ! end do
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_quad_layer_all

 subroutine create_nbhs_face_quad_1v(q_face,q_send,q_recv,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
 
   use mod_initial, only: nvar
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, &
                        gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, &
                           num_nbh, num_send_recv
 
   use mod_ref, only: nmessage
 
   implicit none
 
   !global arrays
   real, intent(inout) :: q_face(2,nq,nface)
   real,intent(in):: q_send(nq,nboun)
   real,intent(in):: q_recv(nq,nboun)
 
   !local variables
 
   real :: wq, a_constant, iflux
   integer ::iface,k
   integer :: iel, ier, ilocl, ilocr
   integer :: ifaceb, nq_i, nq_j, plane_ij
 
   integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
   integer :: ftype, pface, subface, imulti
   integer :: multirate
 
   !Constants
 
   jj=1
   kk=1
   imm = 0
 
   do inbh = 1, num_nbh
      do ib=1,num_send_recv(inbh)
         iface = nbh_send_recv(jj)
         imulti = nbh_send_recv_multi(jj)
 
         ftype = face_type(iface)

         !print*, 'iface', iface, 'ftype', ftype, 'imulti', imulti
 
         !do im = 1,imulti
 
            if (ftype==2) then
               ilocl=face(5,iface)
               iel=face(7,iface)
 
            !    if (multirate==1 .and. lev_list(iel)==0) then
            !       kk=kk+1
            !       cycle
            !    end if
            ! end if
 
            !-------------------------------------
            !Store Left Side Variables
            !-------------------------------------
            call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)

            !Left Element
            q_face(1,:,iface)=q_send(:,kk)

            !Right Element
            q_face(2,:,iface)=q_recv(:,jj)
 
            kk=kk+1
            end if
         !end do
         jj=jj+1
      end do
 
   end do !iface
 
 end subroutine create_nbhs_face_quad_1v
