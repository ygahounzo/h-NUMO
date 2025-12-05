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
   use mod_input, only: nlayers
   use mod_constants, only: gravity
   use mod_variables, only: H_face_ave,ope_face_ave,btp_mass_flux_face_ave, &
                                Qu_face_ave, Qv_face_ave, one_plus_eta_edge_2_ave, &
                                Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, &
                                uvb_face_ave, H_bcl_edge, ope2_face_ave
 
   implicit none
 
   !global arrays
   real, intent(inout) :: rhs(3,npoin)
   real,intent(in):: q_send(nbtp_var,ngl,nboun)
   real,intent(in):: q_recv(nbtp_var,ngl,nboun)
   integer, intent(in) :: nvarb
 
   !local variables
 
   real :: wq
   integer :: iface, k
 
   integer :: jj, imm, inbh, ib, kk, ivar, ii
   integer :: imulti
   integer :: multirate
   real, dimension(4,nq) :: qbl, qbr
   real, dimension(nq) :: pbl, pbr, H_face_tmp
   real :: nxl, nyl, nxr, nyr
   real :: pbpert_edge, c_minus, c_plus
   real :: c_pb_L, c_pb_R, c_pbub_LR, c_pbub_L, c_pbub_R
   real, dimension(nq) :: one_plus_eta_edge, c_pb_LR
   real, dimension(nq) :: flux_edge_x, flux_edge_y
   integer :: el, iquad, n, I, il, jl, kl, er, itype
   real, dimension(nq) :: quu, quv, qvu, qvv, ul, ur, vl, vr
   real :: pU_L, pU_R, lamb, dispu, dispv, flux_x, flux_y, flux
   real :: hi, H_kx, H_ky
   real, dimension(nlayers+1) :: pprime_l, pprime_r
   real, dimension(nq) :: Q_uu_q, Q_uv_q, Q_vu_q, Q_vv_q, H_bcl_q
   real :: ppl, ppr, upl, upr, vpl, vpr

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

         qbl = 0.0 ; qbr = 0.0
         pbl = 0.0 ; pbr = 0.0
         Q_uu_q = 0.0; Q_uv_q = 0.0; Q_vu_q = 0.0; Q_vv_q = 0.0 ; H_bcl_q = 0.0

         do iquad = 1, nq

            nxl = normal_vector_q(1,iquad,1,iface)
            nyl = normal_vector_q(2,iquad,1,iface)

            nxr = -nxl; nyr = -nyl

            do n = 1,ngl
               hi = psiq(n,iquad)
               !Left Element
               qbl(1:4,iquad) = qbl(1:4,iquad) + hi*q_send(1:4,n,kk)
               pbl(iquad) = pbl(iquad) + hi*q_send(5,n,kk)
               !Right Element
               qbr(1:4,iquad) = qbr(1:4,iquad) + hi*q_recv(1:4,n,kk)
               pbr(iquad) = pbr(iquad) + hi*q_recv(5,n,kk)
            end do

            pprime_l(1) = 0.0; pprime_r(1) = 0.0
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

               Q_uu_q(iquad) = Q_uu_q(iquad) + 0.5*alpha_mlswe(k)*((upl*(upl*ppl)) + (upr*(upr*ppr)))
               Q_uv_q(iquad) = Q_uv_q(iquad) + 0.5*alpha_mlswe(k)*((vpl*(upl*ppl)) + (vpr*(upr*ppr)))
               Q_vu_q(iquad) = Q_vu_q(iquad) + 0.5*alpha_mlswe(k)*((upl*(vpl*ppl)) + (upr*(vpr*ppr)))
               Q_vv_q(iquad) = Q_vv_q(iquad) + 0.5*alpha_mlswe(k)*((vpl*(vpl*ppl)) + (vpr*(vpr*ppr)))

               pprime_l(k+1) = pprime_l(k) + ppl
               pprime_r(k+1) = pprime_r(k) + ppr

               H_bcl_q(iquad) = H_bcl_q(iquad) + 0.25*alpha_mlswe(k)*((pprime_l(k+1)**2 - pprime_l(k)**2) + (pprime_r(k+1)**2 - pprime_r(k)**2))
            end do

            pU_L = nxl * qbl(3,iquad) + nyl * qbl(4,iquad)
            pU_R = nxr * qbr(3,iquad) + nyr * qbr(4,iquad)

            c_minus = sqrt(alpha_mlswe(nlayers) * pbr(iquad))
            c_plus  = sqrt(alpha_mlswe(nlayers) * pbl(iquad))
            c_pb_L = c_minus / (c_minus + c_plus)
            c_pb_R = c_plus / (c_minus + c_plus)
            c_pbub_LR  = 1.0 / (c_minus + c_plus)

            pbpert_edge = c_pb_L * qbl(2,iquad) &
                                       + c_pb_R * qbr(2,iquad) &
                                       + c_pbub_LR * (pU_L + pU_R)

            one_plus_eta_edge(iquad) = 1.0 + (pbpert_edge/pbl(iquad))

            ! Compute mass fluxes at each element face.

            c_pbub_L = c_minus / (c_minus + c_plus)
            c_pbub_R = c_plus / (c_minus + c_plus)
            c_pb_LR(iquad) = c_minus * c_plus / (c_minus + c_plus)

            flux_edge_x(iquad) = c_pbub_L * qbl(3,iquad) &
                                       + c_pbub_R * qbr(3,iquad) &
                                       + c_pb_LR(iquad) * &
                                       (nxl * qbl(2,iquad) + nxr * qbr(2,iquad))

            flux_edge_y(iquad) = c_pbub_L * qbl(4,iquad) &
                                       + c_pbub_R * qbr(4,iquad) &
                                       + c_pb_LR(iquad) * &
                                       (nyl * qbl(2,iquad) + nyr * qbr(2,iquad))

         end do !iquad

         ul = qbl(3,:)/qbl(1,:); ur = qbr(3,:)/qbr(1,:)
         vl = qbl(4,:)/qbl(1,:); vr = qbr(4,:)/qbr(1,:)

         quu(:) = 0.5*(ul*qbl(3,:) + ur*qbr(3,:)) + one_plus_eta_edge(:) * Q_uu_q(:)
         quv(:) = 0.5*(vl*qbl(3,:) + vr*qbr(3,:)) + one_plus_eta_edge(:) * Q_uv_q(:)
         qvu(:) = 0.5*(ul*qbl(4,:) + ur*qbr(4,:)) + one_plus_eta_edge(:) * Q_vu_q(:)
         qvv(:) = 0.5*(vl*qbl(4,:) + vr*qbr(4,:)) + one_plus_eta_edge(:) * Q_vv_q(:)

         ! Compute pressure forcing H_face at each element face.
         H_face_tmp(:) = (one_plus_eta_edge(:)**2) * H_bcl_q(:)

         ! Accumulate sums for time averaging

         btp_mass_flux_face_ave(1,:,iface) = btp_mass_flux_face_ave(1,:,iface) + flux_edge_x(:)
         btp_mass_flux_face_ave(2,:,iface) = btp_mass_flux_face_ave(2,:,iface) + flux_edge_y(:)
         H_face_ave(:,iface) = H_face_ave(:,iface) + H_face_tmp(:)
         Qu_face_ave(1,:,iface) = Qu_face_ave(1,:,iface) + quu(:)
         Qu_face_ave(2,:,iface) = Qu_face_ave(2,:,iface) + quv(:)
         Qv_face_ave(1,:,iface) = Qv_face_ave(1,:,iface) + qvu(:)
         Qv_face_ave(2,:,iface) = Qv_face_ave(2,:,iface) + qvv(:)
         ope_face_ave(1,:,iface) = ope_face_ave(1,:,iface) + (1.0 + (qbl(2,:)/pbl))
         ope_face_ave(2,:,iface) = ope_face_ave(2,:,iface) + (1.0 + (qbr(2,:)/pbr))
         ope2_face_ave(1,:,iface) = ope2_face_ave(1,:,iface) + (1.0 + (qbl(2,:)/pbl))**2
         ope2_face_ave(2,:,iface) = ope2_face_ave(2,:,iface) + (1.0 + (qbr(2,:)/pbr))**2
         one_plus_eta_edge_2_ave(:,iface) = one_plus_eta_edge_2_ave(:,iface) &
                                             + one_plus_eta_edge(:)**2
         uvb_face_ave(1,1,:,iface) = uvb_face_ave(1,1,:,iface) + ul
         uvb_face_ave(1,2,:,iface) = uvb_face_ave(1,2,:,iface) + ur
         uvb_face_ave(2,1,:,iface) = uvb_face_ave(2,1,:,iface) + vl
         uvb_face_ave(2,2,:,iface) = uvb_face_ave(2,2,:,iface) + vr

         do iquad = 1, nq

            wq = jac_faceq(iquad,1,iface)

            nxl = normal_vector_q(1,iquad,1,iface)
            nyl = normal_vector_q(2,iquad,1,iface)

            H_kx = nxl*H_face_tmp(iquad)
            H_ky = nyl*H_face_tmp(iquad)

            lamb = c_pb_LR(iquad)

            dispu = 0.5*lamb*(qbr(3,iquad) - qbl(3,iquad))
            dispv = 0.5*lamb*(qbr(4,iquad) - qbl(4,iquad))
            flux_x = nxl*quu(iquad) + nyl*quv(iquad) - dispu
            flux_y = nxl*qvu(iquad) + nyl*qvv(iquad) - dispv

            flux = nxl*flux_edge_x(iquad) + nyl*flux_edge_y(iquad) 

            do n = 1, ngl

               hi = psiq(n,iquad)
               il = imapl(1,n,1,iface)
               jl = imapl(2,n,1,iface)
               kl = imapl(3,n,1,iface)
               I = intma(il,jl,kl,el)

               rhs(1,I) = rhs(1,I) - wq*hi*flux
               rhs(2,I) = rhs(2,I) - wq*hi*(H_kx + flux_x)
               rhs(3,I) = rhs(3,I) - wq*hi*(H_ky + flux_y)
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

            nxl = normal_vector(1,iquad,1,iface)
            nyl = normal_vector(2,iquad,1,iface)

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
