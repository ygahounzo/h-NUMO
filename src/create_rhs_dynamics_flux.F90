

!----------------------------------------------------------------------!
!>@brief Constructs the Inter-Processor Fluxes that lie on the boundary of each processor element
!>@details Requires calling send_receive_boundary, which gets the data lying off-processor and
!> populates q_send and q_recv
!>@author James F. Kelly
!>@date 16 November 2010
!>@date September 7, 2015 by F.X. Giraldo to reuse new ELEMENTAL_FLUX and RUSANOV 
!>routines that use long vectors
!>@date January 2016 by M.A. Kopera - accounts for non-conforming faces 
!----------------------------------------------------------------------!
subroutine create_rhs_dynamics_flux_ip(rhs,q_send,q_recv,multirate)

  use mod_basis, only: ngl, FACE_CHILDREN

  use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send

  use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_ngl, face_type

  use mod_initial, only: nvar

  use mod_input, only: space_method, cgdg_method, llimit, limit_threshold, is_shallow, form_method

  use mod_metrics, only: jac

  use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list

  use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv

  use mod_ref, only: nmessage

  implicit none

  !global arrays
  real, intent(inout) :: rhs(nvar,npoin)
  real,intent(in):: q_send(nmessage,ngl,ngl,nboun)
  real,intent(in):: q_recv(nmessage,ngl,ngl,nboun)

  !local variables
  real, dimension(nvar,ngl,ngl) :: flux_q, flux_ql, flux_qr, flux_lq
  real, dimension(nvar+2,ngl,ngl) :: ql, qr, qlp
  real, dimension(ngl,ngl) :: nxl, nyl, nzl, nxr, nyr, nzr
  real, dimension(ngl,ngl,4) :: nxrp, nyrp, nzrp
  real, dimension(ngl,ngl) :: jac_face_l

  real :: wq, a_constant, iflux
  integer ::iface, i, j, k, l, ip, il, jl, kl, ir, jr, kr
  integer :: iel, ier, ilocl, ilocr
  integer :: ifaceb, ngl_i, ngl_j, plane_ij

  integer :: isub, ic, jj, im, imm, inbh, ib, kk, ivar
  integer :: ftype, pface, subface, imulti
  logical :: compute_mixed_flux
  integer :: multirate

  compute_mixed_flux = (space_method == 'cgd' .and. (cgdg_method == 'mixed'))
  if (compute_mixed_flux) then
     return
  end if

  !Constants
  a_constant=1
  iflux=1.0 !strong form
  if (is_shallow) then
     iflux=0.0 !weak form
  end if
  if (form_method(1:4) == 'weak') iflux=0.0 !weak-form FLUX_IP (only called by DG)

  pface = 0
  isub = 1

  !Loop over all boundary faces
  !  do ifaceb =1,nboun
  !     iface = face_send(ifaceb)  ! Get Local Face

  jj=1
  kk=1
  imm = 0

  do inbh = 1, num_nbh
     do ib=1,num_send_recv(inbh)
        iface = nbh_send_recv(jj)
        imulti = nbh_send_recv_multi(jj)

        ftype = face_type(iface)

        do im = 1,imulti

           if (ftype==2) then
              ilocl=face(5,iface)
              iel=face(7,iface)

              if (multirate==1 .and. lev_list(iel)==0) then
                 kk=kk+1
                 cycle
              end if

           else if (ftype==21) then
              ilocl=face(6,iface)
              i=1
              ic=0
              do while(ic<im)
                 if (face(7+i,iface)>0) then
                    ic=ic+1
                    iel=face(7+i,iface)
                    subface = i
                 end if
                 i=i+1
              end do
           else if (ftype==12) then
              ilocl=face(5,iface)
              i=1
              ic=0
              do while(ic<im)
                 if (nbh_send_recv_half((jj-1)*FACE_CHILDREN+i)==1) then
                    ic=ic+1
                    subface=i
                 end if
                 i=i+1
              end do
              iel=face(7,iface)
           end if

           !-------------------------------------
           !Store Left Side Variables
           !-------------------------------------
           call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

           !Left Element
           do k=1,nvar+2
              ql(k,:,:)=q_send(k,:,:,kk)
           end do

           !Right Element
           do k=1,nvar+2
              qr(k,:,:)=q_recv(k,:,:,kk)
           end do

           !project parent face
           if (ftype==12) then
              qlp=ql
              do ivar=1,nvar+2
                 call scatter_element_2d_subface(qlp(ivar,:,:),ql(ivar,:,:),subface, ngl_i, ngl_j, plane_ij)
              end do
           else if (ftype==21) then
              qlp=qr
              do ivar=1,nvar+2
                 call scatter_element_2d_subface(qlp(ivar,:,:),qr(ivar,:,:),subface, ngl_i, ngl_j, plane_ij)
              end do
           end if

           !------------------
           ! Store normals
           !------------------
           if (ftype==2) then
              do j=1,ngl_j
                 do i=1,ngl_i
                    !Left Normals
                    nxl(i,j)=normal_vector(1,i,j,iface)
                    nyl(i,j)=normal_vector(2,i,j,iface)
                    nzl(i,j)=normal_vector(3,i,j,iface)
                    jac_face_l(i,j) = jac_face(i,j,iface)
                 end do !i
              end do !j
           else if (ftype==21) then
              call compute_normals_element(nxl,nyl,nzl,jac_face_l,iel,ilocl)
           else if (ftype==12) then
              ! get normals and project them
              do i=1,ngl_i
                 do j=1,ngl_j
                    !Left Normals
                    nxl(i,j)=normal_vector(1,i,j,iface)
                    nyl(i,j)=normal_vector(2,i,j,iface)
                    nzl(i,j)=normal_vector(3,i,j,iface)
                    jac_face_l(i,j) = jac_face(i,j,iface)
                 end do !j
              end do !i
              !------------------------------------------
              ! Project normals
              !------------------------------------------
              call scatter_element_2d(nxl,nxrp, ngl_i, ngl_j, plane_ij)
              call scatter_element_2d(nyl,nyrp, ngl_i, ngl_j, plane_ij)
              call scatter_element_2d(nzl,nzrp, ngl_i, ngl_j, plane_ij)

              nxl = nxrp(:,:,subface)
              nyl = nyrp(:,:,subface)
              nzl = nzrp(:,:,subface)

           end if
           nxr=-nxl
           nyr=-nyl
           nzr=-nzl

           !----------------------------------
           !  Compute face and element fluxes
           !----------------------------------
           call elemental_flux(flux_ql,ql,nxl,nyl,nzl,ngl_i,ngl_j)
           call elemental_flux(flux_qr,qr,nxr,nyr,nzr,ngl_i,ngl_j)
           call rusanov_flux(flux_q,flux_ql,flux_qr,ql,qr,nxl,nyl,nzl,a_constant,ngl_i,ngl_j)

           !---------------------------------
           !  Do Gauss-Lobatto Integration
           !---------------------------------

           ! do j=1,ngl
           !    write(*,'(4e16.8,XXX,4e16.8,XXX,4e16.8,XXX,4e16.8)') ql(1,:,j), qr(1,:,j), flux_q(1,:,j), flux_1l(1,:,j)
           ! end do


           if (ftype==12) then
              !need to backward-project the flux
              do ivar=1,nvar
                 call gather_element_2d_subface(flux_q(ivar,:,:),flux_lq(ivar,:,:),subface, ngl_i, ngl_j, plane_ij)
              end do
              flux_q=flux_lq

              !elemental flux needs to be computed differently for face 12 - will do it in the regular nc subroutine
              flux_ql = 0
           end if

           do j=1,ngl_j
              do i=1,ngl_i
                 !Pointers
                 if (ftype==21) then
                    il=imapr(1,i,j,iface)
                    jl=imapr(2,i,j,iface)
                    kl=imapr(3,i,j,iface)
                 else
                    il=imapl(1,i,j,iface)
                    jl=imapl(2,i,j,iface)
                    kl=imapl(3,i,j,iface)
                 end if
                 wq=jac_face_l(i,j)
                 ip=intma(il,jl,kl,iel)

                 !Store RHS
                 do k=1,nvar
                    rhs(k,ip)=rhs(k,ip) + wq*(iflux*flux_ql(k,i,j) - flux_q(k,i,j))
                 end do

              end do !i
           end do !j

           kk=kk+1
        end do
        jj=jj+1
     end do

  end do !iface

end subroutine create_rhs_dynamics_flux_ip


subroutine create_nbhs_face_quad(q_face,q_send,q_recv,nvarb,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
 
   use mod_initial, only: nvar
 
   use mod_input, only: space_method, cgdg_method, llimit, limit_threshold, is_shallow, form_method
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv
 
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

 subroutine create_nbhs_face_df(q_face,q_send,q_recv,nvarb,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
 
   use mod_initial, only: nvar
 
   use mod_input, only: space_method, cgdg_method, llimit, limit_threshold, is_shallow, form_method
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv
 
   use mod_ref, only: nmessage
 
   implicit none
 
   !global arrays
   real, intent(inout) :: q_face(nvarb,2,ngl,nface)
   real,intent(in):: q_send(nvarb,ngl,nboun)
   real,intent(in):: q_recv(nvarb,ngl,nboun)
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
 
            !-------------------------------------
            !Store Left Side Variables
            !-------------------------------------
            
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
 
 end subroutine create_nbhs_face_df

 subroutine create_nbhs_face_quad_all(q_face,grad_uvdp_face,q_send,q_recv,nvarb,multirate)

   use mod_basis, only: ngl, FACE_CHILDREN,nq
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
 
   use mod_initial, only: nvar
 
   use mod_input, only: space_method, cgdg_method, llimit, limit_threshold, is_shallow, form_method
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv
 
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
 
   use mod_input, only: space_method, cgdg_method, llimit, limit_threshold, is_shallow, form_method
 
   !use mod_metrics, only: jac_faceq
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv
 
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
 
   use mod_input, only: space_method, cgdg_method, llimit, limit_threshold, is_shallow, form_method
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv
 
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
 
   use mod_input, only: space_method, cgdg_method, llimit, limit_threshold, is_shallow, form_method
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv
 
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

         !print*, 'iface', iface, 'ftype', ftype, 'imulti', imulti
 
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

  subroutine create_nbhs_face_quad_layer_all(q_face,qprime_face,q_send,q_recv,nvarb,nlayers)

   use mod_basis, only: ngl, FACE_CHILDREN,nq
 
   use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send
 
   use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_nq, face_type,nface
 
   use mod_initial, only: nvar
 
   use mod_input, only: space_method, cgdg_method, llimit, limit_threshold, is_shallow, form_method
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv
 
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
 
   use mod_input, only: space_method, cgdg_method, llimit, limit_threshold, is_shallow, form_method
 
   use mod_metrics, only: jac
 
   use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list
 
   use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv
 
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

!----------------------------------------------------------------------!
!>@brief This subroutine builds the RUSANOV flux for both WEAK and STRONG as well 
!> as for CONFORMING or NON-CONFORMING grids.
!>@author M.A. Kopera to do weak form.
!>@date F.X. Giraldo to do both strong and weak.
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@date November 07, 2014 Daniel S. Abdi
!>@date Sepetember 07, 2015 F.X. Giraldo, rewritten to use long vectors and to 
!>reuse ELEMENTAL_FLUX routines that computes FLUX_QL and FLUX_QR. All we need is the
!>JUMPS in this routine => DISS_Q
!----------------------------------------------------------------------!
subroutine rusanov_flux(flux_q,flux_ql,flux_qr,ql,qr,nx,ny,nz,a_constant,ngl_i,ngl_j)

  use mod_basis, only: ngl

  use mod_constants, only: gamma, gravity

  use mod_initial, only: nvar, rho_layers

  use mod_input, only: is_shallow, eqn_set, is_swe_layers

  implicit none

  !global
  real, dimension(nvar,ngl,ngl), intent(out) :: flux_q
  real, dimension(nvar,ngl,ngl), intent(in)  :: flux_ql, flux_qr
  real, dimension(nvar+2,ngl,ngl), intent(in) :: ql, qr
  real, dimension(ngl,ngl), intent(in) :: nx, ny, nz
  real, intent(in) :: a_constant
  integer, intent(in) :: ngl_i, ngl_j

  !local
  real :: rl_k, ul_k, vl_k, wl_k, tl_k, plt_k, al_k
  real :: rul_k, rvl_k, rwl_k, rtl_k, rql_k
  real :: rr_k, ur_k, vr_k, wr_k, tr_k, prt_k, ar_k
  real :: rur_k, rvr_k, rwr_k, rtr_k, rqr_k
  real :: unl, unr, claml, clamr, clam, clamu, diss_q
  real :: nxl, nyl, nzl
  integer :: i, j, k
  logical :: is_set2nc

  !Equation Set Used
  is_set2nc = (eqn_set == 'set2nc')

  !Do Gauss-Lobatto Integration
  do j=1,ngl_j
     do i=1,ngl_i

        !Store Normal Vectors
        nxl=nx(i,j)
        nyl=ny(i,j)
        nzl=nz(i,j)

        !Store Left and Right States
        rl_k=ql(1,i,j)
        ul_k=ql(2,i,j)
        vl_k=ql(3,i,j)
        wl_k=ql(4,i,j)
        tl_k=ql(5,i,j)
        plt_k=ql(nvar+2,i,j)

        rr_k=qr(1,i,j)
        ur_k=qr(2,i,j)
        vr_k=qr(3,i,j)
        wr_k=qr(4,i,j)
        tr_k=qr(5,i,j)
        prt_k=qr(nvar+2,i,j)

        !Construct Conservation Variables
        if (is_set2nc) then
           rul_k=ul_k
           rvl_k=vl_k
           rwl_k=wl_k
           rtl_k=tl_k
           rur_k=ur_k
           rvr_k=vr_k
           rwr_k=wr_k
           rtr_k=tr_k
        else !set2c or set3c
           rul_k=rl_k*ul_k
           rvl_k=rl_k*vl_k
           rwl_k=rl_k*wl_k
           rtl_k=rl_k*tl_k
           rur_k=rr_k*ur_k
           rvr_k=rr_k*vr_k
           rwr_k=rr_k*wr_k
           rtr_k=rr_k*tr_k
        end if

        !Compute Rusanov flux Constant
        unl=nxl*ul_k + nyl*vl_k + nzl*wl_k
        unr=nxl*ur_k + nyl*vr_k + nzl*wr_k
        al_k = sqrt(gamma*abs(plt_k/rl_k))
        ar_k = sqrt(gamma*abs(prt_k/rr_k))

        claml=abs(unl) + a_constant*al_k
        clamr=abs(unr) + a_constant*ar_k
        clam=max(claml,clamr)
        clamu=max(abs(unl),abs(unr))
       

        !----Mass Equation-----!
        diss_q=clam*(rr_k - rl_k)
        flux_q(1,i,j)=0.5*( flux_ql(1,i,j) - flux_qr(1,i,j) - diss_q)

        !----U-Momentum Equation-----!
        diss_q=clam*(rur_k - rul_k)
        flux_q(2,i,j)=0.5*( flux_ql(2,i,j) - flux_qr(2,i,j) - diss_q)

        !----V-Momentum Equation-----!
        diss_q=clam*(rvr_k - rvl_k)
        flux_q(3,i,j)=0.5*( flux_ql(3,i,j) - flux_qr(3,i,j) - diss_q)

        !----W-Momentum Equation-----!
        diss_q=clam*(rwr_k - rwl_k)
        flux_q(4,i,j)=0.5*( flux_ql(4,i,j) - flux_qr(4,i,j) - diss_q)

        !----Temperature Equation-----!
        !Flux Variables
        diss_q=clam*(rtr_k - rtl_k)
        flux_q(5,i,j)=0.5*( flux_ql(5,i,j) - flux_qr(5,i,j) - diss_q )

        !----Tracer Equations-----!
        do k=6,nvar
           rql_k=ql(1,i,j)*ql(k,i,j)
           rqr_k=qr(1,i,j)*qr(k,i,j)
           diss_q=clamu*(rqr_k - rql_k)
           flux_q(k,i,j)=0.5*( flux_ql(k,i,j) - flux_qr(k,i,j) - diss_q )
        end do !k

     end do !i
  end do !j

end subroutine rusanov_flux

!----------------------------------------------------------------------!
!>@brief This subroutine builds the Local Elemental flux 
!>@author F.X. Giraldo 
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@date November 07, 2014 Daniel S. Abdi
!----------------------------------------------------------------------!
subroutine elemental_flux(flux_ql,ql,nx,ny,nz,ngl_i,ngl_j)

  use mod_basis, only: ngl

  use mod_constants, only: gravity
  
  use mod_initial, only: nvar, rho_layers

  use mod_input, only: eqn_set, is_shallow, nonlinear_swe, is_swe_layers
  
  implicit none

  !global
  real, dimension(nvar,ngl,ngl), intent(out) :: flux_ql
  real, dimension(nvar+2,ngl,ngl), intent(in) :: ql
  real, dimension(ngl,ngl), intent(in) :: nx, ny, nz
  integer, intent(in) :: ngl_i, ngl_j

  !local
  real :: nxl, nyl, nzl
  real :: rl_k, ul_k, vl_k, wl_k, tl_k, pl_k, plt_k
  real :: rul_k, rvl_k, rwl_k, rtl_k, rql_k
  real :: fxl, fyl, fzl
  integer :: i, j, k
  logical :: is_set2c, is_set2nc, is_set3c

  is_set2nc = (eqn_set == 'set2nc')
  is_set2c =  (eqn_set == 'set2c')
  is_set3c =  (eqn_set == 'set3c')

  !Do Gauss-Lobatto Integration
  do j=1,ngl_j
     do i=1,ngl_i

        !Store Normal Vectors
        nxl=nx(i,j)
        nyl=ny(i,j)
        nzl=nz(i,j)

        !Store Left and Right States
        rl_k=ql(1,i,j)
        ul_k=ql(2,i,j)
        vl_k=ql(3,i,j)
        wl_k=ql(4,i,j)
        tl_k=ql(5,i,j)
        pl_k=ql(nvar+1,i,j)
        plt_k = ql(nvar+2,i,j)

        !Construct Conservation Variables
        if (is_set2nc) then
           rul_k=ul_k
           rvl_k=vl_k
           rwl_k=wl_k
!           rtl_k=tl_k
!           pl_k=pl_k/rl_k !FXG SET2NC: doesn't matter for CG
           pl_k=0
           rtl_k=0 !FXG SET2NC: without this, the temperature has large overshoots
        else !set2c or set3c
           rul_k=rl_k*ul_k
           rvl_k=rl_k*vl_k
           rwl_k=rl_k*wl_k
           rtl_k=rl_k*tl_k
        end if

        !----Mass Equation-----!
        fxl=rl_k*ul_k
        fyl=rl_k*vl_k
        fzl=rl_k*wl_k
        flux_ql(1,i,j)=( nxl*fxl + nyl*fyl + nzl*fzl )

        !----U-Momentum Equation-----!
        fxl=nonlinear_swe*rul_k*ul_k + pl_k
        fyl=nonlinear_swe*rul_k*vl_k
        fzl=nonlinear_swe*rul_k*wl_k
        flux_ql(2,i,j)=( nxl*fxl + nyl*fyl + nzl*fzl )

        !----V-Momentum Equation-----!
        fxl=nonlinear_swe*rvl_k*ul_k
        fyl=nonlinear_swe*rvl_k*vl_k + pl_k
        fzl=nonlinear_swe*rvl_k*wl_k
        flux_ql(3,i,j)=( nxl*fxl + nyl*fyl + nzl*fzl )

        !----W-Momentum Equation-----!
        fxl=nonlinear_swe*rwl_k*ul_k
        fyl=nonlinear_swe*rwl_k*vl_k
        fzl=nonlinear_swe*rwl_k*wl_k + pl_k
        flux_ql(4,i,j)=( nxl*fxl + nyl*fyl + nzl*fzl )

        !----Temperature Equation-----!
        if (is_set3c) then
           fxl=(rtl_k+plt_k)*ul_k
           fyl=(rtl_k+plt_k)*vl_k
           fzl=(rtl_k+plt_k)*wl_k
        else !set2c or set2nc
           fxl=rtl_k*ul_k
           fyl=rtl_k*vl_k
           fzl=rtl_k*wl_k
        endif
        flux_ql(5,i,j)=( nxl*fxl + nyl*fyl + nzl*fzl )

        !----Tracer Equations-----!
        do k=6,nvar
           rql_k=ql(1,i,j)*ql(k,i,j)
           fxl=rql_k*ul_k
           fyl=rql_k*vl_k
           fzl=rql_k*wl_k
           flux_ql(k,i,j)=( nxl*fxl + nyl*fyl + nzl*fzl )
        end do !k

     end do !i
  end do !j

end subroutine elemental_flux


