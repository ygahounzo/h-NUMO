!----------------------------------------------------------------------!
!>@brief Communication using Non-Blocking Send/Receives
!>Comes from NSEAM code.
!>@author  F.X. Giraldo
!>      Department of Applied Mathematics
!>      Naval Postgraduate School
!>      Monterey, CA 93943
!----------------------------------------------------------------------!

subroutine unpack_data_dg_general_quad(G, b, par, q_send, q_recv, send_data, recv_data, nvarb)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb

   real, dimension(nvarb,b%nq,par%num_send_recv_total), intent(out) :: q_send, q_recv
   real, dimension(nvarb*b%nq*par%num_send_recv_total), intent(in)  :: send_data, recv_data

   integer :: ii, jj, kk, i, inbh, ib, ifaces, inode, jnode, ivar, ilocl, ilocr
   integer :: nq_i, nq_j, plane_ij, iface, imulti, ftype

   ii = 0
   jj = 1
   kk = 1

   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)
         ftype  = G%face_type(iface)

         if (ftype == 2 .and. imulti > 0) then
            do inode = 1, b%nq
               do ivar = 1, nvarb
                  ii = ii + 1
                  q_send(ivar,inode,kk) = send_data(ii)
                  q_recv(ivar,inode,kk) = recv_data(ii)
               end do
            end do
            kk = kk + 1
         end if
         jj = jj + 1
      end do
   end do

end subroutine unpack_data_dg_general_quad

subroutine unpack_data_dg_general_df(G, b, par, q_send, q_recv, send_data, recv_data, nvarb, nboun_valid)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb
   integer,           intent(in) :: nboun_valid

   real, dimension(nvarb,b%ngl,G%nboun), intent(out) :: q_send, q_recv
   real, dimension(nvarb*b%ngl*G%nboun), intent(in)  :: send_data, recv_data

   integer :: ii, kk, inode, ivar
   integer :: ngl_f, nvarb_f

   ngl_f   = b%ngl
   nvarb_f = nvarb

   ! The send/recv buffers are packed in compact order (kk=1..nboun_valid).
   ! nboun_valid is pre-computed in ref and passed by the caller — no pre-scan needed.
   !$acc parallel loop gang vector collapse(3) &
   !$acc    present(q_send, q_recv, send_data, recv_data) &
   !$acc    firstprivate(ngl_f, nvarb_f, nboun_valid)
   do kk = 1, nboun_valid
      do inode = 1, ngl_f
         do ivar = 1, nvarb_f
            ii = (kk-1)*ngl_f*nvarb_f + (inode-1)*nvarb_f + ivar
            q_send(ivar, inode, kk) = send_data(ii)
            q_recv(ivar, inode, kk) = recv_data(ii)
         end do
      end do
   end do
   !$acc end parallel loop

end subroutine unpack_data_dg_general_df

subroutine unpack_data_dg_general_df_v0(G, b, par, q_send, q_recv, send_data, recv_data, nvarb)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb

   real, dimension(nvarb,b%ngl,G%nboun), intent(out) :: q_send, q_recv
   real, dimension(nvarb*b%ngl*G%nboun), intent(in)  :: send_data, recv_data

   integer :: jj, inode, ivar, ii, ngl_f, nvarb_f, nboun_f

   ! All num_send_recv_total faces have face_type==2 and imulti==1 (no AMR),
   ! so the flat offset (jj-1)*ngl*nvarb mirrors the GPU pack stride exactly.
   ngl_f   = b%ngl
   nvarb_f = nvarb
   nboun_f = G%nboun

   !$acc parallel loop gang vector collapse(3) &
   !$acc    present(q_send, q_recv, send_data, recv_data) &
   !$acc    firstprivate(ngl_f, nvarb_f, nboun_f)
   do jj = 1, nboun_f
      do inode = 1, ngl_f
         do ivar = 1, nvarb_f
            ii = (jj-1)*ngl_f*nvarb_f + (inode-1)*nvarb_f + ivar
            q_send(ivar, inode, jj) = send_data(ii)
            q_recv(ivar, inode, jj) = recv_data(ii)
         end do
      end do
   end do
   !$acc end parallel loop

end subroutine unpack_data_dg_general_df_v0

subroutine unpack_data_dg_general_lap(G, b, par, q_send, q_recv, send_data, recv_data, nvarb)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb

   real, dimension(10,b%ngl,G%nboun), intent(out) :: q_send, q_recv
   real, dimension(10*b%ngl*G%nboun), intent(in)  :: send_data, recv_data

   integer :: ii, kk, inode, ivar, ngl_f, nboun_f

   ngl_f   = b%ngl
   nboun_f = G%nboun

   !$acc parallel loop gang vector collapse(3) &
   !$acc    present(q_send, q_recv, send_data, recv_data) &
   !$acc    firstprivate(ngl_f, nboun_f)
   do kk = 1, nboun_f
      do inode = 1, ngl_f
         do ivar = 1, 10
            ii = (kk-1)*ngl_f*10 + (inode-1)*10 + ivar
            q_send(ivar, inode, kk) = send_data(ii)
            q_recv(ivar, inode, kk) = recv_data(ii)
         end do
      end do
   end do
   !$acc end parallel loop

end subroutine unpack_data_dg_general_lap

subroutine unpack_data_dg_general_bcl(G, b, inp, par, q_send, q_recv, send_data, recv_data)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS
   use mod_input,    only: input

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(input),       intent(in) :: inp
   type(parallel_CS), intent(in) :: par

   real, dimension(3*inp%nlayers,b%ngl,par%num_send_recv_total), intent(out) :: q_send, q_recv
   real, dimension(3*inp%nlayers*b%ngl*par%num_send_recv_total), intent(in)  :: send_data, recv_data

   integer :: ii, jj, kk, i, inbh, ib, ifaces, inode, jnode, ivar, ilocl, ilocr
   integer :: nq_i, nq_j, plane_ij, iface, imulti, ftype, index, ll

   ii = 0
   jj = 1
   kk = 1

   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then
            do ll = 1, inp%nlayers
               index = (ll-1)*3
               do inode = 1, b%ngl
                  do ivar = 1, 3
                     ii = ii + 1
                     q_send(index+ivar,inode,kk) = send_data(ii)
                     q_recv(index+ivar,inode,kk) = recv_data(ii)
                  end do
               end do
            end do
            kk = kk + 1
         end if
         jj = jj + 1
      end do
   end do

end subroutine unpack_data_dg_general_bcl

subroutine unpack_data_dg_general_lap_bcl(G, b, par, q_send, q_recv, send_data, recv_data, nlayers)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nlayers

   real, dimension(5*nlayers,b%ngl,par%num_send_recv_total), intent(out) :: q_send, q_recv
   real, dimension(5*nlayers*b%ngl*par%num_send_recv_total), intent(in)  :: send_data, recv_data

   integer :: ii, jj, kk, i, inbh, ib, ifaces, inode, jnode, ivar, ilocl, ilocr
   integer :: nq_i, nq_j, plane_ij, iface, imulti, ftype, ll, index

   ii = 0
   jj = 1
   kk = 1

   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then
            do ll = 1, nlayers
               index = (ll-1)*5
               do inode = 1, b%ngl
                  do ivar = 1, 5
                     ii = ii + 1
                     q_send(index+ivar,inode,kk) = send_data(ii)
                     q_recv(index+ivar,inode,kk) = recv_data(ii)
                  end do
               end do
            end do
            kk = kk + 1
         end if
         jj = jj + 1
      end do
   end do

end subroutine unpack_data_dg_general_lap_bcl

subroutine unpack_data_dg_general_consistency(G, b, par, q_send, q_recv, send_data, recv_data, nlayers)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nlayers

   real, dimension(nlayers+1,b%ngl,par%num_send_recv_total), intent(out) :: q_send, q_recv
   real, dimension((nlayers+1)*b%ngl*par%num_send_recv_total), intent(in)  :: send_data, recv_data

   integer :: ii, jj, kk, i, inbh, ib, ifaces, inode, jnode, ivar, ilocl, ilocr
   integer :: nq_i, nq_j, plane_ij, iface, imulti, ftype

   ii = 0
   jj = 1
   kk = 1

   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then
            do inode = 1, b%ngl
               do ivar = 1, nlayers+1
                  ii = ii + 1
                  q_send(ivar,inode,kk) = send_data(ii)
                  q_recv(ivar,inode,kk) = recv_data(ii)
               end do
            end do
            kk = kk + 1
         end if
         jj = jj + 1
      end do
   end do

end subroutine unpack_data_dg_general_consistency

subroutine unpack_data_dg_general_quad_layer(G, par, q_send, q_recv, send_data, recv_data, nvarb, nlayers, nq)

   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb, nlayers, nq

   real, dimension(nvarb,nq,par%num_send_recv_total,nlayers), intent(out) :: q_send, q_recv
   real, dimension(nvarb*nq*par%num_send_recv_total*nlayers), intent(in)  :: send_data, recv_data

   integer :: ii, jj, kk, i, inbh, ib, ifaces, inode, jnode, ivar, ilocl, ilocr
   integer :: nq_i, nq_j, plane_ij, iface, imulti, ftype, ll

   ii = 0
   jj = 1
   kk = 1

   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then
            do ll = 1, nlayers
               do inode = 1, nq
                  do ivar = 1, nvarb
                     ii = ii + 1
                     q_send(ivar,inode,kk,ll) = send_data(ii)
                     q_recv(ivar,inode,kk,ll) = recv_data(ii)
                  end do
               end do
            end do
            kk = kk + 1
         end if
         jj = jj + 1

      end do
   end do

end subroutine unpack_data_dg_general_quad_layer

subroutine unpack_data_dg_general_quad_1v(G, b, par, q_send, q_recv, send_data, recv_data)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par

   real, dimension(b%nq,par%num_send_recv_total), intent(out) :: q_send, q_recv
   real, dimension(b%nq*par%num_send_recv_total), intent(in)  :: send_data, recv_data

   integer :: ii, jj, kk, i, inbh, ib, ifaces, inode, jnode, ivar, ilocl, ilocr
   integer :: nq_i, nq_j, plane_ij, iface, imulti, ftype

   ii = 0
   jj = 1
   kk = 1

   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then
            do inode = 1, b%nq
               ii = ii + 1
               q_send(inode,kk) = send_data(ii)
               q_recv(inode,kk) = recv_data(ii)
            end do
            kk = kk + 1
         end if
         jj = jj + 1

      end do
   end do

end subroutine unpack_data_dg_general_quad_1v


subroutine pack_data_dg_quad(G, b, par, q_send, q_face, nvarb)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb

   real, intent(out) :: q_send(nvarb*b%nq*par%num_send_recv_total)
   real, intent(in)  :: q_face(nvarb,2,b%nq,G%nface)

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr
   integer :: iface_type

   ii = 0
   jj = 1
   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then

            ilocl = G%face(5,iface)
            el    = G%face(7,iface)

            do inode = 1, b%nq
               do ivar = 1, nvarb
                  ii = ii + 1
                  q_send(ii) = q_face(ivar,1,inode,iface)
               end do
            end do
         end if
         jj = jj + 1
      end do
   end do

end subroutine pack_data_dg_quad

subroutine pack_data_dg_df(G, b, mf, init, par, q_send, q_face, nvarb)

   use mod_basis,    only: basis
   use mod_face,     only: face_CS
   use mod_grid,     only: grid
   use mod_initial,  only: initial
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(face_CS),     intent(in) :: mf
   type(initial),     intent(in) :: init
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb

   real, intent(out) :: q_send(5*b%ngl*par%num_send_recv_total)
   real, intent(in)  :: q_face(nvarb,2,b%ngl,G%nface)

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr
   integer :: iface_type

   ii = 0
   jj = 1
   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then

            ilocl = G%face(5,iface)
            el    = G%face(7,iface)

            do inode = 1, b%ngl

               il = mf%imapl(1,inode,1,iface)
               jl = mf%imapl(2,inode,1,iface)
               kl = mf%imapl(3,inode,1,iface)

               ip = G%intma(il,jl,kl,el)

               do ivar = 1, nvarb
                  ii = ii + 1
                  q_send(ii) = q_face(ivar,1,inode,iface)
               end do
               ii = ii + 1
               q_send(ii) = init%pbprime_df(ip)
            end do
         end if
         jj = jj + 1
      end do
   end do

end subroutine pack_data_dg_df

subroutine pack_data_dg_df_btp(G, inp, b, mf, init, par, ref, q_send, q, qprime_df, nvarb)

   use mod_basis,    only: basis
   use mod_face,     only: face_CS
   use mod_grid,     only: grid
   use mod_initial,  only: initial
   use mod_parallel, only: parallel_CS
   use mod_input,    only: input
   use mod_ref,      only: mref

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(face_CS),     intent(in) :: mf
   type(initial),     intent(in) :: init
   type(input),       intent(in) :: inp
   type(parallel_CS), intent(in) :: par
   type(mref),        intent(in) :: ref
   integer,           intent(in) :: nvarb

   real, intent(out) :: q_send(ref%nbtp_var*b%ngl*par%num_send_recv_total)
   real, intent(in)  :: q(nvarb,G%npoin)
   real, intent(in)  :: qprime_df(3,G%npoin,inp%nlayers)

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr
   integer :: iface_type, ll

   ii = 0
   jj = 1
   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then

            ilocl = G%face(5,iface)
            el    = G%face(7,iface)

            do inode = 1, b%ngl

               il = mf%imapl(1,inode,1,iface)
               jl = mf%imapl(2,inode,1,iface)
               kl = mf%imapl(3,inode,1,iface)

               ip = G%intma(il,jl,kl,el)

               do ivar = 1, nvarb
                  ii = ii + 1
                  q_send(ii) = q(ivar,ip)
               end do
               ii = ii + 1
               q_send(ii) = init%pbprime_df(ip)

               do ll = 1, inp%nlayers
                  do ivar = 1, 3
                     ii = ii + 1
                     q_send(ii) = qprime_df(ivar,ip,ll)
                  end do
               end do
            end do
         end if
         jj = jj + 1
      end do
   end do

end subroutine pack_data_dg_df_btp

subroutine pack_data_dg_df_btp_lap(G, b, mf, init, par, q_send, q, btp_dpp_graduv, pbprime_visc, nvarb)

   use mod_basis,    only: basis
   use mod_face,     only: face_CS
   use mod_grid,     only: grid
   use mod_initial,  only: initial
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(face_CS),     intent(in) :: mf
   type(initial),     intent(in) :: init
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb

   real, intent(out) :: q_send(10*b%ngl*par%num_send_recv_total)
   real, intent(in)  :: q(nvarb,G%npoin)
   real, intent(in)  :: btp_dpp_graduv(4,G%npoin)
   real, intent(in)  :: pbprime_visc(G%npoin)

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr
   integer :: iface_type

   ii = 0
   jj = 1
   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then

            ilocl = G%face(5,iface)
            el    = G%face(7,iface)

            do inode = 1, b%ngl

               il = mf%imapl(1,inode,1,iface)
               jl = mf%imapl(2,inode,1,iface)
               kl = mf%imapl(3,inode,1,iface)

               ip = G%intma(il,jl,kl,el)

               do ivar = 1, nvarb
                  ii = ii + 1
                  q_send(ii) = q(ivar,ip)
               end do
               ii = ii + 1
               q_send(ii) = init%pbprime_df(ip)
               do ivar = 1, 4
                  ii = ii + 1
                  q_send(ii) = btp_dpp_graduv(ivar,ip)
               end do
               ii = ii + 1
               q_send(ii) = pbprime_visc(ip)

            end do
         end if
         jj = jj + 1
      end do
   end do

end subroutine pack_data_dg_df_btp_lap

subroutine pack_data_dg_df_bcl(G, inp, b, mf, par, q_send, q)

   use mod_basis,    only: basis
   use mod_face,     only: face_CS
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS
   use mod_input,    only: input

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(face_CS),     intent(in) :: mf
   type(input),       intent(in) :: inp
   type(parallel_CS), intent(in) :: par

   real, intent(out) :: q_send(3*inp%nlayers*b%ngl*par%num_send_recv_total)
   real, intent(in)  :: q(3,G%npoin,inp%nlayers)

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr, ll
   integer :: iface_type

   ii = 0
   jj = 1
   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then

            ilocl = G%face(5,iface)
            el    = G%face(7,iface)

            do ll = 1, inp%nlayers
               do inode = 1, b%ngl

                  il = mf%imapl(1,inode,1,iface)
                  jl = mf%imapl(2,inode,1,iface)
                  kl = mf%imapl(3,inode,1,iface)

                  ip = G%intma(il,jl,kl,el)
                  do ivar = 1, 3
                     ii = ii + 1
                     q_send(ii) = q(ivar,ip,ll)
                  end do
               end do
            end do
         end if
         jj = jj + 1
      end do
   end do

end subroutine pack_data_dg_df_bcl

subroutine pack_data_dg_df_bcl_lap(G, b, mf, par, q_send, dpp_graduv, dpprime_visc, nlayers)

   use mod_basis,    only: basis
   use mod_face,     only: face_CS
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(face_CS),     intent(in) :: mf
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nlayers

   real, intent(out) :: q_send(5*nlayers*b%ngl*par%num_send_recv_total)
   real, intent(in)  :: dpp_graduv(4,G%npoin,nlayers)
   real, intent(in)  :: dpprime_visc(G%npoin,nlayers)

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr
   integer :: iface_type, ll

   ii = 0
   jj = 1
   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then

            ilocl = G%face(5,iface)
            el    = G%face(7,iface)

            do ll = 1, nlayers

               do inode = 1, b%ngl

                  il = mf%imapl(1,inode,1,iface)
                  jl = mf%imapl(2,inode,1,iface)
                  kl = mf%imapl(3,inode,1,iface)

                  ip = G%intma(il,jl,kl,el)

                  do ivar = 1, 4
                     ii = ii + 1
                     q_send(ii) = dpp_graduv(ivar,ip,ll)
                  end do
                  ii = ii + 1
                  q_send(ii) = dpprime_visc(ip,ll)
               end do

            end do
         end if
         jj = jj + 1
      end do
   end do

end subroutine pack_data_dg_df_bcl_lap

subroutine pack_data_dg_consistency(G, b, mf, init, par, q_send, dprime_df, nlayers)

   use mod_basis,    only: basis
   use mod_face,     only: face_CS
   use mod_grid,     only: grid
   use mod_initial,  only: initial
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(face_CS),     intent(in) :: mf
   type(initial),     intent(in) :: init
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nlayers

   real, intent(out) :: q_send((nlayers+1)*b%ngl*par%num_send_recv_total)
   real, intent(in)  :: dprime_df(G%npoin,nlayers)

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr
   integer :: iface_type, ll

   ii = 0
   jj = 1
   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then

            ilocl = G%face(5,iface)
            el    = G%face(7,iface)

            do inode = 1, b%ngl

               il = mf%imapl(1,inode,1,iface)
               jl = mf%imapl(2,inode,1,iface)
               kl = mf%imapl(3,inode,1,iface)
               ip = G%intma(il,jl,kl,el)

               do ll = 1, nlayers
                  ii = ii + 1
                  q_send(ii) = dprime_df(ip,ll)
               end do
               ii = ii + 1
               q_send(ii) = init%pbprime_df(ip)
            end do
         end if
         jj = jj + 1
      end do
   end do

end subroutine pack_data_dg_consistency

subroutine pack_data_dg_quad_all(G, b, par, q_send, q_face, grad_face, nvarb)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb

   real, intent(out) :: q_send(2*nvarb*b%nq*par%num_send_recv_total)
   real, intent(in)  :: q_face(nvarb,2,b%nq,G%nface)
   real, intent(in)  :: grad_face(nvarb,2,b%nq,G%nface)

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr
   integer :: iface_type

   ii = 0
   jj = 1
   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then

            ilocl = G%face(5,iface)
            el    = G%face(7,iface)

            do inode = 1, b%nq

               do ivar = 1, nvarb
                  ii = ii + 1
                  q_send(ii) = q_face(ivar,1,inode,iface)
               end do

               do ivar = 1, nvarb
                  ii = ii + 1
                  q_send(ii) = grad_face(ivar,1,inode,iface)
               end do
            end do
         end if
         jj = jj + 1
      end do
   end do

end subroutine pack_data_dg_quad_all

subroutine pack_data_dg_quad_lap(G, b, mf, par, q_send, grad_uvdp, nvarb)

   use mod_basis,    only: basis
   use mod_face,     only: face_CS
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(face_CS),     intent(in) :: mf
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb

   real, intent(out) :: q_send(nvarb*b%nq*par%num_send_recv_total)
   real, intent(in)  :: grad_uvdp(nvarb,G%npoin_q)

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr
   integer :: iface_type

   ii = 0
   jj = 1
   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         !if(face_type(iface) == 2 .and. imulti > 0) then

         ilocl = G%face(5,iface)
         el    = G%face(7,iface)

         do inode = 1, b%nq

            il = mf%imapr_q(1,inode,1,iface)
            jl = mf%imapr_q(2,inode,1,iface)
            kl = mf%imapr_q(3,inode,1,iface)

            ip = G%intma_dg_quad(il,jl,kl,el)

            do ivar = 1, nvarb
               ii = ii + 1
               q_send(ii) = grad_uvdp(ivar,ip)
            end do
         end do
         !end if
         jj = jj + 1
      end do
   end do

end subroutine pack_data_dg_quad_lap

subroutine pack_data_dg_quad_layer(G, par, q_send, q_face, nvarb, nlayers, nq)

   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb, nlayers, nq

   real, intent(out) :: q_send(nvarb*nq*par%num_send_recv_total*nlayers)
   real, intent(in)  :: q_face(nvarb,2,nq,G%nface,nlayers)

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij, ll
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr
   integer :: iface_type

   ii = 0
   jj = 1
   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then

            ilocl = G%face(5,iface)
            el    = G%face(7,iface)

            do ll = 1, nlayers
               do inode = 1, nq
                  do ivar = 1, nvarb
                     ii = ii + 1
                     q_send(ii) = q_face(ivar,1,inode,iface,ll)
                  end do
               end do
            end do
         end if
         jj = jj + 1
      end do
   end do

end subroutine pack_data_dg_quad_layer

subroutine pack_data_dg_quad_layer_all(G, b, par, q_send, q_face, qprime_face, nvarb, nlayers)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nvarb, nlayers

   real, intent(out) :: q_send(2*nvarb*b%nq*par%num_send_recv_total*nlayers)
   real, dimension(nvarb,2,b%nq,G%nface,nlayers), intent(in) :: q_face, qprime_face

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij, ll
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr
   integer :: iface_type

   ii = 0
   jj = 1
   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         if (G%face_type(iface) == 2 .and. imulti > 0) then

            ilocl = G%face(5,iface)
            el    = G%face(7,iface)

            do ll = 1, nlayers
               do inode = 1, b%nq

                  do ivar = 1, nvarb
                     ii = ii + 1
                     q_send(ii) = q_face(ivar,1,inode,iface,ll)
                  end do

                  do ivar = 1, nvarb
                     ii = ii + 1
                     q_send(ii) = qprime_face(ivar,1,inode,iface,ll)
                  end do

               end do
            end do
         end if
         jj = jj + 1
      end do
   end do

end subroutine pack_data_dg_quad_layer_all

subroutine pack_data_dg_quad_1v(G, b, par, q_send, q_face)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par

   real, intent(out) :: q_send(b%nq*par%num_send_recv_total)
   real, intent(in)  :: q_face(2,b%nq,G%nface)

   integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
   integer :: nq_i, nq_j, plane_ij
   real :: h, qu, qv
   integer :: inode, jnode, ip, ilocl, ilocr
   integer :: iface_type

   ii = 0
   jj = 1

   do inbh = 1, par%num_nbh
      do ib = 1, par%num_send_recv(inbh)
         iface  = par%nbh_send_recv(jj)
         imulti = par%nbh_send_recv_multi(jj)

         !if(face_type(iface) == 2 .and. imulti > 0) then

         ilocl = G%face(5,iface)
         el    = G%face(7,iface)

         do inode = 1, b%nq
            ii = ii + 1
            q_send(ii) = q_face(1,inode,iface)
         end do
         !end if

         jj = jj + 1

      end do
   end do

end subroutine pack_data_dg_quad_1v

!----------------------------------------------------------------------!
!>@brief Send boundary data
!>@author James F. Kelly
!>        Department of Applied Mathematics
!>        Naval Postgraduate School
!>        Monterey, CA 93943
!>@date Feb 2015, Daniel S. Abdi, Extended data size for total pressure
!>@date October 9, 2015, F.X. Giraldo, modified to handle LHS operators
!>@date January 7, 2016, M.A. Kopera, modified to handle AMR
!>@date December 25, 2016, F.X. Giraldo, modified to be a general DG communicator of size NSIZE
!----------------------------------------------------------------------!
subroutine send_bound_dg_general(G, b, par, send_data, recv_data, nsize, nreq, ireq, status)

   use mod_basis,    only: basis
   use mod_grid,     only: grid, mod_grid_get_face_ngl
   use mod_parallel, only: parallel_CS
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(parallel_CS), intent(in)  :: par
   integer,           intent(in)  :: nsize

   real, intent(in)  :: send_data(nsize*b%ngl*b%ngl*par%num_send_recv_total)
   real, intent(out) :: recv_data(nsize*b%ngl*b%ngl*par%num_send_recv_total)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size,2*par%num_nbh)

   integer :: inbh, idest, istart, iend, ierr, nq
   integer :: ngl_i, ngl_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

   nreq = 0
   iend = 0
   jj = 1
   status = 0

   do inbh = 1, par%num_nbh

      nq = 0
      do ib = 1, par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)

         ilocl = G%face(5,iface)

         do i = 1, par%nbh_send_recv_multi(jj)
            call mod_grid_get_face_ngl(b, ilocl, ngl_i, ngl_j, plane_ij)
            nq = nq + ngl_i * ngl_j * nsize
         end do

         jj = jj + 1

      end do

      idest = par%nbh_proc(inbh)

      nreq = nreq + 1
      istart = iend + 1
      iend = istart + nq - 1

      if (nq > 0) then
         call mpi_irecv(recv_data(istart:iend), nq, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq), ierr)

         call mpi_isend(send_data(istart:iend), nq, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      endif

      nreq = nreq + 1
   end do

end subroutine send_bound_dg_general

subroutine send_bound_dg_general_quad(G, b, par, send_data, recv_data, nvarb, nreq, ireq, status)

   use mod_basis,    only: basis
   use mod_grid,     only: grid, mod_grid_get_face_nq
   use mod_parallel, only: parallel_CS
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(parallel_CS), intent(in)  :: par
   integer,           intent(in)  :: nvarb

   real, intent(in)  :: send_data(nvarb*b%nq*par%num_send_recv_total)
   real, intent(out) :: recv_data(nvarb*b%nq*par%num_send_recv_total)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size,2*par%num_nbh)

   integer :: inbh, idest, istart, iend, ierr, nqp
   integer :: nq_i, nq_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

   nreq = 0
   iend = 0
   jj = 1
   status = 0

   do inbh = 1, par%num_nbh

      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)

         ilocl = G%face(5,iface)

         do i = 1, par%nbh_send_recv_multi(jj)
            call mod_grid_get_face_nq(b, ilocl, nq_i, nq_j, plane_ij)
            nqp = nqp + b%nq * nvarb
         end do
         jj = jj + 1
      end do

      idest = par%nbh_proc(inbh)

      nreq = nreq + 1
      istart = iend + 1
      iend = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq), ierr)

         call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      endif

      nreq = nreq + 1
   end do

end subroutine send_bound_dg_general_quad

subroutine send_bound_dg_general_df(G, b, par, send_data, recv_data, nvarb, nreq, ireq, status)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(parallel_CS), intent(in)  :: par
   integer,           intent(in)  :: nvarb

   real, intent(in)  :: send_data(nvarb*b%ngl*par%num_send_recv_total)
   real, intent(out) :: recv_data(nvarb*b%ngl*par%num_send_recv_total)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size,2*par%num_nbh)

   integer :: inbh, idest, istart, iend, ierr, nqp
   integer :: nq_i, nq_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

   nreq = 0
   iend = 0
   jj = 1
   status = 0

   do inbh = 1, par%num_nbh

      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)

         ilocl = G%face(5,iface)

         do i = 1, par%nbh_send_recv_multi(jj)
            nqp = nqp + b%ngl * nvarb
         end do
         jj = jj + 1
      end do

      idest = par%nbh_proc(inbh)

      nreq = nreq + 1
      istart = iend + 1
      iend = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq), ierr)

         call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      endif

      nreq = nreq + 1
   end do

end subroutine send_bound_dg_general_df

subroutine pack_and_send_df_btp(G, inp, b, mf, init, par, ref, send_data_dg, recv_data_dg, q, qprime_df, nvarb, nreq, ireq, status)

   use mod_basis,         only: basis
   use mod_face,          only: face_CS
   use mod_grid,          only: grid
   use mod_initial,       only: initial
   use mod_parallel,      only: parallel_CS
   use mod_input,         only: input
   use mod_ref,           only: mref
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(face_CS),     intent(in)  :: mf
   type(initial),     intent(in)  :: init
   type(input),       intent(in)  :: inp
   type(parallel_CS), intent(in)  :: par
   type(mref),        intent(in)  :: ref
   integer,           intent(in)  :: nvarb
   real,              intent(out) :: send_data_dg(ref%nbtp_var*b%ngl*G%nboun)
   real,              intent(out) :: recv_data_dg(ref%nbtp_var*b%ngl*G%nboun)
   real,              intent(in)  :: q(nvarb, G%npoin)
   real,              intent(in)  :: qprime_df(3, G%npoin, inp%nlayers)
   integer,           intent(out) :: nreq
   integer,           intent(out) :: ireq(2*par%num_nbh)
   integer,           intent(out) :: status(mpi_status_size, 2*par%num_nbh)

   integer :: jj, kk, i, inbh, ib, iface, el, ivar
   integer :: inode, ip, ll, ioff
   integer :: nqp, istart, iend, idest, ierr
   integer :: ngl_f, nbtp_var_f, nlayers_f, nboun_valid

   ngl_f        = b%ngl
   nbtp_var_f   = ref%nbtp_var
   nlayers_f    = inp%nlayers
   nboun_valid  = ref%nboun_valid

   !$acc parallel loop gang private(iface, el) &
   !$acc    present(G%face, G%intma, mf%imapl, init%pbprime_df, q, qprime_df, &
   !$acc            send_data_dg, ref%face_pack_list) &
   !$acc    firstprivate(nvarb, ngl_f, nbtp_var_f, nlayers_f, nboun_valid)
   do kk = 1, nboun_valid
      iface = ref%face_pack_list(kk)
      el    = G%face(7, iface)
      !$acc loop vector private(ip, ioff)
      do inode = 1, ngl_f
         ip   = G%intma(mf%imapl(1,inode,1,iface), &
            mf%imapl(2,inode,1,iface), &
            mf%imapl(3,inode,1,iface), el)
         ioff = (kk-1)*ngl_f*nbtp_var_f + (inode-1)*nbtp_var_f
         !$acc loop seq
         do ivar = 1, nvarb
            send_data_dg(ioff+ivar) = q(ivar, ip)
         end do
         send_data_dg(ioff+nvarb+1) = init%pbprime_df(ip)
         !$acc loop seq
         do ll = 1, nlayers_f
            !$acc loop seq
            do ivar = 1, 3
               send_data_dg(ioff+nvarb+1+(ll-1)*3+ivar) = qprime_df(ivar, ip, ll)
            end do
         end do
      end do
   end do
   !$acc end parallel loop

   ! GPU-direct non-blocking send/recv per neighbor.
   nreq   = 0
   iend   = 0
   jj     = 1
   status = 0

   !$acc host_data use_device(send_data_dg, recv_data_dg)
   do inbh = 1, par%num_nbh
      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         do i = 1, par%nbh_send_recv_multi(jj)
            nqp = nqp + b%ngl * ref%nbtp_var
         end do
         jj = jj + 1
      end do

      idest  = par%nbh_proc(inbh)
      nreq   = nreq + 1
      istart = iend + 1
      iend   = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(recv_data_dg(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, ireq(nreq), ierr)
         call mpi_isend(send_data_dg(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      end if
      nreq = nreq + 1
   end do
   !$acc end host_data

end subroutine pack_and_send_df_btp

subroutine pack_and_send_df_btp_lap(G, b, mf, init, par, ref,              &
   q, btp_dpp_graduv, pbprime_visc, nvarb, &
   nreq, ireq, status)

   use mod_basis,         only: basis
   use mod_face,          only: face_CS
   use mod_grid,          only: grid
   use mod_initial,       only: initial
   use mod_parallel,      only: parallel_CS
   use mod_ref,           only: mref
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)    :: G
   type(basis),       intent(in)    :: b
   type(face_CS),     intent(in)    :: mf
   type(initial),     intent(in)    :: init
   type(parallel_CS), intent(in)    :: par
   type(mref),        intent(inout) :: ref
   integer,           intent(in)    :: nvarb

   real,    intent(in)  :: q(nvarb, G%npoin)
   real,    intent(in)  :: btp_dpp_graduv(4, G%npoin)
   real,    intent(in)  :: pbprime_visc(G%npoin)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size, 2*par%num_nbh)

   integer :: kk, inode, iface, el, ip, ioff, ivar
   integer :: jj, ib, inbh, nqp, istart, iend, idest, ierr, i
   integer :: ngl_f, nboun_valid_f

   ngl_f         = b%ngl
   nboun_valid_f = ref%nboun_valid

   ! GPU packing: one gang per boundary face, one vector lane per node.
   !$acc parallel loop gang private(iface, el)                                      &
   !$acc    present(G%face, G%intma, mf%imapl,                                     &
   !$acc            init%pbprime_df, q, btp_dpp_graduv, pbprime_visc,              &
   !$acc            ref%send_data_dg_lap, ref%face_pack_list)                       &
   !$acc    firstprivate(nvarb, ngl_f, nboun_valid_f)
   do kk = 1, nboun_valid_f
      iface = ref%face_pack_list(kk)
      el    = G%face(7, iface)
      !$acc loop vector private(ip, ioff)
      do inode = 1, ngl_f
         ip   = G%intma(mf%imapl(1,inode,1,iface), &
            mf%imapl(2,inode,1,iface), &
            mf%imapl(3,inode,1,iface), el)
         ioff = (kk-1)*ngl_f*10 + (inode-1)*10
         !$acc loop seq
         do ivar = 1, nvarb
            ref%send_data_dg_lap(ioff+ivar) = q(ivar, ip)
         end do
         ref%send_data_dg_lap(ioff+nvarb+1) = init%pbprime_df(ip)
         !$acc loop seq
         do ivar = 1, 4
            ref%send_data_dg_lap(ioff+nvarb+1+ivar) = btp_dpp_graduv(ivar, ip)
         end do
         ref%send_data_dg_lap(ioff+nvarb+6) = pbprime_visc(ip)
      end do
   end do
   !$acc end parallel loop

   ! GPU-direct non-blocking send/recv per neighbor.
   nreq   = 0
   iend   = 0
   jj     = 1
   status = 0

   !$acc host_data use_device(ref%send_data_dg_lap, ref%recv_data_dg_lap)
   do inbh = 1, par%num_nbh
      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         do i = 1, par%nbh_send_recv_multi(jj)
            nqp = nqp + b%ngl * 10
         end do
         jj = jj + 1
      end do

      idest  = par%nbh_proc(inbh)
      nreq   = nreq + 1
      istart = iend + 1
      iend   = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(ref%recv_data_dg_lap(istart:iend), nqp, MPI_PRECISION, &
            idest-1, 99, mpi_comm_world, ireq(nreq),   ierr)
         call mpi_isend(ref%send_data_dg_lap(istart:iend), nqp, MPI_PRECISION, &
            idest-1, 99, mpi_comm_world, ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      end if
      nreq = nreq + 1
   end do
   !$acc end host_data

end subroutine pack_and_send_df_btp_lap

subroutine pack_and_send_df_btp_v0(G, inp, b, mf, init, par, ref, send_data_dg, recv_data_dg, q, qprime_df, nvarb, nreq, ireq, status)

   use mod_basis,         only: basis
   use mod_face,          only: face_CS
   use mod_grid,          only: grid
   use mod_initial,       only: initial
   use mod_parallel,      only: parallel_CS
   use mod_input,         only: input
   use mod_ref,           only: mref
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(face_CS),     intent(in)  :: mf
   type(initial),     intent(in)  :: init
   type(input),       intent(in)  :: inp
   type(parallel_CS), intent(in)  :: par
   type(mref),        intent(in)  :: ref
   integer,           intent(in)  :: nvarb
   real,              intent(out) :: send_data_dg(ref%nbtp_var*b%ngl*G%nboun)
   real,              intent(out) :: recv_data_dg(ref%nbtp_var*b%ngl*G%nboun)
   real,              intent(in)  :: q(nvarb, G%npoin)
   real,              intent(in)  :: qprime_df(3, G%npoin, inp%nlayers)
   integer,           intent(out) :: nreq
   integer,           intent(out) :: ireq(2*par%num_nbh)
   integer,           intent(out) :: status(mpi_status_size, 2*par%num_nbh)

   integer :: jj, inbh, ib, iface, imulti, el, ivar, i
   integer :: inode, ip, ll, ioff
   integer :: nqp, istart, iend, idest, ierr
   integer :: ngl_f, nbtp_var_f, nlayers_f, nboun_f

   ngl_f      = b%ngl
   nbtp_var_f = ref%nbtp_var
   nlayers_f  = inp%nlayers
   nboun_f    = G%nboun

   ! GPU pack — one gang per MPI boundary face.
   !$acc parallel loop gang private(iface, imulti, el) &
   !$acc    present(G%face_type, G%face, G%intma, par%nbh_send_recv, &
   !$acc            par%nbh_send_recv_multi, mf%imapl, init%pbprime_df, &
   !$acc            q, qprime_df, send_data_dg) &
   !$acc    firstprivate(nvarb, ngl_f, nbtp_var_f, nlayers_f)
   do jj = 1, nboun_f
      iface  = par%nbh_send_recv(jj)
      imulti = par%nbh_send_recv_multi(jj)
      if (G%face_type(iface) == 2 .and. imulti > 0) then
         el = G%face(7, iface)
         !$acc loop seq
         do inode = 1, ngl_f
            ip   = G%intma(mf%imapl(1,inode,1,iface), &
               mf%imapl(2,inode,1,iface), &
               mf%imapl(3,inode,1,iface), el)
            ioff = (jj-1)*ngl_f*nbtp_var_f + (inode-1)*nbtp_var_f
            !$acc loop seq
            do ivar = 1, nvarb
               send_data_dg(ioff+ivar) = q(ivar, ip)
            end do
            send_data_dg(ioff+nvarb+1) = init%pbprime_df(ip)
            !$acc loop seq
            do ll = 1, nlayers_f
               !$acc loop seq
               do ivar = 1, 3
                  send_data_dg(ioff+nvarb+1+(ll-1)*3+ivar) = qprime_df(ivar, ip, ll)
               end do
            end do
         end do
      end if
   end do
   !$acc end parallel loop

   ! GPU-direct non-blocking send/recv per neighbor.
   nreq  = 0
   iend  = 0
   jj    = 1
   status = 0

   !$acc host_data use_device(send_data_dg, recv_data_dg)
   do inbh = 1, par%num_nbh
      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         do i = 1, par%nbh_send_recv_multi(jj)
            nqp = nqp + b%ngl * ref%nbtp_var
         end do
         jj = jj + 1
      end do

      idest  = par%nbh_proc(inbh)
      nreq   = nreq + 1
      istart = iend + 1
      iend   = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(recv_data_dg(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, ireq(nreq), ierr)
         call mpi_isend(send_data_dg(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      end if
      nreq = nreq + 1
   end do
   !$acc end host_data

end subroutine pack_and_send_df_btp_v0

subroutine send_bound_dg_general_lap(G, b, par, send_data, recv_data, nvarb, nreq, ireq, status)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(parallel_CS), intent(in)  :: par
   integer,           intent(in)  :: nvarb

   real, intent(in)  :: send_data(10*b%ngl*par%num_send_recv_total)
   real, intent(out) :: recv_data(10*b%ngl*par%num_send_recv_total)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size,2*par%num_nbh)

   integer :: inbh, idest, istart, iend, ierr, nqp
   integer :: nq_i, nq_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

   nreq = 0
   iend = 0
   jj = 1
   status = 0

   do inbh = 1, par%num_nbh

      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)

         ilocl = G%face(5,iface)

         do i = 1, par%nbh_send_recv_multi(jj)
            nqp = nqp + b%ngl * 10
         end do
         jj = jj + 1
      end do

      idest = par%nbh_proc(inbh)

      nreq = nreq + 1
      istart = iend + 1
      iend = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq), ierr)

         call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      endif

      nreq = nreq + 1
   end do

end subroutine send_bound_dg_general_lap

subroutine send_bound_dg_general_bcl(G, b, inp, par, send_data, recv_data, nreq, ireq, status)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION
   use mod_input,    only: input

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(input),       intent(in)  :: inp
   type(parallel_CS), intent(in)  :: par

   real, intent(in)  :: send_data(3*inp%nlayers*b%ngl*par%num_send_recv_total)
   real, intent(out) :: recv_data(3*inp%nlayers*b%ngl*par%num_send_recv_total)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size,2*par%num_nbh)

   integer :: inbh, idest, istart, iend, ierr, nqp
   integer :: nq_i, nq_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

   nreq = 0
   iend = 0
   jj = 1
   status = 0

   do inbh = 1, par%num_nbh

      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)

         ilocl = G%face(5,iface)

         do i = 1, par%nbh_send_recv_multi(jj)
            nqp = nqp + b%ngl * (3*inp%nlayers)
         end do
         jj = jj + 1
      end do

      idest = par%nbh_proc(inbh)

      nreq = nreq + 1
      istart = iend + 1
      iend = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq), ierr)

         call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      endif

      nreq = nreq + 1
   end do

end subroutine send_bound_dg_general_bcl

subroutine send_bound_dg_general_lap_bcl(G, b, par, send_data, recv_data, nlayers, nreq, ireq, status)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(parallel_CS), intent(in)  :: par
   integer,           intent(in)  :: nlayers

   real, intent(in)  :: send_data(5*nlayers*b%ngl*par%num_send_recv_total)
   real, intent(out) :: recv_data(5*nlayers*b%ngl*par%num_send_recv_total)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size,2*par%num_nbh)

   integer :: inbh, idest, istart, iend, ierr, nqp
   integer :: nq_i, nq_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

   nreq = 0
   iend = 0
   jj = 1
   status = 0

   do inbh = 1, par%num_nbh

      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)

         ilocl = G%face(5,iface)

         do i = 1, par%nbh_send_recv_multi(jj)
            nqp = nqp + b%ngl * 5*nlayers
         end do
         jj = jj + 1
      end do

      idest = par%nbh_proc(inbh)

      nreq = nreq + 1
      istart = iend + 1
      iend = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq), ierr)

         call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      endif

      nreq = nreq + 1
   end do

end subroutine send_bound_dg_general_lap_bcl

subroutine send_bound_dg_general_consistency(G, b, par, send_data, recv_data, nlayers, nreq, ireq, status)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(parallel_CS), intent(in)  :: par
   integer,           intent(in)  :: nlayers

   real, intent(in)  :: send_data((nlayers+1)*b%ngl*par%num_send_recv_total)
   real, intent(out) :: recv_data((nlayers+1)*b%ngl*par%num_send_recv_total)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size,2*par%num_nbh)

   integer :: inbh, idest, istart, iend, ierr, nqp
   integer :: nq_i, nq_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

   nreq = 0
   iend = 0
   jj = 1
   status = 0

   do inbh = 1, par%num_nbh

      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)

         ilocl = G%face(5,iface)

         do i = 1, par%nbh_send_recv_multi(jj)
            nqp = nqp + b%ngl * (nlayers+1)
         end do
         jj = jj + 1
      end do

      idest = par%nbh_proc(inbh)

      nreq = nreq + 1
      istart = iend + 1
      iend = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq), ierr)

         call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      endif

      nreq = nreq + 1
   end do

end subroutine send_bound_dg_general_consistency

subroutine send_bound_dg_general_quad_layer(G, b, par, send_data, recv_data, nvarb, nlayers, nq, nreq, ireq, status)

   use mod_basis,    only: basis
   use mod_grid,     only: grid, mod_grid_get_face_nq
   use mod_parallel, only: parallel_CS
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(parallel_CS), intent(in)  :: par
   integer,           intent(in)  :: nvarb, nlayers, nq

   real, intent(in)  :: send_data(nvarb*nq*par%num_send_recv_total*nlayers)
   real, intent(out) :: recv_data(nvarb*nq*par%num_send_recv_total*nlayers)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size,2*par%num_nbh)

   integer :: inbh, idest, istart, iend, ierr, nqp
   integer :: nq_i, nq_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

   nreq = 0
   iend = 0
   jj = 1
   status = 0

   do inbh = 1, par%num_nbh

      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         iface = par%nbh_send_recv(jj)

         ilocl = G%face(5,iface)

         do i = 1, par%nbh_send_recv_multi(jj)
            call mod_grid_get_face_nq(b, ilocl, nq_i, nq_j, plane_ij)
            nqp = nqp + nq * nvarb*nlayers
         end do

         jj = jj + 1

      end do

      idest = par%nbh_proc(inbh)

      nreq = nreq + 1
      istart = iend + 1
      iend = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq), ierr)

         call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, &
            ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      endif

      nreq = nreq + 1
   end do

end subroutine send_bound_dg_general_quad_layer
