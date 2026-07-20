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

subroutine unpack_data_dg_general_lap(G, b, par, ref, q_send, q_recv, send_data, recv_data, nvarb)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS
   use mod_ref,      only: mref

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   type(mref),        intent(in) :: ref
   integer,           intent(in) :: nvarb

   real, dimension(ref%nbtp_var_lap,b%ngl,G%nboun), intent(out) :: q_send, q_recv
   real, dimension(ref%nbtp_var_lap*b%ngl*G%nboun), intent(in)  :: send_data, recv_data

   integer :: ii, kk, inode, ivar, ngl_f, nboun_f, nvarb_lap_f

   ngl_f       = b%ngl
   nboun_f     = G%nboun
   nvarb_lap_f = ref%nbtp_var_lap

   !$acc parallel loop gang vector collapse(3) &
   !$acc    present(q_send, q_recv, send_data, recv_data) &
   !$acc    firstprivate(ngl_f, nboun_f, nvarb_lap_f)
   do kk = 1, nboun_f
      do inode = 1, ngl_f
         do ivar = 1, nvarb_lap_f
            ii = (kk-1)*ngl_f*nvarb_lap_f + (inode-1)*nvarb_lap_f + ivar
            q_send(ivar, inode, kk) = send_data(ii)
            q_recv(ivar, inode, kk) = recv_data(ii)
         end do
      end do
   end do
   !$acc end parallel loop

end subroutine unpack_data_dg_general_lap

subroutine unpack_data_dg_general_bcl(G, b, inp, par, q_send, q_recv, send_data, recv_data, nboun_valid)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS
   use mod_input,    only: input

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(input),       intent(in) :: inp
   type(parallel_CS), intent(in) :: par
   integer,           intent(in) :: nboun_valid

   real, dimension(inp%nvar_bcl*inp%nlayers,b%ngl,par%num_send_recv_total), intent(out) :: q_send, q_recv
   real, dimension(inp%nvar_bcl*inp%nlayers*b%ngl*par%num_send_recv_total), intent(in)  :: send_data, recv_data

   integer :: kk, inode, ivar, ll
   integer :: ngl_f, nlayers_f, nvarb_f
   integer :: ii

   ngl_f     = b%ngl
   nlayers_f = inp%nlayers
   nvarb_f   = inp%nvar_bcl

   ! Each gang unpacks one boundary face using the same flat-buffer index
   ! formula as pack_and_send_df_bcl: (kk-1)*nlayers*ngl*nvarb + (ll-1)*ngl*nvarb + (inode-1)*nvarb + ivar
   !$acc parallel loop gang &
   !$acc    present(send_data, recv_data, q_send, q_recv) &
   !$acc    firstprivate(ngl_f, nlayers_f, nvarb_f, nboun_valid)
   do kk = 1, nboun_valid
      !$acc loop seq
      do ll = 1, nlayers_f
         !$acc loop vector private(ii)
         do inode = 1, ngl_f
            !$acc loop seq
            do ivar = 1, nvarb_f
               ii = (kk-1)*nlayers_f*ngl_f*nvarb_f + (ll-1)*ngl_f*nvarb_f + (inode-1)*nvarb_f + ivar
               q_send((ll-1)*nvarb_f+ivar, inode, kk) = send_data(ii)
               q_recv((ll-1)*nvarb_f+ivar, inode, kk) = recv_data(ii)
            end do
         end do
      end do
   end do
   !$acc end parallel loop

end subroutine unpack_data_dg_general_bcl

subroutine unpack_data_dg_general_lap_bcl(G, b, par, ref, q_send, q_recv, send_data, recv_data, nlayers, nboun_valid)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS
   use mod_ref,      only: mref

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   type(mref),        intent(in) :: ref
   integer,           intent(in) :: nlayers, nboun_valid

   real, dimension(ref%nbcl_var_lap,b%ngl,par%num_send_recv_total), intent(out) :: q_send, q_recv
   real, dimension(ref%nbcl_var_lap*b%ngl*par%num_send_recv_total), intent(in)  :: send_data, recv_data

   integer :: kk, ll, inode, ivar
   integer :: ngl_f, nlayers_f, nboun_valid_f, nvarb_lap_per_layer

   ngl_f               = b%ngl
   nlayers_f           = nlayers
   nboun_valid_f       = nboun_valid
   nvarb_lap_per_layer = ref%nbcl_var_lap / nlayers

   !$acc parallel loop gang &
   !$acc    present(send_data, recv_data, q_send, q_recv) &
   !$acc    firstprivate(ngl_f, nlayers_f, nboun_valid_f, nvarb_lap_per_layer)
   do kk = 1, nboun_valid_f
      !$acc loop seq
      do ll = 1, nlayers_f
         !$acc loop vector
         do inode = 1, ngl_f
            !$acc loop seq
            do ivar = 1, nvarb_lap_per_layer
               q_send((ll-1)*nvarb_lap_per_layer+ivar, inode, kk) = send_data((kk-1)*nlayers_f*ngl_f*nvarb_lap_per_layer + (ll-1)*ngl_f*nvarb_lap_per_layer + (inode-1)*nvarb_lap_per_layer + ivar)
               q_recv((ll-1)*nvarb_lap_per_layer+ivar, inode, kk) = recv_data((kk-1)*nlayers_f*ngl_f*nvarb_lap_per_layer + (ll-1)*ngl_f*nvarb_lap_per_layer + (inode-1)*nvarb_lap_per_layer + ivar)
            end do
         end do
      end do
   end do
   !$acc end parallel loop

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
   real, intent(in)  :: qprime_df(inp%nvar_bcl,G%npoin,inp%nlayers)

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
                  do ivar = 1, inp%nvar_bcl
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

   real, intent(out) :: q_send(inp%nvar_bcl*inp%nlayers*b%ngl*par%num_send_recv_total)
   real, intent(in)  :: q(inp%nvar_bcl,G%npoin,inp%nlayers)

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
                  do ivar = 1, inp%nvar_bcl
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
   real,              intent(in)  :: qprime_df(inp%nvar_bcl, G%npoin, inp%nlayers)
   integer,           intent(out) :: nreq
   integer,           intent(out) :: ireq(2*par%num_nbh)
   integer,           intent(out) :: status(mpi_status_size, 2*par%num_nbh)

   integer :: jj, kk, i, inbh, ib, iface, el, ivar
   integer :: inode, ip, ll, ioff
   integer :: nqp, istart, iend, idest, ierr
   integer :: ngl_f, nbtp_var_f, nlayers_f, nvarb_bcl_f, nboun_valid

   ngl_f        = b%ngl
   nbtp_var_f   = ref%nbtp_var
   nlayers_f    = inp%nlayers
   nvarb_bcl_f  = inp%nvar_bcl
   nboun_valid  = ref%nboun_valid

   !$acc parallel loop gang private(iface, el) &
   !$acc    present(G%face, G%intma, mf%imapl, init%pbprime_df, q, qprime_df, &
   !$acc            send_data_dg, ref%face_pack_list) &
   !$acc    firstprivate(nvarb, ngl_f, nbtp_var_f, nlayers_f, nvarb_bcl_f, nboun_valid)
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
            do ivar = 1, nvarb_bcl_f
               send_data_dg(ioff+nvarb+1+(ll-1)*nvarb_bcl_f+ivar) = qprime_df(ivar, ip, ll)
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
   q, btp_dpp_graduv, pbprime_visc, Uk, nvarb, nvel, &
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
   integer,           intent(in)    :: nvel

   real,    intent(in)  :: q(nvarb, G%npoin)
   real,    intent(in)  :: btp_dpp_graduv(nvarb, G%npoin)
   real,    intent(in)  :: pbprime_visc(G%npoin)
   real,    intent(in)  :: Uk(nvel, G%npoin)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size, 2*par%num_nbh)

   integer :: kk, inode, iface, el, ip, ioff, ivar
   integer :: jj, ib, inbh, nqp, istart, iend, idest, ierr, i
   integer :: ngl_f, nboun_valid_f, nvarb_lap_f
   integer :: off_pbprime, off_graduv, off_visc, off_vel

   ngl_f         = b%ngl
   nboun_valid_f = ref%nboun_valid
   nvarb_lap_f   = ref%nbtp_var_lap

   ! Packed layout per node (nvarb_lap_f = 2*nvarb + 2 + nvel slots):
   !   [1..nvarb]           = q (graduv, local ∇ū)
   !   [nvarb+1]            = pbprime_df (init reference thickness)
   !   [nvarb+2..2*nvarb+1] = btp_dpp_graduv (frozen p̄_b * ∇ū_BTP)
   !   [2*nvarb+2]          = pbprime_visc (current p̄_b, for SIP dp averaging)
   !   [2*nvarb+3..]        = Uk (BTP velocity ū, for SIP [u] jump)
   off_pbprime = nvarb + 1
   off_graduv  = nvarb + 1
   off_visc    = 2*nvarb + 2
   off_vel     = 2*nvarb + 2

   ! GPU packing: one gang per boundary face, one vector lane per node.
   !$acc parallel loop gang private(iface, el)                                      &
   !$acc    present(G%face, G%intma, mf%imapl,                                     &
   !$acc            init%pbprime_df, q, btp_dpp_graduv, pbprime_visc, Uk,          &
   !$acc            ref%send_data_dg_lap, ref%face_pack_list)                       &
   !$acc    firstprivate(nvarb, nvel, ngl_f, nboun_valid_f, nvarb_lap_f,           &
   !$acc                 off_pbprime, off_graduv, off_visc, off_vel)
   do kk = 1, nboun_valid_f
      iface = ref%face_pack_list(kk)
      el    = G%face(7, iface)
      !$acc loop vector private(ip, ioff)
      do inode = 1, ngl_f
         ip   = G%intma(mf%imapl(1,inode,1,iface), &
            mf%imapl(2,inode,1,iface), &
            mf%imapl(3,inode,1,iface), el)
         ioff = (kk-1)*ngl_f*nvarb_lap_f + (inode-1)*nvarb_lap_f
         !$acc loop seq
         do ivar = 1, nvarb
            ref%send_data_dg_lap(ioff+ivar) = q(ivar, ip)
         end do
         ref%send_data_dg_lap(ioff+off_pbprime) = init%pbprime_df(ip)
         !$acc loop seq
         do ivar = 1, nvarb
            ref%send_data_dg_lap(ioff+off_graduv+ivar) = btp_dpp_graduv(ivar, ip)
         end do
         ref%send_data_dg_lap(ioff+off_visc) = pbprime_visc(ip)
         !$acc loop seq
         do ivar = 1, nvel
            ref%send_data_dg_lap(ioff+off_vel+ivar) = Uk(ivar, ip)
         end do
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
            nqp = nqp + b%ngl * nvarb_lap_f
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
   real,              intent(in)  :: qprime_df(inp%nvar_bcl, G%npoin, inp%nlayers)
   integer,           intent(out) :: nreq
   integer,           intent(out) :: ireq(2*par%num_nbh)
   integer,           intent(out) :: status(mpi_status_size, 2*par%num_nbh)

   integer :: jj, inbh, ib, iface, imulti, el, ivar, i
   integer :: inode, ip, ll, ioff
   integer :: nqp, istart, iend, idest, ierr
   integer :: ngl_f, nbtp_var_f, nlayers_f, nvarb_bcl_f, nboun_f

   ngl_f       = b%ngl
   nbtp_var_f  = ref%nbtp_var
   nlayers_f   = inp%nlayers
   nvarb_bcl_f = inp%nvar_bcl
   nboun_f     = G%nboun

   ! GPU pack — one gang per MPI boundary face.
   !$acc parallel loop gang private(iface, imulti, el) &
   !$acc    present(G%face_type, G%face, G%intma, par%nbh_send_recv, &
   !$acc            par%nbh_send_recv_multi, mf%imapl, init%pbprime_df, &
   !$acc            q, qprime_df, send_data_dg) &
   !$acc    firstprivate(nvarb, ngl_f, nbtp_var_f, nlayers_f, nvarb_bcl_f)
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
               do ivar = 1, nvarb_bcl_f
                  send_data_dg(ioff+nvarb+1+(ll-1)*nvarb_bcl_f+ivar) = qprime_df(ivar, ip, ll)
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

   real, intent(in)  :: send_data(inp%nvar_bcl*inp%nlayers*b%ngl*par%num_send_recv_total)
   real, intent(out) :: recv_data(inp%nvar_bcl*inp%nlayers*b%ngl*par%num_send_recv_total)
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
            nqp = nqp + b%ngl * (inp%nvar_bcl*inp%nlayers)
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

subroutine pack_and_send_df_bcl(G, inp, b, mf, par, ref, send_data, recv_data, q, nreq, ireq, status)

   use mod_basis,         only: basis
   use mod_face,          only: face_CS
   use mod_grid,          only: grid
   use mod_parallel,      only: parallel_CS
   use mod_input,         only: input
   use mod_ref,           only: mref
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(face_CS),     intent(in)  :: mf
   type(input),       intent(in)  :: inp
   type(parallel_CS), intent(in)  :: par
   type(mref),        intent(in)  :: ref

   real, intent(out) :: send_data(inp%nvar_bcl*inp%nlayers*b%ngl*par%num_send_recv_total)
   real, intent(out) :: recv_data(inp%nvar_bcl*inp%nlayers*b%ngl*par%num_send_recv_total)
   real, intent(in)  :: q(inp%nvar_bcl, G%npoin, inp%nlayers)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size, 2*par%num_nbh)

   integer :: kk, jj, i, inbh, ib, iface, el, ivar, inode, ip, ll
   integer :: ngl_f, nlayers_f, nvarb_f, nboun_valid_f
   integer :: nqp, istart, iend, idest, ierr

   ngl_f         = b%ngl
   nlayers_f     = inp%nlayers
   nvarb_f       = inp%nvar_bcl
   nboun_valid_f = ref%nboun_valid

   ! GPU pack: each gang handles one boundary face from the precomputed list.
   !$acc parallel loop gang &
   !$acc    present(G%face, G%intma, mf%imapl, q, send_data, ref%face_pack_list) &
   !$acc    firstprivate(ngl_f, nlayers_f, nvarb_f, nboun_valid_f)
   do kk = 1, nboun_valid_f
      iface = ref%face_pack_list(kk)
      el    = G%face(7, iface)
      !$acc loop seq
      do ll = 1, nlayers_f
         !$acc loop vector private(ip)
         do inode = 1, ngl_f
            ip = G%intma(mf%imapl(1,inode,1,iface), &
                         mf%imapl(2,inode,1,iface), &
                         mf%imapl(3,inode,1,iface), el)
            !$acc loop seq
            do ivar = 1, nvarb_f
               send_data((kk-1)*nlayers_f*ngl_f*nvarb_f + (ll-1)*ngl_f*nvarb_f + (inode-1)*nvarb_f + ivar) = q(ivar,ip,ll)
            end do
         end do
      end do
   end do
   !$acc end parallel loop

   ! GPU-direct non-blocking send/receive per neighbor.
   nreq   = 0
   iend   = 0
   jj     = 1
   status = 0

   !$acc host_data use_device(send_data, recv_data)
   do inbh = 1, par%num_nbh
      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         do i = 1, par%nbh_send_recv_multi(jj)
            nqp = nqp + ngl_f * (nvarb_f*nlayers_f)
         end do
         jj = jj + 1
      end do

      idest  = par%nbh_proc(inbh)
      nreq   = nreq + 1
      istart = iend + 1
      iend   = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, ireq(nreq), ierr)
         call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      end if
      nreq = nreq + 1
   end do
   !$acc end host_data

end subroutine pack_and_send_df_bcl

subroutine pack_and_send_df_bcl_lap(G, inp, b, mf, par, ref, send_data, recv_data, &
                                     dpp_graduv, dpprime_visc, qprime_df, nreq, ireq, status)

   use mod_basis,         only: basis
   use mod_face,          only: face_CS
   use mod_grid,          only: grid
   use mod_parallel,      only: parallel_CS
   use mod_input,         only: input
   use mod_ref,           only: mref
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(face_CS),     intent(in)  :: mf
   type(input),       intent(in)  :: inp
   type(parallel_CS), intent(in)  :: par
   type(mref),        intent(in)  :: ref

   real, intent(out) :: send_data(ref%nbcl_var_lap*b%ngl*par%num_send_recv_total)
   real, intent(out) :: recv_data(ref%nbcl_var_lap*b%ngl*par%num_send_recv_total)
   real, intent(in)  :: dpp_graduv(inp%ngraduvw_var, G%npoin, inp%nlayers)
   real, intent(in)  :: dpprime_visc(G%npoin, inp%nlayers)
   real, intent(in)  :: qprime_df(inp%nvar_bcl, G%npoin, inp%nlayers)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size, 2*par%num_nbh)

   integer :: kk, jj, i, inbh, ib, iface, el, ivar, ip, ll, inode
   integer :: ngl_f, nlayers_f, nboun_valid_f, nvarb_lap_per_layer, nw_f, nvel_f
   integer :: nqp, istart, iend, idest, ierr, base

   ngl_f               = b%ngl
   nlayers_f           = inp%nlayers
   nboun_valid_f       = ref%nboun_valid
   nvarb_lap_per_layer = ref%nbcl_var_lap / inp%nlayers
   nw_f                = inp%ngraduvw_var
   nvel_f              = inp%nvar_bcl - 1   ! velocity components (nvar_bcl-1: cartesian=2, sphere_ico/sphere_hex=3)

   ! Buffer layout per layer per quad node (nvarb_lap_per_layer = nw_f+1+nvel_f slots):
   !   1..nw_f      : dpp_graduvw (dp'_k * grad u'_k)
   !   nw_f+1       : dpprime_visc (dp'_k for penalty averaging)
   !   nw_f+2..end  : qprime_df(2..nvel_f+1) = velocity deviation u'_k
   !$acc parallel loop gang &
   !$acc    present(G%face, G%intma, mf%imapl, dpp_graduv, dpprime_visc, qprime_df, send_data, ref%face_pack_list) &
   !$acc    firstprivate(ngl_f, nlayers_f, nboun_valid_f, nvarb_lap_per_layer, nw_f, nvel_f)
   do kk = 1, nboun_valid_f
      iface = ref%face_pack_list(kk)
      el    = G%face(7, iface)
      !$acc loop seq
      do ll = 1, nlayers_f
         !$acc loop vector private(ip, base)
         do inode = 1, ngl_f
            ip   = G%intma(mf%imapl(1,inode,1,iface), &
                           mf%imapl(2,inode,1,iface), &
                           mf%imapl(3,inode,1,iface), el)
            base = (kk-1)*nlayers_f*ngl_f*nvarb_lap_per_layer &
                 + (ll-1)*ngl_f*nvarb_lap_per_layer            &
                 + (inode-1)*nvarb_lap_per_layer
            !$acc loop seq
            do ivar = 1, nw_f
               send_data(base + ivar) = dpp_graduv(ivar,ip,ll)
            end do
            send_data(base + nw_f + 1) = dpprime_visc(ip,ll)
            !$acc loop seq
            do ivar = 1, nvel_f
               send_data(base + nw_f + 1 + ivar) = qprime_df(ivar+1,ip,ll)
            end do
         end do
      end do
   end do
   !$acc end parallel loop

   nreq   = 0
   iend   = 0
   jj     = 1
   status = 0

   !$acc host_data use_device(send_data, recv_data)
   do inbh = 1, par%num_nbh
      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         do i = 1, par%nbh_send_recv_multi(jj)
            nqp = nqp + ngl_f * ref%nbcl_var_lap
         end do
         jj = jj + 1
      end do

      idest  = par%nbh_proc(inbh)
      nreq   = nreq + 1
      istart = iend + 1
      iend   = istart + nqp - 1

      if (nqp > 0) then
         call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, ireq(nreq), ierr)
         call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      end if
      nreq = nreq + 1
   end do
   !$acc end host_data

end subroutine pack_and_send_df_bcl_lap

! ---------------------------------------------------------------------------
! pack_and_send_df_bcl_gradz
!   CPU-only MPI pack+send of grad_z_df (3 components per DOF per layer) at
!   MPI boundary face DOFs.  Mirrors pack_and_send_df_bcl_lap but carries
!   only the 3-component LDG gradient of z_k (no velocity gradient, no visc).
!   Buffer layout per face (kk), per layer (ll), per face node (inode):
!     3 values: grad_z_df(1:3, ip, ll)
! ---------------------------------------------------------------------------
subroutine pack_and_send_df_bcl_gradz(G, inp, b, mf, par, ref, send_data, recv_data, &
                                       grad_z_df, nreq, ireq, status)

   use mod_basis,         only: basis
   use mod_face,          only: face_CS
   use mod_grid,          only: grid
   use mod_parallel,      only: parallel_CS
   use mod_input,         only: input
   use mod_ref,           only: mref
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   type(grid),        intent(in)  :: G
   type(basis),       intent(in)  :: b
   type(face_CS),     intent(in)  :: mf
   type(input),       intent(in)  :: inp
   type(parallel_CS), intent(in)  :: par
   type(mref),        intent(in)  :: ref

   real, intent(out) :: send_data(ref%nbcl_var_gradz*b%ngl*par%num_send_recv_total)
   real, intent(out) :: recv_data(ref%nbcl_var_gradz*b%ngl*par%num_send_recv_total)
   real, intent(in)  :: grad_z_df(3, G%npoin, inp%nlayers)
   integer, intent(out) :: nreq
   integer, intent(out) :: ireq(2*par%num_nbh)
   integer, intent(out) :: status(mpi_status_size, 2*par%num_nbh)

   integer :: kk, jj, ll, i, inode, ip, inbh, ib, nqp, istart, iend, idest, ierr, base, iface, el
   integer :: ngl_f, nlayers_f, nboun_valid_f

   ngl_f         = b%ngl
   nlayers_f     = inp%nlayers
   nboun_valid_f = ref%nboun_valid

   ! Pack: grad_z_df(1:3, ip, ll) for each MPI face / layer / face-node.
   do kk = 1, nboun_valid_f
      iface = ref%face_pack_list(kk)
      el    = G%face(7, iface)
      do ll = 1, nlayers_f
         do inode = 1, ngl_f
            ip   = G%intma(mf%imapl(1,inode,1,iface), &
                           mf%imapl(2,inode,1,iface), &
                           mf%imapl(3,inode,1,iface), el)
            base = (kk-1)*nlayers_f*ngl_f*3 + (ll-1)*ngl_f*3 + (inode-1)*3
            send_data(base+1) = grad_z_df(1, ip, ll)
            send_data(base+2) = grad_z_df(2, ip, ll)
            send_data(base+3) = grad_z_df(3, ip, ll)
         end do
      end do
   end do

   nreq   = 0
   iend   = 0
   jj     = 1
   status = 0

   do inbh = 1, par%num_nbh
      nqp = 0
      do ib = 1, par%num_send_recv(inbh)
         do i = 1, par%nbh_send_recv_multi(jj)
            nqp = nqp + ngl_f * ref%nbcl_var_gradz
         end do
         jj = jj + 1
      end do
      idest  = par%nbh_proc(inbh)
      nreq   = nreq + 1
      istart = iend + 1
      iend   = istart + nqp - 1
      if (nqp > 0) then
         call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, ireq(nreq), ierr)
         call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION, idest-1, 99, mpi_comm_world, ireq(nreq+1), ierr)
      else
         ireq(nreq)   = MPI_REQUEST_NULL
         ireq(nreq+1) = MPI_REQUEST_NULL
      end if
      nreq = nreq + 1
   end do

end subroutine pack_and_send_df_bcl_gradz

! ---------------------------------------------------------------------------
! unpack_data_dg_general_gradz_bcl
!   Unpack flat MPI buffer into q_send/q_recv(3*nlayers, ngl, nboun).
!   Mirrors unpack_data_dg_general_lap_bcl.
! ---------------------------------------------------------------------------
subroutine unpack_data_dg_general_gradz_bcl(G, b, par, ref, q_send, q_recv, &
                                              send_data, recv_data, nlayers, nboun_valid)

   use mod_basis,    only: basis
   use mod_grid,     only: grid
   use mod_parallel, only: parallel_CS
   use mod_ref,      only: mref

   implicit none

   type(grid),        intent(in) :: G
   type(basis),       intent(in) :: b
   type(parallel_CS), intent(in) :: par
   type(mref),        intent(in) :: ref
   integer,           intent(in) :: nlayers, nboun_valid

   real, dimension(3*nlayers, b%ngl, par%num_send_recv_total), intent(out) :: q_send, q_recv
   real, dimension(3*nlayers*b%ngl*par%num_send_recv_total),   intent(in)  :: send_data, recv_data

   integer :: kk, ll, inode, ivar, ngl_f, nlayers_f, nboun_valid_f

   ngl_f         = b%ngl
   nlayers_f     = nlayers
   nboun_valid_f = nboun_valid

   do kk = 1, nboun_valid_f
      do ll = 1, nlayers_f
         do inode = 1, ngl_f
            do ivar = 1, 3
               q_send((ll-1)*3+ivar, inode, kk) = send_data((kk-1)*nlayers_f*ngl_f*3 + (ll-1)*ngl_f*3 + (inode-1)*3 + ivar)
               q_recv((ll-1)*3+ivar, inode, kk) = recv_data((kk-1)*nlayers_f*ngl_f*3 + (ll-1)*ngl_f*3 + (inode-1)*3 + ivar)
            end do
         end do
      end do
   end do

end subroutine unpack_data_dg_general_gradz_bcl

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
