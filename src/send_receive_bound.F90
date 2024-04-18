!----------------------------------------------------------------------!
!>@brief Communication using Non-Blocking Send/Receives
!>Comes from NSEAM code.
!>@author  F.X. Giraldo
!>      Department of Applied Mathematics
!>      Naval Postgraduate School
!>      Monterey, CA 93943
!----------------------------------------------------------------------!

subroutine unpack_data_dg_general_quad(q_send,q_recv,send_data,recv_data,nvarb)

    use mod_basis, only: nq, FACE_CHILDREN

    use mod_grid, only: nboun, face, mod_grid_get_face_nq, face_type

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi
    
    !use mod_ref, only: nfields => nmessage

    implicit none

    !Global Variables
    real, dimension(nvarb,nq,nboun), intent(out) :: q_send,q_recv
    real, dimension(nvarb*nq*nboun), intent(in)  :: send_data, recv_data
    integer, intent(in) :: nvarb

    !Local Variables
    integer ii, jj, kk, i, inbh, ib, ifaces, inode, jnode, ivar, ilocl, ilocr
    integer nq_i, nq_j, plane_ij, iface, imulti, ftype

    ii = 0
    jj = 1
    kk = 1
    
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
            ftype = face_type(iface)

            ! if (ftype == 2 .and. imulti>0) then
                ilocl = face(5,iface)

                call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)

                do inode = 1,nq
                    do ivar = 1,nvarb
                        ii = ii + 1
                        q_send(ivar,inode,kk) = send_data(ii)
                        q_recv(ivar,inode,kk) = recv_data(ii)
                    end do
                end do
                kk=kk+1
            ! end if
            jj = jj + 1
        end do
    end do

end subroutine unpack_data_dg_general_quad

subroutine unpack_data_dg_general_df(q_send,q_recv,send_data,recv_data,nvarb)

    use mod_basis, only: ngl, FACE_CHILDREN

    use mod_grid, only: nboun, face, mod_grid_get_face_nq, face_type

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi
    
    !use mod_ref, only: nfields => nmessage

    implicit none

    !Global Variables
    real, dimension(nvarb,ngl,nboun), intent(out) :: q_send,q_recv
    real, dimension(nvarb*ngl*nboun), intent(in)  :: send_data, recv_data
    integer, intent(in) :: nvarb

    !Local Variables
    integer ii, jj, kk, i, inbh, ib, ifaces, inode, jnode, ivar, ilocl, ilocr
    integer nq_i, nq_j, plane_ij, iface, imulti, ftype

    ii = 0
    jj = 1
    kk = 1
    
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
            ftype = face_type(iface)

            ! if (ftype == 2 .and. imulti>0) then
                ilocl = face(5,iface)

                do inode = 1,ngl
                    do ivar = 1,nvarb
                        ii = ii + 1
                        q_send(ivar,inode,kk) = send_data(ii)
                        q_recv(ivar,inode,kk) = recv_data(ii)
                    end do
                end do
                kk=kk+1
            ! end if
            jj = jj + 1
        end do
    end do

end subroutine unpack_data_dg_general_df

subroutine unpack_data_dg_general_quad_layer(q_send,q_recv,send_data,recv_data,nvarb,nlayers,nq)

    use mod_basis, only: FACE_CHILDREN

    use mod_grid, only: nboun, face, mod_grid_get_face_nq, face_type

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi
    
    !use mod_ref, only: nfields => nmessage

    implicit none

    !Global Variables
    real, dimension(nvarb,nq,nboun,nlayers), intent(out) :: q_send,q_recv
    real, dimension(nvarb*nq*nboun*nlayers), intent(in)  :: send_data, recv_data
    integer, intent(in) :: nvarb,nlayers,nq

    !Local Variables
    integer ii, jj, kk, i, inbh, ib, ifaces, inode, jnode, ivar, ilocl, ilocr
    integer nq_i, nq_j, plane_ij, iface, imulti, ftype, ll

    ii = 0
    jj = 1
    kk = 1
    
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
            ftype = face_type(iface)

            if (ftype == 2 .and. imulti>0) then
                ilocl = face(5,iface)

                call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)

                do ll = 1,nlayers

                    do inode = 1,nq
                        do ivar = 1,nvarb
                            ii = ii + 1
                            q_send(ivar,inode,kk,ll) = send_data(ii)
                            q_recv(ivar,inode,kk,ll) = recv_data(ii)
                        end do
                    end do
                end do
                kk=kk+1
            end if
            jj = jj + 1
        
        end do
    end do

end subroutine unpack_data_dg_general_quad_layer

subroutine unpack_data_dg_general_quad_1v(q_send,q_recv,send_data,recv_data)

    use mod_basis, only: nq, FACE_CHILDREN

    use mod_grid, only: nboun, face, mod_grid_get_face_nq, face_type

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi
    
    !use mod_ref, only: nfields => nmessage

    implicit none

    !Global Variables
    real, dimension(nq,nboun), intent(out) :: q_send,q_recv
    real, dimension(nq*nboun), intent(in)  :: send_data, recv_data

    !Local Variables
    integer ii, jj, kk, i, inbh, ib, ifaces, inode, jnode, ivar, ilocl, ilocr
    integer nq_i, nq_j, plane_ij, iface, imulti, ftype

    ii = 0
    jj = 1
    kk = 1
    
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
            ftype = face_type(iface)

            if (ftype == 2 .and. imulti > 0) then
                ilocl = face(5,iface)

                call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)

                do inode = 1,nq
                    ii = ii + 1
                    q_send(inode,kk) = send_data(ii)
                    q_recv(inode,kk) = recv_data(ii)
                end do
                kk=kk+1
            end if
            jj = jj + 1
        
        end do
    end do

end subroutine unpack_data_dg_general_quad_1v


subroutine pack_data_dg_quad(q_send,q_face,nvarb)
  
    use mod_basis, only: ngl, FACE_CHILDREN, nq

    use mod_face, only: face_send, imapl_q, imapr_q, normal_vector_q, jac_faceq

    use mod_grid, only: nelem, npoin, intma_dg_quad, face_type, nboun, face, mod_grid_get_face_nq, nface

    use mod_initial, only: nvar, q_ref, q_ref_layers, rho_layers

    use mod_input, only: llimit, limit_threshold, is_shallow, eqn_set, is_swe_layers

    use mod_metrics, only: jac

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi
  
    use mod_viscosity, only: grad_q

    implicit none
  
    !Global Variables
    real, intent(out) :: q_send(nvarb*nq*nboun)
    real, intent(in) :: q_face(nvarb,2,nq,nface)
    integer, intent(in) :: nvarb

    !Local Variables
    integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
    integer :: nq_i, nq_j, plane_ij
    real :: h, qu, qv
    integer :: inode, jnode, ip, ilocl, ilocr
    integer :: iface_type
    logical :: is_set2nc

    !check for Set2NC which is Non-Conservative
    is_set2nc = (eqn_set == 'set2nc')

    ii = 0
    jj = 1
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
            ! jj = jj + 1

            if(face_type(iface) == 2 .and. imulti > 0) then

                ilocl=face(5,iface)
                el = face(7,iface)      ! Get Element

                do inode = 1,nq
                    
                    ! Load primitive variables
                    do ivar=1,nvarb
                        ii = ii + 1
                        q_send(ii)=q_face(ivar,1,inode,iface) 
                    end do
                end do
            end if
            jj = jj + 1
        end do
    end do
  
end subroutine pack_data_dg_quad

subroutine pack_data_dg_df(q_send,q_face,nvarb)
  
    use mod_basis, only: ngl, FACE_CHILDREN

    use mod_face, only: face_send, imapl_q, imapr_q, normal_vector_q, jac_faceq

    use mod_grid, only: nelem, npoin, intma_dg_quad, face_type, nboun, face, mod_grid_get_face_nq, nface

    use mod_initial, only: nvar, q_ref, q_ref_layers, rho_layers

    use mod_input, only: llimit, limit_threshold, is_shallow, eqn_set, is_swe_layers

    use mod_metrics, only: jac

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi
  
    use mod_viscosity, only: grad_q

    implicit none
  
    !Global Variables
    real, intent(out) :: q_send(nvarb*ngl*nboun)
    real, intent(in) :: q_face(nvarb,2,ngl,nface)
    integer, intent(in) :: nvarb

    !Local Variables
    integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
    integer :: nq_i, nq_j, plane_ij
    real :: h, qu, qv
    integer :: inode, jnode, ip, ilocl, ilocr
    integer :: iface_type
    logical :: is_set2nc

    !check for Set2NC which is Non-Conservative
    is_set2nc = (eqn_set == 'set2nc')

    ii = 0
    jj = 1
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
            ! jj = jj + 1

            if(face_type(iface) == 2 .and. imulti > 0) then

                ilocl=face(5,iface)
                el = face(7,iface)      ! Get Element

                do inode = 1,ngl
                    
                    ! Load primitive variables
                    do ivar=1,nvarb
                        ii = ii + 1
                        q_send(ii)=q_face(ivar,1,inode,iface) 
                    end do
                end do
            end if
            jj = jj + 1
        end do
    end do
  
end subroutine pack_data_dg_df

subroutine pack_data_dg_quad_all(q_send,q_face,grad_face,nvarb)
  
    use mod_basis, only: ngl, FACE_CHILDREN, nq

    use mod_face, only: face_send, imapl_q, imapr_q, normal_vector_q, jac_faceq

    use mod_grid, only: nelem, npoin, intma_dg_quad, face_type, nboun, face, mod_grid_get_face_nq, nface

    use mod_initial, only: nvar, q_ref, q_ref_layers, rho_layers

    use mod_input, only: llimit, limit_threshold, is_shallow, eqn_set, is_swe_layers

    use mod_metrics, only: jac

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi
  
    use mod_viscosity, only: grad_q

    implicit none
  
    !Global Variables
    real, intent(out) :: q_send(2*nvarb*nq*nboun)
    real, intent(in) :: q_face(nvarb,2,nq,nface)
    real, intent(in) :: grad_face(nvarb,2,nq,nface)
    integer, intent(in) :: nvarb

    !Local Variables
    integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
    integer :: nq_i, nq_j, plane_ij
    real :: h, qu, qv
    integer :: inode, jnode, ip, ilocl, ilocr
    integer :: iface_type
    logical :: is_set2nc

    !check for Set2NC which is Non-Conservative
    is_set2nc = (eqn_set == 'set2nc')

    ii = 0
    jj = 1
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
            ! jj = jj + 1

            if(face_type(iface) == 2 .and. imulti > 0) then

                ilocl=face(5,iface)
                el = face(7,iface)      ! Get Element

                do inode = 1,nq
                    
                    ! Load primitive variables
                    do ivar=1,nvarb
                        ii = ii + 1
                        q_send(ii)=q_face(ivar,1,inode,iface) 
                    end do

                    do ivar=1,nvarb
                        ii = ii + 1
                        q_send(ii)=grad_face(ivar,1,inode,iface) 
                    end do
                end do
            end if
            jj = jj + 1
        end do
    end do
  
end subroutine pack_data_dg_quad_all

subroutine pack_data_dg_quad_lap(q_send,grad_uvdp,nvarb)
  
    use mod_basis, only: ngl, FACE_CHILDREN, nq

    use mod_face, only: face_send, imapl_q, imapr_q, normal_vector_q, jac_faceq

    use mod_grid, only: nelem, npoin, intma_dg_quad, face_type, npoin_q, nboun, face, mod_grid_get_face_nq, nface

    use mod_initial, only: nvar, q_ref, q_ref_layers, rho_layers

    use mod_input, only: llimit, limit_threshold, is_shallow, eqn_set, is_swe_layers

    use mod_metrics, only: jac

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi
  
    use mod_viscosity, only: grad_q

    implicit none
  
    !Global Variables
    real, intent(out) :: q_send(nvarb*nq*nboun)
    real, intent(in) :: grad_uvdp(nvarb,npoin_q)
    integer, intent(in) :: nvarb

    !Local Variables
    integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
    integer :: nq_i, nq_j, plane_ij
    real :: h, qu, qv
    integer :: inode, jnode, ip, ilocl, ilocr
    integer :: iface_type
    logical :: is_set2nc

    !check for Set2NC which is Non-Conservative
    is_set2nc = (eqn_set == 'set2nc')

    ii = 0
    jj = 1
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
            ! jj = jj + 1

            !if(face_type(iface) == 2 .and. imulti > 0) then

                ilocl=face(5,iface)
                el = face(7,iface)      ! Get Element

                do inode = 1,nq

                    il = imapr_q(1,inode,1,iface);
                    jl = imapr_q(2,inode,1,iface)
                    kl = imapr_q(3,inode,1,iface)
        
                    ip=intma_dg_quad(il,jl,kl,el)
                    
                    ! Load primitive variables
                    do ivar=1,nvarb
                        ii = ii + 1
                        q_send(ii)=grad_uvdp(ivar,ip) 
                    end do
                end do
            !end if
            jj = jj + 1
        end do
    end do
  
end subroutine pack_data_dg_quad_lap

subroutine pack_data_dg_quad_layer(q_send,q_face,nvarb,nlayers,nq)
  
    use mod_basis, only: FACE_CHILDREN

    use mod_grid, only: nelem, npoin, intma_dg_quad, face_type, nboun, face, mod_grid_get_face_nq, nface

    use mod_initial, only: nvar, q_ref, q_ref_layers, rho_layers

    use mod_input, only: llimit, limit_threshold, is_shallow, eqn_set, is_swe_layers

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi

    implicit none
  
    !Global Variables
    real, intent(out) :: q_send(nvarb*nq*nboun*nlayers)
    real, intent(in) :: q_face(nvarb,2,nq,nface,nlayers)
    integer, intent(in) :: nvarb, nlayers, nq

    !Local Variables
    integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
    integer :: nq_i, nq_j, plane_ij, ll
    real :: h, qu, qv
    integer :: inode, jnode, ip, ilocl, ilocr
    integer :: iface_type
    logical :: is_set2nc

    !check for Set2NC which is Non-Conservative
    is_set2nc = (eqn_set == 'set2nc')

    ii = 0
    jj = 1
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
        
            if(face_type(iface)== 2 .and. imulti > 0) then

                ilocl=face(5,iface)
                el = face(7,iface)      ! Get Element

                do ll = 1,nlayers

                    do inode = 1,nq
                        
                        ! Load primitive variables
                        do ivar=1,nvarb
                            ii = ii + 1
                            q_send(ii)=q_face(ivar,1,inode,iface,ll) 
                        end do

                    end do
                end do
            end if
            jj = jj + 1
        end do
    end do
  
end subroutine pack_data_dg_quad_layer

subroutine pack_data_dg_quad_layer_all(q_send,q_face,qprime_face,nvarb,nlayers)
  
    use mod_basis, only: ngl, FACE_CHILDREN, nq

    use mod_grid, only: nelem, npoin, intma_dg_quad, face_type, nboun, face, mod_grid_get_face_nq, nface

    use mod_initial, only: nvar, q_ref, q_ref_layers, rho_layers

    use mod_input, only: llimit, limit_threshold, is_shallow, eqn_set, is_swe_layers

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi

    implicit none
  
    !Global Variables
    real, intent(out) :: q_send(2*nvarb*nq*nboun*nlayers)
    real, dimension(nvarb,2,nq,nface,nlayers), intent(in) :: q_face, qprime_face
    integer, intent(in) :: nvarb, nlayers

    !Local Variables
    integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
    integer :: nq_i, nq_j, plane_ij, ll
    real :: h, qu, qv
    integer :: inode, jnode, ip, ilocl, ilocr
    integer :: iface_type
    logical :: is_set2nc

    !check for Set2NC which is Non-Conservative
    is_set2nc = (eqn_set == 'set2nc')

    ii = 0
    jj = 1
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
        
            if(face_type(iface)== 2 .and. imulti > 0) then

                ilocl=face(5,iface)
                el = face(7,iface)      ! Get Element

                do ll = 1,nlayers

                    do inode = 1,nq
                        
                        ! Load primitive variables
                        do ivar = 1,nvarb
                            ii = ii + 1
                            q_send(ii)=q_face(ivar,1,inode,iface,ll) 
                        end do

                        do ivar = 1, nvarb
                            ii = ii + 1
                            q_send(ii)=qprime_face(ivar,1,inode,iface,ll) 
                        end do

                    end do
                end do
            end if
            jj = jj + 1
        end do
    end do
  
end subroutine pack_data_dg_quad_layer_all

subroutine pack_data_dg_quad_1v(q_send,q_face)
  
    use mod_basis, only: ngl, FACE_CHILDREN, nq

    use mod_face, only: face_send, imapl_q, imapr_q, normal_vector_q, jac_faceq

    use mod_grid, only: nelem, npoin, intma_dg_quad, face_type, nboun, face, mod_grid_get_face_nq, nface

    use mod_initial, only: nvar, q_ref, q_ref_layers, rho_layers

    use mod_input, only: llimit, limit_threshold, is_shallow, eqn_set, is_swe_layers

    use mod_metrics, only: jac

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi
  
    use mod_viscosity, only: grad_q

    implicit none
  
    !Global Variables
    real, intent(out) :: q_send(nq*nboun)
    real, intent(in) :: q_face(2,nq,nface)

    !Local Variables
    integer :: ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
    integer :: nq_i, nq_j, plane_ij
    real :: h, qu, qv
    integer :: inode, jnode, ip, ilocl, ilocr
    integer :: iface_type
    logical :: is_set2nc

    !check for Set2NC which is Non-Conservative
    is_set2nc = (eqn_set == 'set2nc')

    ii = 0
    jj = 1

    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
            ! jj = jj + 1
        
            if(face_type(iface) == 2 .and. imulti > 0) then

                ilocl=face(5,iface)
                el = face(7,iface)      ! Get Element

                do inode = 1,nq
                    
                    ! Load primitive variables
                    ii = ii + 1
                    q_send(ii)=q_face(1,inode,iface) 

                end do
            end if

            jj = jj + 1
        
        end do
    end do

end subroutine pack_data_dg_quad_1v

subroutine pack_data_dg_laplacian_mlswe(q_send,q,nfields)
  
    use mod_basis, only: ngl, FACE_CHILDREN

    use mod_face, only: face_send, imapl, imapr, normal_vector, jac_face

    use mod_grid, only: nelem, npoin, intma, face_type, nboun, face, mod_grid_get_face_ngl

    use mod_initial, only: nvar

    use mod_metrics, only: jac

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi
  
    use mod_viscosity, only: grad_q_mlswe, nvar_visc

    implicit none
  
    !Global Variables
    real, intent(out) :: q_send(nfields*ngl*ngl*nboun)
    real, intent(in) :: q(nvar_visc,npoin)
    integer, intent(in) :: nfields

    !Local Variables
    integer ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
    integer ngl_i, ngl_j, plane_ij
    real rho, u, v, w, theta, theta_prime, p, ptot, nx, ny, nz, rhoinv
    integer inode, jnode, ip, ilocl, ilocr
    integer iface_type
    integer illoc(ngl,ngl), jlloc(ngl,ngl), klloc(ngl,ngl)
    real q_t(nfields,ngl,ngl) !nfields=nmessage

    ii = 0
    jj = 1
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
            jj = jj + 1
        
            if(face_type(iface)==2.and.imulti>0) then
                ilocl=face(5,iface)
                el = face(7,iface)      ! Get Element

                call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

                do inode = 1,ngl_i
                    do jnode = 1,ngl_j
                        il = imapl(1,inode,jnode,iface)
                        jl = imapl(2,inode,jnode,iface)
                        kl = imapl(3,inode,jnode,iface)

                        ip=intma(il,jl,kl,el)

                        !Load Variables
                        do ivar=1,nvar_visc
                            ii = ii + 1
                            q_send(ii)=q(ivar,ip) !1 -> nvar_visc
                        end do

                        !--Flux of Gradient
                        nx=normal_vector(1,inode,jnode,iface)
                        ny=normal_vector(2,inode,jnode,iface)
                        nz=normal_vector(3,inode,jnode,iface)
                        do ivar=1,nvar_visc !nvar_visc,...,nvar_visc + 2
                            ii=ii + 1
                            q_send(ii)=nx*grad_q_mlswe(1,ivar,ip) + ny*grad_q_mlswe(2,ivar,ip) + nz*grad_q_mlswe(3,ivar,ip)
                        end do !ivar

                        !SIPG Constant
                        ii = ii + 1
                        q_send(ii)=jac_face(inode,jnode,iface)/jac(il,jl,kl,el) !nvar+3

                    end do
                end do
            end if 
        
        end do
    end do
  
end subroutine pack_data_dg_laplacian_mlswe



!-----------------------------------------------------------------!
!>@brief For (hyper) Laplacian, need 2,...,NVAR variables and N*Gradients.
!>@author James F. Kelly, 12 November 2010
!>@date Feb 2015, Daniel S. Abdi, Extended data size for total pressure
!>@date Sept 2015, F.X. Giraldo, Extended data size for Tracers
!>@date Jan 2016, M.A. Kopera, Modified to handle non-conforming faces
!>@date August 2016, F.X. Giraldo, Updated for the Hyper-diffusion operators for DG
!>@date December 27 2016, F.X. Giraldo, Updated for PISO_FXG
!-----------------------------------------------------------------!
subroutine pack_data_lhs_PISO_fxg(q_send,q,nfields)
  
    use mod_basis, only: ngl, FACE_CHILDREN

    use mod_face, only: face_send, imapl, imapr, normal_vector, jac_face

    use mod_grid, only: nelem, npoin, intma, face_type, nboun, face, mod_grid_get_face_ngl

    use mod_initial, only: nvar

    use mod_metrics, only: jac

    use mod_parallel, only: num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi
  
    use mod_viscosity, only: grad_q

    implicit none
  
    !Global Variables
    real, intent(out) :: q_send(nfields*ngl*ngl*nboun)
    real, intent(in) :: q(npoin)
    integer, intent(in) :: nfields

    !Local Variables
    integer ii, jj, i, inbh, ib, iface, imulti, el, il, jl, kl, ivar
    integer ngl_i, ngl_j, plane_ij
    real rho, u, v, w, theta, theta_prime, p, ptot, nx, ny, nz, rhoinv
    integer inode, jnode, ip, ilocl, ilocr
    integer iface_type
    integer illoc(ngl,ngl), jlloc(ngl,ngl), klloc(ngl,ngl)
    real q_t(3,ngl,ngl) 

    ii = 0
    jj = 1
    do inbh = 1,num_nbh
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)
            jj = jj + 1
        
            if(face_type(iface)==2.and.imulti>0) then
                ilocl=face(5,iface)
                el = face(7,iface)      ! Get Element

                call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

                do inode = 1,ngl_i
                    do jnode = 1,ngl_j
                        il = imapl(1,inode,jnode,iface)
                        jl = imapl(2,inode,jnode,iface)
                        kl = imapl(3,inode,jnode,iface)

                        ip=intma(il,jl,kl,el)

                        !Load Variables
                        ii = ii + 1
                        q_send(ii)=q(ip) !1 -> nvar
                 
                        !--Flux of Gradient
                        nx=normal_vector(1,inode,jnode,iface)
                        ny=normal_vector(2,inode,jnode,iface)
                        nz=normal_vector(3,inode,jnode,iface)
                        do ivar=1,1 !nvar+4,...,2*nvar + 2
                           ii=ii + 1
                           q_send(ii)=nx*grad_q(1,ivar,ip) + ny*grad_q(2,ivar,ip) + nz*grad_q(3,ivar,ip)
                        end do !ivar

                        !SIPG Constant
                        ii = ii + 1
                        q_send(ii)=jac_face(inode,jnode,iface)/jac(il,jl,kl,el) !2*nvar+3

                    end do
                end do

            else if (face_type(iface) == 12.and.imulti>0) then

                !project parent face

                el = face(7,iface)      ! Get Element
                ilocl=face(5,iface)
              
                call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)
      
                do inode = 1,ngl_i
                    do jnode = 1,ngl_j
                        il = imapl(1,inode,jnode,iface)
                        jl = imapl(2,inode,jnode,iface)
                        kl = imapl(3,inode,jnode,iface)
                 
                        ip=intma(il,jl,kl,el)
                       
                        !Variables
                        q_t(1,inode,jnode)=q(ip) 

                        !--Flux of Gradient
                        nx=normal_vector(1,inode,jnode,iface)
                        ny=normal_vector(2,inode,jnode,iface)
                        nz=normal_vector(3,inode,jnode,iface)
                        do ivar=1,1 !nvar+4,...,2*nvar + 2
                            q_t(1+ivar-1,inode,jnode)=nx*grad_q(1,ivar,ip) + ny*grad_q(2,ivar,ip) + nz*grad_q(3,ivar,ip)
                        end do !ivar

                        !SIPG Constant
                        q_t(3,inode,jnode)=jac_face(inode,jnode,iface)/jac(il,jl,kl,el) !2*nvar+3

                    end do
                end do

                do i=1,imulti !send multiple copies of myself to children
                    do inode=1,ngl_i
                        do jnode=1,ngl_j

                            do ivar=1,3
                                ! Load physical variables
                                ii = ii + 1
                                q_send(ii) = q_t(ivar,inode,jnode)
                            end do
                    
                        end do
                    end do
              
                end do

            else if (face_type(iface) == 21.and.imulti>0) then

                do i=1,FACE_CHILDREN !go over four possible children and pack the ones that need to be sent
                    el = face(7+i,iface)      ! Get Element
                    if(el>0)then
                        ilocr=face(6,iface)
                 
                        call mod_grid_get_face_ngl(ilocr, ngl_i, ngl_j, plane_ij)

                        do inode = 1,ngl_i
                            do jnode = 1,ngl_j
                                il = imapr(1,inode,jnode,iface);
                                jl = imapr(2,inode,jnode,iface)
                                kl = imapr(3,inode,jnode,iface)
                 
                                ip=intma(il,jl,kl,el)

                                !Variables
                                ii = ii + 1
                                q_send(ii)=q(ip) 

                                !--Flux of Gradient
                                nx=normal_vector(1,inode,jnode,iface)
                                ny=normal_vector(2,inode,jnode,iface)
                                nz=normal_vector(3,inode,jnode,iface)
                                do ivar=1,1
                                    ii=ii + 1
                                    q_send(ii)=nx*grad_q(1,ivar,ip) + ny*grad_q(2,ivar,ip) + nz*grad_q(3,ivar,ip)
                                end do !ivar
                       
                                !SIPG Constant
                                ii = ii + 1
                                q_send(ii)=jac_face(inode,jnode,iface)/jac(il,jl,kl,el) !2*nvar+3

                            end do
                        end do

                    end if
                end do

            end if
        
        end do
    end do
  
end subroutine pack_data_lhs_PISO_fxg

!----------------------------------------------------------------------!
!>@brief Send boundary data: Not USED => replaced by SEND_BOUND_DG_GENERAL
!>@author James F. Kelly
!>        Department of Applied Mathematics
!>        Naval Postgraduate School
!>        Monterey, CA 93943
!>@date Feb 2015, Daniel S. Abdi, Extended data size for total pressure
!>@date October 9, 2015, F.X. Giraldo, modified to handle LHS operators
!>@date January 7, 2016, M.A. Kopera, modified to handle AMR
!----------------------------------------------------------------------!
subroutine send_boundary_lhs(send_data,recv_data,nreq,ireq,status)

    use mod_basis, only: ngl

    use mod_face, only: face_send

    use mod_grid, only: nboun, face, mod_grid_get_face_ngl, face_type

    use mod_initial, only: nvar

    use mod_parallel, only: nbh_proc, num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi

    use mod_ref, only: nmessage

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    !global variables
    real, intent(in)  :: send_data(nmessage*ngl*ngl*nboun)
    real, intent(out) :: recv_data(nmessage*ngl*ngl*nboun)
    integer, intent(out) :: nreq
    integer, intent(out) :: ireq(2*num_nbh)
    integer, intent(out) :: status(mpi_status_size,2*num_nbh)

    !local variables
    integer inbh, idest, istart, iend, ierr, nq
    integer ngl_i, ngl_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

    nreq = 0
    iend = 0
    jj = 1
    status=0

    do inbh = 1, num_nbh

        !Determine size
        nq = 0
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            ftype = face_type(iface)

            ilocl = face(5,iface)
            if(ftype==21) ilocl = face(6,iface)

            do i=1,nbh_send_recv_multi(jj)
                call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)
                nq = nq + ngl_i * ngl_j * (5 + 2) !5 prognostics + F0 and G0
            end do

            jj = jj + 1

        end do

        !Get Neighboring Processor
        idest = nbh_proc(inbh)

        !Send to NBHs
        nreq = nreq + 1
        istart = iend + 1
        iend = istart + nq - 1

        if(nq > 0) then
          call mpi_irecv(recv_data(istart:iend), nq, &
            MPI_PRECISION,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)

          call mpi_isend(send_data(istart:iend), nq, &
            MPI_PRECISION,idest-1,99,mpi_comm_world, &
            ireq(nreq+1),ierr)
        else
          ireq(nreq) = MPI_REQUEST_NULL

          ireq(nreq+1) = MPI_REQUEST_NULL
        endif

        nreq = nreq + 1
    end do

end subroutine send_boundary_lhs

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
subroutine send_bound_dg_general(send_data,recv_data,nsize,nreq,ireq,status)

    use mod_basis, only: ngl

    use mod_face, only: face_send

    use mod_grid, only: nboun, face, mod_grid_get_face_ngl, face_type

    use mod_parallel, only: nbh_proc, num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi

    use mod_ref, only: nmessage

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    !global variables
    real, intent(in)  :: send_data(nmessage*ngl*ngl*nboun)
    real, intent(out) :: recv_data(nmessage*ngl*ngl*nboun)
    integer, intent(out) :: nreq
    integer, intent(out) :: ireq(2*num_nbh)
    integer, intent(out) :: status(mpi_status_size,2*num_nbh)
    integer, intent(in) :: nsize

    !local variables
    integer inbh, idest, istart, iend, ierr, nq
    integer ngl_i, ngl_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

    nreq = 0
    iend = 0
    jj = 1
    status=0

    do inbh = 1, num_nbh

        !Determine size
        nq = 0
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            ftype = face_type(iface)

            ilocl = face(5,iface)
            if(ftype==21) ilocl = face(6,iface)

            do i=1,nbh_send_recv_multi(jj)
                call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)
                nq = nq + ngl_i * ngl_j * nsize !NSIZE is the size of the message
            end do

            jj = jj + 1

        end do

        !Get Neighboring Processor
        idest = nbh_proc(inbh)

        !Send to NBHs
        nreq = nreq + 1
        istart = iend + 1
        iend = istart + nq - 1

        if(nq > 0) then
          call mpi_irecv(recv_data(istart:iend), nq, &
            MPI_PRECISION,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)

          call mpi_isend(send_data(istart:iend), nq, &
            MPI_PRECISION,idest-1,99,mpi_comm_world, &
            ireq(nreq+1),ierr)
        else
          ireq(nreq) = MPI_REQUEST_NULL

          ireq(nreq+1) = MPI_REQUEST_NULL
        endif

        nreq = nreq + 1
    end do

end subroutine send_bound_dg_general

subroutine send_bound_dg_general_quad(send_data,recv_data,nvarb,nreq,ireq,status)

    use mod_basis, only: nq

    use mod_face, only: face_send

    use mod_grid, only: nboun, face, mod_grid_get_face_nq, face_type

    use mod_parallel, only: nbh_proc, num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi

    use mod_ref, only: nmessage

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    !global variables
    real, intent(in)  :: send_data(nvarb*nq*nboun)
    real, intent(out) :: recv_data(nvarb*nq*nboun)
    integer, intent(out) :: nreq
    integer, intent(out) :: ireq(2*num_nbh)
    integer, intent(out) :: status(mpi_status_size,2*num_nbh)
    integer, intent(in) :: nvarb

    !local variables
    integer inbh, idest, istart, iend, ierr, nqp
    integer nq_i, nq_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

    !recv_data=0.0

    nreq = 0
    iend = 0
    jj = 1
    status=0

    do inbh = 1, num_nbh

        !Determine size
        nqp = 0
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            ftype = face_type(iface)

            ilocl = face(5,iface)
            ! if(ftype==21) ilocl = face(6,iface)

            do i=1,nbh_send_recv_multi(jj)
                call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)
                nqp = nqp + nq *  nvarb!NSIZE is the size of the message
            end do
            jj = jj + 1
        end do

        !Get Neighboring Processor
        idest = nbh_proc(inbh)

        !Send to NBHs
        nreq = nreq + 1
        istart = iend + 1
        iend = istart + nqp - 1

        if(nqp > 0) then
          call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)

          call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION,idest-1,99,mpi_comm_world, &
            ireq(nreq+1),ierr)
        else
          ireq(nreq) = MPI_REQUEST_NULL

          ireq(nreq+1) = MPI_REQUEST_NULL
        endif

        nreq = nreq + 1
    end do

    

end subroutine send_bound_dg_general_quad

subroutine send_bound_dg_general_df(send_data,recv_data,nvarb,nreq,ireq,status)

    use mod_basis, only: ngl

    use mod_face, only: face_send

    use mod_grid, only: nboun, face, mod_grid_get_face_nq, face_type

    use mod_parallel, only: nbh_proc, num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi

    use mod_ref, only: nmessage

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    !global variables
    real, intent(in)  :: send_data(nvarb*ngl*nboun)
    real, intent(out) :: recv_data(nvarb*ngl*nboun)
    integer, intent(out) :: nreq
    integer, intent(out) :: ireq(2*num_nbh)
    integer, intent(out) :: status(mpi_status_size,2*num_nbh)
    integer, intent(in) :: nvarb

    !local variables
    integer inbh, idest, istart, iend, ierr, nqp
    integer nq_i, nq_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

    !recv_data=0.0

    nreq = 0
    iend = 0
    jj = 1
    status=0

    do inbh = 1, num_nbh

        !Determine size
        nqp = 0
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            ftype = face_type(iface)

            ilocl = face(5,iface)
            ! if(ftype==21) ilocl = face(6,iface)

            do i=1,nbh_send_recv_multi(jj)
                nqp = nqp + ngl *  nvarb !NSIZE is the size of the message
            end do
            jj = jj + 1
        end do

        !Get Neighboring Processor
        idest = nbh_proc(inbh)

        !Send to NBHs
        nreq = nreq + 1
        istart = iend + 1
        iend = istart + nqp - 1

        if(nqp > 0) then
          call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)

          call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION,idest-1,99,mpi_comm_world, &
            ireq(nreq+1),ierr)
        else
          ireq(nreq) = MPI_REQUEST_NULL

          ireq(nreq+1) = MPI_REQUEST_NULL
        endif

        nreq = nreq + 1
    end do

end subroutine send_bound_dg_general_df

subroutine send_bound_dg_general_quad_layer(send_data,recv_data,nvarb,nlayers,nq,nreq,ireq,status)

    use mod_face, only: face_send

    use mod_grid, only: nboun, face, mod_grid_get_face_nq, face_type

    use mod_parallel, only: nbh_proc, num_nbh, num_send_recv, nbh_send_recv, nbh_send_recv_multi

    use mod_ref, only: nmessage

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    !global variables
    real, intent(in)  :: send_data(nvarb*nq*nboun*nlayers)
    real, intent(out) :: recv_data(nvarb*nq*nboun*nlayers)
    integer, intent(out) :: nreq
    integer, intent(out) :: ireq(2*num_nbh)
    integer, intent(out) :: status(mpi_status_size,2*num_nbh)
    integer, intent(in) :: nvarb, nlayers, nq

    !local variables
    integer inbh, idest, istart, iend, ierr, nqp
    integer nq_i, nq_j, plane_ij, jj, ilocl, ib, iface, i, ftype, iel

    !recv_data=0.0

    nreq = 0
    iend = 0
    jj = 1
    status=0

    do inbh = 1, num_nbh

        !Determine size
        nqp = 0
        do ib = 1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            ftype = face_type(iface)

            ilocl = face(5,iface)
            if(ftype==21) ilocl = face(6,iface)

            do i=1,nbh_send_recv_multi(jj)
                call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)
                nqp = nqp + nq *  nvarb*nlayers!NSIZE is the size of the message
            end do

            jj = jj + 1

        end do

        !Get Neighboring Processor
        idest = nbh_proc(inbh)

        !Send to NBHs
        nreq = nreq + 1
        istart = iend + 1
        iend = istart + nqp - 1

        if(nqp > 0) then
          call mpi_irecv(recv_data(istart:iend), nqp, &
            MPI_PRECISION,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)

          call mpi_isend(send_data(istart:iend), nqp, &
            MPI_PRECISION,idest-1,99,mpi_comm_world, &
            ireq(nreq+1),ierr)
        else
          ireq(nreq) = MPI_REQUEST_NULL

          ireq(nreq+1) = MPI_REQUEST_NULL
        endif

        nreq = nreq + 1
    end do

end subroutine send_bound_dg_general_quad_layer


