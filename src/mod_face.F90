!----------------------------------------------------------------------!
!>@brief This module builds the Faces and Normals
!>@author  J.F. Kelly and F.X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!----------------------------------------------------------------------!
module mod_face

    use mod_basis, only: ngl, nglx, ngly, nglz, nq, nqx, nqy, nqz

    use mod_grid, only: face, nface

    use mod_input, only: is_mlswe

    public :: mod_face_create, mod_face_create_boundary, &
        normal_vector, jac_face, penalty, &
        imapl, imapr, ip_bound, num_bound, face_send, &
        face_nc, nface_nc, mod_face_create_nc_list, &
        normal_vector_q, jac_faceq, imapl_q,imapr_q

    private
    !-------------------------------------------------------------

    !module variables and parameters
    real,    dimension(:,:,:,:),  allocatable :: normal_vector, normal_vector_q
    real,    dimension(:,:,:),    allocatable :: jac_face, penalty, jac_faceq
    integer, dimension(:,:,:,:),  allocatable :: imapl, imapr, imapl_q, imapr_q
    integer, dimension(:), allocatable :: ip_bound, face_send, face_bound
    integer, dimension(:), allocatable :: face_nc
    integer nboun_max, num_bound, nface_nc

  !-----------------------------------------------------------------------

contains

    !-----------------------------------------------------------------------

    subroutine mod_face_create()

        use mpi

        !use mod_input, only: is_mlswe

        implicit none

        integer :: AllocateStatus, ierr  ,j , iface

        if(allocated(normal_vector)) deallocate(normal_vector,jac_face,imapl,imapr, penalty)
    
        allocate (normal_vector(3,ngl,ngl,nface), jac_face(ngl,ngl,nface), &
            imapl(3,ngl,ngl,nface), imapr(3,ngl,ngl,nface), penalty(ngl,ngl,nface),  stat=AllocateStatus)
        if(is_mlswe) then
            allocate (normal_vector_q(3,nq,nq,nface), jac_faceq(nq,nq,nface), imapl_q(3,nq,nq,nface), imapr_q(3,nq,nq,nface), stat=AllocateStatus)
        end if
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Face_Create **"

        !Create Normals
        call create_normals(normal_vector,jac_face,face,nface)

        if(is_mlswe) then
            call create_normals_quad(normal_vector_q,jac_faceq,face,nface)
            call create_imaplr_quad(imapl_q,imapr_q,normal_vector_q,jac_faceq,face,nface)
            !call create_imaplr_quad1(imapl_q,imapr_q,face,nface)
        endif

        !Create IMAP-DG for the sphere / non-conforming grid
        call create_imaplr(imapl,imapr,normal_vector,jac_face,face,nface)

    end subroutine mod_face_create

    !-----------------------------------------------------------------------
    subroutine mod_face_create_boundary(nboun)

        use mod_input,only: space_method

        implicit none

        integer nboun

        if (space_method /= 'dg') then
            if(allocated(ip_bound)) deallocate(ip_bound)
            allocate(ip_bound(nboun*ngl*ngl))
            call create_boundary_cg(ip_bound,num_bound,0,nboun*ngl*ngl)
        else
            if(allocated(face_send)) deallocate(face_send, face_bound)
            allocate( face_send(nboun), face_bound(nface))
        end if

        !Create jump penalty weight
!!        call compute_jump_penalty_weight()  !FXG: is making AMR hang up.
        
    end subroutine mod_face_create_boundary

    !-----------------------------------------------------------------------
    subroutine mod_face_create_nc_list()

        use mod_grid, only: face, nface, face_type

        implicit none

        integer :: i,j

        j=0

        nface_nc=0
        do i=1,nface
            if(face(9,i)>0 .or. face_type(i)==12) nface_nc = nface_nc + 1
        end do

        if(allocated(face_nc)) deallocate(face_nc)
        allocate(face_nc(nface_nc))
        face_nc=0

        if(nface_nc>0) then
            do i=1,nface
                if(face(9,i)>0 .or. face_type(i)==12) then
                    j=j+1
                    face_nc(j)=i
                end if
            end do

            if (nface_nc.ne.j) then
                print*,"mod_face_create_nc_list: number of non-conforming faces mismatch"
                stop
            end if
        end if

    end subroutine mod_face_create_nc_list

end module mod_face
