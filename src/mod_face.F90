!----------------------------------------------------------------------!
!>@brief This module builds the Faces and Normals
!>@author  J.F. Kelly and F.X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: March 27, 2023
!       Added routines to build numormals at quadrature points
!----------------------------------------------------------------------!
module mod_face

    use mod_basis
    use mod_grid
    use mod_input

    type face_CS
        real,    dimension(:,:,:,:),  allocatable :: normal_vector, normal_vector_q
        real,    dimension(:,:,:),    allocatable :: jac_face, penalty, jac_faceq
        integer, dimension(:,:,:,:),  allocatable :: imapl, imapr, imapl_q, imapr_q
        integer, dimension(:), allocatable :: ip_bound, face_send, face_bound
        integer, dimension(:), allocatable :: face_nc
        integer nboun_max, num_bound, nface_nc
    end type face_CS

    public :: mod_face_create, mod_face_create_boundary, face_CS, mod_face_create_nc_list

    private
    !-------------------------------------------------------------

    !module variables and parameters
    

  !-----------------------------------------------------------------------

contains

    !-----------------------------------------------------------------------

    subroutine mod_face_create(G, inp, b, mf)

        use mpi

        implicit none

        type(grid), intent(in) :: G
        type(input), intent(in) :: inp
        type(basis), intent(in) :: b
        type(face_CS), intent(inout) :: mf


        integer :: AllocateStatus, ierr  ,j , iface

        if(allocated(mf%normal_vector)) deallocate(mf%normal_vector,mf%jac_face,mf%imapl,mf%imapr, mf%penalty)
    
        allocate (mf%normal_vector(3,b%ngl,b%ngl,G%nface), mf%jac_face(b%ngl,b%ngl,G%nface), &
            mf%imapl(3,b%ngl,b%ngl,G%nface), mf%imapr(3,b%ngl,b%ngl,G%nface), mf%penalty(b%ngl,b%ngl,G%nface), &
            stat=AllocateStatus)
        if(inp%is_mlswe) then
            allocate (mf%normal_vector_q(3,b%nq,b%nq,G%nface), mf%jac_faceq(b%nq,b%nq,G%nface), &
                mf%imapl_q(3,b%nq,b%nq,G%nface), mf%imapr_q(3,b%nq,b%nq,G%nface), stat=AllocateStatus)
        end if
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Face_Create **"

        !Create Normals
        call create_normals(mf%normal_vector,mf%jac_face,G%face,G%nface, b, G)
        !Create IMAP-DG for the sphere / non-conforming grid
        call create_imaplr(mf%imapl,mf%imapr,mf%normal_vector,mf%jac_face,G%face,G%nface, b, G)

        if(inp%is_mlswe) then
            call create_normals_quad(mf%normal_vector_q,mf%jac_faceq,G%face,G%nface, b, G)
            call create_imaplr_quad(mf%imapl_q,mf%imapr_q,mf%normal_vector_q,mf%jac_faceq,G%face,G%nface, b, G)
        endif

    end subroutine mod_face_create

    !-----------------------------------------------------------------------
    subroutine mod_face_create_boundary(G, mf, nboun)

        implicit none

        type(grid), intent(in) :: G
        type(face_CS), intent(inout) :: mf

        integer nboun

        if(allocated(mf%face_send)) deallocate(mf%face_send, mf%face_bound)
        allocate( mf%face_send(nboun), mf%face_bound(G%nface))
        
    end subroutine mod_face_create_boundary

    !-----------------------------------------------------------------------
    subroutine mod_face_create_nc_list(mf, G)

        implicit none
        
        type(face_CS), intent(inout) :: mf
        type(grid), intent(in) :: G

        integer :: i,j

        j = 0

        mf%nface_nc = 0
        do i=1,G%nface
            if(G%face(9,i)>0 .or. G%face_type(i)==12) mf%nface_nc = mf%nface_nc + 1
        end do

        if(allocated(mf%face_nc)) deallocate(mf%face_nc)
        allocate(mf%face_nc(mf%nface_nc))
        mf%face_nc=0

        if(mf%nface_nc>0) then
            do i=1,G%nface
                if(G%face(9,i)>0 .or. G%face_type(i)==12) then
                    j=j+1
                    mf%face_nc(j)=i
                end if
            end do

            if (mf%nface_nc.ne.j) then
                print*,"mod_face_create_nc_list: number of non-conforming faces mismatch"
                stop
            end if
        end if

    end subroutine mod_face_create_nc_list

end module mod_face
