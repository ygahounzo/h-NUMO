!----------------------------------------------------------------------!
!> @brief This module builds the Faces and Normals
!> @author J.F. Kelly on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!> @detailss Constructs a face data structure: face(nface,8), where
!> nface = 3*nelem + nelex*neley + nelex*nelez + neley*nelez
!> face(1:4,iface) = 4 vertices of the face
!> face(5,iface) = local face relative to left element
!> face(6,iface) = local face relative to right element
!> face(7,iface) = element to the left
!> face(8,iface) = element to the right or Boundary Condition
!> In the fourth dimension, return i, j, or k variable 
!> @date 17 September 2009 J.F. Kelly
!> @date 11 October 2009 J.F. Kelly
!> @date 17 November 2009 J.F. Kelly
!> @date 18 November 2009 J.F. Kelly
!> @date 3 August 2010
!----------------------------------------------------------------------!
subroutine create_face(face,nface, FACE_LEN)

    use mod_basis, only: nglx, ngly, nglz

    use mod_global_grid, only: intma_g, bsido_g, nelem_g, nboun_g

    use mod_input, only: decomp_type

    use mod_mpi_utilities, only: irank

    implicit none

    !global arrays
    integer FACE_LEN
    integer, intent(in)  :: nface
    integer, intent(out) :: face(FACE_LEN,nface)

    ! local arrays: these are the local verticies of each hexahedral element
    integer in, ip, ie
    integer j1, k1, k2, k3, i1, i2, j2
    integer ips(4)
    integer node_loc(8,3)
    integer face_loc(6,4)
    integer isnew, iface, iboun, temp
 
    logical issame

    !Create 8 local vertices for each hexahedron
    node_loc(1,1) = 1
    node_loc(1,2) = 1
    node_loc(1,3) = 1

    node_loc(2,1) = nglx
    node_loc(2,2) = 1
    node_loc(2,3) = 1

    node_loc(3,1) = 1
    node_loc(3,2) = ngly
    node_loc(3,3) = 1

    node_loc(4,1) = nglx
    node_loc(4,2) = ngly
    node_loc(4,3) = 1

    node_loc(5,1) = 1
    node_loc(5,2) = 1
    node_loc(5,3) = nglz

    node_loc(6,1) = nglx
    node_loc(6,2) = 1
    node_loc(6,3) = nglz

    node_loc(7,1) = 1
    node_loc(7,2) = ngly
    node_loc(7,3) = nglz

    node_loc(8,1) = nglx
    node_loc(8,2) = ngly
    node_loc(8,3) = nglz

    ! Assign local vertices to local faces
    ! Face 1: z = -1 : xy
    face_loc(1,1) = 1
    face_loc(1,2) = 3
    face_loc(1,3) = 4
    face_loc(1,4) = 2

    ! Face 2: z = +1 : xy
    face_loc(2,1) = 5
    face_loc(2,2) = 6
    face_loc(2,3) = 8
    face_loc(2,4) = 7

    ! Face 3: y = -1 : xz
    face_loc(3,1) = 1
    face_loc(3,2) = 2
    face_loc(3,3) = 6
    face_loc(3,4) = 5

    ! Face 4: y = +1 : xz
    face_loc(4,1) = 3
    face_loc(4,2) = 7
    face_loc(4,3) = 8
    face_loc(4,4) = 4

    ! Face 5: x = -1 : yz
    face_loc(5,1) = 1
    face_loc(5,2) = 5
    face_loc(5,3) = 7
    face_loc(5,4) = 3

    ! Face 6 : x =+1 : yz
    face_loc(6,1) = 2
    face_loc(6,2) = 4
    face_loc(6,3) = 8
    face_loc(6,4) = 6

    ! initialize the face data structure
    face = 0
    iface = 0

    ! STEP I: Create the Faces

    !Create faces for first element
    ie = 1
    do i1 = 1,6
        if(nglx == 1 .and. (i1 == 5 .or. i1 == 6)) cycle
        if(ngly == 1 .and. (i1 == 3 .or. i1 == 4)) cycle
        if(nglz == 1 .and. (i1 == 1 .or. i1 == 2)) cycle
        do i2 = 1,4
            !Get global nodes for each face
            j1 = face_loc(i1,i2)
            k1 = node_loc(j1,1)
            k2 = node_loc(j1,2)
            k3 = node_loc(j1,3)
            ips(i2) = intma_g(k1,k2,k3,ie)
        end do
        iface = iface + 1
        face(1:4,iface) = ips
        face(5,iface) = i1
        face(7,iface) = ie
    end do

    do ie = 2,nelem_g
        do i1 = 1,6
            if(nglx == 1 .and. (i1 == 5 .or. i1 == 6)) cycle
            if(ngly == 1 .and. (i1 == 3 .or. i1 == 4)) cycle
            if(nglz == 1 .and. (i1 == 1 .or. i1 == 2)) cycle
            do i2 = 1,4
                !Get global nodes for each face
                j1 = face_loc(i1,i2)
                k1 = node_loc(j1,1)
                k2 = node_loc(j1,2)
                k3 = node_loc(j1,3)
                ips(i2) = intma_g(k1,k2,k3,ie)
            end do
            isnew = 1
            do j2 = 1,iface
                if (issame(face(1:4,j2),ips) ) then
                    if(ie == face(7,j2)) cycle
                    face(6,j2) = i1
                    face(8,j2) = ie
                    isnew = 0
                    exit
                end if
            end do

            if (isnew == 1) then
                !print *, "Creating Face ", iface
                iface = iface + 1
                face(1:4,iface) = ips
                face(5,iface) = i1
                face(7,iface) = ie
            end if
        end do
    end do

    if (iface /= nface) then
        print*,' Error in CREATE_FACE'
        print*,' irank iface nface = ',irank,iface,nface
        stop
    end if

    ! STEP II: Use BSIDO to Assign Info to Faces on the boundary
    do iboun = 1,nboun_g
        if (bsido_g(6,iboun) /= 3 .or. decomp_type(1:4) == 'geom') then
            !FXG: Skip Periodic Faces since they are treated as regular faces now
            do iface = 1,nface
                if (issame(bsido_g(1:4,iboun),face(1:4,iface))) then
                    face(8,iface) = -bsido_g(6,iboun)
                end if
            end do
        end if
    end do

    face(1:4,:) = 0

end subroutine create_face

!----------------------------------------------------------------------!
!>@brief This function returns one if ip1 and ip2 correspond to the same face
!> Returns zero otherwise
!>@todo This is brute force and not elegant...there is a better way
!----------------------------------------------------------------------!
logical function issame(ip1, ip2)

    integer, intent(in) :: ip1(4), ip2(4)
    integer :: i1, i2
    logical :: found

    issame = .true.
  
    do i1 = 1,4
        !vertex from ip1 exists in ip2?
        found = .false.
        do i2 = 1,4
            if (ip1(i1) == ip2(i2)) then
                found = .true.
                exit
            end if
        end do
        if(.not. found) then
            issame = .false.
            return
        endif
     
        !vertex from ip2 exists in ip1?
        found = .false.
        do i2 = 1,4
            if (ip2(i1) == ip1(i2)) then
                found = .true.
                exit
            end if
        end do
        if(.not. found) then
            issame = .false.
            return
        endif
    end do

end function issame
