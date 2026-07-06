!-----------------------------------------------------------------!
! Builds the icosahedral quadrilateral base mesh at an arbitrary
! resolution "nico" (elements per icosahedron-edge, same role as
! nelx for sphere_hex — NOT restricted to a power of 2) and writes
! it to an ABAQUS-format .inp file so that p4est can read it via
! p4est_connectivity_read_inp() (read_external_grid_flg path) and
! derive tree_to_tree/tree_to_face connectivity itself via
! p4est_connectivity_complete().
!
! Construction (building only the coarse
! linear mesh — p4est's existing pipeline does the ngl x ngl
! GLL promotion and the final radial sphere projection):
!   1. The classic 12-vertex, 20-triangle icosahedron.
!   2. Each triangle is uniformly subdivided into nico^2 smaller
!      triangles via integer barycentric coordinates (a,b,c),
!      a+b+c=nico — a linear count, not a quadtree level, so nico
!      can be any positive integer.  Points on a shared icosahedron
!      edge are looked up/created once and reused by both adjacent
!      triangles (get_edge_lattice_point), keyed on the edge's two
!      endpoint vertices plus a canonical (direction-independent)
!      position along it.
!   3. Every small triangle is split into 3 quads (corner ->
!      edge-midpoint -> triangle centroid -> other edge-midpoint);
!      edge midpoints shared between two small triangles (whether
!      within one big triangle or across two) are looked up/created
!      once via get_edge_midpoint, keyed on the two endpoint vertices.
!
! Vertices are left un-normalized (chord positions); mod_p4est.F90
! radially projects every DOF onto the sphere afterward, the same
! way it does for the cubed-sphere's flat cube-corner vertices.
!-----------------------------------------------------------------!
subroutine create_grid_sphere_ico_inp(filename, nico)

    implicit none

    character(len=*), intent(in) :: filename
    integer,           intent(in) :: nico

    integer, parameter :: nvert0 = 12
    integer, parameter :: nface0 = 20

    integer :: nvert_max, nquad, nedge1_max, nedge2_max

    real,    allocatable :: vert(:,:)
    integer :: face(3, nface0)
    integer, allocatable :: quad(:,:)
    integer, allocatable :: lat(:,:,:)          ! lat(a,b,iface) -> global vertex id (c=nico-a-b implicit)

    ! canonical-edge dedup table for the coarse per-face lattice
    integer, allocatable :: e1_lo(:), e1_hi(:), e1_k(:), e1_mid(:)
    integer :: ne1

    ! canonical-edge dedup table for small-triangle edge midpoints
    integer, allocatable :: e2_v1(:), e2_v2(:), e2_mid(:)
    integer :: ne2

    integer :: nvert, iquad
    integer :: iface, v1, v2, v3, a, b, c, p1, p2, p3
    real    :: t, cx, cy, cz, nx, ny, nz
    integer :: fio, i

    nquad      = 60 * nico * nico
    nvert_max  = 12 + 30*max(nico-1,0) + 10*max(nico-1,0)*max(nico-2,0) + 20*nico*nico + 30*nico*nico + 8
    nedge1_max = 30*nico + 8
    nedge2_max = 30*nico*nico + 8

    allocate(vert(3, nvert_max))
    allocate(quad(4, nquad))
    allocate(lat(0:nico, 0:nico, nface0))
    allocate(e1_lo(nedge1_max), e1_hi(nedge1_max), e1_k(nedge1_max), e1_mid(nedge1_max))
    allocate(e2_v1(nedge2_max), e2_v2(nedge2_max), e2_mid(nedge2_max))

    !----------------------------------------------------------------
    ! 12 base icosahedron vertices (golden-ratio construction)
    !----------------------------------------------------------------
    t = (1.0 + sqrt(5.0)) / 2.0

    vert(:,1)  = (/ -1.0,  t,  0.0 /)
    vert(:,2)  = (/  1.0,  t,  0.0 /)
    vert(:,3)  = (/ -1.0, -t,  0.0 /)
    vert(:,4)  = (/  1.0, -t,  0.0 /)
    vert(:,5)  = (/  0.0, -1.0,  t /)
    vert(:,6)  = (/  0.0,  1.0,  t /)
    vert(:,7)  = (/  0.0, -1.0, -t /)
    vert(:,8)  = (/  0.0,  1.0, -t /)
    vert(:,9)  = (/  t,  0.0, -1.0 /)
    vert(:,10) = (/  t,  0.0,  1.0 /)
    vert(:,11) = (/ -t,  0.0, -1.0 /)
    vert(:,12) = (/ -t,  0.0,  1.0 /)
    nvert = nvert0

    !----------------------------------------------------------------
    ! 20 triangular faces (1-indexed; standard icosahedron connectivity)
    !----------------------------------------------------------------
    face(:,1)  = (/ 1, 12,  6 /);  face(:,2)  = (/ 1,  6,  2 /)
    face(:,3)  = (/ 1,  2,  8 /);  face(:,4)  = (/ 1,  8, 11 /)
    face(:,5)  = (/ 1, 11, 12 /);  face(:,6)  = (/ 2,  6, 10 /)
    face(:,7)  = (/ 6, 12,  5 /);  face(:,8)  = (/ 12, 11,  3 /)
    face(:,9)  = (/ 11,  8,  7 /); face(:,10) = (/ 8,  2,  9 /)
    face(:,11) = (/ 4, 10,  5 /);  face(:,12) = (/ 4,  5,  3 /)
    face(:,13) = (/ 4,  3,  7 /);  face(:,14) = (/ 4,  7,  9 /)
    face(:,15) = (/ 4,  9, 10 /);  face(:,16) = (/ 5, 10,  6 /)
    face(:,17) = (/ 3,  5, 12 /);  face(:,18) = (/ 7,  3, 11 /)
    face(:,19) = (/ 9,  7,  8 /);  face(:,20) = (/ 10, 9,  2 /)

    ! Ensure every face normal points outward (away from the origin);
    ! flip winding (swap v2,v3) if it doesn't.
    do iface = 1, nface0
        v1 = face(1,iface); v2 = face(2,iface); v3 = face(3,iface)
        call outward_normal(vert(:,v1), vert(:,v2), vert(:,v3), nx, ny, nz)
        cx = (vert(1,v1)+vert(1,v2)+vert(1,v3))/3.0
        cy = (vert(2,v1)+vert(2,v2)+vert(2,v3))/3.0
        cz = (vert(3,v1)+vert(3,v2)+vert(3,v3))/3.0
        if (nx*cx + ny*cy + nz*cz < 0.0) then
            face(2,iface) = v3
            face(3,iface) = v2
        end if
    end do

    !----------------------------------------------------------------
    ! Build the per-face (a,b,c) barycentric lattice, a+b+c=nico.
    ! Corners: (nico,0,0)->v1, (0,nico,0)->v2, (0,0,nico)->v3.
    ! Edge points (exactly one of a,b,c is 0) are shared with the
    ! neighboring face across that icosahedron edge, so they go
    ! through get_edge_lattice_point; strictly-interior points
    ! (a,b,c all > 0) are unique to this face.
    !----------------------------------------------------------------
    ne1 = 0
    lat = 0

    do iface = 1, nface0
        v1 = face(1,iface); v2 = face(2,iface); v3 = face(3,iface)

        do a = 0, nico
            do b = 0, nico - a
                c = nico - a - b

                if (a == nico) then
                    lat(a,b,iface) = v1                                  ! (nico,0,0)
                else if (b == nico) then
                    lat(a,b,iface) = v2                                  ! (0,nico,0)
                else if (c == nico) then
                    lat(a,b,iface) = v3                                  ! (0,0,nico)
                else if (c == 0) then
                    ! edge v1(a=nico) -- v2(b=nico), local param b=0..nico
                    call get_edge_lattice_point(v1, v2, b, nico, vert, nvert_max, nvert, &
                        e1_lo, e1_hi, e1_k, e1_mid, nedge1_max, ne1, lat(a,b,iface))
                else if (b == 0) then
                    ! edge v1(a=nico) -- v3(c=nico), local param c=0..nico
                    call get_edge_lattice_point(v1, v3, c, nico, vert, nvert_max, nvert, &
                        e1_lo, e1_hi, e1_k, e1_mid, nedge1_max, ne1, lat(a,b,iface))
                else if (a == 0) then
                    ! edge v2(b=nico) -- v3(c=nico), local param c=0..nico
                    call get_edge_lattice_point(v2, v3, c, nico, vert, nvert_max, nvert, &
                        e1_lo, e1_hi, e1_k, e1_mid, nedge1_max, ne1, lat(a,b,iface))
                else
                    ! strictly interior: always new
                    nvert = nvert + 1
                    vert(1,nvert) = (a*vert(1,v1) + b*vert(1,v2) + c*vert(1,v3)) / real(nico)
                    vert(2,nvert) = (a*vert(2,v1) + b*vert(2,v2) + c*vert(2,v3)) / real(nico)
                    vert(3,nvert) = (a*vert(3,v1) + b*vert(3,v2) + c*vert(3,v3)) / real(nico)
                    lat(a,b,iface) = nvert
                end if
            end do
        end do
    end do

    !----------------------------------------------------------------
    ! Split every small triangle (upward and downward) of every
    ! face's lattice into 3 quads.
    !----------------------------------------------------------------
    ne2 = 0
    iquad = 0

    do iface = 1, nface0

        ! "Upward" small triangles: barycentric loop index (a,b),
        ! a+b<=nico-1; corners (a+1,b,c),(a,b+1,c),(a,b,c+1).
        do a = 0, nico-1
            do b = 0, nico-1-a
                p1 = lat(a+1, b,   iface)
                p2 = lat(a,   b+1, iface)
                p3 = lat(a,   b,   iface)
                call split_triangle_into_quads(p1, p2, p3, vert, nvert_max, nvert, &
                    e2_v1, e2_v2, e2_mid, nedge2_max, ne2, quad, nquad, iquad)
            end do
        end do

        ! "Downward" small triangles: loop index (a,b), a+b<=nico-2;
        ! corners (a+1,b+1,c),(a+1,b,c+1),(a,b+1,c+1).
        do a = 0, nico-2
            do b = 0, nico-2-a
                p1 = lat(a+1, b+1, iface)
                p2 = lat(a+1, b,   iface)
                p3 = lat(a,   b+1, iface)
                call split_triangle_into_quads(p1, p2, p3, vert, nvert_max, nvert, &
                    e2_v1, e2_v2, e2_mid, nedge2_max, ne2, quad, nquad, iquad)
            end do
        end do

    end do

    if (iquad /= nquad) then
        print *, "ERROR create_grid_sphere_ico_inp: built ", iquad, " quads, expected ", nquad
        stop 1
    end if

    !----------------------------------------------------------------
    ! Write ABAQUS .inp file (p4est_connectivity_read_inp format)
    !----------------------------------------------------------------
    open(newunit=fio, file=trim(filename), status='replace', action='write')
    write(fio,'(A)') '*Heading'
    write(fio,'(A)') ' icosahedral_base_mesh.inp'
    write(fio,'(A)') '*Node'
    do i = 1, nvert
        write(fio,'(I0,3(",",1x,E23.15))') i, vert(1,i), vert(2,i), vert(3,i)
    end do
    write(fio,'(A)') '*Element, type=CPS4, ELSET=Surface1'
    do i = 1, nquad
        write(fio,'(I0,4(",",1x,I0))') i, quad(1,i), quad(2,i), quad(3,i), quad(4,i)
    end do
    close(fio)

contains

    subroutine outward_normal(a, b, c, nx, ny, nz)
        real, intent(in)  :: a(3), b(3), c(3)
        real, intent(out) :: nx, ny, nz
        real :: e1(3), e2(3)
        e1 = b - a
        e2 = c - a
        nx = e1(2)*e2(3) - e1(3)*e2(2)
        ny = e1(3)*e2(1) - e1(1)*e2(3)
        nz = e1(1)*e2(2) - e1(2)*e2(1)
    end subroutine outward_normal

end subroutine create_grid_sphere_ico_inp

!-----------------------------------------------------------------!
! Returns the global vertex id of the point a fraction tA/nico of
! the way from vA to vB along a coarse icosahedron edge.  tA=0 or
! tA=nico return the endpoints directly.  Otherwise the point is
! looked up (or created, on first visit) in a table keyed on the
! edge's two endpoints plus a canonical position measured from the
! lower-indexed endpoint — so the neighboring face, which walks the
! same edge from the opposite end, resolves to the same vertex.
!-----------------------------------------------------------------!
subroutine get_edge_lattice_point(vA, vB, tA, nico, vert, nvert_max, nvert, &
                                   e_lo, e_hi, e_k, e_mid, nedge_max, nedge, id)

    implicit none

    integer, intent(in)    :: vA, vB, tA, nico, nvert_max, nedge_max
    real,    intent(inout) :: vert(3, nvert_max)
    integer, intent(inout) :: nvert
    integer, intent(inout) :: e_lo(nedge_max), e_hi(nedge_max), e_k(nedge_max), e_mid(nedge_max)
    integer, intent(inout) :: nedge
    integer, intent(out)   :: id

    integer :: vlo, vhi, kcanon, m
    real    :: frac

    if (tA == 0) then
        id = vA; return
    else if (tA == nico) then
        id = vB; return
    end if

    if (vA < vB) then
        vlo = vA; vhi = vB; kcanon = tA
    else
        vlo = vB; vhi = vA; kcanon = nico - tA
    end if

    do m = 1, nedge
        if (e_lo(m) == vlo .and. e_hi(m) == vhi .and. e_k(m) == kcanon) then
            id = e_mid(m)
            return
        end if
    end do

    nvert = nvert + 1
    id = nvert
    frac = real(kcanon) / real(nico)
    vert(:,id) = (1.0 - frac)*vert(:,vlo) + frac*vert(:,vhi)

    nedge = nedge + 1
    e_lo(nedge) = vlo; e_hi(nedge) = vhi; e_k(nedge) = kcanon; e_mid(nedge) = id

end subroutine get_edge_lattice_point

!-----------------------------------------------------------------!
! Looks up the midpoint vertex already created for edge (i,j) (in
! either order), or creates a new one if this is the edge's first
! visit.  Edge identity is keyed on the unordered pair of endpoint
! indices — used both for small-triangle edges internal to one
! coarse face and for small-triangle edges lying on a coarse-face
! boundary (shared with the neighboring face's own subdivision).
!-----------------------------------------------------------------!
subroutine get_edge_midpoint2(i, j, vert, nvert_max, nvert, &
                               e_v1, e_v2, e_mid, nedge_max, nedge, mid)

    implicit none

    integer, intent(in)    :: i, j, nvert_max, nedge_max
    real,    intent(inout) :: vert(3, nvert_max)
    integer, intent(inout) :: nvert
    integer, intent(inout) :: e_v1(nedge_max), e_v2(nedge_max), e_mid(nedge_max)
    integer, intent(inout) :: nedge
    integer, intent(out)   :: mid

    integer :: lo, hi, k

    lo = min(i,j); hi = max(i,j)

    do k = 1, nedge
        if (e_v1(k) == lo .and. e_v2(k) == hi) then
            mid = e_mid(k)
            return
        end if
    end do

    nvert = nvert + 1
    mid = nvert
    vert(:,mid) = 0.5 * (vert(:,i) + vert(:,j))

    nedge = nedge + 1
    e_v1(nedge) = lo
    e_v2(nedge) = hi
    e_mid(nedge) = mid

end subroutine get_edge_midpoint2

!-----------------------------------------------------------------!
! Splits one small triangle (corners p1,p2,p3, already-existing
! global vertex ids) into 3 quads: corner -> edge-midpoint ->
! centroid -> other edge-midpoint.  Appends the 3 quads to quad(:,:)
! starting at iquad+1 and advances iquad.
!-----------------------------------------------------------------!
subroutine split_triangle_into_quads(p1, p2, p3, vert, nvert_max, nvert, &
                                      e_v1, e_v2, e_mid, nedge_max, nedge, quad, nquad, iquad)

    implicit none

    integer, intent(in)    :: p1, p2, p3, nvert_max, nedge_max, nquad
    real,    intent(inout) :: vert(3, nvert_max)
    integer, intent(inout) :: nvert
    integer, intent(inout) :: e_v1(nedge_max), e_v2(nedge_max), e_mid(nedge_max)
    integer, intent(inout) :: nedge
    integer, intent(inout) :: quad(4, nquad)
    integer, intent(inout) :: iquad

    integer :: m12, m23, m31, icent

    call get_edge_midpoint2(p1, p2, vert, nvert_max, nvert, e_v1, e_v2, e_mid, nedge_max, nedge, m12)
    call get_edge_midpoint2(p2, p3, vert, nvert_max, nvert, e_v1, e_v2, e_mid, nedge_max, nedge, m23)
    call get_edge_midpoint2(p3, p1, vert, nvert_max, nvert, e_v1, e_v2, e_mid, nedge_max, nedge, m31)

    nvert = nvert + 1
    icent = nvert
    vert(1,icent) = (vert(1,p1)+vert(1,p2)+vert(1,p3))/3.0
    vert(2,icent) = (vert(2,p1)+vert(2,p2)+vert(2,p3))/3.0
    vert(3,icent) = (vert(3,p1)+vert(3,p2)+vert(3,p3))/3.0

    iquad = iquad + 1
    quad(:,iquad) = (/ p1, m12, icent, m31 /)
    iquad = iquad + 1
    quad(:,iquad) = (/ p2, m23, icent, m12 /)
    iquad = iquad + 1
    quad(:,iquad) = (/ p3, m31, icent, m23 /)

end subroutine split_triangle_into_quads
