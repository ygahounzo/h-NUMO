!-----------------------------------------------------------------!
! Cubed-sphere (sphere_hex) surface mesh with LGL nodes, extruded
! into nelz vertical layers.  Ported from hex_quad / init_hex in
! swe_sphere/src/create_grid_sphere.f90 (FX Giraldo, NPS).
!
! Surface mesh:
!   6-face cubed-sphere, nel x nel elements per face, nglx x ngly
!   LGL nodes per element.  Shared edge nodes are counted once via
!   nodef(nx_s, nx_s, 6).
!
! 3D layout:
!   Vertical levels 1..nz are stacked; 3D node ip3d = (kk-1)*npoin_s + ip_s.
!   intma_table(i,j,k,ie3d) = ip3d for element ie3d at LGL (i,j,k).
!   ele_col_g(ie3d) = 2D surface element index (column).
!   index_g(1,ip3d) = surface node ip_s; index_g(2,ip3d) = vertical level kk.
!-----------------------------------------------------------------!
subroutine create_grid_sphere(coord_cg, index_g, intma_table, ele_col_g, &
    npoin_cg, nelem, &
    xglx, xgly, xglz, nglx, ngly, nglz, &
    nel, nelz, zmin, zmax, nx_s, nz)

    use mod_constants, only: earth_radius

    implicit none

    integer, intent(in)  :: npoin_cg, nelem
    integer, intent(in)  :: nglx, ngly, nglz
    integer, intent(in)  :: nel, nelz, nx_s, nz
    real,    intent(in)  :: xglx(nglx), xgly(ngly), xglz(nglz)
    real,    intent(in)  :: zmin, zmax
    real,    intent(out) :: coord_cg(3, npoin_cg)
    integer, intent(out) :: index_g(2, npoin_cg)
    integer, intent(out) :: intma_table(nglx, ngly, nglz, nelem)
    integer, intent(out) :: ele_col_g(nelem)

    integer, parameter :: nface = 6

    ! surface-mesh work arrays (allocatable to avoid large stack frames)
    real,    allocatable :: coord0(:,:,:)     ! (2, nx_s, nx_s) ref face planar coords
    real,    allocatable :: coors0(:,:,:)     ! (2, nx_s, nx_s) ref face lon/lat
    real,    allocatable :: coorsf(:,:,:,:)   ! (2, nx_s, nx_s, nface) per-face lon/lat
    real,    allocatable :: coors_s(:,:)      ! (2, npoin_s) assembled lon/lat
    integer, allocatable :: nodef(:,:,:)      ! (nx_s, nx_s, nface) surface node map
    integer, allocatable :: inode2d(:,:)      ! (nx_s, nx_s) ref-face node map
    integer, allocatable :: intma2d(:,:,:)    ! (nglx, ngly, nelem_s) 2D connectivity
    real,    allocatable :: z_lev(:)          ! (nz) z-coord at each vertical level

    integer :: npoin_s, nelem_s
    integer :: i, j, k, l, mp, ip, ie, ii, jj, kk
    integer :: ip_s, ip3d, ie2d, ie3d, ielez, ilglz, ilglx, ilgly
    integer :: i1, j1
    real    :: pi, twopi, s, dx, dy, x0, y0, x, y
    real    :: rlon, rlat, alon, alat, clon, clat
    real    :: Lz, zstart

    ! face-centre longitudes and latitudes (degrees converted below)
    real :: lon_f(nface), lat_f(nface)
    data lon_f /  0.0,  90.0, 180.0, 270.0,   0.0,   0.0 /
    data lat_f /  0.0,   0.0,   0.0,   0.0,  90.0, -90.0 /

    pi    = 4.0 * atan(1.0)
    twopi = 2.0 * pi

    do i = 1, nface
        lon_f(i) = lon_f(i) * pi / 180.0
        lat_f(i) = lat_f(i) * pi / 180.0
    end do

    ! exact unique surface node count (cubed-sphere topology)
    npoin_s = 6 * nx_s * nx_s - 12 * nx_s + 8
    nelem_s = nface * nel * nel

    allocate( coord0(2, nx_s, nx_s), coors0(2, nx_s, nx_s),            &
              coorsf(2, nx_s, nx_s, nface), coors_s(2, npoin_s),        &
              nodef(nx_s, nx_s, nface), inode2d(nx_s, nx_s),            &
              intma2d(nglx, ngly, nelem_s), z_lev(nz) )

    !-------------------------------------------------------------------
    ! Step 1 : Build reference face planar LGL grid on [-pi/4, pi/4]^2
    !-------------------------------------------------------------------
    s  = pi / 4.0
    dx = (2.0 * s) / nel
    dy = (2.0 * s) / nel

    ip = 0
    jj = 0
    do k = 1, nel
        y0 = -s + real(k - 1) * dy
        j1 = 2;  if (k == 1) j1 = 1

        do l = j1, ngly
            y = y0 + 0.5 * dy * (xgly(l) + 1.0)
            jj = jj + 1
            ii = 0

            do i = 1, nel
                x0 = -s + real(i - 1) * dx
                i1 = 2;  if (i == 1) i1 = 1

                do j = i1, nglx
                    x = x0 + 0.5 * dx * (xglx(j) + 1.0)
                    ii = ii + 1
                    ip = ip + 1
                    coord0(1, ii, jj) = x
                    coord0(2, ii, jj) = y
                    inode2d(ii, jj)   = ip
                end do
            end do
        end do
    end do

    !-------------------------------------------------------------------
    ! Step 2 : Backward equi-angular gnomonic → spherical for ref face
    !-------------------------------------------------------------------
    do j = 1, nx_s
        do i = 1, nx_s
            call gnomonic_sph(rlon, rlat, coord0(1,i,j), coord0(2,i,j))
            coors0(1, i, j) = rlon
            coors0(2, i, j) = rlat
        end do
    end do

    !-------------------------------------------------------------------
    ! Step 3 : Rotate reference face to all 6 cube faces
    !-------------------------------------------------------------------
    do mp = 1, nface
        clon = lon_f(mp)
        clat = lat_f(mp)
        do j = 1, nx_s
            do i = 1, nx_s
                rlon = coors0(1, i, j)
                rlat = coors0(2, i, j)
                alon = b_rot1_sph(rlon, rlat, clon, clat)
                alat = b_rot2_sph(rlon, rlat, clon, clat)
                if (alon < 0.0)   alon = alon + twopi
                if (alon > twopi) alon = alon - twopi
                coorsf(1, i, j, mp) = alon
                coorsf(2, i, j, mp) = alat
            end do
        end do
    end do

    !-------------------------------------------------------------------
    ! Step 4 : Assemble 6 faces into global surface node numbering
    !-------------------------------------------------------------------
    ip = 0

    ! Face 1 – all nodes are new
    do j = 1, nx_s
        do i = 1, nx_s
            ip = ip + 1
            nodef(i, j, 1)   = ip
            coors_s(1, ip)   = coorsf(1, i, j, 1)
            coors_s(2, ip)   = coorsf(2, i, j, 1)
        end do
    end do

    ! Face 2 – left column i=1 shared with right column i=nx_s of face 1
    do j = 1, nx_s
        nodef(1, j, 2) = nodef(nx_s, j, 1)
    end do
    do j = 1, nx_s
        do i = 2, nx_s
            ip = ip + 1
            nodef(i, j, 2) = ip
            coors_s(1, ip) = coorsf(1, i, j, 2)
            coors_s(2, ip) = coorsf(2, i, j, 2)
        end do
    end do

    ! Face 3 – left column shared with right column of face 2
    do j = 1, nx_s
        nodef(1, j, 3) = nodef(nx_s, j, 2)
    end do
    do j = 1, nx_s
        do i = 2, nx_s
            ip = ip + 1
            nodef(i, j, 3) = ip
            coors_s(1, ip) = coorsf(1, i, j, 3)
            coors_s(2, ip) = coorsf(2, i, j, 3)
        end do
    end do

    ! Face 4 – left shared with face 3, right wraps back to face 1
    do j = 1, nx_s
        nodef(1,    j, 4) = nodef(nx_s, j, 3)
        nodef(nx_s, j, 4) = nodef(1,    j, 1)
    end do
    do j = 1, nx_s
        do i = 2, nx_s - 1
            ip = ip + 1
            nodef(i, j, 4) = ip
            coors_s(1, ip) = coorsf(1, i, j, 4)
            coors_s(2, ip) = coorsf(2, i, j, 4)
        end do
    end do

    ! Face 5 – north pole cap; all 4 boundary rows/cols already assigned
    do i = 1, nx_s
        nodef(i,    1,    5) = nodef(i,          nx_s, 1)
        nodef(i,    nx_s, 5) = nodef(nx_s+1-i,   nx_s, 3)
    end do
    do j = 1, nx_s
        nodef(1,    j, 5) = nodef(nx_s+1-j, nx_s, 4)
        nodef(nx_s, j, 5) = nodef(j,         nx_s, 2)
    end do
    do j = 2, nx_s - 1
        do i = 2, nx_s - 1
            ip = ip + 1
            nodef(i, j, 5) = ip
            coors_s(1, ip) = coorsf(1, i, j, 5)
            coors_s(2, ip) = coorsf(2, i, j, 5)
        end do
    end do

    ! Face 6 – south pole cap; all 4 boundary rows/cols already assigned
    do i = 1, nx_s
        nodef(i,    1,    6) = nodef(nx_s+1-i, 1, 3)
        nodef(i,    nx_s, 6) = nodef(i,        1, 1)
    end do
    do j = 1, nx_s
        nodef(1,    j, 6) = nodef(j,         1, 4)
        nodef(nx_s, j, 6) = nodef(nx_s+1-j, 1, 2)
    end do
    do j = 2, nx_s - 1
        do i = 2, nx_s - 1
            ip = ip + 1
            nodef(i, j, 6) = ip
            coors_s(1, ip) = coorsf(1, i, j, 6)
            coors_s(2, ip) = coorsf(2, i, j, 6)
        end do
    end do

    if (ip /= npoin_s) then
        write(*,'(a,2i8)') 'ERROR create_grid_sphere: ip /= npoin_s: ', ip, npoin_s
        stop
    end if

    !-------------------------------------------------------------------
    ! Step 5 : Build 2D surface intma
    !-------------------------------------------------------------------
    ie = 0
    do mp = 1, nface
        do l = 1, nel
            do k = 1, nel
                ie = ie + 1
                do j = 1, ngly
                    jj = (ngly - 1) * (l - 1) + j
                    do i = 1, nglx
                        ii = (nglx - 1) * (k - 1) + i
                        intma2d(i, j, ie) = nodef(ii, jj, mp)
                    end do
                end do
            end do
        end do
    end do

    !-------------------------------------------------------------------
    ! Step 6 : Build vertical z-level array
    !-------------------------------------------------------------------
    Lz = (zmax - zmin) / nelz
    kk = 0
    do ielez = 1, nelz
        zstart = zmin + real(ielez - 1) * Lz
        i1 = 2;  if (ielez == 1) i1 = 1
        do ilglz = i1, nglz
            kk = kk + 1
            z_lev(kk) = zstart + 0.5 * Lz * (xglz(ilglz) + 1.0)
        end do
    end do

    !-------------------------------------------------------------------
    ! Step 7 : Fill 3D coord_cg and index_g
    !-------------------------------------------------------------------
    do kk = 1, nz
        do ip_s = 1, npoin_s
            ip3d = (kk - 1) * npoin_s + ip_s
            alon = coors_s(1, ip_s)
            alat = coors_s(2, ip_s)
            coord_cg(1, ip3d) = earth_radius * cos(alat) * cos(alon)
            coord_cg(2, ip3d) = earth_radius * cos(alat) * sin(alon)
            coord_cg(3, ip3d) = z_lev(kk)
            index_g(1, ip3d)  = ip_s
            index_g(2, ip3d)  = kk
        end do
    end do

    !-------------------------------------------------------------------
    ! Step 8 : Build 3D intma_table and ele_col_g
    !-------------------------------------------------------------------
    do ielez = 1, nelz
        do ie2d = 1, nelem_s
            ie3d = (ielez - 1) * nelem_s + ie2d
            ele_col_g(ie3d) = ie2d
            do ilglz = 1, nglz
                kk = (ielez - 1) * (nglz - 1) + ilglz
                do ilgly = 1, ngly
                    do ilglx = 1, nglx
                        ip_s = intma2d(ilglx, ilgly, ie2d)
                        ip3d = (kk - 1) * npoin_s + ip_s
                        intma_table(ilglx, ilgly, ilglz, ie3d) = ip3d
                    end do
                end do
            end do
        end do
    end do

    deallocate(coord0, coors0, coorsf, coors_s, nodef, inode2d, intma2d, z_lev)

contains

    !--- Backward equi-angular gnomonic: cube-face (xg,yg) → sphere (olon,olat) ---
    subroutine gnomonic_sph(olon, olat, xg, yg)
        real, intent(out) :: olon, olat
        real, intent(in)  :: xg, yg
        real :: zz
        zz   = 1.0 / sqrt(1.0 + tan(xg)**2 + tan(yg)**2)
        olon = xg
        olat = asin(zz * tan(yg))
    end subroutine gnomonic_sph

    !--- Backward longitude rotation: rotated (lon_p,lat_p) → original lon ---
    real function b_rot1_sph(alon_p, alat_p, alono, alato)
        real, intent(in) :: alon_p, alat_p, alono, alato
        real :: anum, aden
        anum = cos(alat_p) * sin(alon_p)
        aden = cos(alat_p) * cos(alon_p) * cos(alato) - sin(alat_p) * sin(alato)
        b_rot1_sph = alono + atan2(anum, aden)
    end function b_rot1_sph

    !--- Backward latitude rotation: rotated (lon_p,lat_p) → original lat ---
    real function b_rot2_sph(alon_p, alat_p, alono, alato)
        real, intent(in) :: alon_p, alat_p, alono, alato
        real :: anum
        anum = sin(alat_p) * cos(alato) + cos(alat_p) * cos(alon_p) * sin(alato)
        b_rot2_sph = asin(anum)
    end function b_rot2_sph

end subroutine create_grid_sphere
