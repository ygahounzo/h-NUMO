!-----------------------------------------------------------------!
!>@brief Creates a 3D spherical shell grid using spectral elements
!>@details Output:
!> coord_cg(3,npoin_cg) = (x,y,z) coords
!> indexg(2,npoin_cg) = horizontal, vertical indices for physics
!> intma(ngl,ngl,ngl,nelem) = interconnectivity matrix (standard 3D Format)
!> intma_lev(ngl,ngl,nelem_s,npoin_r) = horizontal level-based interconnectivity matrix
!> intma_vertical(npoin_s,ngl,nelem_r) = vertical level-based interconnectivity matrix
!> ele_col(nelem) = level-element that owns global element
!> bsido(6,nboun) = boundary matrix
!> Input: 
!> nelem_r = number of elements in the radial direction
!> nelem_s = number of elements in each concentric sphere
!> ngl = number of LGL points
!> nelem = total number of elements = nelem_r*nelem_s
!> npoin_cg = total number of points
!> nboun = number of boundary elements = 2*nelem_s
!> npoin_r = number of grid points in the radial direction = nop*nelem_r + 1
!> npoin_s = number of grid points on each concentric sphere
!> rmin = radius of the inner (ground) sphere
!> rmax = radius of the outer (top) sphere
!>@author James F. Kelly 26 April 2010 
!>@date 9 May 2010
!>@date 29 June 2010
!-----------------------------------------------------------------!
subroutine create_grid_sphere_hex(coord_cg,indexg,intma,intma_lev,intma_s,ele_col,bsido,&
    iboun,nelem_r,nelem_s,nel,nglx,ngly,nglz,nelem,npoin_cg,nboun,npoin_r,npoin_s,rmin,rmax)
  
    implicit none

    !Global variables
    real coord_cg(3,npoin_cg)
    integer indexg(2,npoin_cg)
    integer intma(nglx,ngly,nglz,nelem)
    integer intma_lev(nglx,ngly,nelem_s,npoin_r)
    integer ele_col(nelem)
    integer bsido(6,nboun), iboun(2)
    integer nboun, npoin_r, npoin_s,nel
    integer nglx,ngly,nglz, nelem,nelem_r,nelem_s,npoin_cg
    real rmin, rmax
  
    integer ier, ipr, ies, ips, ie, ip, ib
    integer i, j, k
    real dr
    real coord_s(3,npoin_s)
    integer intma_s(nglx,ngly,nelem_s)
    real r(npoin_r)
    integer intma_r(nglz,nelem_r)
    integer nface, nelem0, nx, ny, npoin0
    real r2
  
    nface=6
    nx=nel*(nglx-1)+1
    ny=nel*(ngly-1)+1
    npoin0=nx*ny
    nelem0 = nel*nel
  
    !Create 1D Spectral Element Grid in R
    call create1d_grid(r,intma_r,nglz,nelem_r,npoin_r,rmin,rmax)
  
    !Create 2D Spherical Grid: assumes nely=nelx
    call hex_quad(coord_s,intma_s,nel,nglx,nface,nelem_s,npoin_s,nelem0,npoin0,nx,ny)
  
    !Construct Cartesian Product [rmin,rmax] x S^2 to create
    !a spherical shell
    do ier = 1,nelem_r
        do k = 1,nglz
        
            ipr = intma_r(k,ier)
            do ies = 1,nelem_s
                do i = 1,nglx
                    do j = 1,ngly
                        ips = intma_s(i,j,ies)
                        ie = (ier - 1)*nelem_s + ies
                        ip = (ipr - 1)*npoin_s + ips
                        intma(i,j,k,ie) = ip                !3D Intma
                        intma_lev(i,j,ies,ipr) = ip         !Level-Based Intma
                        ele_col(ie) = ies                   !Level-Based element that belongs to ie
                        coord_cg(1,ip) = r(ipr)*coord_s(1,ips)
                        coord_cg(2,ip) = r(ipr)*coord_s(2,ips)
                        coord_cg(3,ip) = r(ipr)*coord_s(3,ips)
                        indexg(1,ip) = ips
                        indexg(2,ip) = ipr
                    end do
                end do
            end do
        end do
    end do
    print *, "Coords Created"
  
    !Create Bsido
    !NFBC on the Ground
    if(nglz > 1) then
        ier = 1
        ib = 1
        do ies =1,nelem_s
            ie = (ier - 1)*nelem_s + ies
            k = 1
     
            i = 1; j = 1
            ip = intma(i,j,k,ie)
            bsido(1,ib) = ip
     
            i = 1; j = ngly
            ip = intma(i,j,k,ie)
            bsido(2,ib) = ip
     
            i = nglx; j = ngly
            ip = intma(i,j,k,ie)
            bsido(3,ib) = ip
     
            i = nglx; j = 1
            ip = intma(i,j,k,ie)
            bsido(4,ib) = ip
            bsido(5,ib) = ie
            bsido(6,ib) = iboun(1)
     
            ib = ib + 1
        end do
  
        !NRBC at the top of the atmosphere
        ier = nelem_r
        do ies =1,nelem_s
            ie = (ier - 1)*nelem_s + ies
            k = nglz
     
            i = 1; j = 1
            ip = intma(i,j,k,ie)
            bsido(1,ib) = ip
     
            i = nglx; j = 1
            ip = intma(i,j,k,ie)
            bsido(2,ib) = ip
     
            i = nglx; j = ngly
            ip = intma(i,j,k,ie)
            bsido(3,ib) = ip
     
            i = 1; j = ngly
            ip = intma(i,j,k,ie)
            bsido(4,ib) = ip
            bsido(5,ib) = ie
            bsido(6,ib) = iboun(2)
     
            ib = ib + 1
     
        end do
    endif
end subroutine create_grid_sphere_hex

!-----------------------------------------------------------------!
subroutine create1d_grid(r,intmar,ngl,nelr,npoinr,rmin,rmax)
 
    use mod_input, only: icase
    use mod_legendre, only : legendre_gauss_lobatto
 
    implicit none
    integer ngl, nelr, npoinr
    real r(npoinr)
    integer intmar(ngl,nelr)
    real xgl(ngl), wgl(ngl)
    real z(nelr+1), dr(nelr), rmin, rmax, num, den
    integer ier, ipr, i
  
    z = 0

    do i = 2,nelr+1
        if(icase == 13) then
            num = (i-1)*(ngl-1)+1
            den = nelr*(ngl-1)+1
            z(i) = (rmax-rmin)*(sqrt(15.0*(num/den)**2+1)-1)/3.0
            dr(i-1) = z(i)-z(i-1)
        else
            dr(i-1) = (rmax-rmin)/nelr
        endif

    end do !i
     
    call legendre_gauss_lobatto(ngl,xgl,wgl)
    ier = 1
    ipr = 0
    do i = 1,ngl
        ipr = ipr + 1
        r(ipr) = rmin +  dr(ier)*(xgl(i) + 1)/2.0
        intmar(i,ier) = ipr
    end do

    z(1) = rmin
    do ier = 2,nelr
        z(ier) = z(ier-1) + dr(ier)
        intmar(1,ier) = ipr
        do i = 2,ngl
            ipr = ipr + 1
            r(ipr) = z(ier) +  dr(ier)*(xgl(i) +1)/2.0
            intmar(i,ier) = ipr
        end do
    end do
  
end subroutine create1d_grid

!----------------------------------------------------------------------!
!>@brief This subroutine constructs the HEXAHEDRAL GRID
!>@author  Francis X. Giraldo on 11/2014
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943
!----------------------------------------------------------------------!
subroutine hex_quad(coord_cg,intma,nel,ngl,nface,nelem,npoin_cg,nelem0,npoin0,nx,ny)
  
    use mod_input, only: geometry_type

    implicit none
  
    !global arrays
    real, intent(out) :: coord_cg(3,npoin_cg)
    integer, intent(out) :: intma(ngl,ngl,nelem)
    integer, intent(in) :: nel, ngl, nface, nelem, npoin_cg
    integer, intent(in) :: nelem0, npoin0, nx, ny

    !local arrays
    real, allocatable :: coors0(:,:,:), coord0(:,:,:), coorsf(:,:,:,:)
    real, allocatable :: coordf(:,:,:,:), coors(:,:)
    integer, allocatable :: intma0(:,:,:), nodef(:,:,:)

    !Face Information
    real lon(6), lat(6)

    !Local Variables
    real :: pi, twopi, pio2
    real :: clon, clat, alon, alat, rlon, rlat, r
    real :: x, y, z, b_gnom1, b_gnom2, b_rot1, b_rot2, s
    integer :: i, j, k, l, m, ip, np, ie, ii, jj
    integer :: AllocateStatus

    data lon / 0, 90, 180, 270,  0,   0 /
    data lat / 0,  0,   0,   0, 90, -90 /

    !local arrays
    allocate( coors0(2,nx,ny), coord0(2,nx,ny), coorsf(2,nx,ny,nface), &
        coordf(3,nx,ny,nface), coors(2,npoin_cg), &
        intma0(ngl,ngl,nelem0), nodef(nx,ny,nface), &
        stat=AllocateStatus )
    if (AllocateStatus /= 0) stop "** Not Enough Memory - CREATE_GRID_SPHERE_HEX/HEX_QUAD **"

    !Constants
    pi=4.0*atan(1.0)
    twopi=2*pi
    pio2=0.5*pi

    do i=1,6
        lon(i)=lon(i)*pi/180
        lat(i)=lat(i)*pi/180
    end do

    !Define domain size based on Projection type
    if (geometry_type(1:13) == 'sphere_hex_es') then
        s=1.0 !Equi-spaced
    else if (geometry_type(1:13) == 'sphere_hex_ea') then
        s=pi/4.0 !Equi-Angular
    end if

    !Construct the 2D Planar Grid
    call init_hex(coord0,intma0,nelem0,nel,nx,ny,ngl,s)

    !Get Rotated Spherical Coordinates
    if (geometry_type(1:13) == 'sphere_hex_es') then
        do j=1,ny
            do i=1,nx
                x=coord0(1,i,j)
                y=coord0(2,i,j)

                !Get Backward Gnomonic Projection - Equi-Spaced
                rlon=b_gnom1(x); rlat=b_gnom2(y,rlon)

                !Store Spherical coordinates
                coors0(1,i,j)=rlon
                coors0(2,i,j)=rlat

            end do !i
        end do !j
    else if (geometry_type(1:13) == 'sphere_hex_ea') then
        do j=1,ny
            do i=1,nx
                x=coord0(1,i,j)
                y=coord0(2,i,j)

                !Get Backward Gnomonic Projection - Equi-Angular
                call gnomonic(rlon,rlat,x,y)

                !Store Spherical coordinates
                coors0(1,i,j)=rlon
                coors0(2,i,j)=rlat

            end do !i
        end do !j
    end if

    !Loop through the 6 Faces
    ip=0
    do l=1,nface
        clon=lon(l)
        clat=lat(l)

        !loop through the Planar Grid and get the Face Coordinates
        do j=1,ny
            do i=1,nx

                !Get Rotated Spherical Coordinates
                rlon=coors0(1,i,j)
                rlat=coors0(2,i,j)

                !Get Backward Rotation Transformation
                alon=b_rot1(rlon,rlat,clon,clat)
                alat=b_rot2(rlon,rlat,clon,clat)

                if (alon < 0 ) alon=alon + twopi
                if (alon > twopi ) alon=alon - twopi

                !Store Spherical coordinates
                coorsf(1,i,j,l)=alon
                coorsf(2,i,j,l)=alat

            end do !i
        end do !j
    end do !l

    !!--------Unite Faces into Global Grid---------!!
    np=0

    !FACE 1
    do j=1,ny
        do i=1,nx
            np=np + 1
            nodef(i,j,1)=np
            do k=1,2
                coors(k,np)=coorsf(k,i,j,1)
            end do !k
        end do !i
    end do !j

    !FACE 2
    do j=1,ny
        nodef(1,j,2)=nodef(nx,j,1)
    end do !j
    do j=1,ny
        do i=2,nx
            np=np + 1
            nodef(i,j,2)=np
            do k=1,2
                coors(k,np)=coorsf(k,i,j,2)
            end do !k
        end do !i
    end do !j

    !FACE 3
    do j=1,ny
        nodef(1,j,3)=nodef(nx,j,2)
    end do !j
    do j=1,ny
        do i=2,nx
            np=np + 1
            nodef(i,j,3)=np
            do k=1,2
                coors(k,np)=coorsf(k,i,j,3)
            end do !k
        end do !i
    end do !j

    !FACE 4
    do j=1,ny
        nodef(1,j,4)=nodef(nx,j,3)
    end do !j
    do j=1,ny
        nodef(nx,j,4)=nodef(1,j,1)
    end do !j
    do j=1,ny
        do i=2,nx-1
            np=np + 1
            nodef(i,j,4)=np
            do k=1,2
                coors(k,np)=coorsf(k,i,j,4)
            end do !k
        end do !i
    end do !j

    !FACE 5
    do i=1,nx
        nodef(i,1,5)=nodef(i,ny,1)
        nodef(i,ny,5)=nodef(nx+1-i,ny,3)
    end do !j
    do j=1,ny
        nodef(1,j,5)=nodef(nx+1-j,ny,4)
        nodef(nx,j,5)=nodef(j,ny,2)
    end do !j
    do j=2,ny-1
        do i=2,nx-1
            np=np + 1
            nodef(i,j,5)=np
            do k=1,2
                coors(k,np)=coorsf(k,i,j,5)
            end do !k
        end do !i
    end do !j

    !FACE 6
    do i=1,nx
        nodef(i,1,6)=nodef(nx+1-i,1,3)
        nodef(i,ny,6)=nodef(i,1,1)
    end do !j
    do j=1,ny
        nodef(1,j,6)=nodef(j,1,4)
        nodef(nx,j,6)=nodef(nx+1-j,1,2)
    end do !j
    do j=2,ny-1
        do i=2,nx-1
            np=np + 1
            nodef(i,j,6)=np
            do k=1,2
                coors(k,np)=coorsf(k,i,j,6)
            end do !k
        end do !i
    end do !j

    if (np /= npoin_cg) then
        print*,' npoin_cg np = ',npoin_cg,np
        stop
    end if

    !GENERATE INTMA
    ie=0
    do m=1,nface
        do l=1,nel
            do k=1,nel
                ie=ie+1
                do j=1,ngl
                    jj=(ngl-1)*(l-1) + j
                    do i=1,ngl
                        ii=(ngl-1)*(k-1) + i
                        ip=nodef(ii,jj,m)
                        intma(i,j,ie)=ip
                    end do !i
                end do !j
            end do !k
        end do !l
    end do !m

    !Convert from Spherical to Cartesian
    do i=1,npoin_cg
        alon=coors(1,i)
        alat=coors(2,i)
        x=cos(alat)*cos(alon)
        y=cos(alat)*sin(alon)
        z=sin(alat)
        r=sqrt(x*x + y*y + z*z)
        coord_cg(1,i)=x/r
        coord_cg(2,i)=y/r
        coord_cg(3,i)=z/r

    end do !i

    !deallocate local arrays
    deallocate( coors0, coord0, coorsf, coordf, coors, &
        intma0, nodef)

end subroutine hex_quad

!----------------------------------------------------------------------!
!>@brief This subroutine constructs the INITIAL HEXAHEDRON
!>@author  Francis X. Giraldo on 11/2014
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943
!----------------------------------------------------------------------!
subroutine init_hex(coord,intma,nelem,nel,nx,ny,ngl,s)
  
    use mod_legendre, only : legendre_gauss_lobatto
  
    !global arrays
    integer nelem, nel,nx,ny,ngl
    real coord(2,nx,ny), s
    integer intma(ngl,ngl,nelem), node(nx,ny)
    character text*72, fname*150
  
    !local arrays
    real xgl(ngl), wgl(ngl)
  
    !set some constants
    pi=4.0*atan(1.0)
    xmin=-s
    xmax=+s
    ymin=-s
    ymax=+s
    dx=(xmax-xmin)/(nel)
    dy=(ymax-ymin)/(nel)
    xl=xmax-xmin
    yl=ymax-ymin
  
    !Generate Gauss-Lobatto Points for NPTS
    call legendre_gauss_lobatto(ngl,xgl,wgl)
  
    !GENERATE COORD
    ip=0
    jj=0
    do k=1,nel
        y0=ymin + real(k-1)*dy
     
        if (k == 1) then
            l1=1
        else
            l1=2
        endif
     
        do l=l1,ngl
            y=( xgl(l)+1 )*dy/2 + y0
        
            jj=jj+1
            ii=0
        
            do i=1,nel
                x0=xmin + real(i-1)*dx
           
                if (i == 1) then
                    j1=1
                else
                    j1=2
                endif
           
                do j=j1,ngl
                    ii=ii + 1
                    ip=ip + 1
                    x=( xgl(j)+1 )*dx/2 + x0
                    coord(1,ii,jj)=x
                    coord(2,ii,jj)=y
                    node(ii,jj)=ip
                end do
            end do
        end do
    end do
  
    !GENERATE INTMA
    ie=0
    do k=1,nel
        do i=1,nel
            ie=ie+1
            do l=1,ngl
                jj=(ngl-1)*(k-1) + l
                do j=1,ngl
                    ii=(ngl-1)*(i-1) + j
                    ip=node(ii,jj)
                    intma(j,l,ie)=ip
                end do
            end do
        end do
    end do
  
end subroutine init_hex

!------------------------------------------------------------------------!
!>@breif This is the Longitude Rotation Transformation - Forward
!------------------------------------------------------------------------!
function f_rot1(alon,alat,alono,alato)
  
    real f_rot1
    real alon, alat, alono, alato
  
    anum=cos(alat)*sin(alon-alono)
    aden=cos(alat)*cos(alon-alono)*cos(alato) + sin(alat)*sin(alato)
    f_rot1=atan2(anum,aden)
  
end function f_rot1
!------------------------------------------------------------------------!
!>@brief This is the Latitude Rotation Transformation - Forward
!------------------------------------------------------------------------!
function f_rot2(alon,alat,alono,alato)
  
    real f_rot2
    real alon, alat, alono, alato
  
    anum=sin(alat)*cos(alato) - cos(alat)*cos(alon-alono)*sin(alato)
    f_rot2=asin(anum)
  
end function f_rot2
!------------------------------------------------------------------------!
!>@brief This is the Longitude Rotation Transformation - Backward
!------------------------------------------------------------------------!
function b_rot1(alon_p,alat_p,alono,alato)
  
    real b_rot1
    real alon_p, alat_p, alono, alato
  
    anum=cos(alat_p)*sin(alon_p)
    aden=cos(alat_p)*cos(alon_p)*cos(alato) - sin(alat_p)*sin(alato)
    b_rot1=alono + atan2(anum,aden)
  
end function b_rot1
!------------------------------------------------------------------------!
!>@brief This is the Latitude Rotation Transformation - Backward
!------------------------------------------------------------------------!
function b_rot2(alon_p,alat_p,alono,alato)
  
    real b_rot2
    real alon_p, alat_p, alono, alato
  
    anum=sin(alat_p)*cos(alato) + cos(alat_p)*cos(alon_p)*sin(alato)
    b_rot2=asin(anum)
  
end function b_rot2
!------------------------------------------------------------------------!
!>@brief This is the Rotation Jacobian
!------------------------------------------------------------------------!
subroutine jac_rotation(jr,alat_p,alon,alat,alono,alato)
    real jr(2,2)
    real alat_p, alon, alat, alono, alato
  
    c=1.0/cos(alat_p)
  
    jr(1,1)=c*( cos(alat)*cos(alato) + sin(alat)*cos(alon-alono)*sin(alato) )
    jr(1,2)=c*(-sin(alon-alono)*sin(alato))
    jr(2,1)=c*(+sin(alon-alono)*sin(alato))
    jr(2,2)=c*( cos(alat)*cos(alato) + sin(alat)*cos(alon-alono)*sin(alato) )
  
end subroutine jac_rotation
!------------------------------------------------------------------------!
!>@brief This is the X Gnomonic Projection - Forward
!------------------------------------------------------------------------!
function f_gnom1(alon_p,alat_p)
  
    real f_gnom1
    real alon_p, alat_p
  
    f_gnom1=tan(alon_p)
  
end function f_gnom1
!------------------------------------------------------------------------!
!>@brief This is the Y Gnomonic Projection - Forward
!------------------------------------------------------------------------!
function f_gnom2(alon_p,alat_p)
  
    real f_gnom2
    real alon_p, alat_p
  
    f_gnom2=tan(alat_p)/cos(alon_p)
  
end function f_gnom2
!------------------------------------------------------------------------!
!>@brief This is the X Gnomonic Projection - Backward
!------------------------------------------------------------------------!
function b_gnom1(x)
  
    real b_gnom1
    real x
  
    b_gnom1=atan(x)
  
end function b_gnom1
!------------------------------------------------------------------------!
!>@brief This is the Y Gnomonic Projection - Backward
!------------------------------------------------------------------------!
function b_gnom2(y,alon_p)
  
    real b_gnom2
    real y, alon_p
  
    b_gnom2=atan(y*cos(alon_p))
  
end function b_gnom2

!------------------------------------------------------------------------!
!>@brief Cartesian to Spherical coordinate transformation
!------------------------------------------------------------------------!
subroutine cart2sph(x,y,z,r,phi,lambda)
  
    use mod_constants, only: tol
    implicit none
    real x,y,z,r,phi,lambda
    real r2
  
    ! Convert to Spherical Coords
    r2 = x*x + y*y + z*z
    r = sqrt(r2)
    phi = asin(z/r)
    lambda = atan2(y,(x + tol))
  
end subroutine cart2sph

!------------------------------------------------------------------------!
!>@brief Spherical to Cartesian coordinate transformation
!------------------------------------------------------------------------!
subroutine sph2cart(x,y,z,r,phi,lambda)
  
    use mod_constants, only: tol
    implicit none
    real x,y,z,r,phi,lambda
    real r2
  
    x = r*cos(phi)*cos(lambda)
    y = r*cos(phi)*sin(lambda)
    z = r*sin(phi)
  
end subroutine sph2cart
!------------------------------------------------------------------------!
!> @brief Maps points on inscribed cubes to the cubed-sphere spherical grid.
!> Given (x,y,z) in a cube we get (x,y,z) in a cubed-sphere.
!> MAK: r_correction introduced to account for a bug in p6est
!>
!> @author F.X. Giraldo
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943
!>
!> @date November 20, 2014
!------------------------------------------------------------------------!
subroutine cube2sphere(x,y,z,r_correction,face)
  
    use mod_constants, only: earth_radius

    use mod_input, only: ztop, geometry_type

    implicit none

    !global variables
    real, intent(inout) :: x, y, z
    real, intent(in) :: r_correction
    integer, intent(out) :: face

    !local variables
    real :: lon_f, lat_f !face space
    real :: x_g, y_g, lon_g, lat_g !gnomonic space
    real :: lon_r, lat_r !rotated gnomonic space
    real :: r, rad, a, b, s, pi, factor
    real :: b_gnom1, b_gnom2, b_rot1, b_rot2
  
    !Constants
    pi=4.0*atan(1.0) !FXG_P4est: move PI to MOD_CONSTANTS once code is verified
    !a=1.0; b=1.0*r_correction !Debugging
    a=earth_radius; b=ztop*r_correction !Real Earth

    !Determinne which mapping to use
    if (geometry_type(1:13) == 'sphere_hex_es') then
        s=1.0 !Equi-spaced
    else if (geometry_type(1:13) == 'sphere_hex_ea') then
        s=pi/4.0 !Equi-spaced
    end if

    !determine which face of cubed-sphere (x,y,z) is on
    call find_panel_face(face,lon_f,lat_f,x,y,z)

    !Map (x,y,z) to (X_G,Y_G) - the gnomonic coordinates on the root cubed-sphere panel in space [-pi/4,+pi/4]
    select case (face)
        case (1) !(1,0,0)
            rad=abs(x)
            factor=s/rad
            x_g=+y*factor; y_g=+z*factor
        case (2) !(0,1,0)
            rad=abs(y)
            factor=s/rad
            x_g=-x*factor; y_g=+z*factor
        case (3) !(-1,0,0)
            rad=abs(x)
            factor=s/rad
            x_g=-y*factor; y_g=+z*factor
        case (4) !(0,-1,0)
            rad=abs(y)
            factor=s/rad
            x_g=+x*factor; y_g=+z*factor
        case (5) !(0,0,1)
            rad=abs(z)
            factor=s/rad
            x_g=+y*factor; y_g=-x*factor
        case (6) !(0,0,-1)
            rad=abs(z)
            factor=s/rad
            x_g=+y*factor; y_g=+x*factor
    end select
    r=a + b*(rad-1.0);

    !Map (X_G,Y_G) to (Lon_G,Lat_G)
    if (geometry_type(1:13) == 'sphere_hex_es') then
        lon_g=b_gnom1(x_g); lat_g=b_gnom2(y_g,lon_g) !Equi-spaced
    else if (geometry_type(1:13) == 'sphere_hex_ea') then
        call gnomonic(lon_g,lat_g,x_g,y_g) !Equi-angular
    end if

    !Rotate point based on which face of Cube it's on
    lon_r=b_rot1(lon_g,lat_g,lon_f,lat_f)
    lat_r=b_rot2(lon_g,lat_g,lon_f,lat_f)

    !using R, expand to the proper vertical height
    call sph2cart(x,y,z,r,lat_r,lon_r)

end subroutine cube2sphere

!------------------------------------------------------------------------!
!> @brief Finds which panel of the cubed sphere a point is on.
!> Given (x,y,z) we get the FACE of the cubed-sphere and the LON and LAT at the 
!> center of the cubed-sphere.
!>
!> @author F.X. Giraldo
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943
!>
!> @date November 20, 2014
!------------------------------------------------------------------------!
subroutine find_panel_face(face,lon_f,lat_f,x,y,z)
  
    use mod_constants, only: pi

    implicit none

    !global variables
    integer, intent(out) :: face
    real,    intent(out) :: lon_f, lat_f
    real,    intent(in)  :: x, y, z

    !local variables
    real, dimension(6) :: xf, yf, zf, lon, lat
    real :: dist, dist_min, dx, dy, dz
    integer :: i
    !Store Panel Face Centers
    data xf / 1, 0, -1, 0, 0, 0 /
    data yf / 0, 1,  0,-1, 0, 0 /
    data zf / 0, 0,  0, 0, 1, -1 /
    data lon / 0, 90, 180, 270,  0,   0 /
    data lat / 0,  0,   0,   0, 90, -90 /
  
    !determine which face (x,y,z) is on
    dist_min=+1e10
    face=0
    do i=1,6
        dx=x-xf(i); dy=y-yf(i); dz=z-zf(i)
        dist=dx*dx + dy*dy + dz*dz
        if (dist < dist_min) then
            dist_min=dist
            face=i
        end if
    end do

    lon_f=lon(face)*pi/180
    lat_f=lat(face)*pi/180

end subroutine find_panel_face

!------------------------------------------------------------------------!
!This is the Equi-Angular Gnomonic Projection
! See FXG Notebook 16 p. 11.25.02.1
!------------------------------------------------------------------------!
subroutine gnomonic(olon,olat,xg,yg)

    implicit none

    real :: olon, olat
    real :: xg, yg, xx, yy, zz

    xx=tan(xg)
    yy=tan(yg)
    zz=1.0/sqrt(1.0+xx*xx+yy*yy)

    olon=xg
    olat=asin( zz*tan(yg) )

end subroutine gnomonic

!------------------------------------------------------------------------!
!This is the Equi-Angular Gnomonic Projection
!------------------------------------------------------------------------!
subroutine backward_gnomonic(olon,olat,xg,yg)

    implicit none

    real :: olon, olat
    real :: x, y, z
    real :: xg, yg

    !Cartesian coords
    x=cos(olat)*cos(olon)
    y=cos(olat)*sin(olon)
    z=sin(olat)

    xg=atan2(y,x)
    yg=atan2(z,x)

end subroutine backward_gnomonic

