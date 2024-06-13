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
