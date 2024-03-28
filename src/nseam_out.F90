!---------------------------------------------------------------------!
!>@brief nseam_out
!---------------------------------------------------------------------!
subroutine nseam_out(q,coord_s,npoin_s)
  
    use mod_constants, only: cv, cp, rgas, p00, earth_radius
  
    implicit none
  
    real eps
    integer npoin_s, i
    real q(5,npoin_s)
    real coord_s(3,npoin_s)
    real x, y, z, phi, lambda
    real sinphi, cosphi, sinlambda, coslambda
    real r
    real rho, u, v, w, theta, temp, press
    real us, vs, ws, pi
  
    pi = 4.0*atan(1.0)
    eps = 1e-6
  
    do i = 1,npoin_s
        x=coord_s(1,i)
        y=coord_s(2,i)
        z=coord_s(3,i)
        r = sqrt(x*x + y*y + z*z)
        lambda=atan2(y,x+eps)
        phi =asin(z/r)
        if ( (x*x + y*y)/r < eps ) lambda =pi/2
     
        sinphi = sin(phi)
        cosphi = cos(phi)
        sinlambda = sin(lambda)
        coslambda = cos(lambda)
     
        rho = q(1,i)
        u=q(2,i)
        v=q(3,i)
        w=q(4,i)
        theta = q(5,i)
        temp = theta*(rho*rgas*theta/p00)**(rgas/cv)
        press = p00 * (rho*rgas*theta/p00)**(cp/cv)
     
        vs = -sinphi*coslambda*u - sinphi*sinlambda*v + cosphi*w
        us = -sinlambda*u + coslambda*v
        ws = cosphi*coslambda*u + cosphi*sinlambda*v + sinphi*w
     
        !vs=w/cos(rlat+eps)
        !if (abs(sin(rlon)) > 0.001 ) then
        !   us=-( u + vs*sin(rlat)*cos(rlon) )/sin(rlon)
        !else
        !   us=+( v + vs*sin(rlat)*sin(rlon) )/cos(rlon)
        !end if
     
        !Construct pressure, zonal vel., meridional vel., and temperature
        q(1,i) = press
        q(2,i) = us
        q(3,i) = vs
        q(4,i) = ws
        q(5,i) = temp
     
    end do
  
end subroutine nseam_out
