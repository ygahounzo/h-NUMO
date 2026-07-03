!-----------------------------------------------------------------!
!>@brief Initial conditions for sphere test cases (MLSWE solver).
!>@author Yao Gahounzo
!-----------------------------------------------------------------!

subroutine initial_conditions_sphere(q_df, pbprime_df, qb_df, alpha_mlswe, &
                                    pbprime_df_face, zbot_df, kvector, nlayers, &
                                    b, G, inp, mf)

   use mod_grid,          only: grid
   use mod_basis,         only: basis
   use mod_input,         only: input
   use mod_face,          only: face_CS
   use mod_constants,     only: gravity, pi, omega, earth_radius
   use mod_initial_mlswe, only: interpolate_pbprime_init

   implicit none

   type(basis),   intent(in) :: b
   type(grid),    intent(in) :: G
   type(input),   intent(in) :: inp
   type(face_CS), intent(in) :: mf

   integer, intent(in) :: nlayers
   real, dimension(inp%nvar_bcl, G%npoin, nlayers), intent(out) :: q_df
   real, dimension(G%npoin),                        intent(out) :: pbprime_df
   real, dimension(inp%nvar_btp, G%npoin),          intent(out) :: qb_df
   real, dimension(nlayers),                        intent(out) :: alpha_mlswe
   real, dimension(2, b%ngl, G%nface),              intent(out) :: pbprime_df_face
   real, dimension(G%npoin),                        intent(out) :: zbot_df
   real, dimension(3, G%npoin),                     intent(in)  :: kvector

   integer :: npoin
   real, dimension(G%npoin, nlayers)   :: u_df, v_df, w_df
   real, dimension(G%npoin, nlayers+1) :: z_interface
   integer :: I1, k
   real :: x, y, z, olon, olat, r, u, v
   real :: twopi, pio2, hc, h0, p, w_rot, rc, tt
   real :: oloni, olati, olonc, olatc, rlon, rlat, beta
   real :: tolon, tolat, talon, talat, alon, alat, alpha
   real :: hs0, rs, u0, atol

   ! External rotation functions (defined at bottom of this file).
   real :: f_rot1, f_rot2, b_rot1, b_rot2

   npoin = G%npoin
   atol  = 1.0e-6   ! guard for atan2 at x=0

   z_interface     = 0.0
   pbprime_df      = 0.0
   q_df            = 0.0
   u_df            = 0.0
   v_df            = 0.0
   w_df            = 0.0
   pbprime_df_face = 0.0

   ! Rotation angle of the initial flow (radians); 0 = no rotation.
   alpha = 0.0

   select case (trim(inp%test_case))

   case ('solid_body') ! Solid body rotation

      beta  = 0.0
      twopi = 2.0*pi
      pio2  = pi/2.0
      hc    = 1000.0
      h0    = 1000.0
      p     = 24.0*3600.0
      w_rot = twopi*earth_radius/(12.0*p)
      oloni = 0.0
      olati = 0.0
      olonc = 3.0*pi/2.0
      olatc = 0.0
      rc    = earth_radius/3.0
      tt    = w_rot*0.0/earth_radius

      rlon = 0.0
      rlat = alpha

      zbot_df(:) = -h0
      do k = 1, nlayers+1
         z_interface(:,k) = -(k-1)*h0/real(nlayers)
      end do
      z_interface(:,nlayers+1) = zbot_df(:)
      alpha_mlswe(1) = 9.7370e-04
      if (nlayers >= 2) alpha_mlswe(2) = 9.7350e-04

      do I1 = 1, npoin
         x = G%coord(1,I1)
         y = G%coord(2,I1)
         z = G%coord(3,I1)

         olon = atan2(y, x+atol)
         olat = asin(z/earth_radius)

         tolon = f_rot1(olon, olat, rlon, rlat)
         tolat = f_rot2(olon, olat, rlon, rlat)

         talon = tolon - 0.0*(tt + beta*sin(tt))
         talat = tolat

         alon = b_rot1(talon, talat, rlon, rlat)
         alat = b_rot2(talon, talat, rlon, rlat)

         if (alon < 0.0)    alon = twopi - (0.0 - alon)
         if (alon > twopi)  alon = 0.0   + (alon - twopi)

         if (alat > pio2) then
            alat = pio2 - (alat - pio2)
            if (alon <= pi) then; alon = alon + pi; else; alon = alon - pi; end if
         end if
         if (alat < -pio2) then
            alat = -pio2 - (alat + pio2)
            if (alon <= pi) then; alon = alon + pi; else; alon = alon - pi; end if
         end if

         if ((alon < 0.0 .or. alon > twopi) .or. (alat < -pio2 .or. alat > pio2)) then
            print *, ' Error in DEPART'
            print *, 'ALON out of Range = ', I1, alon, alat
            stop
         end if

         r = earth_radius*acos( sin(olati)*sin(alat) + &
             cos(olati)*cos(alat)*cos(alon-oloni) )
         if (r < rc) then
            z_interface(I1,1) = z_interface(I1,1) + 0.5*h0*(1.0 + cos(pi*r/rc))
         end if

         u = +w_rot*(cos(olat)*cos(alpha) + sin(olat)*cos(olon+olonc)*sin(alpha))
         v = -w_rot*sin(olon+olonc)*sin(alpha)

         u_df(I1,:) = (-u*sin(olon) - v*sin(olat)*cos(olon))
         v_df(I1,:) = (+u*cos(olon) - v*sin(olat)*sin(olon))
         w_df(I1,:) = (+v*cos(olat))
      end do

   case ('geoFlow') ! Zonal geostrophic flow — solid body rotation

      twopi = 2.0*pi
      pio2  = pi/2.0
      h0    = 2.9981e3
      p     = 24.0*3600.0
      w_rot = twopi*earth_radius/(12.0*p)
      oloni = 0.0
      olati = 0.0
      rc    = earth_radius/3.0
      rlon  = 0.0
      rlat  = alpha

      zbot_df(:) = -h0

      z_interface(:,2) = -h0 + 800.0
      if (nlayers >= 3) z_interface(:,3) = -h0 + 400.0

      alpha_mlswe(1) = 1.0/1000.0
      if (nlayers >= 2) alpha_mlswe(2) = 1.0/1005.0
      if (nlayers >= 3) alpha_mlswe(3) = 1.0/1010.0

      do I1 = 1, npoin
         x = G%coord(1,I1)
         y = G%coord(2,I1)
         z = G%coord(3,I1)
         olon = atan2(y, x+atol)
         olat = asin(z/earth_radius)

         if (olon < 0.0)    olon = olon + twopi
         if (olon > twopi)  olon = olon - twopi
         if (olat < -pi)    olat = olat + pi
         if (olat >  pi)    olat = olat - pi

         z_interface(I1,1) = -((earth_radius*omega*w_rot + 0.5*w_rot**2) * &
              (-cos(olon)*cos(olat)*sin(alpha) + sin(olat)*cos(alpha))**2) / gravity

         u = +w_rot*(cos(olat)*cos(alpha) + sin(olat)*cos(olon)*sin(alpha))
         v = -w_rot*sin(olon)*sin(alpha)

         u_df(I1,:) = -u*sin(olon) - v*sin(olat)*cos(olon)
         v_df(I1,:) = +u*cos(olon) - v*sin(olat)*sin(olon)
         w_df(I1,:) = +v*cos(olat)
      end do

      z_interface(:,nlayers+1) = zbot_df(:)

   case ('sphere_mountain') ! Flow over mountain

      twopi = 2.0*pi
      pio2  = pi/2.0
      h0    = 5960.0
      p     = 24.0*3600.0
      u0    = 20.0
      oloni = 0.0
      olati = 0.0
      rc    = earth_radius
      hs0   = 2000.0
      rs    = pi/9.0
      olonc = 3.0*pi/2.0
      olatc = pi/6.0
      rlon  = 0.0
      rlat  = alpha

      z_interface(:,2) = -3500.0
      alpha_mlswe(:)   = 1.0

      do I1 = 1, npoin
         x = G%coord(1,I1)
         y = G%coord(2,I1)
         z = G%coord(3,I1)
         olon = atan2(y, x+atol)
         olat = asin(z/earth_radius)

         if (olon < 0.0)    olon = olon + twopi
         if (olon > twopi)  olon = olon - twopi
         if (olat < -pi)    olat = olat + pi
         if (olat >  pi)    olat = olat - pi

         r = sqrt(min(rs*rs, (olon-olonc)**2 + (olat-olatc)**2))
         zbot_df(I1) = -h0 + hs0*(1.0 - r/rs)

         z_interface(I1,1) = -((earth_radius*omega*u0 + 0.5*u0**2) * &
              (-cos(olon)*cos(olat)*sin(alpha) + sin(olat)*cos(alpha))**2) / gravity

         u = +u0*(cos(olat)*cos(alpha) + sin(olat)*cos(olon)*sin(alpha))
         v = -u0*sin(olon)*sin(alpha)

         u_df(I1,:) = -u*sin(olon) - v*sin(olat)*cos(olon)
         v_df(I1,:) = +u*cos(olon) - v*sin(olat)*sin(olon)
         w_df(I1,:) = +v*cos(olat)
      end do

      z_interface(:,nlayers+1) = zbot_df(:)

   case default
      print *, "Unknown sphere test case: ", trim(inp%test_case)
      stop

   end select

   ! Clamp interfaces to bathymetry; then enforce monotonicity (no negative layer thickness).
   do I1 = 1, npoin
      do k = 1, nlayers+1
         z_interface(I1,k) = max(zbot_df(I1), z_interface(I1,k))
      end do
      do k = 2, nlayers
         z_interface(I1,k) = min(z_interface(I1,k), z_interface(I1,k-1))
      end do
   end do

   ! Reference pressure pbprime_df = vertical sum of g/alpha * dz.
   do k = 1, nlayers
      if (trim(inp%test_case) == 'solid_body') then
         pbprime_df(:) = pbprime_df(:) + (gravity/alpha_mlswe(1))*(z_interface(:,k) - z_interface(:,k+1))
      else
         pbprime_df(:) = pbprime_df(:) + (gravity/alpha_mlswe(k))*(z_interface(:,k) - z_interface(:,k+1))
      end if
   end do

   call interpolate_pbprime_init(pbprime_df_face, pbprime_df, b, G, inp, mf)

   ! Layer pressure thicknesses dp = g/alpha * dz.
   do k = 1, nlayers
      if (trim(inp%test_case) == 'solid_body') then
         q_df(1,:,k) = (gravity/alpha_mlswe(1))*(z_interface(:,k) - z_interface(:,k+1))
      else
         q_df(1,:,k) = (gravity/alpha_mlswe(k))*(z_interface(:,k) - z_interface(:,k+1))
      end if
   end do

   ! Momentum components: u*dp, v*dp, w*dp.
   do k = 1, nlayers
      q_df(2,:,k) = u_df(:,k)*q_df(1,:,k)
      q_df(3,:,k) = v_df(:,k)*q_df(1,:,k)
      q_df(4,:,k) = w_df(:,k)*q_df(1,:,k)
   end do

   ! Barotropic state: vertical sums.
   qb_df = 0.0
   do k = 1, nlayers
      qb_df(1,:) = qb_df(1,:) + q_df(1,:,k)
      qb_df(3,:) = qb_df(3,:) + q_df(2,:,k)
      qb_df(4,:) = qb_df(4,:) + q_df(3,:,k)
      qb_df(5,:) = qb_df(5,:) + q_df(4,:,k)
   end do
   qb_df(2,:) = qb_df(1,:) - pbprime_df(:)

end subroutine initial_conditions_sphere

!-----------------------------------------------------------------!
! Sphere rotation helpers                                          !
! Forward rotation: (lon,lat) in original frame → (lon_p,lat_p)  !
! in frame rotated about pole (alono,alato).                      !
!-----------------------------------------------------------------!

real function f_rot1(alon, alat, alono, alato)
   implicit none
   real, intent(in) :: alon, alat, alono, alato
   real :: anum, aden
   anum   = cos(alat)*sin(alon - alono)
   aden   = cos(alat)*cos(alon-alono)*cos(alato) + sin(alat)*sin(alato)
   f_rot1 = atan2(anum, aden)
end function f_rot1

real function f_rot2(alon, alat, alono, alato)
   implicit none
   real, intent(in) :: alon, alat, alono, alato
   f_rot2 = asin(sin(alat)*cos(alato) - cos(alat)*cos(alon-alono)*sin(alato))
end function f_rot2

! Backward rotation: (lon_p,lat_p) in rotated frame → (lon,lat). !

real function b_rot1(alon_p, alat_p, alono, alato)
   implicit none
   real, intent(in) :: alon_p, alat_p, alono, alato
   real :: anum, aden
   anum   = cos(alat_p)*sin(alon_p)
   aden   = cos(alat_p)*cos(alon_p)*cos(alato) - sin(alat_p)*sin(alato)
   b_rot1 = alono + atan2(anum, aden)
end function b_rot1

real function b_rot2(alon_p, alat_p, alono, alato)
   implicit none
   real, intent(in) :: alon_p, alat_p, alono, alato
   b_rot2 = asin(sin(alat_p)*cos(alato) + cos(alat_p)*cos(alon_p)*sin(alato))
end function b_rot2
