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

   case ('solid_body') ! Williamson et al. 1992, Test 1: advection of a cosine
                        ! bell by solid-body rotation. The momentum restore in
                        ! layer_mom_boundary_df/btp_mom_boundary_df plays the
                        ! role of "overwrite the predicted wind field every
                        ! step with the analytic advecting wind" (per the
                        ! paper's own recipe for codes solving the full
                        ! shallow water equations), so the bell is advected
                        ! passively rather than dynamically generating its
                        ! own geostrophic adjustment.

      beta  = 0.0
      twopi = 2.0*pi
      pio2  = pi/2.0
      hc    = 1000.0
      h0    = 1000.0
      p     = 24.0*3600.0
      w_rot = twopi*earth_radius/p   ! flow speed for a full circumnavigation in 12 h
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

   case ('geoFlow') ! Williamson et al. 1992, Test 2: steady zonal geostrophic
                    ! flow. Unlike 'solid_body', the elevation field below
                    ! includes the geostrophic-balance term that exactly
                    ! cancels the Coriolis/centrifugal force of the imposed
                    ! wind, so this is a genuine steady state of the full
                    ! shallow water equations (no momentum restore needed).

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

   case ('gravity_wave') ! Chen (2025, QJRMS, 10.1002/qj.4994) Sec 5.2: internal
                         ! gravity-wave propagation test. Two layers, densities
                         ! 1000/1020 kg/m^3 (paper's first configuration; 1040
                         ! is its second), flow initially AT REST with a
                         ! cosine-bell hump on the layer interface, non-rotating
                         ! (f=0, set in wind_stress_coriolis). Excites a fast
                         ! barotropic and a slower baroclinic gravity wave;
                         ! validated by comparing their propagation speeds
                         ! against linear theory, not by a steady-state check.
                         ! Requires nlayers=2.

      hc    = 100.0                 ! hump amplitude (m)
      rc    = earth_radius*pi/18.0  ! hump radius (10 deg -> 20 deg wide)
      oloni = 0.0                   ! hump center: equator, lon=0
      olati = 0.0

      zbot_df(:)     = -2000.0      ! total resting depth: 1000 + 1000 m
      z_interface(:,1) = 0.0        ! flat top surface (at rest)

      alpha_mlswe(1) = 1.0/1000.0
      alpha_mlswe(2) = 1.0/1020.0   ! paper's other tested value: 1.0/1040.0

      do I1 = 1, npoin
         x = G%coord(1,I1)
         y = G%coord(2,I1)
         z = G%coord(3,I1)
         olon = atan2(y, x+atol)
         olat = asin(z/earth_radius)

         r = earth_radius*acos( sin(olati)*sin(olat) + &
             cos(olati)*cos(olat)*cos(olon-oloni) )

         ! phi2 = bottom-layer thickness: 1000 m background + cosine-bell hump.
         if (r < rc) then
            z_interface(I1,2) = 1000.0 + 0.5*hc*(1.0 + cos(pi*r/rc)) - 2000.0
         else
            z_interface(I1,2) = 1000.0 - 2000.0
         end if
      end do
      ! u_df, v_df, w_df stay 0 (flow at rest); Coriolis is set to 0 for this
      ! test case in wind_stress_coriolis (mod_initial_mlswe.F90).

      z_interface(:,nlayers+1) = zbot_df(:)

   case ('sphere_mountain') ! Williamson et al. (1992), Test 5: single-layer zonal
                            ! flow over an isolated mountain. Standard, widely
                            ! validated shallow-water benchmark, used here as a
                            ! single-layer baseline for the two-layer
                            ! 'mountain_2layer' case below (same wind field,
                            ! height balance, and mountain). Requires nlayers=1.

      twopi = 2.0*pi
      pio2  = pi/2.0
      h0    = 5960.0
      u0    = 20.0
      hs0   = 2000.0
      rs    = pi/9.0
      olonc = 3.0*pi/2.0
      olatc = pi/6.0
      rlat  = alpha  ! alpha=0: standard (non-rotated) test 5 configuration

      alpha_mlswe(1) = 1.0/1000.0

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

   case ('mountain_2layer') ! Chen (2025, QJRMS, 10.1002/qj.4994) Sec 5.3: two-layer
                             ! zonal flow over mountain topography, adapting
                             ! Williamson et al. (1992) test #5 to two density
                             ! layers (1000/1010 kg/m^3). Same balanced zonal
                             ! wind and top-surface height field as the existing
                             ! single-layer 'sphere_mountain' case above; the
                             ! water column is split at a flat interface, 3500 m
                             ! above the flat sea floor, into a light top layer
                             ! and a denser bottom layer that submerges the
                             ! mountain. Chen's own results show this case
                             ! develops negative bottom-layer thickness
                             ! (interface outcropping) near the mountain within
                             ! ~15h without an added artificial-potential-energy
                             ! stabilization term, which h-NUMO does not
                             ! currently implement. Requires nlayers=2.

      twopi = 2.0*pi
      pio2  = pi/2.0
      h0    = 5960.0
      u0    = 20.0
      hs0   = 2000.0
      ! hs0   = 0.0
      rs    = pi/9.0
      olonc = 3.0*pi/2.0
      olatc = pi/6.0
      rlat  = alpha  ! alpha=0: standard (non-rotated) test 5 configuration

      alpha_mlswe(1) = 1.0/1000.0
      alpha_mlswe(2) = 1.0/1010.0

      ! Flat layer interface at 3500 m above the flat sea floor, expressed in
      ! h-NUMO's zbot_df-relative height convention (sea floor at -h0).
      z_interface(:,2) = 3500.0 - h0
      ! z_interface(:,2) = 4000.0 - h0

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
