!-----------------------------------------------------------------!
!>@brief This subroutine builds the Initial Conditions for the multilayer Shallow Water Equations
!>@author Written by Yao Gahounzo
!>@date March 27, 2023
!-----------------------------------------------------------------!

subroutine initial_conditions(q_df, pbprime_df, qb_df, alpha, pbprime_df_face, zbot_df, &
   tau_wind_df, z_interface)
    

   use mod_grid, only: nelem, nface, npoin_q, npoin, coord, intma_dg_quad, face, intma
   use mod_constants, only: gravity, pi, tol, omega, earth_radius
   use mod_initial_mlswe, only: interpolate_from_dof_to_quad_uv_init, interpolate_pbprime_init, poslimiter
   use mod_basis, only: ngl, nq
   use mod_initial, only: kvector
   use mod_input, only: gravity_in, &
                        nelx, nelz, eqn_set, &
                        xdims, ydims,nlayers, test_case, dry_cutoff
   use mod_types, only : r8
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION
   use mod_face, only: imapl_q
   use mod_variables, only: z_interface_initial, z_init_flag, z_init_flag_elem

   implicit none
   
   real, dimension(3,npoin, nlayers) :: q_df
   real, dimension(npoin) :: pbprime_df
   real, dimension(2,ngl,nface) :: pbprime_df_face
   real, dimension(npoin) :: zbot_df
   real, dimension(4,npoin) :: qb_df
   real, dimension(nlayers) :: alpha
   real, dimension(2,npoin) :: tau_wind_df
   real, dimension(npoin, nlayers) :: u_df, v_df(npoin, nlayers)
   real, dimension(npoin) :: one_plus_eta_temp
   real, dimension(npoin_q) :: pb_temp, one_plus_eta_temp1
   real, dimension(npoin, nlayers+1) :: z_interface
   real, dimension(npoin,nlayers) :: height_layer

   integer :: e, m, n, I1, iquad, Iq, I2, i, j
   real :: x, y, xmid, ymid, L, amp, r, rho_0, H_bot
   integer :: ierr, iface, ilr, k,ip, il, jl, kl, el
   real :: xmax, xmin, ymax, ymin, zmin, zmax, yl, xm
   real :: xmax_l, xmin_l, ymax_l, ymin_l, zmin_l, zmax_l, Ly
   real :: delta, dp_dry
   real :: c_jump, indep(nlayers+1)
   real, dimension(npoin, nlayers+1) :: z_init
   real :: z_interface_equil(nlayers+1), layer_dz_eq(nlayers)
   real :: xl, xr, ztmp, zl, zr, slope_edge_adjust, slope_init_l, slope_init_r
   real :: dLdx, dzbot_dx, dx, dzbot_dy, dy, yr, dLdy
   
   xmin_l=minval(coord(1,:)); xmax_l=maxval(coord(1,:))
   ymin_l=minval(coord(2,:)); ymax_l=maxval(coord(2,:))
   zmin_l=minval(coord(3,:)); zmax_l=maxval(coord(3,:))

   call mpi_allreduce(xmin_l,xmin,1,MPI_PRECISION,mpi_min,mpi_comm_world,ierr)
   call mpi_allreduce(xmax_l,xmax,1,MPI_PRECISION,mpi_max,mpi_comm_world,ierr)
   call mpi_allreduce(ymin_l,ymin,1,MPI_PRECISION,mpi_min,mpi_comm_world,ierr)
   call mpi_allreduce(ymax_l,ymax,1,MPI_PRECISION,mpi_max,mpi_comm_world,ierr)
   call mpi_allreduce(zmin_l,zmin,1,MPI_PRECISION,mpi_min,mpi_comm_world,ierr)
   call mpi_allreduce(zmax_l,zmax,1,MPI_PRECISION,mpi_max,mpi_comm_world,ierr)

   z_init = 0.0
   z_interface = 0.0
   pbprime_df = 0.0
   q_df = 0.0
   u_df = 0.0
   v_df = 0.0
   tau_wind_df = 0.0
   pbprime_df_face = 0.0
   z_interface_equil = 0.0

   kvector(1,:)=0.0
   kvector(2,:)=0.0
   kvector(3,:)=1.0

   Ly = ydims(2)-ydims(1)

   ! The following is for a two-layer configuration.
   ! The free surface is level, the interface is one cosine (MAK test)
   ! hump centered in the interval, and u = v = 0.
   ! The motion is mostly internal.

   select case (trim(test_case))
      
   case ("bump") ! bump test (wave propagation)

      gravity = 9.806
      H_bot = 40.0   ! total depth

      ! Bottom topography (flat)
      zbot_df(:) = - H_bot
      layer_dz_eq(:) = H_bot/real(nlayers) ! layer thicknesses (m) at global equilibrium

      ! Bump centered at (xm,yl) with radius L and amplitude amp at the second layer interface
      xm=0.5*(xmax+xmin)
      yl=0.5*(ymax+ymin)

      L = 250.0
      amp = 1.0

      do I1 = 1, npoin

         x = coord(1,I1)
         y = coord(2,I1)
         r = sqrt((x-xm)**2 + (y-yl)**2)

         if (r < L) then
            z_init(I1,2) =  0.5*amp*(1.0 + cos(pi*r/L))
         end if
      end do

      ! Layer densities reciprocal (1/rho)
      alpha(1) = 0.9737e-3
      alpha(2) = 0.9735e-3

   case("lakeAtrest") ! lake at rest (well-balanced) test 

      gravity = 9.806
      H_bot = 40.0   ! total depth

      ! Bottom topography (non-flat)
      zbot_df(:) = - H_bot

      xm=0.5*(xdims(1)+xdims(2))
      yl=0.5*(ydims(1)+ydims(2))
      L = 250.0

      do I1 = 1,npoin 

         x = coord(1,I1)
         y = coord(2,I1)
         r = sqrt((x-xm)**2 + (y-yl)**2)

         if (r < L) then
            zbot_df(I1) = zbot_df(I1) + 3.0*(1.0 + cos(pi*r/L))
         end if
      end do  

      ! layer thicknesses (m) at global equilibrium
      if(nlayers < 5) then
         layer_dz_eq(:) = H_bot/real(nlayers)
      else
         layer_dz_eq(1:nlayers-1) = 32.0/real(nlayers-1)
         layer_dz_eq(nlayers) = H_bot - 32.0
      endif

      ! Layer densities reciprocal (1/rho)
      rho_0 = 1027.01037
      alpha(1) = 1.0/rho_0

      do k = 2,nlayers
         alpha(k) = 1.0/(rho_0 + k*0.2110/real(nlayers))
      end do

   case ("double-gyre") ! double-gyre 

      gravity = 9.806
      H_bot = 9928.0  ! total depth

      ! Bottom topography (flat)
      zbot_df(:) = - H_bot

      ! layer thicknesses (m) at global equilibrium
      layer_dz_eq(1) = 1489.5
      layer_dz_eq(2) = H_bot - 1489.5

      ! Layer densities reciprocal (1/rho)
      alpha(1) = 9.7370e-04
      alpha(2) = 9.7350e-04

      ! wind stress
      do I1 = 1, npoin
         y = coord(2,I1)
         tau_wind_df(1,I1) = -0.1*cos(2.0*pi*y/Ly)
      end do

   case ("dam") ! dam-break test

      gravity = 9.806
      H_bot = 2000.0 ! total depth

      ! Bottom topography
      do I1 = 1, npoin
        
         x = coord(1,I1)/1.0e3
         y = coord(2,I1)/1.0e3

         ! if(y <= 600.0) then

         !    if(y <= 300.0) then
         !       zbot_df(I1) = H_bot
         !    else
         !       zbot_df(I1) = H_bot - 10.0*(y - 300.0)
         !    end if
         ! else
         !    if(400.0 <= x .and. x <= 500.0) then
         !       zbot_df(I1) = 600.0
         !    end if
         ! end if

         zbot_df(I1) = 500.0 + 0.5*(H_bot - 500.0)*(1.0 + tanh((x - 401.0)/60.0))

         zbot_df(I1) = - zbot_df(I1)
      end do
        
      ! layer thicknesses (m) at global equilibrium
      do k = 1, nlayers
        layer_dz_eq(k) = H_bot*(real(k)-0.5)/real(nlayers-1)
      enddo

      do k = 2, nlayers
         do I1 = 1,npoin

            x = coord(1,I1)/1.0e3
            y = coord(2,I1)/1.0e3

            ! if((650.0 <= y .and. y <= Ly/1.0e3) .and. (400.0 <= x .and. x <= 500.0)) then
            if(x <= 200.0) then
               ! print*, x
               z_init(I1,k) = 200.0 !max(-100.0, z_interface(I1,k))
            end if

         end do
      end do 

      ! Layer densities reciprocal (1/rho)
      rho_0 = 1027.01037
      alpha(1) = 1.0/rho_0

      do k = 2,nlayers
         alpha(k) = 1.0/(rho_0 + k*0.2110/real(nlayers))
      end do

   case("seamount") ! lake at rest (well-balanced) test 

      gravity = 9.806
      H_bot = 4000.0   ! total depth

      ! Bottom topography (non-flat)
      zbot_df(:) = - H_bot

      xm=0.5*(xdims(1)+xdims(2))
      yl=0.5*(ydims(1)+ydims(2))
      L = 1.0/20.0e3
      delta = 0.4998

      do I1 = 1,npoin

         x = coord(1,I1)
         y = coord(2,I1)
         r = (L*(x-xm))**2 !- (L*(y-ym))**2

         zbot_df(I1) = zbot_df(I1)*(1.0 - delta * exp(-r))
      end do

      ! layer thicknesses (m) at global equilibrium
      do k = 1, nlayers
          layer_dz_eq(k) = (k)*H_bot/real(nlayers)
      end do

      ! Layer densities reciprocal (1/rho)
      rho_0 = 1027.01037
      alpha(1) = 1.0/rho_0

      do k = 2,nlayers
         alpha(k) = 1.0/(rho_0 + k*0.2110/real(nlayers))
      end do
   case default
      print*, "Unknown test case in cube initialization ", test_case
      stop
      
   end select 

   ! Layer interface at the equilibrium state
   do k = 2, nlayers+1
      z_interface_equil(k) = z_interface_equil(k-1) - layer_dz_eq(k-1)
   end do

   ! Making sure the layer interface is not beyond the bottom depth
   do I1 = 1,npoin
      do k = 1, nlayers+1
         z_interface(I1,k) = max(zbot_df(I1), z_interface_equil(k))
      end do
   end do 

   do I1 = 1,npoin
      do k = nlayers, 1, -1
         z_interface(I1,k) = max(z_interface(I1,k), z_interface(I1,k+1) + dry_cutoff)
      end do
   end do

   ! If a horizontal interface intersects sloping bottom topography
   ! in the interior of a grid element, then modify the interface
   ! in that element so that it intersects the topography at 
   ! the element edge that is above  z_interface_equil(k).  
   ! The value  z_init_flag(e,k) = 0  indicates that this action
   ! is taken in element  e  for interface  k,  0 <= k <= nk-1.  
   ! At such locations, the initial perturbation  z_init  should be zero,
   ! otherwise,  z_init_flag(e,k) = 1.  

   z_init_flag(:,:) = 1
   z_init_flag_elem(:,:) = 1

   do e = 1, nelem

      ! x-direction edges of element e
      do j = 1,ngl
         
         I1 = intma(1,j,1,e)
         I2 = intma(ngl,j,1,e)
         dx = coord(1,I2) - coord(1,I1)
         dzbot_dx = (zbot_df(I2) - zbot_df(I1))/dx

         do k = 2, nlayers
            ztmp = z_interface_equil(k)
            if ((zbot_df(I1)- ztmp) * (zbot_df(I2) - ztmp) < 0.0) then
               xl = coord(1,I1)
               xr = coord(1,I2)
               xmid = 0.5*(xl + xr)
               zl = max(zbot_df(I1), ztmp)
               zr = max(zbot_df(I2), ztmp)

               z_init_flag_elem(e,k) = 0

               do i = 1, ngl
                  ip = intma(i,j,1,e)
                  z_init_flag(ip,k) = 0

                  x = coord(1,ip)
                  L = zl*(xr-x)/dx + zr*(x-xl)/dx
                  dLdx = (zr - zl)/dx

                  if (zbot_df(I1) > ztmp) then
                     slope_init_l = dzbot_dx
                     slope_init_r = 0.0
                  else
                     slope_init_l = 0.0
                     slope_init_r = dzbot_dx
                  end if
                  slope_edge_adjust = min(dLdx - slope_init_l, slope_init_r - dLdx)
                  z_interface(ip,k) = L + slope_edge_adjust*((x - xmid)**2 - 0.25*(dx**2))/(dx)
                  ! z_interface(ip,k) = max(z_interface(ip,k), z_interface(ip,k-1))
               end do
            end if
         end do
      end do

      ! y-direction edges of element e
      do i = 1,ngl
         I1 = intma(i,1,1,e)
         I2 = intma(i,ngl,1,e)
         dy = coord(2,I2) - coord(2,I1)
         dzbot_dy = (zbot_df(I2) - zbot_df(I1))/dy

         do k = 2, nlayers
            ztmp = z_interface_equil(k)
            if ((zbot_df(I1)- ztmp) * (zbot_df(I2) - ztmp) < 0.0) then
               yl = coord(2,I1)
               yr = coord(2,I2)
               ymid = 0.5*(yl + yr)
               zl = max(zbot_df(I1), ztmp)
               zr = max(zbot_df(I2), ztmp)

               do j = 1, ngl
                  ip = intma(i,j,1,e)
                  z_init_flag(ip,k) = 0

                  y = coord(2,ip)
                  L = zl*(yr-y)/dy + zr*(y-yl)/dy
                  dLdy = (zr - zl)/dy

                  if (zbot_df(I1) > ztmp) then
                     slope_init_l = dzbot_dy
                     slope_init_r = 0.0
                  else
                     slope_init_l = 0.0
                     slope_init_r = dzbot_dy
                  end if
                  slope_edge_adjust = min(dLdy - slope_init_l, slope_init_r - dLdy)
                  z_interface(ip,k) = L + slope_edge_adjust*((y - ymid)**2 - 0.25*(dy**2))/(dy)
                  ! z_interface(ip,k) = max(z_interface(ip,k), z_interface(ip,k-1))
               end do
            end if
         end do

      end do
   end do

   z_interface_initial = z_interface


   ! Compute the degrees of freedom for  p'_b  in each cell.
   ! In each cell, the computation of degrees of freedom uses
   ! orthogonal projections onto the space of polynomials being used.

   do k = 1, nlayers
      pbprime_df(:) = pbprime_df(:) + (gravity/alpha(k))*(z_interface(:,k) - z_interface(:,k+1))
   end do

   ! Compute  p'_b  at the faces and quadrature points 
   ! in each cell.  This is the value of  
   ! pb = vertical sum of  Delta p  over all layers
   ! at the reference state.

   call interpolate_pbprime_init(pbprime_df_face, pbprime_df)

   ! Use the initial perturbations in array  z_init  
   ! to modify the elevations in array  z_interface  so that 
   ! array  z_interface  then refers to the specified initial state.  
   ! The entries in array  z_init_flag  are either 1 or 0, and they
   ! are defined above;  value  0  is used in situations where the
   ! initial perturbation should be automatically zero.

   do I1 = 1, npoin
      do k = 1, nlayers
         z_interface(I1,k) = z_interface(I1,k) + z_init(I1,k)*z_init_flag(I1,k)
      end do
   end do

   ! Compute pointwise values of  Delta p  at the quadrature points
   ! in each cell, for each layer.
   do k = 1, nlayers
      q_df(1,:,k) = (gravity/alpha(k))*(z_interface(:,k) - z_interface(:,k+1))
   end do

   ! print*, 'z_interface(1,k) = ', z_interface(1,1) - z_interface(1,2), z_interface(1,1), z_interface(1,2)
   ! print*, 'z_interface(1,k) = ', z_interface(1,2) - z_interface(1,3), z_interface(1,2), z_interface(1,3)
   ! print*, q_df(1,1,:)
   ! stop

   ! Compute dofs for u*(Delta p) and v*(Delta p)
   do k = 1, nlayers
      q_df(2,:,k) = u_df(:,k)*q_df(1,:,k)
      q_df(3,:,k) = v_df(:,k)*q_df(1,:,k)
   end do

   call poslimiter(q_df,alpha)

   ! print*, "q_df:", q_df(2,:,2)
   ! stop

   ! === Barotropic variables ===

   qb_df = 0.0
   
   ! Compute degrees of freedom for the barotropic mass and momentum
   ! dependent variables. These are the vertical sums of the
   ! degrees of freedom for the corresponding layer variables.
   ! Also compute degrees of freedom for pbpert = pb - pbprime_init.
   
   do k = 1, nlayers
      qb_df(1,:) = qb_df(1,:) + q_df(1,:,k)
      qb_df(3,:) = qb_df(3,:) + q_df(2,:,k)
      qb_df(4,:) = qb_df(4,:) + q_df(3,:,k)
   end do
   qb_df(2,:) = qb_df(1,:) - pbprime_df(:)
   
end subroutine initial_conditions

!--------------------------------------------------
!>brief Initial grid dimensions and BCs for cubes
!--------------------------------------------------
subroutine initial_grid_coord()

    use mod_input, only: xdims, ydims, x_boundary, y_boundary, z_boundary

    implicit none
    
    integer iboundary(6)
    real :: xmin, xmax, ymin, ymax

    iboundary(1:2)=y_boundary(1:2)
    iboundary(3:4)=x_boundary(1:2)

    xmin=xdims(1) ;  xmax=xdims(2)
    ymin=ydims(1) ;  ymax=ydims(2)

end subroutine initial_grid_coord
