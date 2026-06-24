!-----------------------------------------------------------------!
!>@brief This subroutine builds the Initial Conditions for the multilayer Shallow Water Equations
!>@author Written by Yao Gahounzo
!>@date March 27, 2023
!-----------------------------------------------------------------!

subroutine initial_conditions(inp, G, b, mf, init, mt)

   use mod_input,         only: input
   use mod_grid,          only: grid
   use mod_basis,         only: basis
   use mod_initial,       only: initial
   use mod_constants,     only: gravity, pi
   use mod_initial_mlswe, only: interpolate_pbprime_init, poslimiter
   use mod_metrics,       only: metrics
   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION
   use mod_face,    only: face_CS

   implicit none

   type(input),   intent(in)           :: inp
   type(grid),    intent(in)           :: G
   type(basis),   intent(in)           :: b
   type(face_CS), intent(in)           :: mf
   type(initial), intent(inout)        :: init
   type(metrics), intent(in), optional :: mt
   
   real, dimension(G%npoin, inp%nlayers) :: u_df, v_df

   integer :: e, m, n, I1, iquad, Iq, I2, i, j
   real :: x, y, xmid, ymid, L, amp, r, rho_0, H_bot
   integer :: ierr, iface, ilr, k,ip, il, jl, kl, el
   real :: xmax, xmin, ymax, ymin, zmin, zmax, yl, xm
   real :: xmax_l, xmin_l, ymax_l, ymin_l, zmin_l, zmax_l, Ly
   real :: delta, dp_dry
   real :: c_jump, indep(inp%nlayers+1)
   real, dimension(G%npoin, inp%nlayers+1) :: z_init
   real :: z_interface_equil(inp%nlayers+1), layer_dz_eq(inp%nlayers)
   real :: xl, xr, ztmp, zl, zr, slope_edge_adjust, slope_init_l, slope_init_r
   real :: dLdx, dzbot_dx, dx, dzbot_dy, dy, yr, dLdy
   
   xmin_l=minval(G%coord(1,:)); xmax_l=maxval(G%coord(1,:))
   ymin_l=minval(G%coord(2,:)); ymax_l=maxval(G%coord(2,:))
   zmin_l=minval(G%coord(3,:)); zmax_l=maxval(G%coord(3,:))

   call mpi_allreduce(xmin_l,xmin,1,MPI_PRECISION,mpi_min,mpi_comm_world,ierr)
   call mpi_allreduce(xmax_l,xmax,1,MPI_PRECISION,mpi_max,mpi_comm_world,ierr)
   call mpi_allreduce(ymin_l,ymin,1,MPI_PRECISION,mpi_min,mpi_comm_world,ierr)
   call mpi_allreduce(ymax_l,ymax,1,MPI_PRECISION,mpi_max,mpi_comm_world,ierr)
   call mpi_allreduce(zmin_l,zmin,1,MPI_PRECISION,mpi_min,mpi_comm_world,ierr)
   call mpi_allreduce(zmax_l,zmax,1,MPI_PRECISION,mpi_max,mpi_comm_world,ierr)

   z_init = 0.0
   init%z_interface = 0.0
   init%pbprime_df = 0.0
   init%q_df = 0.0
   u_df = 0.0
   v_df = 0.0
   init%tau_wind_df = 0.0
   init%pbprime_df_face = 0.0
   z_interface_equil = 0.0

   init%kvector(1,:)=0.0
   init%kvector(2,:)=0.0
   init%kvector(3,:)=1.0

   Ly = inp%ydims(2)-inp%ydims(1)

   ! The following is for a two-layer configuration.
   ! The free surface is level, the interface is one cosine (MAK test)
   ! hump centered in the interval, and u = v = 0.
   ! The motion is mostly internal.

   select case (trim(inp%test_case))
      
   case ("bump") ! bump test (wave propagation)

      gravity = 9.806
      H_bot = 40.0   ! total depth

      ! Bottom topography (flat)
      init%zbot_df(:) = - H_bot
      layer_dz_eq(:) = H_bot/real(inp%nlayers) ! layer thicknesses (m) at global equilibrium

      ! Bump centered at (xm,yl) with radius L and amplitude amp at the second layer interface
      xm=0.5*(xmax+xmin)
      yl=0.5*(ymax+ymin)

      L = 250.0
      amp = 1.0

      do I1 = 1, G%npoin

         x = G%coord(1,I1)
         y = G%coord(2,I1)
         r = sqrt((x-xm)**2 + (y-yl)**2)

         if (r < L) then
            z_init(I1,2) =  0.5*amp*(1.0 + cos(pi*r/L))
         end if
      end do

      ! Layer densities reciprocal (1/rho)
      init%alpha_mlswe(1) = 0.9737e-3
      init%alpha_mlswe(2) = 0.9735e-3

   case("lakeAtrest") ! lake at rest (well-balanced) test 

      gravity = 9.806
      H_bot = 40.0   ! total depth

      ! Bottom topography (non-flat)
      init%zbot_df(:) = - H_bot

      xm=0.5*(inp%xdims(1)+inp%xdims(2))
      yl=0.5*(inp%ydims(1)+inp%ydims(2))
      L = 250.0

      do I1 = 1,G%npoin

         x = G%coord(1,I1)
         y = G%coord(2,I1)
         r = sqrt((x-xm)**2 + (y-yl)**2)

         if (r < L) then
            init%zbot_df(I1) = init%zbot_df(I1) + 3.0*(1.0 + cos(pi*r/L))
         end if
      end do  

      ! layer thicknesses (m) at global equilibrium
      if(inp%nlayers < 5) then
         layer_dz_eq(:) = H_bot/real(inp%nlayers)
      else
         layer_dz_eq(1:inp%nlayers-1) = 32.0/real(inp%nlayers-1)
         layer_dz_eq(inp%nlayers) = H_bot - 32.0
      endif

      ! Layer densities reciprocal (1/rho)
      rho_0 = 1027.01037
      init%alpha_mlswe(1) = 1.0/rho_0

      do k = 2,inp%nlayers
         init%alpha_mlswe(k) = 1.0/(rho_0 + k*0.2110/real(inp%nlayers))
      end do

   case ("double-gyre") ! double-gyre 

      gravity = 9.806
      H_bot = 9928.0  ! total depth

      ! Bottom topography (flat)
      init%zbot_df(:) = - H_bot

      ! layer thicknesses (m) at global equilibrium
      layer_dz_eq(1) = 1489.5
      layer_dz_eq(2) = H_bot - 1489.5

      ! Layer densities reciprocal (1/rho)
      init%alpha_mlswe(1) = 9.7370e-04
      init%alpha_mlswe(2) = 9.7350e-04

      ! wind stress
      do I1 = 1, G%npoin
         y = G%coord(2,I1)
         init%tau_wind_df(1,I1) = -0.1*cos(2.0*pi*y/Ly)
      end do

   case ("dam") ! dam-break test

      gravity = 9.806
      H_bot = 2000.0 ! total depth

      ! Bottom topography
      do I1 = 1, G%npoin
        
         x = G%coord(1,I1)/1.0e3
         y = G%coord(2,I1)/1.0e3

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

         ! init%zbot_df(I1) = 500.0 + 0.5*(H_bot - 500.0)*(1.0 + tanh((x - 401.0)/60.0))
         init%zbot_df(I1) = 500.0 + 0.5*(H_bot - 500.0)*(1.0 + tanh((x - 60.0)/20.0))

         init%zbot_df(I1) = - init%zbot_df(I1)
      end do
        
      ! layer thicknesses (m) at global equilibrium: equal-thickness layers
      layer_dz_eq(:) = H_bot / real(inp%nlayers)

      do k = 2, inp%nlayers
         do I1 = 1,G%npoin

            x = G%coord(1,I1)/1.0e3
            y = G%coord(2,I1)/1.0e3

            ! if((650.0 <= y .and. y <= Ly/1.0e3) .and. (400.0 <= x .and. x <= 500.0)) then
            if(x <= 20.0) then
               ! print*, x
               z_init(I1,k) = 400.0 !max(-100.0, z_interface(I1,k))
            end if

         end do
      end do 

      ! Layer densities reciprocal (1/rho)
      rho_0 = 1027.01037
      init%alpha_mlswe(1) = 1.0/rho_0

      do k = 2,inp%nlayers
         init%alpha_mlswe(k) = 1.0/(rho_0 + k*0.2110/real(inp%nlayers))
      end do

   case("seamount") ! lake at rest (well-balanced) test 

      gravity = 9.806
      H_bot = 4000.0   ! total depth

      ! Bottom topography (non-flat)
      init%zbot_df(:) = - H_bot

      xm=0.5*(inp%xdims(1)+inp%xdims(2))
      yl=0.5*(inp%ydims(1)+inp%ydims(2))
      L = 1.0/20.0e3
      delta = 0.4998

      do I1 = 1,G%npoin

         x = G%coord(1,I1)
         y = G%coord(2,I1)
         r = (L*(x-xm))**2 !- (L*(y-ym))**2

         init%zbot_df(I1) = init%zbot_df(I1)*(1.0 - delta * exp(-r))
      end do

      ! layer thicknesses (m) at global equilibrium
      do k = 1, inp%nlayers
          layer_dz_eq(k) = (k)*H_bot/real(inp%nlayers)
      end do

      ! Layer densities reciprocal (1/rho)
      rho_0 = 1027.01037
      init%alpha_mlswe(1) = 1.0/rho_0

      do k = 2,inp%nlayers
         init%alpha_mlswe(k) = 1.0/(rho_0 + k*0.2110/real(inp%nlayers))
      end do
   case default
      print*, "Unknown test case in cube initialization ", inp%test_case
      stop
      
   end select 

   ! Layer interface at the equilibrium state
   do k = 2, inp%nlayers+1
      z_interface_equil(k) = z_interface_equil(k-1) - layer_dz_eq(k-1)
   end do

   ! Making sure the layer interface is not beyond the bottom depth
   do I1 = 1,G%npoin
      do k = 1, inp%nlayers+1
         init%z_interface(I1,k) = max(init%zbot_df(I1), z_interface_equil(k))
      end do
   end do 

   do I1 = 1,G%npoin
      do k = inp%nlayers, 1, -1
         init%z_interface(I1,k) = max(init%z_interface(I1,k), init%z_interface(I1,k+1))
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

   init%z_init_flag(:,:) = 1
   init%z_init_flag_elem(:,:) = 1

   ! do e = 1, G%nelem

   !    ! x-direction edges of element e
   !    do j = 1, b%ngl

   !       I1 = G%intma(1,j,1,e)
   !       I2 = G%intma(b%ngl,j,1,e)
   !       dx = G%coord(1,I2) - G%coord(1,I1)
   !       dzbot_dx = (init%zbot_df(I2) - init%zbot_df(I1))/dx

   !       do k = 2, inp%nlayers
   !          ztmp = z_interface_equil(k)
   !          if ((init%zbot_df(I1) - ztmp) * (init%zbot_df(I2) - ztmp) < 0.0) then
   !             xl = G%coord(1,I1)
   !             xr = G%coord(1,I2)
   !             xmid = 0.5*(xl + xr)
   !             zl = max(init%zbot_df(I1), ztmp)
   !             zr = max(init%zbot_df(I2), ztmp)

   !             init%z_init_flag_elem(e,k) = 0

   !             do i = 1, b%ngl
   !                ip = G%intma(i,j,1,e)
   !                init%z_init_flag(ip,k) = 0

   !                x = G%coord(1,ip)
   !                L = zl*(xr-x)/dx + zr*(x-xl)/dx
   !                dLdx = (zr - zl)/dx

   !                if (init%zbot_df(I1) < ztmp) then
   !                   slope_init_l = 0.0
   !                   slope_init_r = dzbot_dx
   !                else
   !                   slope_init_l = dzbot_dx
   !                   slope_init_r = 0.0
   !                end if
   !                slope_edge_adjust = min(dLdx - slope_init_l, slope_init_r - dLdx)
   !                init%z_interface(ip,k) = L + slope_edge_adjust*((x - xmid)**2/dx - 0.25*dx)
   !                init%z_interface(ip,k) = min(init%z_interface(ip,k), init%z_interface(ip,k-1))
   !             end do
   !          end if
   !       end do
   !    end do

   !    ! y-direction edges of element e — skip nodes already set by the x-sweep
   !    do i = 1, b%ngl
   !       I1 = G%intma(i,1,1,e)
   !       I2 = G%intma(i,b%ngl,1,e)
   !       dy = G%coord(2,I2) - G%coord(2,I1)
   !       dzbot_dy = (init%zbot_df(I2) - init%zbot_df(I1))/dy

   !       do k = 2, inp%nlayers
   !          ztmp = z_interface_equil(k)
   !          if ((init%zbot_df(I1) - ztmp) * (init%zbot_df(I2) - ztmp) < 0.0) then
   !             yl = G%coord(2,I1)
   !             yr = G%coord(2,I2)
   !             ymid = 0.5*(yl + yr)
   !             zl = max(init%zbot_df(I1), ztmp)
   !             zr = max(init%zbot_df(I2), ztmp)

   !             init%z_init_flag_elem(e,k) = 0

   !             do j = 1, b%ngl
   !                ip = G%intma(i,j,1,e)
   !                if (init%z_init_flag(ip,k) == 1) then
   !                   init%z_init_flag(ip,k) = 0

   !                   y = G%coord(2,ip)
   !                   L = zl*(yr-y)/dy + zr*(y-yl)/dy
   !                   dLdy = (zr - zl)/dy

   !                   if (init%zbot_df(I1) < ztmp) then
   !                      slope_init_l = 0.0
   !                      slope_init_r = dzbot_dy
   !                   else
   !                      slope_init_l = dzbot_dy
   !                      slope_init_r = 0.0
   !                   end if
   !                   slope_edge_adjust = min(dLdy - slope_init_l, slope_init_r - dLdy)
   !                   init%z_interface(ip,k) = L + slope_edge_adjust*((y - ymid)**2/dy - 0.25*dy)
   !                   init%z_interface(ip,k) = min(init%z_interface(ip,k), init%z_interface(ip,k-1))
   !                end if
   !             end do
   !          end if
   !       end do

   !    end do
   ! end do

   init%z_interface_initial = init%z_interface


   ! Compute the degrees of freedom for  p'_b  in each cell.
   ! In each cell, the computation of degrees of freedom uses
   ! orthogonal projections onto the space of polynomials being used.

   do k = 1, inp%nlayers
      init%pbprime_df(:) = init%pbprime_df(:) + (gravity/init%alpha_mlswe(k))*(init%z_interface(:,k) - init%z_interface(:,k+1))
   end do

   ! Compute  p'_b  at the faces and quadrature points 
   ! in each cell.  This is the value of  
   ! pb = vertical sum of  Delta p  over all layers
   ! at the reference state.

   call interpolate_pbprime_init(init%pbprime_df_face, init%pbprime_df, b, G, inp, mf)

   ! Use the initial perturbations in array  z_init  
   ! to modify the elevations in array  z_interface  so that 
   ! array  z_interface  then refers to the specified initial state.  
   ! The entries in array  z_init_flag  are either 1 or 0, and they
   ! are defined above;  value  0  is used in situations where the
   ! initial perturbation should be automatically zero.

   do I1 = 1, G%npoin
      do k = 1, inp%nlayers
         init%z_interface(I1,k) = init%z_interface(I1,k) + z_init(I1,k)*init%z_init_flag(I1,k)
      end do
   end do

   ! Compute pointwise values of  Delta p  at the quadrature points
   ! in each cell, for each layer.
   do k = 1, inp%nlayers
      init%q_df(1,:,k) = (gravity/init%alpha_mlswe(k))*(init%z_interface(:,k) - init%z_interface(:,k+1))
   end do

   ! print*, 'z_interface(1,k) = ', z_interface(1,1) - z_interface(1,2), z_interface(1,1), z_interface(1,2)
   ! print*, 'z_interface(1,k) = ', z_interface(1,2) - z_interface(1,3), z_interface(1,2), z_interface(1,3)
   ! print*, q_df(1,1,:)
   ! stop

   ! Compute dofs for u*(Delta p) and v*(Delta p)
   do k = 1, inp%nlayers
      init%q_df(2,:,k) = u_df(:,k)*init%q_df(1,:,k)
      init%q_df(3,:,k) = v_df(:,k)*init%q_df(1,:,k)
   end do

   if (inp%lposlimiter .and. present(mt)) call poslimiter(b, G, inp, mt, init%q_df, init%alpha_mlswe)

   ! print*, "q_df:", q_df(2,:,2)
   ! stop

   ! === Barotropic variables ===

   init%qb_df = 0.0
   
   ! Compute degrees of freedom for the barotropic mass and momentum
   ! dependent variables. These are the vertical sums of the
   ! degrees of freedom for the corresponding layer variables.
   ! Also compute degrees of freedom for pbpert = pb - pbprime_init.
   
   do k = 1, inp%nlayers
      init%qb_df(1,:) = init%qb_df(1,:) + init%q_df(1,:,k)
      init%qb_df(3,:) = init%qb_df(3,:) + init%q_df(2,:,k)
      init%qb_df(4,:) = init%qb_df(4,:) + init%q_df(3,:,k)
   end do
   init%qb_df(2,:) = init%qb_df(1,:) - init%pbprime_df(:)
   
end subroutine initial_conditions

!--------------------------------------------------
!>brief Initial grid dimensions and BCs for cubes
!--------------------------------------------------
subroutine initial_grid_coord(inp)

    use mod_input, only: input

    implicit none
    
    type(input), intent(in) :: inp
    integer iboundary(6)
    real :: xmin, xmax, ymin, ymax

    iboundary(1:2)=inp%y_boundary(1:2)
    iboundary(3:4)=inp%x_boundary(1:2)

    xmin=inp%xdims(1) ;  xmax=inp%xdims(2)
    ymin=inp%ydims(1) ;  ymax=inp%ydims(2)

end subroutine initial_grid_coord
