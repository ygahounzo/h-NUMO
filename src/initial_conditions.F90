!-----------------------------------------------------------------!
!>@brief This subroutine builds the Initial Conditions for the multilayer Shallow Water Equations
!>@author Written by Yao Gahounzo
!>@date March 27, 2023
!-----------------------------------------------------------------!

subroutine initial_conditions(q_df, pbprime_df, qb_df, alpha, pbprime_df_face, zbot_df, &
   tau_wind_df, z_interface)
    

   use mod_grid, only: nelem, nface, npoin_q, npoin, coord, intma_dg_quad, face
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
   real, dimension(npoin, nlayers+1) :: z_init, z_interface
   real, dimension(npoin,nlayers) :: height_layer
   real :: c_jump, indep(nlayers+1)

   integer :: e, m, n, I1, iquad, Iq
   real :: x, y, xmid, ymid, L, amp, r, rho_0, H_bot
   integer :: ierr, iface, ilr, k,ip, il, jl, kl, el
   real :: xmax, xmin, ymax, ymin, zmin, zmax, yl, xm
   real :: xmax_l, xmin_l, ymax_l, ymin_l, zmin_l, zmax_l, Ly
   real :: delta, dp_dry
   
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

      ! Layer interface
      do k = 1, nlayers+1
         z_interface(:,k) = -(k-1)*H_bot/real(nlayers)
      end do

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
            z_interface(I1,2) = z_interface(I1,2) + 0.5*amp*(1.0 + cos(pi*r/L))
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

      ! Layer interface
      do k = 1, nlayers+1
         if(nlayers < 5) then
            z_interface(:,k) = -(k-1)*H_bot/real(nlayers)
         else
            z_interface(:,k) = -(k-1)*32/real(nlayers-1)
            z_interface(:,nlayers+1) = -H_bot
         endif
      end do

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

      ! Layer interface
      z_interface(:,2) = -1489.5
      z_interface(:,3) = - H_bot

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

         zbot_df(I1) = 500.0 + 0.5*(H_bot - 500.0)*(1.0 + tanh((x - 201.0)/60.0))

         zbot_df(I1) = - zbot_df(I1)
      end do
        
      do k = 2, nlayers
        indep(k) = H_bot*(real(k-1)-0.5)/real(nlayers-1)
        !print*, 'k = ', k, indep(k)
        !indep(k) = 200.0
      enddo

      ! Layer interface
      do k = 1, nlayers
        z_interface(:,k) = -indep(k)
      enddo
      z_interface(:,nlayers+1) = zbot_df(:)

      ! do k = 1, nlayers
      !     do I1 = 1,npoin
      !        !if(abs(z_interface(I1,k)) >= abs(zbot_df(I1))) then
      !        !   z_interface(I1,k) = zbot_df(I1)
      !        !endif
      !        z_interface(I1,k) = max(zbot_df(I1), z_interface(I1,k))
      !     end do
      ! end do

      do k = 2, nlayers
         do I1 = 1,npoin

            x = coord(1,I1)/1.0e3
            y = coord(2,I1)/1.0e3

            ! if((650.0 <= y .and. y <= Ly/1.0e3) .and. (400.0 <= x .and. x <= 500.0)) then
            if(x <= 20.0) then
               ! print*, x
               z_interface(I1,k) = -200.0 !max(-100.0, z_interface(I1,k))
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

      ! Layer interface
      do k = 1, nlayers+1
          z_interface(:,k) = -(k-1)*H_bot/real(nlayers)
      end do
      z_interface(:,nlayers+1) = zbot_df(:)

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

   ! Making sure the layer interface is not beyond the bottom depth

   do I1 = 1,npoin
      do k = 1, nlayers+1
         !if(abs(z_interface(I1,k)) >= abs(zbot_df(I1))) then
         !    z_interface(I1,k) = zbot_df(I1)
         !endif 
         z_interface(I1,k) = max(zbot_df(I1), z_interface(I1,k))
      end do
   end do 

   !print*, minval(zbot_df(:)), maxval(zbot_df(:))
   !do k = 1, nlayers+1
   !     print*, minval(z_interface(:,k)), maxval(z_interface(:,k))
   !enddo

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

   ! Compute pointwise values of  Delta p  at the quadrature points
   ! in each cell, for each layer.
   do k = 1, nlayers
      q_df(1,:,k) = (gravity/alpha(k))*(z_interface(:,k) - z_interface(:,k+1))
   end do

   ! Compute dofs for u*(Delta p) and v*(Delta p)
   do k = 1, nlayers
      q_df(2,:,k) = u_df(:,k)*q_df(1,:,k)
      q_df(3,:,k) = v_df(:,k)*q_df(1,:,k)
   end do

   call poslimiter(q_df,alpha)

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
