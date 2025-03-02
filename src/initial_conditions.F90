!-----------------------------------------------------------------!
!>@brief This subroutine builds the Initial Conditions for the multilayer Shallow Water Equations
!>@author Written by Yao Gahounzo
!>@date March 27, 2023
!
!-----------------------------------------------------------------!

subroutine initial_conditions(q, qprime, q_df, pbprime_init, pbprime_df, q_face, qprime_face, pbprime_face, &
   one_over_pbprime, one_over_pbprime_face, pbprime_edge, one_over_pbprime_edge, dpprime_df, one_over_pbprime_df, &
   layer_dz_eq, qb, qb_face, qb_df, qprime_df, alpha, one_over_pbprime_df_face, zbot_df)
    

   use mod_grid, only: nelem, nface, npoin_q, npoin, coord

   use mod_constants, only: gravity, pi, tol, omega, earth_radius

   use mod_initial_mlswe, only: interpolate_from_dof_to_quad_uv_init, interpolate_pbprime_init

   use mod_basis, only: ngl, nq

   use mod_layer_terms, only: evaluate_dp, evaluate_dp_face, evaluate_mom, evaluate_mom_face

   use mod_initial, only: kvector

   use mod_input, only: gravity_in, &
      nelx, nelz, eqn_set, &
      xdims, ydims,nlayers, test_case

   use mod_types, only : r8
   
   use mpi
   
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none
   
   real, dimension(3,npoin_q, nlayers) :: q
   real, dimension(3,npoin_q, nlayers):: qprime
   real, dimension(3,npoin, nlayers) :: q_df
   real, dimension(3,npoin, nlayers) :: qprime_df
   real, dimension(npoin_q) :: pbprime_init
   real, dimension(npoin) :: pbprime_df
   real, dimension(3,2,nq, nface, nlayers) :: q_face
   real, dimension(3,2,nq, nface, nlayers) :: qprime_face
   real, dimension(2,nq,nface) :: pbprime_face
   real, dimension(npoin_q) :: one_over_pbprime
   real, dimension(2,nq,nface) :: one_over_pbprime_face
   real, dimension(2,ngl,nface) :: one_over_pbprime_df_face
   real, dimension(nq, nface) :: pbprime_edge
   real, dimension(nq, nface) :: one_over_pbprime_edge
   real, dimension(npoin,nlayers) :: dpprime_df
   real, dimension(npoin) :: one_over_pbprime_df, zbot_df
   real, dimension(nlayers) :: layer_dz_eq
   real, dimension(4,npoin_q) :: qb
   real, dimension(4,2,nq,nface) :: qb_face
   real, dimension(4,npoin) :: qb_df
   real, dimension(nlayers) :: alpha

   real, dimension(npoin, nlayers) :: u_df, v_df(npoin, nlayers)
   real, dimension(npoin) :: one_plus_eta_temp
   real, dimension(npoin_q) :: pb_temp, one_plus_eta_temp1

   real, dimension(npoin, nlayers+1) :: z_init, z_interface
   real :: c_jump

   integer :: e, m, n, I1, iquad, Iq
   real :: x, y, xmid, ymid, L, amp, r, rho_0
   integer :: ierr, iface, ilr, k,ip

   real xmax, xmin, ymax, ymin, zmin, zmax, yl, xm
   real xmax_l, xmin_l, ymax_l, ymin_l, zmin_l, zmax_l, Ly
   
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
   q = 0.0
   q_df = 0.0
   qprime = 0.0
   q_face = 0.0
   qprime_face = 0.0
   one_over_pbprime = 0.0
   one_over_pbprime_face = 0.0
   pbprime_edge = 0.0
   one_over_pbprime_edge = 0.0
   one_over_pbprime_df = 0.0
   u_df = 0.0
   v_df = 0.0

   kvector(1,:)=0.0
   kvector(2,:)=0.0
   kvector(3,:)=1.0

   ! The following is for a two-layer configuration.
   ! The free surface is level, the interface is one cosine (MAK test)
   ! hump centered in the interval, and u = v = 0.
   ! The motion is mostly internal.

   select case (trim(test_case))
      
   case ("bump") ! bump test (wave propagation)

      gravity = 9.806
       
      layer_dz_eq(1) = 20.0
      layer_dz_eq(2) = 20.0

      alpha(1) = 0.9737e-3
      alpha(2) = 0.9735e-3

      xm=0.5*(xmax+xmin)
      yl=0.5*(ymax+ymin)

      Ly = ymax - ymin

      L = 250.0
      amp = 1.0

      do I1 = 1, npoin

         x = coord(1,I1)
         y = coord(2,I1)

         r = sqrt((x-xm)**2 + (y-yl)**2)

         if (r < L) then
            z_init(I1,2) = 0.5*amp*(1.0 + cos(pi*r/L))
            ! z_init(I1,2) = -10.0*cos(2*pi*x/xmax)
         end if

         ! Bottom topography (flat)

         zbot_df(I1) = - sum(layer_dz_eq)

      end do

   case("lakeAtrest") ! lake at rest (well-balanced) test 

      gravity = 9.806

      rho_0 = 1027.01037
      alpha(1) = 1.0/rho_0

      if(nlayers < 5) then 
         layer_dz_eq(:) = 40.0/real(nlayers)
      else
         layer_dz_eq(:) = 32.0/real(nlayers-1)
         layer_dz_eq(nlayers) = 8.0
      endif

      do k = 2,nlayers
         alpha(k) = 1.0/(rho_0 + k*0.2110/real(nlayers))
      end do

      xm=0.5*(xdims(1)+xdims(2))
      yl=0.5*(ydims(1)+ydims(2))

      L = 250.0

      do I1 = 1,npoin 

         x = coord(1,I1)
         y = coord(2,I1)

         r = sqrt((x-xm)**2 + (y-yl)**2)

         ! Bottom topography (non-flat)
         zbot_df(I1) = - sum(layer_dz_eq)

         if (r < L) then
            zbot_df(I1) = zbot_df(I1) + 3.0*(1.0 + cos(pi*r/L))
         end if

      end do  
   case ("double-gyre") ! double-gyre 

      gravity = 9.806
       
      layer_dz_eq(1) = 1489.5
      layer_dz_eq(2) = 8438.5

      alpha(1) = 9.7370e-04
      alpha(2) = 9.7350e-04

      ! Bottom topography (flat)
      zbot_df(:) = - sum(layer_dz_eq)

   case default
      print*, "Unknown test case in cube initialization ", test_case
      stop
      
   end select 

   ! call interpolate_from_dof_to_quad_uv_init(q, q_df)

   do I1 = 1,npoin
      z_interface(I1,1) = 0.0
      do k = 1, nlayers
         z_interface(I1,k+1) = max(zbot_df(I1), z_interface(I1,k) - layer_dz_eq(k))
      end do
   end do 

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

   call interpolate_pbprime_init(pbprime_init,pbprime_face,one_over_pbprime_df_face, pbprime_df)

   ! Compute  1/p'_b  at the faces and quadrature points

   do iface = 1, nface

      do iquad = 1, nq

         pbprime_edge(iquad,iface) = pbprime_face(1,iquad,iface)

         if(pbprime_edge(iquad,iface) > 0.0) then
            one_over_pbprime_edge(iquad,iface) = 1.0/pbprime_edge(iquad,iface)
         end if

         do ilr = 1,2
            if(pbprime_face(ilr,iquad,iface) > 0.0) then
               one_over_pbprime_face(ilr,iquad,iface) = 1.0/pbprime_face(ilr,iquad,iface)
            end if
         end do
         
      end do
   end do

   do I1 = 1, npoin

      if(pbprime_df(I1) > 0.0) then
         one_over_pbprime_df(I1) = 1.0/pbprime_df(I1)
      end if
   end do

   do Iq = 1, npoin_q

      if(pbprime_init(Iq) > 0.0) then
         one_over_pbprime(Iq) = 1.0/pbprime_init(Iq)
      end if
   end do

   ! Use the initial perturbations in array  z_init  
   ! to modify the elevations in array  z_interface  so that 
   ! array  z_interface  then refers to the specified initial state.  

   do k = 1, nlayers+1
      z_interface(:,k) = z_interface(:,k) + z_init(:,k)
   end do

   ! Compute pointwise values of  Delta p  at the quadrature points
   ! in each cell, for each layer.

   do k = 1, nlayers
      q_df(1,:,k) = (gravity/alpha(k))*(z_interface(:,k) - z_interface(:,k+1))
      one_plus_eta_temp(:) = one_plus_eta_temp(:) + q_df(1,:,k)/pbprime_df(:)
   end do
   
   ! Compute pointwise values of  Delta p  at the dofs

   do k = 1, nlayers
      dpprime_df(:,k) = q_df(1,:,k)/one_plus_eta_temp(:) 
      qprime_df(1,:,k) = q_df(1,:,k)/one_plus_eta_temp(:) 
   end do

   ! Compute pointwise values of  Delta p  at the faces and quadrature points

   call evaluate_dp(q,qprime,q_df,pbprime_init)
   call evaluate_dp_face(q_face, qprime_face,q,qprime)

   ! Compute dofs for u*(Delta p) and v*(Delta p)

   do k = 1, nlayers
      q_df(2,:,k) = u_df(:,k)*q_df(1,:,k)
      q_df(3,:,k) = v_df(:,k)*q_df(1,:,k)
   end do

   ! Compute pointwise values of  u*(Delta p)  and  v*(Delta p) at the faces and quadrature points

   call evaluate_mom(q,q_df)
   call evaluate_mom_face(q_face, q)

   ! === Barotropic variables ===

   qb = 0.0
   qb_face = 0.0
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
   
   ! Compute pointwise values of the barotropic mass and momentum
   ! dependent variables at quadrature points and endpoints of each
   ! grid cell. Use vertical sums of the pointwise values of the
   ! layer variables.
   ! Also compute pointwise values of pbpert = pb - pbprime_init.

   do k = 1, nlayers
      qb(1,:) = qb(1,:) + q(1,:,k)
      qb(3,:) = qb(3,:) + q(2,:,k)
      qb(4,:) = qb(4,:) + q(3,:,k)
   end do
   qb(2,:) = qb(1,:) - pbprime_init(:)

   do iface = 1, nface
      do ilr = 1, 2 ! 1 = left, 2 = right
         do k = 1, nlayers
            qb_face(1,ilr,:,iface) = qb_face(1,ilr,:,iface) + q_face(1,ilr,:,iface,k)
            qb_face(3,ilr,:,iface) = qb_face(3,ilr,:,iface) + q_face(2,ilr,:,iface,k)
            qb_face(4,ilr,:,iface) = qb_face(4,ilr,:,iface) + q_face(3,ilr,:,iface,k)
         end do
         qb_face(2,ilr,:,iface) = qb_face(1,ilr,:,iface) - pbprime_face(ilr,:,iface)
      end do
   end do

   ! Compute values of u' and v' at quadrature points and face of each grid element

   do k = 1,nlayers
      qprime(2,:,k) = q(2,:,k)/q(1,:,k) - qb(3,:)
      qprime(3,:,k) = q(3,:,k)/q(1,:,k) - qb(4,:)
      qprime_df(2,:,k) = q_df(2,:,k)/q_df(1,:,k) - qb_df(3,:)/qb_df(1,:)
      qprime_df(3,:,k) = q_df(3,:,k)/q_df(1,:,k) - qb_df(4,:)/qb_df(1,:)
   end do

   do iface = 1, nface
      do ilr = 1,2
         do k = 1,nlayers
            qprime_face(2,ilr,:,iface,k) = q_face(2,ilr,:,iface,k)/q_face(1,ilr,:,iface,k) - qb_face(3,ilr,:,iface)/qb_face(1,ilr,:,iface)
            qprime_face(3,ilr,:,iface,k) = q_face(3,ilr,:,iface,k)/q_face(1,ilr,:,iface,k) - qb_face(4,ilr,:,iface)/qb_face(1,ilr,:,iface)
         end do
      end do
   end do

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
