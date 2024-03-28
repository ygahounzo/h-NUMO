!-----------------------------------------------------------------!
!>@brief This subroutine builds the Initial Conditions for the multilayer Shallow Water Equations
!>@author Written by Yao Gahounzo
!>@date 2018-03-01
!
!-----------------------------------------------------------------!

subroutine initial_conditions_mlswe(q_init, qprime, q_df, pbprime_init, pbprime_df, q_face, qprime_face, pbprime_face, &
   one_over_pbprime, one_over_pbprime_face, pbprime_edge, one_over_pbprime_edge, dpprime_df, one_over_pbprime_df, &
   layer_dz_eq, qb, qb_face, qb_df, alpha, coord)
    

   use mod_grid, only: nelem, nface, npoin_q, npoin

   use mod_constants, only: gravity, pi, tol, omega, earth_radius

   use mod_initial_mlswe, only: interpolate_from_dof_to_quad_uv_init, interpolate_pbprime_init

   use mod_basis, only: ngl, nq

   use mod_layer_terms, only: evaluate_dp, evaluate_dp_face, evaluate_udp_vdp, evaluate_udp_vdp_face, evaluate_dp_1l

   use mod_initial_mlswe, only: evaluate_dp2

   use mod_initial, only: kvector

   use mod_input, only: gravity_in, &
      nelx, nelz, eqn_set, &
      xdims, ydims,nlayers, icase

   !use mod_input, only: 

   use mod_ti, only: current_layer
   
   use mod_types, only : r8
   
   use mod_constants, only : gamma
   
   use mpi
   
   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   
   
   ! real, dimension(5,npoin_q, nlayers) :: q_init
   ! real, dimension(3,npoin_q, nlayers):: qprime
   ! real, dimension(3,npoin, nlayers) :: q_df
   ! real, dimension(npoin_q) :: pbprime_init
   ! real, dimension(npoin) :: pbprime_df
   ! real, dimension(5,2,nq, nface, nlayers) :: q_face
   ! real, dimension(3,2,nq, nface, nlayers) :: qprime_face
   ! real, dimension(2,nq,nface) :: pbprime_face
   ! real, dimension(npoin_q) :: one_over_pbprime
   ! real, dimension(2,nq,nface) :: one_over_pbprime_face
   ! real, dimension(nq, nface) :: pbprime_edge
   ! real, dimension(nq, nface) :: one_over_pbprime_edge
   ! real, dimension(npoin) :: dpprime_df
   ! real, dimension(npoin) :: one_over_pbprime_df
   ! real, dimension(nlayers) :: layer_dz_eq
   ! real, dimension(6,npoin_q) :: qb
   ! real, dimension(6,2,nq,nface) :: qb_face
   ! real, dimension(3,npoin) :: qb_df
   ! real, dimension(nlayers) :: alpha
   ! real, dimension(3,npoin) :: coord

   
   real(kind=r8) q_init(5,npoin_q,nlayers), q_df(3,npoin, nlayers), qprime(3,npoin_q, nlayers)
   real(kind=r8) q_face(5,2,nq, nface, nlayers), qprime_face(3,2,nq, nface, nlayers)
   real(kind=r8) one_over_pbprime(npoin_q), one_over_pbprime_face(2,nq,nface)
   real(kind=r8) pbprime_edge(nq, nface), one_over_pbprime_edge(nq, nface)
   real(kind=r8) one_over_pbprime_df(npoin)
   real(kind=r8) layer_dz_eq(nlayers)
   real(kind=r8) alpha(nlayers), pbprime_init(npoin_q), pbprime_df(npoin)
   real(kind=r8) dpprime_df(npoin), pbprime_face(2,nq,nface)

   real(kind=r8) qb(6,npoin_q), qb_face(6,2,nq,nface), qb_df(3,npoin)
   real(kind=r8) coord(3,npoin)
   ! real, dimension(6,2,nq,nface), intent(out) :: qb_face
   ! real, dimension(3,npoin), intent(out)      :: qb_df
   ! real, dimension(3,npoin), intent(in)      :: coord

   real :: u_df(npoin, nlayers), v_df(npoin, nlayers), one_plus_eta_temp(npoin)
   real, dimension(npoin_q) :: pb_temp, one_plus_eta_temp1

   real :: z_init(npoin, nlayers+1), z_interface(npoin, nlayers+1)
   integer :: e, m, n, I1, iquad, Iq
   real :: x, y, xmid, ymid, L, amp, r
   integer :: ierr, iface, ilr, k,ip

   real xmax, xmin, ymax, ymin, zmin, zmax, yl, xm
   real xmax_l, xmin_l, ymax_l, ymin_l, zmin_l, zmax_l, Ly
   real dp_temp(npoin_q,nlayers), qprime_temp(npoin_q,nlayers), quv(2,npoin_q,nlayers)
   


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
   q_init = 0.0
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

   kvector(1,:)=0
   kvector(2,:)=0
   kvector(3,:)=1

   ! The following is for a two-layer configuration.
   ! The free surface is level, the interface is one cosine (MAK test)
   ! hump centered in the interval, and u = v = 0.
   ! The motion is mostly internal.

   select case (icase)
      
   case (2023) !as 500, but two layers

      gravity = 9.806
       
      layer_dz_eq(1) = 20.0
      layer_dz_eq(2) = 20.0

      alpha(1) = 0.976e-3
      alpha(2) = 0.972e-3

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
         end if

      end do

   case (1000) !as 500, but two layers

      gravity = 9.806
       
      layer_dz_eq(1) = 100.0
      layer_dz_eq(2) = 900.0

      alpha(1) = 0.976e-3
      alpha(2) = 0.972e-3

      do I1 = 1, npoin

         z_init = 0.0
         u_df = 0.0
         v_df = 0.0

      end do

   case default
      print*, "Unknown case number in cube initialization ", icase
      stop
      
   end select


   call interpolate_from_dof_to_quad_uv_init(q_init, u_df, v_df)


   z_interface(:,1) = 0.0

   do k = 1, nlayers
      z_interface(:,k+1) = z_interface(:,k) - layer_dz_eq(k)
   end do

   ! Compute the degrees of freedom for  p'_b  in each cell.
   ! In each cell, the computation of degrees of freedom uses
   ! orthogonal projections onto the space of polynomials being used.

   do k = 1, nlayers

      pbprime_df = pbprime_df + (gravity/alpha(k))*(z_interface(:,k) - z_interface(:,k+1))
   end do

   ! Compute  p'_b  at the faces and quadrature points 
   ! in each cell.  This is the value of  
   ! pb = vertical sum of  Delta p  over all layers
   ! at the reference state.

   call interpolate_pbprime_init(pbprime_init,pbprime_face,pbprime_df)

   ! open(10,file='pbrime_face.txt')
   ! do iface = 1, nface
   !    write(10,*) iface, pbprime_face(1,:,iface)
   ! end do

   ! Compute  1/p'_b  at the faces and quadrature points

   do iface = 1, nface

      pbprime_edge(:,iface) = pbprime_face(1,:,iface)

      do iquad = 1, nq

         if(pbprime_edge(iquad,iface) > 0.0) then
            one_over_pbprime_edge(iquad,iface) = 1.0/pbprime_edge(iquad,iface)
         end if
      end do

      do ilr = 1,2
         do iquad = 1, nq

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

      do I1 = 1, npoin

         q_df(1,I1,k) = (gravity/alpha(k))*(z_interface(I1,k) - z_interface(I1,k+1))

         one_plus_eta_temp(I1) = one_plus_eta_temp(I1) + q_df(1,I1,k)/pbprime_df(I1)

      end do
      
   end do

   
   ! Compute pointwise values of  Delta p  at the dofs

   do k = 1, nlayers

      do I1 = 1, npoin

         dpprime_df(I1) = q_df(1,I1,k)/one_plus_eta_temp(I1)

      end do
      
   end do


   ! Compute pointwise values of  Delta p  at the faces and quadrature points

   ! do k = 1, nlayers

   !    call evaluate_dp_1l(q_init(1,:,k), q_df(1,:,k))

   !    pb_temp(:) = pb_temp(:) + q_init(1,:,k)

   ! end do

   ! open(11,file='ubar.dat')

   ! do I1 = 1, npoin_q
   !    write(11,*) I1, q_init(1,I1,:)
   ! end do
   ! close(11)
   
   ! one_plus_eta_temp1(:) = pb_temp(:) / pbprime_init(:)

   ! do k = 1, nlayers
   !    do I1 = 1, npoin_q
   !       qprime(1,I1,k) = q_init(1,I1,k) / one_plus_eta_temp1(I1)
   !    end do
   ! end do

   ! open(12,file='q_init.dat')

   ! do I1 = 1, npoin_q
   !    write(12,*) I1, q_init(1,I1,:)
   ! end do
   ! close(12)

   dp_temp = 0.0
   qprime_temp = 0.0
   call evaluate_dp(dp_temp,qprime_temp,q_df,pbprime_init)

   do k = 1, nlayers
      do I1 = 1, npoin_q
         q_init(1,I1,k) = dp_temp(I1,k)
         qprime(1,I1,k) = qprime_temp(I1,k)
      end do
   end do

   ! q_init(1,:,:) = dp_temp
   ! qprime(1,:,:) = qprime_temp
   !call evaluate_dp(q_init(1,:,:),qprime(1,:,:),q_df,pbprime_init)
   call evaluate_dp_face(q_face, qprime_face,q_init,qprime)

   ! print*, 'stop here intial mlswe'

   ! open(11,file='q_init.dat')
   ! do iface = 1, nface
   !    do I1 = 1, nq

   !       write(11,*) iface,q_face(1,1,I1,iface,1)

   !    end do
   ! end do

   ! close(11)

   ! print*, q_init(Iq,k)

   ! open(11,file='q_init.dat')

   ! do I1 = 1, npoin_q
   !    write(11,fmt='(E8.16)') I1, q_init(1,I1,:)
   ! end do
   ! close(11)

   ! open(15,file='ubar.dat')

   ! do I1 = 1, npoin_q
   !    write(15,fmt='(E8.16)') I1, qprime(1,I1,:)
   ! end do
   ! close(15)

   ! print*, 'stop here intial mlswe'

   !stop


   ! Compute dofs for u*(Delta p) and v*(Delta p)

   ! do k = 1, nlayers

   !    do ip = 1, npoin

   !       q_df(2,ip,k) = u_df(ip,k)*q_df(1,ip,k)
   !       q_df(3,ip,k) = v_df(ip,k)*q_df(1,ip,k)

   !    end do
      
   ! end do

   
   ! Compute pointwise values of  u*(Delta p)  and  v*(Delta p) at the faces and quadrature points

   ! call evaluate_udp_vdp(quv,q_df)
   ! q_init(4:5,:,:) = quv(1:2,:,:)
   ! call evaluate_udp_vdp_face(q_face,q_init)

   !open(11,file='q_init.dat')

   
   


   ! === Barotropic variables ===

   qb = 0.0
   qb_face = 0.0
   qb_df = 0.0
   
   ! Compute degrees of freedom for the barotropic mass and momentum
   ! dependent variables. These are the vertical sums of the
   ! degrees of freedom for the corresponding layer variables.
   ! Also compute degrees of freedom for pbpert = pb - pbprime_init.
   
   do I1 = 1, npoin
      do k = 1, nlayers
         qb_df(1,I1) = qb_df(1,I1) + q_df(1,I1,k)
         qb_df(2,I1) = qb_df(2,I1) + q_df(2,I1,k)
         qb_df(3,I1) = qb_df(3,I1) + q_df(3,I1,k)
      end do
      qb_df(2,I1) = qb_df(1,I1) - pbprime_df(I1)
   end do
   
   ! Compute pointwise values of the barotropic mass and momentum
   ! dependent variables at quadrature points and endpoints of each
   ! grid cell. Use vertical sums of the pointwise values of the
   ! layer variables.
   ! Also compute pointwise values of pbpert = pb - pbprime_init.


   !qb(1,:) = sum(q_init(1,:,:),dim=2)
   
   do I1 = 1, npoin_q
      do k = 1, nlayers
         qb(1,I1) = qb(1,I1) + q_init(1,I1,k)
         qb(5,I1) = qb(5,I1) + q_init(4,I1,k)
         qb(6,I1) = qb(6,I1) + q_init(5,I1,k)
      end do
      qb(2,I1) = qb(1,I1) - pbprime_init(I1)
   end do

   ! do I1 = 1, npoin_q
   !    !write(11,*) I1, q_init(1,I1,:)
   !    print*, I1, q_init(1,I1,:)
   !    !print*, I1, qb(2,I1)
   ! end do
   !close(11)
   
   

   do ilr = 1, 2 ! 1 = left, 2 = right
      do iface = 1, nface
         do k = 1, nlayers
            qb_face(1,ilr,:,iface) = qb_face(1,ilr,:,iface) + q_face(1,ilr,:,iface,k)
            qb_face(5,ilr,:,iface) = qb_face(5,ilr,:,iface) + q_face(4,ilr,:,iface,k)
            qb_face(6,ilr,:,iface) = qb_face(6,ilr,:,iface) + q_face(5,ilr,:,iface,k)
         end do
         qb_face(2,ilr,:,iface) = qb_face(1,ilr,:,iface) - pbprime_face(ilr,:,iface)
      end do
   end do
    
    ! Compute the derived quantities ubar, vbar, u', v', q_init'_u, q_init'_v

   qb(3,:) = qb(3,:)/qb(1,:)
   qb(4,:) = qb(4,:)/qb(1,:)

   

   do ilr = 1, 2
      do iface = 1, nface
         qb_face(3,ilr,:,iface) = qb_face(3,ilr,:,iface)/qb_face(1,ilr,:,iface)
         qb_face(4,ilr,:,iface) = qb_face(4,ilr,:,iface)/qb_face(1,ilr,:,iface)
      end do
   end do

   ! Compute values of u' and v' at quadrature points and face of each grid element

   do k = 1, nlayers
      qprime(2,:,k) = q_init(2,:,k) - qb(3,:) ! u' = u - ubar
      qprime(3,:,k) = q_init(3,:,k) - qb(4,:) ! v' = v - vbar
   end do

   do ilr = 1, 2
      do iface = 1, nface
         do k = 1, nlayers
            qprime_face(2,ilr,:,iface,k) = q_face(2,ilr,:,iface,k) - qb_face(3,ilr,:,iface)
            qprime_face(3,ilr,:,iface,k) = q_face(3,ilr,:,iface,k) - qb_face(4,ilr,:,iface)
         end do
      end do
   end do


end subroutine initial_conditions_mlswe



! subroutine initialize_barotropic(qprime_face,qb,qb_df,qb_face,qprime,q_init,q_df,pbprime_df,pbprime_init,pbprime_face,q_face, nlayers)

!    use mod_grid, only: npoin, npoin_q, nface
!    use mod_basis, only: nq

!    implicit none

!    real, dimension(6,npoin_q), intent(out) :: qb
!    real, dimension(6,2,nq,nface), intent(out) :: qb_face
!    real, dimension(3,2,nq,nface), intent(inout) :: qprime_face
!    real, dimension(3,npoin), intent(out) :: qb_df
!    real, dimension(2,npoin_q,nlayers), intent(inout) :: qprime

!    real, dimension(5,npoin_q,nlayers), intent(in) :: q_init
!    real, dimension(3,npoin,nlayers), intent(in) :: q_df
!    real, dimension(npoin), intent(in) :: pbprime_df
!    real, dimension(npoin_q), intent(in) :: pbprime_init
!    real, dimension(2,nq,nface), intent(in) :: pbprime_face
!    real, dimension(5,2,nq,nface,nlayers), intent(in) :: q_face

!    integer, intent(in) :: nlayers
   
!    integer :: ilr, k, i, j, iface, iquad, I1, Iq
    
!    qb = 0.0
!    qb_face = 0.0
!    qb_df = 0.0
   
!    ! Compute degrees of freedom for the barotropic mass and momentum
!    ! dependent variables. These are the vertical sums of the
!    ! degrees of freedom for the corresponding layer variables.
!    ! Also compute degrees of freedom for pbpert = pb - pbprime_init.
   
!    do I1 = 1, npoin
!       do k = 1, nlayers
!          qb_df(1,I1) = qb_df(1,I1) + q_df(1,I1,k)
!          qb_df(2,I1) = qb_df(2,I1) + q_df(2,I1,k)
!          qb_df(3,I1) = qb_df(3,I1) + q_df(3,I1,k)
!       end do
!       qb_df(2,I1) = qb_df(1,I1) - pbprime_df(I1)
!    end do
   
!    ! Compute pointwise values of the barotropic mass and momentum
!    ! dependent variables at quadrature points and endpoints of each
!    ! grid cell. Use vertical sums of the pointwise values of the
!    ! layer variables.
!    ! Also compute pointwise values of pbpert = pb - pbprime_init.
   
!    do I1 = 1, npoin_q
!       do k = 1, nlayers
!          qb(1,I1) = qb(1,I1) + q_init(1,I1,k)
!          qb(5,I1) = qb(5,I1) + q_init(4,I1,k)
!          qb(6,I1) = qb(6,I1) + q_init(5,I1,k)
!       end do
!       qb(2,I1) = qb(1,I1) - pbprime_init(I1)
!    end do
   
!    do ilr = 1, 2 ! 1 = left, 2 = right
!       do iface = 1, nface
!          do k = 1, nlayers
!             qb_face(1,ilr,:,iface) = qb_face(1,ilr,:,iface) + q_face(1,ilr,:,iface,k)
!             qb_face(5,ilr,:,iface) = qb_face(5,ilr,:,iface) + q_face(5,ilr,:,iface,k)
!             qb_face(6,ilr,:,iface) = qb_face(6,ilr,:,iface) + q_face(6,ilr,:,iface,k)
!          end do
!          qb_face(2,ilr,:,iface) = qb_face(1,ilr,:,iface) - pbprime_face(ilr,:,iface)
!       end do
!    end do
    
!     ! Compute the derived quantities ubar, vbar, u', v', q_init'_u, q_init'_v

!    qb(3,:) = qb(3,:)/qb(1,:)
!    qb(4,:) = qb(4,:)/qb(1,:)

!    do ilr = 1, 2
!       do iface = 1, nface
!          qb_face(3,ilr,:,iface) = qb_face(3,ilr,:,iface)/qb_face(1,ilr,:,iface)
!          qb_face(4,ilr,:,iface) = qb_face(4,ilr,:,iface)/qb_face(1,ilr,:,iface)
!       end do
!    end do

!    ! Compute values of u' and v' at quadrature points and face of each grid element

!    do k = 1, nlayers
!       qprime(2,:,k) = q_init(2,:,k) - qb(3,:) ! u' = u - ubar
!       qprime(3,:,k) = q_init(3,:,k) - qb(4,:) ! v' = v - vbar
!    end do

!    do ilr = 1, 2
!       do iface = 1, nface
!          do k = 1, nlayers
!             qprime_face(2,ilr,:,iface,k) = q_face(2,ilr,:,iface,k) - qb_face(3,ilr,:,iface)
!             qprime_face(3,ilr,:,iface,k) = q_face(3,ilr,:,iface,k) - qb_face(4,ilr,:,iface)
!          end do
!       end do
!    end do

! end subroutine initialize_barotropic

