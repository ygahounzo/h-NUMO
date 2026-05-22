module mod_tensor

   ! Provides the tensor_CS derived type and the Tensor_product type-bound
   ! procedure, which pre-computes basis functions and metric terms on both
   ! the quadrature-point grid (npoin_q) and the DOF grid (npoin).

   implicit none

   !  Tensor-product control structure
   type :: tensor_CS

      ! Quadrature-point grid (npoin_q)
      real,    dimension(:,:), allocatable :: psih       ! (npts, npoin_q)
      real,    dimension(:,:), allocatable :: dpsidx     ! (npts, npoin_q)
      real,    dimension(:,:), allocatable :: dpsidy     ! (npts, npoin_q)
      integer, dimension(:,:), allocatable :: indexq     ! (npts, npoin_q)
      integer, dimension(:,:), allocatable :: indexq_e     ! (npts_quad, nelem)
      real,    dimension(:),   allocatable :: wjac       ! (npoin_q)

      ! DOF grid (npoin)
      real,    dimension(:,:), allocatable :: psih_df    ! (npts, npoin)
      real,    dimension(:,:), allocatable :: dpsidx_df  ! (npts, npoin)
      real,    dimension(:,:), allocatable :: dpsidy_df  ! (npts, npoin)
      integer, dimension(:,:), allocatable :: index_df   ! (npts, npoin)
      real,    dimension(:),   allocatable :: wjac_df    ! (npoin)

   end type tensor_CS

   public :: compute_tensor_product, tensor_CS, &
      compute_gradient_quad, compute_gradient_df, interpolate_layer_from_quad_to_node_1d

   private

contains

   subroutine compute_tensor_product(G, b, mt, tsp)

      use mod_grid,    only: grid
      use mod_basis,   only: basis
      use mod_metrics, only: metrics

      implicit none

      type(grid),      intent(in)  :: G
      type(basis),     intent(in)  :: b
      type(metrics),   intent(in)  :: mt
      type(tensor_CS), intent(out) :: tsp

      integer :: e, jquad, iquad, Iq, ip, m, n, I, ipq
      integer :: npts, npts_quad
      real    :: h_e, h_n
      real    :: e_x, e_y, n_x, n_y

      ! npts = number of DOF basis functions per element (nglx*ngly)
      npts = b%nglx * b%ngly
      npts_quad = b%nqx * b%nqy

      ! Allocate quadrature-point arrays
      allocate(tsp%wjac(G%npoin_q))
      allocate(tsp%psih(npts, G%npoin_q))
      allocate(tsp%dpsidx(npts, G%npoin_q))
      allocate(tsp%dpsidy(npts, G%npoin_q))
      allocate(tsp%indexq(npts, G%npoin_q))
      allocate(tsp%indexq_e(npts_quad, G%nelem))

      ! Allocate DOF arrays
      allocate(tsp%wjac_df(G%npoin))
      allocate(tsp%psih_df(npts, G%npoin))
      allocate(tsp%dpsidx_df(npts, G%npoin))
      allocate(tsp%dpsidy_df(npts, G%npoin))
      allocate(tsp%index_df(npts, G%npoin))

      ! Initialise quadrature arrays
      tsp%wjac   = 0.0
      tsp%psih   = 0.0
      tsp%dpsidx = 0.0
      tsp%dpsidy = 0.0
      tsp%indexq = 0
      tsp%indexq_e = 0

      ! Initialise DOF arrays
      tsp%wjac_df   = 0.0
      tsp%psih_df   = 0.0
      tsp%dpsidx_df = 0.0
      tsp%dpsidy_df = 0.0
      tsp%index_df  = 0

      !  Quadrature-point tensor product  (npoin_q grid)
      do e = 1, G%nelem

         ipq = 0
         do jquad = 1, b%nqy
            do iquad = 1, b%nqx

               Iq = G%intma_dg_quad(iquad, jquad, 1, e)

               ipq = ipq + 1
               tsp%indexq_e(ipq,e) = Iq

               tsp%wjac(Iq) = mt%jacq(iquad, jquad, 1, e)

               e_x = mt%ksiq_x(iquad, jquad, 1, e);  e_y = mt%ksiq_y(iquad, jquad, 1, e)
               n_x = mt%etaq_x(iquad, jquad, 1, e);  n_y = mt%etaq_y(iquad, jquad, 1, e)

               ip = 0
               do m = 1, b%ngly
                  do n = 1, b%nglx

                     I  = G%intma(n, m, 1, e)
                     ip = ip + 1

                     tsp%indexq(ip, Iq) = I
                     tsp%psih(ip, Iq)   = b%psiqx(n, iquad) * b%psiqy(m, jquad)

                     h_e = b%dpsiqx(n, iquad) * b%psiqy(m, jquad)   ! Xi  derivative
                     h_n = b%psiqx(n, iquad)  * b%dpsiqy(m, jquad)  ! Eta derivative

                     tsp%dpsidx(ip, Iq) = h_e * e_x + h_n * n_x
                     tsp%dpsidy(ip, Iq) = h_e * e_y + h_n * n_y

                  end do ! n
               end do ! m
            enddo
         enddo

      end do ! e, jquad, iquad

      !  DOF tensor product  (npoin grid)
      do concurrent(e = 1:G%nelem, jquad = 1:b%ngly, iquad = 1:b%nglx)

         Iq = G%intma(iquad, jquad, 1, e)

         tsp%wjac_df(Iq) = mt%jac(iquad, jquad, 1, e)

         e_x = mt%ksi_x(iquad, jquad, 1, e);  e_y = mt%ksi_y(iquad, jquad, 1, e)
         n_x = mt%eta_x(iquad, jquad, 1, e);  n_y = mt%eta_y(iquad, jquad, 1, e)

         ip = 0
         do m = 1, b%ngly
            do n = 1, b%nglx

               I  = G%intma(n, m, 1, e)
               ip = ip + 1

               tsp%index_df(ip, Iq) = I
               tsp%psih_df(ip, Iq)  = b%psix(n, iquad) * b%psiy(m, jquad)

               h_e = b%dpsix(n, iquad) * b%psiy(m, jquad)   ! Xi  derivative
               h_n = b%psix(n, iquad)  * b%dpsiy(m, jquad)  ! Eta derivative

               tsp%dpsidx_df(ip, Iq) = h_e * e_x + h_n * n_x
               tsp%dpsidy_df(ip, Iq) = h_e * e_y + h_n * n_y

            end do ! n
         end do ! m

      end do ! e, jquad, iquad

   end subroutine compute_tensor_product

   subroutine compute_gradient_quad(G, b, mt, grad_q_quad, q)

      use mod_grid,    only: grid
      use mod_basis,   only: basis
      use mod_metrics, only: metrics

      implicit none

      type(grid),    intent(in) :: G
      type(basis),   intent(in) :: b
      type(metrics), intent(in) :: mt

      real, dimension(G%npoin),   intent(in)  :: q
      real, dimension(2,G%npoin_q), intent(out) :: grad_q_quad

      integer :: e, iquad, jquad, kquad, l, m, n, Iq, I
      real :: e_x, e_y, n_x, n_y, h_e, h_n, h_c, c_x, c_y, c_z, e_z, n_z

      grad_q_quad = 0.0

      !Construct Volume Integral Contribution
      do e = 1, G%nelem

         !Loop Integration Points
         !do kquad = 1, nqz
         do jquad = 1, b%nqy
            do iquad = 1, b%nqx


               e_x = mt%ksiq_x(iquad,jquad,1,e); e_y = mt%ksiq_y(iquad,jquad,1,e)
               n_x = mt%etaq_x(iquad,jquad,1,e); n_y = mt%etaq_y(iquad,jquad,1,e)

               Iq = G%intma_dg_quad(iquad,jquad,1,e)

               !Interpolate at Integration Points
               !do l = 1, nglz
               do m = 1, b%ngly
                  do n = 1, b%nglx

                     I = G%intma(n,m,1,e)

                     !Xi derivatives
                     h_e = b%dpsiqx(n,iquad)*b%psiqy(m,jquad)
                     !Eta derivatives
                     h_n = b%psiqx(n,iquad)*b%dpsiqy(m,jquad)

                     grad_q_quad(1,Iq) = grad_q_quad(1,Iq) + (h_e*e_x + h_n*n_x)*q(I)
                     grad_q_quad(2,Iq) = grad_q_quad(2,Iq) + (h_e*e_y + h_n*n_y)*q(I)

                  end do
               end do
               !end do

            end do !iquad
         end do !jquad
         !end do !kquad
      end do !e

   end subroutine compute_gradient_quad

   subroutine compute_gradient_df(G, b, mt, grad_q_df, q)

      use mod_grid,    only: grid
      use mod_basis,   only: basis
      use mod_metrics, only: metrics

      implicit none

      type(grid),    intent(in) :: G
      type(basis),   intent(in) :: b
      type(metrics), intent(in) :: mt

      real, dimension(G%npoin),   intent(in)  :: q
      real, dimension(2,G%npoin), intent(out) :: grad_q_df

      integer :: e, iquad, jquad, m, n, Iq, I
      real :: e_x, e_y, n_x, n_y, h_e, h_n

      grad_q_df = 0.0

      !Construct Volume Integral Contribution
      do e = 1, G%nelem

         !Loop over DOF points
         do jquad = 1, b%ngly
            do iquad = 1, b%nglx

               e_x = mt%ksi_x(iquad,jquad,1,e); e_y = mt%ksi_y(iquad,jquad,1,e)
               n_x = mt%eta_x(iquad,jquad,1,e); n_y = mt%eta_y(iquad,jquad,1,e)

               Iq = G%intma(iquad,jquad,1,e)

               do m = 1, b%ngly
                  do n = 1, b%nglx

                     I = G%intma(n,m,1,e)

                     !Xi derivatives
                     h_e = b%dpsix(n,iquad)*b%psiy(m,jquad)
                     !Eta derivatives
                     h_n = b%psix(n,iquad)*b%dpsiy(m,jquad)

                     grad_q_df(1,Iq) = grad_q_df(1,Iq) + (h_e*e_x + h_n*n_x)*q(I)
                     grad_q_df(2,Iq) = grad_q_df(2,Iq) + (h_e*e_y + h_n*n_y)*q(I)

                  end do
               end do

            end do !iquad
         end do !jquad
      end do !e

   end subroutine compute_gradient_df

   subroutine interpolate_layer_from_quad_to_node_1d(G, b, mt, q_df, q)

      use mod_grid,    only: grid
      use mod_basis,   only: basis
      use mod_metrics, only: metrics

      implicit none

      type(grid),    intent(in) :: G
      type(basis),   intent(in) :: b
      type(metrics), intent(in) :: mt

      real, dimension(G%npoin),   intent(out) :: q_df
      real, dimension(G%npoin_q), intent(in)  :: q

      real :: wq, var_u, var_v, hi
      integer :: iquad, jquad, kquad, e, k, m, n, l, I, Iq, ip

      q_df = 0.0

      do e = 1, G%nelem
         do kquad = 1, b%nqz
            do jquad = 1, b%nqy
               do iquad = 1, b%nqx

                  Iq = G%intma_dg_quad(iquad, jquad, kquad, e)

                  wq = mt%jacq(iquad,jquad,kquad,e)
                  var_u = q(Iq)

                  do l = 1, b%nglz
                     do m = 1, b%ngly
                        do n = 1, b%nglx

                           I = G%intma(n, m, l, e)

                           hi = b%psiqx(n, iquad) * b%psiqy(m, jquad)
                           q_df(I) = q_df(I) + wq*var_u*hi

                        end do
                     end do
                  end do

               end do
            end do
         end do
      end do

      q_df(:) = mt%massinv(:)*q_df(:)

   end subroutine interpolate_layer_from_quad_to_node_1d

end module mod_tensor
