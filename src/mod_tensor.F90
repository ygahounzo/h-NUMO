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

      ! Quadrature-point z-direction arrays (non-zero on sphere)
      real,    dimension(:,:), allocatable :: dpsidz     ! (npts, npoin_q): h_e*e_z + h_n*n_z
      real,    dimension(:,:), allocatable :: dpsidz_x   ! (npts, npoin_q): hi*zetaq_x (outward normal x)
      real,    dimension(:,:), allocatable :: dpsidz_y   ! (npts, npoin_q): hi*zetaq_y
      real,    dimension(:,:), allocatable :: dpsidz_z   ! (npts, npoin_q): hi*zetaq_z (= hi*sin(lat) on sphere)

      ! DOF grid (npoin)
      real,    dimension(:,:), allocatable :: psih_df    ! (npts, npoin)
      real,    dimension(:,:), allocatable :: dpsidx_df  ! (npts, npoin)
      real,    dimension(:,:), allocatable :: dpsidy_df  ! (npts, npoin)
      real,    dimension(:,:), allocatable :: dpsidz_df  ! (npts, npoin): h_e*e_z + h_n*n_z (DOF)
      real,    dimension(:,:), allocatable :: dpsidz_df_x ! (npts, npoin): hi*zetaq_x (DOF)
      real,    dimension(:,:), allocatable :: dpsidz_df_y ! (npts, npoin): hi*zetaq_y (DOF)
      real,    dimension(:,:), allocatable :: dpsidz_df_z ! (npts, npoin): hi*zetaq_z (DOF)
      integer, dimension(:,:), allocatable :: index_df   ! (npts, npoin)
      real,    dimension(:),   allocatable :: wjac_df    ! (npoin)
      integer, dimension(:,:), allocatable :: index_df_elt ! (npts, nelem)

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
      real    :: h_e, h_n, hi
      real    :: e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z

      ! npts = number of DOF basis functions per element (nglx*ngly)
      npts = b%nglx * b%ngly
      npts_quad = b%nqx * b%nqy

      ! Allocate quadrature-point arrays
      allocate(tsp%wjac(G%npoin_q))
      allocate(tsp%psih(npts, G%npoin_q))
      allocate(tsp%dpsidx(npts, G%npoin_q))
      allocate(tsp%dpsidy(npts, G%npoin_q))
      allocate(tsp%dpsidz(npts, G%npoin_q))
      allocate(tsp%dpsidz_x(npts, G%npoin_q))
      allocate(tsp%dpsidz_y(npts, G%npoin_q))
      allocate(tsp%dpsidz_z(npts, G%npoin_q))
      allocate(tsp%indexq(npts, G%npoin_q))
      allocate(tsp%indexq_e(npts_quad, G%nelem))

      ! Allocate DOF arrays
      allocate(tsp%wjac_df(G%npoin))
      allocate(tsp%psih_df(npts, G%npoin))
      allocate(tsp%dpsidx_df(npts, G%npoin))
      allocate(tsp%dpsidy_df(npts, G%npoin))
      allocate(tsp%dpsidz_df(npts, G%npoin))
      allocate(tsp%dpsidz_df_x(npts, G%npoin))
      allocate(tsp%dpsidz_df_y(npts, G%npoin))
      allocate(tsp%dpsidz_df_z(npts, G%npoin))
      allocate(tsp%index_df(npts, G%npoin))
      allocate(tsp%index_df_elt(npts, G%nelem))

      ! Initialise quadrature arrays
      tsp%wjac     = 0.0
      tsp%psih     = 0.0
      tsp%dpsidx   = 0.0
      tsp%dpsidy   = 0.0
      tsp%dpsidz   = 0.0
      tsp%dpsidz_x = 0.0
      tsp%dpsidz_y = 0.0
      tsp%dpsidz_z = 0.0
      tsp%indexq   = 0
      tsp%indexq_e = 0

      ! Initialise DOF arrays
      tsp%wjac_df     = 0.0
      tsp%psih_df     = 0.0
      tsp%dpsidx_df   = 0.0
      tsp%dpsidy_df   = 0.0
      tsp%dpsidz_df   = 0.0
      tsp%dpsidz_df_x = 0.0
      tsp%dpsidz_df_y = 0.0
      tsp%dpsidz_df_z = 0.0
      tsp%index_df    = 0
      tsp%index_df_elt = 0

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
               e_z = mt%ksiq_z(iquad, jquad, 1, e)
               n_x = mt%etaq_x(iquad, jquad, 1, e);  n_y = mt%etaq_y(iquad, jquad, 1, e)
               n_z = mt%etaq_z(iquad, jquad, 1, e)
               c_x = mt%zetaq_x(iquad, jquad, 1, e)
               c_y = mt%zetaq_y(iquad, jquad, 1, e)
               c_z = mt%zetaq_z(iquad, jquad, 1, e)

               ip = 0
               do m = 1, b%ngly
                  do n = 1, b%nglx

                     I  = G%intma(n, m, 1, e)
                     ip = ip + 1
                     hi = b%psiqx(n, iquad) * b%psiqy(m, jquad)

                     tsp%indexq(ip, Iq) = I
                     tsp%psih(ip, Iq)   = hi

                     h_e = b%dpsiqx(n, iquad) * b%psiqy(m, jquad)   ! Xi  derivative
                     h_n = b%psiqx(n, iquad)  * b%dpsiqy(m, jquad)  ! Eta derivative

                     tsp%dpsidx(ip, Iq)   = h_e * e_x + h_n * n_x
                     tsp%dpsidy(ip, Iq)   = h_e * e_y + h_n * n_y
                     tsp%dpsidz(ip, Iq)   = h_e * e_z + h_n * n_z
                     tsp%dpsidz_x(ip, Iq) = hi * c_x
                     tsp%dpsidz_y(ip, Iq) = hi * c_y
                     tsp%dpsidz_z(ip, Iq) = hi * c_z

                  end do ! n
               end do ! m
            enddo
         enddo

      end do ! e, jquad, iquad

      !  DOF tensor product  (npoin grid)
      do e = 1, G%nelem
         ipq = 0
         do jquad = 1, b%ngly
            do iquad = 1, b%nglx

               Iq = G%intma(iquad, jquad, 1, e)
               ipq = ipq + 1
               tsp%index_df_elt(ipq,e) = Iq

               tsp%wjac_df(Iq) = mt%jac(iquad, jquad, 1, e)

               e_x = mt%ksi_x(iquad, jquad, 1, e);  e_y = mt%ksi_y(iquad, jquad, 1, e)
               e_z = mt%ksi_z(iquad, jquad, 1, e)
               n_x = mt%eta_x(iquad, jquad, 1, e);  n_y = mt%eta_y(iquad, jquad, 1, e)
               n_z = mt%eta_z(iquad, jquad, 1, e)
               c_x = mt%zeta_x(iquad, jquad, 1, e)
               c_y = mt%zeta_y(iquad, jquad, 1, e)
               c_z = mt%zeta_z(iquad, jquad, 1, e)

               ip = 0
               do m = 1, b%ngly
                  do n = 1, b%nglx

                     I  = G%intma(n, m, 1, e)
                     ip = ip + 1
                     hi = b%psix(n, iquad) * b%psiy(m, jquad)

                     tsp%index_df(ip, Iq) = I
                     tsp%psih_df(ip, Iq)  = hi

                     h_e = b%dpsix(n, iquad) * b%psiy(m, jquad)   ! Xi  derivative
                     h_n = b%psix(n, iquad)  * b%dpsiy(m, jquad)  ! Eta derivative

                     tsp%dpsidx_df(ip, Iq)   = h_e * e_x + h_n * n_x
                     tsp%dpsidy_df(ip, Iq)   = h_e * e_y + h_n * n_y
                     tsp%dpsidz_df(ip, Iq)   = h_e * e_z + h_n * n_z
                     tsp%dpsidz_df_x(ip, Iq) = hi * c_x
                     tsp%dpsidz_df_y(ip, Iq) = hi * c_y
                     tsp%dpsidz_df_z(ip, Iq) = hi * c_z

                  end do ! n
               end do ! m
            enddo
         enddo

      end do ! e, jquad, iquad

   end subroutine compute_tensor_product

   subroutine compute_gradient_quad(G, inp, b, mt, grad_q_quad, q)

      use mod_grid,    only: grid
      use mod_input,   only: input
      use mod_basis,   only: basis
      use mod_metrics, only: metrics

      implicit none

      type(grid),    intent(in) :: G
      type(input),   intent(in) :: inp
      type(basis),   intent(in) :: b
      type(metrics), intent(in) :: mt

      real, dimension(G%npoin),             intent(in)  :: q
      real, dimension(inp%ngrd_var,G%npoin_q), intent(out) :: grad_q_quad

      integer :: e, iquad, jquad, m, n, Iq, I
      real :: e_x, e_y, e_z, n_x, n_y, n_z, h_e, h_n
      logical :: has_z

      grad_q_quad = 0.0
      has_z = (inp%ngrd_var == 3) ! 3rd (z) component only meaningful on the sphere

      do e = 1, G%nelem
         do jquad = 1, b%nqy
            do iquad = 1, b%nqx

               e_x = mt%ksiq_x(iquad,jquad,1,e); e_y = mt%ksiq_y(iquad,jquad,1,e); e_z = mt%ksiq_z(iquad,jquad,1,e)
               n_x = mt%etaq_x(iquad,jquad,1,e); n_y = mt%etaq_y(iquad,jquad,1,e); n_z = mt%etaq_z(iquad,jquad,1,e)

               Iq = G%intma_dg_quad(iquad,jquad,1,e)

               do m = 1, b%ngly
                  do n = 1, b%nglx

                     I = G%intma(n,m,1,e)

                     h_e = b%dpsiqx(n,iquad)*b%psiqy(m,jquad)   ! Xi  derivative
                     h_n = b%psiqx(n,iquad)*b%dpsiqy(m,jquad)   ! Eta derivative

                     grad_q_quad(1,Iq) = grad_q_quad(1,Iq) + (h_e*e_x + h_n*n_x)*q(I)
                     grad_q_quad(2,Iq) = grad_q_quad(2,Iq) + (h_e*e_y + h_n*n_y)*q(I)
                     if (has_z) grad_q_quad(3,Iq) = grad_q_quad(3,Iq) + (h_e*e_z + h_n*n_z)*q(I)

                  end do
               end do

            end do !iquad
         end do !jquad
      end do !e

   end subroutine compute_gradient_quad

   subroutine compute_gradient_df(G, inp, b, mt, grad_q_df, q)

      use mod_grid,    only: grid
      use mod_input,   only: input
      use mod_basis,   only: basis
      use mod_metrics, only: metrics

      implicit none

      type(grid),    intent(in) :: G
      type(input),   intent(in) :: inp
      type(basis),   intent(in) :: b
      type(metrics), intent(in) :: mt

      real, dimension(G%npoin),           intent(in)  :: q
      real, dimension(inp%ngrd_var,G%npoin), intent(out) :: grad_q_df

      integer :: e, iquad, jquad, m, n, Iq, I
      real :: e_x, e_y, e_z, n_x, n_y, n_z, h_e, h_n
      logical :: has_z

      grad_q_df = 0.0
      has_z = (inp%ngrd_var == 3) ! 3rd (z) component only meaningful on the sphere

      !Construct Volume Integral Contribution
      do e = 1, G%nelem

         !Loop over DOF points
         do jquad = 1, b%ngly
            do iquad = 1, b%nglx

               e_x = mt%ksi_x(iquad,jquad,1,e); e_y = mt%ksi_y(iquad,jquad,1,e); e_z = mt%ksi_z(iquad,jquad,1,e)
               n_x = mt%eta_x(iquad,jquad,1,e); n_y = mt%eta_y(iquad,jquad,1,e); n_z = mt%eta_z(iquad,jquad,1,e)

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
                     if (has_z) grad_q_df(3,Iq) = grad_q_df(3,Iq) + (h_e*e_z + h_n*n_z)*q(I)

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
