!----------------------------------------------------------------------!
!>@brief This subroutine builds the METRIC TERMS see Giraldo 2001 (IJNMF)
!>@author  Francis X. Giraldo on 7/08
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@date 11 September 2009  James F. Kelly
!> Generalized to 3D
!> Combined to do both Cube and Sphere domains.
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: April 03, 2023
!       remove SEMI_ANALYTIC_METRICS and is_sphere flag
!----------------------------------------------------------------------!
subroutine metrics(ksi_x,ksi_y,ksi_z,eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,jac,xjac)

    use mod_basis, only: nglx, ngly, nglz, wglx, wgly, wglz, nqx,nqy,nqz,wnqx,wnqy,wnqz
 
    use mod_grid, only: intma, coord, nelem, intma_dg_quad

    use mod_interface, only: compute_local_gradient_v3, compute_local_gradient_quad_v3

    implicit none

    !global arrays
    real, dimension(nglx,ngly,nglz,nelem), intent(out) :: ksi_x,  ksi_y,  ksi_z
    real, dimension(nglx,ngly,nglz,nelem), intent(out) :: eta_x,  eta_y,  eta_z
    real, dimension(nglx,ngly,nglz,nelem), intent(out) :: zeta_x, zeta_y, zeta_z, jac, xjac

    !local arrays
    real, dimension(nglx,ngly,nglz) :: x, y, z
    real, dimension(nglx,ngly,nglz) :: x_ksi, x_eta, x_zeta
    real, dimension(nglx,ngly,nglz) :: y_ksi, y_eta, y_zeta
    real, dimension(nglx,ngly,nglz) :: z_ksi, z_eta, z_zeta

    real xj
    integer ie, i, j, k
    integer ip, ndim

    !Define Dimension of Matrix going into MXM routine
    ndim=1

    !initialize the global matrix
    ksi_x=0
    ksi_y=0
    ksi_z=0
    eta_x=0
    eta_y=0
    eta_z=0
    zeta_x=0
    zeta_y=0
    zeta_z=0
    jac=0
  
    !loop thru the elements
    do ie=1,nelem

        !Store Element Variables
        do k=1,nglz
            do j=1,ngly
                do i=1,nglx
                    ip=intma(i,j,k,ie)
                    x(i,j,k)=coord(1,ip)
                    y(i,j,k)=coord(2,ip)
                    z(i,j,k)=coord(3,ip)
                end do
            end do
        end do

        !Construct Mapping Derivatives: dx/dksi, dx/deta,dx/dzeta, dy/dksi, dy/deta, dy/dzeta, dz/dksi, dz/deta, dz/dzeta

        call compute_local_gradient_v3(x_ksi,x_eta,x_zeta,x,nglx,ngly,nglz,ndim)
        call compute_local_gradient_v3(y_ksi,y_eta,y_zeta,y,nglx,ngly,nglz,ndim)
        call compute_local_gradient_v3(z_ksi,z_eta,z_zeta,z,nglx,ngly,nglz,ndim)

        !Construct Inverse Mapping
        do k=1,nglz
            do j=1,ngly
                do i=1,nglx
               
                    !set jacobian matrix J for 2D

                    if(nglx == 1) then
                        x_ksi(i,j,k) = 1;  y_ksi(i,j,k) = 0;  z_ksi(i,j,k) = 0;
                        x_eta(i,j,k) = 0;  x_zeta(i,j,k) = 0;
                    endif
                    if(ngly == 1) then
                        x_eta(i,j,k) = 0;  y_eta(i,j,k) = 1;  z_eta(i,j,k) = 0;
                        y_ksi(i,j,k) = 0; y_zeta(i,j,k) = 0;
                    endif
                    if(nglz == 1) then
                        x_zeta(i,j,k) = 0; y_zeta(i,j,k) = 0; z_zeta(i,j,k) = 1;
                        z_ksi(i,j,k) = 0; z_eta(i,j,k) = 0;
                    endif
              
                    !compute inverse of J

                    xj = &
                        (x_ksi(i,j,k)*y_eta(i,j,k)*z_zeta(i,j,k) - x_ksi(i,j,k)*y_zeta(i,j,k)*z_eta(i,j,k))  &
                    - (y_ksi(i,j,k)*x_eta(i,j,k)*z_zeta(i,j,k) - y_ksi(i,j,k)*x_zeta(i,j,k)*z_eta(i,j,k))  &
                    + (z_ksi(i,j,k)*x_eta(i,j,k)*y_zeta(i,j,k) - z_ksi(i,j,k)*x_zeta(i,j,k)*y_eta(i,j,k))

                    ksi_x(i,j,k,ie)=  (y_eta(i,j,k)*z_zeta(i,j,k)-y_zeta(i,j,k)*z_eta(i,j,k))/xj
                    ksi_y(i,j,k,ie)= -(x_eta(i,j,k)*z_zeta(i,j,k)-x_zeta(i,j,k)*z_eta(i,j,k))/xj
                    ksi_z(i,j,k,ie)=  (x_eta(i,j,k)*y_zeta(i,j,k)-x_zeta(i,j,k)*y_eta(i,j,k))/xj
                    eta_x(i,j,k,ie)= -(y_ksi(i,j,k)*z_zeta(i,j,k)-y_zeta(i,j,k)*z_ksi(i,j,k))/xj
                    eta_y(i,j,k,ie)=  (x_ksi(i,j,k)*z_zeta(i,j,k)-x_zeta(i,j,k)*z_ksi(i,j,k))/xj
                    eta_z(i,j,k,ie)= -(x_ksi(i,j,k)*y_zeta(i,j,k)-x_zeta(i,j,k)*y_ksi(i,j,k))/xj
                    zeta_x(i,j,k,ie)= (y_ksi(i,j,k)*z_eta(i,j,k) -y_eta(i,j,k)*z_ksi(i,j,k) )/xj
                    zeta_y(i,j,k,ie)=-(x_ksi(i,j,k)*z_eta(i,j,k) -x_eta(i,j,k)*z_ksi(i,j,k) )/xj
                    zeta_z(i,j,k,ie)= (x_ksi(i,j,k)*y_eta(i,j,k) -x_eta(i,j,k)*y_ksi(i,j,k) )/xj

                    jac(i,j,k,ie)=wglx(i)*wgly(j)*wglz(k)*abs(xj)
                    xjac(i,j,k,ie)=xj
              
                    !fix inverse jacobian matrix for 2D

                    if(nglx == 1)  ksi_x(i,j,k,ie) = 0
                    if(ngly == 1)  eta_y(i,j,k,ie) = 0
                    if(nglz == 1)  zeta_z(i,j,k,ie) = 0
              
                end do
            end do
        end do

    end do !ie
end subroutine metrics
