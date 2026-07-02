!----------------------------------------------------------------------!
!>@brief This subroutine builds the METRIC TERMS see Giraldo 2001 (IJNMF)
!>@author  Francis X. Giraldo on 7/08
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: April 03, 2023
!----------------------------------------------------------------------!
subroutine compute_metrics(ksi_x,ksi_y,ksi_z,eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,jac,xjac, b, G)

    use mod_basis
    use mod_grid
    use mod_gradient, only: compute_local_gradient_v3

    implicit none

    type(basis), intent(in) :: b
    type(grid), intent(in) :: G

    !global arrays
    real, dimension(b%nglx,b%ngly,b%nglz,G%nelem), intent(out) :: ksi_x,  ksi_y,  ksi_z
    real, dimension(b%nglx,b%ngly,b%nglz,G%nelem), intent(out) :: eta_x,  eta_y,  eta_z
    real, dimension(b%nglx,b%ngly,b%nglz,G%nelem), intent(out) :: zeta_x, zeta_y, zeta_z, jac, xjac

    !local arrays
    real, dimension(b%nglx,b%ngly,b%nglz) :: x, y, z
    real, dimension(b%nglx,b%ngly,b%nglz) :: x_ksi, x_eta, x_zeta
    real, dimension(b%nglx,b%ngly,b%nglz) :: y_ksi, y_eta, y_zeta
    real, dimension(b%nglx,b%ngly,b%nglz) :: z_ksi, z_eta, z_zeta

    real xj
    real :: xn, yn, zn, rnorm
    integer ie, i, j, k
    integer ip, ndim

    !Define Dimension of Matrix going into MXM routine
    ndim=1

    !initialize the global matrix
    ksi_x = 0.0
    ksi_y = 0.0
    ksi_z = 0.0
    eta_x = 0.0
    eta_y = 0.0
    eta_z = 0.0
    zeta_x = 0.0
    zeta_y = 0.0
    zeta_z = 0.0
    jac = 0.0
    xjac = 0.0
  
    !loop thru the elements
    do ie = 1, G%nelem

        !Store Element Variables
        do k = 1, b%nglz
            do j = 1, b%ngly
                do i = 1, b%nglx
                    ip = G%intma(i,j,k,ie)
                    x(i,j,k) = G%coord(1,ip)
                    y(i,j,k) = G%coord(2,ip)
                    z(i,j,k) = G%coord(3,ip)
                end do
            end do
        end do

        ! Construct Mapping Derivatives: dx/dksi, dx/deta,dx/dzeta, dy/dksi, dy/deta, 
        ! dy/dzeta, dz/dksi, dz/deta, dz/dzeta

        call compute_local_gradient_v3(x_ksi,x_eta,x_zeta,x,b,ndim)
        call compute_local_gradient_v3(y_ksi,y_eta,y_zeta,y,b,ndim)
        call compute_local_gradient_v3(z_ksi,z_eta,z_zeta,z,b,ndim)

        !Construct Inverse Mapping
        do k = 1, b%nglz
            do j = 1, b%ngly
                do i = 1, b%nglx
               
                    !set jacobian matrix J for 2D

                    if(b%nglx == 1) then
                        x_ksi(i,j,k) = 1;  y_ksi(i,j,k) = 0;  z_ksi(i,j,k) = 0;
                        x_eta(i,j,k) = 0;  x_zeta(i,j,k) = 0;
                    endif
                    if(b%ngly == 1) then
                        x_eta(i,j,k) = 0;  y_eta(i,j,k) = 1;  z_eta(i,j,k) = 0;
                        y_ksi(i,j,k) = 0; y_zeta(i,j,k) = 0;
                    endif
                    if(b%nglz == 1) then
                        xn = y_ksi(i,j,k)*z_eta(i,j,k) - z_ksi(i,j,k)*y_eta(i,j,k)
                        yn = z_ksi(i,j,k)*x_eta(i,j,k) - x_ksi(i,j,k)*z_eta(i,j,k)
                        zn = x_ksi(i,j,k)*y_eta(i,j,k) - y_ksi(i,j,k)*x_eta(i,j,k)
                        rnorm = sqrt(xn*xn + yn*yn + zn*zn)
                        if (rnorm > 0.0) then
                            x_zeta(i,j,k) = xn / rnorm
                            y_zeta(i,j,k) = yn / rnorm
                            z_zeta(i,j,k) = zn / rnorm
                        else
                            x_zeta(i,j,k) = 0.0; y_zeta(i,j,k) = 0.0; z_zeta(i,j,k) = 1.0
                        end if
                    endif
              
                    !compute inverse of J

                    xj = (x_ksi(i,j,k)*y_eta(i,j,k)*z_zeta(i,j,k) &
                            - x_ksi(i,j,k)*y_zeta(i,j,k)*z_eta(i,j,k)) &
                            - (y_ksi(i,j,k)*x_eta(i,j,k)*z_zeta(i,j,k) &
                            - y_ksi(i,j,k)*x_zeta(i,j,k)*z_eta(i,j,k)) &
                            + (z_ksi(i,j,k)*x_eta(i,j,k)*y_zeta(i,j,k) &
                            - z_ksi(i,j,k)*x_zeta(i,j,k)*y_eta(i,j,k))

                    ksi_x(i,j,k,ie)=  (y_eta(i,j,k)*z_zeta(i,j,k)-y_zeta(i,j,k)*z_eta(i,j,k))/xj
                    ksi_y(i,j,k,ie)= -(x_eta(i,j,k)*z_zeta(i,j,k)-x_zeta(i,j,k)*z_eta(i,j,k))/xj
                    ksi_z(i,j,k,ie)=  (x_eta(i,j,k)*y_zeta(i,j,k)-x_zeta(i,j,k)*y_eta(i,j,k))/xj
                    eta_x(i,j,k,ie)= -(y_ksi(i,j,k)*z_zeta(i,j,k)-y_zeta(i,j,k)*z_ksi(i,j,k))/xj
                    eta_y(i,j,k,ie)=  (x_ksi(i,j,k)*z_zeta(i,j,k)-x_zeta(i,j,k)*z_ksi(i,j,k))/xj
                    eta_z(i,j,k,ie)= -(x_ksi(i,j,k)*y_zeta(i,j,k)-x_zeta(i,j,k)*y_ksi(i,j,k))/xj
                    zeta_x(i,j,k,ie)= (y_ksi(i,j,k)*z_eta(i,j,k) -y_eta(i,j,k)*z_ksi(i,j,k) )/xj
                    zeta_y(i,j,k,ie)=-(x_ksi(i,j,k)*z_eta(i,j,k) -x_eta(i,j,k)*z_ksi(i,j,k) )/xj
                    zeta_z(i,j,k,ie)= (x_ksi(i,j,k)*y_eta(i,j,k) -x_eta(i,j,k)*y_ksi(i,j,k) )/xj

                    jac(i,j,k,ie) = b%wglx(i)*b%wgly(j)*b%wglz(k)*abs(xj)
                    xjac(i,j,k,ie)=xj
              
                    !fix inverse jacobian matrix for 2D

                    if(b%nglx == 1)  ksi_x(i,j,k,ie) = 0
                    if(b%ngly == 1)  eta_y(i,j,k,ie) = 0
                    if(b%nglz == 1 .and. .not. G%is_sphere)  zeta_z(i,j,k,ie) = 0
              
                end do
            end do
        end do
    end do !ie
end subroutine compute_metrics
