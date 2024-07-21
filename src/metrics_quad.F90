!----------------------------------------------------------------------!
!>@brief This subroutine builds the METRIC TERMS
!>@ author by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: July 02, 2023
!----------------------------------------------------------------------!
subroutine metrics_quad(ksiq_x,ksiq_y,ksiq_z,etaq_x,etaq_y,etaq_z,zetaq_x,zetaq_y,zetaq_z,jacq,xjacq)

    use mod_basis, only: nglx, ngly, nglz,nqx,nqy,nqz,wnqx,wnqy,wnqz
 
    use mod_grid, only: intma, coord, nelem, intma_dg_quad

    use mod_gradient, only: compute_local_gradient_quad_v3

    implicit none

    !global arrays
    real, dimension(nqx,nqy,nqz,nelem), intent(out) :: ksiq_x,  ksiq_y,  ksiq_z
    real, dimension(nqx,nqy,nqz,nelem), intent(out) :: etaq_x,  etaq_y,  etaq_z
    real, dimension(nqx,nqy,nqz,nelem), intent(out) :: zetaq_x, zetaq_y, zetaq_z, jacq, xjacq

    !local arrays
    real, dimension(nglx,ngly,nglz) :: x, y, z
    real, dimension(nqx,nqy,nqz) :: x_ksiq, x_etaq, x_zetaq
    real, dimension(nqx,nqy,nqz) :: y_ksiq, y_etaq, y_zetaq
    real, dimension(nqx,nqy,nqz) :: z_ksiq, z_etaq, z_zetaq

    real xj
    integer ie, i, j, k
    integer ip, ndim

    !Define Dimension of Matrix going into MXM routine

    !initialize the global matrix
    ksiq_x=0.0
    ksiq_y=0.0
    ksiq_z=0.0
    etaq_x=0.0
    etaq_y=0.0
    etaq_z=0.0
    zetaq_x=0.0
    zetaq_y=0.0
    zetaq_z=0.0
    jacq=0.0
  
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
        call compute_local_gradient_quad_v3(x_ksiq,x_etaq,x_zetaq,x,nglx,ngly,nglz,nqx,nqy,nqz)
        call compute_local_gradient_quad_v3(y_ksiq,y_etaq,y_zetaq,y,nglx,ngly,nglz,nqx,nqy,nqz)
        call compute_local_gradient_quad_v3(z_ksiq,z_etaq,z_zetaq,z,nglx,ngly,nglz,nqx,nqy,nqz)

        !Construct Inverse Mapping
        do k=1,nqz
            do j=1,nqy
                do i=1,nqx
               
                    !set jacobian matrix J for 2D
                    if(nqx == 1) then
                        x_ksiq(i,j,k) = 1.0;  y_ksiq(i,j,k) = 0.0;  z_ksiq(i,j,k) = 0.0;
                        x_etaq(i,j,k) = 0.0;  x_zetaq(i,j,k) = 0.0;
                    endif
                    if(nqy == 1) then
                        x_etaq(i,j,k) = 0.0;  y_etaq(i,j,k) = 1.0;  z_etaq(i,j,k) = 0.0;
                        y_ksiq(i,j,k) = 0.0; y_zetaq(i,j,k) = 0.0;
                    endif
                    if(nqz == 1) then
                        x_zetaq(i,j,k) = 0.0; y_zetaq(i,j,k) = 0.0; z_zetaq(i,j,k) = 1.0;
                        z_ksiq(i,j,k) = 0.0; z_etaq(i,j,k) = 0.0;
                    endif
              
                    !compute inverse of J
                    xj = &
                        (x_ksiq(i,j,k)*y_etaq(i,j,k)*z_zetaq(i,j,k) - x_ksiq(i,j,k)*y_zetaq(i,j,k)*z_etaq(i,j,k))  &
                        - (y_ksiq(i,j,k)*x_etaq(i,j,k)*z_zetaq(i,j,k) - y_ksiq(i,j,k)*x_zetaq(i,j,k)*z_etaq(i,j,k))  &
                        + (z_ksiq(i,j,k)*x_etaq(i,j,k)*y_zetaq(i,j,k) - z_ksiq(i,j,k)*x_zetaq(i,j,k)*y_etaq(i,j,k))

                    ksiq_x(i,j,k,ie)=  (y_etaq(i,j,k)*z_zetaq(i,j,k)-y_zetaq(i,j,k)*z_etaq(i,j,k))/xj
                    ksiq_y(i,j,k,ie)= -(x_etaq(i,j,k)*z_zetaq(i,j,k)-x_zetaq(i,j,k)*z_etaq(i,j,k))/xj
                    ksiq_z(i,j,k,ie)=  (x_etaq(i,j,k)*y_zetaq(i,j,k)-x_zetaq(i,j,k)*y_etaq(i,j,k))/xj
                    etaq_x(i,j,k,ie)= -(y_ksiq(i,j,k)*z_zetaq(i,j,k)-y_zetaq(i,j,k)*z_ksiq(i,j,k))/xj
                    etaq_y(i,j,k,ie)=  (x_ksiq(i,j,k)*z_zetaq(i,j,k)-x_zetaq(i,j,k)*z_ksiq(i,j,k))/xj
                    etaq_z(i,j,k,ie)= -(x_ksiq(i,j,k)*y_zetaq(i,j,k)-x_zetaq(i,j,k)*y_ksiq(i,j,k))/xj
                    zetaq_x(i,j,k,ie)= (y_ksiq(i,j,k)*z_etaq(i,j,k) -y_etaq(i,j,k)*z_ksiq(i,j,k) )/xj
                    zetaq_y(i,j,k,ie)= -(x_ksiq(i,j,k)*z_etaq(i,j,k) -x_etaq(i,j,k)*z_ksiq(i,j,k) )/xj
                    zetaq_z(i,j,k,ie)= (x_ksiq(i,j,k)*y_etaq(i,j,k) -x_etaq(i,j,k)*y_ksiq(i,j,k) )/xj

                    jacq(i,j,k,ie)=wnqx(i)*wnqy(j)*wnqz(k)*abs(xj)
                    xjacq(i,j,k,ie)=xj
              
                    !fix inverse jacobian matrix for 2Dn
                    if(nqx == 1)  ksiq_x(i,j,k,ie) = 0.0
                    if(nqy == 1)  etaq_y(i,j,k,ie) = 0.0
                    if(nqz == 1)  zetaq_z(i,j,k,ie) = 0.0
              
                end do
            end do
        end do

    end do !ie
end subroutine metrics_quad
