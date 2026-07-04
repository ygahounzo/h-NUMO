!----------------------------------------------------------------------!
!>@brief This subroutine builds the METRIC TERMS
!>@ author by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: July 02, 2023
!----------------------------------------------------------------------!
subroutine metrics_quad(ksiq_x,ksiq_y,ksiq_z,etaq_x,etaq_y,etaq_z,zetaq_x,zetaq_y, &
            zetaq_z,jacq,xjacq, b, G)

    use mod_basis
    use mod_grid
    use mod_gradient, only: compute_local_gradient_quad_v3
    use mod_constants, only: earth_radius

    implicit none

    type(basis), intent(in) :: b
    type(grid), intent(in) :: G

    !global arrays
    real, dimension(b%nqx,b%nqy,b%nqz,G%nelem), intent(out) :: ksiq_x,  ksiq_y,  ksiq_z
    real, dimension(b%nqx,b%nqy,b%nqz,G%nelem), intent(out) :: etaq_x,  etaq_y,  etaq_z
    real, dimension(b%nqx,b%nqy,b%nqz,G%nelem), intent(out) :: zetaq_x, zetaq_y, zetaq_z, jacq, xjacq

    !local arrays
    real, dimension(b%nglx,b%ngly,b%nglz) :: x, y, z
    real, dimension(b%nqx,b%nqy,b%nqz) :: x_ksiq, x_etaq, x_zetaq
    real, dimension(b%nqx,b%nqy,b%nqz) :: y_ksiq, y_etaq, y_zetaq
    real, dimension(b%nqx,b%nqy,b%nqz) :: z_ksiq, z_etaq, z_zetaq

    real xj
    real xn, yn, zn, rnorm   ! surface normal (cross product) for nqz==1 case
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
    xjacq=0.0
  
    !loop thru the elements
    do ie=1,G%nelem

        !Store Element Variables
        do k=1,b%nglz
            do j=1,b%ngly
                do i=1,b%nglx
                    ip=G%intma(i,j,k,ie)
                    x(i,j,k)=G%coord(1,ip)
                    y(i,j,k)=G%coord(2,ip)
                    z(i,j,k)=G%coord(3,ip)
                end do
            end do
        end do

        ! Construct Mapping Derivatives: dx/dksi, dx/deta,dx/dzeta, dy/dksi, dy/deta, 
        ! dy/dzeta, dz/dksi, dz/deta, dz/dzeta
        call compute_local_gradient_quad_v3(x_ksiq,x_etaq,x_zetaq,x,b)
        call compute_local_gradient_quad_v3(y_ksiq,y_etaq,y_zetaq,y,b)
        call compute_local_gradient_quad_v3(z_ksiq,z_etaq,z_zetaq,z,b)

        !Construct Inverse Mapping
        do k=1,b%nqz
            do j=1,b%nqy
                do i=1,b%nqx
               
                    !set jacobian matrix J for 2D
                    if(b%nqx == 1) then
                        x_ksiq(i,j,k) = 1.0;  y_ksiq(i,j,k) = 0.0;  z_ksiq(i,j,k) = 0.0;
                        x_etaq(i,j,k) = 0.0;  x_zetaq(i,j,k) = 0.0;
                    endif
                    if(b%nqy == 1) then
                        x_etaq(i,j,k) = 0.0;  y_etaq(i,j,k) = 1.0;  z_etaq(i,j,k) = 0.0;
                        y_ksiq(i,j,k) = 0.0; y_zetaq(i,j,k) = 0.0;
                    endif
                    if(b%nqz == 1) then
                        ! Set the zeta direction to the unit outward normal of the element
                        ! surface (cross product of the two tangent vectors).
                        ! For flat Cartesian meshes z_ksiq=z_etaq=0 so this gives (0,0,1),
                        ! recovering the previous Cartesian behaviour.
                        ! For a curved surface such as the sphere z_ksiq and z_etaq are
                        ! non-zero, so the cross product gives the true outward unit normal
                        ! and xj = |t_ksi x t_eta| is the correct surface area element.
                        xn = y_ksiq(i,j,k)*z_etaq(i,j,k) - z_ksiq(i,j,k)*y_etaq(i,j,k)
                        yn = z_ksiq(i,j,k)*x_etaq(i,j,k) - x_ksiq(i,j,k)*z_etaq(i,j,k)
                        zn = x_ksiq(i,j,k)*y_etaq(i,j,k) - y_ksiq(i,j,k)*x_etaq(i,j,k)
                        rnorm = sqrt(xn*xn + yn*yn + zn*zn)
                        if (rnorm > 0.0) then
                            x_zetaq(i,j,k) = xn / rnorm
                            y_zetaq(i,j,k) = yn / rnorm
                            z_zetaq(i,j,k) = zn / rnorm
                        else
                            x_zetaq(i,j,k) = 0.0; y_zetaq(i,j,k) = 0.0; z_zetaq(i,j,k) = 1.0
                        end if
                    endif
              
                    !compute inverse of J
                    xj = (x_ksiq(i,j,k)*y_etaq(i,j,k)*z_zetaq(i,j,k) &
                        - x_ksiq(i,j,k)*y_zetaq(i,j,k)*z_etaq(i,j,k)) &
                        - (y_ksiq(i,j,k)*x_etaq(i,j,k)*z_zetaq(i,j,k) &
                        - y_ksiq(i,j,k)*x_zetaq(i,j,k)*z_etaq(i,j,k)) &
                        + (z_ksiq(i,j,k)*x_etaq(i,j,k)*y_zetaq(i,j,k) &
                        - z_ksiq(i,j,k)*x_zetaq(i,j,k)*y_etaq(i,j,k))

                    ksiq_x(i,j,k,ie)=  (y_etaq(i,j,k)*z_zetaq(i,j,k) &
                                        -y_zetaq(i,j,k)*z_etaq(i,j,k))/xj
                    ksiq_y(i,j,k,ie)= -(x_etaq(i,j,k)*z_zetaq(i,j,k) &
                                        -x_zetaq(i,j,k)*z_etaq(i,j,k))/xj
                    ksiq_z(i,j,k,ie)=  (x_etaq(i,j,k)*y_zetaq(i,j,k) &
                                        -x_zetaq(i,j,k)*y_etaq(i,j,k))/xj
                    etaq_x(i,j,k,ie)= -(y_ksiq(i,j,k)*z_zetaq(i,j,k) &
                                        -y_zetaq(i,j,k)*z_ksiq(i,j,k))/xj
                    etaq_y(i,j,k,ie)=  (x_ksiq(i,j,k)*z_zetaq(i,j,k) &
                                        -x_zetaq(i,j,k)*z_ksiq(i,j,k))/xj
                    etaq_z(i,j,k,ie)= -(x_ksiq(i,j,k)*y_zetaq(i,j,k) &
                                        -x_zetaq(i,j,k)*y_ksiq(i,j,k))/xj
                    zetaq_x(i,j,k,ie)= (y_ksiq(i,j,k)*z_etaq(i,j,k) &
                                        -y_etaq(i,j,k)*z_ksiq(i,j,k) )/xj
                    zetaq_y(i,j,k,ie)= -(x_ksiq(i,j,k)*z_etaq(i,j,k) &
                                        -x_etaq(i,j,k)*z_ksiq(i,j,k) )/xj
                    zetaq_z(i,j,k,ie)= (x_ksiq(i,j,k)*y_etaq(i,j,k) &
                                        -x_etaq(i,j,k)*y_ksiq(i,j,k) )/xj

                    ! Same rescaling as metrics.F90: x_zetaq/y_zetaq/z_zetaq above is the
                    ! *unit* surface normal (see rnorm block), so zetaq_x/y/z as computed is
                    ! dimensionless. The genuine grad(zeta) uses the true radial position
                    ! (magnitude earth_radius) as the zeta-direction tangent, rescaling by
                    ! 1/earth_radius. ksiq_x/etaq_x/jacq are unaffected (direction-only dependence).
                    if (b%nqz == 1 .and. G%is_sphere) then
                        zetaq_x(i,j,k,ie) = zetaq_x(i,j,k,ie) / earth_radius
                        zetaq_y(i,j,k,ie) = zetaq_y(i,j,k,ie) / earth_radius
                        zetaq_z(i,j,k,ie) = zetaq_z(i,j,k,ie) / earth_radius
                    end if

                    jacq(i,j,k,ie) = b%wnqx(i)*b%wnqy(j)*b%wnqz(k)*abs(xj)
                    xjacq(i,j,k,ie) = xj
              
                    !fix inverse jacobian matrix for 2D
                    if(b%nqx == 1)  ksiq_x(i,j,k,ie) = 0.0
                    if(b%nqy == 1)  etaq_y(i,j,k,ie) = 0.0
                    if(b%nqz == 1 .and. .not. G%is_sphere)  zetaq_z(i,j,k,ie) = 0.0
              
                end do
            end do
        end do

    end do !ie
end subroutine metrics_quad
