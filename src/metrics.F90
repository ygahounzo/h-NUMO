!----------------------------------------------------------------------!
!>@brief This subroutine builds the METRIC TERMS see Giraldo 2001 (IJNMF)
!>@author  Francis X. Giraldo on 7/08
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@date 11 September 2009  James F. Kelly
!> Generalized to 3D
!> Combined to do both Cube and Sphere domains.
!----------------------------------------------------------------------!
subroutine metrics(ksi_x,ksi_y,ksi_z,eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,jac,xjac)

    use mod_basis, only: nglx, ngly, nglz, wglx, wgly, wglz, nqx,nqy,nqz,wnqx,wnqy,wnqz
 
    use mod_grid, only: intma, coord, nelem, intma_dg_quad

    use mod_input, only: geometry_type, lp4est, lp6est

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
    logical is_sphere

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

    is_sphere = .false.
    if(geometry_type(1:6) == 'sphere') is_sphere = .true.
  
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
        if (.not. is_sphere) then
            call compute_local_gradient_v3(x_ksi,x_eta,x_zeta,x,nglx,ngly,nglz,ndim)
            call compute_local_gradient_v3(y_ksi,y_eta,y_zeta,y,nglx,ngly,nglz,ndim)
            call compute_local_gradient_v3(z_ksi,z_eta,z_zeta,z,nglx,ngly,nglz,ndim)
        else
           if ( (lp4est .or. lp6est) .and. geometry_type(1:10) == 'sphere_hex') then
              call semi_analytic_metrics_nep(ie,x,y,z,nglx,ngly,nglz,x_ksi,x_eta,x_zeta,y_ksi,y_eta,y_zeta,z_ksi,z_eta,z_zeta)
!              call semi_analytic_metrics(x,y,z,nglx,ngly,nglz,x_ksi,x_eta,x_zeta,y_ksi,y_eta,y_zeta,z_ksi,z_eta,z_zeta)
           else 
              call semi_analytic_metrics(x,y,z,nglx,ngly,nglz,x_ksi,x_eta,x_zeta,y_ksi,y_eta,y_zeta,z_ksi,z_eta,z_zeta)
!            call compute_local_gradient_v3(x_ksi,x_eta,x_zeta,x,nglx,ngly,nglz,ndim)
!            call compute_local_gradient_v3(y_ksi,y_eta,y_zeta,y,nglx,ngly,nglz,ndim)
!            call compute_local_gradient_v3(z_ksi,z_eta,z_zeta,z,nglx,ngly,nglz,ndim)
           end if
        end if

        !Construct Inverse Mapping
        do k=1,nglz
            do j=1,ngly
                do i=1,nglx
               
                    !set jacobian matrix J for 2D
                    if(.not. is_sphere) then
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
                    endif
              
                    !compute inverse of J
                    !xj = x_ksi(i,j,k)*y_eta(i,j,k)*z_zeta(i,j,k) - x_ksi(i,j,k)*y_zeta(i,j,k)*z_eta(i,j,k) - &
                    !     y_ksi(i,j,k)*x_eta(i,j,k)*z_zeta(i,j,k) + y_ksi(i,j,k)*x_zeta(i,j,k)*z_eta(i,j,k) + &
                    !     z_ksi(i,j,k)*x_eta(i,j,k)*y_zeta(i,j,k) - z_ksi(i,j,k)*x_zeta(i,j,k)*y_eta(i,j,k)
                    !ksi_x(i,j,k,ie)=  (y_eta(i,j,k)*z_zeta(i,j,k)-y_zeta(i,j,k)*z_eta(i,j,k))/xj
                    !ksi_y(i,j,k,ie)= -(x_eta(i,j,k)*z_zeta(i,j,k)-x_zeta(i,j,k)*z_eta(i,j,k))/xj
                    !ksi_z(i,j,k,ie)=  (x_eta(i,j,k)*y_zeta(i,j,k)-x_zeta(i,j,k)*y_eta(i,j,k))/xj
                    !eta_x(i,j,k,ie)= -(y_ksi(i,j,k)*z_zeta(i,j,k)-y_zeta(i,j,k)*z_ksi(i,j,k))/xj
                    !eta_y(i,j,k,ie)=  (x_ksi(i,j,k)*z_zeta(i,j,k)-x_zeta(i,j,k)*z_ksi(i,j,k))/xj
                    !eta_z(i,j,k,ie)= -(x_ksi(i,j,k)*y_zeta(i,j,k)-x_zeta(i,j,k)*y_ksi(i,j,k))/xj
                    !zeta_x(i,j,k,ie)= (y_ksi(i,j,k)*z_eta(i,j,k) -y_eta(i,j,k)*z_ksi(i,j,k) )/xj
                    !zeta_y(i,j,k,ie)=-(x_ksi(i,j,k)*z_eta(i,j,k) -x_eta(i,j,k)*z_ksi(i,j,k) )/xj
                    !zeta_z(i,j,k,ie)= (x_ksi(i,j,k)*y_eta(i,j,k) -x_eta(i,j,k)*y_ksi(i,j,k) )/xj
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
                    if(.not. is_sphere) then
                        if(nglx == 1)  ksi_x(i,j,k,ie) = 0
                        if(ngly == 1)  eta_y(i,j,k,ie) = 0
                        if(nglz == 1)  zeta_z(i,j,k,ie) = 0
                    endif
              
                end do
            end do
        end do

    end do !ie
end subroutine metrics

subroutine metrics1(ksi_x,ksi_y,eta_x,eta_y,jac,xjac)

    use mod_basis, only: nglx, ngly, nglz, wglx, wgly, wglz, nqx,nqy,nqz,wnqx,wnqy,wnqz
 
    use mod_grid, only: intma, coord, nelem, intma_dg_quad

    use mod_input, only: geometry_type, lp4est, lp6est

    use mod_interface, only: compute_local_gradient_v3, compute_local_gradient_quad_v3

    implicit none

    !global arrays
    real, dimension(nglx,ngly,nglz,nelem), intent(out) :: ksi_x,  ksi_y
    real, dimension(nglx,ngly,nglz,nelem), intent(out) :: eta_x,  eta_y
    real, dimension(nglx,ngly,nglz,nelem), intent(out) :: jac, xjac

    !local arrays
    real, dimension(nglx,ngly,nglz) :: x, y, z
    real, dimension(nglx,ngly,nglz) :: x_ksi, x_eta, x_zeta
    real, dimension(nglx,ngly,nglz) :: y_ksi, y_eta, y_zeta
    real, dimension(nglx,ngly,nglz) :: z_ksi, z_eta, z_zeta

    real xj
    integer ie, i, j, k
    integer ip, ndim
    logical is_sphere

    !Define Dimension of Matrix going into MXM routine
    ndim=1

    !initialize the global matrix
    ksi_x=0
    ksi_y=0
    eta_x=0
    eta_y=0
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
                    if(.not. is_sphere) then
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
                    endif
              
                    !compute inverse of J

                    xj = x_ksi(i,j,k)*y_eta(i,j,k) - y_ksi(i,j,k)*x_eta(i,j,k)

                    ksi_x(i,j,k,ie)=  y_eta(i,j,k)/xj
                    ksi_y(i,j,k,ie)= -x_eta(i,j,k)/xj
                    eta_x(i,j,k,ie)= -y_ksi(i,j,k)/xj
                    eta_y(i,j,k,ie)=  x_ksi(i,j,k)/xj

                    jac(i,j,k,ie)=wglx(i)*wgly(j)*abs(xj)
                    xjac(i,j,k,ie)=xj
              
                end do
            end do
        end do

    end do !ie
end subroutine metrics1

!----------------------------------------------------------------------!
!@brief Semi analytic metrics for the sphere
!----------------------------------------------------------------------!
subroutine semi_analytic_metrics(x,y,z,nglx,ngly,nglz,x_ksi,x_eta,x_zeta,&
    y_ksi,y_eta,y_zeta,z_ksi,z_eta,z_zeta)

    use mod_basis, only: dpsix, dpsiy, dpsiz

    use mod_constants, only: pi

    use mod_interface, only: compute_local_gradient_v3

    implicit none

    integer, intent(in) :: nglx, ngly, nglz
    real, intent(inout), dimension(nglx,ngly,nglz) :: x, y, z

    real, dimension(nglx,ngly,nglz),  intent(out) ::  &
        x_ksi,  y_ksi,  z_ksi,  &
        x_eta,  y_eta,  z_eta,  &
        x_zeta, y_zeta, z_zeta

    integer i, j, k, rot, ndim
    real, dimension(nglx,ngly,nglz) :: &
        r, x_r, y_r, z_r, r_ksi, &
        r_eta, r_zeta, phi, lambda, &
        x_p, y_p, z_p, x_l, y_l, z_l, &
        p_ksi, p_eta, p_zeta, l_ksi, l_eta, l_zeta, &
        tmp1, tmp2, tmp3

    !Define Dimension of Matrix going into MXM routine
    ndim=1

    do k=1,nglz
        do j=1,ngly
            do i=1,nglx
                call cart2sph(x(i,j,k),y(i,j,k),z(i,j,k), &
                    r(i,j,k),phi(i,j,k),lambda(i,j,k))
            end do !i
        end do !j
    end do !k

    rot=0
    if(any(phi>(pi/4.0))) then
        rot = 1
    elseif(any(phi<(-pi/4.0))) then
        rot = 2
    elseif(any(abs(lambda) > (3.0 * pi/4.0))) then
        rot = 3
    end if

    do k=1,nglz
        do j=1,ngly
            do i=1,nglx

                if(rot==1) then
                    tmp1(i,j,k)=x(i,j,k)
                    x(i,j,k)=z(i,j,k)
                    z(i,j,k)=-tmp1(i,j,k)
                    call cart2sph(x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k),phi(i,j,k)&
                        ,lambda(i,j,k))
                elseif(rot==2) then
                    tmp1(i,j,k)=x(i,j,k)
                    x(i,j,k)=-z(i,j,k)
                    z(i,j,k)=tmp1(i,j,k)
                    call cart2sph(x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k),phi(i,j,k)&
                        ,lambda(i,j,k))
                elseif(rot==3) then
                    x(i,j,k)=-x(i,j,k)
                    y(i,j,k)=-y(i,j,k)
                    call cart2sph(x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k),phi(i,j,k)&
                        ,lambda(i,j,k))
                endif

                x_r(i,j,k)=cos(phi(i,j,k))*cos(lambda(i,j,k))
                y_r(i,j,k)=cos(phi(i,j,k))*sin(lambda(i,j,k))
                z_r(i,j,k)=sin(phi(i,j,k))
                x_p(i,j,k)=-r(i,j,k)*sin(phi(i,j,k))*cos(lambda(i,j,k))
                y_p(i,j,k)=-r(i,j,k)*sin(phi(i,j,k))*sin(lambda(i,j,k))
                z_p(i,j,k)=r(i,j,k)*cos(phi(i,j,k))
                x_l(i,j,k)=-r(i,j,k)*cos(phi(i,j,k))*sin(lambda(i,j,k))
                y_l(i,j,k)=r(i,j,k)*cos(phi(i,j,k))*cos(lambda(i,j,k))
                z_l(i,j,k)=0.0

            end do !i
        end do !j
    end do !k

    !Construct Mapping Derivatives: dx/dksi, dx/deta,dx/dzeta, dy/dksi, dy/deta, dy/dzeta, dz/dksi, dz/deta, dz/dzeta
    call compute_local_gradient_v3(r_ksi,r_eta,r_zeta,r,     nglx,ngly,nglz,ndim)
    call compute_local_gradient_v3(p_ksi,p_eta,p_zeta,phi,   nglx,ngly,nglz,ndim)
    call compute_local_gradient_v3(l_ksi,l_eta,l_zeta,lambda,nglx,ngly,nglz,ndim)

    do k=1,nglz
        do j=1,ngly
            do i=1,nglx
              
                !set jacobian matrix J for 2D
                if(nglz == 1) then
                    p_zeta(i,j,k) = 0; l_zeta(i,j,k) = 0; r_zeta(i,j,k) = 1e13;
                    r_ksi(i,j,k) = 0; r_eta(i,j,k) = 0;
                endif
              
                x_ksi(i,j,k)=x_r(i,j,k)*r_ksi(i,j,k)+x_p(i,j,k)*p_ksi(i,j,k)&
                    +x_l(i,j,k)*l_ksi(i,j,k)
                x_eta(i,j,k)=x_r(i,j,k)*r_eta(i,j,k)+x_p(i,j,k)*p_eta(i,j,k)&
                    +x_l(i,j,k)*l_eta(i,j,k)
                x_zeta(i,j,k)=x_r(i,j,k)*r_zeta(i,j,k)+x_p(i,j,k)*p_zeta(i,j,k)&
                    +x_l(i,j,k)*l_zeta(i,j,k)

                y_ksi(i,j,k)=y_r(i,j,k)*r_ksi(i,j,k)+y_p(i,j,k)*p_ksi(i,j,k)&
                    +y_l(i,j,k)*l_ksi(i,j,k)
                y_eta(i,j,k)=y_r(i,j,k)*r_eta(i,j,k)+y_p(i,j,k)*p_eta(i,j,k)&
                    +y_l(i,j,k)*l_eta(i,j,k)
                y_zeta(i,j,k)=y_r(i,j,k)*r_zeta(i,j,k)+y_p(i,j,k)*p_zeta(i,j,k)&
                    +y_l(i,j,k)*l_zeta(i,j,k)

                z_ksi(i,j,k)=z_r(i,j,k)*r_ksi(i,j,k)+z_p(i,j,k)*p_ksi(i,j,k)&
                    +z_l(i,j,k)*l_ksi(i,j,k)
                z_eta(i,j,k)=z_r(i,j,k)*r_eta(i,j,k)+z_p(i,j,k)*p_eta(i,j,k)&
                    +z_l(i,j,k)*l_eta(i,j,k)
                z_zeta(i,j,k)=z_r(i,j,k)*r_zeta(i,j,k)+z_p(i,j,k)*p_zeta(i,j,k)&
                    +z_l(i,j,k)*l_zeta(i,j,k)
                
            end do !i
        end do !j
    end do !k

    if(rot==1)then

        tmp1=x_ksi
        tmp2=x_eta
        tmp3=x_zeta

        x_ksi=-z_ksi
        x_eta=-z_eta
        x_zeta=-z_zeta

        z_ksi=tmp1
        z_eta=tmp2
        z_zeta=tmp3

    elseif(rot==2)then

        tmp1=x_ksi
        tmp2=x_eta
        tmp3=x_zeta

        x_ksi=z_ksi
        x_eta=z_eta
        x_zeta=z_zeta

        z_ksi=-tmp1
        z_eta=-tmp2
        z_zeta=-tmp3

    elseif(rot==3)then

        x_ksi=-x_ksi
        x_eta=-x_eta
        x_zeta=-x_zeta

        y_ksi=-y_ksi
        y_eta=-y_eta
        y_zeta=-y_zeta

    endif

    return
end subroutine semi_analytic_metrics

!----------------------------------------------------------------------!
subroutine semi_analytic_metrics_nep(ie,x,y,z,nglx,ngly,nglz, &
                                 x_ksi,x_eta,x_zeta,      &
                                 y_ksi,y_eta,y_zeta,      &
                                 z_ksi,z_eta,z_zeta)

  use mod_constants,   only: pi, earth_radius
  use mod_input,       only: geometry_type
  use mod_interface,   only: compute_local_gradient_v3
  use mod_p4est,       only: cs_face

  implicit none

  logical, parameter :: shallow_atmos = .false.

  integer, intent(in) :: ie
  real, intent(inout), dimension(nglx,ngly,nglz) :: x, y, z
  integer, intent(in) :: nglx, ngly, nglz

  real, dimension(nglx,ngly,nglz),  intent(out) ::  &
       x_ksi,  y_ksi,  z_ksi,  &
       x_eta,  y_eta,  z_eta,  &
       x_zeta, y_zeta, z_zeta

  integer :: i, j, k, l, info, ipiv(3)
  real :: lon(6), lat(6)
  real :: rot(3,3), vec(3), theta
  real :: irot(3,3), wrk(9)
  real :: x_r, y_r, z_r, x_p, y_p, z_p, x_l, y_l, z_l, &
                   pp_ksi, pp_eta, pp_zeta, lp_ksi, lp_eta, lp_zeta, &
                   tmpx, tmpy, pp_xg, pp_yg, lp_xg, lp_yg, &
                   l_lp, l_pp, p_lp, p_pp, &
                   x_lp, x_pp, y_lp, y_pp, z_lp, z_pp, &
                   xg_xgs, xg_ygs, yg_xgs, yg_ygs
  real, dimension(nglx,ngly,nglz) :: &
       r, phi, lambda, xg, yg, xg_ksi, xg_eta, xg_zeta, &
       yg_ksi, yg_eta, yg_zeta, r_ksi, r_eta, r_zeta
  real, allocatable, dimension(:,:,:) :: xgs, ygs
  integer :: ndim

  integer :: grid_proj
  real    :: f_gnom1, f_gnom2, b_gnom1
  real    :: xx, yy, zz, lon_f, lat_f

  if(geometry_type=='sphere_hex_es') then
    grid_proj=0
  else
    grid_proj=1
  end if

  if(grid_proj==1) allocate(xgs(nglx,ngly,nglz),ygs(nglx,ngly,nglz))

  !Initialize MXM dimension
  ndim = 1

  ! Determine face of cubed sphere on which element lies, rotate to face
  ! centered at (lat,lon)=(0,0)

  lon(:) = 0.0
  lon(2) = 90.0
  lon(3) = 180.0
  lon(4) = 270.0

  lat(:) = 0.0
  lat(5) = 90.0
  lat(6) = -90.0

  do i = 1,6
    lon(i) = lon(i)*pi/180.0
    lat(i) = lat(i)*pi/180.0
  end do !i

  l = cs_face(ie)
  rot(:,:) = 0.0
  if(l<5)then
    theta = lon(l)
    rot(1,1) = cos(theta)
    rot(1,2) = sin(theta)
    rot(2,1) = -sin(theta)
    rot(2,2) = cos(theta)
    rot(3,3) = 1.0
  else
    theta = lat(l)
    rot(1,1) = cos(theta)
    rot(1,3) = sin(theta)
    rot(2,2) = 1.0
    rot(3,1) = -sin(theta)
    rot(3,3) = cos(theta)
  endif

  irot = rot
  !if(r8==r8_fix) then
    call dgetrf(3,3,irot,3,ipiv,info)
    call dgetri(3,irot,3,ipiv,wrk,9,info)
  !else
  !  call sgetrf(3,3,irot,3,ipiv,info)
  !  call sgetri(3,irot,3,ipiv,wrk,9,info)
  !end if

  do k = 1,nglz
    do j = 1,ngly
      do i = 1,nglx

        vec(1) = x(i,j,k)
        vec(2) = y(i,j,k)
        vec(3) = z(i,j,k)

! Rotate all points to Face #1 of the cubed sphere

        x(i,j,k) = dot_product(rot(1,:),vec(:))
        y(i,j,k) = dot_product(rot(2,:),vec(:))
        z(i,j,k) = dot_product(rot(3,:),vec(:))

! Convert from Cartesian to spherical coordinates

        call cart2sph(x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k),phi(i,j,k)&
                      ,lambda(i,j,k))

! Forward gnomonic projection from spherical coordinates to 2D local
! Cartesian plane

        xg(i,j,k) = f_gnom1(lambda(i,j,k),phi(i,j,k))
        yg(i,j,k) = f_gnom2(lambda(i,j,k),phi(i,j,k))

        if(grid_proj==1)then

! Backwards projection to equiangular 2D plane

          xgs(i,j,k) = b_gnom1(xg(i,j,k))
          ygs(i,j,k) = b_gnom1(yg(i,j,k))

        endif
      end do !i
    end do !j
  end do !k

! Compute canonical derivatives of original 2D plane (xg,yg)

  call compute_local_gradient_v3(r_ksi,r_eta,r_zeta,r,nglx,ngly,nglz,ndim)

  if(shallow_atmos) r(:,:,:) = earth_radius

  if(grid_proj==1)then
    call compute_local_gradient_v3(xg_ksi,xg_eta,xg_zeta,xgs,nglx,ngly,nglz,ndim)
    call compute_local_gradient_v3(yg_ksi,yg_eta,yg_zeta,ygs,nglx,ngly,nglz,ndim)
  else
    call compute_local_gradient_v3(xg_ksi,xg_eta,xg_zeta,xg,nglx,ngly,nglz,ndim)
    call compute_local_gradient_v3(yg_ksi,yg_eta,yg_zeta,yg,nglx,ngly,nglz,ndim)
  endif

  do k = 1,nglz
    do j = 1,ngly
      do i = 1,nglx

! Analytic Jacobian from Cartesian (x,y,z) to spherical (l,p,r) on Face #1

        x_r = cos(phi(i,j,k))*cos(lambda(i,j,k))
        y_r = cos(phi(i,j,k))*sin(lambda(i,j,k))
        z_r = sin(phi(i,j,k))

        x_p = -r(i,j,k)*sin(phi(i,j,k))*cos(lambda(i,j,k))
        y_p = -r(i,j,k)*sin(phi(i,j,k))*sin(lambda(i,j,k))
        z_p = r(i,j,k)*cos(phi(i,j,k))

        x_l = -r(i,j,k)*cos(phi(i,j,k))*sin(lambda(i,j,k))
        y_l = r(i,j,k)*cos(phi(i,j,k))*cos(lambda(i,j,k))
        z_l = 0.0

! Backwards rotation of analytic derivatives from Face #1
! to Face containing this element

        vec(1)=x_l
        vec(2)=y_l
        vec(3)=z_l
        x_lp=dot_product(irot(1,:),vec(:))
        y_lp=dot_product(irot(2,:),vec(:))
        z_lp=dot_product(irot(3,:),vec(:))

        vec(1)=x_p
        vec(2)=y_p
        vec(3)=z_p
        x_pp=dot_product(irot(1,:),vec(:))
        y_pp=dot_product(irot(2,:),vec(:))
        z_pp=dot_product(irot(3,:),vec(:))

        vec(1)=x_r
        vec(2)=y_r
        vec(3)=z_r
        x_r=dot_product(irot(1,:),vec(:))
        y_r=dot_product(irot(2,:),vec(:))
        z_r=dot_product(irot(3,:),vec(:))

        if(grid_proj==1)then

! Analytic Jacobian of forward projection from 2D equiangular plane to local
! Cartesian plane

          xg_xgs=1.0/(cos(xgs(i,j,k))**2)
          xg_ygs=0.0
          yg_xgs=0.0
          yg_ygs=1.0/(cos(ygs(i,j,k))**2)

          tmpx = xg_xgs*xg_ksi(i,j,k)+xg_ygs*yg_ksi(i,j,k)
          tmpy = yg_xgs*xg_ksi(i,j,k)+yg_ygs*yg_ksi(i,j,k)
          xg_ksi(i,j,k) = tmpx
          yg_ksi(i,j,k) = tmpy

          tmpx = xg_xgs*xg_eta(i,j,k)+xg_ygs*yg_eta(i,j,k)
          tmpy = yg_xgs*xg_eta(i,j,k)+yg_ygs*yg_eta(i,j,k)
          xg_eta(i,j,k) = tmpx
          yg_eta(i,j,k) = tmpy

          tmpx = xg_xgs*xg_zeta(i,j,k)+xg_ygs*yg_zeta(i,j,k)
          tmpy = yg_xgs*xg_zeta(i,j,k)+yg_ygs*yg_zeta(i,j,k)
          xg_zeta(i,j,k) = tmpx
          yg_zeta(i,j,k) = tmpy

        endif

! Analytic Jacobian of backwards gnomonic projection from local Cartesian
! plane to spherical coordinates on Face #1

        lp_xg=1.0/(1.0+xg(i,j,k)**2)
        lp_yg=0.0
        pp_xg=-yg(i,j,k)*sin(lambda(i,j,k))*lp_xg/(1.0+(yg(i,j,k)*cos(lambda(i,j,k)))**2)
        pp_yg=cos(lambda(i,j,k))/(1.0+(yg(i,j,k)*cos(lambda(i,j,k)))**2)

! Chain rule to transform from planar to Cartesian canonical derivatives

        lp_ksi  = lp_xg*xg_ksi(i,j,k) + lp_yg*yg_ksi(i,j,k)
        lp_eta  = lp_xg*xg_eta(i,j,k) + lp_yg*yg_eta(i,j,k)
        lp_zeta = lp_xg*xg_zeta(i,j,k)+ lp_yg*yg_zeta(i,j,k)

        pp_ksi  = pp_xg*xg_ksi(i,j,k) + pp_yg*yg_ksi(i,j,k)
        pp_eta  = pp_xg*xg_eta(i,j,k) + pp_yg*yg_eta(i,j,k)
        pp_zeta = pp_xg*xg_zeta(i,j,k)+ pp_yg*yg_zeta(i,j,k)

        x_ksi(i,j,k) = x_r*r_ksi(i,j,k) + x_pp*pp_ksi + x_lp*lp_ksi
        x_eta(i,j,k) = x_r*r_eta(i,j,k) + x_pp*pp_eta + x_lp*lp_eta
        x_zeta(i,j,k)= x_r*r_zeta(i,j,k)+ x_pp*pp_zeta+ x_lp*lp_zeta

        y_ksi(i,j,k) = y_r*r_ksi(i,j,k) + y_pp*pp_ksi + y_lp*lp_ksi
        y_eta(i,j,k) = y_r*r_eta(i,j,k) + y_pp*pp_eta + y_lp*lp_eta
        y_zeta(i,j,k)= y_r*r_zeta(i,j,k)+ y_pp*pp_zeta+ y_lp*lp_zeta

        z_ksi(i,j,k) = z_r*r_ksi(i,j,k) + z_pp*pp_ksi + z_lp*lp_ksi
        z_eta(i,j,k) = z_r*r_eta(i,j,k) + z_pp*pp_eta + z_lp*lp_eta
        z_zeta(i,j,k)= z_r*r_zeta(i,j,k)+ z_pp*pp_zeta+ z_lp*lp_zeta

      end do !i
    end do !j
  end do !k

end subroutine semi_analytic_metrics_nep
