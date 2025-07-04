!-------------------------------------------------------------
!> @brief Module for basis functions and interpolation points
!-------------------------------------------------------------
module mod_legendre

    use mod_types, only : r8
    use mod_types, only : rQ=>r8

    implicit none

    real(kind=rQ), parameter :: thres = epsilon(1.0_rQ)

    logical, parameter :: reduce_round_off = .true.

    !---------------------------------
    !> @brief Legendre polynomial type
    !----------------------------------
    type legendre_poly_type
        real(rQ) :: p0, p0_1, p0_2, &
            p1, p1_1, p1_2, &
            p2, p2_1, p2_2, &
            p00,p00_1,p00_2
    end type legendre_poly_type

    public :: legendre_gauss_lobatto, &
        legendre_gauss, &
        legendre_poly,  &
        legendre_basis, &
        lagrange_basis, &
        lagrange_basis_filter,&
        lagrange_basis2,&
        legendre_poly_type

contains

    !--------------------------------------
    !> @brief Compute legendre polynomial
    !--------------------------------------
    subroutine legendre_poly(n,x,lpoly)
        integer,  intent(in) :: n
        real(r8), intent(in) :: x
        type(legendre_poly_type), intent(out) :: lpoly

        call legendre_poly_loc(n,real(x,rQ),lpoly)
    end subroutine legendre_poly

    !---------------------------------------------------------------------!
    !>@brief This subroutine finds the Legendre-Gauss-Lobatto Roots.
    !>@author  Francis X. Giraldo on 7/08
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !---------------------------------------------------------------------!
    subroutine legendre_gauss_lobatto(ngl,xgl,wgl)
        integer, intent(in) :: ngl
        real(r8), intent(out) :: xgl(ngl), wgl(ngl)

        integer, parameter :: kmax=20

        integer :: n, nh,k, i

        real(rQ) :: pi, x, dx
        real(rQ) :: p0, p0_1, p0_2
        type(legendre_poly_type) :: lpoly

        !nop=0
        if(ngl == 1) then
            xgl(1) = 0
            wgl(1) = 2
            return
        endif
    
        pi=4.0_rQ*atan(1.0_rQ)
        n=ngl-1
        nh=(n+1)/2

        !First Find Half of the Roots
        do i=1,nh
            x=cos( (2.0_rQ*i - 1.)/(2.0_rQ*n + 1.0_rQ)*pi )

            do k=1,kmax
                !Construct Legendre Polynomial and Derivatives
                call legendre_poly_loc(n,x,lpoly)
                p0=lpoly%p0 ; p0_1=lpoly%p0_1 ; p0_2=lpoly%p0_2

                !Get next Newton Iterative
                dx = -(1.0_rQ-x**2)*p0_1/(-2.0_rQ*x*p0_1 + (1.0_rQ-x**2)*p0_2)
                x  = x + dx
                if (abs(dx) < thres) exit
            end do !K loop

            xgl(n+2-i)=real(x,kind=r8)
            wgl(n+2-i)=real(2.0_rQ/( real(n*(n+1),kind=rQ)*p0**2 ),kind=r8)
        end do !I loop

        !Check for Zero
        if (n+1 /= 2*nh) then
            x=0.0_rQ
            call legendre_poly_loc(n,x,lpoly)
            p0=lpoly%p0

            xgl(nh+1)=real(x,r8)
            wgl(nh+1)=real(2.0_rQ/( real(n*(n+1),kind=rQ)*p0**2 ),r8)
        endif

        !Find Remainder of Roots via Symmetry
        do i=1,nh
            xgl(i)=-xgl(n+2-i)
            wgl(i)=+wgl(n+2-i)
        end do
    end subroutine legendre_gauss_lobatto

    !---------------------------------------------------------------------!
    !>@brief This subroutine finds the Legendre-Gauss Roots.
    !>@author  Francis X. Giraldo on 7/08
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !---------------------------------------------------------------------!
    subroutine legendre_gauss(ngl,xgl,wgl)
        integer, intent(in) :: ngl
        real(r8), intent(out) :: xgl(ngl), wgl(ngl)

        integer, parameter :: kmax=20

        integer :: n, nh,k, i

        real(rQ) :: pi, x, dx
        real(rQ) :: p00, p00_1

        type(legendre_poly_type) :: lpoly

        !nop=0
        if(ngl == 1) then
            xgl(1) = 0
            wgl(1) = 2
            return
        endif
  
        pi=4.0_rQ*atan(1.0_rQ)
        n=ngl-1
        nh=(n+1)/2

        !First Find Half of the Roots
        do i=1,nh
            x=cos( (2.0_rQ*i - 1.)/(2.0_rQ*n + 1.0_rQ)*pi )
            do k=1,kmax
                !Construct Legendre Polynomial and Derivatives
                call legendre_poly_loc(n,x,lpoly)
                p00=lpoly%p00 ; p00_1=lpoly%p00_1

                !Get next Newton Iterative
                dx = -p00/p00_1
                x  = x + dx
                if (abs(dx) < thres) exit
            end do !K loop

            xgl(n+2-i) = real(x,kind=r8)
            wgl(n+2-i) = real(2.0_rQ/( (1.0_rQ-x**2)*p00_1**2 ),kind=r8)

        end do !I loop

        !Check for Zero
        if (n+1 /= 2*nh) then
            x=0.0_rQ
            call legendre_poly_loc(n,x,lpoly)
            p00_1=lpoly%p00_1

            xgl(nh+1) = real(x,kind=r8)
            wgl(nh+1) = real(2.0_rQ/( (1.0_rQ-x**2)*p00_1**2 ),kind=r8)
        endif

        !Find Remainder of Roots via Symmetry
        do i=1,nh
            xgl(i) = -xgl(n+2-i)
            wgl(i) = +wgl(n+2-i)
        end do
   
    end subroutine legendre_gauss

    !---------------------------------------------------------------------!
    !>@brief This subroutine finds the Legendre Cardinal Basis Functions and
    !>their Derivatives.
    !>@author  Francis X. Giraldo on 7/08
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !---------------------------------------------------------------------!
    subroutine legendre_poly_loc(n,x,lpoly)
        integer,  intent(in) :: n
        real(rQ), intent(in) :: x
        type(legendre_poly_type), intent(out) :: lpoly

        real(rQ) :: p0, p0_1, p0_2, &
            p1, p1_1, p1_2, &
            p2, p2_1, p2_2, &
            p00,p00_1,p00_2
        real(rQ) :: a, b
        integer  :: j

        !Construct Legendre Polynomials 00=N+1, 0=N, 1=N-1, 2=N-2, and derivs
        p2=0.0_rQ  ; p2_1=0.0_rQ  ; p2_2=0.0_rQ
        p1=0.0_rQ  ; p1_1=0.0_rQ  ; p1_2=0.0_rQ
        p0=1.0_rQ  ; p0_1=0.0_rQ  ; p0_2=0.0_rQ
        p00=1.0_rQ ; p00_1=0.0_rQ ; p00_2=0.0_rQ

        !Construct Nth Order Legendre Polynomial
        do j=1,n
            p2=p1
            p2_1=p1_1
            p2_2=p1_2

            p1=p0
            p1_1=p0_1
            p1_2=p0_2

            a=(2.0_rQ*real(j,rQ)-1.0_rQ)/real(j,rQ)
            b=(real(j,rQ)-1.0_rQ)/real(j,rQ)
            p0=a*x*p1 - b*p2
            p0_1=a*( p1 + x*p1_1 ) - b*p2_1
            p0_2=a*( 2.0_rQ*p1_1 + x*p1_2 ) - b*p2_2

            a=(2.0_rQ*real(j,rQ)+1.0_rQ)/(real(j,rQ)+1.0_rQ)
            b=(real(j,rQ))/(real(j,rQ)+1.0_rQ)
            p00=a*x*p0 - b*p1
            p00_1=a*( p0 + x*p0_1 ) - b*p1_1
            p00_2=a*( 2.0_rQ*p0_1 + x*p0_2 ) - b*p1_2
        end do !J loop

        lpoly%p0   = p0 ; lpoly%p0_1 = p0_1 ; lpoly%p0_2 = p0_2
        lpoly%p1   = p1 ; lpoly%p1_1 = p1_1 ; lpoly%p1_2 = p1_2
        lpoly%p2   = p2
        lpoly%p2_1 = p2_1
        lpoly%p2_2 = p2_2
        lpoly%p00  = p00; lpoly%p00_1= p00_1; lpoly%p00_2= p00_2

    end subroutine legendre_poly_loc

    !---------------------------------------------------------------------!
    !>@brief This subroutine finds the Legendre Cardinal Basis Functions and
    !>their Derivatives.
    !>@author  Francis X. Giraldo on 7/08
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !---------------------------------------------------------------------!

    subroutine legendre_basis(ngl,xgl,psi,dpsi)
        integer,  intent(in)  :: ngl
        real(r8), intent(in)  :: xgl(ngl)
        real(r8), intent(out) :: psi(ngl,ngl), dpsi(ngl,ngl)

        integer  :: n, nh, i,j
        real(rQ) :: xj, ksi
        real(rQ) :: bb(ngl), cc(ngl)

        type(legendre_poly_type) :: lpoly_i, lpoly_j

        n=ngl-1

        if(.not.reduce_round_off) then
            do j=1,ngl
                xj=real(xgl(j),rQ)
                call legendre_poly_loc(n,xj,lpoly_j)
                do i=1,ngl
                    ksi=real(xgl(i),rQ)
                    call legendre_poly_loc(n,ksi,lpoly_i)
                    if (i == j) then
                        psi(i,j)=1.0_r8
                    else
                        psi(i,j)=0.0_r8
                    endif
            
                    if (i == j) then
                        if (i /= 1 .and. i /= ngl) then
                            dpsi(i,j)=0.0_r8
                        else if (i == 1) then
                            dpsi(i,j)=-real(n*(n+1),r8)/4.0_r8
                        else if (i == ngl) then
                            dpsi(i,j)=+real(n*(n+1),r8)/4.0_r8
                        endif
                    else if (i /= j) then
                        dpsi(i,j)=real(lpoly_j%p0/(lpoly_i%p0*(xj-ksi)),r8)
                    endif

                end do
            end do
        else
            bb(:) = 0.0_rQ
            cc(:) = 0.0_rQ
            do j=1,ngl
                xj=real(xgl(j),rQ)
                do i=1,ngl
                    ksi=real(xgl(i),rQ)
                    if (i==j) then
                        psi(i,j) = 1.0_r8
                    else
                        bb(j) = bb(j) + real(log(abs(xj - ksi)),rQ)
                        psi(i,j) = 0.0_r8
                    end if
                end do
            end do

            do j=1,ngl
                xj=real(xgl(j),rQ)
                do i=1,ngl
                    ksi=real(xgl(i),rQ)
                    if (i /= j) then
                        dpsi(i,j) = (-1)**(j+i)*exp(bb(j)-bb(i))/(xj-ksi)
                        cc(j) = cc(j) + dpsi(i,j)
                    end if
                end do
            end do

            do j=1,ngl
                do i=1,ngl
                    if(i==j) dpsi(i,j) = real(-cc(j),r8)
                end do
            end do
        end if
   
    end subroutine legendre_basis

    !---------------------------------------------------------------------!
    !>@brief This subroutine finds the Legendre Cardinal Basis Functions and
    !>their Derivatives in in Lagrange Polynomial Form.
    !>@author  Francis X. Giraldo on 7/08
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !---------------------------------------------------------------------!
    subroutine lagrange_basis_filter(ngl,xgl,nq,xnq,wnq,psi,dpsi)
        integer, intent(in)  :: ngl, nq
        real,    intent(in)  :: xgl(ngl)
        real,    intent(in) :: xnq(nq), wnq(nq)
        real,    intent(out) :: psi(ngl,nq), dpsi(ngl,nq)

        real :: ksi, xl, xi, xj, xk, ddpsi
        integer :: i, j, k, l

        !Get High Order Roots
        !>  call legendre_gauss_lobatto(nq,xnq,wnq)

        !Do Quadrature
        do l=1,nq
            xl=xnq(l)

            !Construct Bases
            do i=1,ngl
                ksi=xgl(i)
                psi(i,l)=1
                dpsi(i,l)=0

                do j=1,ngl
                    xj=xgl(j)

                    !Basis Functions
                    if (j /= i) then
                        psi(i,l)=psi(i,l)*( xl - xj )/( ksi - xj )
                    endif

                    ddpsi=1
               
                    !Derivative of Basis Functions
                    if (j /= i) then
                        do k=1,ngl
                            xk=xgl(k)
                            if (k /= i .and. k /= j) then
                                ddpsi=ddpsi*( xl - xk )/( ksi - xk )
                            endif
                        end do
                        dpsi(i,l)=dpsi(i,l) + ddpsi/( ksi - xj )
                    endif
                end do
            end do
        end do
    end subroutine lagrange_basis_filter

    !---------------------------------------------------------------------!
    !>@brief This subroutine finds the Legendre Cardinal Basis Functions and
    !>their Derivatives in in Lagrange Polynomial Form.
    !>@author  Francis X. Giraldo on 7/08
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !---------------------------------------------------------------------!
    subroutine lagrange_basis(ngl,xgl,nq,xnq,wnq,psiq,dpsiq)
        integer, intent(in)  :: ngl, nq
        real,    intent(in)  :: xgl(ngl)
        real,    intent(out) :: xnq(nq), wnq(nq)
        real,    intent(out) :: psiq(ngl,nq), dpsiq(ngl,nq)

        real :: ksi, xl, xi, xj, xk, ddpsi
        integer :: i, j, k, l

        !Get High Order Roots
        call legendre_gauss_lobatto(nq,xnq,wnq)

        !Do Quadrature
        do l=1,nq
            xl=xnq(l)


            !Construct Bases
            do i=1,ngl
                ksi=xgl(i)
                psiq(i,l)=1.0
                dpsiq(i,l)=0.0

                do j=1,ngl
                    xj=xgl(j)

                    !Basis Functions
                    if (j /= i) then
                        psiq(i,l)=psiq(i,l)*( xl - xj )/( ksi - xj )
                    endif

                    ddpsi=1.0
               
                    !Derivative of Basis Functions
                    if (j /= i) then
                        do k=1,ngl
                            xk=xgl(k)
                            if (k /= i .and. k /= j) then
                                ddpsi=ddpsi*( xl - xk )/( ksi - xk )
                            endif
                        end do
                        dpsiq(i,l)=dpsiq(i,l) + ddpsi/( ksi - xj )
                    endif
                 end do
              end do
           end do
         end subroutine lagrange_basis

    !---------------------------------------------------------------------!
    !>@brief This subroutine finds the Legendre Cardinal Basis Functions and
    !>their Derivatives in Lagrange Polynomial Form.
    !>@author  Francis X. Giraldo on 7/08
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !---------------------------------------------------------------------!
    subroutine lagrange_basis2(ngl,xgl,xd,psi,dpsi)
        integer, intent(in)  :: ngl
        real,    intent(in)  :: xd, xgl(ngl)
        real,    intent(out) :: psi(ngl), dpsi(ngl)

        real :: xl, xi, xj, xk, dpsii
        integer :: i, j, k
      
        xl=xd

        !Construct Bases
        do i=1,ngl
            xi=xgl(i)
            psi(i)=1
            dpsi(i)=0

            do j=1,ngl
                xj=xgl(j)

                !Basis Functions
                if (j /= i) then
                    psi(i)=psi(i)*( xl - xj )/( xi - xj )
                endif

                dpsii=1
               
                !Derivative of Basis Functions
                if (j /= i) then
                    do k=1,ngl
                        xk=xgl(k)
                        if (k /= i .and. k /= j) then
                            dpsii=dpsii*( xl - xk )/( xi - xk )
                        endif
                    end do
                    dpsi(i)=dpsi(i) + dpsii/( xi - xj )
                endif
            end do
        end do
   
    end subroutine lagrange_basis2

end module mod_legendre
