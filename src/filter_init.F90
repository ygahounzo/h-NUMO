!----------------------------------------------------------------------!
!>@brief This subroutine builds the Filter Matrices using the MODAL FILTER
!>which is the BOYD-VANDEVEN FILTER
!>@author  Francis X. Giraldo on 7/08
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@date FXG to F90 on May 7, 2014
!----------------------------------------------------------------------!
subroutine filter_init(f,wgl,xgl,ngl,nfl,xmu,filter_weight_type, &
    filter_basis_type)

    use mod_legendre, only : legendre_poly, legendre_gauss_lobatto, &
        lagrange_basis_filter, legendre_poly_type

    use mod_types, only : r16

    implicit none

    !global variables
    real, intent(out) :: f(ngl,ngl)
    real, intent(in)  :: wgl(ngl), xgl(ngl)
    integer, intent(in) :: ngl, nfl
    real, intent(in) :: xmu
    character, intent(in) :: filter_basis_type*8, filter_weight_type*4

    !local variables
    real(r16) :: leg(ngl,ngl), leg2(ngl,ngl), leg_inv(ngl,ngl)
    real :: weight(ngl)
    real :: vandeven_modal
    real :: gam(ngl)
    real :: xfl(nfl), wfl(ngl)
    real :: psi_l(nfl,ngl), dpsi_l(nfl,ngl)
    real :: psi_h(ngl,nfl), dpsi_h(ngl,nfl)
    type(legendre_poly_type) :: lpoly
    integer :: nop
    real :: exp_alpha, exp_order, quad_order, quad_alpha, erf_alpha, erf_order
    real :: xi, xmode2, amp, sum
    integer :: i, j, k, jj, mode_filter, k0, ierr

    !Initialize
    f=0
  
    !Constants
    nop=ngl-1
    exp_alpha=17
    exp_order=18
    quad_alpha=1.0
    quad_order=ngl/3 !1/3 Truncation
    erf_alpha=0.0
    erf_order=12

    !Compute the Legendre Polynomial Matrix
    leg=0
    do i=1,ngl
        xi=xgl(i)
        do j=1,ngl
            jj=j-1
            call legendre_poly(jj,xi,lpoly)
            leg(i,j)=real(lpoly%p0)
        end do !j
    end do !i

    !Construct the Hierarchical Modal Legendre Basis (Szabo)
    leg2=leg
    if (filter_basis_type == 'modal') then
        do i=1,ngl
            xi=xgl(i)
            leg2(i,1)=0.5*(1 - xi)
            if(ngl > 1) then
                leg2(i,2)=0.5*(1 + xi)
                do j=3,ngl
                    leg2(i,j)=leg(i,j) - leg(i,j-2)
                end do              !j
            endif
        end do                 !i
    else if (filter_basis_type == 'legendre') then
        leg2=leg
    end if

    !Compute Inverse Matrix
    leg_inv=leg2
    call gaujordf_r16(leg_inv,ngl,ierr)
    if (ierr /= 0) then
        print*,' Error in GAUJORDF in FILTER_INIT '
        print*,' ierr = ',ierr
        stop
    end if

    !Compute the Boyd-Vandeven (ERF-LOG) Transfer Function
    if (filter_weight_type == 'erf') then
        do k=1,ngl
            weight(k)=vandeven_modal(k,ngl,erf_order)  !Boyd Filter
           !            weight(k)=vandeven_modal(k-1,ngl,erf_order) !Modified Boyd Filter
        end do                 !k

       !Compute Paul Fischer's Transfer Function
    else if (filter_weight_type == 'quad') then
        mode_filter=int(quad_order)
        k0 = ngl-mode_filter
        xmode2=mode_filter*mode_filter
        weight=1
        do k=k0+1,ngl
            amp = quad_alpha*(k-k0)*(k-k0)/(xmode2) ! quadratic growth
            weight(k) = 1.- amp
        end do

       !Compute Tim Warburton's Exponential Smoothing filter
    else if (filter_weight_type == 'exp') then
        do k=1,ngl
            weight(k)=exp( -exp_alpha*( real(k-1)/nop )**exp_order )
        end do                 !i
    end if !filter_basis_type

    !Construct 1D Filter Matrix
    do i=1,ngl
        do j=1,ngl
            sum=0
            do k=1,ngl
                sum=sum + leg2(i,k)*weight(k)*leg_inv(k,j)
            end do              !k
            f(i,j)=xmu*sum
        end do                 !j
        f(i,i)=f(i,i) + (1.0-xmu)
    end do                    !i

    if (filter_basis_type == 'nodal') then
        !Construct N-1 - > N Basis
        call legendre_gauss_lobatto(nfl,xfl,wfl)
        call lagrange_basis_filter(nfl,xfl,ngl,xgl,wgl,psi_l,dpsi_l)

        !Construct N -> N-1 Basis
        call  lagrange_basis_filter(ngl,xgl,nfl,xfl,wfl,psi_h,dpsi_h)

        !Construct Filter matrix
        f=0
        do i=1,ngl
            do j=1,ngl
                sum=0
                do k=1,nfl
                    sum=sum + psi_l(k,i)*psi_h(j,k)
                end do           !k
                f(i,j)=xmu*sum
            end do              !i
            f(i,i)=f(i,i) + (1.0-xmu)
        end do                 !j
    end if


end subroutine filter_init
!----------------------------------------------------------------------!
!>@brief This is the Boyd-Vandeven filter
!----------------------------------------------------------------------!
real function vandeven_modal(kk,ngl,p)

    implicit none

    !global variables
    integer, intent(in) :: kk, ngl
    real, intent(in) :: p

    !local variables
    real :: x, omega, z
    integer :: k, n, i
    real :: pe, a1, a2, a3, a4, a5, eps
    real :: xlog, c, diff, square_root, zc, t

    !Constants - ERF
    pe=0.3275911
    a1=0.254829592
    a2=-0.284496736
    a3=1.421413741
    a4=-1.453152027
    a5=1.061405429

    !Constants - Vandeven
    n=ngl-1
    k=kk-1
    i=2*n/3
    eps=1.0e-10

    if (k <= i) then
        x=0
        vandeven_modal=1.0
    else if (k > i .and. k < n) then
        x=real(k-i)/real(n-i)
        omega=abs(x) - 0.5
        xlog=log(1.0-4.0*omega**2)
        c=4.0*omega**2
        diff=abs(x-0.5)
        if (diff < eps) then
            square_root=1
        else
            square_root=sqrt(-xlog/c)
        end if

        z=2.0*sqrt(p)*omega*square_root
        zc=abs(z)

        !ERF
        t=1.0/(1.0 + pe*zc)
        c=1.0 - (a1*t + a2*t**2 + a3*t**3 + a4*t**4 + a5*t**5)*exp(-zc*zc)
        if (zc < eps) then
            c=0
        else
            c=c*z/zc
        end if
        vandeven_modal=0.5*(1.0-c)
    else if (k == n) then
        x=1
        vandeven_modal=0
    end if

end function vandeven_modal

!-----------------------------------------------------------------------
!>@brief     Gauss-Jordan matrix inversion with full pivoting
!>     Num. Rec. p. 30, 2nd Ed., Fortran
!>     a     is an n x n matrix
!>     rmult is a  work array of dimension n
!-----------------------------------------------------------------------
subroutine gaujordf(a,n,ierr)

    implicit none

    !global variables
    real, intent(inout) :: a(n,n)
    integer, intent(in) :: n
    integer, intent(out) :: ierr

    !local variables
    integer :: indr(n), indc(n), ipiv(n)
    integer :: i, j, k, l, ll, irow, icol
    real :: eps, big, dum, piv

    !Initialize
    ierr = 0
    eps = 1.e-9
    ipiv=0

    do i=1,n
        big=0

        !Pivot search
        do j=1,n
            if (ipiv(j) /= 1) then
                do k=1,n
                    if (ipiv(k) == 0) then
                        if (abs(a(j,k)) >= big) then
                            big=abs(a(j,k))
                            irow=j
                            icol=k
                        endif
                    elseif (ipiv(k) > 1) then
                        ierr=-ipiv(k)
                        return
                    endif
                enddo            !k
            endif
        enddo                  !j
        ipiv(icol)=ipiv(icol) + 1

        !Swap rows
        if (irow /= icol) then
            do l=1,n
                dum=a(irow,l)
                a(irow,l)=a(icol,l)
                a(icol,l)=dum
            enddo               !l
        endif
        indr(i)=irow
        indc(i)=icol
        if (abs(a(icol,icol)) < eps) then
            write(6 ,*) 'small Gauss Jordan Piv:',icol,a(icol,icol)
            ierr = icol
            return
        endif
        piv = 1./a(icol,icol)
        a(icol,icol)=1.
        do l=1,n
            a(icol,l)=a(icol,l)*piv
        enddo                  !l

        do ll=1,n
            if (ll /= icol) then
                dum=a(ll,icol)
                a(ll,icol)=0
                do l=1,n
                    a(ll,l)=a(ll,l)-a(icol,l)*dum
                end do           !l
            end if
        end do                 !ll
    end do !i

    !Unscramble Matrix
    do l=n,1,-1
        if (indr(l) /= indc(l)) then
            do k=1,n
                dum=a(k,indr(l))
                a(k,indr(l))=a(k,indc(l))
                a(k,indc(l))=dum
            end do !k
        end if
    end do !l

end subroutine gaujordf


!-----------------------------------------------------------------------
!>@brief     Gauss-Jordan matrix inversion with full pivoting for Quad Precision
!>     Num. Rec. p. 30, 2nd Ed., Fortran
!>     a     is an n x n matrix
!>     rmult is a  work array of dimension n
!-----------------------------------------------------------------------
subroutine gaujordf_r16(a,n,ierr)

    use mod_types, only: r16

    implicit none

    !global variables
    real(r16), intent(inout) :: a(n,n)
    integer, intent(in) :: n
    integer, intent(out) :: ierr

    !local variables
    integer :: indr(n), indc(n), ipiv(n)
    integer :: i, j, k, l, ll, irow, icol
    real(r16) :: eps, big, dum, piv

    !Initialize
    ierr = 0
    eps = 1.e-9
    ipiv=0

    do i=1,n
        big=0

        !Pivot search
        do j=1,n
            if (ipiv(j) /= 1) then
                do k=1,n
                    if (ipiv(k) == 0) then
                        if (abs(a(j,k)) >= big) then
                            big=abs(a(j,k))
                            irow=j
                            icol=k
                        endif
                    elseif (ipiv(k) > 1) then
                        ierr=-ipiv(k)
                        return
                    endif
                enddo            !k
            endif
        enddo                  !j
        ipiv(icol)=ipiv(icol) + 1

        !Swap rows
        if (irow /= icol) then
            do l=1,n
                dum=a(irow,l)
                a(irow,l)=a(icol,l)
                a(icol,l)=dum
            enddo               !l
        endif
        indr(i)=irow
        indc(i)=icol
        if (abs(a(icol,icol)) < eps) then
            write(6 ,*) 'small Gauss Jordan Piv:',icol,a(icol,icol)
            ierr = icol
            return
        endif
        piv = 1./a(icol,icol)
        a(icol,icol)=1.
        do l=1,n
            a(icol,l)=a(icol,l)*piv
        enddo                  !l

        do ll=1,n
            if (ll /= icol) then
                dum=a(ll,icol)
                a(ll,icol)=0
                do l=1,n
                    a(ll,l)=a(ll,l)-a(icol,l)*dum
                end do           !l
            end if
        end do                 !ll
    end do !i

    !Unscramble Matrix
    do l=n,1,-1
        if (indr(l) /= indc(l)) then
            do k=1,n
                dum=a(k,indr(l))
                a(k,indr(l))=a(k,indc(l))
                a(k,indc(l))=dum
            end do !k
        end if
    end do !l

end subroutine gaujordf_r16
