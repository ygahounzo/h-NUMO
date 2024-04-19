!----------------------------------------------------------------------!
!>@brief This module contains basic operators; derivatives (1st and 2nd),
!> curl, divergence and laplacian
!>@author  Sasa Gabersek on 11/2011
!>            Naval Research Laboratory
!>            Monterey, CA 93943-5216
!>@update Sasa put the routines into the module. The subroutines were written
!>by FX Giraldo
!----------------------------------------------------------------------!
module mod_interface

    use mod_basis, only: ngl, dpsi, &
        nglx, ngly, nglz, npts, &
        dpsix, dpsiy, dpsiz, &
        dpsix_tr, dpsiy_tr, dpsiz_tr, &
        wgl, wglx, wgly, wglz, &
        fx, fy, fz, fx_tr, fy_tr, fz_tr, &
        psiqx, psiqy, psiqz, dpsiqx, dpsiqy, dpsiqz

    use mod_grid, only: intma, intma_1d, npoin, nelem, nz

    use mod_initial, only: nvar

    use mod_input, only: nelz, eqn_set, lvisc_anisotropic, lincompressible, &
         space_method, cgdg_method, nlaplacian_flg

    use mod_metrics, only: &
        ksi_x, ksi_y, ksi_z, &
        eta_x, eta_y, eta_z, &
        zeta_x, zeta_y, zeta_z, &
        jac, massinv, massinv_1d

    use mod_parallel, only : num_send_recv_total

    public :: &
        compute_local_gradient, &
        compute_local_gradient_v2, &
        compute_local_gradient_v3, &
        compute_local_gradient_transpose_v3, &
        compute_local_gradient_filter_v3, &
        compute_local_gradient_quad_v3, &  ! added by Yao Gahounzo

    private
  !----------------------------------------------------------------------!


  !----------------------------------------------------------------------!

contains

    !----------------------------------------------------------------------!
    !>@brief This subroutine computes the canonical derivative of the vector Q as one NxNDIMxN system.
    !>That is, it builds Q_e, Q_n, and Q_c.
    !>@author  Francis X. Giraldo on 12/2013
    !>                  Department of Applied Mathematics
    !>                  Naval Postgraduate School
    !>                  Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_local_gradient(q_e,q_n,q_c,q,nglx,ngly,nglz,ndim)

        implicit none

        !global arrays
        real, intent(out) :: q_e(nglx,ndim,ngly,nglz)
        real, intent(out) :: q_n(nglx,ndim,ngly,nglz)
        real, intent(out) :: q_c(nglx,ndim,ngly,nglz)
        real, intent(in)  :: q(nglx,ndim,ngly,nglz)
        integer, intent(in) :: nglx, ngly, nglz, ndim

        !local arrays
        integer :: m1, m2, m3, k

        if(nglx > 1) then
            m1=nglx !nglx
            m2=nglx !nglx
            m3=ndim*ngly*nglz !nsize*ngly*nglz
            call mxm(dpsix_tr,m1,q,m2,q_e,m3)
        else
            q_e = 0
        endif

        if(ngly > 1) then
            m1=nglx*ndim !nglx*nsize
            m2=ngly !ngly
            m3=ngly !ngly
            do k=1,nglz !nglz
                call mxm(q(1,1,1,k),m1,dpsiy,m2,q_n(1,1,1,k),m3)
            enddo !k
        else
            q_n = 0
        endif

        if(nglz > 1) then
            m1=nglx*ndim*ngly !nglx*nsize*ngly
            m2=nglz !nglz
            m3=nglz !nglz
            call mxm(q,m1,dpsiz,m2,q_c,m3)
        else
            q_c = 0
        endif

    end subroutine compute_local_gradient

    !----------------------------------------------------------------------!
    !>@brief This subroutine computes the canonical derivative of the vector Q as (N,N,N,NDIM).
    !>That is, it builds Q_e, Q_n, and Q_c.
    !>@author  Francis X. Giraldo on 12/2013
    !>                  Department of Applied Mathematics
    !>                  Naval Postgraduate School
    !>                  Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_local_gradient_v2(q_e,q_n,q_c,q,nglx,ngly,nglz,ndim)

        implicit none

        !global arrays
        real, dimension(nglx,ngly,nglz,ndim), intent(out) :: q_e, q_n, q_c
        real, dimension(nglx,ngly,nglz,ndim), intent(in)  :: q
        integer, intent(in) :: nglx, ngly, nglz, ndim

        !local arrays
        integer :: m1, m2, m3, m, k
        real, dimension(nglx,ngly,nglz) :: qt, qt_e, qt_n, qt_c

        !Construct Local Derivatives
        do m=1,ndim
            !Store Input
            qt=q(:,:,:,m)

            !KSI Derivative
            if(nglx > 1) then
                m1=nglx !nglx
                m2=nglx !nglx
                m3=ngly*nglz !nsize*ngly*nglz
                call mxm(dpsix_tr,m1,qt,m2,qt_e,m3)
            else
                qt_e = 0
            endif

            !ETA Derivative
            if(ngly > 1) then
                m1=nglx !nglx
                m2=ngly !ngly
                m3=ngly !ngly
                do k=1,ngl
                    call mxm(qt(1,1,k),m1,dpsiy,m2,qt_n(1,1,k),m3)
                enddo !k
            else
                qt_n = 0
            endif

            !Zeta Derivative
            if(nglz > 1) then
                m1=nglx*ngly !nglx*nsize*ngly
                m2=nglz !nglz
                m3=nglz !nglz
                call mxm(qt,m1,dpsiz,m2,qt_c,m3)
            else
                qt_c = 0
            endif

            !Store Output
            q_e(:,:,:,m)=qt_e
            q_n(:,:,:,m)=qt_n
            q_c(:,:,:,m)=qt_c
        end do !m

    end subroutine compute_local_gradient_v2

    !----------------------------------------------------------------------!
    !>@brief This subroutine computes the canonical derivative of the vector Q as (NDIM,N,N,N)
    !>That is, it builds Q_e, Q_n, and Q_c.
    !>@author  Francis X. Giraldo on 12/2013
    !>                  Department of Applied Mathematics
    !>                  Naval Postgraduate School
    !>                  Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_local_gradient_v3(q_e,q_n,q_c,q,nglx,ngly,nglz,ndim)

        implicit none

        !global arrays
        real, dimension(ndim,nglx,ngly,nglz), intent(out) :: q_e, q_n, q_c
        real, dimension(ndim,nglx,ngly,nglz), intent(in)  :: q
        integer, intent(in) :: nglx, ngly, nglz, ndim

        !local arrays
        integer :: m1, m2, m3, m, k
        real, dimension(nglx,ngly,nglz) :: qt, qt_e, qt_n, qt_c

        !Construct Local Derivatives
        do m=1,ndim
            !Store Input
            qt=q(m,:,:,:)

            !KSI Derivative
            if(nglx > 1) then
                m1=nglx !nglx
                m2=nglx !nglx
                m3=ngly*nglz !nsize*ngly*nglz
                call mxm(dpsix_tr,m1,qt,m2,qt_e,m3)
            else
                qt_e = 0
            endif

            !ETA Derivative
            if(ngly > 1) then
                m1=nglx !nglx
                m2=ngly !ngly
                m3=ngly !ngly
                do k=1,nglz
                    call mxm(qt(1,1,k),m1,dpsiy,m2,qt_n(1,1,k),m3)
                enddo !k
            else
                qt_n = 0
            endif

            !Zeta Derivative
            if(nglz > 1) then
                m1=nglx*ngly !nglx*nsize*ngly
                m2=nglz !nglz
                m3=nglz !nglz
                call mxm(qt,m1,dpsiz,m2,qt_c,m3)
            else
                qt_c = 0
            endif

            !Store Output
            q_e(m,:,:,:)=qt_e
            q_n(m,:,:,:)=qt_n
            q_c(m,:,:,:)=qt_c
        end do !m

    end subroutine compute_local_gradient_v3

    !----------------------------------------------------------------------!
    !>@brief This subroutine computes the canonical derivative of the vector Q as (Nx,Ny,Nz)
    !>That is, it builds Q_e, Q_n, and Q_c.
    !>@author  Yao Gahounzo 04/03/2023
    !>                  Computing PhD
    !>                  Boise State University
    !----------------------------------------------------------------------!
    subroutine compute_local_gradient_quad_v3(q_e,q_n,q_c,q,nglx,ngly,nglz,nqx,nqy,nqz)

        use mod_basis, only: psiqx, psiqy, dpsiqx, dpsiqy, dpsiqz

        implicit none

        !global arrays
        real, dimension(nqx,nqy,nqz), intent(out) :: q_e, q_n, q_c
        real, dimension(nglx,ngly,nglz), intent(in)  :: q
        integer, intent(in) :: nqx, nqy, nqz, nglx,ngly,nglz

        !local arrays
        integer :: m1, m2, m3, m, k, l, kquad, jquad, iquad,n
        real :: hix, hiy, temp

        q_e = 0
        q_n = 0
        q_c = 0

        !Construct Local Derivatives

        !KSI Derivative

        !open(1,file='hix.dat',status='unknown')

        do kquad = 1,nqz
            do jquad = 1,nqy
                do iquad = 1,nqx

                    !temp = 0.0
                    
                    do l = 1,nglz
                        do m = 1,ngly
                            do n = 1,nglx

                                q_e(iquad,jquad,kquad) = q_e(iquad,jquad,kquad) + dpsiqx(n,iquad)*psiqy(m,jquad)*q(n,m,l)
                                q_n(iquad,jquad,kquad) = q_n(iquad,jquad,kquad) + psiqx(n,iquad)*dpsiqy(m,jquad)*q(n,m,l)

                                if(nglz > 1) then
                                    q_c(iquad,jquad,kquad) = q_c(iquad,jquad,kquad) + psiqx(n,iquad)*psiqy(m,jquad)*dpsiz(l,kquad)*q(n,m,l)
                                endif
                                
                            enddo !n
                        enddo !m
                    enddo !l

                enddo !iquad
            enddo !jquad
        enddo !kquad

    end subroutine compute_local_gradient_quad_v3


    !----------------------------------------------------------------------!
    !>@brief This subroutine computes the inner product of Grad Psi .dot. Grad Q
    !>of the vector Q as (NDIM,N,N,N)
    !>That is, it builds the Laplacian Q_xx, Q_yy, and Q_zz assuming that:
    !>q_e=g1*q_e + g2*q_n + g3_g_c,
    !>q_n=g2*q_e + g4*q_n + g5*q_c,
    !>q_c=g3*q_e + g5*q_n + g6*q_c.
    !>@author  Francis X. Giraldo on 1/2014
    !>                  Department of Applied Mathematics
    !>                  Naval Postgraduate School
    !>                  Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_local_gradient_transpose_v3(q,q_e,q_n,q_c,nglx,ngly,nglz,ndim)

        implicit none

        !global arrays
        real, dimension(ndim,nglx,ngly,nglz), intent(out)  :: q
        real, dimension(ndim,nglx,ngly,nglz), intent(in)   :: q_e, q_n, q_c
        integer, intent(in) :: nglx, ngly, nglz, ndim

        !local arrays
        integer :: m1, m2, m3, m, k
        real, dimension(nglx,ngly,nglz)      :: qt_e, qt_n, qt_c
        real, dimension(nglx,ngly,nglz) :: q1, q2, q3

        !Construct Local Derivatives
        do m=1,ndim
            !Store Input
            qt_e=q_e(m,:,:,:)
            qt_n=q_n(m,:,:,:)
            qt_c=q_c(m,:,:,:)

            !KSI Derivative
            if(nglx > 1) then
                m1=nglx !nglx
                m2=nglx !nglx
                m3=ngly*nglz !nsize*ngly*nglz
                call mxm(dpsix,m1,qt_e,m2,q1,m3)
            else
                q1 = 0
            endif

            !ETA Derivative
            if(ngly > 1) then
                m1=nglx !nglx
                m2=ngly !ngly
                m3=ngly !ngly
                do k=1,nglz
                    call mxm(qt_n(1,1,k),m1,dpsiy_tr,m2,q2(1,1,k),m3)
                enddo !k
            else
                q2 = 0
            endif

            !Zeta Derivative
            if(nglz > 1) then
                m1=nglx*ngly !nglx*nsize*ngly
                m2=nglz !nglz
                m3=nglz !nglz
                call mxm(qt_c,m1,dpsiz_tr,m2,q3,m3)
            else
                q3 = 0
            endif

            !Store Output
            q(m,:,:,:)=q1(:,:,:) + q2(:,:,:) + q3(:,:,:)
        end do !m

    end subroutine compute_local_gradient_transpose_v3

    !----------------------------------------------------------------------!
    !>@brief This subroutine computes the canonical filter of the vector Q as (NDIM,N,N,N)
    !>That is, it builds Q_e, Q_n, and Q_c.
    !>@author  Francis X. Giraldo on 12/2013
    !>                  Department of Applied Mathematics
    !>                  Naval Postgraduate School
    !>                  Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_local_gradient_filter_v3(fqf,q,nglx,ngly,nglz,ndim)

        implicit none

        !global arrays
        real, dimension(ndim,nglx,ngly,nglz), intent(out) :: fqf
        real, dimension(ndim,nglx,ngly,nglz), intent(in)  :: q
        integer, intent(in) :: nglx, ngly, nglz, ndim

        !local arrays
        integer :: m1, m2, m3, m, k
        real, dimension(nglx,ngly,nglz) :: qt, qt_i, qt_ij, qt_ijk

        !Construct Local Derivatives
        do m=1,ndim
            !Store Input
            qt=q(m,:,:,:)

            !KSI Derivative
            if(nglx > 1) then
                m1=nglx !nglx
                m2=nglx !nglx
                m3=ngly*nglz !nsize*ngly*nglz
                call mxm(fx,m1,qt,m2,qt_i,m3)
            else
                qt_i = qt
            endif

            !ETA Derivative
            if(ngly > 1) then
                m1=nglx !nglx
                m2=ngly !ngly
                m3=ngly !ngly
                do k=1,nglz
                    call mxm(qt_i(1,1,k),m1,fy_tr,m2,qt_ij(1,1,k),m3)
                enddo !k
            else
                qt_ij = qt_i
            endif

            !Zeta Derivative
            if(nglz > 1) then
                m1=nglx*ngly !nglx*nsize*ngly
                m2=nglz !nglz
                m3=nglz !nglz
                call mxm(qt_ij,m1,fz_tr,m2,qt_ijk,m3)
            else
                qt_ijk = qt_ij
            endif

            !Store Output
            fqf(m,:,:,:)=qt_ijk
        end do !m

    end subroutine compute_local_gradient_filter_v3

 end module mod_interface
