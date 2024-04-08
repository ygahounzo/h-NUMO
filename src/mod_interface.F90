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

    use mod_initial, only: nvar, q_ref

    use mod_input, only: nelz, eqn_set, lvisc_anisotropic, lincompressible, &
         space_method, cgdg_method, nlaplacian_flg

    use mod_metrics, only: &
        ksi_x, ksi_y, ksi_z, &
        eta_x, eta_y, eta_z, &
        zeta_x, zeta_y, zeta_z, &
        jac, massinv, massinv_1d

    use mod_viscosity, only: visc_elem, metrics_visc, q_visc, grad_q, nmetrics_visc, nvar_visc, grad_q_mlswe, q_visc_mlswe

    use mod_parallel, only : num_send_recv_total

    public :: &
        compute_vorticity, &
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
    !>@brief This subroutine computes the Vorticity in the Radial Direction of a 3-vector Q
    !>@author  Francis X. Giraldo on 9/9/2011
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_vorticity(vorticity_q,q,imass)

        implicit none

        !global arrays
        real, intent(out) :: vorticity_q(3,npoin)
        real, intent(in)  :: q(3,npoin)
        integer, intent(in)  :: imass

        !local arrays
        real :: u(ngl,ngl,ngl)
        real :: v(ngl,ngl,ngl)
        real :: w(ngl,ngl,ngl)
        integer :: inode(ngl,ngl,ngl)
        real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z
        real :: u_e, u_n, u_c
        real :: v_e, v_n, v_c
        real :: w_e, w_n, w_c
        real :: h_e, h_n, h_c
        real :: u_x, u_y, u_z
        real :: v_x, v_y, v_z
        real :: w_x, w_y, w_z
        integer :: i, j, k, l, ie, ip
        integer :: AllocateStatus

        real, dimension(3,num_send_recv_total) :: vorticity_q_nbh

        !initialize the global matrix
        vorticity_q=0

        !loop thru the elements
        do ie=1,nelem

            !Store Element Variables
            do k=1,ngl
                do j=1,ngl
                    do i=1,ngl
                        ip=intma(i,j,k,ie)
                        inode(i,j,k)=ip
                        ip=inode(i,j,k)
                        u(i,j,k)  = q(1,ip)
                        v(i,j,k)  = q(2,ip)
                        w(i,j,k)  = q(3,ip)
                    end do
                end do
            end do

            !Loop through I (and L, LGL Points)
            do k=1,ngl
                do j=1,ngl
                    do i=1,ngl
                        ip=inode(i,j,k)

                        !Gauss-Lobatto Weight and Jacobian
                        wq=jac(i,j,k,ie)

                        !Store Metric Terms
                        e_x=ksi_x(i,j,k,ie);  e_y=ksi_y(i,j,k,ie);  e_z=ksi_z(i,j,k,ie)
                        n_x=eta_x(i,j,k,ie);  n_y=eta_y(i,j,k,ie);  n_z=eta_z(i,j,k,ie)
                        c_x=zeta_x(i,j,k,ie); c_y=zeta_y(i,j,k,ie); c_z=zeta_z(i,j,k,ie)

                        !construct Derivatives in Computational Space
                        u_e=0; u_n=0; u_c=0
                        v_e=0; v_n=0; v_c=0
                        w_e=0; w_n=0; w_c=0

                        do l = 1,ngl
                            ! Derivatives of Basis functions
                            h_e = dpsi(l,i)
                            h_n = dpsi(l,j)
                            h_c = dpsi(l,k)

                            !KSI Derivatives
                            u_e = u_e + u(l,j,k)*h_e
                            v_e = v_e + v(l,j,k)*h_e
                            w_e = w_e + w(l,j,k)*h_e

                            !ETA Derivatives
                            u_n = u_n + u(i,l,k)*h_n
                            v_n = v_n + v(i,l,k)*h_n
                            w_n = w_n + w(i,l,k)*h_n

                            !ZETA Derivatives
                            u_c = u_c + u(i,j,l)*h_c
                            v_c = v_c + v(i,j,l)*h_c
                            w_c = w_c + w(i,j,l)*h_c
                        end do !l

                        !Construct Derivatives in Physical Space (via the Chain Rule)
                        u_x=u_e*e_x + u_n*n_x + u_c*c_x
                        u_y=u_e*e_y + u_n*n_y + u_c*c_y
                        u_z=u_e*e_z + u_n*n_z + u_c*c_z

                        v_x=v_e*e_x + v_n*n_x + v_c*c_x
                        v_y=v_e*e_y + v_n*n_y + v_c*c_y
                        v_z=v_e*e_z + v_n*n_z + v_c*c_z

                        w_x=w_e*e_x + w_n*n_x + w_c*c_x
                        w_y=w_e*e_y + w_n*n_y + w_c*c_y
                        w_z=w_e*e_z + w_n*n_z + w_c*c_z

                        !Store the solution
                        vorticity_q(1,ip) = vorticity_q(1,ip) + wq*(w_y - v_z)
                        vorticity_q(2,ip) = vorticity_q(2,ip) + wq*(u_z - w_x)
                        vorticity_q(3,ip) = vorticity_q(3,ip) + wq*(v_x - u_y)
                    end do !i
                end do !j
            end do !k
        end do !ie

        !Apply DSS
        call create_global_rhs(vorticity_q,vorticity_q_nbh,3,imass)

    end subroutine compute_vorticity

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
