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

    use mod_limiter, only: flim_x, flim_y, flim_z

    use mod_metrics, only: &
        ksi_x, ksi_y, ksi_z, &
        eta_x, eta_y, eta_z, &
        zeta_x, zeta_y, zeta_z, &
        jac, massinv, massinv_1d

    use mod_viscosity, only: visc_elem, metrics_visc, q_visc, grad_q, nmetrics_visc, nvar_visc, grad_q_mlswe, q_visc_mlswe

    use mod_parallel, only : num_send_recv_total

    public :: &
        compute_gradient, &
        compute_divergence, &
        compute_vorticity, &
        compute_curl, &  !New Operator to replace COMPUTE_VORTICITY (FXG: 12/2013)
        compute_laplacian_v0, &
        compute_laplacian_IBP, &
        compute_laplacian_IBP_anisotropic, &
        compute_laplacian_IBP_horizontal, &
        compute_laplacian_IBP_vertical, &
        compute_laplacian_IBP_set2nc, &
        compute_laplacian_IBP_set2nc_anisotropic, &
        compute_laplacian_IBP_set2nc_with_gradient, &
        compute_laplacian_IBP_set2c, &
        compute_laplacian_IBP_set2c_v1, &
        compute_laplacian_IBP_set3c, &
        compute_local_gradient, &
        compute_local_gradient_v2, &
        compute_local_gradient_v3, &
        compute_local_gradient_transpose_v3, &
        compute_local_gradient_filter_v3, &
        compute_local_gradient_limiter_dss, &
        compute_strain_derivative_v0, &
        compute_strain_derivative_IBP, &
        compute_local_gradient_quad_v3, &  ! added by Yao Gahounzo
        compute_laplacian_IBP_with_gradient_mlswe ! added by Yao Gahounzo

    private
  !----------------------------------------------------------------------!


  !----------------------------------------------------------------------!

contains

    !----------------------------------------------------------------------!
    !>@brief This subroutine constructs the GRADIENT of a 1-vector Q
    !>@author  Francis X. Giraldo on 8/2006
    !>Extended to 3D by Francis X. Giraldo on 11/2009
    !>Modified to use MXM calls by F.X. Giraldo on 12/2013
    !>                  Department of Applied Mathematics
    !>                  Naval Postgraduate School
    !>                  Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_gradient(grad_q,q,imass,penalize)

        implicit none

        ! global arrays
        real, intent(out) :: grad_q(3,npoin)
        real, intent(in)  :: q(npoin)
        integer, intent(in)  :: imass
        logical, intent(in), optional :: penalize

        !local arrays
        real, dimension(npts) :: p, p_e, p_n, p_c
        integer :: inode(npts)
        real :: wq
        real, dimension(npts) :: e_x, e_y, e_z
        real, dimension(npts) :: n_x, n_y, n_z
        real, dimension(npts) :: c_x, c_y, c_z, jac_e
        real :: p_x, p_y, p_z
        integer :: i, j, k, e, ip, ii
        integer :: ndim
        logical :: pen

        real, dimension(3,num_send_recv_total) :: grad_q_nbh

        !penalize jump
        pen = .false.
        if(present(penalize)) pen = penalize
        
        !Constants for MXM
        ndim=1

        !initialize the global matrix
        grad_q=0.0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii+1
                        inode(ii)=ip
                        p(ii)  = q(ip)

                        !Metrics
                        e_x(ii)=ksi_x(i,j,k,e);  e_y(ii)=ksi_y(i,j,k,e);  e_z(ii)=ksi_z(i,j,k,e)
                        n_x(ii)=eta_x(i,j,k,e);  n_y(ii)=eta_y(i,j,k,e);  n_z(ii)=eta_z(i,j,k,e)
                        c_x(ii)=zeta_x(i,j,k,e); c_y(ii)=zeta_y(i,j,k,e); c_z(ii)=zeta_z(i,j,k,e)
                        jac_e(ii)=jac(i,j,k,e)
                    end do
                end do
            end do

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(p_e,p_n,p_c,p,nglx,ngly,nglz,ndim)

            !Loop through I (1-Vector of LGL Points)
            do i=1,npts
                ip=inode(i)

                !Gauss-Lobatto Weight and Jacobian
                wq=jac_e(i)

                !Construct Derivatives in Physical Space (via the Chain Rule)
                p_x=p_e(i)*e_x(i) + p_n(i)*n_x(i) + p_c(i)*c_x(i)
                p_y=p_e(i)*e_y(i) + p_n(i)*n_y(i) + p_c(i)*c_y(i)
                p_z=p_e(i)*e_z(i) + p_n(i)*n_z(i) + p_c(i)*c_z(i)

                !Store RHS
                grad_q(1,ip) = grad_q(1,ip) + wq*p_x
                grad_q(2,ip) = grad_q(2,ip) + wq*p_y
                grad_q(3,ip) = grad_q(3,ip) + wq*p_z
            end do !i
        end do !e

        !compute fluxes
        call compute_gradient_penalty_flux(grad_q,q,pen)
        
        !Apply DSS
        call create_global_rhs(grad_q,grad_q_nbh,3,imass)

    end subroutine compute_gradient

    !----------------------------------------------------------------------!
    !>@brief This subroutine constructs the Divergence of a 3-vector Q
    !>@author  Francis X. Giraldo on 8/2006
    !>Extended to 3D by Francis X. Giraldo on 11/2009
    !>Rewritten to use MXM routines by F.X. Giraldo on 12/2013
    !>                  Department of Applied Mathematics
    !>                  Naval Postgraduate School
    !>                  Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_divergence(div_q,q,imass,penalize)

        implicit none

        !global arrays
        real, intent(out) :: div_q(npoin)
        real, intent(in)  :: q(3,npoin)
        integer, intent(in)  :: imass
        logical, intent(in), optional :: penalize

        !local arrays
        real, dimension(3,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(1,num_send_recv_total) :: div_q_nbh
        integer :: inode(npts)
        real :: wq
        real, dimension(npts) :: e_x, e_y, e_z
        real, dimension(npts) :: n_x, n_y, n_z
        real, dimension(npts) :: c_x, c_y, c_z, jac_e
        real :: u_x, v_y, w_z
        integer :: i, j, k, e, ip, ii
        integer :: ndim
        logical pen

        !penalize jump
        pen = .false.
        if(present(penalize)) pen = penalize
        
        !Constants for MXM
        ndim=3

        !initialize the global matrix
        div_q=0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii+1
                        inode(ii)=ip
                        qq(1,ii)= q(1,ip)
                        qq(2,ii)= q(2,ip)
                        qq(3,ii)= q(3,ip)

                        !Metrics
                        e_x(ii)=ksi_x(i,j,k,e);  e_y(ii)=ksi_y(i,j,k,e);  e_z(ii)=ksi_z(i,j,k,e)
                        n_x(ii)=eta_x(i,j,k,e);  n_y(ii)=eta_y(i,j,k,e);  n_z(ii)=eta_z(i,j,k,e)
                        c_x(ii)=zeta_x(i,j,k,e); c_y(ii)=zeta_y(i,j,k,e); c_z(ii)=zeta_z(i,j,k,e)
                        jac_e(ii)=jac(i,j,k,e)
                    end do
                end do
            end do

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Loop through I (1-Vector of LGL Points)
            do i=1,npts
                ip=inode(i)

                !Gauss-Lobatto Weight and Jacobian
                wq=jac_e(i)

                !Construct Derivatives in Physical Space (via the Chain Rule)
                u_x=qq_e(1,i)*e_x(i) + qq_n(1,i)*n_x(i) + qq_c(1,i)*c_x(i)
                v_y=qq_e(2,i)*e_y(i) + qq_n(2,i)*n_y(i) + qq_c(2,i)*c_y(i)
                w_z=qq_e(3,i)*e_z(i) + qq_n(3,i)*n_z(i) + qq_c(3,i)*c_z(i)

                !Store RHS
                div_q(ip) = div_q(ip) + wq*(u_x + v_y + w_z)
            end do !i

        end do !e

        !compute fluxes
        call compute_divergence_penalty_flux(div_q,q,pen)
        
        !Apply DSS
        if(imass.ne.2) then !if we need to dss
           call create_global_rhs(div_q,div_q_nbh,1,imass)
        end if

    end subroutine compute_divergence

    !----------------------------------------------------------------------!
    !>@brief This subroutine computes the Curl of the 3-vector Q
    !>@author  Francis X. Giraldo on 9/9/2011
    !>Rewritten to use MXM routines by F.X. Giraldo on 12/2013
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_curl(curl_q,q,imass)

        implicit none

        !global arrays
        real, intent(out) :: curl_q(3,npoin)
        real, intent(in)  :: q(3,npoin)
        integer, intent(in)  :: imass

        !local arrays
        real, dimension(3,npts) :: qq, qq_e, qq_n, qq_c
        integer :: inode(npts)
        real :: wq
        real, dimension(npts) :: e_x, e_y, e_z
        real, dimension(npts) :: n_x, n_y, n_z
        real, dimension(npts) :: c_x, c_y, c_z, jac_e
        real :: u_y, u_z
        real :: v_x, v_z
        real :: w_x, w_y
        integer :: i, j, k, m, e, ip, ii
        integer :: m1, m2, m3, ndim

        real, dimension(3,num_send_recv_total) :: curl_q_nbh

        !initialize the global matrix
        curl_q=0
        ndim=3

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii+1
                        inode(ii)=ip
                        qq(1,ii)= q(1,ip)
                        qq(2,ii)= q(2,ip)
                        qq(3,ii)= q(3,ip)

                        !Metrics
                        e_x(ii)=ksi_x(i,j,k,e);  e_y(ii)=ksi_y(i,j,k,e);  e_z(ii)=ksi_z(i,j,k,e)
                        n_x(ii)=eta_x(i,j,k,e);  n_y(ii)=eta_y(i,j,k,e);  n_z(ii)=eta_z(i,j,k,e)
                        c_x(ii)=zeta_x(i,j,k,e); c_y(ii)=zeta_y(i,j,k,e); c_z(ii)=zeta_z(i,j,k,e)
                        jac_e(ii)=jac(i,j,k,e)
                    end do
                end do
            end do

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Loop through I (and L, LGL Points)
            do i=1,npts
                ip=inode(i)

                !Gauss-Lobatto Weight and Jacobian
                wq=jac_e(i)

                !Construct Derivatives in Physical Space (via the Chain Rule)
                u_y=qq_e(1,i)*e_y(i) + qq_n(1,i)*n_y(i) + qq_c(1,i)*c_y(i)
                u_z=qq_e(1,i)*e_z(i) + qq_n(1,i)*n_z(i) + qq_c(1,i)*c_z(i)
                v_x=qq_e(2,i)*e_x(i) + qq_n(2,i)*n_x(i) + qq_c(2,i)*c_x(i)
                v_z=qq_e(2,i)*e_z(i) + qq_n(2,i)*n_z(i) + qq_c(2,i)*c_z(i)
                w_x=qq_e(3,i)*e_x(i) + qq_n(3,i)*n_x(i) + qq_c(3,i)*c_x(i)
                w_y=qq_e(3,i)*e_y(i) + qq_n(3,i)*n_y(i) + qq_c(3,i)*c_y(i)

                !Store the solution
                curl_q(1,ip) = curl_q(1,ip) + wq*(w_y - v_z)
                curl_q(2,ip) = curl_q(2,ip) + wq*(u_z - w_x)
                curl_q(3,ip) = curl_q(3,ip) + wq*(v_x - u_y)
            end do !i
        end do !e

        !Apply DSS
        call create_global_rhs(curl_q,curl_q_nbh,3,imass)

    end subroutine compute_curl

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
    !>@brief This subroutine builds the Laplacian Operator of a NDIM-Vector Q
    !>Using the storage (NDIM,N,N,N)
    !>This is different from COMPUTE_LAPLACIAN_V1 in that we apply viscosity different to each dimension.
    !>@author  Francis X. Giraldo on 1/2014
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_IBP(lap_q,q,ndim)

      implicit none

      !Global Arrays
      real, intent(out) :: lap_q(ndim,npoin)
      real, intent(in)  :: q(ndim,npoin)
      integer, intent(in) :: ndim

      !Note: No Communication is done in any of these routines
      if (lincompressible) then
         call compute_laplacian_IBP_PISO(lap_q,q,ndim)
      else if (eqn_set(1:6) == 'set2nc') then
         call compute_laplacian_IBP_set2nc(lap_q,q,ndim)
      else if (eqn_set(1:5) == 'set2c') then
            !!       call compute_laplacian_IBP_set2c(lap_q,q,ndim)
         call compute_laplacian_IBP_set2nc_with_gradient(lap_q,q,ndim)
      else if (eqn_set(1:5) == 'set3c') then
         call compute_laplacian_IBP_set3c(lap_q,q,ndim)
      end if

    end subroutine compute_laplacian_IBP

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator of a NDIM-Vector Q
    !>Using the storage (NDIM,N,N,N)
    !>This is different from COMPUTE_LAPLACIAN_V1 in that we apply viscosity different to each dimension.
    !>@author  Francis X. Giraldo on 1/2014
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_IBP_anisotropic(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Note: No Communication is done in any of these routines
        if (eqn_set(1:6) == 'set2nc') then
            call compute_laplacian_IBP_set2nc_anisotropic(lap_q,q,ndim)
        else if (eqn_set(1:5) == 'set2c') then
            !!       call compute_laplacian_IBP_set2c(lap_q,q,ndim)
            call compute_laplacian_IBP_set2nc_with_gradient(lap_q,q,ndim)
        else if (eqn_set(1:5) == 'set3c') then
            call compute_laplacian_IBP_set3c(lap_q,q,ndim)
        end if

    end subroutine compute_laplacian_IBP_anisotropic

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator of a NDIM-Vector Q
    !>Using the storage (NDIM,N,N,N)
    !>This is different from COMPUTE_LAPLACIAN_IBP_set2NC in that we separate out the
    !>Horizontal and Vertical directions and apply the viscosity differently.
    !>@author  Francis X. Giraldo on 10/31/2016 - Happy Halloween!!
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_IBP_set2nc(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(npts)
        integer :: e, i, j, k, l, m, ip, ii
        real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim,npts) :: qq_1, qq_2, qq_3
        real, dimension(6,npts) :: g
        real :: gg(6)

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii + 1
                        inode(ii)=ip
                        do m=1,ndim
                            qq(m,ii)=q(m,ip)
                        end do
                        do m=1,6
                            g(m,ii)=metrics_visc(m,i,j,k,e)
                        end do
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Add Metric Terms (they are already pre-multiplied by wgl^3*Jac
            do i=1,npts
                gg(:)=g(:,i)
                do m=1,ndim
                    qq_1(m,i)=gg(1)*qq_e(m,i) + gg(2)*qq_n(m,i) + gg(3)*qq_c(m,i)
                    qq_2(m,i)=gg(2)*qq_e(m,i) + gg(4)*qq_n(m,i) + gg(5)*qq_c(m,i)
                    qq_3(m,i)=gg(3)*qq_e(m,i) + gg(5)*qq_n(m,i) + gg(6)*qq_c(m,i)
                end do !m
            end do !i

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_transpose_v3(qq,qq_1,qq_2,qq_3,nglx,ngly,nglz,ndim)

            !Store Values
            do i=1,npts
                ip=inode(i)
                do m=2,ndim
                    lap_q(m,ip)=lap_q(m,ip) - qq(m,i)
                end do !m
            end do !i

        end do !e

    end subroutine compute_laplacian_IBP_set2nc

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator of a NDIM-Vector Q
    !>Using the storage (NDIM,N,N,N)
    !>This is different from COMPUTE_LAPLACIAN_IBP_set2NC in that we separate out the
    !>Horizontal and Vertical directions and apply the viscosity differently.
    !>@author  Francis X. Giraldo on 10/31/2016 - Happy Halloween!!
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_IBP_set2nc_anisotropic(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(npts)
        integer :: e, i, j, k, l, m, ip, ii
        real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim,npts) :: qq_1, qq_2, qq_3
        real, dimension(nmetrics_visc,npts) :: g
        real :: gg(nmetrics_visc)

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii + 1
                        inode(ii)=ip
                        do m=1,ndim
                            qq(m,ii)=q(m,ip)
                        end do
                        do m=1,nmetrics_visc
                            g(m,ii)=metrics_visc(m,i,j,k,e)
                        end do
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Add Metric Terms (they are already pre-multiplied by wgl^3*Jac
            do i=1,npts
                gg(:)=g(:,i)
                do m=1,nvar
                    qq_1(m,i)=gg(1)*qq_e(m,i) + gg(2)*qq_n(m,i) + gg(3)*qq_c(m,i)
                    qq_2(m,i)=gg(2)*qq_e(m,i) + gg(4)*qq_n(m,i) + gg(5)*qq_c(m,i)
                    qq_3(m,i)=gg(3)*qq_e(m,i) + gg(5)*qq_n(m,i) + gg(6)*qq_c(m,i)
                end do !m
                do m=nvar+1,nvar_visc
                    qq_1(m,i)=gg(7)*qq_e(m,i) + gg(8)*qq_n(m,i)  + gg(9)*qq_c(m,i)
                    qq_2(m,i)=gg(8)*qq_e(m,i) + gg(10)*qq_n(m,i) + gg(11)*qq_c(m,i)
                    qq_3(m,i)=gg(9)*qq_e(m,i) + gg(11)*qq_n(m,i) + gg(12)*qq_c(m,i)
                end do !m
            end do !i

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_transpose_v3(qq,qq_1,qq_2,qq_3,nglx,ngly,nglz,ndim)

            !Store Values
            do i=1,npts
                ip=inode(i)
                do m=1,ndim
                    lap_q(m,ip)=lap_q(m,ip) - qq(m,i)
                end do !m
            end do !i

        end do !e

    end subroutine compute_laplacian_IBP_set2nc_anisotropic

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator of a NDIM-Vector Q.
    !> Note that it is not mathematically equivalent to Set2NC but
    !> simplified in order to generalize to hyper-viscosity.
    !> It is in fact identical to IBP_SET2NC except we store GRAD_Q
    !>@author  Francis X. Giraldo on 8/14/2016
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_IBP_set2nc_with_gradient(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(npts)
        integer :: e, i, j, k, l, m, ip, ii
        real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim,npts) :: qq_1, qq_2, qq_3
        real, dimension(ndim)      ::  qq_x, qq_y, qq_z
        real, dimension(npts) :: rho
        real, dimension(6,npts) :: g
        real :: g1, g2, g3, g4, g5, g6

        !Metric Variables and Constants
        real :: wq
        real :: e_x, e_y, e_z
        real :: n_x, n_y, n_z
        real :: c_x, c_y, c_z
        real :: h_e, h_n, h_c
        real, dimension(npts) :: ksi_x_e, ksi_y_e, ksi_z_e, &
            eta_x_e, eta_y_e, eta_z_e, &
            zeta_x_e,zeta_y_e,zeta_z_e, jac_e

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii + 1
                        inode(ii)=ip
                        do m=1,ndim
                            qq(m,ii)=q(m,ip)
                        end do
                        do m=1,6
                            g(m,ii)=metrics_visc(m,i,j,k,e)
                        end do

                        !Metrics
                        ksi_x_e(ii)=ksi_x(i,j,k,e);   ksi_y_e(ii)=ksi_y(i,j,k,e);   ksi_z_e(ii)=ksi_z(i,j,k,e)
                        eta_x_e(ii)=eta_x(i,j,k,e);   eta_y_e(ii)=eta_y(i,j,k,e);   eta_z_e(ii)=eta_z(i,j,k,e)
                        zeta_x_e(ii)=zeta_x(i,j,k,e); zeta_y_e(ii)=zeta_y(i,j,k,e); zeta_z_e(ii)=zeta_z(i,j,k,e)
                        jac_e(ii)=jac(i,j,k,e)
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Add Metric Terms
            do i=1,npts
                g1=g(1,i); g2=g(2,i); g3=g(3,i); g4=g(4,i); g5=g(5,i); g6=g(6,i)
                do m=1,ndim
                    qq_1(m,i)=g1*qq_e(m,i) + g2*qq_n(m,i) + g3*qq_c(m,i)
                    qq_2(m,i)=g2*qq_e(m,i) + g4*qq_n(m,i) + g5*qq_c(m,i)
                    qq_3(m,i)=g3*qq_e(m,i) + g5*qq_n(m,i) + g6*qq_c(m,i)
                end do !m

                !Store Metric Terms
                e_x=ksi_x_e(i);  e_y=ksi_y_e(i);  e_z=ksi_z_e(i)
                n_x=eta_x_e(i);  n_y=eta_y_e(i);  n_z=eta_z_e(i)
                c_x=zeta_x_e(i); c_y=zeta_y_e(i); c_z=zeta_z_e(i)

                !Construct Derivatives in Physical Space (via the Chain Rule)
                !----Perturbation Variables------!
                do m=1,ndim
                    qq_x(m)=qq_e(m,i)*e_x + qq_n(m,i)*n_x + qq_c(m,i)*c_x
                    qq_y(m)=qq_e(m,i)*e_y + qq_n(m,i)*n_y + qq_c(m,i)*c_y
                    qq_z(m)=qq_e(m,i)*e_z + qq_n(m,i)*n_z + qq_c(m,i)*c_z
                end do !m

                !Store Derivatives and Pass them viad MOD_VISCOSITY
                do m=1,ndim
                    grad_q(1,m,ip)=qq_x(m)
                    grad_q(2,m,ip)=qq_y(m)
                    grad_q(3,m,ip)=qq_z(m)
                end do

            end do !i

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_transpose_v3(qq,qq_1,qq_2,qq_3,nglx,ngly,nglz,ndim)

            !Store Values
            do i=1,npts
                ip=inode(i)
                do m=1,ndim
                    lap_q(m,ip)=lap_q(m,ip) - qq(m,i)
                end do !m
            end do !i
        end do !e

    end subroutine compute_laplacian_IBP_set2nc_with_gradient

    subroutine compute_laplacian_IBP_with_gradient_mlswe(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(npts)
        integer :: e, i, j, k, l, m, ip, ii
        real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim,npts) :: qq_1, qq_2, qq_3
        real, dimension(ndim)      ::  qq_x, qq_y, qq_z
        real, dimension(npts) :: rho
        real, dimension(6,npts) :: g
        real :: g1, g2, g3, g4, g5, g6

        !Metric Variables and Constants
        real :: wq
        real :: e_x, e_y, e_z
        real :: n_x, n_y, n_z
        real :: c_x, c_y, c_z
        real :: h_e, h_n, h_c
        real, dimension(npts) :: ksi_x_e, ksi_y_e, ksi_z_e, &
            eta_x_e, eta_y_e, eta_z_e, &
            zeta_x_e,zeta_y_e,zeta_z_e, jac_e

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii + 1
                        inode(ii)=ip
                        do m=1,ndim
                            qq(m,ii)=q(m,ip)
                        end do
                        do m=1,6
                            g(m,ii)=metrics_visc(m,i,j,k,e)
                        end do

                        !Metrics
                        ksi_x_e(ii)=ksi_x(i,j,k,e);   ksi_y_e(ii)=ksi_y(i,j,k,e);   ksi_z_e(ii)=ksi_z(i,j,k,e)
                        eta_x_e(ii)=eta_x(i,j,k,e);   eta_y_e(ii)=eta_y(i,j,k,e);   eta_z_e(ii)=eta_z(i,j,k,e)
                        zeta_x_e(ii)=zeta_x(i,j,k,e); zeta_y_e(ii)=zeta_y(i,j,k,e); zeta_z_e(ii)=zeta_z(i,j,k,e)
                        jac_e(ii)=jac(i,j,k,e)
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Add Metric Terms
            do i=1,npts
                g1=g(1,i); g2=g(2,i); g3=g(3,i); g4=g(4,i); g5=g(5,i); g6=g(6,i)
                do m=1,ndim
                    qq_1(m,i)=g1*qq_e(m,i) + g2*qq_n(m,i) + g3*qq_c(m,i)
                    qq_2(m,i)=g2*qq_e(m,i) + g4*qq_n(m,i) + g5*qq_c(m,i)
                    qq_3(m,i)=g3*qq_e(m,i) + g5*qq_n(m,i) + g6*qq_c(m,i)
                end do !m

                !Store Metric Terms
                e_x=ksi_x_e(i);  e_y=ksi_y_e(i);  e_z=ksi_z_e(i)
                n_x=eta_x_e(i);  n_y=eta_y_e(i);  n_z=eta_z_e(i)
                c_x=zeta_x_e(i); c_y=zeta_y_e(i); c_z=zeta_z_e(i)

                !Construct Derivatives in Physical Space (via the Chain Rule)
                !----Perturbation Variables------!
                do m=1,ndim
                    qq_x(m)=qq_e(m,i)*e_x + qq_n(m,i)*n_x + qq_c(m,i)*c_x
                    qq_y(m)=qq_e(m,i)*e_y + qq_n(m,i)*n_y + qq_c(m,i)*c_y
                    qq_z(m)=qq_e(m,i)*e_z + qq_n(m,i)*n_z + qq_c(m,i)*c_z
                end do !m

                !Store Derivatives and Pass them viad MOD_VISCOSITY
                do m=1,ndim
                    grad_q_mlswe(1,m,ip)=qq_x(m)
                    grad_q_mlswe(2,m,ip)=qq_y(m)
                    grad_q_mlswe(3,m,ip)=qq_z(m)
                end do

            end do !i

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_transpose_v3(qq,qq_1,qq_2,qq_3,nglx,ngly,nglz,ndim)

            !Store Values
            do i=1,npts
                ip=inode(i)
                do m=1,ndim
                    lap_q(m,ip)=lap_q(m,ip) - qq(m,i)
                end do !m
            end do !i
        end do !e

    end subroutine compute_laplacian_IBP_with_gradient_mlswe

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator for the Poisson 
    !>for Pressure. Based on COMPUTE_LAPLACIAN_IBP_SET2NC_with_Gradient
    !>@author  Francis X. Giraldo on 12/17/2016
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_IBP_PISO(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(npts)
        integer :: e, i, j, k, l, m, ip, ii
        real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim,npts) :: qq_1, qq_2, qq_3
!        real, dimension(ndim)      ::  qq_x, qq_y, qq_z
!        real, dimension(npts) :: rho
        real, dimension(6,npts) :: g
        real :: g1, g2, g3, g4, g5, g6

!         !Metric Variables and Constants
!         real :: wq
!         real :: e_x, e_y, e_z
!         real :: n_x, n_y, n_z
!         real :: c_x, c_y, c_z
!         real :: h_e, h_n, h_c
!         real, dimension(npts) :: ksi_x_e, ksi_y_e, ksi_z_e, &
!             eta_x_e, eta_y_e, eta_z_e, &
!             zeta_x_e,zeta_y_e,zeta_z_e, jac_e

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii + 1
                        inode(ii)=ip
                        do m=1,ndim
                            qq(m,ii)=q(m,ip)
                        end do
                        do m=1,6
                            g(m,ii)=metrics_visc(m,i,j,k,e)
                        end do

!                         !Metrics
!                         ksi_x_e(ii)=ksi_x(i,j,k,e);   ksi_y_e(ii)=ksi_y(i,j,k,e);   ksi_z_e(ii)=ksi_z(i,j,k,e)
!                         eta_x_e(ii)=eta_x(i,j,k,e);   eta_y_e(ii)=eta_y(i,j,k,e);   eta_z_e(ii)=eta_z(i,j,k,e)
!                         zeta_x_e(ii)=zeta_x(i,j,k,e); zeta_y_e(ii)=zeta_y(i,j,k,e); zeta_z_e(ii)=zeta_z(i,j,k,e)
!                         jac_e(ii)=jac(i,j,k,e)
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Add Metric Terms
            do i=1,npts
                g1=g(1,i); g2=g(2,i); g3=g(3,i); g4=g(4,i); g5=g(5,i); g6=g(6,i)
                do m=1,ndim
                    qq_1(m,i)=g1*qq_e(m,i) + g2*qq_n(m,i) + g3*qq_c(m,i)
                    qq_2(m,i)=g2*qq_e(m,i) + g4*qq_n(m,i) + g5*qq_c(m,i)
                    qq_3(m,i)=g3*qq_e(m,i) + g5*qq_n(m,i) + g6*qq_c(m,i)
                end do !m

!                 !Store Metric Terms
!                 e_x=ksi_x_e(i);  e_y=ksi_y_e(i);  e_z=ksi_z_e(i)
!                 n_x=eta_x_e(i);  n_y=eta_y_e(i);  n_z=eta_z_e(i)
!                 c_x=zeta_x_e(i); c_y=zeta_y_e(i); c_z=zeta_z_e(i)

!                 !Construct Derivatives in Physical Space (via the Chain Rule)
!                 !----Perturbation Variables------!
!                 do m=1,ndim
!                     qq_x(m)=qq_e(m,i)*e_x + qq_n(m,i)*n_x + qq_c(m,i)*c_x
!                     qq_y(m)=qq_e(m,i)*e_y + qq_n(m,i)*n_y + qq_c(m,i)*c_y
!                     qq_z(m)=qq_e(m,i)*e_z + qq_n(m,i)*n_z + qq_c(m,i)*c_z
!                 end do !m

!                 !Store Derivatives and Pass them via MOD_VISCOSITY
!                 do m=1,ndim
!                     grad_q(1,m,ip)=qq_x(m)
!                     grad_q(2,m,ip)=qq_y(m)
!                     grad_q(3,m,ip)=qq_z(m)
!                 end do
            end do !i

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_transpose_v3(qq,qq_1,qq_2,qq_3,nglx,ngly,nglz,ndim)

            !Store Values
            do i=1,npts
                ip=inode(i)
                do m=1,ndim
                    lap_q(m,ip)=lap_q(m,ip) - qq(m,i)
                end do !m
            end do !i
        end do !e

      end subroutine compute_laplacian_IBP_PISO

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator 
    !>Using the storage (NDIM,N,N,N) in horizontal direction only.
    !>Based on compute_laplacian_IBP_set2nc_anisotropic by FXG.
    !>@author  Michal A. Kopera on 07/24/2017 
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_IBP_horizontal(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(npts)
        integer :: e, i, j, k, l, m, ip, ii
        real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim,npts) :: qq_1, qq_2, qq_3
        real, dimension(nmetrics_visc,npts) :: g
        real :: gg(nmetrics_visc)

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii + 1
                        inode(ii)=ip
                        do m=1,ndim
                            qq(m,ii)=q(m,ip)
                        end do
                        do m=1,nmetrics_visc
                            g(m,ii)=metrics_visc(m,i,j,k,e)
                        end do
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Add Metric Terms (they are already pre-multiplied by wgl^3*Jac
            do i=1,npts
                gg(:)=g(:,i)
                !horizontal
                do m=1,ndim
                    qq_1(m,i)=gg(1)*qq_e(m,i) + gg(2)*qq_n(m,i) + gg(3)*qq_c(m,i)
                    qq_2(m,i)=gg(2)*qq_e(m,i) + gg(4)*qq_n(m,i) + gg(5)*qq_c(m,i)
                    qq_3(m,i)=gg(3)*qq_e(m,i) + gg(5)*qq_n(m,i) + gg(6)*qq_c(m,i)
                end do !m

                !vertical
!                 do m=1,ndim
!                     qq_1(m,i)=gg(7)*qq_e(m,i) + gg(8)*qq_n(m,i)  + gg(9)*qq_c(m,i)
!                     qq_2(m,i)=gg(8)*qq_e(m,i) + gg(10)*qq_n(m,i) + gg(11)*qq_c(m,i)
!                     qq_3(m,i)=gg(9)*qq_e(m,i) + gg(11)*qq_n(m,i) + gg(12)*qq_c(m,i)
!                 end do !m
            end do !i

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_transpose_v3(qq,qq_1,qq_2,qq_3,nglx,ngly,nglz,ndim)

            !Store Values
            do i=1,npts
                ip=inode(i)
                do m=1,ndim
                    lap_q(m,ip)=lap_q(m,ip) - qq(m,i)
                end do !m
            end do !i

        end do !e

      end subroutine compute_laplacian_IBP_horizontal

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator of a single variable
    !>Using the storage (N,N,N) in a vertical direction only
    !>Based on compute_laplacian_IBP_set2nc_anisotropic by FXG.
    !>@author  Michal A. Kopera on 07/24/2017 
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_IBP_vertical(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(npts)
        integer :: e, i, j, k, l, m, ip, ii
        real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim,npts) :: qq_1, qq_2, qq_3
        real, dimension(nmetrics_visc,npts) :: g
        real :: gg(nmetrics_visc)

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii + 1
                        inode(ii)=ip
                        do m=1,ndim
                            qq(m,ii)=q(m,ip)
                        end do
                        do m=1,nmetrics_visc
                            g(m,ii)=metrics_visc(m,i,j,k,e)
                        end do
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Add Metric Terms (they are already pre-multiplied by wgl^3*Jac
            do i=1,npts
                gg(:)=g(:,i)
                !horizontal
!                 do m=1,ndim
!                     qq_1(m,i)=gg(1)*qq_e(m,i) + gg(2)*qq_n(m,i) + gg(3)*qq_c(m,i)
!                     qq_2(m,i)=gg(2)*qq_e(m,i) + gg(4)*qq_n(m,i) + gg(5)*qq_c(m,i)
!                     qq_3(m,i)=gg(3)*qq_e(m,i) + gg(5)*qq_n(m,i) + gg(6)*qq_c(m,i)
!                 end do !m

                !vertical
                do m=1,ndim                   
                    qq_1(m,i)=gg(7)*qq_e(m,i) + gg(8)*qq_n(m,i)  + gg(9)*qq_c(m,i)
                    qq_2(m,i)=gg(8)*qq_e(m,i) + gg(10)*qq_n(m,i) + gg(11)*qq_c(m,i)
                    qq_3(m,i)=gg(9)*qq_e(m,i) + gg(11)*qq_n(m,i) + gg(12)*qq_c(m,i)
                end do !m
            end do !i

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_transpose_v3(qq,qq_1,qq_2,qq_3,nglx,ngly,nglz,ndim)

            !Store Values
            do i=1,npts
                ip=inode(i)
                do m=1,ndim
                    lap_q(m,ip)=lap_q(m,ip) - qq(m,i)
                end do !m
            end do !i

        end do !e

      end subroutine compute_laplacian_IBP_vertical

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator of a NDIM-Vector Q.
    !> Note that it is not mathematically equivalent to Set2NC but
    !> simplified in order to generalize to hyper-viscosity.
    !> It is in fact identical to IBP_SET2NC except we store GRAD_Q
    !>@author  Francis X. Giraldo on 8/14/2016
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_IBP_set2c_v1(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(npts)
        integer :: e, i, j, k, l, m, ip, ii
        real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim,npts) :: qq_1, qq_2, qq_3
        real, dimension(ndim)      ::  qq_x, qq_y, qq_z
        real, dimension(npts) :: rho
        real, dimension(6,npts) :: g
        real :: g1, g2, g3, g4, g5, g6

        !Metric Variables and Constants
        real :: wq
        real :: e_x, e_y, e_z
        real :: n_x, n_y, n_z
        real :: c_x, c_y, c_z
        real :: h_e, h_n, h_c
        real, dimension(npts) :: ksi_x_e, ksi_y_e, ksi_z_e, &
            eta_x_e, eta_y_e, eta_z_e, &
            zeta_x_e,zeta_y_e,zeta_z_e, jac_e

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii + 1
                        inode(ii)=ip
                        do m=1,ndim
                            qq(m,ii)=q(m,ip)
                        end do
                        do m=1,6
                            g(m,ii)=metrics_visc(m,i,j,k,e)
                        end do

                        !Metrics
                        ksi_x_e(ii)=ksi_x(i,j,k,e);   ksi_y_e(ii)=ksi_y(i,j,k,e);   ksi_z_e(ii)=ksi_z(i,j,k,e)
                        eta_x_e(ii)=eta_x(i,j,k,e);   eta_y_e(ii)=eta_y(i,j,k,e);   eta_z_e(ii)=eta_z(i,j,k,e)
                        zeta_x_e(ii)=zeta_x(i,j,k,e); zeta_y_e(ii)=zeta_y(i,j,k,e); zeta_z_e(ii)=zeta_z(i,j,k,e)
                        jac_e(ii)=jac(i,j,k,e)
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Add Metric Terms
            do i=1,npts
                g1=g(1,i); g2=g(2,i); g3=g(3,i); g4=g(4,i); g5=g(5,i); g6=g(6,i)
                do m=1,ndim
                    qq_1(m,i)=g1*qq_e(m,i) + g2*qq_n(m,i) + g3*qq_c(m,i)
                    qq_2(m,i)=g2*qq_e(m,i) + g4*qq_n(m,i) + g5*qq_c(m,i)
                    qq_3(m,i)=g3*qq_e(m,i) + g5*qq_n(m,i) + g6*qq_c(m,i)
                end do !m

                !Store Metric Terms
                e_x=ksi_x_e(i);  e_y=ksi_y_e(i);  e_z=ksi_z_e(i)
                n_x=eta_x_e(i);  n_y=eta_y_e(i);  n_z=eta_z_e(i)
                c_x=zeta_x_e(i); c_y=zeta_y_e(i); c_z=zeta_z_e(i)

                !Construct Derivatives in Physical Space (via the Chain Rule)
                !----Perturbation Variables------!
                do m=1,ndim
                    qq_x(m)=qq_e(m,i)*e_x + qq_n(m,i)*n_x + qq_c(m,i)*c_x
                    qq_y(m)=qq_e(m,i)*e_y + qq_n(m,i)*n_y + qq_c(m,i)*c_y
                    qq_z(m)=qq_e(m,i)*e_z + qq_n(m,i)*n_z + qq_c(m,i)*c_z
                end do !m

                !Store Derivatives and Pass them viad MOD_VISCOSITY
                do m=1,ndim
                    grad_q(1,m,ip)=qq_x(m)
                    grad_q(2,m,ip)=qq_y(m)
                    grad_q(3,m,ip)=qq_z(m)
                end do

            end do !i

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_transpose_v3(qq,qq_1,qq_2,qq_3,nglx,ngly,nglz,ndim)

            !Store Values
            do i=1,npts
                ip=inode(i)
                do m=2,ndim
                    lap_q(m,ip)=lap_q(m,ip) - qq(m,i)
                end do !m
            end do !i
        end do !e

    end subroutine compute_laplacian_IBP_set2c_v1

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator of a NDIM-Vector Q
    !> for SET2C. It is mathematically equivalent to IBP_SET2NC
    !>@author  Francis X. Giraldo on 1/2014
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_IBP_set2c(lap_q,q,ndim)
        !!  subroutine compute_laplacian_IBP_set2c_v0(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(npts)
        integer :: e, i, j, k, l, m, ip, ii
        real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim,npts) :: qq_1, qq_2, qq_3
        real, dimension(ndim)      :: qq_x, qq_y, qq_z
        real, dimension(npts) :: rho
        real, dimension(6,npts) :: g
        real :: g1, g2, g3, g4, g5, g6
        real :: r_k, t_k

        real, dimension(ndim,num_send_recv_total) :: lap_q_nbh

        !Metric Variables and Constants
        real :: wq
        real :: e_x, e_y, e_z
        real :: n_x, n_y, n_z
        real :: c_x, c_y, c_z
        real :: h_e, h_n, h_c
        real, dimension(npts) :: ksi_x_e, ksi_y_e, ksi_z_e, &
            eta_x_e, eta_y_e, eta_z_e, &
            zeta_x_e,zeta_y_e,zeta_z_e, jac_e

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii + 1
                        inode(ii)=ip
                        r_k=q(1,ip) + q_ref(1,ip) !Total Density
                        t_k=( q(5,ip) + q_ref(5,ip) )/r_k - q_ref(5,ip)/q_ref(1,ip) !Perturbation Potential Temperature
                        rho(ii)=r_k
                        do m=2,ndim
                            qq(m,ii)=q(m,ip)/r_k
                        end do
                        qq(1,ii)=r_k
                        qq(5,ii)=t_k
                        do m=1,6
                            g(m,ii)=metrics_visc(m,i,j,k,e)
                        end do

                        !Metrics
                        ksi_x_e(ii)=ksi_x(i,j,k,e);   ksi_y_e(ii)=ksi_y(i,j,k,e);   ksi_z_e(ii)=ksi_z(i,j,k,e)
                        eta_x_e(ii)=eta_x(i,j,k,e);   eta_y_e(ii)=eta_y(i,j,k,e);   eta_z_e(ii)=eta_z(i,j,k,e)
                        zeta_x_e(ii)=zeta_x(i,j,k,e); zeta_y_e(ii)=zeta_y(i,j,k,e); zeta_z_e(ii)=zeta_z(i,j,k,e)
                        jac_e(ii)=jac(i,j,k,e)
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Add Metric Terms
            do i=1,npts
                g1=g(1,i); g2=g(2,i); g3=g(3,i); g4=g(4,i); g5=g(5,i); g6=g(6,i)
                do m=1,ndim
                    qq_1(m,i)=g1*qq_e(m,i) + g2*qq_n(m,i) + g3*qq_c(m,i)
                    qq_2(m,i)=g2*qq_e(m,i) + g4*qq_n(m,i) + g5*qq_c(m,i)
                    qq_3(m,i)=g3*qq_e(m,i) + g5*qq_n(m,i) + g6*qq_c(m,i)
                end do !m

                !Store Metric Terms
                e_x=ksi_x_e(i);  e_y=ksi_y_e(i);  e_z=ksi_z_e(i)
                n_x=eta_x_e(i);  n_y=eta_y_e(i);  n_z=eta_z_e(i)
                c_x=zeta_x_e(i); c_y=zeta_y_e(i); c_z=zeta_z_e(i)

                !Construct Derivatives in Physical Space (via the Chain Rule)
                !----Perturbation Variables------!
                do m=1,ndim
                    qq_x(m)=qq_e(m,i)*e_x + qq_n(m,i)*n_x + qq_c(m,i)*c_x
                    qq_y(m)=qq_e(m,i)*e_y + qq_n(m,i)*n_y + qq_c(m,i)*c_y
                    qq_z(m)=qq_e(m,i)*e_z + qq_n(m,i)*n_z + qq_c(m,i)*c_z
                end do !m

                !Store Derivatives and Pass them viad MOD_VISCOSITY
                do m=1,ndim
                    grad_q(1,m,ip)=qq_x(m)
                    grad_q(2,m,ip)=qq_y(m)
                    grad_q(3,m,ip)=qq_z(m)
                end do

            end do !i

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_transpose_v3(qq,qq_1,qq_2,qq_3,nglx,ngly,nglz,ndim)

            !Store Values
            do i=1,npts
                ip=inode(i)
                do m=2,ndim
                    lap_q(m,ip)=lap_q(m,ip) - rho(i)*qq(m,i) &
                        - 0*( qq_e(1,i)*qq_1(m,i) + qq_n(1,i)*qq_2(m,i) + qq_c(1,i)*qq_3(m,i) ) !2nd term should be included to make it similar to Integral psi rho Laplacian Q but result is closer when this term is neglected.
                end do !m
            end do !i
        end do !e

      !Apply DSS

    end subroutine compute_laplacian_IBP_set2c

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator of a NDIM-Vector Q
    !>@author  Simone Marras on 2/2014
    !>from compute_laplacian_set2c(lap_q,q,ndim) above
    !>
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_IBP_set3c(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(npts)
        integer :: e, i, j, k, l, m, ip, ii
        real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim,npts) :: qq_1, qq_2, qq_3
        real, dimension(npts) :: rho
        real, dimension(6,npts) :: g
        real :: g1, g2, g3, g4, g5, g6
        real :: r_k, t_k

        !!    real, dimension(ndim,num_send_recv_total) :: lap_q_nbh

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii + 1
                        inode(ii)=ip
                        r_k=q(1,ip) + q_ref(1,ip) !rho total
                        t_k=(q(5,ip) + q_ref(5,ip))/r_k - q_ref(5,ip)/q_ref(1,ip) !Energy Perturbation
                        rho(ii)=r_k
                        do m=2,ndim
                            qq(m,ii)=q(m,ip)/r_k
                        end do
                        qq(1,ii)=r_k !rho total
                        qq(5,ii)=t_k !energy perturbation: not really the correct term here
                        do m=1,6
                            g(m,ii)=metrics_visc(m,i,j,k,e)
                        end do
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Add Metric Terms
            do i=1,npts
                g1=g(1,i); g2=g(2,i); g3=g(3,i); g4=g(4,i); g5=g(5,i); g6=g(6,i)
                do m=1,ndim
                    qq_1(m,i)=g1*qq_e(m,i) + g2*qq_n(m,i) + g3*qq_c(m,i)
                    qq_2(m,i)=g2*qq_e(m,i) + g4*qq_n(m,i) + g5*qq_c(m,i)
                    qq_3(m,i)=g3*qq_e(m,i) + g5*qq_n(m,i) + g6*qq_c(m,i)
                end do !m
            end do !i

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_transpose_v3(qq,qq_1,qq_2,qq_3,nglx,ngly,nglz,ndim)

            !Store Values
            do i=1,npts
                ip=inode(i)
                do m=2,ndim
                    lap_q(m,ip)=lap_q(m,ip) - rho(i)*qq(m,i) &
                        - 0*( qq_e(1,i)*qq_1(m,i) + qq_n(1,i)*qq_2(m,i) + qq_c(1,i)*qq_3(m,i) ) !2nd term should be included to make it similar to Integral psi rho Laplacian Q but result is closer when this term is neglected.
                end do !m
            end do !i
        end do !e

      !Apply DSS

    end subroutine compute_laplacian_IBP_set3c

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator of a NDIM-Vector Q
    !>Using the storage (NDIM,N,N,N)
    !>@author  Francis X. Giraldo on 1/2014
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_v1(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(npts)
        integer :: e, i, j, k, l, m, ip, ii
        real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim,npts) :: qq_1, qq_2, qq_3
        real, dimension(6,npts) :: g
        real :: g1, g2, g3, g4, g5, g6

        real, dimension(ndim,num_send_recv_total) :: lap_q_nbh

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            ii=0
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        ii=ii + 1
                        inode(ii)=ip
                        do m=1,ndim
                            qq(m,ii)=q(m,ip)
                        end do
                        do m=1,6
                            g(m,ii)=metrics_visc(m,i,j,k,e)
                        end do
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Add Metric Terms (they are already pre-multiplied by wgl^3*Jac
            do i=1,npts
                g1=g(1,i); g2=g(2,i); g3=g(3,i); g4=g(4,i); g5=g(5,i); g6=g(6,i)
                do m=1,ndim
                    qq_1(m,i)=g1*qq_e(m,i) + g2*qq_n(m,i) + g3*qq_c(m,i)
                    qq_2(m,i)=g2*qq_e(m,i) + g4*qq_n(m,i) + g5*qq_c(m,i)
                    qq_3(m,i)=g3*qq_e(m,i) + g5*qq_n(m,i) + g6*qq_c(m,i)
                end do !m
            end do !i

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_transpose_v3(qq,qq_1,qq_2,qq_3,nglx,ngly,nglz,ndim)

            !Store Values
            do i=1,npts
                ip=inode(i)
                do m=1,ndim
                    lap_q(m,ip)=lap_q(m,ip) - qq(m,i)
                end do !m
            end do !i
        end do !e

        !Apply DSS
        call create_global_rhs(lap_q,lap_q_nbh,ndim,1)

    end subroutine compute_laplacian_v1

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator of a NDIM-Vector Q
    !>Using the storage (N,N,N,NDIM)
    !>@author  Francis X. Giraldo on 11/2009
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_v0(lap_q,q,ndim)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(ndim,npoin)
        real, intent(in)  :: q(ndim,npoin)
        integer, intent(in) :: ndim

        !Local Arrays
        integer :: inode(nglx,ngly,nglz)
        integer :: e, i, j, k, l, m, ip, ipe, ipn, ipc
        !Perturbation Variables
        real, dimension(nglx,ngly,nglz,ndim) :: qq, qq_e, qq_n, qq_c
        real, dimension(ndim) :: qq_x, qq_y, qq_z
        !RHS and Metric Variables
        real :: wq
        real :: e_x, e_y, e_z
        real :: n_x, n_y, n_z
        real :: c_x, c_y, c_z
        real :: h_e, h_n, h_c
        real :: dhqde, dhqdn, dhqdc
        integer :: m1, m2, m3

        real, dimension(ndim,num_send_recv_total) :: lap_q_nbh

        !initialize
        lap_q = 0

        !loop thru the elements
        do e=1,nelem

            !Store Element Variables
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,e)
                        inode(i,j,k)=ip
                        do m=1,ndim
                            qq(i,j,k,m)=q(m,ip)
                        end do
                    end do !i
                end do !j
            end do !k

            !Construct Local Derivatives for Prognostics Variables
            call compute_local_gradient_v2(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

            !Loop through I (and L, LGL Points)
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=inode(i,j,k)

                        !Gauss-Lobatto Weight and Jacobian
                        wq=jac(i,j,k,e)

                        !Store Metric Terms
                        e_x=ksi_x(i,j,k,e);  e_y=ksi_y(i,j,k,e);  e_z=ksi_z(i,j,k,e)
                        n_x=eta_x(i,j,k,e);  n_y=eta_y(i,j,k,e);  n_z=eta_z(i,j,k,e)
                        c_x=zeta_x(i,j,k,e); c_y=zeta_y(i,j,k,e); c_z=zeta_z(i,j,k,e)

                        !Construct Derivatives in Physical Space (via the Chain Rule)
                        !----Perturbation Variables------!
                        do m=1,ndim
                            qq_x(m)=qq_e(i,j,k,m)*e_x + qq_n(i,j,k,m)*n_x + qq_c(i,j,k,m)*c_x
                            qq_y(m)=qq_e(i,j,k,m)*e_y + qq_n(i,j,k,m)*n_y + qq_c(i,j,k,m)*c_y
                            qq_z(m)=qq_e(i,j,k,m)*e_z + qq_n(i,j,k,m)*n_z + qq_c(i,j,k,m)*c_z
                        end do !m

                        !Diffusion Operator
                        do l=1,ngl
                            h_e = dpsi(l,i)
                            ipe=inode(l,j,k)
                            h_n = dpsi(l,j)
                            ipn=inode(i,l,k)
                            h_c = dpsi(l,k)
                            ipc=inode(i,j,l)
                            do m=1,ndim
                                !Derivatives along KSI
                                dhqde=h_e*( e_x*qq_x(m) + e_y*qq_y(m) + e_z*qq_z(m) )
                                lap_q(m,ipe)=lap_q(m,ipe) - wq*dhqde
                                !Derivatives along ETA
                                dhqdn=h_n*( n_x*qq_x(m) + n_y*qq_y(m) + n_z*qq_z(m) )
                                lap_q(m,ipn)=lap_q(m,ipn) - wq*dhqdn
                                !Derivatives along ZETA
                                dhqdc=h_c*( c_x*qq_x(m) + c_y*qq_y(m) + c_z*qq_z(m) )
                                lap_q(m,ipc)=lap_q(m,ipc) - wq*dhqdc
                            end do !m
                        end do !l

                    end do !i
                end do !j
            end do !k
        end do !e

        !Apply DSS
        call create_global_rhs(lap_q,lap_q_nbh,ndim,1)

    end subroutine compute_laplacian_v0

    !----------------------------------------------------------------------!
    !>@brief This subroutine builds the Laplacian Operator for a 1-Vector Q
    !>@author  Francis X. Giraldo on 11/2009
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_laplacian_1vector(lap_q,q)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(npoin)
        real, intent(in)  :: q(npoin)

        !Local Arrays
        integer inode(ngl,ngl,ngl)
        integer ie, i, j, k, l, m, ip, ipe, ipn, ipc
        !Perturbation Variables
        real qq(ngl,ngl,ngl), bq
        real qq_e, qq_n, qq_c
        real qq_x, qq_y, qq_z
        !RHS and Metric Variables
        real wq
        real e_x, e_y, e_z
        real n_x, n_y, n_z
        real c_x, c_y, c_z
        real h_e, h_n, h_c
        real dhqde, dhqdn, dhqdc

        real, dimension(1,num_send_recv_total) :: lap_q_nbh

        !initialize
        lap_q = 0

        !loop thru the elements
        do ie=1,nelem

            !Store Element Variables
            do k=1,ngl
                do j=1,ngl
                    do i=1,ngl
                        ip=intma(i,j,k,ie)
                        inode(i,j,k)=ip
                        qq(i,j,k)=q(ip)
                    end do !i
                end do !j
            end do !k

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
                        qq_e=0; qq_n=0; qq_c=0

                        do l = 1,ngl
                            ! Derivatives of Basis functions
                            h_e = dpsi(l,i)
                            h_n = dpsi(l,j)
                            h_c = dpsi(l,k)

                            !KSI, ETA, ZETA Derivatives
                            qq_e = qq_e + qq(l,j,k)*h_e
                            qq_n = qq_n + qq(i,l,k)*h_n
                            qq_c = qq_c + qq(i,j,l)*h_c
                        end do !l

                        !Construct Derivatives in Physical Space (via the Chain Rule)
                        !----Perturbation Variables------!
                        qq_x=qq_e*e_x + qq_n*n_x + qq_c*c_x
                        qq_y=qq_e*e_y + qq_n*n_y + qq_c*c_y
                        qq_z=qq_e*e_z + qq_n*n_z + qq_c*c_z

                        !Diffusion Operator
                        do l=1,ngl
                            h_e = dpsi(l,i)
                            ipe=inode(l,j,k)
                            h_n = dpsi(l,j)
                            ipn=inode(i,l,k)
                            h_c = dpsi(l,k)
                            ipc=inode(i,j,l)
                            !Derivatives along KSI
                            dhqde=h_e*( e_x*qq_x + e_y*qq_y + e_z*qq_z )
                            lap_q(ipe)=lap_q(ipe) + wq*dhqde
                            !Derivatives along ETA
                            dhqdn=h_n*( n_x*qq_x + n_y*qq_y + n_z*qq_z )
                            lap_q(ipn)=lap_q(ipn) + wq*dhqdn
                            !Derivatives along ZETA
                            dhqdc=h_c*( c_x*qq_x + c_y*qq_y + c_z*qq_z )
                            lap_q(ipc)=lap_q(ipc) + wq*dhqdc
                        end do !l

                    end do !i
                end do !j
            end do !k
        end do !ie

        !Apply DSS
        call create_global_rhs(lap_q,lap_q_nbh,1,1)

    end subroutine compute_laplacian_1vector


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

    ! subroutine compute_local_gradient_quad(q_e,q_n,q_c,q,nglx,ngly,nglz,ndim)

    !     use mod_basis, only: psiqx, psiqy, dpsiqx, dpsiqy, dpsiqz, nqx, nqy, nqz


    !     implicit none

    !     !global arrays
    !     real, dimension(ndim,nqx,nqy,nqz), intent(out) :: q_e, q_n, q_c
    !     real, dimension(ndim,nqx,nqy,nqz), intent(in)  :: q
    !     integer, intent(in) :: nglx, ngly, nglz, ndim

    !     !local arrays
    !     integer :: m1, m2, m3, m, k
    !     real, dimension(nqx,nqy,nqz) :: qt, qt_e, qt_n, qt_c

    !     !Construct Local Derivatives
    !     do m=1,ndim
    !         !Store Input
    !         qt=q(m,:,:,:)

    !         !KSI Derivative
    !         if(nqx > 1) then
    !             m1=nqx !nglx
    !             m2=nglx !nglx
    !             m3=nqy*nqz !nsize*ngly*nglz
    !             call mxm(dpsiqx_tr,m1,qt,m2,qt_e,m3)
    !         else
    !             qt_e = 0
    !         endif

    !         !ETA Derivative
    !         if(nqy > 1) then
    !             m1=nqx !nglx
    !             m2=ngly !ngly
    !             m3=nqy !ngly
    !             do k=1,nglz
    !                 call mxm(qt(1,1,k),m1,dpsiqy,m2,qt_n(1,1,k),m3)
    !             enddo !k
    !         else
    !             qt_n = 0
    !         endif

    !         !Zeta Derivative
    !         if(nqz > 1) then
    !             m1=nglx*ngly !nglx*nsize*ngly
    !             m2=nglz !nglz
    !             m3=nglz !nglz
    !             call mxm(qt,m1,dpsiqz,m2,qt_c,m3)
    !         else
    !             qt_c = 0
    !         endif

    !         !Store Output
    !         q_e(m,:,:,:)=qt_e
    !         q_n(m,:,:,:)=qt_n
    !         q_c(m,:,:,:)=qt_c
    !     end do !m

    ! end subroutine compute_local_gradient_quad

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

    !----------------------------------------------------------------------!
    !>@brief For CG, after limiting via an element-wise limiter, we need to enforce 
    !>C^0 continuity via a DSS operator. However, we need to perform integration before 
    !>applying DSS which will divide by the mass matrix.
    !>@author  Francis X. Giraldo on 2/2017
    !>                  Department of Applied Mathematics
    !>                  Naval Postgraduate School
    !>                  Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_local_gradient_limiter_dss(fqf,q,nglx,ngly,nglz,ndim)

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
                call mxm(flim_x,m1,qt,m2,qt_i,m3)
            else
                qt_i = qt
            endif

            !ETA Derivative
            if(ngly > 1) then
                m1=nglx !nglx
                m2=ngly !ngly
                m3=ngly !ngly
                do k=1,nglz
                    call mxm(qt_i(1,1,k),m1,flim_y,m2,qt_ij(1,1,k),m3)
                enddo !k
            else
                qt_ij = qt_i
            endif

            !Zeta Derivative
            if(nglz > 1) then
                m1=nglx*ngly !nglx*nsize*ngly
                m2=nglz !nglz
                m3=nglz !nglz
                call mxm(qt_ij,m1,flim_z,m2,qt_ijk,m3)
            else
                qt_ijk = qt_ij
            endif

            !Store Output
            fqf(m,:,:,:)=qt_ijk
        end do !m

   end subroutine compute_local_gradient_limiter_dss
   
   
   !----------------------------------------------------------------------!
   !>@brief This subroutine computes the inner product of Grad Psi .dot. Grad Q
   !>for the strain tensor
   !>
   !> This is done WITHOUT integration by parts. 
   !>
   !>@author  Simone Marras on 11/2018
   !>                  Department of Mechanical Engineering
   !>                  New Jersey Institute of Technology
   !>
   !----------------------------------------------------------------------!
   subroutine compute_strain_derivative_v0(dtau_dxij,q,ndim)

     implicit none

     !Global Arrays
     real,    intent(out) :: dtau_dxij(ndim,npoin)
     real,    intent(in)  :: q(ndim,npoin)
     integer, intent(in)  :: ndim

     !Local Arrays
     integer :: inode(nglx,ngly,nglz)
     integer :: e, i, j, k, l, m, ip, ipe, ipn, ipc
     !Perturbation Variables
     real, dimension(nglx,ngly,nglz,ndim) :: dqdx, dqdy, dqdz
     real, dimension(nglx,ngly,nglz,ndim) :: dqdx_e, dqdx_n, dqdx_c
     real, dimension(nglx,ngly,nglz,ndim) :: dqdy_e, dqdy_n, dqdy_c
     real, dimension(nglx,ngly,nglz,ndim) :: dqdz_e, dqdz_n, dqdz_c
    
     real, dimension(ndim) :: qq_xx, qq_yx, qq_zx
     real, dimension(ndim) :: qq_xy, qq_yy, qq_zy
     real, dimension(ndim) :: qq_xz, qq_yz, qq_zz
     
     !RHS and Metric Variables
     real :: wq
     real :: e_x, e_y, e_z
     real :: n_x, n_y, n_z
     real :: c_x, c_y, c_z
     real :: h_e, h_n, h_c
     real :: dhqde, dhqdn, dhqdc
     integer :: m1, m2, m3

     !real, dimension(ndim,num_send_recv_total) :: lap_q_nbh

     !initialize
     dtau_dxij = 0.0
     
     !loop thru the elements
     do e=1,nelem

        !Store Element Variables
        do k=1,nglz
           do j=1,ngly
              do i=1,nglx
                 ip=intma(i,j,k,e)
                 inode(i,j,k)=ip
                 do m=1,ndim
                    dqdx(i,j,k,m)=grad_q(1,m,ip)
                    dqdy(i,j,k,m)=grad_q(2,m,ip)
                    dqdz(i,j,k,m)=grad_q(3,m,ip)
                 end do
              end do !i
           end do !j
        end do !k

        !Construct Local Derivatives for Prognostics Variables
        call compute_local_gradient_v2(dqdx_e,dqdx_n,dqdx_c,dqdx,nglx,ngly,nglz,ndim)
        call compute_local_gradient_v2(dqdy_e,dqdy_n,dqdy_c,dqdy,nglx,ngly,nglz,ndim)
        call compute_local_gradient_v2(dqdz_e,dqdz_n,dqdz_c,dqdz,nglx,ngly,nglz,ndim)

        !Loop through I (and L, LGL Points)
        do k=1,nglz
           do j=1,ngly
              do i=1,nglx
                 ip=inode(i,j,k)

                 !Gauss-Lobatto Weight and Jacobian
                 wq=jac(i,j,k,e)

                 !Store Metric Terms
                 e_x=ksi_x(i,j,k,e);  e_y=ksi_y(i,j,k,e);  e_z=ksi_z(i,j,k,e)
                 n_x=eta_x(i,j,k,e);  n_y=eta_y(i,j,k,e);  n_z=eta_z(i,j,k,e)
                 c_x=zeta_x(i,j,k,e); c_y=zeta_y(i,j,k,e); c_z=zeta_z(i,j,k,e)

                 !Construct Derivatives in Physical Space (via the Chain Rule)
                 !----Perturbation Variables------!
                 
                 qq_xx(1)=dqdx_e(i,j,k,1)*e_x + dqdx_n(i,j,k,1)*n_x + dqdx_c(i,j,k,1)*c_x !d^2u/dx^2
                 qq_xx(2)=dqdx_e(i,j,k,2)*e_x + dqdx_n(i,j,k,2)*n_x + dqdx_c(i,j,k,2)*c_x !d^2v/dx^2
                 qq_xx(2)=dqdx_e(i,j,k,3)*e_x + dqdx_n(i,j,k,3)*n_x + dqdx_c(i,j,k,3)*c_x !d^2w/dx^2
                 
                 qq_yy(1)=dqdx_e(i,j,k,1)*e_x + dqdx_n(i,j,k,1)*n_x + dqdx_c(i,j,k,1)*c_x !d^2u/dy^2
                 qq_yy(2)=dqdx_e(i,j,k,2)*e_x + dqdx_n(i,j,k,2)*n_x + dqdx_c(i,j,k,2)*c_x !d^2v/dy^2
                 qq_yy(2)=dqdx_e(i,j,k,3)*e_x + dqdx_n(i,j,k,3)*n_x + dqdx_c(i,j,k,3)*c_x !d^2w/dy^2
                 
                 qq_zz(1)=dqdz_e(i,j,k,1)*e_z + dqdz_n(i,j,k,1)*n_z + dqdz_c(i,j,k,1)*c_z !d^2u/dz^2
                 qq_zz(2)=dqdz_e(i,j,k,2)*e_z + dqdz_n(i,j,k,2)*n_z + dqdz_c(i,j,k,2)*c_z !d^2v/dz^2
                 qq_zz(2)=dqdz_e(i,j,k,3)*e_z + dqdz_n(i,j,k,3)*n_z + dqdz_c(i,j,k,3)*c_z !d^2w/dz^2

                 !
                 qq_xy(1)=dqdx_e(i,j,k,1)*e_y + dqdx_n(i,j,k,1)*n_y + dqdx_c(i,j,k,1)*c_y !d^2u/dxdy
                 qq_xy(2)=dqdx_e(i,j,k,2)*e_y + dqdx_n(i,j,k,2)*n_y + dqdx_c(i,j,k,2)*c_y !d^2v/dxdy
                 qq_xy(3)=dqdx_e(i,j,k,3)*e_y + dqdx_n(i,j,k,3)*n_y + dqdx_c(i,j,k,3)*c_y !d^2w/dxdy
                 
                 qq_yx(1)=dqdy_e(i,j,k,1)*e_x + dqdy_n(i,j,k,1)*n_x + dqdy_c(i,j,k,1)*c_x !d^2u/dxdy
                 qq_yx(2)=dqdy_e(i,j,k,2)*e_x + dqdy_n(i,j,k,2)*n_x + dqdy_c(i,j,k,2)*c_x !d^2v/dxdy
                 qq_yx(3)=dqdy_e(i,j,k,3)*e_x + dqdy_n(i,j,k,3)*n_x + dqdy_c(i,j,k,3)*c_x !d^2w/dxdy

                 !
                 qq_xz(1)=dqdx_e(i,j,k,1)*e_z + dqdx_n(i,j,k,1)*n_z + dqdx_c(i,j,k,1)*c_z !d^2u/dxdz
                 qq_xz(2)=dqdx_e(i,j,k,2)*e_z + dqdx_n(i,j,k,2)*n_z + dqdx_c(i,j,k,2)*c_z !d^2v/dxdz
                 qq_xz(3)=dqdx_e(i,j,k,3)*e_z + dqdx_n(i,j,k,3)*n_z + dqdx_c(i,j,k,3)*c_z !d^2w/dxdz
                  
                 qq_zx(1)=dqdz_e(i,j,k,1)*e_x + dqdz_n(i,j,k,1)*n_x + dqdz_c(i,j,k,1)*c_x !d^2u/dxdz
                 qq_zx(2)=dqdz_e(i,j,k,2)*e_x + dqdz_n(i,j,k,2)*n_x + dqdz_c(i,j,k,2)*c_x !d^2v/dxdz
                 qq_zx(3)=dqdz_e(i,j,k,3)*e_x + dqdz_n(i,j,k,3)*n_x + dqdz_c(i,j,k,3)*c_x !d^2w/dxdz

                 !
                 qq_yz(1)=dqdy_e(i,j,k,1)*e_z + dqdy_n(i,j,k,1)*n_z + dqdy_c(i,j,k,1)*c_z !d^2u/dzdx
                 qq_yz(2)=dqdy_e(i,j,k,2)*e_z + dqdy_n(i,j,k,2)*n_z + dqdy_c(i,j,k,2)*c_z !d^2v/dzdx
                 qq_yz(3)=dqdy_e(i,j,k,3)*e_z + dqdy_n(i,j,k,3)*n_z + dqdy_c(i,j,k,3)*c_z !d^2w/dzdx
                 
                 qq_zy(1)=dqdz_e(i,j,k,1)*e_y + dqdz_n(i,j,k,1)*n_y + dqdz_c(i,j,k,1)*c_y !d^2u/dzdx
                 qq_zy(2)=dqdz_e(i,j,k,2)*e_y + dqdz_n(i,j,k,2)*n_y + dqdz_c(i,j,k,2)*c_y !d^2v/dzdx
                 qq_zy(3)=dqdz_e(i,j,k,3)*e_y + dqdz_n(i,j,k,3)*n_y + dqdz_c(i,j,k,3)*c_y !d^2w/dzdx

                 dtau_dxij(1,npoin) = qq_xx(1) + qq_xy(2) + qq_xz(3)
                 dtau_dxij(2,npoin) = qq_yx(1) + qq_yy(2) + qq_yz(3)
                 dtau_dxij(3,npoin) = qq_zx(1) + qq_zy(2) + qq_zz(3)
                 
              end do !i
           end do !j
        end do !k

     end do !e

     !Apply DSS
     !call create_global_rhs(lap_q,lap_q_nbh,ndim,1)

   end subroutine compute_strain_derivative_v0


 end module mod_interface
