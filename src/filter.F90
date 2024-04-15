!----------------------------------------------------------------------!
!>@brief This subroutine applies the Filter
!>@author  Francis X. Giraldo on 8/99
!>@date 11/2009 Extended to 3D
!----------------------------------------------------------------------!
subroutine filter(q)

    use mod_basis, only: nglx, ngly, nglz, npts, ngl, f, fx, fy, fz

    use mod_filter, only: b, b_data

    use mod_grid, only: intma, npoin, nelem

    use mod_initial, only: nvar, nvart

    use mod_input, only: filter_tracers_flg, space_method, equations

    use mod_interface, only: compute_local_gradient_filter_v3

    use mod_metrics, only: jac, massinv
  
    implicit none
  
    !global arrays
    real, intent(inout) :: q(nvar,npoin)
  
    !local arrays
    real, dimension(nvar,npts) :: qq, fqf, qq_i, qq_ij, qq_ijk
    real, dimension(npts) :: jac_e
    real :: wq
    integer :: inode(npts)
    integer :: i, j, k, l, m, e, ip, ip1, ii, ie
    integer :: ndim, ndim2
    real r_k, u_k, v_k, w_k, t_k
  
    !Store dimensions of MxM object
    ndim  = nvar
    ndim2 = nvart + (nvar-nvart)*filter_tracers_flg
     
    !Initialize
    b=0

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
                    do m=1,nvar
                        qq(m,ii) = q(m,ip)
                    end do
                    jac_e(ii)=jac(i,j,k,e)
                end do !i
            end do !j
        end do !k
     
        !Construct Local Derivatives for Prognostics Variables
        call compute_local_gradient_filter_v3(fqf,qq,nglx,ngly,nglz,ndim)

        !Do Numerical Integration
        do i=1,npts
            ip=inode(i)

            !Gauss-Lobatto Weight and Jacobian
            wq=jac_e(i)
              
            !Store RHS
            do m=1,ndim
                b(m,ip)=b(m,ip) + wq*fqf(m,i)
            end do !m
        
        end do !i
  
    end do !e
  
    !Global DSS
    if (space_method /= 'dg') then
        call create_global_rhs(b,b_data,nvar,1) !Periodicity is applied inside
    else
        do m=1,ndim
            b(m,:) = b(m,:)*massinv(:)
        end do
    end if

    !q=b
    q(1:ndim2,:) = b(1:ndim2,:)
  
end subroutine filter


!----------------------------------------------------------------------!
!>@brief This subroutine applies the Filter without the Filter Matrix
!>in order to perform DSS with Mass matrix weighting.
!>@author  Francis X. Giraldo on 8/30/2017
!----------------------------------------------------------------------!
subroutine filter_dss(q)

    use mod_basis, only: nglx, ngly, nglz, npts, ngl

    use mod_filter, only: b, b_data

    use mod_grid, only: intma, npoin, nelem

    use mod_initial, only: nvar, nvart

    use mod_input, only: filter_tracers_flg, space_method, equations

    use mod_interface, only: compute_local_gradient_filter_v3

    use mod_metrics, only: jac, massinv
  
    implicit none
  
    !global arrays
    real, intent(inout) :: q(nvar,npoin)
  
    !local arrays
    real, dimension(nvar,npts) :: qq, fqf, qq_i, qq_ij, qq_ijk
    real, dimension(npts) :: jac_e
    real :: wq
    integer :: inode(npts)
    integer :: i, j, k, l, m, e, ip, ip1, ii, ie
    integer :: ndim, ndim2
    real r_k, u_k, v_k, w_k, t_k
  
    !Store dimensions of MxM object
    ndim  = nvar
    ndim2 = nvart + (nvar-nvart)*filter_tracers_flg

    !Initialize
    b=0

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
                    do m=1,nvar
                        qq(m,ii) = q(m,ip)
                    end do
                    jac_e(ii)=jac(i,j,k,e)
                end do !i
            end do !j
        end do !k
     
        !Do Numerical Integration
        do i=1,npts
            ip=inode(i)

            !Gauss-Lobatto Weight and Jacobian
            wq=jac_e(i)
              
            !Store RHS
            do m=1,ndim
                b(m,ip)=b(m,ip) + wq*qq(m,i)
            end do !m
        
        end do !i
  
    end do !e
  
    !Global DSS
    if (space_method /= 'dg') then
        call create_global_rhs(b,b_data,nvar,1) !Periodicity is applied inside
    else
        do m=1,ndim
            b(m,:) = b(m,:)*massinv(:)
        end do
    end if
  
    !q=b
    q(1:ndim2,:) = b(1:ndim2,:)
  
end subroutine filter_dss
