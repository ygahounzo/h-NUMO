!----------------------------------------------------------------------!
!>@brief This module contains the subroutines that create 
!> what is needed for the Limiters
!>@author Francis X. Giraldo
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@date November 16, 2016.  Originally from NUMA2dCGDG_AMR
!----------------------------------------------------------------------!
module mod_limiter

  use mod_basis, only: nglx, ngly, nglz, npts, wglx, wgly, wglz

  use mod_convert_variables, only: mod_convert_variables_eqnset_to_primitive_perturbation, &
       mod_convert_variables_primitive_perturbation_to_eqnset

  use mod_grid, only: npoin, nelem, intma

  use mod_initial, only: nvar, q_ref, q_ref_layers, rho_layers

  use mod_input, only: limiter_type, space_method, limit_threshold, limiter_qoi, is_shallow, is_swe_layers
  
  public :: &
       mod_limiter_create, &
       mod_limiter_apply, &
       q_limiter, &
       flim_x, flim_y, flim_z

  private

  !-----------------------------------------------------------------------
  !module variables and parameters
  real, dimension(:,:), allocatable :: q_limiter
  real, dimension(:,:), allocatable :: flim_x, flim_y, flim_z
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine mod_limiter_create()

    implicit none

    integer :: AllocateStatus
    integer :: e, i, j, k, l, ip

    if (allocated(q_limiter)) then
       deallocate(q_limiter, flim_x, flim_y, flim_z)
    endif
    allocate( q_limiter(nvar,npoin), flim_x(nglx,nglx), flim_y(ngly,ngly), flim_z(nglz,nglz), &
         stat=AllocateStatus )
    if (AllocateStatus /= 0) stop "** Not Enough Memory - MOD_LIMITER_CREATE **"

    !Initialize
    q_limiter=0
    flim_x=0
    flim_y=0
    flim_z=0

    !Form identity matrix for FLIM to apply DSS in MOD_LIMITER_APPLY
    do i=1,nglx
       flim_x(i,i)=1.0
    end do

    do i=1,ngly
       flim_y(i,i)=1.0
    end do

    do i=1,nglz
       flim_z(i,i)=1.0
    end do

  end subroutine mod_limiter_create


  !----------------------------------------------------------------------!
  !> @brief This subroutine calls the limiter. Comes from NUMA2dCGDG_AMR
  !> @author Francis X. Giraldo on 9/2013
  !>           Department of Applied Mathematics
  !>           Naval Postgraduate School
  !>           Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine mod_limiter_apply(q)

    implicit none

    !global arrays
    real, dimension(nvar,npoin), intent(inout) :: q

    !Convert Variables to Primitive Perturbation form
    q_limiter=q

    !Apply Limiter
    if (limiter_type == "qoi") then
       call limiter_Positivity_Preserving_QOI(q_limiter)
    else if (limiter_type == "tracers") then
       call limiter_Positivity_Preserving_Tracers(q_limiter)
    else if (limiter_type == "inundation") then
       call limiter_Positivity_Preserving_Inundation(q_limiter)
    else if (limiter_type == "bound") then
!!       call mod_limiter_Bound_Preserving(q_limiter)
    else 
       print*,'----------Warning from MOD_LIMITER.F90------------'
       print*,' limiter_type = ',limiter_type,' not supported '
       print*,' No Limiter will be Applied '
    end if

  end subroutine mod_limiter_apply

  !----------------------------------------------------------------------!
  !> @brief This subroutine applies the Positivity Preserving Limiter of Zhang and Shu JCP 2010. 
  !> @details Works for Tracers AND Prognostic Variables. However, should not apply to 
  !> Density perturbation since this is sometimes negative only (e.g., RTB case)
  !> @author Francis X. Giraldo on 9/2012 using 1D Matlab Version by Haley Lane
  !>           Department of Applied Mathematics
  !>           Naval Postgraduate School
  !>           Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine limiter_Positivity_Preserving_QOI(qp)

    implicit none

    !global arrays
    real, intent(inout) :: qp(nvar,npoin)

    !local arrays
    real, dimension(nvar,npts) :: qlocal, qlim
    real, dimension(nvar) :: theta, qmean, qmin
    real :: eps, rho
    integer :: e, i, j, k, l, m, inodes(npts), ii, ip

    !Define Constants
    eps=1e-16

    do e=1,nelem
       ii=0
       qmean=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                ii=ii+1
                inodes(ii)=ip
                !Store Primitive Variables
                qlocal(limiter_qoi,ii)=qp(limiter_qoi,ip)
                !Compute Mean Value
                qmean(limiter_qoi)=qmean(limiter_qoi) + wglx(i)*wgly(j)*wglz(k)*qlocal(limiter_qoi,ii)
             end do !i
          end do !j
       end do !k
       qmean=qmean/8 !area of canonical element in 3D (4 in 2D)

       !Compute Weight
       qmin(limiter_qoi)=minval(qlocal(limiter_qoi,:))
       theta(limiter_qoi)=min( 1.0, abs((qmean(limiter_qoi)-eps)/( qmean(limiter_qoi) - qmin(limiter_qoi)  + eps )) )

       !Apply Weighted combination of Mean Value and High-Order Solution
       ii=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ii=ii+1
                qlim(limiter_qoi,ii)=qmean(limiter_qoi) + theta(limiter_qoi)*( qlocal(limiter_qoi,ii) - qmean(limiter_qoi) )
                if (qmin(limiter_qoi) == 0 .and. qmean(limiter_qoi) == 0) qlim(limiter_qoi,ii) = 0
             end do !i
          end do !j
       end do !k
       do ii=1,npts
          ip=inodes(ii)
          qp(limiter_qoi,ip)=qlim(limiter_qoi,ii)
       end do !ii

    end do !e

  end subroutine limiter_Positivity_Preserving_QOI

  !----------------------------------------------------------------------!
  !> @brief This subroutine applies the Positivity Preserving Limiter of Zhang and Shu JCP 2010. 
  !> @details Works for Tracers AND Prognostic Variables. However, should not apply to 
  !> Density perturbation since this is sometimes negative only (e.g., RTB case)
  !> @author Francis X. Giraldo on 9/2012 using 1D Matlab Version by Haley Lane
  !>           Department of Applied Mathematics
  !>           Naval Postgraduate School
  !>           Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine limiter_Positivity_Preserving_Tracers(qp)

    implicit none

    !global arrays
    real, intent(inout) :: qp(nvar,npoin)

    !local arrays
    real, dimension(nvar,npts) :: qlocal, qlim
    real, dimension(nvar) :: theta, qmean, qmin
    real :: eps, rho
    integer :: e, i, j, k, l, m, inodes(npts), ii, ip

    !Define Constants
    eps=1e-16

    do e=1,nelem

       !Store local DOFs and compute Mean
       ii=0
       qmean=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                ii=ii+1
                inodes(ii)=ip
                !Store Primitive Variables
                do m=6,nvar
                   qlocal(m,ii)=qp(m,ip)
                end do !m
                !Compute Mean Value
                do m=6,nvar
                   qmean(m)=qmean(m) + wglx(i)*wgly(j)*wglz(k)*qlocal(m,ii)
                end do !m
             end do !i
          end do !j
       end do !k
       qmean=qmean/8 !area of canonical element in 3D (4 in 2D)

       !Compute Weight
       do m=6,nvar
          qmin(m)=minval(qlocal(m,:))
          theta(m)=min( 1.0, abs((qmean(m)-eps)/( qmean(m) - qmin(m)  + eps )) )
       end do !m

       !Apply Weighted combination of Mean Value and HO Solution
       ii=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ii=ii+1
                do m=6,nvar
                   qlim(m,ii)=qmean(m) + theta(m)*( qlocal(m,ii) - qmean(m) )
                   if (qmin(m) == 0 .and. qmean(m) == 0) qlim(m,ii) = 0
                end do !m
             end do !i
          end do !j
       end do !k
       do ii=1,npts
          ip=inodes(ii)
          do m=6,nvar
             qp(m,ip)=qlim(m,ii)
          end do !m
       end do !ii

    end do !e

  end subroutine limiter_Positivity_Preserving_Tracers


  !----------------------------------------------------------------------!
  !> @brief This subroutine applies the Positivity Preserving Limiter of Zhang and Shu JCP 2010. 
  !> @details Works for Tracers AND Prognostic Variables. However, should not apply to 
  !> Density perturbation since this is sometimes negative only (e.g., RTB case)
  !> @author Francis X. Giraldo on 9/2012 using 1D Matlab Version by Haley Lane
  !>           Department of Applied Mathematics
  !>           Naval Postgraduate School
  !>           Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine limiter_Positivity_Preserving_Inundation(qp)

    implicit none

    !global arrays
    real, intent(inout) :: qp(nvar,npoin)

    !local arrays
    real, dimension(nvar,npts) :: qlocal, qlim
    real, dimension(nvar) :: theta, qmean, qmin
    real :: eps, rho
    integer :: e, i, j, k, l, m, inodes(npts), ii, ip

    !Define Constants
    eps=1e-16

    !loop through elements
    do e=1,nelem

       !Store local DOFs and compute Mean
       ii=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                ii=ii+1
                inodes(ii)=ip
                !Store Primitive Variables
                qlocal(1,ii)=qp(1,ip) + q_ref(1,ip)
                do m=2,nvar
                   qlocal(m,ii)=qp(m,ip) 
                end do !m
             end do !i
          end do !j
       end do !k

       qmean=0
       do ii=1,npts
          do m=1,nvar
             qmean(m)=qmean(m) + qlocal(m,ii)
          end do !m
       end do !i
       qmean=qmean/npts !area of canonical element in 3D (4 in 2D)

       !Compute Weight
       do m=1,nvar
          qmin(m)=minval(qlocal(m,:))
          theta(m)=min( 1.0, qmean(m)/( abs(qmean(m) - qmin(m))  + eps ) )
       end do !m

       !Apply Weighted combination of Mean Value and High-Order Solution
       do ii=1,npts
          do m=1,nvar
             qlim(m,ii)=qmean(m) + theta(1)*( qlocal(m,ii) - qmean(m) )
             if (qmin(m) == 0 .and. qmean(m) == 0) qlim(m,ii) = 0
          end do !m
       end do !ii

       !Check for Positivity and Store Values
       do ii=1,npts
          ip=inodes(ii)
          qlim(1,ii)=max(0.0, qlim(1,ii))
          if (qlim(1,ii) <= limit_threshold) then
             qlim(2:nvar,ii)=0
          end if
          qp(1,ip)=qlim(1,ii) - q_ref(1,ip)
          do m=2,nvar
             qp(m,ip)=qlim(m,ii)
          end do !m
       end do !ii

    end do !e

  end subroutine limiter_Positivity_Preserving_Inundation

end module mod_limiter
