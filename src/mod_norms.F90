!-----------------------------------------------------------------------------------------!
!This module builds computes a quantity \infty-norm on the whole domain and element-wise
!Written by Simone Marras on 01/2014
!           Department of Applied Mathematics
!           Naval Postgraduate School
!           Monterey, CA 93943-5216
!-----------------------------------------------------------------------------------------!
module mod_norms
  
    use mod_basis, only: npts

    use mod_grid, only: nelem

    public :: &
        mod_norms_create, &
        mod_norms_get_residuals_infty_norms, &
        res_projected, res_br, res_bu, res_bv, res_bw, res_bt, res_btr1, res_bene, &
        norm_inf_br, norm_inf_bu, norm_inf_bv, norm_inf_bw, norm_inf_bt, norm_inf_btr1, norm_inf_bene
    
    private

    !module variables and parameters
    !
    !Local Residuals of the 5 dynamics equations
    real,    dimension(:),   allocatable :: res_br, res_bu, res_bv, res_bw, res_bt, res_btr1, res_bene
    real,    dimension(:,:,:), allocatable :: res_projected

    !element infty norms of the 5 dynamics equations
    real, dimension(:), allocatable :: norm_inf_br, norm_inf_bu, norm_inf_bv, norm_inf_bw
      real, dimension(:), allocatable :: norm_inf_bt, norm_inf_btr1, norm_inf_bene
    !The infty-norm is defined as the maximum value of the residuals at each element node
    contains

  subroutine mod_norms_create()

      implicit none
  
      integer :: AllocateStatus
  
      !Allocate element norms:
      if(allocated(norm_inf_br)) then
          deallocate(norm_inf_br, norm_inf_bu, norm_inf_bv, norm_inf_bw, &
              norm_inf_bt, norm_inf_btr1, norm_inf_bene, &
              res_projected, res_br, res_bu, res_bv, res_bw, res_bt, res_btr1, res_bene)
      endif
      allocate(norm_inf_br(nelem), norm_inf_bu(nelem), norm_inf_bv(nelem), norm_inf_bw(nelem), &
          norm_inf_bt(nelem), norm_inf_btr1(nelem), norm_inf_bene(nelem), stat=AllocateStatus )
      if (AllocateStatus /= 0) stop "** Not Enough Memory to allocate local element &
                              norms - Mod_norms **"
  
      allocate(res_projected(3,nelem,npts), res_br(npts), res_bu(npts), res_bv(npts), &
            res_bw(npts), res_bt(npts), res_btr1(npts), res_bene(npts), stat=AllocateStatus)
      if (AllocateStatus /= 0) stop "** Not Enough Memory to allocate local element &
                              residuals - Mod_norms **"
  
  
      !Initialize arrays to zero:
      !Local infty norms:
      norm_inf_br   = 0.0
      norm_inf_bu   = 0.0
      norm_inf_bv   = 0.0
      norm_inf_bw   = 0.0
      norm_inf_bt   = 0.0
      norm_inf_btr1 = 0.0
      norm_inf_bene = 0.0
  
      !Local residuals of each dynamics equation
      res_br   = 0.0
      res_bu   = 0.0
      res_bv   = 0.0
      res_bw   = 0.0
      res_bt   = 0.0
      res_btr1 = 0.0
      res_bene = 0.0

    end subroutine mod_norms_create


    subroutine mod_norms_get_residuals_infty_norms(ie, ndim, nelem, res_br, res_bu, res_bv, &
        res_bw, res_bt, norm_inf_br, norm_inf_bu, norm_inf_bv, norm_inf_bw, norm_inf_bt)
        !
        ! This subroutine finds the maximum (infty norm) 
        ! of the equations residuals that are stored in create_rhs_volume at every element node
        ! S.Marras, 01/2014
        !
        implicit none

        integer, intent(in)  :: ie, ndim, nelem
        real,    intent(in)  :: res_br(ndim), res_bu(ndim), res_bv(ndim), res_bw(ndim), res_bt(ndim)
        real,    intent(out) :: norm_inf_br(nelem), norm_inf_bu(nelem), norm_inf_bv(nelem)
        real,    intent(out) :: norm_inf_bw(nelem), norm_inf_bt(nelem)

        integer :: AllocateStatus
        norm_inf_br(ie) = maxval(res_br(:))
        norm_inf_bu(ie) = maxval(res_bu(:))
        norm_inf_bv(ie) = maxval(res_bv(:))
        norm_inf_bw(ie) = maxval(res_bw(:))
        norm_inf_bt(ie) = maxval(res_bt(:))

    end subroutine mod_norms_get_residuals_infty_norms

end module mod_norms
