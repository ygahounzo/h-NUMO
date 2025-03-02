!----------------------------------------------------------------------!
!>@brief This module defines the Constants that were in PARAM.H
!>@author F.X. Giraldo
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943-5216
!
!>@modified by Yao Gahounzo, Computing PhD, Boise State University
!> Date: April 03, 2024
!----------------------------------------------------------------------!
module mod_constants

    use mod_types, only : r8

    implicit none

    public

    private r8

    !-----------------------------------------------------------------------
    real(kind=r8), parameter :: tol       = 1e-10

    real(kind=r8), parameter :: pi_trig   = 3.1415926535897932346_r8
    real(kind=r8), parameter :: deg2rad   = pi_trig/180.0_r8
    real(kind=r8), parameter :: rad2deg   = 180.0_r8/pi_trig

    !-------------------------------------------------------------------------------
    ! Physical constants
    !-------------------------------------------------------------------------------
    real(kind=r8), parameter :: omega     = 7.29212e-5_r8 !Earth Angular Velocity
    real(kind=r8), parameter :: earth_rad = 6.37122e6_r8  !Average Earth Radius
    real(kind=r8)            :: gravity   = 9.80616_r8
    !-------------------------------------------------------------------------------
    ! LES constants
    !-------------------------------------------------------------------------------
    real(kind=r8), parameter :: Ck        = 0.15_r8       !Smagorinsky constant
  
    !-------------------------------------------------------------------------------
    ! Derived constants
    !-------------------------------------------------------------------------------

    real(kind=r8) :: earth_radius, pi

    real(kind=r8) :: cm,      &
        km_max,  &
        km_min,  &
        prandtl
    integer nnorm
  !-----------------------------------------------------------------------
  
contains
  
    !-----------------------------------------------------------------------
    subroutine mod_constants_create()
    
        implicit none

        !Initialize
        earth_radius=earth_rad

        pi=pi_trig
        nnorm=1000
        km_min=1.e-2 !

    end subroutine mod_constants_create
  
end module mod_constants
