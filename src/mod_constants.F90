!----------------------------------------------------------------------!
!>@brief This module defines the Constants that were in PARAM.H
!>@author F.X. Giraldo
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943-5216
!>@author S. Marras added nu_air_in
!----------------------------------------------------------------------!
module mod_constants
  
    use mod_input, only : icase, lgravity, p00_in, nu_air_in

    use mod_types, only : r8

    implicit none

    public

    private icase, r8

    !-----------------------------------------------------------------------
    real(kind=r8), parameter :: tol       = 1e-10

    real(kind=r8), parameter :: pi_trig   = 3.1415926535897932346_r8
    real(kind=r8), parameter :: deg2rad   = pi_trig/180.0_r8
    real(kind=r8), parameter :: rad2deg   = 180.0_r8/pi_trig

    !-------------------------------------------------------------------------------
    ! Physical constants
    !-------------------------------------------------------------------------------
    real(kind=r8), parameter :: omega     = 7.29212e-5_r8 !Earth Angular Velocity
    real(kind=r8), parameter :: p00_air   = 1.0e5_r8      !Reference pressure
    real(kind=r8), parameter :: earth_rad = 6.37122e6_r8  !Average Earth Radius
    real(kind=r8), parameter :: cp        = 1004.67_r8    !Specific Heat (constant pressure)
    real(kind=r8), parameter :: cv        = 717.5_r8      !Specific Heat (constant volume)
    real(kind=r8), parameter :: rgas_air   = 287.17_r8    !Dry air gas constant
    real(kind=r8), parameter :: rgasocp   = rgas_air/cp   !R/cp
    real(kind=r8), parameter :: rgasocv   = rgas_air/cv   !R/cv
    real(kind=r8), parameter :: rvap      = 461.50_r8     !Water vapor gas constant
    real(kind=r8), parameter :: tfreeze   = 273.16_r8     !H20 freezing point at p00
    real(kind=r8), parameter :: Lvapor0   = 2.50e6_r8     !Latent Heat of Vapor (@ 0 C)
    real(kind=r8), parameter :: Lvapor100 = 2.25e6_r8     !Latent Heat of Vapor (@ 100 C)
    real(kind=r8), parameter :: Lsub      = 2.85e6_r8     !Latent Heat of Sublimation
    real(kind=r8), parameter :: Lfussion  = 3.34e5_r8     !Latent Heat of Fussion
    real(kind=r8), parameter :: von_karman= 0.4_r8        !von Karman constant
    real(kind=r8), parameter :: scon      = 1362.0_r8     !Solar constant (W/m2)
    real(kind=r8), parameter :: nuair     = 2.0e-5_r8     !Kinematic viscosity of dry air
    real(kind=r8)            :: gravity   = 9.80616_r8
    !-------------------------------------------------------------------------------
    ! LES constants
    !-------------------------------------------------------------------------------
    real(kind=r8), parameter :: Ck        = 0.15_r8       !Smagorinsky constant
  
    !-------------------------------------------------------------------------------
    ! Derived constants
    !-------------------------------------------------------------------------------
    real(kind=r8), parameter :: eps0      = rgas_air/rvap
    real(kind=r8), parameter :: zvirtual  = rvap/rgas_air - 1.0_r8
    real(kind=r8), parameter :: kappa     = rgas_air/cp
    real(kind=r8), parameter :: gamma_air = cp/cv
    real(kind=r8), parameter :: gammam1   = gamma_air - 1.0_r8
    real(kind=r8), parameter :: C0        = (rgas_air**gamma_air)/(p00_air)**(gamma_air-1)

    real(kind=r8) :: earth_radius, p00, rgas, gamma, pi, nu_air

    real(kind=r8) :: cm,      &
        km_max,  &
        km_min,  &
        prandtl
    integer nnorm
  !-----------------------------------------------------------------------
  
contains
  
    !-----------------------------------------------------------------------
    subroutine mod_constants_create()

        use mod_input, only: equations, icase
    
        implicit none

        !Initialize
        earth_radius=earth_rad
        

        p00=p00_air
        nu_air=nuair
        rgas=rgas_air
        pi=pi_trig
        gamma=gamma_air
        nnorm=1000
        cm=0.21 !for mixing
        prandtl=1.0 ! ratio of mass versus heat eddy diffusivities
        km_max=25.0 ! following Durran & Klemp formulation in meso00, meso06
        km_min=1.e-2 !

    end subroutine mod_constants_create
  
end module mod_constants
