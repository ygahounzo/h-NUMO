!-------------------------------
!> @brief Data types and constants
!-------------------------------
module mod_types

    implicit none
    private

    public :: i4, i8, r4, c4, r8, c8, r16
    public :: PI, DEG2RAD, RAD2DEG

    integer, parameter :: i4 = SELECTED_INT_KIND(8)
    integer, parameter :: i8 = SELECTED_INT_KIND(13)
    integer, parameter :: r4 = SELECTED_REAL_KIND(6,30)
    integer, parameter :: c4 = SELECTED_REAL_KIND(6,30)
    !> comment in only one of the following lines, not both.   the first is the default.
    !>integer, parameter :: r8 = SELECTED_REAL_KIND(12)   ! real r8
    !>integer, parameter :: r8 = max(SELECTED_REAL_KIND(15),SELECTED_REAL_KIND(32))   ! real r8
#ifdef SINGLE
integer, parameter :: r8 = r4                     ! alias r8 to r4
#else
    integer, parameter :: r8 = SELECTED_REAL_KIND(15,307)   ! real r8
#endif
    integer, parameter :: c8 = SELECTED_REAL_KIND(12)

    integer, parameter :: r16 = max(r8,SELECTED_REAL_KIND(32))   ! real r16

    real(kind=r8), parameter :: PI = 3.1415926535897932346_r8
    real(kind=r8), parameter :: DEG2RAD = PI / 180.0_r8
    real(kind=r8), parameter :: RAD2DEG = 180.0_r8 / PI

end module mod_types


