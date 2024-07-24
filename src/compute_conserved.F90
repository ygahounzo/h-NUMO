!----------------------------------------------------------------------!
!>@brief This subroutine calculates the Mass
!>@author  Yao Gahounzo on 06/2024
!>           Computing PhD
!>           Boise State University
!----------------------------------------------------------------------!
subroutine compute_conserved(mass_conserv,q)

    use mod_basis, only: npts

    use mod_grid, only: npoin

    use mod_input, only: nlayers

    use mod_initial, only: psih_df,wjac_df, index_df

    implicit none

    !global arrays
    real, intent(out) :: mass_conserv
    real, intent(in) :: q(npoin)

    !local
    real :: wq, hi
    integer :: Iq, I, ip

    mass_conserv = 0.0

    do Iq = 1, npoin

        wq = wjac_df(Iq)

        do ip = 1,npts

            I = index_df(ip,Iq)
            hi = psih_df(ip,Iq)
            
            mass_conserv = mass_conserv + wq*hi*q(I)

        end do
    end do

    
end subroutine compute_conserved