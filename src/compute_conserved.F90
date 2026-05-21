!----------------------------------------------------------------------
!>@brief This subroutine calculates the Mass
!>@author  Yao Gahounzo on 06/2024
!>           Computing PhD
!>           Boise State University
!----------------------------------------------------------------------
subroutine compute_conserved(G, b, tsp, mass_conserv, q)

    use mod_grid,    only: grid
    use mod_basis,   only: basis
    use mod_tensor, only: tensor_CS

    implicit none

    type(grid),      intent(in) :: G
    type(basis),     intent(in) :: b
    type(tensor_CS), intent(in) :: tsp

    !global arrays
    real, intent(out) :: mass_conserv
    real, intent(in)  :: q(G%npoin)

    !local
    real :: wq, hi
    integer :: Iq, I, ip

    mass_conserv = 0.0

    do Iq = 1, G%npoin

        wq = tsp%wjac_df(Iq)

        do ip = 1, b%npts

            I  = tsp%index_df(ip,Iq)
            hi = tsp%psih_df(ip,Iq)

            mass_conserv = mass_conserv + wq*hi*q(I)

        end do
    end do

end subroutine compute_conserved
