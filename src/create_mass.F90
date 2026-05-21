!----------------------------------------------------------------------!
!>@brief This subroutine constructs the Mass matrix for the SEM using
!>Quadrilateral Elements for the Euler equations
!----------------------------------------------------------------------!
subroutine create_mass(G, b, mass, jac)

    use mod_basis, only: basis
    use mod_grid,  only: grid

    implicit none

    type(grid),  intent(in) :: G
    type(basis), intent(in) :: b

    !global arrays
    real, intent(out) :: mass(G%npoin)
    real, intent(in)  :: jac(b%nglx, b%ngly, b%nglz, G%nelem)

    !local
    integer :: ie, i, j, k, ip
    real    :: wq

    !initialize the global matrix
    mass = 0.0

    !loop thru the elements
    do ie = 1, G%nelem
       do k = 1, b%nglz
          do j = 1, b%ngly
             do i = 1, b%nglx
                wq = jac(i,j,k,ie)
                ip = G%intma(i,j,k,ie)
                mass(ip) = mass(ip) + wq
             end do
          end do
       end do
    end do

end subroutine create_mass
