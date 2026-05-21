!----------------------------------------------------------------------!
!>@brief Writes Output
!>@author  J.F. Kelly and F.X. Giraldo
!>Department of Applied Maths
!>Naval Postgraduate School
!>Monterey, CA 93943
!> @date March 12 2015, S. Marras added xml output from Michal's subroutine in the CG-only code
!----------------------------------------------------------------------!

subroutine write_output_mlswe(G, inp, b, init, gg, par, qp, qb, fnp1, time, layer)

    use mod_grid,        only: grid
    use mod_input,       only: input
    use mod_basis,       only: basis
    use mod_initial,     only: initial
    use mod_global_grid, only: grid_global
    use mod_parallel,    only: parallel_CS

    implicit none

    type(grid),        intent(in) :: G
    type(input),       intent(in) :: inp
    type(basis),       intent(in) :: b
    type(initial),     intent(in) :: init
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par

    real, intent(in)  :: qp(init%nvar, G%npoin), qb(4, G%npoin)
    integer, intent(in) :: layer
    real, dimension(:,:), allocatable :: q_aux, qq, vorticity
    real :: time

    character :: fnp1*9, fnp*100

    allocate(qq(init%nvar, G%npoin), vorticity(3, G%npoin), q_aux(init%nvar, G%npoin))
    qq = qp

    fnp = trim(inp%fname_root) // '_' // trim(fnp1) // '.vtk'
    call outvtk_g_binary_mlswe(G, inp, b, init, gg, par, qq, qb, fnp, time)

    deallocate(qq, q_aux, vorticity)

end subroutine write_output_mlswe

subroutine write_output_mlswe_global(G, inp, b, gg, par, qp, qb, qprime, fnp1, time, layer)

    use mod_grid,        only: grid
    use mod_input,       only: input
    use mod_basis,       only: basis
    use mod_global_grid, only: grid_global
    use mod_parallel,    only: parallel_CS

    implicit none

    type(grid),        intent(in) :: G
    type(input),       intent(in) :: inp
    type(basis),       intent(in) :: b
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par

    real, intent(in)  :: qp(3, G%npoin), qb(3, G%npoin), qprime(3, G%npoin)
    integer, intent(in) :: layer
    real, dimension(:,:), allocatable :: qq
    real :: time

    character :: fnp1*9, fnp*100

    allocate(qq(3, G%npoin))
    qq = qp

    fnp = trim(inp%fname_root) // '_' // trim(fnp1) // '.vtk'
    call outvtk_g_binary_mlswe_global(G, inp, b, gg, par, qq, qb, qprime, fnp, time)

    deallocate(qq)

end subroutine write_output_mlswe_global

!----------------------------!
!> @brief Timestamp
!----------------------------!
subroutine timestamp(stamp)

    implicit none

    character * (12) stamp
    character * ( 8 ) date
    character * ( 10 ) time
    call date_and_time(date, time)
    time = time(1:4)
    stamp = date // time

end subroutine timestamp
