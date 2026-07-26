!----------------------------------------------------------------------!
!>@brief Writes Output
!>@author  J.F. Kelly and F.X. Giraldo
!>Department of Applied Maths
!>Naval Postgraduate School
!>Monterey, CA 93943
!> @date March 12 2015, S. Marras added xml output from Michal's subroutine in the CG-only code
!----------------------------------------------------------------------!

subroutine write_output_mlswe(G, inp, b, init, gg, par, tsp, q_df, qout, qb, fnp1, time, layer)

    use mod_grid,        only: grid
    use mod_input,       only: input
    use mod_basis,       only: basis
    use mod_initial,     only: initial
    use mod_global_grid, only: grid_global
    use mod_parallel,    only: parallel_CS
    use mod_tensor,      only: tensor_CS
    use mod_gradient,    only: compute_vorticity_mlswe
    use mod_constants,   only: gravity

    implicit none

    type(grid),        intent(in) :: G
    type(input),       intent(in) :: inp
    type(basis),       intent(in) :: b
    type(initial),     intent(in) :: init
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par
    type(tensor_CS),   intent(in) :: tsp

    real, intent(in)    :: q_df(inp%nvar_bcl, G%npoin, inp%nlayers), qb(inp%nvar_btp, G%npoin)
    real, intent(inout) :: qout(init%nvar, G%npoin)
    integer, intent(in) :: layer
    real, dimension(:,:), allocatable :: q_aux
    real, dimension(:,:), allocatable :: mslwe_elevation
    real, dimension(:), allocatable :: abs_vort, pot_vort
    real :: time
    logical :: has_w
    integer :: k

    character :: fnp1*9, fnp*100

    allocate(q_aux(init%nvar, G%npoin))
    allocate(mslwe_elevation(G%npoin, inp%nlayers+1))
    allocate(abs_vort(G%npoin), pot_vort(G%npoin))

    has_w = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

    ! Compute the (h,u,v,w/dp,ssh) diagnostic fields for this layer, mirroring
    ! diagnostics.F90's q_df -> q conversion. Written back into qout so it stays
    ! valid for the caller's later print_diagnostics_mlswe call.
    qout(1,:) = (init%alpha_mlswe(layer)/gravity) * q_df(1,:,layer)
    qout(2,:) = q_df(2,:,layer) / q_df(1,:,layer)
    qout(3,:) = q_df(3,:,layer) / q_df(1,:,layer)
    if (has_w) then
        qout(4,:) = q_df(4,:,layer) / q_df(1,:,layer)
    else
        qout(4,:) = q_df(1,:,layer)
    end if

    mslwe_elevation = 0.0
    mslwe_elevation(:,inp%nlayers+1) = init%zbot_df
    do k = inp%nlayers, 1, -1
        mslwe_elevation(:,k) = mslwe_elevation(:,k+1) + (init%alpha_mlswe(k)/gravity)*q_df(1,:,k)
    end do
    qout(5,:) = mslwe_elevation(:,layer)

    ! Absolute vorticity η = ζ + f and potential vorticity q = η/h for this
    ! layer -- see mod_gradient.F90's compute_vorticity_mlswe for the formula.
    ! Gated by inp%lwrite_vorticity since it costs an extra element loop per
    ! layer per dump; outvtk_g_binary_mlswe skips writing these fields when off.
    if (inp%lwrite_vorticity) then
        call compute_vorticity_mlswe(G, b, tsp, init, has_w, qout(2,:), qout(3,:), qout(4,:), qout(1,:), &
                                      abs_vort, pot_vort)
    else
        abs_vort = 0.0
        pot_vort = 0.0
    end if

    fnp = trim(inp%fname_root) // '_' // trim(fnp1) // '.vtk'
    call outvtk_g_binary_mlswe(G, inp, b, init, gg, par, qout, qb, abs_vort, pot_vort, fnp, time)

    deallocate(q_aux, mslwe_elevation, abs_vort, pot_vort)

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
