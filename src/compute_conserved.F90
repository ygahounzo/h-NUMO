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

! Runtime physical energies, matching energy.py.  Work on the prognostic
! state so diagnostics do not depend on the selected visualization format.
subroutine write_energy_runtime(G, inp, tsp, init, q_df, ape, time, itime, reset)
    use mpi
    use mod_mpi_utilities, only: MPI_PRECISION
    use mod_grid, only: grid
    use mod_input, only: input
    use mod_tensor, only: tensor_CS
    use mod_initial, only: initial
    use mod_constants, only: gravity
    implicit none
    type(grid), intent(in) :: G
    type(input), intent(in) :: inp
    type(tensor_CS), intent(in) :: tsp
    type(initial), intent(in) :: init
    real, intent(in) :: q_df(inp%nvar_bcl,G%npoin,inp%nlayers), time
    ! Chen APE energy (mod_ape%ape_energy, global); used only when c_APE > 0,
    ! adding the columns ape_J, total_aug_J (= total_J + ape_J), relative_change_aug.
    real, intent(in) :: ape
    integer, intent(in) :: itime
    logical, intent(in) :: reset
    real :: energy_global(2), total, relative, relative_aug
    real, save :: initial_total = 0.0, initial_aug = 0.0
    integer, save :: last_step = -1
    integer :: rank, ierr, iu
    logical :: with_ape

    if (reset) last_step = -1
    if (itime == last_step) return ! final diagnostics may repeat the last output step
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    with_ape = (inp%c_APE > 0.0 .and. inp%nlayers > 1)

    call compute_runtime_energy(G, inp, tsp, init, q_df, energy_global)
    if (rank == 0) then
        total = sum(energy_global)
        if (reset) then
            initial_total = total
            initial_aug = total + ape
            open(newunit=iu, file='energy_runtime.csv', status='replace', action='write')
            if (with_ape) then
                write(iu,'(a)') 'step,time_days,kinetic_J,potential_J,total_J,relative_change_total,' // &
                    'ape_J,total_aug_J,relative_change_aug'
            else
                write(iu,'(a)') 'step,time_days,kinetic_J,potential_J,total_J,relative_change_total'
            end if
        else
            open(newunit=iu, file='energy_runtime.csv', status='old', position='append', action='write')
        end if
        relative = 0.0
        if (initial_total /= 0.0) relative = (total-initial_total)/initial_total
        if (with_ape) then
            relative_aug = 0.0
            if (initial_aug /= 0.0) relative_aug = (total+ape-initial_aug)/initial_aug
            write(iu,'(i0,8(",",es24.16e3))') itime, time/86400.0, energy_global, total, relative, &
                ape, total+ape, relative_aug
        else
            write(iu,'(i0,5(",",es24.16e3))') itime, time/86400.0, energy_global, total, relative
        end if
        close(iu)
    end if
    last_step = itime
end subroutine write_energy_runtime

! Globally integrated energies; returned values are defined on MPI rank 0.
subroutine compute_runtime_energy(G, inp, tsp, init, q_df, energy_global)
    use mpi
    use mod_mpi_utilities, only: MPI_PRECISION
    use mod_grid, only: grid
    use mod_input, only: input
    use mod_tensor, only: tensor_CS
    use mod_initial, only: initial
    use mod_constants, only: gravity
    implicit none
    type(grid), intent(in) :: G
    type(input), intent(in) :: inp
    type(tensor_CS), intent(in) :: tsp
    type(initial), intent(in) :: init
    real, intent(in) :: q_df(inp%nvar_bcl,G%npoin,inp%nlayers)
    real, intent(out) :: energy_global(2)
    real :: height(G%npoin), means_local(inp%nlayers+1), means_global(inp%nlayers+1)
    real :: energy_local(2), h, rho, rho_above, drho, speed2
    integer :: k, i, ierr
    means_local = 0.0
    means_local(1) = sum(tsp%wjac_df)
    height = init%zbot_df
    do k = inp%nlayers, 1, -1
        height = height + (init%alpha_mlswe(k)/gravity)*q_df(1,:,k)
        means_local(k+1) = sum(tsp%wjac_df*height)
    end do
    call mpi_allreduce(means_local, means_global, inp%nlayers+1, MPI_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    means_global(2:) = means_global(2:)/means_global(1)

    energy_local = 0.0
    height = init%zbot_df
    do k = inp%nlayers, 1, -1
        rho = 1.0/init%alpha_mlswe(k)
        rho_above = 0.0
        if (k > 1) rho_above = 1.0/init%alpha_mlswe(k-1)
        drho = rho-rho_above
        do i = 1, G%npoin
            h = (init%alpha_mlswe(k)/gravity)*q_df(1,i,k)
            height(i) = height(i)+h
            speed2 = 0.0
            if (q_df(1,i,k) /= 0.0) then
            speed2 = (q_df(2,i,k)/q_df(1,i,k))**2 + (q_df(3,i,k)/q_df(1,i,k))**2
            if (inp%nvar_bcl == 4) speed2 = speed2 + (q_df(4,i,k)/q_df(1,i,k))**2
            end if
            energy_local(1) = energy_local(1)+0.5*tsp%wjac_df(i)*rho*h*speed2
            energy_local(2) = energy_local(2)+0.5*tsp%wjac_df(i)*gravity*drho* &
                             (height(i)-means_global(k+1))**2
        end do
    end do
    call mpi_reduce(energy_local, energy_global, 2, MPI_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
end subroutine compute_runtime_energy

! Rank-0 bookkeeping along the actual filtered trajectory. Reset on each
! invocation (including a restarted run); no counterfactual unfiltered claim.
subroutine write_filter_energy(time, itime, before, after, reset)
    implicit none
    real, intent(in) :: time, before(2), after(2)
    integer, intent(in) :: itime
    logical, intent(in) :: reset
    real, save :: reference_total, previous_total, cumulative_filter, cumulative_solver
    real :: total_before, total_after, delta_filter, delta_solver, filter_rel, solver_rel
    integer :: iu
    if (reset) then
        reference_total = sum(after)
        previous_total = reference_total
        cumulative_filter = 0.0
        cumulative_solver = 0.0
        open(newunit=iu, file='energy_filter.csv', status='replace', action='write')
        write(iu,'(a)') 'step,time_days,kinetic_before_J,potential_before_J,total_before_J,' // &
            'kinetic_after_J,potential_after_J,total_after_J,delta_filter_J,cumulative_delta_filter_J,' // &
            'delta_solver_J,cumulative_delta_solver_J,cumulative_filter_relative_initial,' // &
            'cumulative_solver_relative_initial'
        close(iu)
        return
    end if
    total_before = sum(before)
    total_after = sum(after)
    delta_filter = total_after-total_before
    delta_solver = total_before-previous_total
    cumulative_filter = cumulative_filter+delta_filter
    cumulative_solver = cumulative_solver+delta_solver
    filter_rel = 0.0
    solver_rel = 0.0
    if (reference_total /= 0.0) then
        filter_rel = cumulative_filter/reference_total
        solver_rel = cumulative_solver/reference_total
    end if
    open(newunit=iu, file='energy_filter.csv', status='old', position='append', action='write')
    write(iu,'(i0,13(",",es24.16e3))') itime, time/86400.0, before, total_before, after, total_after, &
        delta_filter, cumulative_filter, delta_solver, cumulative_solver, filter_rel, solver_rel
    close(iu)
    previous_total = total_after
end subroutine write_filter_energy
