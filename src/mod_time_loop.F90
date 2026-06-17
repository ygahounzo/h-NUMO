!----------------------------------------------------------------------!
!> @brief This module contains the entire time loop of the code.
!>@ author by Yao Gahounzo
!>      Computing PhD
!       Boise State University
!       Date: July 02, 2023
!----------------------------------------------------------------------!
module mod_time_loop

    use mod_mpi_utilities,    only: irank, irank0, MPI_PRECISION, wtime, numproc
    use mod_types,            only: r8
    use mod_grid,             only: grid
    use mod_basis,            only: basis
    use mod_input,            only: input
    use mod_initial,          only: initial
    use mod_metrics,          only: metrics
    use mod_face,             only: face_CS
    use mod_tensor,           only: tensor_CS
    use mod_variables,        only: btp_CS, bcl_CS
    use mod_global_grid,      only: grid_global
    use mod_parallel,         only: parallel_CS
    use mod_ref,              only: mref
    use mod_mpi_communicator, only: mpi_communicator
    use mpi
    use mod_openacc_utilities, only: openacc_initialize, openacc_enter_data

    implicit none

    private

    public :: time_loop, rhs_time, write_time, walltime
    real(kind=r8) :: rhs_time, walltime

contains

    subroutine time_loop(G, inp, b, gg, par, mf, btp, bcl, init, ref, mpic, tsp, mt)

        use mod_restart,  only: restart_mlswe

        implicit none

        type(grid),             intent(in)    :: G
        type(input),            intent(in)    :: inp
        type(basis),            intent(in)    :: b
        type(grid_global),      intent(in)    :: gg
        type(parallel_CS),      intent(in)    :: par
        type(face_CS),          intent(in)    :: mf
        type(btp_CS),           intent(inout) :: btp
        type(bcl_CS),           intent(inout) :: bcl
        type(initial),          intent(in)    :: init
        type(mref),             intent(inout) :: ref
        type(mpi_communicator), intent(inout) :: mpic
        type(tensor_CS),        intent(in)    :: tsp
        type(metrics),          intent(in)    :: mt

        real, dimension(:,:,:), allocatable :: qout_mlswe
        integer :: AllocateStatus

        real :: mass_conserv_l, mass_conserv_g
        real :: mass_conserv0_g(inp%nlayers)
        integer :: m, itime
        integer :: irestart, ntime, inorm
        integer :: iloop, i, j, icol, ivar, npoin_bound, ifnp, idone, ii
        character :: fnp1*4, fnp4*4, fnp2*9, fnp*100, fnpm*72, fnpg*72
        character :: s_layers*3, fnp11*18, fnps*100
        integer :: ilevel

        real :: x, y, z, phi, lambda, r, time
        integer :: ip, ik, ipr, ikr, l, unit0, ierr
        real, dimension(G%npoin, inp%nlayers+1) :: mslwe_elevation

        real(kind=r8) :: time1, time2, tick
        real :: t1, t2
        real, dimension(3, G%npoin, inp%nlayers) :: q_df_read
        real, dimension(3, G%npoin)              :: qb_df_read

        call mpi_comm_rank(mpi_comm_world, irank, ierr)

        ! Setup Number of Time-steps
        time  = 0
        ntime = ceiling(inp%time_final / inp%dt)

        ! Allocate Memory
        allocate(qout_mlswe(5, G%npoin, inp%nlayers),   &
                 stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - MOD_TIME_LOOP:Time_Loop **"

        ! Initialize all layers
        bcl%q_df = init%q_df
        btp%qb_df       = init%qb_df

        ! Initialize/Restart
        inorm  = 0
        itime  = 0
        idone  = 0

        irestart = nint(inp%time_restart / inp%dt)
        if(irestart == 0) irestart = 1

        time = inp%time_initial

        if (inp%time_initial == 0) then

            if (irank == irank0) print *, "Begin Writing"

            fnp1 = '0000'

            if(inp%dump_data) then
                if(trim(inp%out_type) == 'txt') then
                    call diagnostics(G, inp, gg, par, init, qout_mlswe, bcl%q_df, btp%qb_df(1:4,:), itime, idone)
                elseif(trim(inp%out_type) == 'vtk') then
                    do l = 1, inp%nlayers
                        ifnp = l
                        write(fnp1,'(i3)') ifnp
                        iloop = 2 - int(log10(real(ifnp)))
                        do j = 1, iloop
                            fnp1(j:j) = '0'
                        end do
                        write(fnp2,'(a1,a3,a5)') "l", fnp1, "_0000"
                        call write_output_mlswe(G, inp, b, init, gg, par, qout_mlswe(:,:,l), btp%qb_df(1:4,:), fnp2, time, l)
                    end do
                end if
            end if

            if (irank == irank0) print *, "End Writing"

        else

            itime = inp%irestart_file_number
            inorm = inp%irestart_file_number

            if(trim(inp%out_type) == 'txt' .or. trim(inp%out_type) == 'nc') then

                ifnp = inp%irestart_file_number

                write(fnp4,'(i4)') ifnp

                iloop = 3 - int(log10(real(ifnp)))
                do j = 1, iloop
                    fnp4(j:j) = '0'
                end do

                fnp = trim('mlswe') // trim(fnp4)
                if(trim(inp%out_type) == 'nc') fnp = trim('mlswe') // trim(fnp4) // trim('.nc')

                call restart_mlswe(G, inp, b, init, bcl%q_df, btp%qb_df, qout_mlswe, fnp)

            end if

            if(irank == 0) print *, ' Done Reading'
        end if

        if(inp%lcheck_conserved) then

            if (.not. inp%dump_data) then
                call diagnostics(G, inp, gg, par, init, qout_mlswe, bcl%q_df, btp%qb_df(1:4,:), itime, 1)
            end if

            do l = 1, inp%nlayers
                call compute_conserved(G, b, tsp, mass_conserv_l, qout_mlswe(1,:,l))

                mass_conserv_g = 0.0
                call mpi_reduce(mass_conserv_l, mass_conserv_g, 1, MPI_PRECISION, &
                                mpi_sum, 0, mpi_comm_world, ierr)

                mass_conserv0_g(l) = mass_conserv_g
            end do

            fnps = trim('mass_mlswe.cons')

            unit0 = 111
            open(unit=unit0, file=fnps)

            write(s_layers,'(i3)') inp%nlayers
            fnp11 = '(i8,' // trim(adjustl(s_layers)) // '(e16.8,1x))'

        end if

        if (inp%lprint_diagnostics) then
            call print_diagnostics_mlswe(G, inp, b, tsp, init, qout_mlswe, btp%qb_df(1:4,:), time, itime, inp%dt, idone, &
                                         mass_conserv0_g, ntime, fnp11, unit0)
        end if

        if (irank == irank0) then
            print *, "--------------"
            print *, "Begin Time Integration: "
        end if

        call openacc_initialize(par)
        call openacc_enter_data(G, b, mf, btp, bcl, init, tsp, mt, par, ref)

        ! Time Loop

        idone   = 0
        rhs_time = 0.0

        do while (time < inp%time_final)

            itime = itime + 1
            time  = time  + inp%dt

            call mpi_barrier(mpi_comm_world, ierr)

            call cpu_time(time1)

            if (trim(inp%bcl_time_method) == 'rk3') then
              call ti_rk3_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, bcl%q_df, btp%qb_df)
            else
              call ti_2levels_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, bcl%q_df, btp%qb_df)
            endif

            call cpu_time(time2)

            rhs_time = rhs_time + time2 - time1

            if (mod(itime, irestart) == 0 .and. inp%dump_data) then
                inorm = inorm + 1

                ifnp = inorm
                write(fnp1,'(i4)') ifnp
                iloop = 3 - int(log10(real(ifnp)))
                do j = 1, iloop
                    fnp1(j:j) = '0'
                end do

                if(trim(inp%out_type) == 'txt') then
                    call diagnostics(G, inp, gg, par, init, qout_mlswe, bcl%q_df, btp%qb_df(1:4,:), inorm, idone)
                elseif(trim(inp%out_type) == 'vtk') then
                    do l = 1, inp%nlayers
                        ifnp = l
                        write(fnp4,'(i3)') ifnp
                        iloop = 2 - int(log10(real(ifnp)))
                        do j = 1, iloop
                            fnp4(j:j) = '0'
                        end do
                        write(fnp2,'(a1,a3,a1,a4)') "l", fnp4, "_", fnp1
                        call write_output_mlswe(G, inp, b, init, gg, par, qout_mlswe(:,:,l), btp%qb_df(1:4,:), fnp2, time, l)
                    end do
                end if

                if (inp%lprint_diagnostics) then
                    call print_diagnostics_mlswe(G, inp, b, tsp, init, qout_mlswe, btp%qb_df(1:4,:), time, itime, inp%dt, idone, &
                                                 mass_conserv0_g, ntime, fnp11, unit0)
                end if
            end if

        end do

        time2 = wtime()

        close(unit0)

        idone = 1

        if (.not. inp%dump_data) then
            call diagnostics(G, inp, gg, par, init, qout_mlswe, bcl%q_df, btp%qb_df(1:4,:), inorm, idone)
        end if
        call print_diagnostics_mlswe(G, inp, b, tsp, init, qout_mlswe, btp%qb_df(1:4,:), time, itime, inp%dt, idone, &
                                     mass_conserv0_g, ntime, fnp11, unit0)

    end subroutine time_loop

    subroutine write_time()
        open(unit=222, file='time.csv')
        write(222,*) rhs_time
        close(222)
    end subroutine write_time

end module mod_time_loop
