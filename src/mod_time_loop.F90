!----------------------------------------------------------------------!
!> @brief This module contains the entire time loop of the code.
!>@ author by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: July 02, 2023
!----------------------------------------------------------------------!
module mod_time_loop

    use mod_mpi_communicator, only: ierr, ireq, nreq, status
    use mod_input, only: dt, lprint_diagnostics, irestart_file_number, out_type
    use mod_mpi_utilities, only : irank, irank0, MPI_PRECISION, wtime, numproc
    use mod_types, only: r8
    use mod_global_grid, only: npoin_g, nelem_g
    use mpi

    implicit none

    private

    public time_loop, rhs_time, write_time, walltime
    real(kind=r8) :: rhs_time, walltime

contains

    subroutine time_loop()

        use mod_grid, only: npoin,  index2d, ncol, coord, npoin_q, nface
        use mod_basis, only: nq
        use mod_initial, only: q_df_mlswe_init, qb_df_mlswe_init, qprime_df_init
        use mod_input, only: dt,time_initial, time_final, time_restart, &
            fname_root, lprint_diagnostics, &
            nlayers, is_mlswe, ti_method_btp, dump_data, lcheck_conserved
        use mod_constants, only: gravity
        use mod_restart, only: restart_mlswe
        use mpi

        implicit none

        real, dimension(:,:,:), allocatable :: q0_df_mlswe, qout_mlswe, qprime0_df
        real, dimension(:,:), allocatable :: qb0_df_mlswe
        integer AllocateStatus
        real mass_conserv_l, mass_conserv_g, mass_conserv0_g(nlayers)
        integer m, itime
        integer irestart, ntime, inorm
        integer iloop, i, j, icol, ivar, npoin_bound, ifnp, idone, ii
        character fnp1*4, fnp4*4, fnp2*9, fnp*100, fnpm*72, fnpg*72, s_layers*3, fnp11*18, fnps*100
        integer ierr, irank, ilevel

        real x, y, z, phi, lambda, r, time
        integer ip, ik, ipr, ikr, l, unit0
        real, dimension(npoin,nlayers+1) :: mslwe_elevation

        real(kind=r8) :: time1, time2, tick
        real :: t1, t2
        real, dimension(3,npoin, nlayers) :: q_df_read
        real, dimension(3,npoin) :: qb_df_read

        call mpi_comm_rank(mpi_comm_world,irank,ierr)

        !Setup Number of Time-steps
        time=0
        ntime=ceiling(time_final/dt)

        ! Allocate Memory
        allocate(q0_df_mlswe(3,npoin,nlayers), qb0_df_mlswe(4,npoin), &
            qout_mlswe(5,npoin, nlayers), qprime0_df(3,npoin,nlayers), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - MOD_TIME_LOOP:Time_Loop:q0_mlswe **"

        !initialize all layers here
        q0_df_mlswe = q_df_mlswe_init
        qb0_df_mlswe = qb_df_mlswe_init
        qprime0_df = qprime_df_init

        !Initialize/Restart
        inorm=0
        itime=0
        idone = 0
    
        irestart=nint(time_restart/dt)
        if(irestart == 0) then
            irestart = 1
        end if

        time = time_initial

        if (time_initial == 0) then


           !write Initial Data
           if (irank == irank0) then
              print *, "Begin Writing"
           end if
           fnp1='0000'

           !write layers output
           if(dump_data) then 
                if(trim(out_type) == 'txt') then
                    call diagnostics(qout_mlswe,q0_df_mlswe, qb0_df_mlswe(1:4,:),itime,idone)
                elseif(trim(out_type) == 'nc') then
                    call diagnostics_nc(qout_mlswe,q0_df_mlswe, qb0_df_mlswe(1:4,:),itime,idone)
                elseif(trim(out_type) == 'vtk') then
                    do l=1,nlayers
                        !Write Snapshot File
                        ifnp= l
                        write(fnp1,'(i3)')ifnp
                        iloop=2 - int(log10(real(ifnp)))
                        do j=1,iloop
                        fnp1(j:j)='0'
                        end do
                        write(fnp2,'(a1,a3,a5)')"l",fnp1,"_0000"
                        call write_output_mlswe(qout_mlswe(:,:,l),qb0_df_mlswe(1:4,:),fnp2,time,l)

                    end do
                end if
            end if

           if (irank == irank0) then
              print *, "End Writing"
           end if

        else 

            itime = irestart_file_number
            inorm = irestart_file_number
           
            !Read Snapshot File
            if(trim(out_type) == 'txt' .or. trim(out_type) == 'nc') then 

                ifnp = irestart_file_number

                write(fnp4,'(i4)')ifnp

                iloop=3 - int(log10(real(ifnp)))
                do j=1,iloop
                fnp4(j:j)='0'
                end do

                fnp=trim('mlswe') // trim(fnp4)
                if(trim(out_type) == 'nc') fnp = trim('mlswe') // trim(fnp4)//trim('.nc')

                call restart_mlswe(q0_df_mlswe,qb0_df_mlswe,qprime0_df, qout_mlswe, fnp)

                qb0_df_mlswe = qb0_df_mlswe
            end if
           
           if(irank==0)  print*,' Done Reading'
        end if

        if(lcheck_conserved) then 

            if (.not. dump_data) then
                call diagnostics(qout_mlswe,q0_df_mlswe, qb0_df_mlswe(1:4,:),itime,1)
            end if

            do l = 1,nlayers
                call compute_conserved(mass_conserv_l,qout_mlswe(1,:,l))

                mass_conserv_g = 0.0
                call mpi_reduce(mass_conserv_l,mass_conserv_g,1,MPI_PRECISION, &
                                mpi_sum,0,mpi_comm_world,ierr)

                mass_conserv0_g(l) = mass_conserv_g
            end do 

            fnps = trim('mass_mlswe.cons')

            unit0 = 111
            open(unit=unit0, file = fnps)
            
            write(s_layers,'(i3)')nlayers

            !fnp4 = trim("(") // trim(nlayers) // trim("e24.12)")
            fnp11 = '(i8,' // trim(adjustl(s_layers)) // '(e16.8,1x))'

            !write(unit0,fnp11) mass_conserv0_g
            !write(111,"(2e24.12)") mass_conserv0_g

            !close(111)

        endif 

        if (lprint_diagnostics) then
            call print_diagnostics_mlswe(qout_mlswe,qb0_df_mlswe(1:4,:),time,itime,dt,idone, &
            mass_conserv0_g,ntime,fnp11,unit0)
        end if

        if (irank == irank0) then
            print *, "--------------"
            print *, "Begin Time Integration: "
        end if

        !     Time Loop

        idone=0
        rhs_time = 0.0
        
        do while (time < time_final)

            !--- Update time counters
            itime = itime + 1
            time = time + dt

            call mpi_barrier(mpi_comm_world,ierr)

            !time1 = wtime()
            call cpu_time(time1)

            call ti_rk_bcl(q0_df_mlswe, qb0_df_mlswe, qprime0_df)

            call cpu_time(time2)

            rhs_time = rhs_time + time2-time1
            !rhs_time = time2-time1

            !rhs_time = rhs_time + (wtime() - time2)
            
            !Write Out Restart File
            if (mod(itime,irestart) == 0 .and. dump_data) then
                inorm=inorm + 1

                !Write Snapshot File
                ifnp= inorm
                write(fnp1,'(i4)')ifnp
                iloop=3 - int(log10(real(ifnp)))
                do j=1,iloop
                    fnp1(j:j)='0'
                end do

                if(trim(out_type)=='txt') then
                    call diagnostics(qout_mlswe,q0_df_mlswe,qb0_df_mlswe(1:4,:),inorm,idone)
                elseif(trim(out_type)=='nc') then
                    call diagnostics_nc(qout_mlswe,q0_df_mlswe,qb0_df_mlswe(1:4,:),inorm,idone)
                elseif(trim(out_type)=='vtk') then
                    do l=1,nlayers
                        !Write Snapshot File
                        ifnp= l
                        write(fnp4,'(i3)')ifnp
                        iloop=2 - int(log10(real(ifnp)))
                        do j=1,iloop
                            fnp4(j:j)='0'
                        end do
                        write(fnp2,'(a1,a3,a1,a4)')"l",fnp4,"_",fnp1

                        call write_output_mlswe(qout_mlswe(:,:,l),qb0_df_mlswe(1:4,:),fnp2,time,l)
                    end do
                end if

                !Append q max and q min to file:
                if (lprint_diagnostics) then 
                    call print_diagnostics_mlswe(qout_mlswe,qb0_df_mlswe(1:4,:),time,itime,dt,idone,&
                    mass_conserv0_g,ntime,fnp11,unit0)
                end if
            end if !mod
        end do !time

        time2 = wtime()

        close(unit0)

        idone = 1

        !if (lprint_diagnostics) then 
            if (.not. dump_data) then
                call diagnostics(qout_mlswe,q0_df_mlswe,qb0_df_mlswe(1:4,:),inorm,idone)
            end if
            call print_diagnostics_mlswe(qout_mlswe,qb0_df_mlswe(1:4,:),time,itime,dt,idone,&
            mass_conserv0_g,ntime,fnp11,unit0)
        !end if

        if(allocated(q0_df_mlswe)) deallocate(q0_df_mlswe)
        if(allocated(qb0_df_mlswe)) deallocate(qb0_df_mlswe)
        if(allocated(qprime0_df)) deallocate(qprime0_df)
        
    end subroutine time_loop

    subroutine write_time()
       open(unit=222,file='time.csv')
       
       write(222,*) rhs_time
       close(222)

    end subroutine write_time

end module mod_time_loop
