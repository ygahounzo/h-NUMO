!---------------------------------------------------------!
!> @brief Write parallel data structures to disk
!---------------------------------------------------------!
subroutine write_parallel(par)

    use mod_parallel, only: parallel_CS

    implicit none

    type(parallel_CS), intent(in) :: par

    integer iproc, iloop, j
    character fnp1*72, fnp*72

    do iproc = 1, par%nproc

        write(fnp1,'(i5)') iproc
        iloop=4 - int(log10(real(iproc)))
        do j=1,iloop
            fnp1(j:j)='0'
        end do

        fnp= 'parallel'// '_' // trim(fnp1) // '.out'
        open(1,file=fnp,form="unformatted")
        write(1) par%nbh_proc_l(1:par%num_nbh_l(iproc),iproc)
        write(1) par%local_global_elem_l(1:par%nelem_l(iproc),iproc)
        write(1) par%local_global_face_l(1:par%nface_l(iproc),iproc)
        write(1) par%local_global_poin_l(1:par%npoin_l(iproc),iproc)
        write(1) par%local_global_poin_periodic_l(1:par%npoin_l(iproc),iproc)
        write(1) par%ipoin_proc_l(1:par%npoin_l(iproc),iproc)
        write(1) par%flag_periodic_l(1:par%npoin_l(iproc),iproc)
        close(1)
    end do

end subroutine write_parallel

!---------------------------------------------------------!
!> @brief Read parallel data structures from disk
!---------------------------------------------------------!
subroutine read_parallel(par)

    use mod_parallel, only: parallel_CS
    use mpi

    implicit none

    type(parallel_CS), intent(inout) :: par

    integer iproc, iloop, j, irank, ierr
    character fnp1*72, fnp*72

    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    iproc = irank + 1
    write(fnp1,'(i5)') iproc
    iloop=4 - int(log10(real(iproc)))
    do j=1,iloop
        fnp1(j:j)='0'
    end do

    fnp= 'parallel'// '_' // trim(fnp1) // '.out'

    open(1,file=fnp,form="unformatted")
    read(1) par%nbh_proc(:)
    read(1) par%local_global_elem(:)
    read(1) par%local_global_face(:)
    read(1) par%local_global_poin(:)
    read(1) par%local_global_poin_periodic(:)
    read(1) par%ipoin_proc(:)
    read(1) par%flag_periodic(:)
    close(1)

end subroutine read_parallel

!---------------------------------------------------------!
!> @brief Read grid from disk
!> @todo Not up to date with communicate stencils
!---------------------------------------------------------!
subroutine read_grid(G)

    use mod_grid, only: grid
    use mpi

    implicit none

    type(grid), intent(inout) :: G

    integer iproc, iloop, j, irank, ierr
    character fnp1*72, fnp*72

    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    iproc = irank + 1
    write(fnp1,'(i5)') iproc
    iloop=4 - int(log10(real(iproc)))
    do j=1,iloop
        fnp1(j:j)='0'
    end do

    fnp= 'grid'// '_' // trim(fnp1) // '.gri'

    open(1,file=fnp,form="unformatted")
    read(1) G%coord_cg(:,:)
    read(1) G%index2d(:,:)
    read(1) G%intma_table(:,:,:,:)
    read(1) G%bsido(:,:)
    close(1)

end subroutine read_grid
