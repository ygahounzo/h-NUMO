!---------------------------------------------------------!
!> @brief Write parallel data structures to disk
!---------------------------------------------------------!
subroutine write_parallel()
  
    use mod_parallel, only: nproc, &
        nbh_proc_l, local_global_elem_l, local_global_face_l, local_global_poin_l, ipoin_proc_l , &
        num_nbh_l, nelem_l, nface_l, npoin_l, flag_periodic_l, local_global_poin_periodic_l
  
    implicit none
  
    integer iproc, iloop, j
    character fnp1*72, fnp*72
  
    do iproc = 1,nproc
     
        write(fnp1,'(i5)') iproc
        iloop=4 - int(log10(real(iproc)))
        do j=1,iloop
            fnp1(j:j)='0'
        end do
     
        fnp= 'parallel'// '_' // trim(fnp1) // '.out'
        open(1,file=fnp,form="unformatted")
        write(1) nbh_proc_l(1:num_nbh_l(iproc),iproc)
        write(1) local_global_elem_l(1:nelem_l(iproc),iproc)
        write(1) local_global_face_l(1:nface_l(iproc),iproc)
        write(1) local_global_poin_l(1:npoin_l(iproc),iproc)
        write(1) local_global_poin_periodic_l(1:npoin_l(iproc),iproc)
        write(1) ipoin_proc_l(1:npoin_l(iproc),iproc)
        write(1) flag_periodic_l(1:npoin_l(iproc),iproc)
        close (1)
    end do
  
end subroutine write_parallel

!---------------------------------------------------------!
!> @brief Read parallel data structures from disk
!---------------------------------------------------------!
subroutine read_parallel()
  
    use mod_parallel, only: nbh_proc, local_global_elem, local_global_poin, ipoin_proc, &
        local_global_face, local_global_poin_periodic, flag_periodic

    use mpi

    implicit none
  
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
    read(1) nbh_proc(:)
    read(1) local_global_elem(:)
    read(1) local_global_face(:)
    read(1) local_global_poin(:)
    read(1) local_global_poin_periodic(:)
    read(1) ipoin_proc(:)
    read(1) flag_periodic(:)
    close(1)
  
end subroutine read_parallel

!---------------------------------------------------------!
!> @brief Write grid to disk
!> @todo Not up to date with communicate stencils
!---------------------------------------------------------!
subroutine read_grid()
  
    use mod_grid, only: coord_cg, index2d, intma_table, bsido

    use mod_parallel, only: nproc

    use mpi
  
    implicit none
  
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
    read(1) coord_cg(:,:)
    read(1) index2d(:,:)
    read(1) intma_table(:,:,:,:)
    read(1) bsido(:,:)
    close (1)
  
end subroutine read_grid
