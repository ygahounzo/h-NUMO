!-----------------------------------------------------------------!
!>@brief Subroutine to exchange stencil data between procs, to replace the
!> file-based communication earlier.
!>@author  Debojyoti Ghosh (ghosh@mcs.anl.gov)
!> Move this to a more-appropriately name file
!-----------------------------------------------------------------!
SUBROUTINE communicate_stencils(irank, root)
  
    USE mod_parallel

    USE mpi
  
    IMPLICIT NONE
  
    !Global Variables
    INTEGER,INTENT(IN) :: irank, root
  
    !Local Variables
    INTEGER :: iproc, i, offset, ierr, total_size
    INTEGER, DIMENSION(:), ALLOCATABLE :: package
    INTEGER :: MPIStatus(MPI_STATUS_SIZE),AllocateStatus
  
  
    !First Executable Statement
    IF (irank == root) THEN ! Root processor
        !   Save its own data
        DO i = 1,num_nbh_l(1)
            nbh_proc(i) = nbh_proc_l(i,1)
        ENDDO
        DO i = 1,nelem_l(1)
            local_global_elem(i) = local_global_elem_l(i,1)
        ENDDO
        DO i = 1,npoin_l(1)
            local_global_poin_periodic(i) = local_global_poin_periodic_l(i,1)
            local_global_poin(i) = local_global_poin_l(i,1)
            ipoin_proc(i) = ipoin_proc_l(i,1)
            flag_periodic(i) = flag_periodic_l(i,1)
        ENDDO
        DO i = 1,npoin_l_cg(1)
            local_global_poin_cg(i) = local_global_poin_l_cg(i,1)
        ENDDO
        DO i = 1,nface_l(1)
            local_global_face(i) = local_global_face_l(i,1)
        ENDDO

        !   Prepare and send the data for the other procs
        DO iproc = 2,nproc
            total_size = num_nbh_l(iproc) + nelem_l(iproc) + 1*nface_l(iproc) + 4*npoin_l(iproc) + npoin_l_cg(iproc)
            ALLOCATE(package(1:total_size),stat=AllocateStatus)
            IF (AllocateStatus /= 0) STOP "** Failed to allocate in subroutine communicate_stencils (1) **"

            offset = 0;
            DO i = 1,num_nbh_l(iproc)
                package(i+offset) = nbh_proc_l(i,iproc)
            ENDDO
            offset = offset + num_nbh_l(iproc)
            DO i = 1,nelem_l(iproc)
                package(i+offset) = local_global_elem_l(i,iproc)
            ENDDO
            offset = offset + nelem_l(iproc)
            DO i = 1,nface_l(iproc)
                package(i+offset) = local_global_face_l(i,iproc)
            ENDDO
            offset = offset + nface_l(iproc)
            DO i = 1,npoin_l(iproc)
                package(i+offset) = local_global_poin_l(i,iproc)
            ENDDO
            offset = offset + npoin_l(iproc)
            DO i = 1,npoin_l(iproc)
                package(i+offset) = local_global_poin_periodic_l(i,iproc)
            ENDDO
            offset = offset + npoin_l(iproc)
            DO i = 1,npoin_l(iproc)
                package(i+offset) = ipoin_proc_l(i,iproc)
            ENDDO
            offset = offset + npoin_l(iproc)
            DO i = 1,npoin_l(iproc)
                package(i+offset) = flag_periodic_l(i,iproc)
            ENDDO
            offset = offset + npoin_l(iproc)
            DO i = 1,npoin_l_cg(iproc)
                package(i+offset) = local_global_poin_l_cg(i,iproc)
            ENDDO
            CALL MPI_Send(package,total_size,MPI_INTEGER,iproc-1,1983,&
                MPI_COMM_WORLD,ierr)
            IF (ierr /= MPI_SUCCESS) THEN
                PRINT*, ' ** ERROR: MPI_Send was not successful in subroutine communicate_stencils. '
                PRINT*, ' ** ERROR: for iproc ', iproc, ' ierr is ', ierr
            ENDIF
            DEALLOCATE(package)
        ENDDO

    ELSE
        !   Other procs: Receive the data from the root proc
        total_size = num_nbh + nelem + 1*nface + 4*npoin + npoin_cg
        ALLOCATE(package(1:total_size),stat=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "** Failed to allocate in subroutine communicate_stencils (2) **"

        CALL MPI_Recv(package,total_size,MPI_INTEGER,0,1983,&
            MPI_COMM_WORLD,MPIStatus,ierr)
        IF (ierr /= MPI_SUCCESS) THEN
            PRINT*, ' ** ERROR: MPI_Recv was not successful in subroutine communicate_stencils. '
            PRINT*, ' ** ERROR: for irank', irank, ' ierr is ', ierr
        ENDIF

        offset = 0;
        DO i = 1,num_nbh
            nbh_proc(i) = package(i+offset)
        ENDDO
        offset = offset + num_nbh
        DO i = 1,nelem
            local_global_elem(i) = package(i+offset)
        ENDDO
        offset = offset + nelem
        DO i = 1,nface
            local_global_face(i) = package(i+offset)
        ENDDO
        offset = offset + nface
        DO i = 1,npoin
            local_global_poin(i) = package(i+offset)
        ENDDO
        offset = offset + npoin
        DO i = 1,npoin
            local_global_poin_periodic(i) = package(i+offset)
        ENDDO
        offset = offset + npoin
        DO i = 1,npoin
            ipoin_proc(i) = package(i+offset)
        ENDDO
        offset = offset + npoin
        DO i = 1,npoin
            flag_periodic(i) = package(i+offset)
        ENDDO
        offset = offset + npoin
        DO i = 1,npoin_cg
            local_global_poin_cg(i) = package(i+offset)
        ENDDO

        DEALLOCATE(package)
    ENDIF
  
END SUBROUTINE communicate_stencils
