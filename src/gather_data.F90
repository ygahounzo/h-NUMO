subroutine gather_data(q_g,q,nvar)

    use mod_global_grid, only: npoin_g, nelem_g

    use mod_grid, only: npoin

    use mod_parallel, only: local_global_poin_l, npoin_l, npoin_l_max, nproc

    use mod_input, only: space_method

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    integer nvar
    real q(nvar,npoin)
    real q_g(nvar,npoin_g)
    integer displs1(nproc)
    integer ierr, irank, iproc, ip, ip_g
    real, dimension(:,:,:), allocatable :: q_l
    integer,dimension(:),allocatable :: recvcount

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    allocate(recvcount(1:nproc))

    if (irank == 0) then
        allocate (q_l(nvar,npoin_l_max,nproc))
        do iproc = 1,nproc
            displs1(iproc) = nvar*(iproc - 1)*npoin_l_max
            recvcount(iproc) = nvar*npoin_l(iproc)
        end do
    else
        do iproc = 1,nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(q,nvar*npoin,MPI_PRECISION, &
        q_l,recvcount, displs1,MPI_PRECISION, 0, mpi_comm_world,ierr)

    if (irank == 0) then
        if(space_method == 'cgc') then
            do iproc = 1,nproc
                do ip = 1,npoin_l(iproc)
                    ip_g = local_global_poin_l(ip,iproc)
                    q_g(:,ip_g) = q_l(:,ip,iproc)
                end do
            end do
        else
            ip_g=0
            do iproc = 1,nproc
                do ip = 1,npoin_l(iproc)
                    ip_g = ip_g+1
                    q_g(:,ip_g) = q_l(:,ip,iproc)
                end do
            end do
        endif
        deallocate(q_l)
    end if

    deallocate(recvcount)

end subroutine gather_data

subroutine gather_data_column(q_g,q)

    use mod_global_grid, only: nelem_g, ncol_g, nz

    use mod_grid, only: npoin_cg, ncol, node_column

    use mod_parallel, only: local_global_poin_l, npoin_l, npoin_l_max, ncol_l_max, nproc

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    real, intent(in)  :: q(ncol)
    real, intent(out) :: q_g(ncol_g)
    integer           :: displs1(nproc)
    integer           :: ncol_l
    integer           :: ierr, irank, iproc, ip, ip_g

    real, dimension(:,:), allocatable :: q_l
    integer,dimension(:),allocatable  :: recvcount

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    allocate(recvcount(1:nproc))

    if (irank == 0) then
        allocate (q_l(ncol_l_max,nproc))
        do iproc = 1,nproc
            displs1(iproc) = (iproc - 1)*ncol_l_max
            recvcount(iproc) =  npoin_l(iproc)/nz !ncol_l(iproc) <--> npoin_l(iproc)/nz
        end do
    else
        do iproc = 1,nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(q,ncol,MPI_PRECISION, &
        q_l,recvcount, displs1,MPI_PRECISION, 0, mpi_comm_world,ierr)

    if (irank ==0) then
        do iproc = 1,nproc
            ncol_l = npoin_l(iproc)/nz
            do ip = 1,ncol_l
                ip_g = node_column(ip,1)

                q_g(ip_g) = q_l(ip,iproc)
            end do
        end do

        deallocate(q_l)
    end if
    deallocate(recvcount)

end subroutine gather_data_column


subroutine gather_data_column_p4est(q_g, q, nvar)
    !
    ! Added by Simone marras.
    ! It's a modification of
    ! gather_data_column_v0. This version uses also nvar.
    !
    ! October 2014
    !
    ! Modified by MA Kopera to not use local_global_poin_l

    use mod_global_grid, only: npoin_g, nelem_g, ncol_g

    use mod_grid, only: npoin, ncol, nz

    use mod_parallel, only: npoin_l, npoin_l_max, ncol_l_max, nproc, ncol_l

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    integer, intent(in)  :: nvar
    real,    intent(in)  :: q(nvar,ncol)
    real,    intent(out) :: q_g(nvar,ncol_g)

    integer              :: displs1(nproc)
    integer              :: ierr, irank, iproc, ip, ip_g, ip_aux

    real, dimension(:,:,:), allocatable :: q_l
    integer,dimension(:),   allocatable :: recvcount

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    allocate(recvcount(1:nproc))
  
    !nz corresponds to npoin_r for the sphere:

    if (irank == 0) then
        allocate (q_l(nvar,ncol_l_max,nproc))
        do iproc = 1,nproc
            displs1(iproc) = nvar*(iproc - 1)*ncol_l_max
            recvcount(iproc) =  nvar*ncol_l(iproc)
        end do
    !     print*,"ncol_l_max",ncol_l_max
    !     print*,"recvcount:",recvcount
    !     print*,"displ:",displs1
     
    else
        do iproc = 1,nproc
            recvcount(iproc) = 0
        end do
    !     print*,"nvar*ncol",nvar*ncol
    end if
  
    ! Gather all the data onto the head node
    !q_l

    call mpi_gatherv(q,nvar*ncol,MPI_PRECISION, &
        q_l,recvcount, displs1,MPI_PRECISION, 0, mpi_comm_world,ierr)
  
    if (irank ==0) then
        ip_g=0
        do iproc = 1,nproc
            do ip = 1,ncol_l(iproc)
                ip_g = ip_g+1
                q_g(:,ip_g) = q_l(:,ip,iproc)
            end do
        end do

        deallocate(q_l)
    end if

    deallocate(recvcount)
  
end subroutine gather_data_column_p4est


!>@author M.A.Kopera
!>Added on 01/06/2015 to gather data without local_global
subroutine gather_data_p4est(q_g,q,nvar)

    use mod_grid, only: npoin

    use mod_global_grid, only: npoin_g

    use mod_parallel, only: npoin_l, npoin_l_max, nproc

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    integer nvar
    real q(nvar,npoin)
    real q_g(nvar,npoin_g)
    integer displs1(nproc)
    integer ierr, irank, iproc, ip, ip_g
    real, dimension(:,:,:), allocatable :: q_l
    integer,dimension(:),allocatable :: recvcount

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    allocate(recvcount(1:nproc))

    if (irank == 0) then
        allocate (q_l(nvar,npoin_l_max,nproc))
        do iproc = 1,nproc
            displs1(iproc) = nvar*(iproc - 1)*npoin_l_max
            recvcount(iproc) = nvar*npoin_l(iproc)
        end do
    else
        do iproc = 1,nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(q,nvar*npoin,MPI_PRECISION, &
        q_l,recvcount, displs1,MPI_PRECISION, 0, mpi_comm_world,ierr)

    if (irank ==0) then
        ip_g=0
        do iproc = 1,nproc
            do ip = 1,npoin_l(iproc)
                ip_g = ip_g+1
                q_g(:,ip_g) = q_l(:,ip,iproc)
            end do
        end do

        deallocate(q_l)
    end if
    deallocate(recvcount)

end subroutine gather_data_p4est

!>@author M.A.Kopera
!>Added on 11/2016 to scatter data without local_global
subroutine scatter_data_p4est(q_g,q,nvar)

  use mod_grid, only: npoin

  use mod_global_grid, only: npoin_g

  use mod_parallel, only: npoin_l, npoin_l_max, nproc

  use mpi

  use mod_mpi_utilities, only: MPI_PRECISION

  implicit none

  integer nvar
  real q(nvar,npoin)
  real q_g(nvar,npoin_g)
  integer displs1(nproc)
  integer ierr, irank, iproc, ip, ip_g
  real, dimension(:,:,:), allocatable :: q_l
  integer,dimension(:),allocatable :: sndcount

  call mpi_comm_rank(mpi_comm_world,irank,ierr)
  allocate(sndcount(1:nproc))

  if (irank == 0) then 
     allocate (q_l(nvar,npoin_l_max,nproc))
     do iproc = 1,nproc
        displs1(iproc) = nvar*(iproc - 1)*npoin_l_max 
        sndcount(iproc) = nvar*npoin_l(iproc)
     end do
  else
     do iproc = 1,nproc
         sndcount(iproc) = 0
     end do
  end if

  if (irank ==0) then
     ip_g=0
     do iproc = 1,nproc
        do ip = 1,npoin_l(iproc)
           ip_g = ip_g+1
           q_l(:,ip,iproc) = q_g(:,ip_g)
        end do
     end do

  end if


  ! Scatter all the data from the head node
  call mpi_scatterv(q_l,sndcount,displs1,MPI_PRECISION, & 
       q,nvar*npoin,MPI_PRECISION, 0, mpi_comm_world,ierr) 

  if(irank==0) deallocate(q_l)

  deallocate(sndcount)

end subroutine scatter_data_p4est


!---------------------------------------------------------------------!  
!> @brief Added on 01/06/2015 to gather data without local_global 
!> from a specified communicator
!> It is a bit rusty for now since it needs a lot of input
!>
!> @author Michal A. Kopera on 01/2015 
!---------------------------------------------------------------------!  
subroutine gather_data_from(q_g,q,nvar,npoin,npoin_g,npoin_l,npoin_l_max,nproc,mpi_comm_custom)

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    integer mpi_comm_custom
    integer npoin_l_max, npoin_l(nproc)
    integer nvar, npoin, npoin_g, nproc
    real q(nvar,npoin)
    real q_g(nvar,npoin_g)
    integer ierr, irank, iproc, ip, ip_g
    real, dimension(:,:,:), allocatable :: q_l
    integer,dimension(:),allocatable :: recvcount, displs1

    call mpi_comm_rank(mpi_comm_custom,irank,ierr)

    allocate(recvcount(1:nproc),displs1(1:nproc))
    if (irank == 0) then
        allocate (q_l(nvar,npoin_l_max,nproc))
        do iproc = 1,nproc
            displs1(iproc) = nvar*(iproc - 1)*npoin_l_max
            recvcount(iproc) = nvar*npoin_l(iproc)
        end do
    else
        do iproc = 1,nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(q,nvar*npoin,MPI_PRECISION, &
        q_l,recvcount, displs1,MPI_PRECISION, 0, mpi_comm_custom,ierr)

    if (irank ==0) then
        ip_g=0
        do iproc = 1,nproc
            do ip = 1,npoin_l(iproc)
                ip_g = ip_g+1
                q_g(:,ip_g) = q_l(:,ip,iproc)
            end do
        end do

        deallocate(q_l)
    end if

    deallocate(recvcount,displs1)

end subroutine gather_data_from

!---------------------------------------------------------------------!  
!> @brief Added on 03/05/2015 to scatter data without local_global 
!> in a specified communicator
!> It is a bit rusty for now since it needs a lot of input
!>
!> @author Michal A. Kopera on 01/2015 
!---------------------------------------------------------------------!  
subroutine scatter_data_to(q_g,q,nvar,npoin,npoin_g,npoin_l,npoin_l_max,nproc,mpi_comm_custom)

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none

    integer mpi_comm_custom
    integer npoin_l_max, npoin_l(nproc)
    integer nvar, npoin, npoin_g, nproc
    real q(nvar,npoin)
    real q_g(nvar,npoin_g)
    integer ierr, irank, iproc, ip, ip_g
    real, dimension(:,:,:), allocatable :: q_l
    integer,dimension(:),allocatable :: sndcount, displs1

    call mpi_comm_rank(mpi_comm_custom,irank,ierr)

    allocate(sndcount(1:nproc),displs1(1:nproc))
    if (irank == 0) then
        allocate (q_l(nvar,npoin_l_max,nproc))
        do iproc = 1,nproc
            displs1(iproc) = nvar*(iproc - 1)*npoin_l_max
            sndcount(iproc) = nvar*npoin_l(iproc)
        end do
    else
        do iproc = 1,nproc
            sndcount(iproc) = 0
        end do
    end if

    if (irank ==0) then
        ip_g=0
        do iproc = 1,nproc
            do ip = 1,npoin_l(iproc)
                ip_g = ip_g+1
                q_l(:,ip,iproc)=q_g(:,ip_g)
            end do
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_scatterv(q_l,sndcount,displs1,MPI_PRECISION, &
        q,nvar*npoin,MPI_PRECISION, 0, mpi_comm_custom,ierr)

    if (irank ==0) then
        deallocate(q_l)
    end if

    deallocate(sndcount,displs1)

end subroutine scatter_data_to



!>@author M.A.Kopera
!>Added on 01/28/2015 to gather integer data without local_global using offset
!>This is used to offset connectivity while assembling data 
subroutine gather_connectivity(conn_g,conn,ncells,ncells_g)

    use mod_basis, only: nglx, ngly, nglz, nsubcells, CELL_CHILDREN

    use mod_grid, only: npoin

    use mod_parallel, only: npoin_l, npoin_l_max, nproc, nelem_l, nelem_l_max

    use mpi

    implicit none

    integer ncells, ncells_g
    integer conn(CELL_CHILDREN,ncells)
    integer conn_g(CELL_CHILDREN,ncells_g)
    integer displs1(nproc)
    integer ierr, irank, iproc, ip, ip_g
    integer, dimension(:,:,:), allocatable :: conn_l
    integer,dimension(:),allocatable :: recvcount
    integer :: offset

    integer :: ncells_l_max
    integer, dimension(:), allocatable :: ncells_l

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    allocate(recvcount(1:nproc), ncells_l(nproc))



    if (irank == 0) then
        ncells_l = nelem_l*nsubcells !number of cells per element = (ngl-1)*(ngl-1)*(ngl-1)
        ncells_l_max = nelem_l_max*nsubcells

        allocate (conn_l(CELL_CHILDREN,ncells_l_max,nproc))

        do iproc = 1,nproc
            displs1(iproc) = CELL_CHILDREN*(iproc - 1)*ncells_l_max
            recvcount(iproc) = CELL_CHILDREN*ncells_l(iproc)
        end do
    else
        do iproc = 1,nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(conn,CELL_CHILDREN*ncells,mpi_integer, &
        conn_l,recvcount, displs1,mpi_integer, 0, mpi_comm_world,ierr)

    if (irank ==0) then
        ip_g=0
        offset=0
        do iproc = 1,nproc
            do ip = 1,ncells_l(iproc)
                ip_g = ip_g+1
                conn_g(:,ip_g) = conn_l(:,ip,iproc)+offset
            end do
            offset = offset+npoin_l(iproc)
        end do

        deallocate(conn_l)
    end if
    deallocate(recvcount,ncells_l)

end subroutine gather_connectivity

!---------------------------------------------------------------------!  
!> @brief Added on 01/30/2015 to gather integer data without 
!> local_global using offset from a specified communicator
!> This is used to offset connectivity while assembling data 
!> It is a bit rusty for now since it needs a lot of input
!>
!> @author Michal A. Kopera on 01/2015 
!---------------------------------------------------------------------!  
subroutine gather_connectivity_from(conn_g,conn,ncells,ncells_g,npoin_l, npoin_l_max, nproc, nelem_l, nelem_l_max, mpi_comm_custom)

    use mod_basis, only: nglx, ngly, nglz, nsubcells, CELL_CHILDREN

    use mpi

    implicit none

    integer ncells, ncells_g
    integer conn(CELL_CHILDREN,ncells)
    integer conn_g(CELL_CHILDREN,ncells_g)
    integer nproc, nelem_l_max, npoin_l_max
    integer nelem_l(nproc), npoin_l(nproc)
    integer mpi_comm_custom

    integer displs1(nproc)
    integer ierr, irank, iproc, ip, ip_g
    integer, dimension(:,:,:), allocatable :: conn_l
    integer,dimension(:),allocatable :: recvcount
    integer :: offset

    integer :: ncells_l_max
    integer, dimension(:), allocatable :: ncells_l

    call mpi_comm_rank(mpi_comm_custom,irank,ierr)
    allocate(recvcount(1:nproc), ncells_l(nproc))



    if (irank == 0) then
        ncells_l = nelem_l*nsubcells !number of cells per element = (ngl-1)*(ngl-1)*(ngl-1)
        ncells_l_max = nelem_l_max*nsubcells

        allocate (conn_l(CELL_CHILDREN,ncells_l_max,nproc))

        do iproc = 1,nproc
            displs1(iproc) = CELL_CHILDREN*(iproc - 1)*ncells_l_max
            recvcount(iproc) = CELL_CHILDREN*ncells_l(iproc)
        end do
    else
        do iproc = 1,nproc
            recvcount(iproc) = 0
        end do
    end if

    ! Gather all the data onto the head node
    call mpi_gatherv(conn,CELL_CHILDREN*ncells,mpi_integer, &
        conn_l,recvcount, displs1,mpi_integer, 0, mpi_comm_custom,ierr)

    if (irank ==0) then
        ip_g=0
        offset=0
        do iproc = 1,nproc
            do ip = 1,ncells_l(iproc)
                ip_g = ip_g+1
                conn_g(:,ip_g) = conn_l(:,ip,iproc)+offset
            end do
            offset = offset+npoin_l(iproc)
        end do

        deallocate(conn_l)
    end if
    deallocate(recvcount,ncells_l)

end subroutine gather_connectivity_from

