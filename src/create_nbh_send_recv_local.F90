!-----------------------------------------------------------------!
!>@brief Creates data structures nbh_send  and
!> num_send_recv, which contain the global grid points and
!> number of global grid points which are sent to neighboring processors
!> This routine create nbh_send in parallel, thus reducing the serial
!> portion of EULER3D-PAR
!>@author James F. Kelly 25 March 2010
!>@date 14 June 2010 Parallel mode
!>@date 28 July 2010 Updated to Unstructured Grids (for sphNUMApar)
!>@date January 2014 F.X. Giraldo Create one sing NBH_SEND_RECV array.
!-----------------------------------------------------------------!
subroutine create_nbh_send_recv_local_cg(nbh_send_recv,num_send_recv,num_bound_max)
  
    use mod_basis, only: ngl, nglx, ngly, nglz
  
    use mod_face, only: ip_bound, num_bound
  
    use mod_global_grid, only: npoin_g_cg, iperiodic_g
  
    use mod_grid, only: npoin_cg, nboun
  
    use mod_parallel, only:  local_global_poin_cg, local_global_poin, local_global_poin_periodic, nproc, nbh_proc, num_nbh, num_nbh_max, nboun_max
  
    use mpi

    implicit none
  
    !global arrays
    integer, intent(in) :: num_bound_max
    integer, dimension(num_bound_max,num_nbh), intent(out) :: nbh_send_recv
    integer, intent(out) :: num_send_recv(num_nbh)

    !local arrays
    integer :: ierr, ireq(2*num_nbh)
    integer :: status(mpi_status_size,2*num_nbh)
    integer :: inbh, iprocr, ip, ip_g, ib, ii, idest, nreq, ipp_g
    integer :: n_a, n_b, n_c
    integer, dimension(num_bound_max) ::  a, b, c, d
    integer :: recv_data(num_nbh*num_bound_max)
    integer :: istart(num_nbh), iend(num_nbh)
    integer :: i_a, i_b, i_c, ibb, i_d, n_d
    integer :: ip1, ip2
    integer :: num_bound_nbh(num_nbh_max)
    integer, dimension(:), allocatable :: global_local, dd, order
    integer :: AllocateStatus

    !Create Global-Local
    allocate(global_local(npoin_g_cg), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop "** Not Enough Memory - Create_Nhb_Send_Local 1**"

    do ip = 1,npoin_cg
        !ip_g = local_global_poin_periodic(ip)
        ip_g = local_global_poin_cg(ip)
        global_local(ip_g) = ip !only store Periodic Point that is Active
    end do
  
    !STEP 1: Send number of boundary points to NBHs
    nreq = 0
    do inbh = 1,num_nbh
     
        !Get Neighboring Processor
        idest = nbh_proc(inbh)
     
        nreq=nreq + 1
        call mpi_isend(num_bound, 1, &
            mpi_integer,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)
     
        !Receive from NBHs
        nreq=nreq + 1
        call mpi_irecv(num_bound_nbh(inbh), 1, &
            mpi_integer,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)
    end do !inbh
  
    call mpi_waitall(nreq,ireq,status,ierr)

    !Form START and END of Send Array
    istart(1) = 1
    iend(1) = istart(1) + num_bound_nbh(1) - 1
    do inbh = 2,num_nbh
        istart(inbh) = iend(inbh-1) + 1
        iend(inbh) =  istart(inbh) + num_bound_nbh(inbh) - 1
    end do
  
    !STEP 2: SEND/RECEIVE BOUNDARY POINTS FOR ALL NBHS
    ! Initialize data structure
    nbh_send_recv = 0; num_send_recv = 0
  
    ! Boundary points on sending processor
    n_a = num_bound
    !  a(1:n_a) = local_global_poin_periodic( ip_bound(1:n_a))
    a(1:n_a) = local_global_poin_cg( ip_bound(1:n_a))
  
    ! Loop over the neighboring processors
    nreq = 0
    do inbh = 1,num_nbh
     
        !Get Neighboring Processor
        idest = nbh_proc(inbh)
     
        ! Boundary point on receiving processor
        nreq=nreq + 1
        call mpi_isend(a, n_a, &
            mpi_integer,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)
     
        !Receive from NBHs
        nreq=nreq + 1
        call mpi_irecv(recv_data(istart(inbh):iend(inbh)), num_bound_nbh(inbh), &
            mpi_integer,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)
    end do
  
    !Wait for everybody to send and receive boundary data
    call mpi_waitall(nreq,ireq,status,ierr)
  
    ! Now loop over neighbors and find intersecton of on-processor boundary points with
    ! neighboring boundary points
    do inbh = 1,num_nbh
        n_b = num_bound_nbh(inbh)
        b(1:n_b) = recv_data(istart(inbh):iend(inbh))
        !Get the intersection of boundary points: neighbors to local processor
        d = 0
        i_d = 0
        do i_b = 1,n_b
            do i_a = 1,n_a
                if (a(i_a) == b(i_b)) then
                    i_d = i_d + 1
                    d(i_d) = b(i_b)
                end if
            end do
        end do
        n_d = i_d

        !Sort the Points
        allocate(dd(n_d), order(n_d), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Create_Nhb_Send_Local 2 **"
        order=0
        do i_d=1,n_d
            dd(i_d)=d(i_d)
        end do
        call quick_sort(dd,order,n_d)

        !Remove redundant Points: from Periodicity, etc.
        i_c=1
        d(i_c)=dd(i_c)
        do i_d=2,n_d
            if ( dd(i_d) /= dd(i_d-1) ) then
                i_c=i_c + 1
                d(i_c)=dd(i_d)
            end if
        end do !i_d
        n_c=i_c

        !Number of Grid points to Send
        num_send_recv(inbh) = n_c
        do ib = 1,n_c
            ip = d(ib)
            nbh_send_recv(ib,inbh) = global_local(ip)
        end do !ib
     
        !deallocate
        deallocate(dd,order)

    end do !inbh

    call mpi_barrier(mpi_comm_world,ierr)

    deallocate(global_local)

end subroutine create_nbh_send_recv_local_cg

!-----------------------------------------------------------------!
!      DG-version
!-----------------------------------------------------------------!
!-----------------------------------------------------------------!
!>@brief Creates data structures nbh_send_recv_dg and
!> num_send, which contain the global faces and
!> number of global faces which are sent to neighboring processors
!> This routine create nbh_send_recv_dg in parallel, thus reducing the serial
!> portion of EULER3D-PAR
!>@author James F. Kelly 25 March 2010
!>@date 14 June 2010 Parallel mode
!>@date 28 July 2010 Updated to Unstructured Grids (for sphNUMApar)
!>@date January 2014 F.X. Giraldo Create one sing NBH_SEND_RECV array.
!>@date March 2015 Daniel S. Abdi Imported from DG-code
!-----------------------------------------------------------------!
subroutine create_nbh_send_recv_local_dg(nbh_send_recv,num_send_recv)
  
    use mod_face, only: face_send
  
    use mod_grid, only: nboun, nface

    use mod_parallel, only: local_global_face, nboun_max
  
    use mod_global_grid, only: nface_g

    use mod_parallel, only:  nproc, nbh_proc, num_nbh, num_nbh_max

    use mpi
  
    implicit none
  
    ! MPI Variables
    integer ierr, ireq(2*num_nbh), irank
    integer status(mpi_status_size,2*num_nbh)
    integer num_send_recv(num_nbh)
    integer nbh_send_recv(nboun,num_nbh)
    integer num_bound_nbh(num_nbh)
    integer inbh, idest, nreq
    integer a(nboun_max), b(nboun_max), c(nboun_max)
    integer recv_data(num_nbh*nboun_max)
    integer istart(num_nbh), iend(num_nbh)
    integer n_a, n_b, n_c
    integer ifaceb , iface, iface_g
    integer ip, ib , i_a, i_b, i_c

    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)

    !STEP 1: Send number of boundary faces to NBHs
    nreq = 0
    do inbh = 1,num_nbh
     
        !Get Neighboring Processor
        idest = nbh_proc(inbh)
   
        nreq=nreq + 1
        call mpi_isend(nboun, 1, &
            mpi_integer,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)
     
        !Receive from NBHs
        nreq=nreq + 1
        call mpi_irecv(num_bound_nbh(inbh), 1, &
            mpi_integer,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)

    end do
 
    ! Wait till num_bound_nbh is populated
    call mpi_waitall(nreq,ireq,status,ierr)

    istart(1) = 1
    iend(1) = istart(1) + num_bound_nbh(1) - 1
    do inbh = 2,num_nbh
        istart(inbh) = iend(inbh-1) + 1
        iend(inbh) =  istart(inbh) + num_bound_nbh(inbh) - 1
    end do

    ! STEP 2: Send Boundary Faces to an Adjacent Processor
    n_a = nboun
    ! Map local faces to Global Faces
    do ifaceb = 1,nboun
        iface = face_send(ifaceb)
        iface_g = local_global_face(iface)
        a(ifaceb) = iface_g
    end do

  
    nreq = 0
    ! Loop over the neighboring processors
    do inbh = 1,num_nbh
     
        !Get Neighboring Processor
        idest = nbh_proc(inbh)
         
        ! Boundary point on receiving processor
        nreq=nreq + 1
        call mpi_isend(a, n_a, &
            mpi_integer,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)
     
        !Receive from NBHs
        nreq=nreq + 1
  

        call mpi_irecv(recv_data(istart(inbh):iend(inbh)), num_bound_nbh(inbh), &
            mpi_integer,idest-1,99,mpi_comm_world, &
            ireq(nreq),ierr)
    end do
  
    !Wait for everybody to send and receive boundary data
    call mpi_waitall(nreq,ireq,status,ierr)
    
   
    do inbh = 1,num_nbh
     
        !Data from a NBH
        n_b = num_bound_nbh(inbh)
        b(1:n_b) = recv_data(istart(inbh):iend(inbh))
     
     
        ! Get the intersection of boundary points on processor
        ! and those of its neighbor
        c = 0
        i_c = 0
        do i_a = 1,n_a
            do i_b = 1,n_b
                if (a(i_a) == b(i_b)) then
                    i_c = i_c + 1
                    c(i_c) = a(i_a)
                end if
            end do
        end do
        n_c = i_c
     
        ! Number of Grid points to Send
        num_send_recv(inbh) = n_c
        do ib = 1,n_c
            ip = c(ib)
            do ifaceb = 1,nboun
                iface = face_send(ifaceb)
                iface_g = local_global_face(iface)
                if (iface_g == ip) then
                    nbh_send_recv(ib,inbh) = iface
                end if
            end do
        end do
     
    end do
  
end subroutine
