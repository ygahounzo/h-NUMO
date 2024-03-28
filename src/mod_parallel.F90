!----------------------------------------------------------------------!
!>@brief This module builds the Parallel (MPI) Communication
!>@author  Jim Kelly in 2010
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>@date F.X. Giraldo 1/2014 to include variable NGL orders and to use INTENT 
!>of passed arrays
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>@date S. Marras on 03/2014 to add ncol_l_max
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>@date S. Marras on 06/2014 to add nboun_poin_l
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>
!----------------------------------------------------------------------!
module mod_parallel
  
    use mod_basis, only: nglx, ngly, nglz
  
    use mod_global_grid, only: nelem_g, npoin_g, npoin_g_cg, nface_g, ncol_g, intma_g, &
        nelem_s, ele_col_g, flag_periodic_g, nz, nboun_poin_g
  
    use mod_input, only: decomp_type, geometry_type, nelx, nely, nelz
  
    public :: &
        ! Constructor Functions
        mod_parallel_create, &
        mod_parallel_create_nbh_proc, &
        mod_parallel_scatter, &
        mod_parallel_create_send_recv, &
        mod_parallel_print_status, &
        mod_parallel_create_ipoin_proc, &
        mod_parallel_destroy, &
        mod_parallel_create_nbh_send, &
        mod_parallel_reorder, &

        ! Variables
        nproc, &
        nelem_l, &
        nelemx_l, &
        nelemy_l, &
        nelemz_l, &
        npoin_l, &
        npoin_l_cg, &
        ncol_l,  &
        nbsido_l, &
        nbound_poin_l, &
        nface_l, &
        nface_boundary_l, &
        nelem_l_max, &
        nelemx_l_max, &
        nelemy_l_max, &
        nelemz_l_max, &
        npoin_l_max, &
        nface_l_max, &
        ncol_l_max, &
        nboun_max, &
        nboun_poin_l_max, &
        num_nbh_max, &
        face_type_l,&
        local_global_face_l, &
        local_global_elem_l, &
        global_proc, &
        local_global_poin_l_cg, &
        local_global_poin_l, &
        local_global_poin_periodic_l, &
        flag_periodic_l, &
        bsido_l, &
        nbh_proc_l, &
        num_nbh_l, &
        ipoin_proc_l, &
        num_nbh, &
       
        local_global_elem, &
        local_global_poin, &
        local_global_poin_cg, &
        local_global_poin_periodic, &
        local_global_face, &
        flag_periodic, &
        num_proc, &
        nbh_proc, &
        ipoin_proc, &
       
        num_send_recv, &
        num_send_recv_total, &
        nbh_send_recv, &
        nbh_send_recv_multi, &
        nbh_send_recv_half, &
        npoin, &
        npoin_cg, &
        nelem, &
        nelemx,&
        nelemy,&
        nelemz,&
        nface
  
    private
    !-----------------------------------------------------------------------
    !Arrays on the head node
    integer nproc, num_nbh, num_send_total, num_send_recv_total
    integer, dimension(:),     allocatable :: nelem_l, npoin_l, npoin_l_cg, nface_l, nface_boundary_l, nbsido_l
    integer, dimension(:),     allocatable :: nboun_poin_l
    integer, dimension(:),     allocatable :: nelemx_l, nelemy_l, nelemz_l
    integer, dimension(:),     allocatable :: ncol_l
    integer, dimension(:,:),   allocatable :: local_global_elem_l, global_local_elem_l, local_global_face_l,face_type_l
    integer, dimension(:),     allocatable :: global_proc
    integer, dimension(:,:),   allocatable :: local_global_poin_l_cg, local_global_poin_l, local_global_poin_periodic_l, global_local_poin_l
    integer, dimension(:,:),   allocatable :: flag_periodic_l
    integer, dimension(:,:),   allocatable :: nbh_proc_l
    integer, dimension(:),     allocatable :: num_nbh_l
    integer, dimension(:,:,:), allocatable :: bsido_l
    integer, dimension(:,:),   allocatable :: ipoin_proc_l
  
    !Arrays on individual processors
    integer, dimension(:),   allocatable :: local_global_elem, global_local_elem
    integer, dimension(:),   allocatable :: local_global_poin_cg, local_global_poin, local_global_poin_periodic,local_global_face
    integer, dimension(:),   allocatable :: flag_periodic
    integer, dimension(:),   allocatable :: nbh_proc
    integer, dimension(:),   allocatable :: ipoin_proc
    integer, dimension(:),   allocatable :: num_send_recv
    integer, dimension(:),   allocatable :: nbh_send_recv
    integer, dimension(:),   allocatable :: nbh_send_recv_multi
    integer, dimension(:),   allocatable :: nbh_send_recv_half
    integer npoin, npoin_cg, nelem, nelemx, nelemy, nelemz
    integer ncol_l_max
  !-----------------------------------------------------------------------
  
contains
  
    !----------------------------------------------------------------------!
    !>@brief This subroutine creates the local/global mappings using METIS
    !----------------------------------------------------------------------!
    subroutine mod_parallel_create(numproc,nproc_z)
    
        implicit none
    
        !global
        integer, intent(in) :: numproc, nproc_z

        !local
        integer :: AllocateStatus
        integer :: ierr, iproc
        integer :: nprocx, nprocy, nprocz
        integer :: ie_g, icol, ie_s, i, ix, iy
    
        ! Get total number of processors
        nproc = numproc
    
        allocate (global_proc(nelem_g) , nelem_l(nproc), &
            nelemx_l(nproc),nelemy_l(nproc),nelemz_l(nproc),&
            npoin_l(nproc), npoin_l_cg(nproc), &
            ncol_l(nproc),  &
            nface_l(nproc), nface_boundary_l(nproc), stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Parallel **"
    
        !Initialize
        nprocx=1; nprocy=1; nprocz=1 !These are not updated for METIS decomposition

        !Domain Decomposition
        if (decomp_type(1:7) == 'metis2d') then
            print *, "Performing 2D Domain Decomposition"
            if (nproc < nelem_s) then
                call domain_decomp_metis_adjacency2d(global_proc, nproc)
            else if (nproc == nelem_s) then
                print *, "Constructing Identity Map"
                do ie_g = 1,nelem_g
                    ie_s = ele_col_g(ie_g)
                    global_proc(ie_g) = ie_s
                end do
            else
                print *, "Number of surface elements must be less than or equal to nproc"
                stop
            end if
       
        else if (decomp_type(1:7) == 'metis3d') then
            print *, "Performing 3D Domain Decomposition"
            if (nproc < nelem_g) then
                call domain_decomp_metis_adjacency3d(global_proc, nproc)
            else if (nproc == nelem_g) then
                do ie_g = 1,nproc
                    global_proc(ie_g) = ie_g
                end do
            else
                print *, "Number of elements must be less than or equal to nproc"
                stop
            end if
       
           ! Geometric Decomposition (NSEAM style...requires a cubed grid)
        else if (decomp_type(1:5) == 'nseam') then
            print *, "Performing NSEAM Domain Decomposition"
            if (nproc <= 96) then
                call nseam_decomp(nprocx,nprocy,nproc)
            else
                nprocx = int(sqrt(nproc/6.0))
                nprocy = nprocx
            end if
            call domain_decomp_nseam(global_proc,nprocx,nprocy,nproc)
       
        else if (decomp_type(1:4) == 'geom') then
            print *, "Performing Geometric Domain Decomposition"
            if (decomp_type(1:5) == 'geom1') then
                nprocx = nproc
                nprocy = 1
                nprocz = 1
            else if (decomp_type(1:5) == 'geom2') then
                nprocx = int(sqrt(nproc/1.0))
                nprocy = int(sqrt(nproc/1.0))
                nprocz = 1
                if (nprocx*nprocy /= nproc) then
                    !check factors
                    do i=2,8
                        if (mod(nproc,i) == 0) then
                            ix=i
                            iy=nproc/ix
                            if ( mod(nelx,ix) == 0 .and. mod(nely,iy) == 0 ) then
                                nprocx=ix
                                nprocy=iy
                            end if
                        end if
                    end do
                end if
            else if (decomp_type(1:5) == 'geom3') then
                nprocx = int(nproc**(1/3))
                nprocy = nprocx
                nprocz = nprocx
            else if (decomp_type(1:5) == 'geom4') then
                nprocz=nproc_z
                nprocx = int(sqrt(1.0*nproc/nprocz))
                nprocy = int(sqrt(1.0*nproc/nprocz))
            else if (decomp_type(1:5) == 'geom5') then
                nprocx = 1
                nprocy = nproc
                nprocz = 1
            else if (decomp_type(1:5) == 'geom6') then
                nprocx = 1
                nprocy = int(sqrt(nproc/1.0))
                nprocz = int(sqrt(nproc/1.0))
            else if (decomp_type(1:5) == 'geom7') then
                nprocx = int(sqrt(nproc/1.0))
                nprocy = 1
                nprocz = int(sqrt(nproc/1.0))
            else if (decomp_type(1:5) == 'geom8') then
                nprocx = 1
                nprocy = 1
                nprocz = nproc
            else if (decomp_type(1:5) == 'geom9') then
                nprocx = int(1.0*nproc/nproc_z)
                nprocy = 1
                nprocz = nproc_z
            end if
            print*,' nprocx nprocy nprocz = ',nprocx,nprocy,nprocz
            call domain_decomp_cube(global_proc,nprocx,nprocy,nprocz,nproc)
       
        else
            print*,' Error in MOD_PARALLEL.F90'
            print*,' Unknown decomp_type = ',decomp_type
            stop
        end if
    
        !Print Global_Proc
        open(1,file='global_proc.out')
        do ie_g=1,nelem_g
            write(1,'(" ie global_proc(ie_g) = ",2(i6,1x))')ie_g,global_proc(ie_g)
        end do
        close(1)

        !Count the number of elements on each processor
        call count_elements(nelem_l)
    
        !Count the number of elements in every cardinal direction, for each processor
        call count_elements_xyz(nelemx_l,nelemy_l,nelemz_l, nprocx, nprocy, nprocz)
    
        ! Count Total Faces
        call count_faces(nface_l,nface_boundary_l)
    
        ! Get maximum number of elements on any processor
        nelem_l_max  = maxval(nelem_l)
        nelemx_l_max = maxval(nelemx_l)
        nelemy_l_max = maxval(nelemy_l)
        nelemz_l_max = maxval(nelemz_l)
        nface_l_max  = maxval(nface_l)
        nboun_max    = maxval(nface_boundary_l)
        npoin_l_max  = nelem_l_max*nglx*ngly*nglz
        ncol_l_max   = (nelem_l_max/nelz)*nglx*ngly

        allocate(local_global_elem_l(nelem_l_max,nproc) , local_global_poin_l(npoin_l_max,nproc), &
            local_global_poin_l_cg(npoin_l_max,nproc), local_global_poin_periodic_l(npoin_l_max,nproc), &
            nbsido_l(nproc), flag_periodic_l(npoin_l_max,nproc), &
            local_global_face_l(nface_l_max,nproc), face_type_l(nface_l_max,nproc),&
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Parallel**"
    
        ! Create Local/Global Mappings (Element-Wise)
        print *, "Creating Element Local/Global Mappings"
        call create_local_global_elem(local_global_elem_l,nelem_l_max, nproc,nelem_l,global_proc)
        ! Create Point-wise Local/Global Mappings (Point-Wise)
        print *, "Create GridPoint Local/Global Mappings"
        call create_local_global_poin(local_global_poin_l,local_global_poin_l_cg, local_global_poin_periodic_l,local_global_elem_l,  &
            flag_periodic_l,npoin_l_max,npoin_l,npoin_l_cg,nelem_l,nelem_l_max,nproc)
        ! Create Local/Global Mappings (Face-Wise)
        print *, "Create Face Local/Global Mappings"
        call create_local_global_face(local_global_face_l,nface_l,face_type_l, &
            nface_l_max,global_proc,nproc)
        ! Create coord_l,intma_l, and bsido_l on the head node
        print *, "Store Domain Decomposition to Local Arrays on Root Proc to be Scattered"
        call domain_decomp_store_local(npoin_l, npoin_l_cg,nbsido_l,npoin_l_max, nface_l,nface_l_max,&
            nelem_l,nelem_l_max,nelemx_l,nelemx_l_max,nelemy_l, nelemy_l_max, &
            nelemz_l,nelemz_l_max,nboun_max,nproc,local_global_elem_l, &
            local_global_poin_l,local_global_poin_l_cg, local_global_poin_periodic_l,local_global_face_l,global_proc)

    end subroutine mod_parallel_create
  
    !-------------------------------------------------!
    !>@brief Destroy parallel data-structures
    !>@todo This subroutine currently not used
    !-------------------------------------------------!
    subroutine mod_parallel_destroy()
    
        deallocate(local_global_elem_l,nbsido_l,bsido_l,local_global_poin_l,&
            local_global_poin_periodic_l, local_global_face_l)
    
    end subroutine mod_parallel_destroy
  
    !----------------------------------------------------------------------!
    !>@brief This subroutine creates:
    !>NUM_NBH_L(NPROC), NUMB_NBH_MAX, NBH_PROC_L(INBH,NPROC)
    !----------------------------------------------------------------------!
    subroutine mod_parallel_create_nbh_proc(adjncy_proc,xadj_proc,nadj_proc)

        implicit none

        !global arrays
        integer, intent(in) :: nadj_proc
        integer, intent(in) :: adjncy_proc(nadj_proc), xadj_proc(nproc + 1)

        !local arrays
        integer :: num_nbh_max
        integer :: iproc, istart, iend, ieb, inbh
    
        ! Allocate memory for the number of processors
        allocate (num_nbh_l(nproc))
    
        ! Construct NUM_NBH = number of processors adjacent to iproc
        if (nproc > 1) then
            do iproc = 1,nproc
                num_nbh_l(iproc) = xadj_proc(iproc + 1) - xadj_proc(iproc)
            end do
            num_nbh_l(nproc) = num_nbh_l(nproc) + 1
       
            ! Maximum number of neighbors
            num_nbh_max = maxval(num_nbh_l)
       
            allocate(nbh_proc_l(num_nbh_max,nproc))
       
            do iproc = 1,nproc
                istart = xadj_proc(iproc)
                iend = istart + num_nbh_l(iproc) - 1
                inbh = 1
                do ieb = istart,iend
                    nbh_proc_l(inbh,iproc) = adjncy_proc(ieb)
                    inbh = inbh + 1
                end do
            end do
       
        else
            num_nbh_l = 0
        end if
    
    end subroutine mod_parallel_create_nbh_proc
  
    !----------------------------------------------------------------------!
    !>@brief This subroutine scatters the global_local data to each processor.
    !----------------------------------------------------------------------!
    subroutine mod_parallel_scatter()
    
        use mod_input, only : bcast_type
    
        use mpi
    
        use mod_mpi_utilities, only : irank0

        implicit none
    
        !local arrays
        integer :: displs1(nproc), displs2(nproc), displs3(nproc)
        integer :: irank, ierr, iproc, AllocateStatus
    
        call mpi_comm_rank(mpi_comm_world,irank,ierr)
    
        ! Broadcast NUM_NBH_MAX
        if (irank == irank0) then
            num_nbh_max = maxval(num_nbh_l)
        end if
        call mpi_bcast(num_nbh_max,1,mpi_integer,0,mpi_comm_world,ierr)
    
        ! Scatter number of neighbors to each processor
        call mpi_scatter(num_nbh_l,1,mpi_integer, &
            num_nbh,1,mpi_integer,0,mpi_comm_world,ierr)
        ! Scatter the number of local elements to each processor
        call mpi_scatter(nelem_l,1,mpi_integer, &
            nelem,1,mpi_integer,0,mpi_comm_world,ierr)
        ! Scatter the number of local grid points to each processor
        call mpi_scatter(npoin_l,1,mpi_integer, &
            npoin,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_scatter(npoin_l_cg,1,mpi_integer, &
            npoin_cg,1,mpi_integer,0,mpi_comm_world,ierr)
        ! Scatter Number of Faces to Each Processor
        call mpi_scatter(nface_l,1,mpi_integer, &
            nface,1,mpi_integer,0,mpi_comm_world,ierr)

        ! Allocate Memory on each node
        allocate(nbh_proc(num_nbh), local_global_elem(nelem), &
            local_global_poin_cg(npoin_cg), local_global_poin(npoin), local_global_poin_periodic(npoin), &
            ipoin_proc(npoin), flag_periodic(npoin), &
            local_global_face(nface), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Parallel_Scatter **"

        if (bcast_type=='mpi') then
            call communicate_stencils(irank,irank0) !communicate_stencils.f90
        elseif (bcast_type=='file') then
            ! Write Local/Global and Communication Stencils to Disk
            if (irank == irank0) then
                call write_parallel() !FXG: Not up-to-date with COMMUNICATE_STENCILS
            end if
            call mpi_barrier(mpi_comm_world,ierr)
            ! Read Mappings from Disk
            call read_parallel() !FXG: Not up-to-date with COMMUNICATE_STENCILS
        endif
    
        ! Broadcast Global Variables
        call mpi_bcast(nboun_max,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(npoin_g_cg,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(ncol_g,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(nelem_g,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(nface_g,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(nproc,1,mpi_integer,0,mpi_comm_world,ierr)
    end subroutine mod_parallel_scatter
  
    !----------------------------------------------------------------------!
    !>@brief This routine creates nbh_send_recv which is used in the MPI SEND/RECEIVES
    !>NBH_SEND/NBH_RECV are redesigned to have size num_send_total, which reduces the
    !>amt. of memory and the number of page faults.
    !>@author 14 September 2010
    !>@date 01/2014  F.X. Giraldo Rewritten to combine all into NBH_SEND_RECV
    !>@date 02/2015  Daniel S. Abdi Added DG
    !----------------------------------------------------------------------!
    subroutine mod_parallel_create_send_recv(nboun)
    
        use mpi
    
        use mod_input, only: space_method
    
        implicit none

        !local arrays
        integer :: nsize, ierr, AllocateStatus, irank, iproc, ib, ibb, nboun
        integer :: num_bound_max, inbh, ngl2_max
        integer, dimension(:,:), allocatable :: nbh_send_recv_local
    
        allocate(num_send_recv(num_nbh))
    
        if (space_method == 'dg') then
            allocate(nbh_send_recv_local(nboun,num_nbh), stat = AllocateStatus)
            if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Parallel_Create_Send_Recv: 0 **"

            ! Create NBH_SEND = local points to send to NBH
            !        NUM_SEND = number of points to send to NBH
            call create_nbh_send_recv_local_dg(nbh_send_recv_local,num_send_recv)
        else
            ngl2_max=max(nglx*ngly, ngly*nglz, nglx*nglz)
            num_bound_max = nboun_max*ngl2_max
	
            ! Allocate Space for NBH_SEND_RECV
            allocate(nbh_send_recv_local(num_bound_max,num_nbh), stat=AllocateStatus )
            if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Parallel_Create_Send_Recv: 1 **"
	
            ! Create NBH_SEND = local points to send to NBH
            !        NUM_SEND = number of points to send to NBH
            call create_nbh_send_recv_local_cg(nbh_send_recv_local,num_send_recv,num_bound_max)
        end if
    
        ! Total number of points (counting duplicates) or faces (DG) sent to ALL NBHs
        num_send_recv_total = 0
        do inbh = 1,num_nbh
            num_send_recv_total = num_send_recv_total + num_send_recv(inbh)
        end do
	
        !global arrays that are passed via MOD_PARALLEL
        allocate(nbh_send_recv(num_send_recv_total), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Parallel_Create_Send_Recv: 2 **"

        allocate(nbh_send_recv_multi(num_send_recv_total), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Parallel_Create_Send_Recv: 2 **"

        nbh_send_recv_multi=1

        ibb = 0
        do inbh = 1,num_nbh
            do ib = 1,num_send_recv(inbh)
                ibb = ibb + 1
                nbh_send_recv(ibb) = nbh_send_recv_local(ib,inbh)
            end do
        end do

        deallocate(nbh_send_recv_local)

    end subroutine mod_parallel_create_send_recv
  
    !-------------------------------------------
    !-  Reorder face_send
    !-------------------------------------------
    subroutine mod_parallel_reorder(face_send)
      use mod_input, only: space_method
      use mod_mpi_utilities, only: irank

      implicit none

      integer, dimension(:):: face_send
      integer:: ind, jnd, imult, inbh, ib, iface

      if (space_method == 'dg') then
        ind = 0
        jnd = 0
        do inbh = 1,num_nbh
          do ib = 1,num_send_recv(inbh)
            ind = ind + 1

            iface = nbh_send_recv(ind)
            do imult = 1,nbh_send_recv_multi(ind)
              jnd = jnd + 1
              face_send(jnd) = iface
            enddo
          end do
        end do
        if(jnd .ne. size(face_send)) then
          print*,"on rank = ", irank
          print*,"face_send mismatch"
          print*,"filled ", jnd, " of ", size(face_send)
          stop "** mod_parallel_reorder problem"
        endif
      endif
    end subroutine

    !----------------------------------------------------------------------!
    !>@brief Prints Domain Decomposition Status
    !----------------------------------------------------------------------!
    subroutine mod_parallel_print_status()
    
        use mpi

        implicit none

        integer :: ierr, numprocs
    
        call mpi_comm_size(mpi_comm_world,numprocs,ierr)
    
        print *, "Domain Decomposition Performed Sucessfully"
        print *, "No. Global Elements         ", nelem_g
        print *, "No. Global Grid Points      ", npoin_g_cg
        print *, "No. Global COlumns          ", ncol_g
        print *, "No. Global Faces            ", nface_g
        print *, "No. Processors              ", numprocs
        print *, "Maximum No Local Elements   ", nelem_l_max
        print *, "Maximum No Local Gridpoints ", npoin_l_max
        print *, "Maximum No Local Faces      ", nface_l_max
    
    end subroutine mod_parallel_print_status
  
    !----------------------------------------------------------------------!
    !>@brief This subroutine creates ipoin_proc
    !> Call this subroutine on the head node
    !----------------------------------------------------------------------!
    subroutine mod_parallel_create_ipoin_proc()
  
        use mod_input, only: space_method
    
        implicit none

        !local arrays
        integer :: AllocateStatus
        integer :: ip_g, iproc, ipp
        integer, dimension(:),  allocatable :: ipoin_proc_g
    
        allocate( ipoin_proc_l(npoin_l_max,nproc), ipoin_proc_g(npoin_g_cg),  stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Parallel_Create_Ipoin_Proc **"
    
        if(space_method /= 'cgc') then
            ipoin_proc_g = 1
            ipoin_proc_l = 1
        else
            ipoin_proc_g = 0
            ipoin_proc_l = 0

            !Count how many Processors claim a Node
            do iproc=1,nproc
                do ipp=1,npoin_l_cg(iproc)
                    ip_g=local_global_poin_periodic_l(ipp,iproc)
                    ipoin_proc_g(ip_g)=ipoin_proc_g(ip_g) + 1
                end do !ipp
            end do !iproc

            do iproc =1,nproc
                do ipp=1,npoin_l(iproc)
                    ip_g=local_global_poin_periodic_l(ipp,iproc)
                    ipoin_proc_l(ipp,iproc)=ipoin_proc_g(ip_g)
                end do !ipp
            end do !n
        endif
    
        deallocate(ipoin_proc_g)
       
    end subroutine mod_parallel_create_ipoin_proc
  
end module mod_parallel




