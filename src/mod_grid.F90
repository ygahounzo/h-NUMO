!----------------------------------------------------------------------!
!>@brief This module builds the Grid Geometry: INTMA, COORD, etc.
!>@author  J.F. Kelly on 11/2010
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>@date F.X. Giraldo on 3/2011 to include 1D Column-wise Semi-implicit
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>@date Daniel S. Abdi on Dec 2014 -- Jun 2015 Added unified interface 
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!----------------------------------------------------------------------!
module mod_grid

    use mod_basis, only: ngl, nglx, ngly, nglz, nopx, nopy, nopz, FACE_LEN,    &
      P4EST_FACES, P8EST_EDGES, &
      nqx, nqy, nqz, nq

    use mod_input, only: ti_method, si_dimension, geometry_type, delta, decomp_type, &
        space_method, cgdg_method, nelz,         &
        is_non_conforming_flg

    use mod_parallel, only:   nproc, &
        npoin_l_cg,      npoin_l_max, &
        nelem_l,      nelem_l_max, &
        nelemx_l,     nelemy_l,     nelemz_l, &
        nelemx_l_max, nelemy_l_max, nelemz_l_max, &
        bsido_l, nbsido_l, nboun_max, &
        nface_l, nface_boundary_l, face_type_l

      use iso_c_binding, only: C_INT

    public :: &
        mod_grid_create, &
        mod_grid_init_unified, &
        mod_grid_init_coord, &
        mod_grid_init_coord_dg_to_cg, &
        mod_grid_get_face_ngl, &
        mod_grid_get_face_nq, &
        mod_grid_create_column, &
        mod_grid_init_node_column, &
        coord_cg, &
        coord, &
        sigma, &
        index2d, &
        intma_table, &
        intma_dg_quad, &  !< DG quadrature points indexing added by Yao
        NC_face, &
        NC_edge, &
        EToNC, nNC, &
        bsido, &
        face, is_non_conforming, &
        face_type, &
        npoin_cg, npoin_massinv, npoin, nelem, nbsido, npoin_q, & !< npoin_q for DG quadrature points added by Yao
        nface, nboun, &
        nelemx, nelemy, nelemz, &
        node_column_table, &
        node_column_table_dg, &
        ncol, ncol_cg, nz, nz_cg, &
        do_1dIMEX, &
        is_cg_coupled, &
        coupling_flg, &
        dump_field, &
        mod_grid_rotate, &
        mod_grid_deform, &
        edge_limit, &
        edge_mask_dg, &
        face_limit, &
        face_mask_dg

    private
    !-----------------------------------------------------------------------

    !module variables and parameters
    real,    dimension(:,:),      allocatable :: coord_cg, coord
    real,    dimension(:),        allocatable :: sigma
    integer, dimension(:,:),      allocatable :: index2d
    integer, dimension(:,:,:,:),  allocatable :: intma_table
    integer(C_INT), dimension(:,:), allocatable :: NC_face
    integer(C_INT), dimension(:,:), allocatable :: NC_edge
    integer(C_INT), dimension(:),   allocatable :: EToNC
    integer(C_INT)                              :: nNC
    integer, dimension(:,:),      allocatable :: bsido, node_column_table, node_column_table_dg
    integer, dimension(:),        allocatable :: face_type, is_non_conforming
    integer, dimension(:),        allocatable :: coupling_flg
    integer, dimension(:,:),      allocatable :: face
    integer                                   :: npoin_cg, npoin_massinv, npoin, nelem, nbsido, ncol, ncol_cg, nz, nz_cg
    integer                                   :: nface, nboun
    integer                                   :: nelemx, nelemy, nelemz
    integer                                   :: nh_i, nz_i
    logical                                   :: do_1dIMEX
    integer, dimension(12)                    :: edge_limit
    integer, dimension(2,6)                   :: face_limit
    !-----------------------------------------------------------------------

    !-----------------------------------------
    !>@brief Interface for various intmas
    !-----------------------------------------
    abstract interface
        function intma_interface(i,j,k,e) result(r)
            integer,intent(in)::i,j,k,e
            integer::r
        end function intma_interface
        function intma_1d_interface(k,e) result(r)
            integer,intent(in)::k,e
            integer::r
        end function intma_1d_interface
        function node_column_interface(ic,iz) result(r)
            integer,intent(in)::ic,iz
            integer::r
        end function node_column_interface
        function intma_dg_to_cg_interface(index) result(r)
            integer,intent(in)::index
            integer::r
        end function intma_dg_to_cg_interface
    end interface

    public:: intma, intma_cg, intma_dg_to_cg
    public:: intma_1d, intma_1d_cg
    public:: node_column, node_column_cg

    procedure(intma_interface), pointer:: intma
    procedure(intma_1d_interface), pointer:: intma_1d
    procedure(node_column_interface), pointer:: node_column

    procedure(intma_dg_to_cg_interface), pointer::intma_dg_to_cg

contains

    function face_mask_dg(n, m, face) result(ind)
      integer, intent(in) :: n, m, face
      integer :: i, j, k
        select case(face)
        case (1)
          i = 0
          j = n-1
          k = m-1
        case(2)
          i = nglx-1
          j = n-1
          k = m-1
        case(3)
          i = n-1
          j = 0
          k = m-1
        case(4)
          i = n-1
          j = ngly-1
          k = m-1
        case(5)
          i = n-1
          j = m-1
          k = 0
        case(6)
          i = n-1
          j = m-1
          k = nglz-1
        case default
          stop "** face_mask_dg: invalid face"
        end select

        ind = i + j * nglx + k * nglx * ngly + 1
    end

    function edge_mask_dg(n, edge) result(ind)
      integer, intent(in) :: n, edge
      integer :: ind, i, j, k

        select case(edge)
        case (1)
          i = n-1
          j = 0
          k = 0
        case (2)
          i = n-1
          j = ngly-1
          k = 0
        case (3)
          i = n-1
          j = 0
          k = nglz-1
        case (4)
          i = n-1
          j = ngly-1
          k = nglz-1
        case (5)
          i = 0
          j = n-1
          k = 0
        case (6)
          i = nglx-1
          j = n-1
          k = 0
        case (7)
          i = 0
          j = n-1
          k = nglz-1
        case (8)
          i = nglx-1
          j = n-1
          k = nglz-1
        case (9)
          i = 0
          j = 0
          k = n-1
        case (10)
          i = nglx-1
          j = 0
          k = n-1
        case (11)
          i = 0
          j = ngly-1
          k = n-1
        case (12)
          i = nglx-1
          j = ngly-1
          k = n-1
        case default
          stop "** edge_mask_dg: invalid edge"
        end select

        ind = i + j * nglx + k * nglx * ngly + 1
    end function

    !------------------------------------------
    !>@brief intma for CG indexing
    !------------------------------------------
    function intma_cg(i,j,k,e)
        integer,intent(in)::i,j,k,e
        intma_cg = intma_table(i,j,k,e)
    end function intma_cg
    !------------------------------------------
    !>@brief intma for DG indexing
    !------------------------------------------
    function intma_dg(i,j,k,e)
        integer,intent(in)::i,j,k,e
        intma_dg = (e-1) * nglz * ngly * nglx + &
            (k-1) * ngly * nglx + &
            (j-1) * nglx + &
            (i-1) + &
            1
    end function intma_dg

    !------------------------------------------
    !>@brief intma for DG indexing on quads
    !------------------------------------------
    function intma_dg_quad(i,j,k,e)
        integer,intent(in)::i,j,k,e
        intma_dg_quad = (e-1) * nqz * nqy * nqx + &
            (k-1) * nqy * nqx + &
            (j-1) * nqx + &
            (i-1) + &
            1
    end function intma_dg_quad
      !------------------------------------------
      !>@brief intma for 1d-CG indexing
      !------------------------------------------
    function intma_1d_cg(k,e)
        integer,intent(in)::k,e
        intma_1d_cg = (e - 1) * (nglz - 1) + k
    end function intma_1d_cg
    !------------------------------------------
    !>@brief intma for 1d-DG indexing
    !------------------------------------------
    function intma_1d_dg(k,e)
        integer,intent(in)::k,e
        intma_1d_dg = (e - 1) * (nglz) + k
    end function intma_1d_dg
    !------------------------------------------
    !>@brief node_column for 1d-3d-CG indexing
    !------------------------------------------
    function node_column_cg(ic,iz)
        integer,intent(in)::ic,iz
        node_column_cg = node_column_table(iz,ic)
    end function node_column_cg
    !------------------------------------------
    !>@brief node_column for 1d-3d-DG indexing
    !------------------------------------------
    function node_column_dg(ic,iz)
        integer,intent(in)::ic,iz
        node_column_dg = node_column_table_dg(iz,ic)
    end function node_column_dg
    !-----------------------------------
    !>@brief DG to CG conversion for CG0
    !-----------------------------------
    function intma_dg_to_cg_1(index)
        integer,intent(in)::index
        intma_dg_to_cg_1 = index
    end function intma_dg_to_cg_1
    !-----------------------------------
    !>@brief DG to CG conversion for CG/DG
    !-----------------------------------
    function intma_dg_to_cg_2(index)
        integer,intent(in)::index
        integer::i,j,k,e, rem
        rem = index - 1
        e = rem / (nglz * ngly * nglx) + 1
        rem = mod(rem, nglz * ngly * nglx)
        k = rem / (ngly * nglx) + 1
        rem = mod(rem, ngly * nglx)
        j = rem / nglx + 1
        i = mod(rem, nglx) + 1
        intma_dg_to_cg_2 = intma_cg(i,j,k,e)
    end function intma_dg_to_cg_2

    !---------------------------------------------------------
    !>@brief Initialize unifiedCGDG tricks (npoin_cg/npoin) and
    !> intma_cg/intma_dg/intma based on selected space method.
    !---------------------------------------------------------
    subroutine mod_grid_init_unified()

        implicit none

        integer:: AllocateStatus

        integer :: npoin_dg

        if(space_method == 'cgc') then
            !3d
            npoin = npoin_cg
            npoin_dg = nelem * (nopx + 1) * (nopy + 1) * (nopz + 1) !just for coord_dg from p4est
            intma => intma_cg
            intma_dg_to_cg => intma_dg_to_cg_1

            !1D
            nz_cg = nelz * nopz + 1
            nz = nz_cg
            ncol_cg = npoin_cg / nz_cg
            ncol = ncol_cg
            intma_1d => intma_1d_cg

            !1D in 3D
            node_column => node_column_cg
        else
            !3D
            npoin = nelem * (nopx + 1) * (nopy + 1) * (nopz + 1)
            npoin_q = nelem * nqx * nqy * nqz
            intma => intma_dg
            intma_dg_to_cg => intma_dg_to_cg_2

            !1D
            nz_cg = nelz * nopz + 1
            nz = nelz * (nopz + 1)
            ncol_cg = npoin_cg / nz_cg
            ncol = npoin / nz
            intma_1d => intma_1d_dg

            !1D in 3D
            node_column => node_column_dg

            if(is_non_conforming_flg > 0 .and. space_method == 'cgd') then
              npoin_massinv = npoin_cg
            else
              npoin_massinv = npoin
            endif
        endif

        edge_limit = (/nglx, nglx, nglx, nglx, &
                       ngly, ngly, ngly, ngly, &
                       nglz, nglz, nglz, nglz/)
        face_limit = reshape((/ ngly, nglz, ngly, nglz, &
                                nglx, nglz, nglx, nglz, &
                                nglx, ngly, nglx, ngly /), &
                              (/2, 6/))

        !-----------------------
        ! Allocate local arrays
        !-----------------------
        if(allocated(coord_cg)) then
            deallocate(coord_cg, coord,  sigma, &
                index2d, &
                intma_table, &
                NC_face, &
                NC_edge, &
                EToNC, &
                face, &
                face_type, &
                bsido, &
                coupling_flg, &
                is_non_conforming, &
                stat=AllocateStatus )
            if (AllocateStatus /= 0) stop "** Deallocate failed - Mod_Grid **"
        endif
        allocate(coord_cg(3,npoin_cg), coord(3,npoin),  sigma(npoin_cg), &
             index2d(2,npoin_cg), &
            intma_table(nglx,ngly,nglz,nelem), &
            NC_face(P4EST_FACES, nNC),  &
            NC_edge(P8EST_EDGES, nNC),  &
            EToNC(nelem),  &
            face(FACE_LEN,nface), &
            face_type(nface), &
            bsido(6,nbsido), &
            coupling_flg(npoin_cg), &
            is_non_conforming(nface), &
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Grid **"
        EToNC = 0

        coord_cg = 0.0
        sigma = 0.0
        index2d = 0
        intma_table = 0

    end subroutine mod_grid_init_unified

    !-------------------------------------------------
    !>@brief Initialize coord which could use CG/DG storage
    !-------------------------------------------------
    subroutine mod_grid_init_coord()

        implicit none

        integer :: i, j, k, e, ip , ip1

        do e = 1,nelem
            do k = 1,nglz
                do j = 1,ngly
                    do i = 1,nglx
                        ip  = intma_table(i,j,k,e);
                        ip1 = intma(i,j,k,e);
                        coord(:,ip1) = coord_cg(:,ip)
                    enddo
                enddo
            enddo
        enddo
        call initialize_coupling()
    end subroutine mod_grid_init_coord

    !-------------------------------------------------
    !>@brief Initialize coord_cg - temporary subroutine to be replaced by L2 projection and DSS
    !-------------------------------------------------
    subroutine mod_grid_init_coord_dg_to_cg()
        integer :: i, j, k, e, ip , ip1
        if(space_method == 'cgc') then
            coord_cg(:,:) = coord(:,:)
        else
            do i = 1,nglx
                do j = 1,ngly
                    do k = 1,nglz
                        do e = 1,nelem
                            ip = intma_table(i,j,k,e);
                            ip1 = intma_dg(i,j,k,e);
                            coord_cg(:,ip) = coord(:,ip1)
                        enddo
                    enddo
                enddo
            enddo
        endif
        call initialize_coupling()
    end subroutine mod_grid_init_coord_dg_to_cg


    !-------------------------------------------------
    !>@brief Get face ngls
    !-------------------------------------------------
    subroutine mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)
        implicit none
        integer :: ilocl, ngl_i, ngl_j, plane_ij
        if (ilocl == 1 .or. ilocl == 2) then !Zeta=+/-1
            ngl_i=nglx
            ngl_j=ngly
            plane_ij=1
        elseif (ilocl == 3 .or. ilocl == 4) then !Eta=+/-1
            ngl_i=nglx
            ngl_j=nglz
            plane_ij=2
        elseif (ilocl == 5 .or. ilocl == 6) then !Ksi=+/-1
            ngl_i=ngly
            ngl_j=nglz
            plane_ij=3
        end if
    end subroutine mod_grid_get_face_ngl

    !-------------------------------------------------
    !>@brief Get face nqs, added by Yao
    !-------------------------------------------------
    subroutine mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)
        implicit none
        integer :: ilocl, nq_i, nq_j, plane_ij
        if (ilocl == 1 .or. ilocl == 2) then !Zeta=+/-1
            nq_i=nqx
            nq_j=nqy
            plane_ij=1
        elseif (ilocl == 3 .or. ilocl == 4) then !Eta=+/-1
            nq_i=nqx
            nq_j=nqz
            plane_ij=2
        elseif (ilocl == 5 .or. ilocl == 6) then !Ksi=+/-1
            nq_i=nqy
            nq_j=nqz
            plane_ij=3
        end if
    end subroutine mod_grid_get_face_nq

    !------------------------------------------
    !>@brief Mixed CG-DG coupling
    !------------------------------------------
    function is_cg_coupled(i,j,k,e)
        integer,intent(in)::i,j,k,e
        integer ip
        ip = intma_cg(i,j,k,e)
        is_cg_coupled = coupling_flg(ip)
    end function is_cg_coupled

    subroutine initialize_coupling()
        integer ip
        if(space_method == 'cgd' .and. cgdg_method == 'mixed') then
            do ip = 1, npoin_cg
                != Change this problem-dependent parameter
                coupling_flg(ip) = merge(1, 0, coord_cg(1,ip) > 500)
            end do
        else !unified/separate
            do ip = 1, npoin_cg
                coupling_flg(ip) = 1
            end do
        endif
    end subroutine initialize_coupling
    !-----------------------------------------------------------------------
    subroutine mod_grid_create()

        use mod_global_grid, only: saved_coord, saved_sigma, saved_index,&
            saved_intma, saved_bsido, saved_face

        use mod_input, only: bcast_type

        use mpi

        implicit none

        integer :: AllocateStatus, ierr, iproc, irank
        integer :: ngl3, nreq
        integer :: status(mpi_status_size)
        integer :: i,j,k,l
        integer :: MPIStatus(MPI_STATUS_SIZE)

        call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr)

        ! Scatter the number of local elements to each processor
        call mpi_scatter(nelem_l,1,mpi_integer, &
            nelem,1,mpi_integer,0,mpi_comm_world,ierr)

        ! Scatter the number of local elements along x to each processor
        call mpi_scatter(nelemx_l,1,mpi_integer, &
            nelemx,1,mpi_integer,0,mpi_comm_world,ierr)
        ! Scatter the number of local elements along y to each processor
        call mpi_scatter(nelemy_l,1,mpi_integer, &
            nelemy,1,mpi_integer,0,mpi_comm_world,ierr)
        ! Scatter the number of local elements along z to each processor
        call mpi_scatter(nelemz_l,1,mpi_integer, &
            nelemz,1,mpi_integer,0,mpi_comm_world,ierr)

        ! Scatter the number of local grid points to each processor
        call mpi_scatter(npoin_l_cg,1,mpi_integer, &
            npoin_cg,1,mpi_integer,0,mpi_comm_world,ierr)
        ! Scatter the number of physical boundary faces to each processor
        call mpi_scatter(nbsido_l,1,mpi_integer, &
            nbsido,1,mpi_integer,0,mpi_comm_world,ierr)

        ! Scatter Number of Faces to Each Processor
        call mpi_scatter(nface_l,1,mpi_integer, &
            nface,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_scatter(nface_boundary_l,1,mpi_integer, &
            nboun,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_allreduce(nboun,nboun_max,1,mpi_integer, &
            mpi_max, mpi_comm_world,ierr)

        ! init unified
        call mod_grid_init_unified()

        !-------------------------------------
        ! Communicate grid via 'mpi' or 'file
        !-------------------------------------
        if (bcast_type=='mpi') then
            !     Transfer data between root and slave processors
            if (irank == 0) then
                !       Root copies its own data first
                !       The arrays saved_* are being allocated
                !       on the room processor only in domain_decomp()
                !       They are deallocated here
                do i = 1,3
                    do j = 1,npoin_cg
                        coord_cg(i,j) = saved_coord(i,j,1)
                    enddo
                enddo
                do j = 1,nface
                    face(:,j)    = saved_face(:,j,1)
                    face_type(j) = face_type_l(j,1)
                enddo

                do j = 1,npoin_cg
                    index2d(1,j) = saved_index(1,j,1)
                    index2d(2,j) = saved_index(2,j,1)
                    sigma(j)     = saved_sigma(j,1)
                enddo

                do l = 1,nelem
                    do k = 1,nglz
                        do j = 1,ngly
                            do i = 1,nglx
                                intma_table(i,j,k,l) = saved_intma(i,j,k,l,1)
                            enddo
                        enddo
                    enddo
                enddo

                do i = 1,6
                    do j = 1,nbsido
                        bsido(i,j) = saved_bsido(i,j,1)
                    enddo
                enddo
                !       Now root will send to other procs
                do iproc = 2,nproc
                    call MPI_Send(saved_coord(:,1:npoin_l_cg(iproc),iproc:iproc),&
                        3*npoin_l_cg(iproc),&
                        MPI_DOUBLE_PRECISION,iproc-1,6170,MPI_COMM_WORLD,ierr)
                    call MPI_Send(saved_index(:,1:npoin_l_cg(iproc),iproc:iproc),&
                        2*npoin_l_cg(iproc),&
                        MPI_INTEGER,iproc-1,6171,MPI_COMM_WORLD,ierr)
                    call MPI_Send(saved_intma(:,:,:,1:nelem_l(iproc),iproc:iproc),&
                        nglx*ngly*nglz*nelem_l(iproc),&
                        MPI_INTEGER,iproc-1,6172,MPI_COMM_WORLD,ierr)
                    call MPI_Send(saved_bsido(:,1:nbsido_l(iproc),iproc:iproc),&
                        6*nbsido_l(iproc),&
                        MPI_INTEGER,iproc-1,6173,MPI_COMM_WORLD,ierr)
                    call MPI_Send(saved_sigma(1:npoin_l_cg(iproc),iproc:iproc),&
                        npoin_l_cg(iproc),&
                        MPI_DOUBLE_PRECISION,iproc-1,6175,MPI_COMM_WORLD,ierr)
                    call MPI_Send(saved_face(:,1:nface_l(iproc),iproc:iproc),&
                        FACE_LEN*nface_l(iproc),&
                        MPI_INTEGER,iproc-1,6177,MPI_COMM_WORLD,ierr)
                    call MPI_Send(face_type_l(1:nface_l(iproc),iproc:iproc),&
                        nface_l(iproc),&
                        MPI_INTEGER,iproc-1,6178,MPI_COMM_WORLD,ierr)
                enddo
                deallocate(saved_coord,saved_sigma,saved_index,saved_intma,saved_bsido)
            else
                !       Receive data using MPI_Recv()
                call MPI_Recv(coord_cg,3*npoin_cg,MPI_DOUBLE_PRECISION,0,6170,&
                    MPI_COMM_WORLD,MPIStatus,ierr)
                call MPI_Recv(index2d,2*npoin_cg,MPI_INTEGER,0,6171,&
                    MPI_COMM_WORLD,MPIStatus,ierr)
                call MPI_Recv(intma_table,nglx*ngly*nglz*nelem,MPI_INTEGER,0,6172,&
                    MPI_COMM_WORLD,MPIStatus,ierr)
                call MPI_Recv(bsido,6*nbsido,MPI_INTEGER,0,6173,&
                    MPI_COMM_WORLD,MPIStatus,ierr)
                call MPI_Recv(sigma,npoin_cg,MPI_DOUBLE_PRECISION,0,6175,&
                    MPI_COMM_WORLD,MPIStatus,ierr)
                call MPI_Recv(face,FACE_LEN*nface,MPI_INTEGER,0,6177,&
                    MPI_COMM_WORLD,MPIStatus,ierr)
                call MPI_Recv(face_type,nface,MPI_INTEGER,0,6178,&
                    MPI_COMM_WORLD,MPIStatus,ierr)
            endif
        elseif (bcast_type=='file') then
            ! All Slave Nodes Read Grid From Disk
            call read_grid()
        endif

    end subroutine mod_grid_create

    !----------------------------------------------------------------------!
    !>This module builds the Column Storage of INTMA for the 1D SI.
    !>@author  J.F. Kelly on 3/2011
    !>           Department of Applied Mathematics
    !>           Naval Postgraduate School
    !>           Monterey, CA 93943-5216
    !>@date Jun 20, 2015 Daniel S. Abdi
    !> Rewrite for unification
    !----------------------------------------------------------------------!
    subroutine mod_grid_create_column()

        implicit none

        !local
        integer icol, iz, ip, iicol, ii, ie, k, i
        real x, y, z, r, eps, x1, y1, z1, x_dot_x1
        logical isnew
        real, dimension(:), allocatable :: xcol, ycol, zcol

        eps = 1e-14 !!!---- 1e-5, 1e-8, 1e-14, 1e-16

        if(allocated(node_column_table)) then
            deallocate(node_column_table, xcol, ycol, zcol)
        endif
        allocate(node_column_table(nz_cg,ncol_cg))
        allocate(xcol(ncol_cg), ycol(ncol_cg), zcol(ncol_cg))

        !First get the unique cols: Assume everything is numbered in xy-planes.
        if (geometry_type(1:4) == 'cube') then
            do i=1,ncol_cg
                xcol(i)=coord_cg(1,i);
                ycol(i)=coord_cg(2,i);
                zcol(i)=0
            end do
        elseif (geometry_type(1:6) == 'sphere') then
            do i=1,ncol_cg
                x=coord_cg(1,i);
                y=coord_cg(2,i);
                z=coord_cg(3,i)
                r=sqrt( x**2 + y**2 + z**2 )
                xcol(i)=x/r; ycol(i)=y/r; zcol(i)=z/r
            end do
        end if

        !------ Now fill up the local node_column table
        do icol = 1,ncol_cg
            x1 = xcol(icol)
            y1 = ycol(icol)
            z1 = zcol(icol)

            iz = 0
            do ip = 1,npoin_cg
                if (geometry_type(1:4) == 'cube') then
                    x = coord_cg(1,ip)
                    y = coord_cg(2,ip)
                    if ( abs(x - x1) < eps .and. abs(y - y1) < eps) then
                        iz = iz + 1
                        node_column_table(iz,icol) = ip
                    end if
                elseif (geometry_type(1:6) == 'sphere') then
                    x = coord_cg(1,ip)
                    y = coord_cg(2,ip)
                    z = coord_cg(3,ip)
                    r = sqrt( x**2 + y**2 + z**2 )
                    x = x/r; y=y/r; z=z/r
                    x_dot_x1 = x*x1 + y*y1 + z*z1
                    if (abs(x_dot_x1 - 1) < eps) then
                        iz = iz + 1
                        node_column_table(iz,icol) = ip
                    end if
                end if
            end do !ip

            if (iz /= nz_cg) then
                print*,'Error in MOD_GRID_CREATE_COLUMN: iz nz_cg =',iz,nz_cg
                print*,' xcol ycol zcol = ',x1,y1,z1
                stop
            end if

        end do !icol

        deallocate (xcol, ycol, zcol)

    end subroutine mod_grid_create_column

    !-----------------------------------------
    !>@brief Initialize node column table for DG
    !-----------------------------------------
    subroutine mod_grid_init_node_column()

        use mod_input, only: lp6est

        implicit none

        integer ic,iz,i,j,k,l,ecol,ie
        integer nelems

        if(allocated(node_column_table_dg)) deallocate(node_column_table_dg)
        allocate(node_column_table_dg(nz,ncol))

        if (space_method == 'cgc') then
            node_column_table_dg = node_column_table
        else
            do l = 1,nelem
                if(lp6est) then
                    !p4est grid column and element numbering dependent code
                    ecol = (l - 1) / nelz + 1
                    ie = l - ((ecol - 1) * nelz + 1) + 1
                else
                    !calculate on-processor number of elements on a shell
                    nelems = nelem / nelz
                    !numa grid column and element numbering dependent code
                    ecol = mod(l - 1 , nelems) + 1
                    ie = (l - 1) / nelems + 1
                endif

                do j = 1,ngly
                    do i = 1,nglx
                        ic = (ecol - 1) * ngly * nglx + (j - 1) * nglx + i
                        do k = 1,nglz
                            iz = intma_1d_dg(k,ie)
                            node_column_table_dg(iz,ic) = intma_dg(i,j,k,l)
                        enddo
                    enddo
                enddo
            enddo
        endif

    end subroutine mod_grid_init_node_column

    !-----------------------------------------
    !>@brief Dump field values to file
    !-----------------------------------------
    subroutine dump_field(rhs, name, ex)

        use mod_mpi_utilities, only: irank
        
        implicit none

        real, dimension(:,:), intent(in) :: rhs
        character(len=*), intent(in) :: name
        character(len=24) :: namer
        character(len=4) :: rankn
        logical, intent(in) :: ex
        character*1, parameter:: tab = ' '

        integer:: j, iloop
        
        !build name
        write(unit=rankn,fmt='(I4)') irank
        if(irank == 0) then
            iloop = 3
        else
            iloop = 3 - int(log10(real(irank)))
        endif
        do j=1,iloop
          rankn(j:j)='0'
        end do
        namer=name//rankn
        
        !write coordinates and rhs to file
        open(unit=2,file=namer)

        do j=1,npoin
            write(2,'(es16.9,a1,es16.9,a1,es16.9,a1,es24.16,a1,es24.16,a1,es24.16,a1,es24.16,a1,es24.16)') &
                coord(1,j), tab, coord(2,j), tab, coord(3,j), tab, &
                rhs(1,j), tab,  rhs(2,j), tab, rhs(3,j), tab, rhs(4,j), tab, rhs(5,j)
        end do

        close(2)

        if(ex) then
            call exit(0)
        endif

    end subroutine dump_field

    !---------------------------------------------------------
    !>@brief Rotates the grid in xy plane around the origin
    !---------------------------------------------------------
    subroutine mod_grid_rotate(beta)

      use mod_constants, only: pi

      implicit none

      real beta, x,y
      integer ip
      

      if(beta>0) then
         do ip = 1,npoin
            x = coord(1,ip)
            y = coord(3,ip) 
               
            coord(1,ip) = x * cos(beta) + y * sin (beta)
            coord(3,ip) = -x * sin(beta) + y * cos(beta)
            
         end do
      end if

      end subroutine mod_grid_rotate

    !---------------------------------------------------------
    !>@brief Deforms the grid in z direction based on prescribed function
    !---------------------------------------------------------
    subroutine mod_grid_deform()

      use mod_input, only: ztop
      use mod_constants, only: pi

      implicit none

      real x,y,z,f,h,L,A
      integer ip
      
      L = 100
      A = 1.0
      h = 15

      do ip = 1,npoin
         x = coord(1,ip)                
         y = coord(2,ip)                
         z = coord(3,ip)                

         f = -2.0*A/L*x + A
!          f = A*cos(2*pi*x/L) + pi*A**2/L* &
!               (1 - 1.0/(4*cosh(2*pi*h/L)*cosh(2*pi*h/L)) &
!               + 3.0/(4*sinh(2*pi*h/L)*sinh(2*pi*h/L))) * &
!               cos(4*pi*x/L)
         coord(3,ip) = (ztop+f)/ztop*z

      end do

    end subroutine mod_grid_deform

end module mod_grid
