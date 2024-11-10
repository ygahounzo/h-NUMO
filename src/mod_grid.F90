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

    use mod_input, only: ti_method, si_dimension, geometry_type, &
        space_method, nelz,         &
        is_non_conforming_flg

      use iso_c_binding, only: C_INT

    public :: &
        mod_grid_init_unified, &
        mod_grid_init_coord, &
        mod_grid_init_coord_dg_to_cg, &
        mod_grid_get_face_ngl, &
        mod_grid_get_face_nq, &
        coord_cg, &
        coord, &
        sigma, &
        index2d, &
        intma_table, &
        D2C_mask, &
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
    integer, dimension(:,:,:,:),  allocatable :: D2C_mask
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
                D2C_mask, &
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
            D2C_mask(nglx,ngly,nglz,nelem), &
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
        D2C_mask = 0

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

        do ip = 1, npoin_cg
            coupling_flg(ip) = 1
        end do

    end subroutine initialize_coupling
    !-----------------------------------------------------------------------

end module mod_grid
