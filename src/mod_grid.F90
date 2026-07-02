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

    use mod_basis
    use mod_input
    use iso_c_binding, only: C_INT

    implicit none

    !-----------------------------------------------------------------------
    ! grid type defined first so abstract interfaces below can reference it.
    !-----------------------------------------------------------------------
    type :: grid
        real,    dimension(:,:),      allocatable :: coord_cg, coord
        real,    dimension(:),        allocatable :: sigma
        integer, dimension(:,:),      allocatable :: index2d
        integer, dimension(:,:,:,:),  allocatable :: intma_table
        integer, dimension(:,:,:,:),  allocatable :: intma
        integer, dimension(:,:,:,:),  allocatable :: intma_dg_quad  !< precomputed DG quad-point index
        integer, dimension(:,:,:,:),  allocatable :: D2C_mask
        integer(C_INT), dimension(:,:), allocatable :: NC_face
        integer(C_INT), dimension(:,:), allocatable :: NC_edge
        integer(C_INT), dimension(:),   allocatable :: EToNC
        integer(C_INT)                              :: nNC = 0
        integer, dimension(:,:),      allocatable :: bsido
        integer, dimension(:,:),      allocatable :: node_column_table, node_column_table_dg
        integer, dimension(:),        allocatable :: face_type, is_non_conforming, coupling_flg
        integer, dimension(:,:),      allocatable :: face
        integer :: npoin_cg = 0, npoin_massinv = 0, npoin = 0, nelem = 0, nbsido = 0
        integer :: ncol = 0, ncol_cg = 0, nz = 0, nz_cg = 0
        integer :: nface = 0, nboun = 0
        integer :: nelemx = 0, nelemy = 0, nelemz = 0
        integer :: npoin_q = 0
        logical :: do_1dIMEX = .false.
        logical :: is_sphere  = .false.
        integer, dimension(12)  :: edge_limit = 0
        integer, dimension(2,6) :: face_limit = 0
    end type grid
    !-----------------------------------------------------------------------

    !-----------------------------------------
    !>@brief Interfaces for procedure pointers.
    !> G and b are explicit arguments so target functions need no module-level state.
    !-----------------------------------------
    abstract interface
        pure function intma_interface(G, i,j,k,e) result(r)
            import :: grid
            type(grid), intent(in) :: G
            integer, intent(in) :: i, j, k, e
            integer :: r
        end function intma_interface
        function intma_1d_interface(b, k, e) result(r)
            import :: basis
            type(basis), intent(in) :: b
            integer, intent(in) :: k, e
            integer :: r
        end function intma_1d_interface
        function node_column_interface(G, ic, iz) result(r)
            import :: grid
            type(grid), intent(in) :: G
            integer, intent(in) :: ic, iz
            integer :: r
        end function node_column_interface
        function intma_dg_to_cg_interface(G, b, index) result(r)
            import :: grid, basis
            type(grid), intent(in) :: G
            type(basis), intent(in) :: b
            integer, intent(in) :: index
            integer :: r
        end function intma_dg_to_cg_interface
    end interface

    integer, private :: nh_i, nz_i

    ! Procedure pointers — set by mod_grid_init_unified based on space_method.
    procedure(intma_1d_interface),       pointer :: intma_1d    
    procedure(node_column_interface),    pointer :: node_column
    procedure(intma_dg_to_cg_interface), pointer :: intma_dg_to_cg

    public :: &
        mod_grid_init_unified,       &
        mod_grid_init_coord,         &
        mod_grid_init_coord_dg_to_cg,&
        mod_grid_get_face_ngl,       &
        mod_grid_get_face_nq,        &
        grid,                        &
        intma_cg,                    &
        intma_dg_to_cg,              &
        intma_1d, intma_1d_cg,       &
        node_column, node_column_cg, &
        is_cg_coupled

    private

contains

    ! ── Indexing functions ────────────────────────────────────────────────
    ! All take explicit G and/or b — no module-level state accessed.

    pure function intma_cg(G, i,j,k,e)
        type(grid), intent(in) :: G
        integer, intent(in) :: i, j, k, e
        integer :: intma_cg
        intma_cg = G%intma_table(i,j,k,e)
    end function intma_cg

    pure function intma_dg(b, i,j,k,e)
        type(basis), intent(in) :: b
        integer, intent(in) :: i, j, k, e
        integer :: intma_dg
        intma_dg = (e-1)*b%nglz*b%ngly*b%nglx + &
                   (k-1)*b%ngly*b%nglx         + &
                   (j-1)*b%nglx                + &
                   (i-1)                        + &
                   1
    end function intma_dg

    !------------------------------------------
    !>@brief intma for 1d-CG indexing
    !------------------------------------------
    function intma_1d_cg(b, k, e)
        type(basis), intent(in) :: b
        integer, intent(in) :: k, e
        integer :: intma_1d_cg
        intma_1d_cg = (e - 1) * (b%nglz - 1) + k
    end function intma_1d_cg

    !------------------------------------------
    !>@brief intma for 1d-DG indexing
    !------------------------------------------
    function intma_1d_dg(b, k, e)
        type(basis), intent(in) :: b
        integer, intent(in) :: k, e
        integer :: intma_1d_dg
        intma_1d_dg = (e - 1) * b%nglz + k
    end function intma_1d_dg

    !------------------------------------------
    !>@brief node_column for CG 
    !------------------------------------------
    function node_column_cg(G, ic, iz)
        type(grid), intent(in) :: G
        integer, intent(in) :: ic, iz
        integer :: node_column_cg
        node_column_cg = G%node_column_table(iz,ic)
    end function node_column_cg

    !------------------------------------------
    !>@brief node_column for DG 
    !------------------------------------------
    function node_column_dg(G, ic, iz)
        type(grid), intent(in) :: G
        integer, intent(in) :: ic, iz
        integer :: node_column_dg
        node_column_dg = G%node_column_table_dg(iz,ic)
    end function node_column_dg

    !-----------------------------------
    !>@brief DG to CG — identity for CG0  
    !-----------------------------------
    function intma_dg_to_cg_1(G, b, index)
        type(grid), intent(in) :: G
        type(basis), intent(in) :: b
        integer, intent(in) :: index
        integer :: intma_dg_to_cg_1
        intma_dg_to_cg_1 = index
    end function intma_dg_to_cg_1

    !-----------------------------------
    !>@brief DG to CG — full decomposition  
    !-----------------------------------
    function intma_dg_to_cg_2(G, b, index)
        type(grid), intent(in) :: G
        type(basis), intent(in) :: b
        integer, intent(in) :: index
        integer :: intma_dg_to_cg_2
        integer :: i, j, k, e, rem
        rem = index - 1
        e   = rem / (b%nglz * b%ngly * b%nglx) + 1
        rem = mod(rem, b%nglz * b%ngly * b%nglx)
        k   = rem / (b%ngly * b%nglx) + 1
        rem = mod(rem, b%ngly * b%nglx)
        j   = rem / b%nglx + 1
        i   = mod(rem, b%nglx) + 1
        intma_dg_to_cg_2 = intma_cg(G, i,j,k,e)
    end function intma_dg_to_cg_2

    ! ── Initialisation subroutines ────────────────────────────────────────

    !---------------------------------------------------------
    !>@brief Initialize unified CG/DG indexing and allocate grid arrays.
    !---------------------------------------------------------
    subroutine mod_grid_init_unified(G, inp, b)

        implicit none

        type(grid),  intent(inout) :: G
        type(input), intent(in)    :: inp
        type(basis), intent(in)    :: b

        integer :: AllocateStatus, i, j, k, e
        integer :: npoin_dg

        if (allocated(G%intma)) then
            !!$acc exit data delete(G%intma)
            deallocate(G%intma)
        endif

        allocate(G%intma(b%nglx, b%ngly, b%nglz, G%nelem))

        if (inp%space_method == 'cgc') then
            !3d
            G%npoin   = G%npoin_cg
            npoin_dg  = G%nelem * (b%nopx + 1) * (b%nopy + 1) * (b%nopz + 1)
            do concurrent (i=1:b%nglx, j=1:b%ngly, k=1:b%nglz, e=1:G%nelem)
                G%intma(i,j,k,e) = intma_cg(G, i,j,k,e)
            end do
            intma_dg_to_cg => intma_dg_to_cg_1

            !1D
            G%nz_cg   = inp%nelz * b%nopz + 1
            G%nz      = G%nz_cg
            G%ncol_cg = G%npoin_cg / G%nz_cg
            G%ncol    = G%ncol_cg
            intma_1d  => intma_1d_cg

            !1D in 3D
            node_column => node_column_cg
        else
            !3D
            G%npoin   = G%nelem * (b%nopx + 1) * (b%nopy + 1) * (b%nopz + 1)
            G%npoin_q = G%nelem * b%nqx * b%nqy * b%nqz
            do concurrent (i=1:b%nglx, j=1:b%ngly, k=1:b%nglz, e=1:G%nelem)
                G%intma(i,j,k,e) = intma_dg(b, i,j,k,e)
            end do
            intma_dg_to_cg => intma_dg_to_cg_2

            !1D
            G%nz_cg   = inp%nelz * b%nopz + 1
            G%nz      = inp%nelz * (b%nopz + 1)
            G%ncol_cg = G%npoin_cg / G%nz_cg
            G%ncol    = G%npoin / G%nz
            intma_1d  => intma_1d_dg

            !1D in 3D
            node_column => node_column_dg

            if (inp%is_non_conforming_flg > 0 .and. inp%space_method == 'cgd') then
                G%npoin_massinv = G%npoin_cg
            else
                G%npoin_massinv = G%npoin
            endif

            ! Precompute DG quadrature-point global indices
            if (allocated(G%intma_dg_quad)) deallocate(G%intma_dg_quad)
            allocate(G%intma_dg_quad(b%nqx, b%nqy, b%nqz, G%nelem))
            do concurrent (i=1:b%nqx, j=1:b%nqy, k=1:b%nqz, e=1:G%nelem)
                G%intma_dg_quad(i,j,k,e) = &
                    (e-1)*b%nqz*b%nqy*b%nqx + (k-1)*b%nqy*b%nqx + (j-1)*b%nqx + (i-1) + 1
            end do
        endif

        G%edge_limit = (/b%nglx, b%nglx, b%nglx, b%nglx, &
                         b%ngly, b%ngly, b%ngly, b%ngly,  &
                         b%nglz, b%nglz, b%nglz, b%nglz/)
        G%face_limit = reshape((/ b%ngly, b%nglz, b%ngly, b%nglz, &
                                  b%nglx, b%nglz, b%nglx, b%nglz, &
                                  b%nglx, b%ngly, b%nglx, b%ngly /), &
                               (/2, 6/))

        !-----------------------
        ! Allocate arrays
        !-----------------------
        if (allocated(G%coord_cg)) then
            deallocate(G%coord_cg, G%coord, G%sigma,   &
                G%index2d,                              &
                G%intma_table,                          &
                G%D2C_mask,                             &
                G%NC_face,                              &
                G%NC_edge,                              &
                G%EToNC,                                &
                G%face,                                 &
                G%face_type,                            &
                G%bsido,                                &
                G%coupling_flg,                         &
                G%is_non_conforming,                    &
                stat=AllocateStatus )
            if (AllocateStatus /= 0) stop "** Deallocate failed - Mod_Grid **"
        endif
        allocate(G%coord_cg(3, G%npoin_cg), G%coord(3, G%npoin),    &
            G%sigma(G%npoin_cg),                                     &
            G%index2d(2, G%npoin_cg),                                &
            G%intma_table(b%nglx, b%ngly, b%nglz, G%nelem),         &
            G%D2C_mask(b%nglx, b%ngly, b%nglz, G%nelem),            &
            G%NC_face(b%P4EST_FACES, G%nNC),                         &
            G%NC_edge(b%P8EST_EDGES, G%nNC),                         &
            G%EToNC(G%nelem),                                        &
            G%face(b%FACE_LEN, G%nface),                             &
            G%face_type(G%nface),                                    &
            G%bsido(6, G%nbsido),                                    &
            G%coupling_flg(G%npoin_cg),                              &
            G%is_non_conforming(G%nface),                            &
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Grid **"

        G%EToNC       = 0
        G%coord_cg    = 0.0
        G%sigma       = 0.0
        G%index2d     = 0
        G%intma_table = 0
        G%D2C_mask    = 0

        call initialize_coupling(G)

    end subroutine mod_grid_init_unified

    !-------------------------------------------------
    !>@brief Initialize coord which could use CG/DG storage
    !-------------------------------------------------
    subroutine mod_grid_init_coord(G, b)

        implicit none

        type(grid),  intent(inout) :: G
        type(basis), intent(in)    :: b

        integer :: i, j, k, e, ip, ip1

        do e = 1, G%nelem
            do k = 1, b%nglz
                do j = 1, b%ngly
                    do i = 1, b%nglx
                        ip  = G%intma_table(i,j,k,e)
                        ip1 = G%intma(i,j,k,e)
                        G%coord(:,ip1) = G%coord_cg(:,ip)
                    enddo
                enddo
            enddo
        enddo
        call initialize_coupling(G)

    end subroutine mod_grid_init_coord

    !-------------------------------------------------
    !>@brief Initialize coord_cg
    !-------------------------------------------------
    subroutine mod_grid_init_coord_dg_to_cg(G, inp, b)

        implicit none

        type(grid),  intent(inout) :: G
        type(input), intent(in)    :: inp
        type(basis), intent(in)    :: b

        integer :: i, j, k, e, ip, ip1

        if (inp%space_method == 'cgc') then
            G%coord_cg(:,:) = G%coord(:,:)
        else
            do i = 1, b%nglx
                do j = 1, b%ngly
                    do k = 1, b%nglz
                        do e = 1, G%nelem
                            ip  = G%intma_table(i,j,k,e)
                            ip1 = intma_dg(b, i,j,k,e)
                            G%coord_cg(:,ip) = G%coord(:,ip1)
                        enddo
                    enddo
                enddo
            enddo
        endif
        call initialize_coupling(G)

    end subroutine mod_grid_init_coord_dg_to_cg

    !-------------------------------------------------
    !>@brief Get face ngls
    !-------------------------------------------------
    subroutine mod_grid_get_face_ngl(b, ilocl, ngl_i, ngl_j, plane_ij)
        implicit none
        type(basis), intent(in)  :: b
        integer,     intent(in)  :: ilocl
        integer,     intent(out) :: ngl_i, ngl_j, plane_ij
        if (ilocl == 1 .or. ilocl == 2) then
            ngl_i = b%nglx; ngl_j = b%ngly; plane_ij = 1
        elseif (ilocl == 3 .or. ilocl == 4) then
            ngl_i = b%nglx; ngl_j = b%nglz; plane_ij = 2
        elseif (ilocl == 5 .or. ilocl == 6) then
            ngl_i = b%ngly; ngl_j = b%nglz; plane_ij = 3
        end if
    end subroutine mod_grid_get_face_ngl

    !-------------------------------------------------
    !>@brief Get face nqs
    !-------------------------------------------------
    subroutine mod_grid_get_face_nq(b, ilocl, nq_i, nq_j, plane_ij)
        implicit none
        type(basis), intent(in)  :: b
        integer,     intent(in)  :: ilocl
        integer,     intent(out) :: nq_i, nq_j, plane_ij
        if (ilocl == 1 .or. ilocl == 2) then
            nq_i = b%nqx; nq_j = b%nqy; plane_ij = 1
        elseif (ilocl == 3 .or. ilocl == 4) then
            nq_i = b%nqx; nq_j = b%nqz; plane_ij = 2
        elseif (ilocl == 5 .or. ilocl == 6) then
            nq_i = b%nqy; nq_j = b%nqz; plane_ij = 3
        end if
    end subroutine mod_grid_get_face_nq

    !------------------------------------------
    !>@brief Mixed CG-DG coupling
    !------------------------------------------
    function is_cg_coupled(G, i,j,k,e)
        type(grid), intent(in) :: G
        integer, intent(in) :: i, j, k, e
        integer :: is_cg_coupled, ip
        ip = intma_cg(G, i,j,k,e)
        is_cg_coupled = G%coupling_flg(ip)
    end function is_cg_coupled

    subroutine initialize_coupling(G)
        type(grid), intent(inout) :: G
        integer :: ip
        do ip = 1, G%npoin_cg
            G%coupling_flg(ip) = 1
        end do
    end subroutine initialize_coupling

end module mod_grid
