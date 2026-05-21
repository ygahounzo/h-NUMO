!----------------------------------------------------------------------!
!>@brief This module builds the Grid Geometry: INTMA, COORD, etc.
!>@author  Francis X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>@date FX Giraldo to add sigma_g          on Jun 2014
!>
!>@date S. Marras  to add bdy_flg_g        on Jun 2014
!>
!----------------------------------------------------------------------!
module mod_global_grid

    use mod_basis,     only: basis
    use mod_input,     only: input
    use mod_constants, only: earth_radius, tol, gravity
    use mod_types,     only: r8

    implicit none
    private

    type, public :: grid_global
        real,    dimension(:,:),       allocatable :: coord_g_cg, coord_g
        real,    dimension(:),         allocatable :: sigma_g
        integer, dimension(:,:),       allocatable :: index_g, face_g
        integer, dimension(:,:,:,:),   allocatable :: intma_g, intma_lev_g
        integer, dimension(:,:,:),     allocatable :: intma_s
        integer, dimension(:,:),       allocatable :: bsido_g
        integer, dimension(:),         allocatable :: ele_col_g
        integer, dimension(:),         allocatable :: iperiodic_g
        logical, dimension(:),         allocatable :: flag_periodic_g
        real    :: xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0
        logical :: xperiodic = .false., yperiodic = .false., zperiodic = .false.
        integer :: iboundary(6)   = 0
        integer :: nface_int      = 0, nface_g      = 0
        integer :: npoin_g_cg     = 0, npoin_g      = 0
        integer :: nelem_g        = 0, nelem_r      = 0
        integer :: npoin_g_cg_q   = 0, npoin_g_q    = 0
        integer :: nboun_g        = 0, nboun_poin_g  = 0, ncol_g = 0
        integer :: nelem_s        = 0, npoin_r      = 0, npoin_s = 0
        integer :: nx             = 0, ny           = 0, nz      = 0
        real,    dimension(:,:,:),       allocatable :: saved_coord
        real,    dimension(:,:),         allocatable :: saved_sigma
        integer, dimension(:,:,:),       allocatable :: saved_index, saved_bsido, saved_face
        integer, dimension(:,:,:,:,:),   allocatable :: saved_intma
    end type grid_global

    public :: &
        mod_global_grid_create,     &
        mod_global_grid_bcs,        &
        mod_global_grid_dimensions, &
        mod_global_grid_destroy,    &
        mod_global_grid_out,        &
        mod_global_grid_init_coord

contains

    !-----------------------------------------------------------------------
    subroutine mod_global_grid_create(gg, inp, b)

        implicit none

        type(grid_global), intent(inout) :: gg
        type(input),       intent(in)    :: inp
        type(basis),       intent(in)    :: b

        call mod_global_grid_cube_create(gg, inp, b)
        call mod_global_grid_init_coord(gg, inp, b)
        call create_face(gg%face_g, gg%nface_g, b%FACE_LEN)

    end subroutine mod_global_grid_create

    !-----------------------------------------------------------------------
    subroutine mod_global_grid_init_coord(gg, inp, b)

        implicit none

        type(grid_global), intent(inout) :: gg
        type(input),       intent(in)    :: inp
        type(basis),       intent(in)    :: b

        integer :: i, j, k, e, ip, ip1
        logical :: is_cgc

        is_cgc = .false.
        if(inp%space_method == 'cgc') is_cgc = .true.

        do i = 1, b%nglx
            do j = 1, b%ngly
                do k = 1, b%nglz
                    do e = 1, gg%nelem_g
                        ip = gg%intma_g(i,j,k,e)
                        if(is_cgc) then
                            ip1 = ip
                        else
                            ip1 = (e-1) * b%nglz * b%ngly * b%nglx + &
                                  (k-1) * b%ngly * b%nglx + &
                                  (j-1) * b%nglx + &
                                  (i-1) &
                                  + 1
                        endif
                        gg%coord_g(:,ip1) = gg%coord_g_cg(:,ip)
                    enddo
                enddo
            enddo
        enddo

    end subroutine mod_global_grid_init_coord

    !-----------------------------------------------------------------------
    subroutine mod_global_grid_cube_create(gg, inp, b)

        implicit none

        type(grid_global), intent(inout) :: gg
        type(input),       intent(in)    :: inp
        type(basis),       intent(in)    :: b

        integer :: AllocateStatus
        integer :: i, j
        real :: xc, yc, zc, ac, hc
        real :: x, y, z, ztop
        real, dimension(:), allocatable :: zsurf

        gg%nx = inp%nelx * b%nopx + 1
        gg%ny = inp%nely * b%nopy + 1
        gg%nz = inp%nelz * b%nopz + 1

        gg%nelem_s = inp%nelx * inp%nely
        gg%npoin_s = gg%nx * gg%ny

        gg%nelem_g    = inp%nelx * inp%nely * inp%nelz
        gg%npoin_g_cg = gg%nx * gg%ny * gg%nz
        gg%npoin_g_cg_q = ((b%nqx-1)*inp%nelx + 1) * &
                          ((b%nqy-1)*inp%nely + 1) * &
                          ((b%nqz-1)*inp%nelz + 1)

        if(inp%space_method == 'cgc') then
            gg%npoin_g = gg%npoin_g_cg
        else
            gg%npoin_g   = gg%nelem_g * (b%nopx + 1) * (b%nopy + 1) * (b%nopz + 1)
            gg%npoin_g_q = gg%nelem_g * b%nqx * b%nqy * b%nqz
        endif
        gg%ncol_g = gg%npoin_g_cg / gg%nz

        gg%nboun_g = 0
        if(b%nopx > 0) gg%nboun_g = gg%nboun_g + 2*inp%nely*inp%nelz
        if(b%nopy > 0) gg%nboun_g = gg%nboun_g + 2*inp%nelx*inp%nelz
        if(b%nopz > 0) gg%nboun_g = gg%nboun_g + 2*inp%nelx*inp%nely

        gg%nboun_poin_g = 8
        if(b%nopx > 0) gg%nboun_poin_g = gg%nboun_poin_g + 2*(gg%ny-2)*(gg%nz-2) + 4*(gg%nx-2)
        if(b%nopy > 0) gg%nboun_poin_g = gg%nboun_poin_g + 2*(gg%nx-2)*(gg%nz-2) + 4*(gg%ny-2)
        if(b%nopz > 0) gg%nboun_poin_g = gg%nboun_poin_g + 2*(gg%nx-2)*(gg%ny-2) + 4*(gg%nz-2)

        gg%nface_g = 0
        if(b%nopx > 0) gg%nface_g = gg%nface_g + (inp%nelx + 1)*inp%nely*inp%nelz
        if(b%nopy > 0) gg%nface_g = gg%nface_g + (inp%nely + 1)*inp%nelx*inp%nelz
        if(b%nopz > 0) gg%nface_g = gg%nface_g + (inp%nelz + 1)*inp%nelx*inp%nely

        gg%iboundary(1:2) = inp%z_boundary(1:2)
        gg%iboundary(3:4) = inp%y_boundary(1:2)
        gg%iboundary(5:6) = inp%x_boundary(1:2)

        gg%xmin = inp%xdims(1)  ;  gg%xmax = inp%xdims(2)
        gg%ymin = inp%ydims(1)  ;  gg%ymax = inp%ydims(2)
        gg%zmin = inp%zbottom   ;  gg%zmax = inp%ztop

        gg%xperiodic = .false.
        gg%yperiodic = .false.
        gg%zperiodic = .false.

        if (gg%iboundary(1) == 3 .and. gg%iboundary(2) == 3) gg%zperiodic = .true.
        if (gg%iboundary(3) == 3 .and. gg%iboundary(4) == 3) gg%yperiodic = .true.
        if (gg%iboundary(5) == 3 .and. gg%iboundary(6) == 3) gg%xperiodic = .true.

        print *, ' Boundary Conditions are: '
        do i = 1, 6
            print *, ' i, iboundary = ', i, gg%iboundary(i)
        end do
        print *, ' xperiodic = ', gg%xperiodic
        print *, ' yperiodic = ', gg%yperiodic
        print *, ' zperiodic = ', gg%zperiodic

        allocate(gg%coord_g_cg(3, gg%npoin_g_cg), gg%coord_g(3, gg%npoin_g), &
                 gg%face_g(b%FACE_LEN, gg%nface_g),                           &
                 gg%sigma_g(gg%npoin_g_cg), gg%index_g(2, gg%npoin_g_cg),    &
                 gg%intma_g(b%nglx, b%ngly, b%nglz, gg%nelem_g),              &
                 gg%bsido_g(6, gg%nboun_g),                                    &
                 zsurf(gg%npoin_g_cg),                                         &
                 gg%ele_col_g(gg%nelem_g),                                     &
                 stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Grid **"

        gg%coord_g_cg = 0.0
        gg%coord_g    = 0.0
        gg%sigma_g    = 0.0
        gg%index_g    = 0
        gg%intma_g    = 0
        gg%bsido_g    = 0
        zsurf         = 0.0
        gg%ele_col_g  = 0

        call create_grid_cube(inp, gg%coord_g_cg, gg%index_g, gg%intma_g, gg%ele_col_g, gg%bsido_g, &
            gg%npoin_g_cg, gg%npoin_g, gg%nelem_g, gg%nboun_g,                                  &
            b%xglx, b%xgly, b%xglz, b%nglx, b%ngly, b%nglz,                                    &
            inp%nelx, inp%nely, inp%nelz,                                                         &
            gg%xmin, gg%xmax, gg%ymin, gg%ymax, gg%zmin, gg%zmax,                              &
            gg%nx, gg%ny, gg%nz,                                                                  &
            gg%iboundary, gg%xperiodic, gg%yperiodic, gg%zperiodic)

    end subroutine mod_global_grid_cube_create

    !-----------------------------------------------------------------------
    subroutine mod_global_grid_bcs(gg)

        use mpi

        implicit none

        type(grid_global), intent(inout) :: gg

        integer :: ierr

        call mpi_bcast(gg%xperiodic, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(gg%yperiodic, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(gg%zperiodic, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(gg%iboundary, 6, mpi_integer, 0, mpi_comm_world, ierr)

    end subroutine mod_global_grid_bcs

    !-----------------------------------------------------------------------
    subroutine mod_global_grid_dimensions(gg)

        use mpi
        use mod_mpi_utilities, only: MPI_PRECISION

        implicit none

        type(grid_global), intent(inout) :: gg

        integer :: ierr

        call mpi_bcast(gg%xmin, 1, MPI_PRECISION, 0, mpi_comm_world, ierr)
        call mpi_bcast(gg%xmax, 1, MPI_PRECISION, 0, mpi_comm_world, ierr)
        call mpi_bcast(gg%ymin, 1, MPI_PRECISION, 0, mpi_comm_world, ierr)
        call mpi_bcast(gg%ymax, 1, MPI_PRECISION, 0, mpi_comm_world, ierr)
        call mpi_bcast(gg%zmin, 1, MPI_PRECISION, 0, mpi_comm_world, ierr)
        call mpi_bcast(gg%zmax, 1, MPI_PRECISION, 0, mpi_comm_world, ierr)

    end subroutine mod_global_grid_dimensions

    !-----------------------------------------------------------------------
    subroutine mod_global_grid_destroy(gg)

        implicit none

        type(grid_global), intent(inout) :: gg

        deallocate(gg%coord_g_cg, gg%coord_g, gg%intma_g, gg%bsido_g)

    end subroutine mod_global_grid_destroy

    !-----------------------------------------------------------------------
    subroutine mod_global_grid_out(gg, b)

        implicit none

        type(grid_global), intent(in) :: gg
        type(basis),       intent(in) :: b

        integer :: ivar, ip, ie, k, mk, j, mj, i, mi
        character :: text*72
        character :: fname*72

        fname = 'globalgrid.gri'
        open(1, file=fname)

        text = " NPOIN  NELEM NBOUN NOPX NOPY NOPZ "
        write(1,'(a)') text
        write(1,'(6(i7,1x))') gg%npoin_g_cg, gg%nelem_g, gg%nboun_g, b%nopx, b%nopy, b%nopz

        text = " COORDINATES "
        write(1,'(a)') text
        do ip = 1, gg%npoin_g_cg
            write(1, *) ip, gg%coord_g_cg(1,ip), gg%coord_g_cg(2,ip), gg%coord_g_cg(3,ip)
        end do

        text = " ELEMENT CONNECTIVITY "
        write(1,'(a)') text
        do ie = 1, gg%nelem_g
            do k = 1, b%nglz
                mk = (k - 1)*b%nglx*b%ngly
                do j = 1, b%ngly
                    mj = (j-1)*b%nglx
                    do i = 1, b%nglx
                        mi = mk + mj + i
                        write(1,'(3(i7,1x))') ie, mi, gg%intma_g(i,j,k,ie)
                    end do
                end do
            end do
        end do

        close(1)

    end subroutine mod_global_grid_out

end module mod_global_grid
