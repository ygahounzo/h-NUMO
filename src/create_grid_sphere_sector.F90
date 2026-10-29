!-----------------------------------------------------------------!
! Builds a longitude-latitude sector of the sphere (a closed ocean
! basin bounded by two meridians and two latitude circles) as nelx
! x nely quads and writes it to an ABAQUS-format .inp file that
! p4est reads via p4est_connectivity_read_inp() (the
! read_external_grid_flg path, as for sphere_ico).
!
!   lon(1:2) = western/eastern meridian   (degrees, inp%xdims)
!   lat(1:2) = southern/northern latitude (degrees, inp%ydims)
!   bc_x(1:2) = west/east wall,  bc_y(1:2) = south/north wall
!               (inp%x_boundary, inp%y_boundary)
!
! The four sides are written as T3D2 line elements in one ELSET per
! side, named <SIDE>:BC_<TYPE>, which p4est_bc_read_inp maps onto
! face(8,:) = -bc. Only free-slip walls (4) are accepted: the face
! flux and nodal wall routines treat -4 (slip) and -2 (no-slip) as
! walls, and the reader's NO_SLIP code (5) is neither.
!
! Vertices are put on the unit sphere at their exact (lon,lat);
! mod_p4est.F90 radially projects every DOF onto the sphere of
! radius earth_radius afterward, as for sphere_ico. Element edges
! are therefore great-circle arcs between vertices: exact on the
! meridian walls, a chord-sized approximation of the latitude walls.
! Quads are counterclockwise seen from outside (east, then north).
!-----------------------------------------------------------------!
subroutine create_grid_sphere_sector_inp(filename, nelx, nely, lon, lat, bc_x, bc_y)

    use mod_constants, only: deg2rad

    implicit none

    character(len=*), intent(in) :: filename
    integer,          intent(in) :: nelx, nely
    real,             intent(in) :: lon(2), lat(2)
    integer,          intent(in) :: bc_x(2), bc_y(2)

    character(len=5), parameter :: side_name(4) = (/ 'SOUTH', 'NORTH', 'WEST ', 'EAST ' /)

    integer :: bc_side(4), nline_side(4)
    integer :: nvert, nquad, iside, i, j, iline, iline0, fio
    real    :: olon, olat

    if (nelx < 1 .or. nely < 1) error stop "sphere_sector requires nelx >= 1 and nely >= 1"
    if (.not. (lon(2) > lon(1) .and. lon(2) - lon(1) < 360.0)) &
        error stop "sphere_sector requires xdims(1) < xdims(2) < xdims(1) + 360 (longitudes, degrees)"
    if (.not. (lat(1) > -90.0 .and. lat(2) > lat(1) .and. lat(2) < 90.0)) &
        error stop "sphere_sector requires -90 < ydims(1) < ydims(2) < 90 (latitudes, degrees)"

    bc_side = (/ bc_y(1), bc_y(2), bc_x(1), bc_x(2) /)
    if (any(bc_side /= 4)) &
        error stop "sphere_sector supports free-slip walls only: set x_boundary = 4, 4 and y_boundary = 4, 4"

    nline_side = (/ nelx, nelx, nely, nely /)
    nvert = (nelx+1)*(nely+1)
    nquad = nelx*nely

    open(newunit=fio, file=trim(filename), status='replace', action='write')
    write(fio,'(A)') '*Heading'
    write(fio,'(A)') ' sphere_sector_base_mesh.inp'

    write(fio,'(A)') '*Node'
    do j = 0, nely
        olat = (lat(1) + (lat(2) - lat(1))*real(j)/real(nely))*deg2rad
        do i = 0, nelx
            olon = (lon(1) + (lon(2) - lon(1))*real(i)/real(nelx))*deg2rad
            write(fio,'(I0,3(",",1x,E23.15))') vid(i,j), cos(olat)*cos(olon), cos(olat)*sin(olon), sin(olat)
        end do
    end do

    write(fio,'(A)') '*Element, type=CPS4, ELSET=Surface1'
    do j = 0, nely-1
        do i = 0, nelx-1
            write(fio,'(I0,4(",",1x,I0))') 1 + i + j*nelx, vid(i,j), vid(i+1,j), vid(i+1,j+1), vid(i,j+1)
        end do
    end do

    ! Boundary line elements, labelled after the quads, side by side
    ! in the order of side_name.
    write(fio,'(A)') '*Element, type=T3D2, ELSET=Line1'
    iline = nquad
    do i = 0, nelx-1                                             ! south
        iline = iline + 1
        write(fio,'(I0,2(",",1x,I0))') iline, vid(i,0), vid(i+1,0)
    end do
    do i = 0, nelx-1                                             ! north
        iline = iline + 1
        write(fio,'(I0,2(",",1x,I0))') iline, vid(i,nely), vid(i+1,nely)
    end do
    do j = 0, nely-1                                             ! west
        iline = iline + 1
        write(fio,'(I0,2(",",1x,I0))') iline, vid(0,j), vid(0,j+1)
    end do
    do j = 0, nely-1                                             ! east
        iline = iline + 1
        write(fio,'(I0,2(",",1x,I0))') iline, vid(nelx,j), vid(nelx,j+1)
    end do

    iline0 = nquad
    do iside = 1, 4
        write(fio,'(A)') '*ELSET, ELSET='//trim(side_name(iside))//':BC_'//trim(bc_label(bc_side(iside)))
        do i = 1, nline_side(iside), 10
            write(fio,'(I0,9(",",1x,I0))') (iline0 + j, j = i, min(i+9, nline_side(iside)))
        end do
        iline0 = iline0 + nline_side(iside)
    end do

    close(fio)

contains

    ! Global (1-based) id of lattice vertex (i,j), i east, j north.
    integer function vid(i, j)
        integer, intent(in) :: i, j
        vid = 1 + i + j*(nelx+1)
    end function vid

    ! ELSET suffix for a boundary code (inverse of which_bc in p4est.c).
    character(len=16) function bc_label(bc)
        integer, intent(in) :: bc
        select case (bc)
        case (4)
            bc_label = 'FREE_SLIP'
        case default
            bc_label = ''
        end select
    end function bc_label

end subroutine create_grid_sphere_sector_inp
