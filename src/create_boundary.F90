!----------------------------------------------------------------------!
!>@brief This subroutine creates the boundary of each processor
!>@detail
!> ip_bound = list of (local) boundary points
!> num_bound = number of boundary points
!> btype = BC type (-7 = bounding another processor)s
!----------------------------------------------------------------------!
subroutine create_boundary_cg(G, b, mf, ip_bound, num_bound, btype, ndim)

    use mod_grid,  only: grid, mod_grid_get_face_ngl
    use mod_basis, only: basis
    use mod_face,  only: face_CS

    implicit none

    type(grid),    intent(in) :: G
    type(basis),   intent(in) :: b
    type(face_CS), intent(in) :: mf

    !global arrays
    integer, intent(out) :: ip_bound(ndim)
    integer, intent(out) :: num_bound  !In general, num_bound is less than or equal to nboun
    integer, intent(in)  :: btype, ndim

    !local arrays
    integer :: iface, iloc, iel, i1, j1, i, j, k, ip, ib, ipr, ngl_i, ngl_j, plane_ij
    logical :: is_new

    ! This is a local face counter
    num_bound = 0
    ip_bound = 0

    do iface = 1, G%nface
        if ( G%face(8,iface) == -btype ) then
            iloc = G%face(5,iface)
            iel = G%face(7,iface)

            call mod_grid_get_face_ngl(b, iloc, ngl_i, ngl_j, plane_ij)

            do i1 = 1,ngl_i
                do j1 = 1,ngl_j
                    !Map Face Coords to Grid Points
                    i = mf%imapl(1,i1,j1,iface)
                    j = mf%imapl(2,i1,j1,iface)
                    k = mf%imapl(3,i1,j1,iface)

                    !Get Grid point
                    ip = G%intma(i,j,k,iel)
                    is_new = .true.
                    do ib = 1,num_bound
                        ipr = ip_bound(ib)
                        if (ipr == ip) then
                            is_new = .false.
                        end if
                    end do
                    if (is_new) then
                        num_bound = num_bound + 1
                        ip_bound(num_bound) = ip
                    end if
                end do !j1
            end do !i1
        end if !face=-btype
    end do !iface

end subroutine create_boundary_cg
