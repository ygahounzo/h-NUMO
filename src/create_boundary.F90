!----------------------------------------------------------------------!
!>@brief This subroutine creates the boundary of each processor
!>@author James F. Kelly 19 May 2010 
!>@detail
!> ip_bound = list of (local) boundary points
!> num_bound = number of boundary points
!> btype = BC type (-7 = bounding another processor)s
!----------------------------------------------------------------------!
subroutine create_boundary_cg(ip_bound,num_bound,btype,ndim)
  
    use mod_basis, only: ngl, nglx, ngly, nglz
    use mod_face, only: imapl
    use mod_grid, only: intma_cg, nelem, nboun, nface, face, mod_grid_get_face_ngl
  
    implicit none

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

    do iface = 1, nface
        if ( face(8,iface) == -btype ) then
            iloc = face(5,iface)
            iel = face(7,iface)

            call mod_grid_get_face_ngl(iloc, ngl_i,ngl_j, plane_ij)

            do i1 = 1,ngl_i
                do j1 = 1,ngl_j
                    !Map Face Coords to Grid Points
                    i = imapl(1,i1,j1,iface)
                    j = imapl(2,i1,j1,iface)
                    k = imapl(3,i1,j1,iface)

                    !Get Grid point
                    ip = intma_cg(i,j,k,iel)
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
  
end subroutine
