!------------------------------------------------------------------------------
!>@brief This subroutine computes the Courant Number
!>@author  Yao Gahounzo
!>           PhD Computing
!>           Boise State University
!>           Date: June 22, 2025
!------------------------------------------------------------------------------

subroutine courant_mlswe(G, b, cfl_vector, q_layers, qb, dt, dt_btp, nlayers, nvarb, min_dx_vec)

    use mod_grid,  only: grid
    use mod_basis, only: basis

    implicit none

    type(grid),  intent(in) :: G
    type(basis), intent(in) :: b

    !global arrays
    real,    intent(in)  :: q_layers(5,G%npoin,nlayers), qb(nvarb,G%npoin)
    integer, intent(in)  :: nlayers, nvarb
    real,    intent(in)  :: dt, dt_btp
    real,    intent(out) :: cfl_vector(2), min_dx_vec(2)

    !local arrays
    real :: min_dx, min_dy, cfl_b, cfl

    call courant_cube_mlswe(G, b, cfl, cfl_b, q_layers, qb, dt, dt_btp, nlayers, nvarb, min_dx, min_dy)

    cfl_vector(1) = cfl_b
    cfl_vector(2) = cfl

    min_dx_vec(1) = min_dx
    min_dx_vec(2) = min_dy

end subroutine courant_mlswe

subroutine courant_cube_mlswe(G, b, cfl, cfl_b, q_layers, qb, dt, dt_btp, nlayers, nvarb, min_dx, min_dy)

    use mod_grid,  only: grid
    use mod_basis, only: basis

    implicit none

    type(grid),  intent(in) :: G
    type(basis), intent(in) :: b

    !global arrays
    real,    intent(in)  :: q_layers(5,G%npoin,nlayers), qb(nvarb,G%npoin)
    integer, intent(in)  :: nlayers, nvarb
    real,    intent(in)  :: dt, dt_btp
    real,    intent(out) :: cfl, cfl_b, min_dx, min_dy

    !local arrays
    real    :: x(4), y(4)
    integer :: inode(4)
    integer :: ie, i, j, k, m, ii, jj, kk, il
    integer :: npoints_per_cell
    real :: dx, dy, ub, vb, uk, vk

    npoints_per_cell = 4 !if always 8, then no issue

    !initialize
    cfl   = -1.0d10
    cfl_b = cfl

    min_dx = 1e16
    min_dy = 1e16

    do ie=1,G%nelem
        do j=1,max(b%ngly-1,1)
            do i=1,max(b%nglx-1,1)

                ii = min(i+1,b%nglx)
                jj = min(j+1,b%ngly)
                kk = 1
                inode(1) = G%intma( i, j,1,ie)
                inode(2) = G%intma(ii, j,1,ie)
                inode(3) = G%intma( i,jj,1,ie)
                inode(4) = G%intma(ii,jj,1,ie)

                ub = 0.0; vb = 0.0

                ! Barotropic
                do m=1,npoints_per_cell

                    x(m) = G%coord(1,inode(m))
                    y(m) = G%coord(2,inode(m))

                    ub = ub + qb(3,inode(m)) / npoints_per_cell
                    vb = vb + qb(4,inode(m)) / npoints_per_cell

                end do !m

                dx = maxval(x(:)) - minval(x(:))
                dy = maxval(y(:)) - minval(y(:))

                min_dx = min(min_dx,dx)
                min_dy = min(min_dy,dy)

                cfl_b = max(cfl_b, abs(ub)*dt_btp/min_dx)
                cfl_b = max(cfl_b, abs(vb)*dt_btp/min_dy)

                ! Baroclinic
                do k = 1,nlayers

                    uk = 0.0; vk = 0.0

                    do m = 1,npoints_per_cell

                        x(m) = G%coord(1,inode(m))
                        y(m) = G%coord(2,inode(m))

                        uk = uk + q_layers(2,inode(m),k) / npoints_per_cell
                        vk = vk + q_layers(3,inode(m),k) / npoints_per_cell

                    end do !m

                    cfl = max(abs(uk)*dt/min_dx, cfl)
                    cfl = max(abs(vk)*dt/min_dy, cfl)

                end do !k
            end do !i
        end do !j
    end do !ie

end subroutine courant_cube_mlswe
