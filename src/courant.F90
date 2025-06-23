!------------------------------------------------------------------------------
!>@brief This subroutine computes the Courant Number
!>@author  Francis X. Giraldo on 7/08
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@ modified by Yao Gahounzo 
!------------------------------------------------------------------------------

subroutine courant_mlswe(cfl_vector,q_layers,qb,dt,dt_btp,nlayers,min_dx_vec)

    use mod_grid, only: npoin

    implicit none

    !global arrays
    real, intent(in) ::  q_layers(5,npoin,nlayers), qb(4,npoin)
    integer, intent(in) ::  nlayers
    real, intent(in) ::  dt, dt_btp
    real, intent(out) ::  cfl_vector(2), min_dx_vec(2)

    !local arrays
    real ::  min_dx, min_dy, cfl_b, cfl

    call courant_cube_mlswe(cfl,cfl_b,q_layers,qb,dt, dt_btp,nlayers,min_dx,min_dy)

    cfl_vector(1) = cfl_b
    cfl_vector(2) = cfl

    min_dx_vec(1) = min_dx
    min_dx_vec(2) = min_dy

end subroutine courant_mlswe

subroutine courant_cube_mlswe(cfl,cfl_b,q_layers,qb,dt,dt_btp,nlayers,min_dx,min_dy)

    use mod_basis, only: nglx, ngly, nglz

    use mod_constants, only: gravity

    use mod_grid, only: coord, intma, npoin, nelem

    use mod_constants, only: gravity

    use mod_initial, only: alpha_mlswe

    implicit none

    !global arrays
    real, intent(in)    ::  q_layers(5,npoin,nlayers), qb(4,npoin)
    integer, intent(in) ::  nlayers
    real, intent(in)    ::  dt, dt_btp
    real, intent(out)   ::  cfl, cfl_b, min_dx, min_dy

    !local arrays
    real    ::  x(4), y(4)
    integer ::  inode(4)
    integer ::  ie, i, j, k, m, ii, jj, kk, il
    integer ::  npoints_per_cell
    real :: ub,vb, uk, vk

    npoints_per_cell = 4 !if always 8, then no issue

    !initialize
    cfl = -1.0d10
    cfl_b = cfl

    min_dx = 1e16
    min_dy = 1e16
    
    do ie=1,nelem
        do j=1,max(ngly-1,1)
            do i=1,max(nglx-1,1)
                
                ii = min(i+1,nglx)
                jj = min(j+1,ngly)
                kk = 1
                inode(1) = intma( i, j, 1,ie)
                inode(2) = intma(ii, j, 1,ie)
                inode(3) = intma( i,jj, 1,ie)
                inode(4) = intma(ii,jj, 1,ie)

                ub = 0.0; vb = 0.0

                ! Barotropic 
                do m=1,npoints_per_cell

                    x(m) = coord(1,inode(m))
                    y(m) = coord(2,inode(m))

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

                    uk = 0.0 ; vk = 0.0

                    do m = 1,npoints_per_cell

                        x(m) = coord(1,inode(m))
                        y(m) = coord(2,inode(m))

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