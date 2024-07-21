!----------------------------------------------------------------------!
!>@brief This subroutine computes the Courant Number
!>@author  Francis X. Giraldo on 7/08
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!----------------------------------------------------------------------!

subroutine courant_mlswe(cfl_vector,q_layers,qb,dt,dt_btp,nlayers,min_dx_vec)

    use mod_grid, only: npoin

    implicit none

    !global arrays
    real ::  q_layers(5,npoin,nlayers), qb(4,npoin)
    integer :: nlayers
    real ::  cfl_vector(5), dt, min_dx_vec(2), dt_btp

    !local arrays
    real ::  cfl_h, cfl_v, cflu_h, cflu_v, min_dx,min_dy
    real :: cfl_b, cflu_b, cfl

    call courant_cube_mlswe(cfl,cfl_h,cfl_b,cflu_h,cflu_b,q_layers,qb,dt,dt_btp,nlayers,min_dx,min_dy)

    cfl_vector(1)=cfl_h
    cfl_vector(2)=cfl_b
    cfl_vector(3)=cflu_h
    cfl_vector(4)=cflu_b
    cfl_vector(5)=cfl

    min_dx_vec(1) = min_dx
    min_dx_vec(2) = min_dy

end subroutine courant_mlswe

subroutine courant_cube_mlswe(cfl,cfl_h,cfl_b,cflu_h,cflu_b,q_layers,qb,dt,dt_btp,nlayers,min_dx,min_dy)

    use mod_basis, only: nglx, ngly, nglz

    use mod_constants, only: gravity

    use mod_grid, only: coord, intma, npoin, nelem

    use mod_constants, only: gravity

    use mod_initial, only: alpha_mlswe

    implicit none

    !global arrays
    real    ::  q_layers(5,npoin,nlayers), qb(4,npoin)
    real    ::  dt, cfl, cfl_h, cflu_h, cfl_b, cflu_b, dt_btp
    real    ::  min_dx, min_dy
    integer ::  nlayers
    
    !local arrays
    real    ::  x(8), y(8), z(8)
    integer ::  inode(8)
    real    ::  c, cflx, cfly, cflz, cflxu, cflyv, cflzw
    real    ::  r, rho, u, v, w, theta, dx, dy, dz, ds, vel, E, p, a, htot
    real    ::  rho0, u0, v0, w0, theta0, rho_tot
    real    ::  sig, sigu, sigx, sigy, sigz, sigxu, sigyv, sigzw
    integer ::  ie, i, j, k, m, ii, jj, kk, il,ill
    real    ::  dxc, dyc, dzc
    integer ::  npoints_per_cell
    real :: u1,u2,v1,v2,ub,vb,hb,h1,h2,r1,r2,vel1,vel2,vel_b,sig1,sig2,sigb,sig1uv,sig2uv,lam,lamb,c_maxB
    real :: gprime, rb,cfl1,cfl2,cfl1uv,cfl2uv
  
    npoints_per_cell=4 !FXG: if always 8, then no issue

    !initialize
    cfl=-1.0d10

    cfl1=cfl
    cfl2=cfl
    cfl_b=cfl
    cfl1uv=cfl
    cfl2uv=cfl
    cflu_b=cfl
    c_maxB = cfl


    min_dx = 1e16
    min_dy = 1e16

    gprime = gravity *alpha_mlswe(1)*(1.0/alpha_mlswe(2) - 1.0/alpha_mlswe(1))

    !Simplified Old method
    do ie=1,nelem
        do k=1,max(nglz-1,1)
            do j=1,max(ngly-1,1)
                do i=1,max(nglx-1,1)
                    ii=min(i+1,nglx)
                    jj=min(j+1,ngly)
                    kk=min(k+1,nglz)
                    inode(1)=intma( i, j, k,ie)
                    inode(2)=intma(ii, j, k,ie)
                    inode(3)=intma( i,jj, k,ie)
                    inode(4)=intma(ii,jj, k,ie)
                    inode(5)=intma( i, j,kk,ie)
                    inode(6)=intma(ii, j,kk,ie)
                    inode(7)=intma( i,jj,kk,ie)
                    inode(8)=intma(ii,jj,kk,ie)

                    u1=0.0; v1=0.0; u2=0.0; v2=0.0; ub=0.0; vb=0.0; hb=0.0; h1=0.0; h2=0.0; r1=0.0; r2=0.0
                    do m=1,npoints_per_cell
                        x(m)=coord(1,inode(m))
                        y(m)=coord(2,inode(m))
                        z(m)=coord(3,inode(m))

                        ! r1 = (alpha_mlswe(1)/gravity)*q_layers(1,inode(m),1)
                        ! r2 = (alpha_mlswe(2)/gravity)*q_layers(1,inode(m),2)

                        r1 = q_layers(1,inode(m),1)
                        r2 = q_layers(1,inode(m),2)

                        rb = (sum(alpha_mlswe(:))/gravity)*qb(1,inode(m))
                        
                        ! u1 = u1 + (q_layers(2,inode(m),1) / q_layers(1,inode(m),1)) / npoints_per_cell
                        ! v1 = v1 + (q_layers(3,inode(m),1) / q_layers(1,inode(m),1)) / npoints_per_cell
                        ! u2 = u2 + (q_layers(2,inode(m),2) / q_layers(1,inode(m),2)) / npoints_per_cell
                        ! v2 = v2 + (q_layers(3,inode(m),2) / q_layers(1,inode(m),2)) / npoints_per_cell

                        u1 = u1 + q_layers(2,inode(m),1) / npoints_per_cell
                        v1 = v1 + q_layers(3,inode(m),1) / npoints_per_cell
                        u2 = u2 + q_layers(2,inode(m),2) / npoints_per_cell
                        v2 = v2 + q_layers(3,inode(m),2) / npoints_per_cell

                        ub = ub + qb(3,inode(m)) / npoints_per_cell
                        vb = vb + qb(4,inode(m)) / npoints_per_cell

                        ! ub = ub + (qb(3,inode(m)) / qb(1,inode(m))) / npoints_per_cell
                        ! vb = vb + (qb(4,inode(m)) / qb(1,inode(m))) / npoints_per_cell

                        hb = hb + rb/npoints_per_cell
                        h1 = h1 + r1/npoints_per_cell
                        h2 = h2 + r2/npoints_per_cell
                        
                    end do !m

                    dx=maxval(x(:))-minval(x(:))
                    dy=maxval(y(:))-minval(y(:))

                    min_dx = min(min_dx,dx)
                    min_dy = min(min_dy,dy)
                    ds = sqrt(dx*dx + dy*dy)

                    lamb = sqrt(gravity*hb)
                    lam = sqrt(gprime*(h1*h2)/(h1+h2))

                    ! Barotropic 

                    vel_b = sqrt(ub*ub + vb*vb)
                    cflu_b = dt_btp*vel_b/ds

                    c_maxB = max(c_maxB,vel_b + lamb)
                    cfl_b = dt_btp*c_maxB/ds

                    sig = dt*c_maxB/ds
                    cfl = max(cfl,sig)

                    ! Baroclinic

                    vel1 = sqrt(u1*u1 + v1*v1)
                    sig1 = dt*(vel1+lam)/ds
                    sig1uv = dt*vel1/ds

                    vel2 = sqrt(u2*u2 + v2*v2)
                    sig2 = dt*(vel2+lam)/ds
                    sig2uv = dt*vel2/ds

                    cfl1 = max(cfl1,sig1)
                    cfl2 = max(cfl2,sig2)
                    cfl1uv = max(cfl1uv,sig1uv)
                    cfl2uv = max(cfl2uv,sig2uv)

                end do !i
            end do !j
        end do !k
    end do !ie

    cfl_h = max(cfl1,cfl2)
    cflu_h = max(cfl1uv,cfl2uv)

end subroutine courant_cube_mlswe