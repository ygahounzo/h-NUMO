!----------------------------------------------------------------------!
!>@brief This subroutine computes size of the element hexa 
!>       and the mean spacing !between the lgl points inside the element.
!>
!>@author Simone Marras
!>@date 11/2013
!>
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>@author Simone Marras
!>@date 12/2014
!>  added the distance computation for elements on the cubed sphere
!>
!> The following local node numbering is assumed (as it is everywhere in NUMA)
!>
!>      7---------- 8
!>    / . Top    /  |
!>   /  .       /   |
!>  5----------6    |
!>  |   .      |    |
!>  |   3------|----4
!>  |  .  Bott |  /
!>  | .        | /
!>  1---------2
!>
!> WARNING: the distances on the sphere assume that nglx=ngly=nglz
!>
!----------------------------------------------------------------------!

subroutine compute_element_size_driver(ie, dx,dy,dz, dx_mean, dy_mean, dz_mean, &
                                       ds_mean, &
                                       Lambda, delta, fcoe)
  
    use mod_basis, only: nglx, ngly, nglz, npts
    
    use mod_input, only: geometry_type, min_length_flg, max_length_flg, ladd_full_stress_tensor

    use mod_types, only: r8
    
    implicit none

    integer, intent(in)  :: ie
    real,    intent(out) :: dx,dy,dz, dx_mean, dy_mean, dz_mean, ds_mean
    real,    intent(out) :: Lambda, delta, fcoe

    real(kind=r8)        :: area_el
    real(kind=r8)        ::  y_flg
    real(kind=r8)        :: delta2
    real(kind=r8)        :: deltai, deltak
    real(kind=r8)        :: a1, a2, fcoe2
    integer              :: nop, nop2
    integer              :: nopx, nopy, nopz
    
    !Auxiliary parameters:
    real(kind=r8), parameter :: oneoverthree = 1.0/3.0
    real(kind=r8), parameter :: fourovertwentyseven = 4.0/27.0
    
    !Begin computations:
    nopx = nglx - 1
    nopy = ngly - 1
    nopz = nglz - 1
    nop  = max(nglx, ngly, nglz) - 1
    nop2 = nop*nop

    call compute_element_size(dx,dy,dz,dx_mean,dy_mean,dz_mean,ie)

    ds_mean = min_length_flg*min(dx_mean, dy_mean*y_flg, dz_mean) + max_length_flg*max(dx_mean, dy_mean*y_flg, dz_mean)

    area_el = dx*dy*dz
    !delta   = (area_el*pi2/nop2*ds_mean)**oneoverthree !delta according to Karmano, Sherwin, and Morrison
    !if(lcylinder) then
    !   delta  = sqrt(dx_mean*dz_mean)
    !   delta = (dx*dz/(nopx*nopz))**oneoverthree
    !else
    !   delta  = (dx_mean*dy_mean*dz_mean)**oneoverthree
    delta = (dx*dy*dz/(nopx*nopy*nopz))**oneoverthree
    !end if
    !delta2 = delta*delta
    delta2 = ds_mean*ds_mean

    !Get the two smaller dimensions of the cell:
    if(dx > dy .and. dx > dz) then
       deltai = dy
       deltak = dz
    else if(dy > dx .and. dy > dz) then
       deltai = dx
       deltak = dz
    else if(dz > dx .and. dz > dy) then
       deltai = dx
       deltak = dy
    end if
    a1 = deltai/max(dx, dy, dz)
    a2 = deltak/max(dx, dy, dz)

    fcoe = cosh(fourovertwentyseven*(log(a1)*log(a1) - log(a1)*log(a2) + log(a2)*log(a2)))
    
    !Compute Lambda:
    Lambda = fcoe*delta
    
end subroutine compute_element_size_driver

  
subroutine compute_element_size(dx,dy,dz, dx_mean, dy_mean, dz_mean,ie)

    use mod_basis, only: nglx, ngly, nglz
    
    use mod_grid, only: coord, npoin, nelem, intma
    !global arrays
    real    :: dx, dy, dz
    real    :: dx_mean, dy_mean, dz_mean
    integer :: ie

    !local arrays
    real    :: x(8), y(8), z(8)
    integer :: inode(8)
    integer :: i, j, k, m

    inode(1)=intma(1,1,1,ie)
    inode(2)=intma(nglx,1,1,ie)
    inode(3)=intma(1,ngly,1,ie)
    inode(4)=intma(nglx,ngly,1,ie)
    inode(5)=intma(1,1,nglz,ie)
    inode(6)=intma(nglx,1,nglz,ie)
    inode(7)=intma(1,ngly,nglz,ie)
    inode(8)=intma(nglx,ngly,nglz,ie)

    do m=1,8
        x(m)=coord(1,inode(m))
        y(m)=coord(2,inode(m))
        z(m)=coord(3,inode(m))
    end do !m

    !Element sizes (as if it were linear)
    dx=maxval(x(:))-minval(x(:))
    dy=maxval(y(:))-minval(y(:))
    dz=maxval(z(:))-minval(z(:))

    !Average distance between LGL points inside the element:
    dx_mean = dx/max(nglx-1,1)
    dy_mean = dy/max(ngly-1,1)
    dz_mean = dz/max(nglz-1,1)

end subroutine compute_element_size

