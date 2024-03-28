!-----------------------------------------------------------------!
!>@brief Creates a 3D spherical shell grid using spectral elements
!>@details Output:
!> coord_cg(3,npoin_cg) = (x,y,z) coords
!> indexg(2,npoin_cg) = horizontal, vertical indices for physics
!> intma(ngl,ngl,ngl,nelem) = interconnectivity matrix (standard 3D Format)
!> intma_lev(ngl,ngl,nelem_s,npoin_r) = horizontal level-based interconnectivity matrix
!> intma_vertical(npoin_s,ngl,nelem_r) = vertical level-based interconnectivity matrix
!> ele_col(nelem) = level-element that owns global element
!> bsido(6,nboun) = boundary matrix
!> Input: 
!> nelem_r = number of elements in the radial direction
!> nelem_s = number of elements in each concentric sphere
!> ngl = number of LGL points
!> nelem = total number of elements = nelem_r*nelem_s
!> npoin_cg = total number of points
!> nboun = number of boundary elements = 2*nelem_s
!> npoin_r = number of grid points in the radial direction = nop*nelem_r + 1
!> npoin_s = number of grid points on each concentric sphere
!> rmin = radius of the inner (ground) sphere
!> rmax = radius of the outer (top) sphere
!>@author James F. Kelly 26 April 2010 
!>@date 9 May 2010
!>@date 29 June 2010
!>@date 4 August 2010 for Icosohedral Grids
!-----------------------------------------------------------------!
subroutine create_grid_sphere_ico(coord_cg,indexg,intma,intma_lev,intma_s,ele_col,bsido,&
    iboun,nelem_r,nelem_s,nel,nglx,ngly,nglz,nelem,npoin_cg,nboun,npoin_r,npoin_s,rmin,rmax)
  
    !Global
    implicit none
    real coord_cg(3,npoin_cg)
    integer indexg(2,npoin_cg)
    integer intma(nglx,ngly,nglz,nelem)
    integer intma_lev(nglx,ngly,nelem_s,npoin_r)
    integer ele_col(nelem)
    integer bsido(6,nboun), iboun(2)
    integer nboun, npoin_r, npoin_s,nel
    integer nglx,ngly,nglz,nelem,nelem_r,nelem_s,nop,npoin_cg
    real rmin, rmax
  
    integer ier, ipr, ies, ips, ie, ip, ib
    integer i, j, k
    real dr
  
    real r(npoin_r)
    integer intma_r(nglz,nelem_r)
    integer nface, nelem0, nx, ny, npoin0
    real r2
  
    !Global arrays
    real coord_s(3,npoin_s), coordh(npoin_s,3)
    integer intma_s(nglx,ngly,nelem_s), intmah(nelem_s,nglx,ngly)
    integer npoint, nelemt, nsidet
  
  
    nface=6
    nop = nglx - 1 !assumes that ngly=nglx
    npoint=10*nel*nel + 2
    nelemt=2*(npoint-2)
  
  
    !Create 1D Spectral Element Grid in R
    call create1d_grid(r,intma_r,nglz,nelem_r,npoin_r,rmin,rmax)
  
    !Create 2D Spherical Grid
    call ico_quad(coordh,intmah,npoint,nelemt,npoin_s,nelem_s, &
        nop,nglx,nel)
    print *, "Icosohedral Grid Constructed"
    !Construct Cartesian Product [rmin,rmax] x S^2 to create
    !a spherical shell
    do ier = 1,nelem_r
        do k = 1,nglz
        
            ipr = intma_r(k,ier)
            do ies = 1,nelem_s
                do i = 1,nglx
                    do j = 1,ngly
                        intma_s(i,j,ies) = intmah(ies,i,j)
                        ips = intma_s(i,j,ies)
                        ie = (ier - 1)*nelem_s + ies
                        ip = (ipr - 1)*npoin_s + ips
                        intma(i,j,k,ie) = ip                !3D Intma
                        intma_lev(i,j,ies,ipr) = ip         !Level-Based Intma
                        ele_col(ie) = ies                   !Level-Based element that belongs to ie
                        coord_s(:,ips) = coordh(ips,:)
                        coord_cg(1,ip) = r(ipr)*coord_s(1,ips)
                        coord_cg(2,ip) = r(ipr)*coord_s(2,ips)
                        coord_cg(3,ip) = r(ipr)*coord_s(3,ips)
                        indexg(1,ip) = ips
                        indexg(2,ip) = ipr
                    end do
                end do
            end do
        end do
    end do
    print *, "Coords Created"
  
    !Create Bsido
    !NFBC on the Ground
    if(nglz > 1) then
        ier = 1
        ib = 1
        do ies =1,nelem_s
            ie = (ier - 1)*nelem_s + ies
            k = 1
     
            i = 1; j = 1
            ip = intma(i,j,k,ie)
            bsido(1,ib) = ip
     
            i = 1; j = ngly
            ip = intma(i,j,k,ie)
            bsido(2,ib) = ip
     
            i = nglx; j = ngly
            ip = intma(i,j,k,ie)
            bsido(3,ib) = ip
     
            i = nglx; j = 1
            ip = intma(i,j,k,ie)
            bsido(4,ib) = ip
            bsido(5,ib) = ie
            bsido(6,ib) = iboun(1)
     
            ib = ib + 1
        end do
  
        !NRBC at the top of the atmosphere
        ier = nelem_r
        do ies =1,nelem_s
            ie = (ier - 1)*nelem_s + ies
            k = nglz
     
            i = 1; j = 1
            ip = intma(i,j,k,ie)
            bsido(1,ib) = ip
     
            i = nglx; j = 1
            ip = intma(i,j,k,ie)
            bsido(2,ib) = ip
     
            i = nglx; j = ngly
            ip = intma(i,j,k,ie)
            bsido(3,ib) = ip
     
            i = 1; j = ngly
            ip = intma(i,j,k,ie)
            bsido(4,ib) = ip
            bsido(5,ib) = ie
            bsido(6,ib) = iboun(2)
     
            ib = ib + 1
     
        end do
    endif
end subroutine create_grid_sphere_ico

