!----------------------------------------------------------------------!
!>@brief This subroutine constructs the Mass matrix for the SEM using 
!>Quadrilateral Elements for the Euler equations
!----------------------------------------------------------------------!
subroutine create_mass(mass,jac)
  
    use mod_basis, only: nglx, ngly, nglz

    use mod_grid, only: intma, npoin, nelem

    implicit none
  
    !global arrays
    real mass(npoin), jac(nglx,ngly,nglz,nelem)
  
    !local
    integer ie, i, j, k, ip
    real wq
  
    !initialize the global matrix
    mass=0
    
    !loop thru the elements
    do ie=1,nelem
       !Do Numerical Integration
       do k = 1,nglz
          do j=1,ngly
             do i=1,nglx
                wq=jac(i,j,k,ie)
                ip = intma(i,j,k,ie) !FIX-ME iperiodic
                mass(ip)=mass(ip) + wq
             end do
          end do
       end do
    end do
  
  !Periodicity

end subroutine create_mass


!----------------------------------------------------------------------!
!>@brief This subroutine constructs the Mass matrix along the Vertical Only (Columns)
!>for the SEM using Quadrilateral Elements for the Euler equations
!>@author  F.X. Giraldo 12/14/2010
!----------------------------------------------------------------------!
subroutine create_mass_column(mass_1d,jac_1d)
  
    use mod_basis, only: nglz, wglz
  
    use mod_grid, only: coord, intma_1d, ncol, nz, node_column

    use mod_input, only: nelz
  
    implicit none
  
    !global arrays
    real mass_1d(ncol,nz), jac_1d(ncol,nglz,nelz)
  
    !local
    integer ie, ix, iez, i, j, k, ip, i0, i1
    real wq, ds, dz
  
    !initialize
    mass_1d=0
    jac_1d=0

    do ix=1,ncol !loop through horizontal (rows that define a vertical column)
        do iez=1,nelz !loop through vertical components of the column
     
            !Find element size along z/r-direction
            i0=intma_1d(1,iez)
            i0=node_column(ix,i0)
            i1=intma_1d(nglz,iez)
            i1=node_column(ix,i1)
            ds=sqrt( (coord(1,i1) - coord(1,i0))**2 &
                + (coord(2,i1) - coord(2,i0))**2  &
                + (coord(3,i1) - coord(3,i0))**2 ) !Euclidean Distance
            dz=0.5*ds

            !Do Numerical Integration
            do i=1,nglz
                ip=intma_1d(i,iez)
                wq=dz*wglz(i)
                jac_1d(ix,i,iez)=wq
                mass_1d(ix,ip)=mass_1d(ix,ip) + wq
            end do !i
        end do !iez
    end do !ix
  
end subroutine create_mass_column
