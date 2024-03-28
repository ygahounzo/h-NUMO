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

subroutine create_mass_element(mass_e, massinv_e, jacq)


    use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, npts, nqx, nqy

    use mod_grid, only: intma, npoin, nelem

    implicit none
    real :: jacq(nqx,nqy,nqz,nelem)
    real :: mass_e(npts, npts)
    real :: massinv_e(npts, npts)

    integer :: iquad, jquad, e, m, n, l, s, jp, ip
    real :: wq, h_j, h_i
    real, dimension(npts) :: work  ! work array for LAPACK
    integer, dimension(npts) :: ipiv   ! pivot indices
    integer :: info

    ! External procedures defined in LAPACK
    ! external DGETRF
    ! external DGETRI

    ! Initialize
    mass_e = 0.0d0

    ! Compute mass matrix for element e
    e = 1

    ! Do LGL Integration
    do jquad = 1, nqy
        do iquad = 1, nqx
            wq = jacq(iquad, jquad, 1, e)

            do m = 1, ngly
                do n = 1, nglx
                    jp = intma(n, m, 1, e)
                    h_j = psiqx(n, iquad) * psiqy(m, jquad)

                    do l = 1, ngly
                        do s = 1, nglx
                            ip = intma(s, l, 1, e)
                            h_i = psiqx(s, iquad) * psiqy(l, jquad)
                            mass_e(ip, jp) = mass_e(ip, jp) + wq * h_i * h_j
                        end do
                    end do
                end do
            end do
        end do
    end do

    ! Here we invert the mass matrix using LAPACK
    ! Store mass_e in massinv_e to prevent it from being overwritten by LAPACK
    massinv_e = mass_e

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(npts, npts, massinv_e, npts, ipiv, info)

    if (info /= 0) then
        stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(npts, massinv_e, npts, ipiv, work, npts, info)

    if (info /= 0) then
        stop 'Matrix inversion failed!'
    end if

    ! do ip = 1, npts
    !     print *, massinv_e(ip, :)
    ! end do

    ! stop

end subroutine create_mass_element


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



