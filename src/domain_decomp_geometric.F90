!-----------------------------------------------------------------!
!>@brief This subroutine creates global_proc, which maps global elements to processors, using naive
!> geometric domain decomposition.  Using this approach, all processors receive the same
!> number of processors.  DOMAIN_DECOMP_CUBE requirs geometry_type == 'cube'
!>@author James F. Kelly 31 July 2010
!-----------------------------------------------------------------!
subroutine domain_decomp_cube(global_proc,nprocx,nprocy,nprocz,nproc)
  
    use mod_global_grid, only:  nelem_g

    use mod_input, only: nelx, nely, nelz

    implicit none
  
    ! Number of processors in the x, y, and z directions and total number of processors, nproc
    integer nprocx, nprocy, nprocz, nproc
  
    ! Data structures that map
    ! global_proc: Maps global elements to processors (that is, what processor owns a particular element)
    ! Using these data structures and the global grid (coord_cg and element connectivity),
    ! local grids may be constructed
    integer global_proc(nelem_g)
  
    integer nelx_local, nely_local, nelz_local
    integer ip_x, ip_y, ip_z, iel_x, iel_y, iel_z
    integer ielex, ieley, ielez, ie_g, ip
  
    ! Get number of elements per processor in each direction
    if ( (mod(nelx,nprocx) > 0) .or. (mod(nely,nprocy) > 0) .or. (mod(nelz,nprocz) > 0) ) then
        print*,' Number of elements is not an integer multiple of the number of processors!>'
        stop
    end if
  
    ! Determine the number of elements in each cardinal direction for each processor space
    ! nelem_local = nelx_local*nely_local*nelz_local
    nelx_local = nelx/nprocx
    nely_local = nely/nprocy
    nelz_local = nelz/nprocz
  
    ! Loop over number of processors in each direction and number of elements per proc.
    do ip_z = 1, nprocz
        do iel_z = 1, nelz_local
            ielez =  (ip_z-1)*nelz_local +  iel_z             !Global z-element
        
            do ip_y = 1, nprocy
                do iel_y = 1, nely_local
                    ieley =  (ip_y-1)*nely_local +  iel_y          !Global y-element
                    do ip_x = 1, nprocx
                        do iel_x = 1, nelx_local
                            ielex =  (ip_x-1)*nelx_local +  iel_x       !Global x-element
                    
                            ie_g = ielex + nelx*(ieley -1) + nelx*nely*(ielez -1)     !Global Element
                            ip   = ip_x + nprocx*(ip_y-1) + nprocx*nprocy*(ip_z - 1)        ! Processor Nummber
                    
                            global_proc(ie_g) = ip
                        end do
                    end do
                end do
            end do
        end do
    end do
  
end subroutine domain_decomp_cube

!-----------------------------------------------------------------!
!>@brief This subroutine constructs NSEAM style domain decomposition, where
!> the total number of processors is nproc = 6 * nprocx * nprocy
!> This domain decomposition requires geometry_type == 'sphere"
!>@author James F. Kelly 31 July 2010
!-----------------------------------------------------------------!
subroutine domain_decomp_nseam(global_proc,nprocx,nprocy,nproc)
  
    use mod_global_grid, only:  nelem_g, nelem_s, nelem_r

    use mod_input, only: nelx

    implicit none
  
    integer nproc, nel2, nel
    integer global_proc(nelem_g)
    integer ier, iface, ie_x, ie_y, ies, ie_g, ip_y, iel_y, ip_x, iel_x
    integer nprocx,nprocy,nelx_local,nely_local, ip
  
    nel = nelx
  
    ! Get number of elements per processor in each direction
    if ( (mod(nel,nprocx) > 0) .or. (mod(nel,nprocy) > 0) ) then
        print*,' Number of elements is not an integer multiple of the number of processors!>'
        stop
    end if
  
    nelx_local = nel/nprocx
    nely_local = nel/nprocy
    nel2 = nel*nel
  
    do ier = 1,nelem_r
        do iface = 1,6
            do ip_y = 1, nprocy
                do iel_y = 1, nely_local
                    ie_y =  (ip_y-1)*nely_local +  iel_y          !Global y-element
                    do ip_x = 1, nprocx
                        do iel_x = 1, nelx_local
                            ie_x =  (ip_x-1)*nelx_local +  iel_x       !Global x-element
                    
                            ! Horz. Element Number
                            ies = nel2*(iface - 1) + (ie_x - 1)*nel + ie_y
                            ! Global Element Number
                            ie_g = (ier - 1)*nelem_s + ies
                            ! Processor Nummber
                            ip   = ip_x + nprocx*(ip_y-1) + nprocx*nprocy*(iface - 1)
                            ! Construct Global Element to Processor Mapping
                            !           print *, ie_g, ip
                            global_proc(ie_g) = ip
                        end do
                    end do
                end do
            end do
        end do
    end do
  
end subroutine domain_decomp_nseam

!-----------------------------------------------------------------!
!>@brief Determines the appropriate number of processers in N/S and E/W
!> directions on the cubed-sphere using so-called NSEAM Decomposition
!>@author James F. Kelly 31 July 2010
!-----------------------------------------------------------------!
subroutine nseam_decomp(nprocx,nprocy,nproc)
  
    integer nprocx, nprocy, nproc
  
    if (nproc == 6) then
        nprocx = 1
        nprocy = 1
    else if (nproc == 12) then
        nprocx = 2
        nprocy = 1
    else if (nproc == 18) then
        nprocx = 3
        nprocy = 1
    else if (nproc == 24) then
        nprocx = 2
        nprocy = 2
    else if (nproc == 30) then
        nprocx = 5
        nprocy = 1
    else if (nproc == 36) then
        nprocx = 3
        nprocy = 2
    else if (nproc == 48) then
        nprocx = 4
        nprocy = 2
    else if (nproc == 54) then
        nprocx = 3
        nprocy = 3
    else if (nproc == 60) then
        nprocx = 5
        nprocy = 2
    else if (nproc == 72) then
        nprocx = 4
        nprocy = 3
    else if (nproc == 96) then
        nprocx = 4
        nprocy = 4
    else
        print *, "Nproc must be 6, 12, 18 , 24, 30, 36, 48, 54, 60, 72, or 96"
    end if
  
end subroutine nseam_decomp
