!------------------------------------------------------------------------------------------------
!>@brief Generates an LGL 3D Mesh
!>@details
!> All arrays below are GLOBAL arrays (should have _g ending)
!> coord_cg = npoin_cg x 3 matrix with (x,y,z) coordinates
!> indexg(2,npoin_cg) = horizontal, vertical indices for physics
!> intma_table = nelem x ngl x ngl x ngl elementy connectivity matrix
!> iperiodic = npoin_cg by 1 vector with periodicity info
!> bsido = nboun by 6 matrix with boundary data
!> bsido(ib,1:4) = verticies of ib-th face st normal points outward
!> bsido(ib,5) = element to which ib-th face belongs
!> bsido(ib,6) = BC enforced on ib-th face
!> npoin_cg = number of grid points
!> xgl = LGL points
!> ngl = number of LGL points 
!> nelem = number of elements
!> nelex = number of elements in x-direction
!> neley = number of elements in y-direction    
!> nelez = number of elements in z-direction   
!> xmin, xmax = x-dimensions of mesh
!> ymin, ymax=  y-dimensions of mesh
!> zmin, zmax = z-dimensions of mesh
!>@author James F. Kelly
!> 			Naval Postgraduate School
!>@date 25 August 2009
!>@date 11 September 2009
!>
!>@date August 2013 Simone Marras 
!> Grid node clustering. This subroutine creates the grid with the option
!> of clustering and stretching the nodes along all directions. If no 
!> stretching option is set in input by the user, a regular grid will 
!> be built by default
!--------------------------------------------------------------------------------------------------
subroutine create_grid_cube(coord_cg,indexg,intma_table,ele_col,bsido,npoin_cg,npoin,nelem,nboun,&
    xglx,xgly,xglz,nglx,ngly,nglz,nelex,neley,nelez,xmin,xmax,ymin,ymax,zmin,zmax,&
    nx,ny,nz,iboun,xperiodic,yperiodic,zperiodic)

    use mod_input, only: xstretch_coe, ystretch_coe, zstretch_coe, &
        lxstretch, lystretch, lzstretch, &
        xstretch_type, ystretch_type, zstretch_type, space_method
  
    implicit none

    !global
    integer :: npoin_cg, npoin, nelem, nboun, ngl, nglx, ngly, nglz, nelex, neley, nelez
    real    :: xmin, xmax, ymin, ymax, zmin, zmax
    real    :: coord_cg(3,npoin_cg)
    integer :: indexg(2,npoin_cg)
    real    :: xglx(nglx), xgly(ngly), xglz(nglz)
    integer :: nx, ny, nz
    integer :: intma_table(nglx,ngly,nglz,nelem)
    integer :: bsido(6,nboun), iboun(6)
    logical :: xperiodic, yperiodic, zperiodic, distort_grid
    integer :: ele_col(nelem)

    !local
    integer :: ix, iy, iz, ip, ib, ip1, ip2, ies
    integer :: ielex, ieley, ielez, icell
    integer :: iglx, igly, iglz
    integer :: ilglz, ilgly, ilglx,ie
    integer :: i1, j1, k1, ii, jj, kk, ih
    integer :: nglm1, i, j, k, e
    real    :: x, y, z, xstart, ystart, zstart, r, ztop
    real    :: Lx, Ly, Lz
    real    :: flag, z0,h
    real    :: xx(nglx*nelex + 1), yy(ngly*neley + 1), zz(nglz*nelez + 1)

    integer, dimension(:,:,:), allocatable :: inode
    integer :: AllocateStatus

    !Local stretching variables:
    real, dimension(:), allocatable :: ksi,     eta,     zeta     !Logical space
    real                            :: ksi_ref, eta_ref, zeta_ref !Logical space

    allocate (inode(nx,ny,nz),   stat=AllocateStatus )
    if (AllocateStatus /= 0) stop "** Not Enough Memory - CREATE_GRID **"
    allocate (ksi(nelex+1),   stat=AllocateStatus )
    if (AllocateStatus /= 0) stop "** Not Enough Memory - CREATE_GRID **"
    allocate (eta(neley+1),   stat=AllocateStatus )
    if (AllocateStatus /= 0) stop "** Not Enough Memory - CREATE_GRID **"
    allocate (zeta(nelez+1),   stat=AllocateStatus )
    if (AllocateStatus /= 0) stop "** Not Enough Memory - CREATE_GRID **"

  
    !-------------------------------------------------------------------
    !Build the logical (uniform) mesh in (0,1)
    !-------------------------------------------------------------------
    !In the ksi direction
    Lx = 1.0/nelex
    ksi(1) = 0.0
    do icell = 1,nelex
        ix = icell + 1
        ksi(ix) = ksi(ix-1) + Lx
    end do
    ksi(nelex+1) = 1.0
    ksi_ref = (ksi(nelex+1) - ksi(1))/2.0

    !In the Eta direction
    Ly = 1.0/neley
    eta(1) = 0.0
    do icell = 1,neley
        iy = icell + 1
        eta(iy) = eta(iy-1) + Ly
    end do
    eta(neley+1) = 1.0
    eta_ref = (eta(neley+1) - eta(1))/2.0

    !In the Zeta direction
    Lz = 1.0/nelez
    zeta(1) = 0.0
    do icell = 1,nelez
        iz = icell + 1
        zeta(iz) = zeta(iz-1) + Lz
    end do
    zeta(nelez+1) = 1.0
    zeta_ref = (zeta(nelez+1) - zeta(1))/2.0

    !Force no clustering when clustering_type_* == 'none'
    !X:
    if( xstretch_type == 'none') xstretch_coe= 0.0
    !Y:
    if( ystretch_type == 'none') ystretch_coe = 0.0
    !Z:
    if( zstretch_type == 'none') zstretch_coe= 0.0

    !>  nglm1 = ngl - 1
    Lx = (xmax - xmin)/nelex
    Ly = (ymax - ymin)/neley
    Lz = (zmax - zmin)/nelez

    ! Generate Coordinates Here
    ip = 1
    ii = 0
    jj = 0
    kk = 0

    zz(1) = zmin
    do ielez = 1, nelez
        iz = ielez + 1

        if (ielez == 1) then
            k1 = 1
        else
            k1 = 2
        end if

        if(.not.lzstretch) then
            !
            ! No clustering:
            !
            zstart = zmin + (ielez - 1)*Lz
            zz(iz-1) = 0.0
        else
            !
            ! Clustering (defined by the function 'h'):
            !
            zstart = 0.0

            if( zstretch_type == 'one_sided' ) then
                !
                ! one_sided
                !
                h = (zmax-zmin)*(exp(zstretch_coe*zeta(iz)) - 1.0)/(exp(zstretch_coe) - 1.0)

            else if( zstretch_type == 'middle' ) then
                !
                ! middle
                !
                if( zeta(iz) >= 0.0 .and. zeta(iz) <= zeta_ref) then
                    h = (zmax-zmin)*zeta_ref*( (exp(zstretch_coe) - exp(zstretch_coe*(1.0 - &
                        zeta(iz)/zeta_ref)))/(exp(zstretch_coe) - 1.0))
                else if( zeta(iz) > zeta_ref .and. zeta(iz) <= 1.0) then
                    h = (zmax-zmin)*zeta_ref + (zmax-zmin)*(1.0 - zeta_ref)*(exp(zstretch_coe* &
                        (zeta(iz)-zeta_ref)/(1.0 - zeta_ref)) - 1.0)/(exp(zstretch_coe) - 1.0)
                end if

            else
                !
                ! symmetric (Eriksson)
                !
                if(zeta(iz) <= zeta_ref ) then
                    h = (zmax-zmin)*zeta_ref*(exp(zstretch_coe*zeta(iz)/zeta_ref) - &
                        1.0)/(exp(zstretch_coe) - 1.0)
                else

                    h = (zmax-zmin) - (zmax-zmin)*(1.0-zeta_ref)*(exp(zstretch_coe* &
                        (1-zeta(iz))/(1-zeta_ref)) - 1.0)/(exp(zstretch_coe) - 1.0)
                end if
            end if

            zz(iz) = zmin + h

            Lz = zz(iz) - zz(iz-1)
        end if

        do ilglz = k1, nglz
            kk = kk + 1

            z = zstart + zz(iz-1) + 0.5*Lz*(1 + xglz(ilglz))

            !y
            jj = 0
            !>        yy(1) = ymin !!!!
            do ieley = 1, neley
                iy = ieley + 1

                if (ieley == 1) then
                    j1 = 1
                else
                    j1 = 2
                end if

                if(.not.lystretch) then
                    !
                    ! No clustering in the x direction
                    !
                    ystart = ymin + (ieley - 1)*Ly
                    yy(iy-1) = 0.0

                else
                    !
                    ! Clustering (defined by the function 'h'):
                    !
                    ystart = 0.0
                    if( ystretch_type == 'one_sided' ) then
                        !
                        ! one_sided
                        !
                        h = (ymax-ymin)*(exp(ystretch_coe*eta(iy)) - 1.0)/(exp(ystretch_coe) - 1.0)
                    else if( ystretch_type == 'middle' ) then
                        !
                        ! middle
                        !
                        if( eta(iy) >= 0.0 .and. eta(iy) <= eta_ref) then
                            h = (ymax-ymin)*eta_ref*( (exp(ystretch_coe) - exp(ystretch_coe* &
                            (1.0 - eta(iy)/eta_ref)))/(exp(ystretch_coe) - 1.0))

                        else if( eta(iy) > eta_ref .and. eta(iy) <= 1.0) then
                            h = (ymax-ymin)*eta_ref + (ymax-ymin)*(1.0 - eta_ref)* &
                                (exp(ystretch_coe*(eta(iy)-eta_ref)/(1.0 - eta_ref)) - &
                                1.0)/(exp(ystretch_coe) - 1.0)
                        end if
                    else 
                        !
                        ! symmetric
                        !
                        if(eta(iy) <= eta_ref ) then
                            h = (ymax-ymin)*eta_ref*(exp(ystretch_coe*eta(iy)/eta_ref) - &
                                1.0)/(exp(ystretch_coe) - 1.0)
                        else
                            h = (ymax-ymin) - (ymax-ymin)*(1.0-eta_ref)*(exp(ystretch_coe* &
                                (1-eta(iy))/(1-eta_ref)) - 1.0)/(exp(ystretch_coe) - 1.0)
                        end if
                    end if

                    yy(iy) = ymin + h
                    Ly = yy(iy) - yy(iy-1)
                end if

                do ilgly = j1, ngly
                    jj = jj + 1

                    y= ystart + yy(iy-1) + 0.5*Ly*(1 + xgly(ilgly))

                    !x
                    ii = 0
                    xx(1) = xmin
                    do ielex = 1,nelex
                        ix = ielex + 1

                        if (ielex == 1) then
                            i1 = 1
                        else
                            i1 = 2
                        end if

                        if(.not.lxstretch) then
                            !
                            ! No clustering in the x direction
                            !
                            xstart = xmin + (ielex - 1)*Lx
                            xx(ix-1) = 0.0

                        else
                            !
                            ! Clustering (defined by the function 'h'):
                            !
                            xstart = 0.0
                            if( xstretch_type == 'one_sided' ) then
                                !
                                ! one_sided
                                !
                                h = (xmax-xmin)*(exp(xstretch_coe*ksi(ix)) - 1.0) &
                                    /(exp(xstretch_coe) - 1.0)
                            else if( xstretch_type == 'middle' ) then
                                !
                                ! middle
                                !
                                if( ksi(ix) >= 0.0 .and. ksi(ix) <= ksi_ref) then
                                    h = (xmax-xmin)*ksi_ref*( (exp(xstretch_coe) - &
                                        exp(xstretch_coe*(1.0 - ksi(ix)/ksi_ref)))&
                                        /(exp(xstretch_coe) - 1.0))
                                else if( ksi(ix) > ksi_ref .and. ksi(ix) <= 1.0) then
                                    h = (xmax-xmin)*ksi_ref + (xmax-xmin)*(1.0 - ksi_ref)* &
                                        (exp(xstretch_coe*(ksi(ix)-ksi_ref)/(1.0 - ksi_ref)) - &
                                        1.0)/(exp(xstretch_coe) - 1.0)
                                end if
                            else
                                !
                                ! symmetric
                                !
                                if(ksi(ix) <= ksi_ref ) then
                                    h = (xmax-xmin)*ksi_ref*(exp(xstretch_coe*ksi(ix)/ksi_ref) - &
                                        1.0)/(exp(xstretch_coe) - 1.0)
                                else
                                    h = (xmax-xmin) - (xmax-xmin)*(1.0-ksi_ref)* &
                                        (exp(xstretch_coe*(1-ksi(ix))/(1-ksi_ref)) - 1.0) &
                                        /(exp(xstretch_coe) - 1.0)
                                end if
                            end if

                            xx(ix) = xmin + h
                            Lx = xx(ix) - xx(ix-1)
                        end if

                        do ilglx = i1, nglx
                            ii = ii + 1

                            x = xstart + xx(ix-1) + 0.5*Lx*(1 + xglx(ilglx))

                            ie = ielex + nelex*(ieley -1) + nelex*neley*(ielez -1)
                            coord_cg(1,ip) = x
                            coord_cg(2,ip) = y
                            coord_cg(3,ip) = z
                            indexg(1,ip) = (jj-1)*nx+ii
                            indexg(2,ip) = kk
                            inode(ii,jj,kk) = ip
                            ip = ip + 1
                        end do
                    end do
                end do
            end do
        end do
    end do

    !Construct Intma
    ie = 0
    do ielez = 1, nelez
        ies = 0
        do ieley = 1,neley
            do ielex = 1,nelex
                ies = ies + 1
                ie = ie + 1

                ele_col(ie) = ies
                do ilglz = 1, nglz
                    kk = (ielez - 1)*(nglz-1) + ilglz
                    do ilgly = 1,ngly
                        jj = (ieley - 1)*(ngly-1) + ilgly
                        do ilglx = 1,nglx
                            ii = (ielex - 1)*(nglx-1) + ilglx
                            !>                    ih=nx*(ieley-1)*(ngly-1) + nx*(ilgly-1) + ii
                            ip=inode(ii,jj,kk)
                            intma_table(ilglx,ilgly,ilglz,ie)=ip
                        end do
                    end do
                end do
            end do
        end do
    end do

    ! Create Boundary Condition Data
    call create_bsido(bsido,intma_table,inode,npoin_cg,nboun,nelem, &
        nelex,neley,nelez,iboun,nglx,ngly,nglz, nx, ny, nz)


    !Free local memory
    deallocate(ksi,eta,zeta)

end subroutine create_grid_cube

!----------------------------------------------------------------------!
!>@brief Creates Boundary condition data in data structure bsido(6,nboun)
!>@details Requires inode(nx,ny,nz) and iboun(6) as input
!>@todo Fix face oreintation
!----------------------------------------------------------------------!
subroutine create_bsido(bsido,intma,inode,npoin_cg,nboun,nelem,nelex,neley,nelez,iboun, &
    nglx, ngly, nglz, nx, ny, nz)

    implicit none

    integer :: nboun, nelem, nelex,neley,nelez, nglx, ngly, nglz, nglm1
    integer :: nx, ny, nz
    integer :: npoin_cg
    integer :: bsido(6,nboun)
    integer :: iboun(6)
    integer :: inode(nx,ny,nz)
    integer :: ip1, ip2, ip3, ip4, ip5, ip6, ip7, ip8
    integer :: ii, jj, kk
    integer :: ielex, ieley, ielez
    integer :: ib, ie, ip
    integer :: ibx, iby, ibz
    integer :: ix, iy, iz
    integer, dimension(nglx,ngly,nglz,nelem)::intma
    !  nglm1 = ngl - 1
  
    !Construct Boundary Conditions using the alphabetical convention
    ib=0
    bsido=0
  
    ! Bottom Boundary xy, z = -1
    if(nglz > 1) then
  
        ielez = 1
        do ielex =1,nelex
            do ieley = 1,neley
                ib=ib+1
                ! Node 1
                ii=(ielex-1)*(nglx-1) + 1
                jj=(ieley-1)*(ngly-1) + 1
                kk=1
                ip1=inode(ii,jj,kk)
                ! Node 3
                ii=(ielex-1)*(nglx-1) + 1
                jj=(ieley-1)*(ngly-1) + ngly
                kk=1
                ip3=inode(ii,jj,kk)
                ! Node 4
                ii=(ielex-1)*(nglx-1) + nglx
                jj=(ieley-1)*(ngly-1) + ngly
                kk=1
                ip4=inode(ii,jj,kk)
                ! Node 2
                ii=(ielex-1)*(nglx-1) + nglx
                jj=(ieley-1)*(ngly-1) + 1
                kk=1
                ip2=inode(ii,jj,kk)

                ! Get element number
                ie = ielex + nelex*(ieley -1) + nelex*neley*(ielez -1)

                ! Populate BSIDO
                bsido(1,ib) = ip1
                bsido(2,ib) = ip3
                bsido(3,ib) = ip4
                bsido(4,ib) = ip2
                bsido(5,ib) = ie
                bsido(6,ib) = iboun(1)
            end do
        end do

        ! Top Boundary xy, z = +1
        !>  if (iboun(2) /= 3) then !don't store for Periodic
        ielez = nelez
        do ielex = 1,nelex
            do ieley = 1,neley
                ib = ib + 1
                ! Node 5
                ii = (ielex - 1)*(nglx-1) + 1
                jj = (ieley - 1)*(ngly-1) + 1
                kk = nz
                ip5 = inode(ii,jj,kk)
                ! Node 6
                ii = (ielex - 1)*(nglx-1) + nglx
                jj = (ieley - 1)*(ngly-1) + 1
                kk = nz
                ip6 = inode(ii,jj,kk)
                ! Node 8
                ii = (ielex - 1)*(nglx-1) + nglx
                jj = (ieley - 1)*(ngly-1) + ngly
                kk = nz
                ip8 = inode(ii,jj,kk)
                ! Node 7
                ii = (ielex - 1)*(nglx-1) + 1
                jj = (ieley - 1)*(ngly-1) + ngly
                kk = nz
                ip7 = inode(ii,jj,kk)

                ! Get element number
                ie = ielex + nelex*(ieley -1) + nelex*neley*(ielez -1)

                ! Populate BSIDO
                bsido(1,ib) = ip5
                bsido(2,ib) = ip6
                bsido(3,ib) = ip8
                bsido(4,ib) = ip7
                bsido(5,ib) = ie
                bsido(6,ib) = iboun(2)
            end do
        end do
    !>  end if !iboun(2) == 3
    endif
   
    if(ngly > 1) then
        ! Left Boundary xz, y = -1
        ieley = 1
        do ielex =1,nelex
            do ielez = 1,nelez
                ib=ib+1
                ! Node 1
                ii=(ielex-1)*(nglx-1) + 1
                jj = 1
                kk = (ielez -1)*(nglz-1) + 1
                ip1 = inode(ii,jj,kk)
                ! Node 2
                ii=(ielex-1)*(nglx-1) + nglx
                jj = 1
                kk = (ielez -1)*(nglz-1) + 1
                ip2 = inode(ii,jj,kk)
                ! Node 6
                ii=(ielex-1)*(nglx-1) + nglx
                jj = 1
                kk = (ielez -1)*(nglz-1) + nglz
                ip6 = inode(ii,jj,kk)
                ! Node 5
                ii=(ielex-1)*(nglx-1) + 1
                jj = 1
                kk = (ielez -1)*(nglz-1) + nglz
                ip5 = inode(ii,jj,kk)

                ! Get element number
                ie = ielex + nelex*(ieley -1) + nelex*neley*(ielez -1)

                ! Populate BSIDO
                bsido(1,ib) = ip1
                bsido(2,ib) = ip2
                bsido(3,ib) = ip6
                bsido(4,ib) = ip5
                bsido(5,ib) = ie
                bsido(6,ib) = iboun(3)
            end do
        end do

        ! Right Boundary xz, y = 1
        !>  if (iboun(4) /= 3) then !don't store if Periodic
        ieley = neley
        do ielex = 1,nelex
            do ielez = 1,nelez
                ib = ib + 1
                ! Node 3
                ii = (ielex - 1)*(nglx-1) + 1
                jj = ny
                kk = (ielez - 1)*(nglz-1) + 1
                ip3 = inode(ii,jj,kk)
                ! Node 7
                ii = (ielex - 1)*(nglx-1) + 1
                jj = ny
                kk = (ielez - 1)*(nglz-1) + nglz
                ip7 = inode(ii,jj,kk)
                ! Node 8
                ii = (ielex - 1)*(nglx-1) + nglx
                jj = ny
                kk = (ielez - 1)*(nglz-1) + nglz
                ip8 = inode(ii,jj,kk)
                ! Node 4
                ii = (ielex - 1)*(nglx-1) + nglx
                jj = ny
                kk = (ielez - 1)*(nglz-1) + 1
                ip4 = inode(ii,jj,kk)

                ! Get element number
                ie = ielex + nelex*(ieley -1) + nelex*neley*(ielez -1)

                ! Populate BSIDO
                bsido(1,ib) = ip3
                bsido(2,ib) = ip7
                bsido(3,ib) = ip8
                bsido(4,ib) = ip4
                bsido(5,ib) = ie
                bsido(6,ib) = iboun(4)
            end do
        end do
    !>  end if
    endif
  
    if(nglx > 1) then
        ! Back Boundary yz, x = -1
        ielex = 1
        do ieley = 1,neley
            do ielez = 1, nelez
                ib = ib + 1
                ! Node 1
                ii = 1
                jj = (ieley - 1)*(ngly-1) + 1
                kk = (ielez - 1)*(nglz-1) + 1
                ip1 = inode(ii,jj,kk)
                ! Node 5
                ii = 1
                jj = (ieley - 1)*(ngly-1) + 1
                kk = (ielez - 1)*(nglz-1) + nglz
                ip5 = inode(ii,jj,kk)
                ! Node 7
                ii = 1
                jj = (ieley - 1)*(ngly-1) + ngly
                kk = (ielez - 1)*(nglz-1) + nglz
                ip7 = inode(ii,jj,kk)
                ! Node 3
                ii = 1
                jj = (ieley - 1)*(ngly-1) + ngly
                kk = (ielez - 1)*(nglz-1) + 1
                ip3 = inode(ii,jj,kk)

                ! Get element number
                ie = ielex + nelex*(ieley -1) + nelex*neley*(ielez -1)

                ! Populate BSIDO
                bsido(1,ib) = ip1
                bsido(2,ib) = ip5
                bsido(3,ib) = ip7
                bsido(4,ib) = ip3
                bsido(5,ib) = ie
                bsido(6,ib) = iboun(5)
            end do
        end do

        ! Front Boundary yz, x = 1
        ! if (iboun(6) .ne. 6) then
        !>  if (iboun(6) /= 3) then !don't store if Periodic
        ielex = nelex
        do ieley = 1,neley
            do ielez = 1, nelez
                ib = ib + 1
                ! Node 2
                ii = nx
                jj = (ieley - 1)*(ngly-1) + 1
                kk = (ielez - 1)*(nglz-1) + 1
                ip2 = inode(ii,jj,kk)
                ! Node 4
                ii = nx
                jj = (ieley - 1)*(ngly-1) + ngly
                kk = (ielez - 1)*(nglz-1) + 1
                ip4 = inode(ii,jj,kk)
                ! Node 8
                ii = nx
                jj = (ieley - 1)*(ngly-1) + ngly
                kk = (ielez - 1)*(nglz-1) + nglz
                ip8 = inode(ii,jj,kk)
                ! Node 6
                ii = nx
                jj = (ieley - 1)*(ngly-1) + 1
                kk = (ielez - 1)*(nglz-1) + nglz
                ip6 = inode(ii,jj,kk)

                ! Get element number
                ie = ielex + nelex*(ieley -1) + nelex*neley*(ielez -1)

                ! Populate BSIDO
                bsido(1,ib) = ip2
                bsido(2,ib) = ip4
                bsido(3,ib) = ip8
                bsido(4,ib) = ip6
                bsido(5,ib) = ie
                bsido(6,ib) = iboun(6)
            end do
        end do
    !>  end if
    endif
  
end subroutine create_bsido
