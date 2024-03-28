!------------------------------------------------------------------------!
!>@brief This subroutine controls the entire Grid Generation for a high order 
!>spectral element/discontinuous Galerkin ICOSAHEDRAL grid.
!>@details Each of the 20 faces of the initial icosahedron is subdivided by NICO*NICO 
!>triangular elements. Then each triangular element is subdivided into 
!>3 quadrilaterals. Then in each of these elements (NOP+1)*(NOP+1)
!>collocation points at the Legendre-Gauss-Lobatto points are constructed.
!>@author  F.X. Giraldo 11/00
!>@date FXG on May 7, 2014 for F90
!------------------------------------------------------------------------!
subroutine ico_quad(coordh,intmah,npoin_cg,nelem,npoinh,nelemh,nop,ngl,nico)

    use mod_constants, only:  tol, pi

    implicit none

    integer, parameter :: npoin0=12, nelem0=2*(npoin0-2), nside0=3*(npoin0-2), nd = 4

    !Global arrays
    integer, intent(in) :: npoin_cg, nelem, npoinh, nelemh, nop, ngl, nico
    real, intent(out) :: coordh(npoinh,3)
    integer, intent(out) :: intmah(nelemh,ngl,ngl)

    !Local arrays == Triangle and Quadrilateral
    real :: coorsh(npoinh,2), coord0(npoin0,3), coord_cg(npoin_cg,3), coors(npoin_cg,2)
    integer :: intma0(nelem0,3), intma(nelem,3), intmaq(nelemh,4)
    real :: twopi, pio2, x, y, z, r, alon, alat
    integer :: i, npl, nelemq, npoinq, nsideh, nside

    !Constants
    twopi=2*pi
    pio2=0.5*pi
    npl=nico + 1
    nelemq=nelemh
    nsideh = 120*nico*nico
    nside = 3*(npoin_cg - 2)

    !!--------------------------------------------------------------!!
    !!-----------------Construct the Triangular Grid----------------!!
    !!--------------------------------------------------------------!!

    !Construct Initial Icosahedron
    call init_ico(coord0,intma0,npoin0,nelem0)

    !Store all points in spherical coordinates
    do i=1,npoin0
        x=coord0(i,1)
        y=coord0(i,2)
        z=coord0(i,3)
        alon=atan2(y,x+tol)
        alat=asin(z)
        coors(i,1)=alon
        coors(i,2)=alat
    end do

    !Generate the Higher Order Points NICO in Spherical Coordinates
    call make_tris(coors,intma0,intma,npoin0,nelem0,nside0,npoin_cg, &
        nelem,nico,npl)
    print *, "Triangles Constructed"
    !Construct the Cartesian Coordinates for Plotting
    do i=1,npoin_cg
        alon=coors(i,1)
        alat=coors(i,2)

        x=cos(alat)*cos(alon)
        y=cos(alat)*sin(alon)
        z=sin(alat)

        coord_cg(i,1)=x
        coord_cg(i,2)=y
        coord_cg(i,3)=z

    end do
    print *, "Cartesian Coords constructed"
    !!--------------------------------------------------------------!!
    !!--------------Construct the Quadrilateral Grid----------------!!
    !!--------------------------------------------------------------!!

    !Split the Triangles into 3 Quadrilaterals
    call make_quadl(coord_cg,intma,coordh,intmaq,npoin_cg,nelem,nside, &
        npoinh,nelemq,npoinq)
    print *, "Quads constructed"
    !Project all Points onto Sphere and store in spherical coordinates
    do i=1,npoinq
        x=coordh(i,1)
        y=coordh(i,2)
        z=coordh(i,3)
        r=sqrt( x*x + y*y + z*z )
        x=x/r
        y=y/r
        z=z/r
        alon=atan2(y,x)
        alat=asin(z)
        coorsh(i,1)=alon
        coorsh(i,2)=alat
    end do

    !Generate the Higher Order Points NOP in Spherical Coordinates
    call make_quadh(coorsh,intmaq,intmah,npoinq,npoinh,nelemq, &
        nsideh,nop,ngl)

    !Construct the Cartesian Coordinates for Plotting
    do i=1,npoinh
        alon=coorsh(i,1)
        alat=coorsh(i,2)

        x=cos(alat)*cos(alon)
        y=cos(alat)*sin(alon)
        z=sin(alat)

        coordh(i,1)=x
        coordh(i,2)=y
        coordh(i,3)=z
    end do

    write(*,'("---------------------------------------------------")')
    write(*,'(" nico nop npoin_cg nelem nside = ",5(i6,1x))')nico,nop,npoinh,nelemq,nsideh
    write(*,'("---------------------------------------------------")')


end subroutine ico_quad

!------------------------------------------------------------------------!
!>@brief This subroutine creates the initial Icosahedron.
!------------------------------------------------------------------------!
subroutine init_ico(coord_cg,intma,npoin_cg,nelem)

    use mod_constants, only:  tol, pi

    implicit none

    !global variables
    integer, intent(in) :: npoin_cg, nelem
    real, intent(out) :: coord_cg(npoin_cg,3)
    integer, intent(out) :: intma(nelem,3)

    !local variables
    real :: alon(12), alat(12), x, y, z, r
    integer :: i, ip, ie

    !construct north and south pole
    alon(1)=0.0
    alat(1)=pi/2.0
    alon(12)=0.0
    alat(12)=-pi/2.0

    !construct the north hemisphere nodes (5 more)
    do i=2,6
        alon(i)=real(2*i-3)/2.0*(2.0*pi/5.0)
        alat(i)=+2.0*asin( 1.0/(2.0*cos(3.0*pi/10.0)) ) - pi/2.0
    end do

    !construct the south hemisphere nodes (5 more)
    do i=7,11
        alon(i)=real(i-7)*(2.0*pi/5.0)
        alat(i)=-2.0*asin( 1.0/(2.0*cos(3.0*pi/10.0)) ) + pi/2.0
    end do

    !map spherical coords to Cartesian coords
    do ip=1,npoin_cg
        x=cos(alat(ip))*cos(alon(ip))
        y=cos(alat(ip))*sin(alon(ip))
        z=sin(alat(ip))
        if (abs(x) < tol) x=0
        if (abs(y) < tol) y=0
        if (abs(z) < tol) z=0
        coord_cg(ip,1)=x
        coord_cg(ip,2)=y
        coord_cg(ip,3)=z
        r=sqrt(coord_cg(ip,1)**2 + coord_cg(ip,2)**2 + coord_cg(ip,3)**2)
    end do

    !Construct Initial Connectivity
    !South
    ie=0
    ie=ie+1
    intma(ie,1)=12
    intma(ie,2)=8
    intma(ie,3)=7
    ie=ie+1
    intma(ie,1)=12
    intma(ie,2)=9
    intma(ie,3)=8
    ie=ie+1
    intma(ie,1)=12
    intma(ie,2)=10
    intma(ie,3)=9
    ie=ie+1
    intma(ie,1)=12
    intma(ie,2)=11
    intma(ie,3)=10
    ie=ie+1
    intma(ie,1)=12
    intma(ie,2)=7
    intma(ie,3)=11

    !Equator
    ie=ie+1
    intma(ie,1)=2
    intma(ie,2)=7
    intma(ie,3)=8
    ie=ie+1
    intma(ie,1)=2
    intma(ie,2)=8
    intma(ie,3)=3
    ie=ie+1
    intma(ie,1)=3
    intma(ie,2)=8
    intma(ie,3)=9
    ie=ie+1
    intma(ie,1)=3
    intma(ie,2)=9
    intma(ie,3)=4
    ie=ie+1
    intma(ie,1)=4
    intma(ie,2)=9
    intma(ie,3)=10
    ie=ie+1
    intma(ie,1)=4
    intma(ie,2)=10
    intma(ie,3)=5
    ie=ie+1
    intma(ie,1)=5
    intma(ie,2)=10
    intma(ie,3)=11
    ie=ie+1
    intma(ie,1)=5
    intma(ie,2)=11
    intma(ie,3)=6
    ie=ie+1
    intma(ie,1)=6
    intma(ie,2)=11
    intma(ie,3)=7
    ie=ie+1
    intma(ie,1)=6
    intma(ie,2)=7
    intma(ie,3)=2

    !North
    ie=ie+1
    intma(ie,1)=1
    intma(ie,2)=2
    intma(ie,3)=3
    ie=ie+1
    intma(ie,1)=1
    intma(ie,2)=3
    intma(ie,3)=4
    ie=ie+1
    intma(ie,1)=1
    intma(ie,2)=4
    intma(ie,3)=5
    ie=ie+1
    intma(ie,1)=1
    intma(ie,2)=5
    intma(ie,3)=6
    ie=ie+1
    intma(ie,1)=1
    intma(ie,2)=6
    intma(ie,3)=2

end subroutine init_ico

!------------------------------------------------------------------------!
!>@brief This subroutine creates the Linear Quadrilateral Grid.
!------------------------------------------------------------------------!
subroutine make_quadl(coord_cg,intma,coordh,intmaq,npoin_cg,nelem,nside,&
    npoinh,nelemh,npoinq)

    implicit none

    real, parameter :: etol=1.0e-10

    !global variables
    integer, intent(in) :: npoin_cg, nelem, nside, npoinh, nelemh
    real, intent(out) :: coordh(npoinh,3)
    integer, intent(out) :: intmaq(nelemh,4)
    integer, intent(out) :: npoinq
    real, intent(in) :: coord_cg(npoin_cg,3)
    integer, intent(in) :: intma(nelem,3)

    !local variables
    integer :: iside(nside,4), lwher(npoin_cg), lhowm(npoin_cg)
    integer :: icone(6*npoin_cg), jesid(nelem,3)
    integer :: ipoin(nside), inode(7)
    real :: dsmin, areat_total, area_min, area_max
    integer :: nelemq, i, j, is, i1, i2, ie, ip1, j1, j2
    real :: x1, y1, z1, x2, y2, z2, xm, ym, zm, x, y, z

    !constants
    dsmin=+1.0e5
    areat_total=0
    area_min=+1.0e5
    area_max=-1.0e5

    !initialize
    nelemq=0
    npoinq=npoin_cg
    coordh=0
    intmaq=0

    !initialize new coordinates
    do i=1,npoin_cg
        do j=1,3
            coordh(i,j)=coord_cg(i,j)
        end do !j
    end do !i

    !form side information
    call side_ico(nelem,npoin_cg,nside,intma,iside,lwher,lhowm,icone, &
        jesid,3)

    !obtain midpoints of sides
    do is=1,nside
        i1=iside(is,1)
        i2=iside(is,2)
        x1=coord_cg(i1,1)
        y1=coord_cg(i1,2)
        z1=coord_cg(i1,3)
        x2=coord_cg(i2,1)
        y2=coord_cg(i2,2)
        z2=coord_cg(i2,3)

        xm=(x1 + x2)/2.0
        ym=(y1 + y2)/2.0
        zm=(z1 + z2)/2.0

        npoinq=npoinq + 1

        if (abs(xm) < etol) xm=0
        if (abs(ym) < etol) ym=0
        if (abs(zm) < etol) zm=0

        coordh(npoinq,1)=xm
        coordh(npoinq,2)=ym
        coordh(npoinq,3)=zm

        ipoin(is)=npoinq
    end do !is

    !loop thru TRIANGULAR elements and create the new QUADRILATERAL ones
    do ie=1,nelem

        !store vertices
        do i=1,3
            inode(i)=intma(ie,i)
        end do !i

        !find the new point for the side
        do i=1,3
            i1=inode(i)
            ip1=i+1
            if (ip1 > 3) ip1=1
            i2=inode(ip1)

            do j=1,3
                is=jesid(ie,j)
                j1=iside(is,1)
                j2=iside(is,2)
                if ( (i1 == j1 .and. i2 == j2) .or. (i1 == j2 .and. i2 == j1) ) then
                inode(i+3)=ipoin(is)
                goto 10
            endif
        end do !j

        print*,' Error in make_mesh. no side found '
        print*,' ie i1 i2 i3 = ',ie,inode(1),inode(2),inode(3)
        stop
10  continue
    end do !i

    !form the centroid
    npoinq=npoinq + 1
    inode(7)=npoinq
    x=1.0/3.0*( coord_cg(inode(1),1) + coord_cg(inode(2),1) + coord_cg(inode(3),1) )
    y=1.0/3.0*( coord_cg(inode(1),2) + coord_cg(inode(2),2) + coord_cg(inode(3),2) )
    z=1.0/3.0*( coord_cg(inode(1),3) + coord_cg(inode(2),3) + coord_cg(inode(3),3) )

    if (abs(x) < etol) x=0
    if (abs(y) < etol) y=0
    if (abs(z) < etol) z=0

    coordh(npoinq,1)=x
    coordh(npoinq,2)=y
    coordh(npoinq,3)=z

    !form the 3 quadrilaterals
    !1st Quad
    intmaq(nelemq+1,1)=inode(1)
    intmaq(nelemq+1,2)=inode(4)
    intmaq(nelemq+1,3)=inode(7)
    intmaq(nelemq+1,4)=inode(6)

    !2nd Quad
    intmaq(nelemq+2,1)=inode(4)
    intmaq(nelemq+2,2)=inode(2)
    intmaq(nelemq+2,3)=inode(5)
    intmaq(nelemq+2,4)=inode(7)

    !3rd Quad
    intmaq(nelemq+3,1)=inode(6)
    intmaq(nelemq+3,2)=inode(7)
    intmaq(nelemq+3,3)=inode(5)
    intmaq(nelemq+3,4)=inode(3)

    !Update the Element Counter
    nelemq=nelemq + 3

end do !ie

if (nelemq /= nelemh) then
    print*,' Error in MAKE_QUADL. nelemq nelem = ',nelemq,nelemh
    stop
end if

end subroutine make_quadl
!------------------------------------------------------------------------!
!>@brief This subroutine creates the higher order quads
!------------------------------------------------------------------------!
subroutine make_quadh(coors,intmaq,intma,npoinq,npoin_cg,nelem,nside,nop,ngl)

    use mod_legendre, only : legendre_gauss_lobatto

    implicit none

    real, parameter :: etol=1.0e-10

    !global variables
    integer, intent(in) :: npoinq, npoin_cg, nelem, nside, nop, ngl
    real, intent(inout) :: coors(npoin_cg,2)
    integer, intent(out) :: intma(nelem,ngl,ngl)
    integer, intent(in) :: intmaq(nelem,4)

    !local variables
    integer :: iside(nside,4), lwher(npoinq), lhowm(npoinq)
    integer :: icone(6*npoinq), jesid(nelem,4)
    integer :: ipoin(nside,ngl)
    real :: lon(4), lat(4)
    real :: xgl(ngl), wgl(ngl), ksi_i(4), eta_j(4)
    data ksi_i/ -1, +1, +1, -1/
    data eta_j/ -1, -1, +1, +1/

    !Looping arrays
    real :: x(4), y(4)
    integer :: ipp(ngl,ngl), inode(4)
    integer :: ii(2,4), inn(4)
    data ii/ 1, 2, 4, 3, 1, 4, 2, 3 /

    integer :: npoinh
    integer :: i, j, is, iel, ier, ip, ie, in, in1, in2, i1, i2, jn, j1, j2
    integer :: ix, iy, k
    real :: xm, ym, zm, xt, yt, zt, alono, alato, alon, alat, alon_p, alat_p, r
    real :: dx, dy, eta, ksi, xsum, ysum, h_i, h_j
    real :: f_rot1, f_rot2, b_rot1, b_rot2
    real :: f_gnom1, f_gnom2, b_gnom1, b_gnom2
    
    !initialize
    npoinh=npoinq

    !Generat Gauss-Lobatto Points
    call legendre_gauss_lobatto(ngl,xgl,wgl)
    do i=1,ngl
        print*,'i xgl wgl = ',i,xgl(i),wgl(i)
    end do
    inn(1)=1
    inn(2)=ngl
    inn(3)=1
    inn(4)=ngl

    !form side information
    call side_ico(nelem,npoinq,nside,intmaq,iside,lwher,lhowm,icone,jesid,4)

    !obtain points at sides
    do is=1,nside
        do i=1,2
            inode(i)=iside(is,i)
        end do
        iel=iside(is,3)
        ier=iside(is,4)

        !Store Corner Nodes of Element IEL and get Midpoint
        xm=0
        ym=0
        zm=0
        do i=1,4
            ip=intmaq(iel,i)
            lon(i)=coors(ip,1)
            lat(i)=coors(ip,2)
            xt=cos(lat(i))*cos(lon(i))
            yt=cos(lat(i))*sin(lon(i))
            zt=sin(lat(i))
            xm=xm + xt
            ym=ym + yt
            zm=zm + zt
        end do !i
        xm=0.25*xm
        ym=0.25*ym
        zm=0.25*zm
        r=sqrt( xm*xm + ym*ym + zm*zm )
        xm=xm/r
        ym=ym/r
        zm=zm/r
        alono=atan2(ym,xm)
        alato=asin(zm)

        !Get Gnomonic coordinates of Side Nodes
        do i=1,2
            alon=coors(inode(i),1)
            alat=coors(inode(i),2)

            !Rotated coords
            alon_p=f_rot1(alon,alat,alono,alato)
            alat_p=f_rot2(alon,alat,alono,alato)

            !Gnomonic coords
            x(i)=f_gnom1(alon_p,alat_p)
            y(i)=f_gnom2(alon_p,alat_p)
        end do !i

        ipoin(is,1)=inode(1)
        ipoin(is,ngl)=inode(2)

        !construct Gauss-Lobatto Points in Gnomonic Space on this Side
        dx=x(2)-x(1)
        dy=y(2)-y(1)

        do i=2,ngl-1
            xt=dx/2.0*(xgl(i)+1.0) + x(1)
            yt=dy/2.0*(xgl(i)+1.0) + y(1)

            npoinh=npoinh + 1

            !Get Spherical coords
            alon_p=b_gnom1(xt)
            alat_p=b_gnom2(yt,alon_p)
            alon=b_rot1(alon_p,alat_p,alono,alato)
            alat=b_rot2(alon_p,alat_p,alono,alato)
            coors(npoinh,1)=alon
            coors(npoinh,2)=alat
            ipoin(is,i)=npoinh
        end do !i

    end do !is

    !loop thru elements and create the new ones
    do ie=1,nelem

        !Store Corner Nodes and get Midpoint
        xm=0
        ym=0
        zm=0
        do i=1,4
            ip=intmaq(ie,i)
            lon(i)=coors(ip,1)
            lat(i)=coors(ip,2)
            xt=cos(lat(i))*cos(lon(i))
            yt=cos(lat(i))*sin(lon(i))
            zt=sin(lat(i))
            xm=xm + xt
            ym=ym + yt
            zm=zm + zt
        end do
        xm=0.25*xm
        ym=0.25*ym
        zm=0.25*zm
        r=sqrt( xm*xm + ym*ym + zm*zm )
        xm=xm/r
        ym=ym/r
        zm=zm/r
        alono=atan2(ym,xm)
        alato=asin(zm)

        !Get Gnomonic coordinates of Corner Nodes of Element IE
        do i=1,4
            alon=lon(i)
            alat=lat(i)

            !Rotated coords
            alon_p=f_rot1(alon,alat,alono,alato)
            alat_p=f_rot2(alon,alat,alono,alato)

            !Gnomonic coords
            x(i)=f_gnom1(alon_p,alat_p)
            y(i)=f_gnom2(alon_p,alat_p)
        end do

        !Loop Twice to get to Opposite Edges for Sweeping 4 -- 3
        !
        !                                                 1 __ 2
        do in=1,4
            in1=ii(1,in)
            in2=ii(2,in)

            !Store Element Edges
            i1=intmaq(ie,in1)
            i2=intmaq(ie,in2)

            !Find which Side corresponds to this Edge
            do jn=1,4
                is=jesid(ie,jn)
                j1=iside(is,1)
                j2=iside(is,2)
                if (i1 == j1 .and. i2 == j2) then
                    do i=1,ngl
                        ip=ipoin(is,i)
                        if (in <= 2) then
                            ix=i
                            iy=inn(in)
                        else if (in > 2) then
                            ix=inn(in)
                            iy=i
                        endif

                        !Store Node
                        ipp(ix,iy)=ip
                    end do !i
                    goto 10
                else if (i1 == j2 .and. i2 == j1) then
                    do i=1,ngl
                        ip=ipoin(is,i)
                        if (in <= 2) then
                            ix=ngl-i+1
                            iy=inn(in)
                        else if (in > 2) then
                            ix=inn(in)
                            iy=ngl-i+1
                        endif

                        !Store Node
                        ipp(ix,iy)=ip
                    end do !i
                    goto 10
                endif
            end do !jn
            print*,' Error in Make_Quad. No Side found!>'
            stop
10      continue
        end do !in

        !Get Internal Gauss-Lobatto Points in Gnomonic Space
        do j=2,ngl-1
            eta=xgl(j)

            do i=2,ngl-1
                ksi=xgl(i)

                xsum=0
                ysum=0
                do k=1,4
                    h_i=0.5*( 1.0 + ksi*ksi_i(k) )
                    h_j=0.5*( 1.0 + eta*eta_j(k) )
                    xsum=xsum + h_i*h_j*x(k)
                    ysum=ysum + h_i*h_j*y(k)
                end do

                xt=xsum
                yt=ysum
                npoinh=npoinh + 1

                !Get Spherical coords
                alon_p=b_gnom1(xt)
                alat_p=b_gnom2(yt,alon_p)
                alon=b_rot1(alon_p,alat_p,alono,alato)
                alat=b_rot2(alon_p,alat_p,alono,alato)
                coors(npoinh,1)=alon
                coors(npoinh,2)=alat
                ipp(i,j)=npoinh
            end do !i
        end do !j

        !Store the Element
        do j=1,ngl
            do i=1,ngl
                intma(ie,i,j)=ipp(i,j)
            end do !i
        end do !j

    end do !IE

    !check grid constants
    if (npoin_cg /= npoinh) then
        print*,' Error in MAKE_QUADH. npoin_cg npoinh = ',npoin_cg,npoinh
        stop
    end if

end subroutine make_quadh
!------------------------------------------------------------------------!
!>@brief This subroutine creates the higher order triangles
!------------------------------------------------------------------------!
subroutine make_tris(coors,intma,intmah,npoin0,nelem0,nside0,npoin_cg,nelem,nop,ngl)

    implicit none

    real, parameter :: etol=1.0e-10

    !global variables
    integer, intent(in) :: npoin0, nelem0, nside0, npoin_cg, nelem, nop, ngl
    real, intent(inout) :: coors(npoin_cg,2)
    integer, intent(in) :: intma(nelem0,3)
    integer, intent(out) :: intmah(nelem,3)

    !local variables
    integer :: iside(nside0,4), lwher(npoin0), lhowm(npoin0)
    integer :: icone(6*npoin0), jesid(nelem0,3)
    integer :: ipoin(nside0,ngl)
    real :: lon(3), lat(3)
    real :: xgl(ngl)

    !Looping arrays
    real :: x(3), y(3)
    integer :: ipp(ngl,ngl), inode(3)
    integer :: ii(2,3)
    data ii/ 1, 2, 1, 3, 2, 3 /

    integer :: np, ne, i, is, iel, ier, ip, ie, in, in1, in2, i1, i2, jn
    integer :: j1, j2, ix, iy, j
    real :: dksi, xc, yc, zc, xt, yt, zt, r, alono, alato, alon, alat
    real :: rlon, rlat, dx, dy, ksi, eta, psi_1, psi_2, psi_3
    real :: f_rot1, f_rot2, b_rot1, b_rot2
    real :: f_gnom1, f_gnom2, b_gnom1, b_gnom2
    
    !Initialize
    np=npoin0
    ne=0

    !Construct Collocation Points
    dksi=1.0/nop
    do i=1,ngl
        xgl(i)=real(i-1)*dksi
    end do

    !form side information
    call side_ico(nelem0,npoin0,nside0,intma,iside,lwher,lhowm,icone,jesid,3)

    !obtain points at sides
    do is=1,nside0
        do i=1,2
            inode(i)=iside(is,i)
        end do
        iel=iside(is,3)
        ier=iside(is,4)

        !Store Vertices of Element IEL and get Centroid
        xc=0
        yc=0
        zc=0
        do i=1,3
            ip=intma(iel,i)
            lon(i)=coors(ip,1)
            lat(i)=coors(ip,2)
            xt=cos(lat(i))*cos(lon(i))
            yt=cos(lat(i))*sin(lon(i))
            zt=sin(lat(i))
            xc=xc + xt
            yc=yc + yt
            zc=zc + zt
        end do !i
        xc=xc/3.0
        yc=yc/3.0
        zc=zc/3.0
        r=sqrt( xc*xc + yc*yc + zc*zc )
        xc=xc/r
        yc=yc/r
        zc=zc/r
        alono=atan2(yc,xc)
        alato=asin(zc)

        !Get Gnomonic coordinates of Side Nodes
        do i=1,2
            alon=coors(inode(i),1)
            alat=coors(inode(i),2)

            !Rotated coords
            rlon=f_rot1(alon,alat,alono,alato)
            rlat=f_rot2(alon,alat,alono,alato)

            !Gnomonic coords
            x(i)=f_gnom1(rlon,rlat)
            y(i)=f_gnom2(rlon,rlat)
        end do !i

        ipoin(is,1)=inode(1)
        ipoin(is,ngl)=inode(2)

        !construct Points in Gnomonic Space on this Side
        dx=x(2)-x(1)
        dy=y(2)-y(1)

        do i=2,ngl-1
            xt=dx*xgl(i) + x(1)
            yt=dy*xgl(i) + y(1)

            np=np + 1

            !Get Spherical coords
            rlon=b_gnom1(xt)
            rlat=b_gnom2(yt,rlon)
            alon=b_rot1(rlon,rlat,alono,alato)
            alat=b_rot2(rlon,rlat,alono,alato)
            coors(np,1)=alon
            coors(np,2)=alat
            ipoin(is,i)=np
        end do !i

    end do !is

    !loop thru elements and create the new ones
    !Note NELEM0=NTREE because it is the original Icosahedron
    !and so INTMA=NTREE
    do ie=1,nelem0

        !Store Corner Nodes and get Midpoint
        xc=0
        yc=0
        zc=0
        do i=1,3
            ip=intma(ie,i)
            lon(i)=coors(ip,1)
            lat(i)=coors(ip,2)
            xt=cos(lat(i))*cos(lon(i))
            yt=cos(lat(i))*sin(lon(i))
            zt=sin(lat(i))
            xc=xc + xt
            yc=yc + yt
            zc=zc + zt
        end do
        xc=xc/3
        yc=yc/3
        zc=zc/3
        r=sqrt( xc*xc + yc*yc + zc*zc )
        xc=xc/r
        yc=yc/r
        zc=zc/r
        alono=atan2(yc,xc)
        alato=asin(zc)

        !Get Gnomonic coordinates of Vertices of Element IE
        do i=1,3
            alon=lon(i)
            alat=lat(i)

            !Rotated coords
            rlon=f_rot1(alon,alat,alono,alato)
            rlat=f_rot2(alon,alat,alono,alato)

            !Gnomonic coords
            x(i)=f_gnom1(rlon,rlat)
            y(i)=f_gnom2(rlon,rlat)
        end do

        !Loop Twice to get to Opposite Edges for Sweeping   3
        !                                                 /   \
        !                                                 1 _ 2
        do in=1,3
            in1=ii(1,in)
            in2=ii(2,in)

            !Store Element Edges
            i1=intma(ie,in1)
            i2=intma(ie,in2)

            !Find which Side corresponds to this Edge
            do jn=1,3
                is=jesid(ie,jn)
                j1=iside(is,1)
                j2=iside(is,2)
                if (i1 == j1 .and. i2 == j2) then
                    do i=1,ngl
                        ip=ipoin(is,i)
                        if (in == 1) then
                            ix=i
                            iy=1
                        else if (in == 2) then
                            ix=1
                            iy=i
                        else if (in == 3) then
                            ix=ngl+1-i
                            iy=i
                        endif

                        !Store Node
                        ipp(ix,iy)=ip
                    end do
                    goto 10
                else if (i1 == j2 .and. i2 == j1) then
                    do i=1,ngl
                        ip=ipoin(is,i)
                        if (in == 1) then
                            ix=ngl+1-i
                            iy=1
                        else if (in == 2) then
                            ix=1
                            iy=ngl+1-i
                        else if (in == 3) then
                            ix=i
                            iy=ngl+1-i
                        endif

                        !Store Node
                        ipp(ix,iy)=ip
                    end do
                    goto 10
                endif
            end do
            print*,' Error in Make_Tris. No Side found!>'
            stop
10      continue
        end do !in

        !Get Internal Points in Gnomonic Space
        do j=2,ngl-1
            eta=xgl(j)

            do i=2,ngl-j
                ksi=xgl(i)

                psi_1=1-ksi-eta
                psi_2=ksi
                psi_3=eta
                xt=psi_1*x(1) + psi_2*x(2) + psi_3*x(3)
                yt=psi_1*y(1) + psi_2*y(2) + psi_3*y(3)
                np=np + 1

                !Get Spherical coords
                rlon=b_gnom1(xt)
                rlat=b_gnom2(yt,rlon)
                alon=b_rot1(rlon,rlat,alono,alato)
                alat=b_rot2(rlon,rlat,alono,alato)
                coors(np,1)=alon
                coors(np,2)=alat
                ipp(i,j)=np
            end do !i
        end do !j

        !Store the Element
        do j=1,ngl-1
            j1=j
            j2=j+1
            do i=1,ngl-j
                i1=i
                i2=i+1

                !1st Triangle
                ne=ne+1
                intmah(ne,1)=ipp(i1,j1)
                intmah(ne,2)=ipp(i2,j1)
                intmah(ne,3)=ipp(i1,j2)

                !2nd Triangle
                if (i < ngl-j) then
                    ne=ne+1
                    intmah(ne,1)=ipp(i2,j2)
                    intmah(ne,2)=ipp(i1,j2)
                    intmah(ne,3)=ipp(i2,j1)
                end if
            end do
        end do

    end do !IE

    if (np /= npoin_cg) then
        print*,' Error in MAKE_TRIS.'
        print*,' np npoin_cg = ',np,npoin_cg
        stop
    end if
    if (ne /= nelem) then
        print*,' Error in MAKE_TRIS.'
        print*,' ne nelem = ',ne,nelem
        stop
    end if

end subroutine make_tris

!------------------------------------------------------------------------!
!>@brief This subroutine creates the array ISIDE which stores all of
!>the information concerning the sides of all the elements.
!------------------------------------------------------------------------!
subroutine side_ico(nelem,npoin_cg,nside,intma,iside,lwher,lhowm,icone,jesid,ncount)

    implicit none

    !global variables
    integer, intent(in) :: nelem, npoin_cg, nside, ncount
    integer, intent(in) :: intma(nelem,ncount)
    integer, intent(out) :: iside(nside,4), jesid(nelem,ncount)
    integer, intent(out) :: lwher(npoin_cg), lhowm(npoin_cg), icone(6*npoin_cg)

    !local variables
    integer :: in, ie, ip, iloca, jloca, iloc1, iele, iwher, ip1, ip2, iel, ier
    integer :: in1, in2, ipt, j, is, inode, jnode, i1, i2, jnod

    !initialize
    iside=0
    lwher=0
    lhowm=0
    icone=0
    jesid=0

    !count how many elements own each node
    do in=1,ncount
        do ie=1,nelem
            ip=intma(ie,in)
            lhowm(ip)=lhowm(ip) + 1
        end do
    end do

    !track elements owning each node
    lwher(1)=0
    do ip=2,npoin_cg
        lwher(ip)=lwher(ip-1) + lhowm(ip-1)
    end do

    !another tracker array
    lhowm=0

    do in=1,ncount
        do ie=1,nelem
            ip=intma(ie,in)
            lhowm(ip)=lhowm(ip) + 1
            jloca=lwher(ip) + lhowm(ip)
            icone(jloca)=ie
        end do
    end do

    !LOOP OVER THE NODES
    iloca=0
    do ip=1,npoin_cg
        iloc1=iloca
        iele=lhowm(ip)
        if (iele == 0 ) goto 1000
        iwher=lwher(ip)

        !LOOP OVER THOSE ELEMENTS SURROUNDING NODE IP

        ip1=ip
        do iel=1,iele
            ie=icone(iwher+iel)

            !find out position of ip in intma
            do in=1,ncount
                in1=in
                ipt=intma(ie,in)
                if (ipt == ip) goto 110
            end do !in
110     continue

        j=0
        do jnod=1,ncount-1,ncount-2
            j=j+1
            in2=in + jnod
            if (in2 > ncount) in2=in2-ncount
            ip2=intma(ie,in2)
            if (ip2 < ip1) goto 400

            !check whether side is old or new
            if (iloca == iloc1) goto 200
            do is=iloc1+1,iloca
                jloca=is
                if (iside(is,2) .eq. ip2) goto 250
            end do !is
200     continue

        !NEW SIDE
        iloca=iloca + 1
        iside(iloca,1)=ip1
        iside(iloca,2)=ip2
        iside(iloca,2+j)=ie
        goto 300

    !OLD SIDE
250 continue
    iside(jloca,2+j)=ie
300 continue

400 continue
    end do !jnod
end do !iel

!Perform some Shifting to order the nodes of a side in CCW direction
do is=iloc1+1,iloca
    if (iside(is,3) /= 0) goto 600
    iside(is,3)=iside(is,4)
    iside(is,4)=0
    iside(is,1)=iside(is,2)
    iside(is,2)=ip1
600 continue
    end do !is

!END LOOP OVER THE NODES
1000 continue
     end do !ip

     if (iloca /= nside) then
         print*,' Error in SIDE. iloca nside = ',iloca,nside
         stop
     end if

     !FORM ELEMENT/SIDE CONNECTIVITY ARRAY
     do is=1,nside
         iel=iside(is,3)
         ier=iside(is,4)
         inode=iside(is,1)
         jnode=iside(is,2)

         !LEFT SIDE
         do in=1,ncount
             i1=intma(iel,in)
             in1=in + 1
             if (in1 > ncount) in1=1
             i2=intma(iel,in1)
             if ((inode == i1) .and. (jnode == i2)) jesid(iel,in)=is
         end do !in

         !RIGHT SIDE
         if (ier <= 0) goto 2000
         do in=1,ncount
             i1=intma(ier,in)
             in1=in + 1
             if (in1 > ncount) in1=1
             i2=intma(ier,in1)
             if ((inode == i2) .and. (jnode == i1)) jesid(ier,in)=is
         end do !in

2000 continue
     end do !is

 end subroutine side_ico
