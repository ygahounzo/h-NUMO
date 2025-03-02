!---------------------------------------------------------------------!
! This subroutine reads in the mesh from .msh file gerated
! by GMSH. The mesh that is read in is 1st order. The read_gmsh_size
! subroutine reads only the size of the mesh.
! The high order populating of the elements is done a-posteriori 
! in the mod_grid part.
!
! Written by Michal Kopera on 08/22/2014
!            Naval Postgraduate School
!            Monterey, CA
!---------------------------------------------------------------------!
subroutine read_gmsh_size(fname,nelem,nnodes,nsides,nbounds,nbc)

    implicit none

    character :: fname*72,text*72
    integer, dimension(:,:),allocatable  :: elements

    integer :: nnodes,nelem,i,ip,j,nelements,nsides,nbounds,nbc

    character(len=80) :: line
    integer           :: success, indx, prev, beginning
    integer           :: value,e
    real              :: x,y,z

    print*,"reading mesh file:", fname

    open(1,file=fname)

    print*,"Reading NODES"
    do while (text.ne."$Nodes")
        read(1,*) text
    end do
    read(1,'(1I10)')nnodes

    do i=1,nnodes
        read(1,*)ip,x,y,z
    end do
    print*,"DONE Reading NODES: nnodes = ",nnodes

    print*,"Reading Elements"
    do while (text.ne."$Elements")
        read(1,*) text
    end do

    read(1,'(1I10)')nelements

    allocate(elements(nelements,9))

    nelem=0
    do e=1,nelements !read elments data - internal reads

        read(1,'(A)',iostat=success) line
        prev      = 1
        beginning = 1
        j=0
        do i=1,len(line)

            indx = index('0123456789', line(i:i))

            if (indx.eq.0 .and. prev.gt.0) then
                read(line(beginning:i-1), *) value
                j=j+1
                elements(e,j) = value
            else if (indx.gt.0 .and. prev.eq.0) then
                beginning = i
            end if
            prev = indx
        end do

        if (elements(e,2)==1) then
            nbounds=nbounds+1
        else if (elements(e,2)==3) then
            nelem=nelem+1
        end if
    end do
    nsides = ((nelem*4)-nbounds)/2+nbounds
    print*,"DONE Reading Elements: nelem nbounds nsides = ",nelem,nbounds,nsides

    print*,"Reading Boundaries"
    do while (text.ne."$BC")
        read(1,*) text
    end do
    read(1,*) nbc
    print*,"Done Reading Boundaries"

    close(1)

end subroutine read_gmsh_size

!---------------------------------------------------------------------!
subroutine read_gmsh(fname,coord0,intma0,nelem,nnodes,boundary,nboun,bc,nbc)

    implicit none

    real, dimension(2,nnodes)    :: coord0
    integer, dimension(nelem,4)  :: intma0
    integer, dimension(nboun,6)  :: boundary
    integer, dimension(2,nbc)   :: bc
    integer, dimension(:,:),allocatable  :: elements

    character :: fname*72,text*72

    integer :: nnodes,nelem,i,ip,j,ie,nelements,ib,nbc,nboun

    character(len=80) :: line
    integer           :: success, indx, prev, beginning
    integer           :: value,e
    real              :: x,y,z

    open(1,file=fname)

    print*,"Reading NODES"
    do while (text.ne."$Nodes")
        read(1,*) text
    end do
    read(1,*)text

    coord0=0.0
    do i=1,nnodes
        read(1,*)ip,x,y,z
        coord0(1,i)=x !correct that
        coord0(2,i)=y
    end do
    print*,"DONE Reading NODES"

    print*,"Reading Elements"
    do while (text.ne."$Elements")
        read(1,*) text
    end do
    read(1,*) nelements

    allocate(elements(nelements,9))

    intma0=0

    ie=0
    ib=0
    do e=1,nelements !read elments data - internal reads

        read(1,'(A)',iostat=success) line
        prev      = 1
        beginning = 1
        j=0
        do i=1,len(line)

            indx = index('0123456789', line(i:i))

            if (indx.eq.0 .and. prev.gt.0) then
                read(line(beginning:i-1), *) value
                j=j+1
                elements(e,j) = value
            else if (indx.gt.0 .and. prev.eq.0) then
                beginning = i
            end if
            prev = indx
        end do
        if(elements(e,2)==1) then
            ib=ib+1
            boundary(ib,1:6)=elements(e,2:7)
        else if(elements(e,2)==3) then
            ie=ie+1
            intma0(ie,:)=elements(e,6:9)
        end if
    end do
    print*,"DONE Reading Elements"

    do while (text.ne."$BC")
        read(1,*) text
    end do
    read(1,*) nbc
    do e=1,nbc
        read(1,*)bc(:,e)
    end do

end subroutine read_gmsh

subroutine read_bathy(fname,bathy0,nnodes)
    !
    ! S. Marras, April 2016
    !
    implicit none

    real, dimension(nnodes) :: bathy0
  
    character :: fname*72,text*72
    integer   :: nnodes,i,ip,j
    real      :: z

    open(1,file=fname)

    print*,"Reading Bathymetry (point-wise from linear grid)"
    do while (text.ne."$Bathy")
        read(1,*) text
    end do
    read(1,*)text

    bathy0=0.0
    do i=1,nnodes
        read(1,*)ip, z
        bathy0(i)=z
    end do
    print*,"DONE Reading Bathymetry data"

    close(1)

end subroutine read_bathy

!---------------------------------------------------------------------!
!program test_gmsh
!
!  implicit none
!
!
!  real, dimension(:,:), allocatable :: coord0
!  integer, dimension(:,:), allocatable :: intma0

!   integer nelem,nnodes,e

!   character :: fname*72

!   fname="untitled2.msh"

!   call read_gmsh_size(fname,nelem,nnodes)

!   allocate(coord0(2,nnodes),intma0(nelem,4))
!   coord0=0
!   intma0=0

!   print*,"nelem=",nelem," nnodes=",nnodes

!   call read_gmsh(fname,coord0,intma0,nelem,nnodes)

!   do e=1,nelem
!      print*,e,intma0(e,1:4)
!   end do

!   print*,"nelem=",nelem," nnodes=",nnodes

!   do e=1,nnodes
!      print*,e,coord0(1:2,e)
!   end do

! end program test_gmsh


!------------------------------------------------------------------------!
!Populate linear quads with high order elements
!
! This subroutine adds the LGL points
! to a linear grid of quads on the plane.
!
! Written by F.X. Giraldo 11/00 for the sphere
!
! Modified by:
! Simone Marras to populate, with LGL points, a linear grid that is read
! from an external file
!
! The routine to read the grid is called read_mesh.
! October 2012, Cambridge, UK
!------------------------------------------------------------------------!

!------------------------------------------------------------------------!
!This subroutine creates the higher order quads
!------------------------------------------------------------------------!
subroutine make_quadh_gmsh(coordh,intmaq,intma,bsido,iboun,npoinq,npoin,nelem,nboun,boundary,nside,nop,ngl,npoin_out,bc,nbc)

    use mod_constants, only: tol
  
    use mod_legendre, only: legendre_gauss_lobatto

    implicit none

    integer npoin, nelem, nboun, nside, npoinh, nelemh, nelemq, npoinq, nop,ngl
    integer i,j,k,ip,i1,i2,ir,is,ie,ip1,ip2,j1,j2,iel,ier,in,in1,in2,jn
    integer ix,iy,iz
    integer npoin_out
    integer nbc, bc(2,nbc) !boundary condition types

    !global arrays
    real coordh(2,npoin)
    integer intmaq(nelem,4), intma(ngl,ngl,nelem), bsido(4,nboun), iboun(4), boundary(nboun,6)

    !local arrays
    integer iside(nside,4), jeside(nelem,4)
    integer ipoin(nside,ngl)
    real xgl(ngl), wgl(ngl), ksi_i(4), eta_j(4)
    real ksi, eta
    real h_i,h_j
    data ksi_i/ -1, +1, +1, -1/
    data eta_j/ -1, -1, +1, +1/

    !Looping arrays
    real x(4), y(4)
    integer ipp(ngl,ngl), inode(4)

    real dsmin, areat_total, area_min, area_max
    real x1, x2, y1, y2, z1, z2
    real xm,ym,zm,xt,yt,zt
    real dx, dy, dz
    real r
    real xsum,ysum,zsum
    real xmax, xmin, ymax, ymin
    real s1, s2
    integer ib, ib1, ib2, ibc,ic

    integer btype

    !initialize
    npoinh=npoinq


    !Generat Gauss-Lobatto Points
    call legendre_gauss_lobatto(ngl,xgl,wgl)
    do i=1,ngl
        print*,'i xgl wgl = ',i,xgl(i),wgl(i)
    end do

    !form side information
    call side_gmsh(iside,jeside,intmaq,bsido,npoinq,nelem,nboun,nside,ngl)

    intma=0

    ! first put corner nodes into elements
    do ie=1,nelem
        intma(1,1,ie)=intmaq(ie,1)
        intma(ngl,1,ie)=intmaq(ie,2)
        intma(ngl,ngl,ie)=intmaq(ie,3)
        intma(1,ngl,ie)=intmaq(ie,4)
    end do

    !then go over sides and create lgl points
    !obtain points at sides
    do is=1,nside
        do i=1,2
            inode(i)=iside(is,i)
        end do
        iel=iside(is,3) !element on the left  side of iside
        ier=iside(is,4) !element on the right side of iside

        !Get coordinates of Side Nodes
        do i=1,2
            xt=coordh(1,inode(i))
            yt=coordh(2,inode(i))

            !Gnomonic coords
            x(i)=xt
            y(i)=yt
        end do !i

        ipoin(is,1)   = inode(1)
        ipoin(is,ngl) = inode(2)

        !construct Gauss-Lobatto Points in Gnomonic Space on this Side
        dx=x(2)-x(1)
        dy=y(2)-y(1)

        do i=2,ngl-1
            xt=dx/2.0*(xgl(i)+1.0) + x(1)
            yt=dy/2.0*(xgl(i)+1.0) + y(1)

            npoinh=npoinh + 1

            coordh(1,npoinh)=xt
            coordh(2,npoinh)=yt
            ipoin(is,i)=npoinh !ipoin stores points at each side
        end do !i

    end do !is

    ic = npoinh !point counter

    !loop thru elements and create the interior points
    do ie=1,nelem

        do i=1,4
            ip=intmaq(ie,i)
            xt=coordh(1,ip)
            yt=coordh(2,ip)
            x(i)=xt
            y(i)=yt
        end do

        !corners
        intma(1,1,ie)=intmaq(ie,1)
        intma(ngl,1,ie)=intmaq(ie,2)
        intma(ngl,ngl,ie)=intmaq(ie,3)
        intma(1,ngl,ie)=intmaq(ie,4)

        !get side information from ipoin
        !south
        is=jeside(ie,1)
        ip1=iside(is,1)
        ip2=iside(is,2)

        if (ip1==intma(1,1,ie)) then
            intma(2:ngl-1,1,ie)=ipoin(is,2:ngl-1)
        else
            !inverse order
            do i=2,ngl-1
                intma(i,1,ie)=ipoin(is,ngl+1-i)
            end do
        end if

        !north
        is=jeside(ie,2)
        ip1=iside(is,1)
        ip2=iside(is,2)

        if (ip1==intma(1,ngl,ie)) then
            intma(2:ngl-1,ngl,ie)=ipoin(is,2:ngl-1)
        else
            !inverse order
            do i=2,ngl-1
                intma(i,ngl,ie)=ipoin(is,ngl+1-i)
            end do
        end if

        !west
        is=jeside(ie,3)
        ip1=iside(is,1)
        ip2=iside(is,2)

        if (ip1==intma(1,1,ie)) then
            intma(1,2:ngl-1,ie)=ipoin(is,2:ngl-1)
        else
            !inverse order
            do i=2,ngl-1
                intma(1,i,ie)=ipoin(is,ngl+1-i)
            end do
        end if

        !east
        is=jeside(ie,4)
        ip1=iside(is,1)
        ip2=iside(is,2)

        if (ip1==intma(ngl,1,ie)) then
            intma(ngl,2:ngl-1,ie)=ipoin(is,2:ngl-1)
        else
            !inverse order
            do i=2,ngl-1
                intma(ngl,i,ie)=ipoin(is,ngl+1-i)
            end do
        end if

        !interior
        do j=2,ngl-1
            do i=2,ngl-1
                ic=ic+1
                intma(i,j,ie)=ic
                call linear_interp(coordh(1,ic),xgl(i),xgl(j),x)
                call linear_interp(coordh(2,ic),xgl(i),xgl(j),y)
            end do
        end do


    end do !IE

    npoin_out=ic

    !Use ISIDE to construct BSIDO
    xmax=maxval(coordh(1,:)); xmin=minval(coordh(1,:))
    ymax=maxval(coordh(2,:)); ymin=minval(coordh(2,:))
    ib=0

    bsido=0
    ib=0
    do is=1,nside
        if (iside(is,4)==0) then
            ib=ib+1
            bsido(1:3,ib)=iside(is,1:3)

            do j=1,nboun
                if(((bsido(1,ib)==boundary(j,5)).and.(bsido(2,ib)==boundary(j,6)))&
                    .or.(((bsido(1,ib)==boundary(j,6)).and.(bsido(2,ib)==boundary(j,5))))) then
                    btype=boundary(j,3)
                    do k=1,nbc
                        if(btype==bc(1,k)) then
                            bsido(4,ib)=bc(2,k)
                        end if
                    end do
                end if
            end do
           !        write(*,'(5I7)')ib,bsido(:,ib)
           !        bsido(4,ib)=4 !assume for now only one bc
        end if

    end do

    if (ib.ne.nboun) then
        print*,"ALL BOUNDARIES NOT ACCOUNTED FOR!!!"
        stop
    end if
    return

end subroutine make_quadh_gmsh

!!! SM
!------------------------------------------------------------------------!
!This subroutine creates the higher order quads
!------------------------------------------------------------------------!
subroutine make_bathyh(bathyh, intmaq,intma,bsido,iboun, &
    npoinq,npoin,nelem,nboun,boundary,nside,nop,ngl,npoin_out,bc,nbc)
  
    use mod_constants, only: tol
  
    use mod_legendre, only: legendre_gauss_lobatto

    implicit none

    integer npoin, nelem, nboun, nside, npoinh, nelemh, nelemq, npoinq, nop,ngl
    integer i,j,k,ip,i1,i2,ir,is,ie,ip1,ip2,j1,j2,iel,ier,in,in1,in2,jn
    integer ix,iy,iz
    integer npoin_out
    integer nbc, bc(2,nbc) !boundary condition types

    !global arrays
    real bathyh(npoin)
    integer intmaq(nelem,4), intma(ngl,ngl,nelem), bsido(4,nboun), iboun(4), boundary(nboun,6)

    !local arrays
    integer iside(nside,4), jeside(nelem,4)
    integer ipoin(nside,ngl)
    real xgl(ngl), wgl(ngl), ksi_i(4), eta_j(4)
    real ksi, eta
    real h_i,h_j
    data ksi_i/ -1, +1, +1, -1/
    data eta_j/ -1, -1, +1, +1/

    !Looping arrays
    real x(4), y(4), z(4)
    integer ipp(ngl,ngl), inode(4)

    real dsmin, areat_total, area_min, area_max
    real x1, x2, y1, y2, z1, z2
    real xm,ym,zm,xt,yt,zt
    real dx, dy, dz
    real r
    real xsum,ysum,zsum
    real xmax, xmin, ymax, ymin
    real s1, s2
    integer ib, ib1, ib2, ibc,ic

    integer btype

    character fname*100
  
    !initialize
    npoinh=npoinq

    !Generat Gauss-Lobatto Points
    call legendre_gauss_lobatto(ngl,xgl,wgl)
    do i=1,ngl
        print*,'i xgl wgl = ',i,xgl(i),wgl(i)
    end do

    !form side information
    call side_gmsh(iside,jeside,intmaq,bsido,npoinq,nelem,nboun,nside,ngl)

    intma=0

    ! first put corner nodes into elements
    do ie=1,nelem
        intma(1,1,ie)=intmaq(ie,1)
        intma(ngl,1,ie)=intmaq(ie,2)
        intma(ngl,ngl,ie)=intmaq(ie,3)
        intma(1,ngl,ie)=intmaq(ie,4)
    end do

    !then go over sides and create lgl points
    !obtain points at sides
    do is=1,nside
        do i=1,2
            inode(i)=iside(is,i)
        end do
        iel=iside(is,3) !element on the left  side of iside
        ier=iside(is,4) !element on the right side of iside

        !Get coordinates of Side Nodes
        do i=1,2
            zt=bathyh(  inode(i))

            !Gnomonic coords
            z(i)=zt
        end do !i

        ipoin(is,1)   = inode(1)
        ipoin(is,ngl) = inode(2)

        !construct Gauss-Lobatto Points in Gnomonic Space on this Side
        dz=z(2)-z(1)

        do i=2,ngl-1
            zt=dz/2.0*(xgl(i)+1.0) + z(1)
        
            npoinh=npoinh + 1

            !coordh(1,npoinh)=xt
            !coordh(2,npoinh)=yt
            bathyh(  npoinh)=zt
        
            ipoin(is,i)=npoinh !ipoin stores points at each side
        end do !i

    end do !is
  
    ic = npoinh !point counter

    !loop thru elements and create the interior points
    do ie=1,nelem

        do i=1,4
            ip=intmaq(ie,i)
            zt=bathyh(  ip)
            z(i)=zt
        end do

        !corners
        intma(1,1,ie)=intmaq(ie,1)
        intma(ngl,1,ie)=intmaq(ie,2)
        intma(ngl,ngl,ie)=intmaq(ie,3)
        intma(1,ngl,ie)=intmaq(ie,4)

        !get side information from ipoin
        !south
        is=jeside(ie,1)
        ip1=iside(is,1)
        ip2=iside(is,2)

        if (ip1==intma(1,1,ie)) then
            intma(2:ngl-1,1,ie)=ipoin(is,2:ngl-1)
        else
            !inverse order
            do i=2,ngl-1
                intma(i,1,ie)=ipoin(is,ngl+1-i)
            end do
        end if

        !north
        is=jeside(ie,2)
        ip1=iside(is,1)
        ip2=iside(is,2)

        if (ip1==intma(1,ngl,ie)) then
            intma(2:ngl-1,ngl,ie)=ipoin(is,2:ngl-1)
        else
            !inverse order
            do i=2,ngl-1
                intma(i,ngl,ie)=ipoin(is,ngl+1-i)
            end do
        end if

        !west
        is=jeside(ie,3)
        ip1=iside(is,1)
        ip2=iside(is,2)

        if (ip1==intma(1,1,ie)) then
            intma(1,2:ngl-1,ie)=ipoin(is,2:ngl-1)
        else
            !inverse order
            do i=2,ngl-1
                intma(1,i,ie)=ipoin(is,ngl+1-i)
            end do
        end if

        !east
        is=jeside(ie,4)
        ip1=iside(is,1)
        ip2=iside(is,2)

        if (ip1==intma(ngl,1,ie)) then
            intma(ngl,2:ngl-1,ie)=ipoin(is,2:ngl-1)
        else
            !inverse order
            do i=2,ngl-1
                intma(ngl,i,ie)=ipoin(is,ngl+1-i)
            end do
        end if

        !interior
        do j=2,ngl-1
            do i=2,ngl-1
                ic=ic+1
                intma(i,j,ie)=ic
                call linear_interp(bathyh(  ic),xgl(i),xgl(j),z)
            end do
        end do

    end do !IE

    npoin_out=ic
  
end subroutine make_bathyh

!Linear interpolation subroutine
subroutine linear_interp(coord,ksi,eta,x)

    implicit none

    real :: coord, ksi, eta
    real, dimension(4) :: x

    coord=0.25*((1-ksi)*(1-eta)*x(1)&
        +(1+ksi)*(1-eta)*x(2)&
        +(1+ksi)*(1+eta)*x(3)&
        +(1-ksi)*(1+eta)*x(4))

end subroutine linear_interp

!----------------------------------------------------------------------!
!This subroutine creates the array ISIDE which stores all of
!the information concerning the sides of all the elements.
!Written by Francis X. Giraldo on 1/01
!           Naval Postgraduate School
!           Department of Applied Mathematics
!           Monterey, CA 93943-5216
!----------------------------------------------------------------------!
subroutine side_gmsh(iside,jeside,intma,bsido,npoin,nelem,nboun,nside,ngl)

    implicit none

    !global arrays
    integer iside(nside,4), jeside(nelem,4)
    integer intma(nelem,4), bsido(4,nboun)
    integer npoin, nelem, nboun, nside, ngl

    !local arrays
    integer lwher(npoin), lhowm(npoin), icone(6*npoin)
    integer temp1, temp2, temp3, temp4
    integer in, ie, ip, e, jloca, iloca, iloc1, iele, iwher, ip1, ip2, ipt
    integer iel, ier, in1, in2, i, i1, i2, j, jnod, iold, is, is1, is2, il, ir
    integer ib, ibe, ibc, ilb, irb

    !Swap Nodes to make CCW
    !  do ie=1,nelem
    !     itemp=intma(ie,3)
    !     intma(ie,3)=intma(ie,4)
    !     intma(ie,4)=itemp
    !  end do

    !initialize
    iside=0
    lwher=0
    lhowm=0
    icone=0
    jeside=0

    !count how many elements own each node
    do in=1,4
        do ie=1,nelem
            ip=intma(ie,in)
            lhowm(ip)=lhowm(ip) + 1
        end do
    end do

    !track elements owning each node
    lwher(1)=0
    do ip=2,npoin
        lwher(ip)=lwher(ip-1) + lhowm(ip-1)
    end do

    !another tracker array
    lhowm=0
    do in=1,4
        do ie=1,nelem
            ip=intma(ie,in)
            lhowm(ip)=lhowm(ip) + 1
            jloca=lwher(ip) + lhowm(ip)
            icone(jloca)=ie
        end do
    end do

    !LOOP OVER THE NODES
    iloca=0
    do ip=1,npoin
        iloc1=iloca
        iele=lhowm(ip)
        if (iele /= 0 ) then
            iwher=lwher(ip)

            !LOOP OVER THOSE ELEMENTS SURROUNDING NODE IP

            ip1=ip
            do iel=1,iele
                ie=icone(iwher+iel)

                !find out position of ip in intma
                do in=1,4
                    in1=in
                    ipt=intma(ie,in)
                    if (ipt == ip) exit
                end do           !in

                j=0
                do jnod=1,3,2
                    iold=0
                    j=j+1
                    in2=in + jnod
                    if (in2 > 4) in2=in2-4
                    ip2=intma(ie,in2)
                    if (ip2 >= ip1) then

                        !check whether side is old or new
                        if (iloca /= iloc1) then
                            do is=iloc1+1,iloca
                                jloca=is
                                if (iside(is,2) == ip2) then
                                    iold=1
                                    exit
                                end if !iside
                            end do  !is
                        end if !iloca

                        if (iold == 0) then
                            !NEW SIDE
                            iloca=iloca + 1
                            iside(iloca,1)=ip1
                            iside(iloca,2)=ip2
                            iside(iloca,2+j)=ie
                        else if (iold == 1) then
                            !OLD SIDE
                            iside(jloca,2+j)=ie
                        end if !iold
                    end if !ip2
                end do           !jnod
            end do              !iel

            !Perform some Shifting to order the nodes of a side in CCW direction
            do is=iloc1+1,iloca
                if (iside(is,3) == 0) then
                    iside(is,3)=iside(is,4)
                    iside(is,4)=0
                    iside(is,1)=iside(is,2)
                    iside(is,2)=ip1
                end if !iside
            end do              !is

           !END LOOP OVER THE NODES
        end if !iele
    end do !ip

    if (iloca /= nside) then
        print*,' Error in SIDE. iloca nside = ',iloca,nside
        stop
    end if

    !RESET THE BOUNDARY MARKERS
    goto 10
    do is=1,nside
        if (iside(is,4) == 0) then
            il=iside(is,1)
            ir=iside(is,2)
            ie=iside(is,3)
            do ib=1,nboun
                ibe=bsido(3,ib)
                ibc=bsido(4,ib)
                if (ibe == ie) then
                    ilb=bsido(1,ib)
                    irb=bsido(2,ib)
                    if ( (ilb == il) .and. (irb == ir) ) then
                        iside(is,4)=-ibc
                        exit
                    end if
                end if
            end do              !ib
        end if !iside
    end do                 !is
10 continue

   !FORM ELEMENT/SIDE CONNECTIVITY ARRAY
   do is=1,nside
       iel=iside(is,3)
       ier=iside(is,4)
       is1=iside(is,1)
       is2=iside(is,2)

       !LEFT SIDE
       do in=1,4
           i1=intma(iel,in)
           in1=in + 1
           if (in1 > 4) in1=1
           i2=intma(iel,in1)
           if ((is1 == i1) .and. (is2 == i2)) jeside(iel,in)=is
       end do !in

       !RIGHT SIDE
       if (ier > 0) then
           do in=1,4
               i1=intma(ier,in)
               in1=in + 1
               if (in1 > 4) in1=1
               i2=intma(ier,in1)
               if ((is1 == i2) .and. (is2 == i1)) jeside(ier,in)=is
           end do              !in
       end if !ier
   end do !is

   !swap local face numbering in jeside to assure S,N,W,E order
   do e=1,nelem
       temp2=jeside(e,2)
       temp3=jeside(e,3)
       temp4=jeside(e,4)
       jeside(e,2)=temp3
       jeside(e,3)=temp4
       jeside(e,4)=temp2
   end do


   end subroutine side_gmsh

   subroutine read_bathymetry_structured_size(fname,Ox,Oy,Dx,Dy,nx,ny)

       implicit none

       character :: fname*72
       real    :: Ox, Oy, Oz !origin coordinates
       real    :: Dx, Dy, Dz !increments
       integer :: nx, ny, nz !number of points in each direction

       open(1,file=fname)

       read(1,*) Ox, Oy, Oz
       read(1,*) Dx, Dy, Dz
       read(1,*) nx, ny, nz

       !  print*,"bathymetry read from file"
       !  print*, Ox, Oy, Oz
       !  print*, Dx, Dy, Dz
       !  print*, nx, ny, nz

       close(1)

   end subroutine read_bathymetry_structured_size


   subroutine read_bathymetry_structured_data(fname,hb_struct,nx,ny)

       implicit none

       character :: fname*72
       real :: Ox, Oy, Oz !origin coordinates
       real    :: Dx, Dy, Dz !increments
       integer :: nx, ny, nz !number of points in each direction
       integer :: i,j,k

       real :: x,y,z,hb
       real, dimension(nx,ny) :: hb_struct

       open(1,file=fname)

       read(1,*) Ox, Oy, Oz
       read(1,*) Dx, Dy, Dz
       read(1,*) nx, ny, nz

       do i=1,nx
           do j=1,ny
               read(1,*)hb
               hb_struct(i,j)=hb
           end do
       end do

       close(1)
   end subroutine read_bathymetry_structured_data

   subroutine interpolate_bathymetry(hb,xb,yb,hb_struct, Ox,Oy,Dx,Dy,nx,ny)
       !get hb value for given xb,yb and hb_struct bathymetry field

       implicit none

       real, dimension(nx,ny) :: hb_struct
       real :: Ox,Oy,Dx,Dy,xb,yb,hb
       integer:: nx,ny

       integer :: i,j
       real,dimension(2)::x,y
       real :: ddx,ddy,a,b,c,d

       i = floor((xb-Ox)/Dx)+1
       j = floor((yb-Oy)/Dy)+1


       if(i==0) then
           i=1
       end if

       if(j==0) then
           j=1
       end if

       if(i>=nx) then
           i=nx-1
       end if

       if(j>=ny) then
           j=ny-1
       end if

       x(1) = Ox+(i-1)*Dx
       x(2) = x(1)+Dx
       y(1) = Oy+(j-1)*Dy
       y(2) = y(1)+Dy

       ddx = x(2)-x(1)
       ddy = y(2)-y(1)

       a = xb-x(1)
       b = x(2)-xb
       c = yb-y(1)
       d = y(2)-yb

       hb = b*d*hb_struct(i,j)
       hb = hb + a*d*hb_struct(i+1,j)
       hb = hb + a*c*hb_struct(i+1,j+1)
       hb = hb + b*c*hb_struct(i,j+1)
       hb = hb/(ddx*ddy)

   end subroutine interpolate_bathymetry

   !read high order mesh
   subroutine read_ho_mesh_size(nelem,ngl,npoin,nside,nboun)

       use mod_input, only: mesh_file
  
       use mod_basis, only: nop

       implicit none

       integer:: nelem, ngl, npoin, nside, nboun

       integer:: i
       real::x,y,z

       open(1,file=mesh_file)
       read(1,*)npoin
       do i=1,npoin
           read(1,*)x,y,z
       end do
       read(1,*)nelem,ngl,nside,nboun
       close(1)

       if(ngl.ne.nop+1) then
           print*,"READ_HO_MESH_SIZE:: ngl in the mesh file does not match the input file"
           stop
       end if

   end subroutine read_ho_mesh_size

   subroutine read_ho_mesh(coordg,npoinT,intma,nelemT,ngl)

       use mod_input, only: mesh_file

       implicit none

       integer :: npoinT,nelemT,ngl
       real, dimension(2,npoinT) :: coordg
       integer, dimension(ngl,ngl,nelemT) :: intma

       integer :: i,j,el,ii,jj,ee,nelem,nside,npoin
       real :: buf
       real :: x,y,z

       open(1,file=mesh_file)
       read(1,*)npoin
       do i=1,npoin
           read(1,*)x,y,z
           coordg(1,i)=x
           coordg(2,i)=y
       end do

       read(1,*)nelem,i,j,ii
       do el=1,nelem
           do j=1,ngl
               do i=1,ngl
                   read(1,*)ee,ii,jj,intma(ii,jj,ee)
               end do
           end do
       end do
       close(1)
       print*,'circular coords read'

   end subroutine read_ho_mesh
