!----------------------------------------------------------------------!
!>This module builds the Faces and Normals
!>@author  J.F. Kelly on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@date 11 October 2009
!> Creates normals for spectral element hexahedra
!>@date 16 November 2009
!>@date December 2014, Daniel S. Abdi, Fixed bug in jac_face
!> and modified for inter-processor faces
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: April 03, 2023
!       remove SEMI_ANALYTIC_METRICS and is_sphere flag
!----------------------------------------------------------------------!
subroutine create_normals(nv,jac_face,face,nface)

    use mod_basis, only: ngl, nglx, ngly, nglz, wglx, wgly, wglz, FACE_LEN

    use mod_grid, only: intma, coord

    use mod_input, only: geometry_type

    use mod_interface, only: compute_local_gradient_v3

    implicit none

    !global
    integer, intent(in) :: nface
    real, dimension(3,ngl,ngl,nface), intent(out) :: nv
    real, dimension(ngl,ngl,nface), intent(out) :: jac_face
    integer, intent(in) :: face(FACE_LEN,nface)

    !local
    integer :: iface, i, j, k, ip, l, m
    integer :: ilocl, ilocr, iel, ier, ndim
    real :: ww, nx, ny, nz, nlen

    !local arrays
    real, dimension(nglx,ngly,nglz) :: x, y, z
    real, dimension(nglx,ngly,nglz) :: x_ksi, x_eta, x_zeta, &
        y_ksi, y_eta, y_zeta, &
        z_ksi, z_eta, z_zeta
  
    !Initialize MXM dimension
    ndim=1

    !initialize the global matrix
    nv=0
    jac_face=0

    !loop thru the faces
    do iface=1,nface

        !Store Left and Right Elements
        ilocl=face(5,iface)
        ilocr=face(6,iface)
        iel= face(7,iface)
        ier= face(8,iface)

        if(iel.ne.0) then
        
            !Store Element Variables
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        ip=intma(i,j,k,iel)
                        x(i,j,k)=coord(1,ip)
                        y(i,j,k)=coord(2,ip)
                        z(i,j,k)=coord(3,ip)
                    end do !i
                end do !j
            end do !k
        
            !Construct Mapping Derivatives
            call compute_local_gradient_v3(x_ksi,x_eta,x_zeta,x,nglx,ngly,nglz,ndim)
            call compute_local_gradient_v3(y_ksi,y_eta,y_zeta,y,nglx,ngly,nglz,ndim)
            call compute_local_gradient_v3(z_ksi,z_eta,z_zeta,z,nglx,ngly,nglz,ndim)
           
            !set jacobian matrix J for 2D

            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        if(nglx == 1) then
                            x_ksi(i,j,k) = 1;  y_ksi(i,j,k) = 0;  z_ksi(i,j,k) = 0;
                            x_eta(i,j,k) = 0;  x_zeta(i,j,k) = 0;
                        endif
                        if(ngly == 1) then
                            x_eta(i,j,k) = 0;  y_eta(i,j,k) = 1;  z_eta(i,j,k) = 0;
                            y_ksi(i,j,k) = 0; y_zeta(i,j,k) = 0;
                        endif
                        if(nglz == 1) then
                            x_zeta(i,j,k) = 0; y_zeta(i,j,k) = 0; z_zeta(i,j,k) = 1;
                            z_ksi(i,j,k) = 0; z_eta(i,j,k) = 0;
                        endif
                    enddo
                enddo
            enddo
              
            !Compute Normals
            select case (ilocl)
           
                case(1)
                    ! Face 1: zeta = -1
                    do l=1,nglx
                        do m = 1,ngly
                            ww=wglx(l)*wgly(m)
                            i = l
                            j = m
                            k = 1
                            nx = + y_eta(i,j,k)*z_ksi(i,j,k) - z_eta(i,j,k)*y_ksi(i,j,k)
                            ny = - x_eta(i,j,k)*z_ksi(i,j,k) + z_eta(i,j,k)*x_ksi(i,j,k)
                            nz = + x_eta(i,j,k)*y_ksi(i,j,k) - y_eta(i,j,k)*x_ksi(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            jac_face(l,m,iface) = ww*nlen
                            nv(1,l,m,iface) = nx/nlen
                            nv(2,l,m,iface) = ny/nlen
                            nv(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l
           
                case(2)
                    ! Face 2: zeta = +1
                    do l=1,nglx
                        do m=1,ngly
                            ww=wglx(l)*wgly(m)
                            i = l
                            j = m
                            k = nglz
                            nx = - y_eta(i,j,k)*z_ksi(i,j,k) + z_eta(i,j,k)*y_ksi(i,j,k)
                            ny = + x_eta(i,j,k)*z_ksi(i,j,k) - z_eta(i,j,k)*x_ksi(i,j,k)
                            nz = - x_eta(i,j,k)*y_ksi(i,j,k) + y_eta(i,j,k)*x_ksi(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            jac_face(l,m,iface) = ww*nlen
                            nv(1,l,m,iface) = nx/nlen
                            nv(2,l,m,iface) = ny/nlen
                            nv(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l

                case(3)
                    ! Face 3: eta = -1
                    do l=1,nglx
                        do m=1,nglz
                            ww=wglx(l)*wglz(m)
                            i = l
                            j = 1
                            k = m
                            nx = + y_ksi(i,j,k)*z_zeta(i,j,k) - z_ksi(i,j,k)*y_zeta(i,j,k)
                            ny = - x_ksi(i,j,k)*z_zeta(i,j,k) + z_ksi(i,j,k)*x_zeta(i,j,k)
                            nz = + x_ksi(i,j,k)*y_zeta(i,j,k) - y_ksi(i,j,k)*x_zeta(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            jac_face(l,m,iface) = ww*nlen
                            nv(1,l,m,iface) = nx/nlen
                            nv(2,l,m,iface) = ny/nlen
                            nv(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l

                case(4)
                    ! Face 4: eta = +1
                    do l=1,nglx
                        do m=1,nglz
                            ww=wglx(l)*wglz(m)
                            i = l
                            j = ngly
                            k = m
                            nx = - y_ksi(i,j,k)*z_zeta(i,j,k) + z_ksi(i,j,k)*y_zeta(i,j,k)
                            ny = + x_ksi(i,j,k)*z_zeta(i,j,k) - z_ksi(i,j,k)*x_zeta(i,j,k)
                            nz = - x_ksi(i,j,k)*y_zeta(i,j,k) + y_ksi(i,j,k)*x_zeta(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            jac_face(l,m,iface) = ww*nlen
                            nv(1,l,m,iface) = nx/nlen
                            nv(2,l,m,iface) = ny/nlen
                            nv(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l

                case(5)
                    ! Face 5: ksi = -1
                    do l=1,ngly
                        do m=1,nglz
                            ww=wgly(l)*wglz(m)
                            i = 1
                            j = l
                            k = m
                            nx = - y_eta(i,j,k)*z_zeta(i,j,k) + z_eta(i,j,k)*y_zeta(i,j,k)
                            ny = + x_eta(i,j,k)*z_zeta(i,j,k) - z_eta(i,j,k)*x_zeta(i,j,k)
                            nz = - x_eta(i,j,k)*y_zeta(i,j,k) + y_eta(i,j,k)*x_zeta(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            jac_face(l,m,iface) = ww*nlen
                            nv(1,l,m,iface) = nx/nlen
                            nv(2,l,m,iface) = ny/nlen
                            nv(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l

                case(6)
                    ! Face 6: ksi = +1
                    do l=1,ngly
                        do m=1,nglz
                            ww=wgly(l)*wglz(m)
                            i = nglx
                            j = l
                            k = m
                            nx = + y_eta(i,j,k)*z_zeta(i,j,k) - z_eta(i,j,k)*y_zeta(i,j,k)
                            ny = - x_eta(i,j,k)*z_zeta(i,j,k) + z_eta(i,j,k)*x_zeta(i,j,k)
                            nz = + x_eta(i,j,k)*y_zeta(i,j,k) - y_eta(i,j,k)*x_zeta(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            jac_face(l,m,iface) = ww*nlen
                            nv(1,l,m,iface) = nx/nlen
                            nv(2,l,m,iface) = ny/nlen
                            nv(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l
            end select
        end if
    end do !iface

end subroutine create_normals
!----------------------------------------------------------------------!
!>@brief Create (global) data structure imapl and imapr, which denote the node
!> points which need to be sent to the left element and the right element
!> in the construction of DG-fluxes
!>@author James F. Kelly
!>@date 11 October 2010
!>@date 28 May 2015, Daniel S. Abdi
!----------------------------------------------------------------------!
subroutine create_imaplr(imapl,imapr,nv,jac_face,face,nface)

    use mod_basis, only: ngl, ngl2, nglx, ngly, nglz, FACE_LEN
  
    use mod_grid, only: intma, coord, mod_grid_get_face_ngl, face_type

    use mod_input, only : lgpu, is_non_conforming_flg
  
    implicit none
  
    integer nface, ngl_i, ngl_j, ngl_ij, plane_ij
    integer face(FACE_LEN,nface), imapl(3,ngl,ngl,nface), imapr(3,ngl,ngl,nface)
    integer iface, ilocl, ilocr, i, j, k, ii, l, m, iel, ier, il, ir,ip
    integer jmapl(3,ngl2), jmapr(3,ngl2),  orderl(ngl2), orderr(ngl2)
    real nodel(ngl2), noder(ngl2)
    real nv(3,ngl,ngl,nface),jac_face(ngl,ngl,nface)
    real nvc(3,ngl,ngl), jac_facec(ngl,ngl)
    real x, y, z

    integer ik,jk,kk

    do iface = 1,nface
        ilocl=face(5,iface)
        ilocr=face(6,iface)
        iel= face(7,iface)
        ier= face(8,iface)

        call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)
        if(face_type(iface)==21) call mod_grid_get_face_ngl(ilocr, ngl_i, ngl_j,ngl_ij)

        ii = 0
        do l = 1,ngl_i
            do m = 1,ngl_j
                ii = ii + 1

                !store Left-Side local GridPoint Pointers
                select case (ilocl)
                    case(1)
                        ! Face 1: zeta = -1
                        i = l
                        j = m
                        k = 1
                    case(2)
                        ! Face 2: zeta = +1
                        i = l
                        j = m
                        k = nglz
                    case(3)
                        ! Face 3: eta = -1
                        i = l
                        j = 1
                        k = m
                    case(4)
                        ! Face 4: eta = +1
                        i = l
                        j = ngly
                        k = m
                    case(5)
                        ! Face 5: ksi = -1
                        i = 1
                        j = l
                        k = m
                    case(6)
                        ! Face 6: ksi = +1
                        i = nglx
                        j = l
                        k = m
                end select
                jmapl(1,ii)=i
                jmapl(2,ii)=j
                jmapl(3,ii)=k

                if(iel>0) then
                    !sort based on coordinate NOT local node id
                    ip = intma(i,j,k,iel);
                    x = coord(1,ip);
                    y = coord(2,ip);
                    z = coord(3,ip);
                    nodel(ii)=sqrt(x*x+y*y+z*z)+x-2*y+3*z
                else
                    nodel(ii)=ii
                end if

                !Now Store Right-Side Local GridPoint Pointers
                select case (ilocr)

                    case(1)
                        ! Face 1: zeta = -1
                        i = l
                        j = m
                        k = 1
                    case(2)
                        ! Face 2: zeta = +1
                        i = l
                        j = m
                        k = nglz
                    case(3)
                        ! Face 3: eta = -1
                        i = l
                        j = 1
                        k = m
                    case(4)
                        ! Face 4: eta = +1
                        i = l
                        j = ngly
                        k = m
                    case(5)
                        ! Face 5: ksi = -1
                        i = 1
                        j = l
                        k = m
                    case(6)
                        ! Face 6: ksi = +1
                        i = nglx
                        j = l
                        k = m
                end select
                jmapr(1,ii)=i
                jmapr(2,ii)=j
                jmapr(3,ii)=k
                if ((ier > 0).and.(face_type(iface).ne.21)) then
                    !sort based on coordinate NOT local node id
                    ip = intma(i,j,k,ier);
                    x = coord(1,ip);
                    y = coord(2,ip);
                    z = coord(3,ip);
                    noder(ii)=sqrt(x*x+y*y+z*z)+x-2*y+3*z
                else
                    noder(ii)=nodel(ii)
                end if
            end do   !m
        end do !l
     
        ngl_ij=ngl_i * ngl_j
        do l=1,ngl_ij
            orderl(l)=l
            orderr(l)=l
        end do

        !store imapl,imapr
        ii = 0
        do i=1,ngl_i
            do j=1,ngl_j
                ii = ii + 1
                il=orderl(ii)
                ir=orderr(ii)
                imapl(:,i,j,iface)=jmapl(:,il)
                imapr(:,i,j,iface)=jmapr(:,ir)
         
                !sort normals
                l = (il - 1)/ngl_j + 1
                m = mod(il - 1,ngl_j) + 1
                nvc(:,i,j)=nv(:,l,m,iface)
                jac_facec(i,j)=jac_face(l,m,iface)
            end do
        end do !l

        !Store sorted normals
        do i=1,ngl_i
            do j=1,ngl_j
                nv(:,i,j,iface)=nvc(:,i,j)
                jac_face(i,j,iface)=jac_facec(i,j)
            end do
        end do !l


    end do !iface

!stop

end subroutine create_imaplr
