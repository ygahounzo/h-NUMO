!----------------------------------------------------------------------!
!>This module builds the Faces and Normals at quadrature points
!>@ author by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: March 27, 2023
!----------------------------------------------------------------------!
subroutine create_normals_quad(nv_q,jac_faceq,face,nface)

    use mod_basis, only: nq, nqx, nqy, nqz, wnqx, wnqy, wnqz, FACE_LEN, nglx,ngly,nglz

    use mod_grid, only: intma_dg_quad, coord, intma

    use mod_gradient, only: compute_local_gradient_quad_v3

    implicit none

    !global
    integer, intent(in) :: nface
    real, dimension(3,nq,nq,nface), intent(out) :: nv_q
    real, dimension(nq,nq,nface), intent(out) :: jac_faceq
    integer, intent(in) :: face(FACE_LEN,nface)

    !local
    integer :: iface, i, j, k, ip, l, m
    integer :: ilocl, ilocr, iel, ier, ndim
    real :: ww, nx, ny, nz, nlen

    !local arrays
    real, dimension(nglx,ngly,nglz) :: x, y, z
    real, dimension(nqx,nqy,nqz) :: x_ksiq, x_etaq, x_zetaq
    real, dimension(nqx,nqy,nqz) :: y_ksiq, y_etaq, y_zetaq
    real, dimension(nqx,nqy,nqz) :: z_ksiq, z_etaq, z_zetaq
  
    !Initialize MXM dimension
    ndim=1

    !initialize the global matrix
    nv_q=0.0
    jac_faceq=0.0

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

                        !print* ,'x = ',x(i,j,k)
                    end do !i
                end do !j
            end do !k

            !Construct Mapping Derivatives

            call compute_local_gradient_quad_v3(x_ksiq,x_etaq,x_zetaq,x,nglx,ngly,nglz,nqx,nqy,nqz)
            call compute_local_gradient_quad_v3(y_ksiq,y_etaq,y_zetaq,y,nglx,ngly,nglz,nqx,nqy,nqz)
            call compute_local_gradient_quad_v3(z_ksiq,z_etaq,z_zetaq,z,nglx,ngly,nglz,nqx,nqy,nqz)

            !set jacobian matrix J for 2D

            do k=1,nqz
                do j=1,nqy
                    do i=1,nqx
                        if(nqx == 1) then
                            x_ksiq(i,j,k) = 1.0;  y_ksiq(i,j,k) = 0.0;  z_ksiq(i,j,k) = 0.0;
                            x_etaq(i,j,k) = 0.0;  x_zetaq(i,j,k) = 0.0;
                        endif
                        if(nqy == 1) then
                            x_etaq(i,j,k) = 0.0;  y_etaq(i,j,k) = 1.0;  z_etaq(i,j,k) = 0.0;
                            y_ksiq(i,j,k) = 0.0; y_zetaq(i,j,k) = 0.0;
                        endif
                        if(nqz == 1) then
                            x_zetaq(i,j,k) = 0.0; y_zetaq(i,j,k) = 0.0; z_zetaq(i,j,k) = 1.0;
                            z_ksiq(i,j,k) = 0.0; z_etaq(i,j,k) = 0.0;
                        endif
                    enddo
                enddo
            enddo
              
            !Compute Normals
            select case (ilocl)
           
                case(1)
                    ! Face 1: zeta = -1
                    do l=1,nqx
                        do m = 1,nqy
                            ww=wnqx(l)*wnqy(m)
                            i = l
                            j = m
                            k = 1
                            nx = + y_etaq(i,j,k)*z_ksiq(i,j,k) - z_etaq(i,j,k)*y_ksiq(i,j,k)
                            ny = - x_etaq(i,j,k)*z_ksiq(i,j,k) + z_etaq(i,j,k)*x_ksiq(i,j,k)
                            nz = + x_etaq(i,j,k)*y_ksiq(i,j,k) - y_etaq(i,j,k)*x_ksiq(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            !print*, 'nlen = ', nlen
                            jac_faceq(l,m,iface) = ww*nlen
                            nv_q(1,l,m,iface) = nx/nlen
                            nv_q(2,l,m,iface) = ny/nlen
                            nv_q(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l
           
                case(2)
                    ! Face 2: zeta = +1
                    do l=1,nqx
                        do m=1,nqy
                            ww=wnqx(l)*wnqy(m)
                            i = l
                            j = m
                            k = nqz
                            nx = - y_etaq(i,j,k)*z_ksiq(i,j,k) + z_etaq(i,j,k)*y_ksiq(i,j,k)
                            ny = + x_etaq(i,j,k)*z_ksiq(i,j,k) - z_etaq(i,j,k)*x_ksiq(i,j,k)
                            nz = - x_etaq(i,j,k)*y_ksiq(i,j,k) + y_etaq(i,j,k)*x_ksiq(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            jac_faceq(l,m,iface) = ww*nlen
                            nv_q(1,l,m,iface) = nx/nlen
                            nv_q(2,l,m,iface) = ny/nlen
                            nv_q(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l

                case(3)
                    ! Face 3: eta = -1
                    do l=1,nqx
                        do m=1,nqz
                            ww=wnqx(l)*wnqz(m)
                            i = l
                            j = 1
                            k = m
                            nx = + y_ksiq(i,j,k)*z_zetaq(i,j,k) - z_ksiq(i,j,k)*y_zetaq(i,j,k)
                            ny = - x_ksiq(i,j,k)*z_zetaq(i,j,k) + z_ksiq(i,j,k)*x_zetaq(i,j,k)
                            nz = + x_ksiq(i,j,k)*y_zetaq(i,j,k) - y_ksiq(i,j,k)*x_zetaq(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            !print*, 'nlen = ', nlen
                            jac_faceq(l,m,iface) = ww*nlen
                            nv_q(1,l,m,iface) = nx/nlen
                            nv_q(2,l,m,iface) = ny/nlen
                            nv_q(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l

                case(4)
                    ! Face 4: eta = +1
                    do l=1,nqx
                        do m=1,nqz
                            ww=wnqx(l)*wnqz(m)
                            i = l
                            j = nqy
                            k = m
                            nx = - y_ksiq(i,j,k)*z_zetaq(i,j,k) + z_ksiq(i,j,k)*y_zetaq(i,j,k)
                            ny = + x_ksiq(i,j,k)*z_zetaq(i,j,k) - z_ksiq(i,j,k)*x_zetaq(i,j,k)
                            nz = - x_ksiq(i,j,k)*y_zetaq(i,j,k) + y_ksiq(i,j,k)*x_zetaq(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            jac_faceq(l,m,iface) = ww*nlen
                            nv_q(1,l,m,iface) = nx/nlen
                            nv_q(2,l,m,iface) = ny/nlen
                            nv_q(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l

                case(5)
                    ! Face 5: ksi = -1
                    do l=1,nqy
                        do m=1,nqz
                            ww=wnqy(l)*wnqz(m)
                            i = 1
                            j = l
                            k = m
                            nx = - y_etaq(i,j,k)*z_zetaq(i,j,k) + z_etaq(i,j,k)*y_zetaq(i,j,k)
                            ny = + x_etaq(i,j,k)*z_zetaq(i,j,k) - z_etaq(i,j,k)*x_zetaq(i,j,k)
                            nz = - x_etaq(i,j,k)*y_zetaq(i,j,k) + y_etaq(i,j,k)*x_zetaq(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            jac_faceq(l,m,iface) = ww*nlen
                            nv_q(1,l,m,iface) = nx/nlen
                            nv_q(2,l,m,iface) = ny/nlen
                            nv_q(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l

                case(6)
                    ! Face 6: ksi = +1
                    do l=1,nqy
                        do m=1,nqz
                            ww=wnqy(l)*wnqz(m)
                            i = nqx
                            j = l
                            k = m
                            nx = + y_etaq(i,j,k)*z_zetaq(i,j,k) - z_etaq(i,j,k)*y_zetaq(i,j,k)
                            ny = - x_etaq(i,j,k)*z_zetaq(i,j,k) + z_etaq(i,j,k)*x_zetaq(i,j,k)
                            nz = + x_etaq(i,j,k)*y_zetaq(i,j,k) - y_etaq(i,j,k)*x_zetaq(i,j,k)
                            nlen = sqrt(nx*nx + ny*ny + nz*nz)
                            jac_faceq(l,m,iface) = ww*nlen
                            nv_q(1,l,m,iface) = nx/nlen
                            nv_q(2,l,m,iface) = ny/nlen
                            nv_q(3,l,m,iface) = nz/nlen
                        end do !m
                    end do !l
            end select
        end if
    end do !iface

end subroutine create_normals_quad


!----------------------------------------------------------------------!
!>@brief Create (global) data structure imapl_q and imapr_q, which denote the node
!> points which need to be sent to the left element and the right element
!>@ author by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: March 27, 2023
!----------------------------------------------------------------------!
subroutine create_imaplr_quad(imapl_q,imapr_q,nv_q,jac_faceq,face,nface)

    use mod_basis, only: nq, nq2, nqx, nqy, nqz, FACE_LEN
  
    use mod_grid, only: intma_dg_quad, mod_grid_get_face_nq, face_type

    use mod_input, only : lgpu, is_non_conforming_flg
  
    implicit none
  
    integer nface, nq_i, nq_j, nq_ij, plane_ij
    integer face(FACE_LEN,nface), imapl_q(3,nq,nq,nface), imapr_q(3,nq,nq,nface)
    integer iface, ilocl, ilocr, i, j, k, ii, l, m, iel, ier, il, ir,ip
    integer jmapl_q(3,nq2), jmapr_q(3,nq2),  orderl(nq2), orderr(nq2)
    real nv_q(3,nq,nq,nface),jac_faceq(nq,nq,nface)
    real nvc_q(3,nq,nq), jac_facecq(nq,nq)
    real x, y, z

    integer ik,jk,kk

    do iface = 1,nface
        ilocl=face(5,iface)
        ilocr=face(6,iface)
        iel= face(7,iface)
        ier= face(8,iface)

        call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)

        ii = 0
        do l = 1,nq_i
            do m = 1,nq_j
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
                        k = nqz
                    case(3)
                        ! Face 3: eta = -1
                        i = l
                        j = 1
                        k = m
                    case(4)
                        ! Face 4: eta = +1
                        i = l
                        j = nqy
                        k = m
                    case(5)
                        ! Face 5: ksi = -1
                        i = 1
                        j = l
                        k = m
                    case(6)
                        ! Face 6: ksi = +1
                        i = nqx
                        j = l
                        k = m
                end select
                jmapl_q(1,ii)=i
                jmapl_q(2,ii)=j
                jmapl_q(3,ii)=k

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
                        k = nqz
                    case(3)
                        ! Face 3: eta = -1
                        i = l
                        j = 1
                        k = m
                    case(4)
                        ! Face 4: eta = +1
                        i = l
                        j = nqy
                        k = m
                    case(5)
                        ! Face 5: ksi = -1
                        i = 1
                        j = l
                        k = m
                    case(6)
                        ! Face 6: ksi = +1
                        i = nqx
                        j = l
                        k = m
                end select
                jmapr_q(1,ii)=i
                jmapr_q(2,ii)=j
                jmapr_q(3,ii)=k
            end do   !m
        end do !l
     
        nq_ij=nq_i * nq_j
        do l=1,nq_ij
            orderl(l)=l
            orderr(l)=l
        end do

        !store imapl_q,imapr_q
        ii = 0
        do i=1,nq_i
            do j=1,nq_j
                ii = ii + 1
                il=orderl(ii)
                ir=orderr(ii)
                imapl_q(:,i,j,iface)=jmapl_q(:,il)
                imapr_q(:,i,j,iface)=jmapr_q(:,ir)
         
                !sort normals
                l = (il - 1)/nq_j + 1
                m = mod(il - 1,nq_j) + 1
                nvc_q(:,i,j)=nv_q(:,l,m,iface)
                jac_facecq(i,j)=jac_faceq(l,m,iface)
            end do
        end do !l

        !Store sorted normals
        do i=1,nq_i
            do j=1,nq_j
                nv_q(:,i,j,iface)=nvc_q(:,i,j)
                jac_faceq(i,j,iface)=jac_facecq(i,j)
            end do
        end do !l


    end do !iface

!stop

end subroutine create_imaplr_quad
