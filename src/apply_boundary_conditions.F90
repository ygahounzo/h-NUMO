!-------------------------------------------------
!>@brief Create NORMALS and PMATRIX which are used for NFBC
!-------------------------------------------------
subroutine create_nfbc_vector(normals,ipoin_bound,npoin_bound,ldirichlet)

  use mod_basis, only: nglx, ngly, nglz, is_2d

  use mod_constants, only: tol

  use mod_face, only: normal_vector, imapl

  use mod_global_grid, only: xmin, xmax, ymin, ymax, zmin, zmax, iboundary

  use mod_grid, only: intma, intma_dg_to_cg, npoin, coord, sigma, face, nface, &
       mod_grid_get_face_ngl

  use mod_initial, only: kvector

  use mod_input, only: space_method

  use mod_parallel, only: num_send_recv_total

  implicit none

  !global arrays
  real, intent(out)    :: normals(3,npoin)
  integer, intent(out) :: ipoin_bound(npoin)
  integer, intent(out) :: npoin_bound
  logical, intent(out) :: ldirichlet

  !local array
  real, dimension(:,:), allocatable :: recv_data_matrix
  real, dimension(:),   allocatable :: recv_data_vector, ipoin
  integer :: AllocateStatus
  real :: rnx, rny, rnz, rnr
  real :: x, y, z, s
  integer :: iface, iel, ier, ilocl, ilocr, ip, jp, ipcg, i, j, il, jl, kl
  integer :: nbb, xflag, yflag, zflag, ngl_i, ngl_j, plane_ij, ip_g, jp_g

  allocate (ipoin(npoin), recv_data_vector(num_send_recv_total),&
       recv_data_matrix(3,num_send_recv_total), stat=AllocateStatus )
  if (AllocateStatus /= 0) stop "** Not Enough Memory - CREATE_PMATRIX **"

  !Initialize
  ipoin=0
  ipoin_bound=0
  normals=0
  ldirichlet=.false.

  !Modify Normal Vectors
  do iface=1,nface
     iel=face(7,iface)
     ier=face(8,iface)
     ilocl=face(5,iface)

     if (ier == -4 .or. ier == -8) then !n*U=0

        call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

        do j=1,ngl_j
           do i=1,ngl_i

              il=imapl(1,i,j,iface)
              jl=imapl(2,i,j,iface)
              kl=imapl(3,i,j,iface)
              ip=intma(il,jl,kl,iel)

              !Store Normals
              rnx=normal_vector(1,i,j,iface)
              rny=normal_vector(2,i,j,iface)
              rnz=normal_vector(3,i,j,iface)
              ipoin(ip)=ipoin(ip) + 1

              !Store Normals
              normals(1,ip)=normals(1,ip) + rnx
              normals(2,ip)=normals(2,ip) + rny
              normals(3,ip)=normals(3,ip) + rnz
           end do !i
        end do !j

     end if !ier

  end do !iface

  !Perform Global Assembly for the Normal Vectors and
  if (space_method == 'cgc') then
     call create_global_rhs(normals,recv_data_matrix,3,0)
     call create_global_rhs(ipoin,recv_data_matrix,1,0)
  end if

  nbb=0
  do i=1,npoin
     if (ipoin(i) /= 0) then
        nbb=nbb + 1
        ipoin_bound(nbb)=i
        rnx=normals(1,i)
        rny=normals(2,i)
        rnz=normals(3,i)
        rnr=sqrt( rnx*rnx + rny*rny + rnz*rnz)
        if(abs(rnr) <= 1e-7) rnr=1.0
        normals(1,i)=rnx/rnr
        normals(2,i)=rny/rnr
        normals(3,i)=rnz/rnr
     end if
  end do
  npoin_bound=nbb

  deallocate (ipoin, recv_data_vector, recv_data_matrix)

end subroutine create_nfbc_vector