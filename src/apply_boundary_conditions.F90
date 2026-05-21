!-------------------------------------------------
!>@brief Create NORMALS and PMATRIX which are used for NFBC
!-------------------------------------------------
subroutine create_nfbc_vector(G, b, mf, normals, ipoin_bound, npoin_bound, ldirichlet)

  use mod_grid,  only: grid, mod_grid_get_face_ngl
  use mod_basis, only: basis
  use mod_face,  only: face_CS

  implicit none

  type(grid),    intent(in) :: G
  type(basis),   intent(in) :: b
  type(face_CS), intent(in) :: mf

  !global arrays
  real,    intent(out) :: normals(3,G%npoin)
  integer, intent(out) :: ipoin_bound(G%npoin)
  integer, intent(out) :: npoin_bound
  logical, intent(out) :: ldirichlet

  !local array
  real, dimension(:), allocatable :: ipoin
  integer :: AllocateStatus
  real :: rnx, rny, rnz, rnr
  real :: x, y, z, s
  integer :: iface, iel, ier, ilocl, ilocr, ip, jp, ipcg, i, j, il, jl, kl
  integer :: nbb, xflag, yflag, zflag, ngl_i, ngl_j, plane_ij, ip_g, jp_g

  allocate(ipoin(G%npoin), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop "** Not Enough Memory - CREATE_PMATRIX **"

  !Initialize
  ipoin=0
  ipoin_bound=0
  normals=0
  ldirichlet=.false.

  !Modify Normal Vectors
  do iface=1,G%nface
     iel=G%face(7,iface)
     ier=G%face(8,iface)
     ilocl=G%face(5,iface)

     if (ier == -4 .or. ier == -8) then !n*U=0

        call mod_grid_get_face_ngl(b, ilocl, ngl_i, ngl_j, plane_ij)

        do j=1,ngl_j
           do i=1,ngl_i

              il=mf%imapl(1,i,j,iface)
              jl=mf%imapl(2,i,j,iface)
              kl=mf%imapl(3,i,j,iface)
              ip=G%intma(il,jl,kl,iel)

              !Store Normals
              rnx=mf%normal_vector(1,i,j,iface)
              rny=mf%normal_vector(2,i,j,iface)
              rnz=mf%normal_vector(3,i,j,iface)
              ipoin(ip)=ipoin(ip) + 1

              !Store Normals
              normals(1,ip)=normals(1,ip) + rnx
              normals(2,ip)=normals(2,ip) + rny
              normals(3,ip)=normals(3,ip) + rnz
           end do !i
        end do !j

     end if !ier

  end do !iface

  nbb=0
  do i=1,G%npoin
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

  deallocate(ipoin)

end subroutine create_nfbc_vector
