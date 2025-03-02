!------------------------------------------------------------------------------!
!>@brief This module contains the call to the library p4est
!
!>@author Written by Simone Marras on Aug 26, 2014
!> Naval Postgraduate School, Monterey, CA
!>
!>@author Modified by Andreas Mueller on October 15, 2014
!> added p4est_initial
!>
!>@author Modified by Michal Kopera on December 15, 2014
!> added p6est capability
!>
!>@author Modified by Simone Marras on Feb 26, 2015
!> removed  use mod_io_grid, only: mod_write_grid because 'call mod_write_grid'
!> was unused
!------------------------------------------------------------------------------!

#define approx_equal(x, y, r, a) \
  abs((x)-(y)) <= (a) + (r) * max(abs((x)), abs((y)))

#define SAFE_DEALLOCATE(x) if(allocated(x)) deallocate(x)

module mod_p4est

  use mod_input, only: nelx, nely, nelz, nproc_z, nopx, nopy, nopz, xdims,     &
    ydims, ztop, zbottom, x_boundary, y_boundary, z_boundary,     &
    x_periodic, y_periodic, z_periodic, nel_root_h,            &
    refinement_levels_h, &
    space_method, read_external_grid_flg, is_non_conforming_flg, xlim_min,     &
    xlim_max, ylim_min, ylim_max, zlim_min, zlim_max, space_method,            &
    luse_hybrid_cpu_gpu, platformWeight, platformWeight2, cpus_per_node,       &
    gpus_per_node, amr_indicator_variables, amr_smoothness_limits,  &
    amr_smoothness_qL2_limit, lrestoring_sponge

  use mod_basis, only: nopz, ngl, nglx, ngly, nglz, xgl, xglx, xgly, xglz,     &
    npts, is_2d, FACE_CHILDREN, FACE_LEN, P4EST_FACES, P8EST_EDGES

  use mod_bc, only: vc_el_type

  use mod_global_grid, only: nelem_g, npoin_g_cg, npoin_g, ncol_g, xmin, xmax, &
    ymin, ymax, zmin, zmax, iboundary, xperiodic, yperiodic, zperiodic

  use mod_face, only: mod_face_create_nc_list

  use mod_grid, only: intma, intma_table, nelemx, nelemy, nelemz,              &
    npoin_cg, npoin, nelem, coord_cg, coord, nbsido, bsido, sigma,             &
    ncol, ncol_cg, nz_cg, node_column_table, face, nface, face_type, nboun,    &
    index2d, mod_grid_init_unified, mod_grid_init_coord,                       &
    mod_grid_init_coord_dg_to_cg, is_non_conforming,                           &
    NC_face, NC_edge, EToNC, nNC, D2C_mask

  use mod_metrics, only: mod_metrics_create_metrics

  use mod_mpi_utilities, only : irank, irank0, MPI_PRECISION

  use mod_parallel, only: nproc, num_nbh, num_send_recv, num_send_recv_total,  &
    nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half,  nbh_proc,         &
    ipoin_proc, npoin_l, npoin_l_max, ncol_l, ncol_l_max, nelem_l,             &
    nelem_l_max, nboun_max

  use mod_types, only : r8

  use iso_c_binding, only: C_INT64_T

  use mpi

  public :: &
       mod_p4est_create, &
       mod_p4est_bcs, &
       mod_p4est_dimensions, &
       mod_p4est_adapt, &
       mod_p4est_out, &
       mod_p4est_finalize, &
       mod_p4est_dump_mesh, &
       scatter_element_2d, scatter_element_2d_subface, &
       gather_element_2d, gather_element_2d_subface, &
       plist, lev_list, partition, npartition, &
       cs_face, interp_x, interp_y, interp_z

  private
  !-----------------------------------------------------------------------
  !module variables and parameters
  type :: p4esttonuma
     private
     integer, pointer :: p
  end type p4esttonuma

  real,    dimension(:, :, :), allocatable :: PsgXY, PsgXZ, PsgYZ
  integer, dimension(:),    allocatable :: lev_list, plist, partition
  integer :: npartition
  integer, dimension(:),    allocatable :: refine_coarsen_elements


  ! 1-D interpolation/projection matrices:
  real, dimension(:, :, :), allocatable :: interp_x, project_x
  real, dimension(:, :, :), allocatable :: interp_y, project_y
  real, dimension(:, :, :), allocatable :: interp_z, project_z

  real, dimension(:, :), allocatable :: Vx_inv, Vy_inv, Vz_inv
  real, dimension(:), allocatable :: bx, by, bz

  integer(C_INT64_T), dimension(:), allocatable :: qid_src, qid_dst

  integer, dimension(:), allocatable :: cs_face, csp_face
  logical :: init_refine = .false.


  !-----------------------------------------------------------------------

contains

  subroutine mod_p4est_finalize()
    implicit none
    call p4est_stop()
    SAFE_DEALLOCATE(PsgXY)
    SAFE_DEALLOCATE(PsgXZ)
    SAFE_DEALLOCATE(PsgYZ)
    SAFE_DEALLOCATE(lev_list)
    SAFE_DEALLOCATE(plist)
    SAFE_DEALLOCATE(partition)
    SAFE_DEALLOCATE(refine_coarsen_elements)
    SAFE_DEALLOCATE(interp_x)
    SAFE_DEALLOCATE(interp_y)
    SAFE_DEALLOCATE(interp_z)
    SAFE_DEALLOCATE(project_x)
    SAFE_DEALLOCATE(project_y)
    SAFE_DEALLOCATE(project_z)
    SAFE_DEALLOCATE(qid_src)
    SAFE_DEALLOCATE(qid_dst)
    SAFE_DEALLOCATE(Vx_inv)
    SAFE_DEALLOCATE(Vy_inv)
    SAFE_DEALLOCATE(Vz_inv)
    SAFE_DEALLOCATE(bx)
    SAFE_DEALLOCATE(by)
    SAFE_DEALLOCATE(bz)
  end subroutine

!-----------------------------------------------------------------------
!@>brief Create p4est grid + static AMR
!-----------------------------------------------------------------------
  ! {{{

  subroutine mod_p4est_create()
    use mod_initial, only: mod_initial_create, q_init
    use mod_input, only: max_mesh_lvl, p4est_log_level
    use iso_c_binding, only: C_INT
    implicit none

    integer:: ierr, temp, k
    integer(C_INT) :: p4est_log_level_c

    p4est_log_level_c = p4est_log_level

    call mpi_comm_size(mpi_comm_world, nproc, ierr)

    allocate(npoin_l(nproc), ncol_l(nproc), nelem_l(nproc))
    allocate(partition(1))
    allocate(refine_coarsen_elements(1))

    !create projection matrices
    if(is_non_conforming_flg > 0) then
        call create_2d_projection_matrices_numa2d(PsgXY,1)
        call create_2d_projection_matrices_numa2d(PsgXZ,2)
        call create_2d_projection_matrices_numa2d(PsgYZ,3)
    endif

    if(.not.allocated(interp_x)) then
      call mod_p4est_build_projection_1d(interp_x, project_x, xglx)
      call mod_p4est_build_projection_1d(interp_y, project_y, xgly)
      call mod_p4est_build_projection_1d(interp_z, project_z, xglz)
    endif

    !start p4est
    call p4est_start(p4est_log_level_c)

    !create initial grid
    temp = is_non_conforming_flg
    is_non_conforming_flg = 0
    call mod_p4est_create_grid(.true.)
    is_non_conforming_flg = temp

    !refine initial grid
    if(is_non_conforming_flg == 1 .or. is_non_conforming_flg == 3) then
      call mod_p4est_refine()
    elseif(is_non_conforming_flg == 4) then
      if(irank == irank0) print*, '----------- Begin iterative initial refinement -------------'
      init_refine = .true.
      do k = 1, max_mesh_lvl
        call mod_initial_create()
        call mod_p4est_adapt(q_init)
      enddo
      init_refine = .false.
      if(irank == irank0) print*, '----------- End iterative initial refinement -------------'
    endif

  end subroutine

!-----------------------------------------------------------------------
!@>brief Refine grid
!-----------------------------------------------------------------------
  subroutine mod_p4est_refine()
    implicit none

    if(irank == irank0) print*, '----------- Begin refinement and coarsening -------------'
    call set_refine_coarsen_flags()
    call mod_p4est_create_grid(.true.)
    call mod_face_create_nc_list()

    if(irank == irank0) print*, '----------- End refinement and coarsening   -------------'

  end subroutine mod_p4est_refine

!-----------------------------------------------------------------------
!@>brief Generate p4est grid
!-----------------------------------------------------------------------
  subroutine mod_p4est_create_grid(init_p4est)

    use mpi

    implicit none

    logical, intent(in) :: init_p4est

    integer :: nop
    integer :: ierr, i, j, k, ip, e
    character(len=72) :: fnp, fnp1, fname
    integer :: iloop, AllocateStatus
    real :: x, y, z
    integer :: lnx, lny, lnz, orient = 0
    integer :: log_root_v
    real :: height_correction

    type (p4esttonuma) :: p2n

    integer :: is_cube=1
    integer :: is_dg=0
    integer :: is_cgc=0

    integer irank

    nop = ngl - 1

    call mpi_comm_rank(mpi_comm_world, irank, ierr)

    !------------------
    !  Boundary flags
    !------------------
    iboundary(1)=z_boundary(1) !bottom (z=-1)
    iboundary(2)=z_boundary(2) !top    (z=+1)
    iboundary(3)=y_boundary(1) !left   (y=-1)
    iboundary(4)=y_boundary(2) !right  (y=+1)
    iboundary(5)=x_boundary(1) !front  (x=-1)
    iboundary(6)=x_boundary(2) !back   (x=+1)

    !-----------------------------------------
    ! Initialize p4est
    !------------------------------------------
    if(space_method == 'dg') is_dg = 1
    if(space_method == 'cgc') is_cgc = 1
    
    !--- p4est initialization

    if(nglx == 1) then
      lnx = nely; lny = nelz
      orient = 0
    else if(ngly == 1) then
      lnx = nelx; lny = nelz
      orient = 1
    else
      orient = 2
      lnx = nelx; lny = nely
    endif
    lnz = 1

    if(init_p4est) &
      call p4esttonuma_init(is_cube, lnx, lny, lnz, &
      refinement_levels_h, read_external_grid_flg, &
      is_non_conforming_flg, refine_coarsen_elements, &
      x_periodic, y_periodic, z_periodic, FACE_LEN, orient, lrestoring_sponge)

    ! TODO: fix to use xglx, xgly, xglz
    call p4esttonuma_fill_data(nop, xgl, is_dg, iboundary, read_external_grid_flg, p2n, lrestoring_sponge)

    call p4esttonuma_get_mesh_scalars(p2n, npoin_cg, nelem, num_nbh, &
      num_send_recv_total, nbsido, nface, nboun, nNC)

    if(is_non_conforming_flg == 0) then
      nelemx = nelx *2**refinement_levels_h
      nelemy = nely *2**refinement_levels_h
      nelemz = nelz *2**refinement_levels_h
    end if

    height_correction = 1

    call mod_grid_init_unified()

    !-----------------------------------------
    ! Allocate parallel data structures
    !------------------------------------------
    if(allocated(nbh_proc)) then
    deallocate(nbh_proc, &
        ipoin_proc, &
        num_send_recv, &
        nbh_send_recv, &
        nbh_send_recv_multi, &
        nbh_send_recv_half, &
        lev_list, plist)
    endif
    allocate(nbh_proc(num_nbh), &
         ipoin_proc(npoin), &
         num_send_recv(num_nbh), &
         nbh_send_recv(num_send_recv_total), &
         nbh_send_recv_multi(num_send_recv_total), &
         nbh_send_recv_half(num_send_recv_total*FACE_CHILDREN), &
         lev_list(nelem), plist(nelem), &
         stat=AllocateStatus )
    if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_p4est **"

    if (allocated(cs_face)) deallocate(cs_face)
    allocate(cs_face(nelem))
    if (allocated(csp_face)) deallocate(csp_face)
    allocate(csp_face(npoin))

    ! Initialize lev_list to 0
    lev_list(:) = 0
    
    !-----------------------------------------
    ! Get grid data from p4est
    !------------------------------------------


    if (irank == irank0) print*, '------------------Entering P4est_Mesh_Arrays-------------------------'
      ! FIXME: Kill the NC_face and NC_edge arrays?
      call p4esttonuma_get_mesh_arrays(p2n, coord, intma_table,           &
          NC_face, NC_edge, EToNC, face, face_type,                      &
          num_send_recv, nbh_proc, nbh_send_recv, nbh_send_recv_multi,   &
          nbh_send_recv_half, bsido, lev_list, plist, is_cgc, vc_el_type, &
          lrestoring_sponge, D2C_mask)
      call p4esttonuma_free(p2n,lrestoring_sponge)

    if (irank == irank0) print*, '------------------Leaving P4est_Mesh_Arrays-------------------------'

    !adjust cube dimensions
    if(is_2d.and.read_external_grid_flg==0) then
        do i=1, npoin
!                if(read_external_grid_flg==0) then !don't adjust if mesh read from file
              if(nglx == 1) then
                  x=0; y=coord(1, i)/nely; z=coord(2, i)/nelz
              else if(ngly == 1) then
                  x=coord(1, i)/nelx; y=0; z=coord(2, i)/nelz
              else
                  x=coord(1, i)/nelx; y=coord(2, i)/nely; z=0
              endif
              coord(1, i)=x*(xdims(2)-xdims(1))+xdims(1)
              coord(2, i)=y*(ydims(2)-ydims(1))+ydims(1)
              coord(3, i)=z*(ztop-zbottom)+zbottom
!                 end if
        end do
    else
        if(read_external_grid_flg==0) then !don't adjust if mesh read from file
          do i=1, npoin
              x=coord(1, i)/nelx; y=coord(2, i)/nely; z=coord(3, i)/nelz
              coord(1, i)=x*(xdims(2)-xdims(1))+xdims(1)
              coord(2, i)=y*(ydims(2)-ydims(1))+ydims(1)
              coord(3, i)=z*(ztop-zbottom)+zbottom
          end do
        end if
        
    endif

    !multirate stuff
    npartition = sum(plist)

    deallocate(partition)
    allocate(partition(npartition))

    j=1
    do i=1, nelem
      if(plist(i)>0) then
         partition(j) = i
         j=j+1
      end if
    end do

    !initialize cg coordinate
    call mod_grid_init_coord_dg_to_cg() !temporary fix - need to change if using non-conforming CG

    !initialize some data
    call mod_p4est_init_data()

    if (irank == irank0) then
       print*,'--------------------------'
       print*,'P4est Boundary Info'
       print*,'--------------------------'       
       print*,' Boundary Conditions are: '
       do i=1,6
          print*,' i, iboundary = ',i,iboundary(i)
       end do
       print*,' xperiodic = ',xperiodic
       print*,' yperiodic = ',yperiodic
       print*,' zperiodic = ',zperiodic
!!$       print*,'--------------------------'
!!$       print*,' P4est Grid Info: '
!!$       print*,'--------------------------'
!!$       print*,' npoin npoin_cg nelem = ',npoin,npoin_cg,nelem
!!$       print*,' nbsido nboun nface   = ',nbsido,nboun,nface
       print*,'--------------------------'
       print*,'--------------------------'
    end if


    deallocate(csp_face)

  end subroutine mod_p4est_create_grid

!--------------------------------------------------------------------!
!>@brief Initialization of some data after extracting p4est grid
!--------------------------------------------------------------------!
  subroutine mod_p4est_init_data

    use mod_input, only: lread_bc

    use mod_bc, only: read_bc

    implicit none

    integer :: ip, i, j, or, iboun, AllocateStatus, ierr
    
    !----------------------------
    ! read boundary conditions
    !----------------------------
    if (lread_bc) then
    
       call read_bc(face, nface, bsido, nbsido)
        
    end if

    !-------------------
    ! init ippoin_proc
    !-------------------
    ipoin_proc = 1
    if(space_method /= 'dg') then
       do i=1, num_send_recv_total
          ipoin_proc(nbh_send_recv(i)) = ipoin_proc(nbh_send_recv(i)) + 1
       end do !i
    end if

    !--------------------------------------------
    ! Initialize some values when using p4est:
    !--------------------------------------------

    !Determine Periodicity
    xperiodic=.false.
    yperiodic=.false.
    zperiodic=.false.

    !THIS Works: Note that iboundary == 30 in order to let COUNT_FACES count
    !the right number that CREATE_FACE counts since Periodic Faces are counted twice
    !since they are seen by CREATE_FACE as boundary faces but in reality are interior faces.
    if (iboundary(1) == 3 .and. iboundary(2) == 3) zperiodic=.true.
    if (iboundary(3) == 3 .and. iboundary(4) == 3) yperiodic=.true.
    if (iboundary(5) == 3 .and. iboundary(6) == 3) xperiodic=.true.

    !-----------------------
    !    init global data
    !-----------------------

    call mpi_allgather(npoin, 1, mpi_integer, npoin_l, 1, mpi_integer, mpi_comm_world, ierr)
    call mpi_allgather(ncol, 1, mpi_integer, ncol_l, 1, mpi_integer, mpi_comm_world, ierr)
    call mpi_allgather(nelem, 1, mpi_integer, nelem_l, 1, mpi_integer, mpi_comm_world, ierr)

    ncol_g=sum(ncol_l)
    npoin_g=sum(npoin_l)
    nelem_g=sum(nelem_l)

    ncol_l_max = maxval(ncol_l)
    npoin_l_max = maxval(npoin_l)
    nelem_l_max = maxval(nelem_l)

  end subroutine mod_p4est_init_data

!--------------------------------------------------------------------!
!>@brief Static Adaptive Mesh Refinement
!--------------------------------------------------------------------!

  subroutine set_refine_coarsen_flags()
    implicit none

    integer :: e, nr, r

    integer, save :: cr = 0

    !find elements to refine
    deallocate(refine_coarsen_elements)
    allocate(refine_coarsen_elements(nelem))

    nr=0
    do e=1, nelem
      r = init_ref_crit(e)
      if(r == 1) nr = nr + 1
      if(cr == 1) r = -1
      refine_coarsen_elements(e) = r
    end do

    nr = nr * 2 !assume 2x more than what we asked for
    nr = nr * 8 + (nelem - nr) !correct for 2D with 4 children

  end subroutine

!--------------------------------------------------------------------!
!>@brief Set refinement criteriea
!--------------------------------------------------------------------!
  integer function init_ref_crit(e)

    use mod_basis, only: nglx, ngly, nglz

    use mod_grid, only: coord, intma
    use mod_input, only: max_mesh_lvl, nc_box_invert

    implicit none

    integer :: e

    real :: x(nglx, ngly, nglz), y(nglx, ngly, nglz), z(nglx, ngly, nglz)

    real :: xmin, xmax, ymin, ymax, zmin, zmax

    integer :: i, j, k, ip, iref

    do k=1, nglz
       do j=1, ngly
          do i=1, nglx
             ip=intma(i, j, k, e)
             x(i, j, k) = coord(1, ip)
             y(i, j, k) = coord(2, ip)
             z(i, j, k) = coord(3, ip)

          end do
       end do
    end do

    xmin = minval(x); xmax = maxval(x)
    ymin = minval(y); ymax = maxval(y)
    zmin = minval(z); zmax = maxval(z)

    iref = 0
    if(.not. nc_box_invert) then
      if(  (xmax >= xlim_min .and. xmin <= xlim_max) .and. &
        (ymax >= ylim_min .and. ymin <= ylim_max) .and. &
        (zmax >= zlim_min .and. zmin <= zlim_max)) then

        if(max_mesh_lvl .eq. 0) then
          iref = 1
        else
          iref = max_mesh_lvl
        endif
      end if
    else
      if(  (xmax >= xlim_min .and. xmin <= xlim_max) .and. &
        (ymax >= ylim_min .and. ymin <= ylim_max) .and. &
        (zmax >= zlim_min .and. zmin <= zlim_max)) then
        iref = 0
      else
        if(max_mesh_lvl .eq. 0) then
          iref = 1
        else
          iref = max_mesh_lvl
        endif
      end if
    end if

    init_ref_crit = iref

 end function init_ref_crit

!--------------------------------------------------------------------!
!>@brief Create 2D projection matrices
!--------------------------------------------------------------------!
 subroutine create_2d_projection_matrices_numa2d(Psg, plane)

    use mod_legendre, only: legendre_gauss_lobatto

    use mod_types, only: r8, r4

    implicit none

    real, dimension(:,:,:), allocatable, intent(out):: Psg
    integer, intent(in):: plane

    integer :: i, j, k, l, mm, nn, p, q, ierr, AllocateStatus
    real    :: s, o1, o2
    real(kind=r8), dimension(:), allocatable :: wq1, wq2, xq1, xq2, xq3, xq4, wgl1, wgl2, xgl1, xgl2
    real(kind=r8), dimension(:, :), allocatable :: L1, L2, L3, L4, L1o, L2o
    real, dimension(:, :), allocatable :: M               !mass matrix

    !scatter and gather matrices
    integer :: ngl1, ngl2, nq1, nq2, nngl
    if (plane == 1) then
       ngl1=nglx
       ngl2=ngly
    elseif (plane == 2) then
       ngl1=nglx
       ngl2=nglz
    elseif (plane == 3) then
       ngl1=ngly
       ngl2=nglz
    end if

    nq1=ngl1+1
    nq2=ngl2+1

    nngl=ngl1*ngl2

    allocate(Psg(nngl,nngl,8), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_AMR_Create **"

    allocate(M(nngl, nngl), &
            xgl1(ngl1), xgl2(ngl2), &
            wgl1(ngl1), wgl2(ngl2), &
            xq1(nq1), xq2(nq2), &
            wq1(nq1), wq2(nq2), &
            xq3(nq1), xq4(nq2), &
            L1o(ngl1, nq1), L2o(ngl2, nq2), &
            L1(ngl1, nq1), L2(ngl2, nq2), &
            L3(ngl1, nq1), L4(ngl2, nq2), &
            stat=AllocateStatus)
    if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_AMR_Create **"

    !generate LGL quadrature points
    call legendre_gauss_lobatto(ngl1, xgl1, wgl1)
    call legendre_gauss_lobatto(ngl2, xgl2, wgl2)
    call legendre_gauss_lobatto(nq1, xq1, wq1)
    call legendre_gauss_lobatto(nq2, xq2, wq2)

    !generate LGL values on quadrature points
    call lagrange_poly(L1o, xq1, xgl1, nq1, ngl1)
    call lagrange_poly(L2o, xq2, xgl2, nq2, ngl2)

    !set offset and scale values for interface
    s = 0.5
    o1 = -0.5
    o2 = 0.5

    !offset xq values for children elements on the interface
    xq3 = xq1
    xq4 = xq2
    do i=1, nq1
       xq1(i) = o1+s*xq1(i)
       xq3(i) = o2+s*xq3(i)
    end do
    do i=1, nq2
       xq2(i) = o1+s*xq2(i)
       xq4(i) = o2+s*xq4(i)
    end do

    call lagrange_poly(L1, xq1, xgl1, nq1, ngl1)
    call lagrange_poly(L2, xq2, xgl2, nq2, ngl2)
    call lagrange_poly(L3, xq3, xgl1, nq1, ngl1)
    call lagrange_poly(L4, xq4, xgl2, nq2, ngl2)

    !initialize integral matrices
    M = 0
    Psg = 0

    !compute integral matrices
    do l=1, ngl1
    do k=1, ngl2
       mm=(l-1)*ngl2+k
       do j=1, ngl1
       do i=1, ngl2
          nn=(j-1)*ngl2+i
          do q=1, nq1
          do p=1, nq2
             M(nn, mm) = M(nn, mm) + L2o(i, p)*L1o(j, q)*L2o(k, p)*L1o(l, q)&
                  *wq2(p)*wq1(q)
             !order of subfaces in (ngl1,ngl2) is: --,-+,+-,++
             Psg(nn, mm, 1) = Psg(nn, mm, 1) +  L2(i, p)*L1(j, q)*L2o(k, p)*L1o(l, q)&
                  *wq2(p)*wq1(q)
             Psg(nn, mm, 2) = Psg(nn, mm, 2) +  L4(i, p)*L1(j, q)*L2o(k, p)*L1o(l, q)&
                  *wq2(p)*wq1(q)
             Psg(nn, mm, 3) = Psg(nn, mm, 3) +  L2(i, p)*L3(j, q)*L2o(k, p)*L1o(l, q)&
                  *wq2(p)*wq1(q)
             Psg(nn, mm, 4) = Psg(nn, mm, 4) +  L4(i, p)*L3(j, q)*L2o(k, p)*L1o(l, q)&
                  *wq2(p)*wq1(q)
             Psg(nn, mm, 5) = Psg(nn, mm, 5) + L2o(i, p)*L1o(j, q)*L2(k, p)*L1(l, q)&
                  *wq2(p)*wq1(q)*s*s
             Psg(nn, mm, 6) = Psg(nn, mm, 6) + L2o(i, p)*L1o(j, q)*L4(k, p)*L1(l, q)&
                  *wq2(p)*wq1(q)*s*s
             Psg(nn, mm, 7) = Psg(nn, mm, 7) + L2o(i, p)*L1o(j, q)*L2(k, p)*L3(l, q)&
                  *wq2(p)*wq1(q)*s*s
             Psg(nn, mm, 8) = Psg(nn, mm, 8) + L2o(i, p)*L1o(j, q)*L4(k, p)*L3(l, q)&
                  *wq2(p)*wq1(q)*s*s
          end do
          end do
       end do
       end do
    end do
    end do

    !invert M
    call gaujordf(M, nngl, ierr)

    !compute projection matrices!!!
    do i = 1, 8
        Psg(:,:,i) = matmul(Psg(:,:,i), M)
    end do

    !fix for 2d
    if(ngl2 == 1) then
        Psg(:,:,1) = 0.5*(Psg(:,:,1)+Psg(:,:,2))
        Psg(:,:,2) = 0.5*(Psg(:,:,3)+Psg(:,:,4))
        Psg(:,:,3) =     (Psg(:,:,5)+Psg(:,:,6))
        Psg(:,:,4) =     (Psg(:,:,7)+Psg(:,:,8))
    else if(ngl1 == 1) then
        Psg(:,:,1) = 0.5*(Psg(:,:,1)+Psg(:,:,3))
        Psg(:,:,2) = 0.5*(Psg(:,:,2)+Psg(:,:,4))
        Psg(:,:,3) =     (Psg(:,:,5)+Psg(:,:,7))
        Psg(:,:,4) =     (Psg(:,:,6)+Psg(:,:,8))
    endif

    !test gather-scatter
!    call test_projection_matrices_numa2d(ngl1,ngl2,plane)

    !free
    deallocate(xq1, wq1, xq2, wq2, xq3, xq4, &
               xgl1, wgl1, xgl2, wgl2, &
               L1o, L2o, L1, L2, L3, L4, M)

  end subroutine create_2d_projection_matrices_numa2d

!--------------------------------------------------------------------!
!>@brief Scatter to one subface
!--------------------------------------------------------------------!
  subroutine scatter_element_2d_subface(qe, qec, ic, ngl1, ngl2, plane)

    implicit none

    integer, intent(in):: ngl1, ngl2, plane
    integer, intent(in) :: ic
    integer j, k, n
    real, dimension(ngl, ngl), intent(in) :: qe
    real, dimension(ngl, ngl), intent(out) :: qec
    real, dimension(1, ngl1*ngl2) :: qa

    do j=1, ngl2
       do k=1, ngl1
          n=(j-1)*ngl1+k
          qa(1, n) = qe(k, j)
       end do
    end do

    if(plane == 1) then
        qa = matmul(qa, PsgXY(:,:,ic))
    else if(plane == 2) then
        qa = matmul(qa, PsgXZ(:,:,ic))
    else
        qa = matmul(qa, PsgYZ(:,:,ic))
    endif

    do j=1, ngl2
       do k=1, ngl1
          n=(j-1)*ngl1+k
          qec(k, j) = qa(1, n)
       end do
    end do


  end subroutine scatter_element_2d_subface

!--------------------------------------------------------------------!
!>@brief Gather from one subface
!--------------------------------------------------------------------!
  subroutine gather_element_2d_subface(qec, qe, ic, ngl1, ngl2, plane)

    implicit none

    integer, intent(in):: ngl1, ngl2, plane
    integer, intent(in) :: ic
    integer i, j, k, n
    real, dimension(ngl, ngl) :: qe
    real, dimension(ngl, ngl) :: qec
    real, dimension(1, ngl1*ngl2) :: qa

    do j=1, ngl2
       do k=1, ngl1
          n=(j-1)*ngl1+k
          qa(1, n) = qec(k, j)
       end do
    end do

    if(plane == 1) then
        qa = matmul(qa, PsgXY(:,:,ic + FACE_CHILDREN))
    else if(plane == 2) then
        qa = matmul(qa, PsgXZ(:,:,ic + FACE_CHILDREN))
    else
        qa = matmul(qa, PsgYZ(:,:,ic + FACE_CHILDREN))
    endif

    do j=1, ngl2
       do k=1, ngl1
          n=(j-1)*ngl1+k
          qe(k, j) = qa(1, n)
       end do
    end do

  end subroutine gather_element_2d_subface

!--------------------------------------------------------------------!
!>@brief Scatter fields to children faces
!--------------------------------------------------------------------!
  subroutine scatter_element_2d(qe, qec, ngl1, ngl2, plane)

    implicit none

    real, dimension(ngl, ngl) :: qe
    real, dimension(ngl, ngl, FACE_CHILDREN) :: qec
    integer, intent(in):: ngl1, ngl2, plane
    integer ic

    do ic=1, FACE_CHILDREN
        call scatter_element_2d_subface(qe, qec(:,:,ic), ic, ngl1, ngl2, plane)
    enddo

  end subroutine scatter_element_2d

!--------------------------------------------------------------------!
!>@brief Gather fields from children faces
!--------------------------------------------------------------------!
  subroutine gather_element_2d(qec, qe, ngl1, ngl2, plane)

    implicit none

    real, dimension(ngl, ngl) :: qe
    real, dimension(ngl, ngl, FACE_CHILDREN) :: qec
    real, dimension(ngl, ngl) :: qe_
    integer, intent(in):: ngl1, ngl2, plane
    integer ic

    qe = 0
    do ic=1, FACE_CHILDREN
        call gather_element_2d_subface(qec(:,:,ic), qe_, ic, ngl1, ngl2, plane)
        qe = qe + qe_
    enddo

  end subroutine gather_element_2d

!--------------------------------------------------------------------!
!>@brief Test
!--------------------------------------------------------------------!
  subroutine test_projection_matrices_numa2d(ngl1, ngl2, plane)

    implicit none
    integer, intent(in):: ngl1, ngl2, plane
    character :: fmt*10, ngl_str*2
    real, dimension(ngl1, ngl2) :: u1, u3
    real, dimension(ngl1, ngl2, FACE_CHILDREN) :: u2
    integer i, j, ic

    u1=1
    u3=0

    do j=1, ngl2
       do i=1, ngl1
          u1(i, j) = 0.25*(i-1.0)/(ngl1) + 0.75*(j-1.0)/(ngl2)
       end do
    end do

    call scatter_element_2d(u1, u2, ngl1, ngl2, plane)
    call gather_element_2d(u2, u3, ngl1, ngl2, plane)

    write(ngl_str, '(I2)')ngl1*ngl2
    fmt='('//trim(ngl_str)//'e16.4)'

    print*, "PsgXY"
    do i = 1,FACE_CHILDREN*2
        do j=1, ngl1*ngl2
           write(*, fmt)PsgXY(j, :, i)
        end do
        print*, ""
    enddo

    write(ngl_str, '(I2)') ngl1
    fmt='('//trim(ngl_str)//'e16.4)'

    print*, "ORIGINAL"
    do j=1, ngl1
       write(*, fmt)u1(j, :)
    end do

    print*, "SCATTER"
    do ic=1, FACE_CHILDREN
       print*, "Block", ic
       do j=1, ngl1
          write(*, fmt)u2(j, :, ic)
       end do
    end do

    print*, "GATHER"
    do j=1, ngl1
       write(*, fmt)u3(j, :)
    end do

!    stop

  end subroutine test_projection_matrices_numa2d

  !----------------------------------------------------------------------!
  !> @brief generates interpolation of LGL polynomials defined by xgl on xq quadrature points
  !> @author M.A.Kopera on 09 Feb 2012
  !----------------------------------------------------------------------!
  subroutine lagrange_poly(L, xq, xgl, nq, ngl)

    use mod_types, only: r8, r4

    integer :: nq, ngl
    integer :: i, j, k
    real(kind=r8), dimension(ngl, nq) :: L
    real(kind=r8), dimension(ngl)    :: xgl
    real(kind=r8), dimension(nq)     :: xq

    L = 1
    do i=1, ngl
       do j=1, ngl
          if(i.ne.j) then
             do k=1, nq
                L(i, k) = L(i, k)*(xq(k)-xgl(j))/(xgl(i)-xgl(j))
             end do
          endif
       end do
    end do

  end subroutine lagrange_poly

  ! }}}

  ! {{{ Adaptivity Routines

  integer function binary_search(table, val)

    use iso_c_binding, only: C_INT64_T

    implicit none

    integer(C_INT64_T), dimension(:), intent(in) :: table
    integer(C_INT64_T), intent(in) :: val
    integer :: l
    integer :: r
    integer :: m

    l = 1
    r = size(table,1)
    do while(r > l+1)
      m = (r+l) / 2
      if(table(m) == val) then
        do while(m < size(table, 1))
          if(table(m+1) .ne. val) exit
          m = m + 1
        enddo
        binary_search = m
        return
      elseif(table(m) < val) then
        l = m
      else
        r = m
      end if
    enddo

    binary_search = l
  end function binary_search

  subroutine mod_p4est_mark_elements(q, lvl)

    use mod_input, only: amr_mark_max_min, amr_mark_modes, amr_mark_threshold, &
      amr_mark_random

    use iso_c_binding, only: C_INT8_T

    implicit none

    real, dimension(:, :), allocatable, intent(in):: q
    integer(C_INT8_T), dimension(:) :: lvl

    if (amr_mark_random) then
      call mod_p4est_mark_elements_random(q, lvl)
    endif
    if (amr_mark_max_min) then
      call mod_p4est_mark_elements_max_min(q, lvl)
    endif
    if (amr_mark_modes) then
      call mod_p4est_mark_elements_modes(q, lvl)
    endif
    if (amr_mark_threshold) then
      call mod_p4est_mark_threshold(q, lvl)
    endif
    ! call mod_p4est_mark_elements_region(q, lvl)

  end subroutine mod_p4est_mark_elements

  subroutine mod_p4est_adapt(q, q1, q2, q3, q4, q5, q6)

    use iso_c_binding, only: C_INT32_T, C_INT8_T, C_CHAR, C_NULL_CHAR

    use mod_input, only: amr_num_neigh_iter, amr_mark_random, amr_mark_set2nc, &
      eqn_set

    use mod_initial, only: nvar, q_ref

    use mod_grid, only: npoin

    implicit none
    real, dimension(:, :), allocatable, intent(inout):: q
    real, dimension(:, :), allocatable, intent(inout), optional :: q1, q2, q3, &
                                                                   q4, q5, q6
    real, dimension(:, :, :), allocatable :: qs, q_recv
    real, dimension(:, :), allocatable :: q_aux
    integer(C_INT32_T) :: num_loc_elem
    integer(C_INT32_T) :: num_neigh_iter

    integer(C_INT32_T), dimension(:), allocatable :: dst_to_src_elem_id,       &
                                                     coarse_dst_to_src_elem_id

    integer(C_INT8_T), dimension(:), allocatable :: lvl_src, lvl_dst

    integer :: sz1, sz2, loc_npoin, rk, ierr, tag=777
    integer :: comm_start, comm_end, comm_size, comm_oset
    integer :: recv_rank_start, recv_rank_end
    integer :: send_rank_start, send_rank_end
    integer, dimension(:), allocatable :: recv_requests, send_requests
    integer, dimension(:,:), allocatable :: recv_status, send_status

    if(.not.allocated(qid_src)) allocate(qid_src(0:nproc))
    if(.not.allocated(qid_dst)) allocate(qid_dst(0:nproc))

    num_neigh_iter = amr_num_neigh_iter

    if (amr_mark_random) then
      num_neigh_iter = 0
    endif

    if(irank == irank0) print*, '----------- Begin mesh adapt -------------'

    sz2 = 1
    if(present(q1)) sz2 = sz2+1
    if(present(q2)) sz2 = sz2+1
    if(present(q3)) sz2 = sz2+1
    if(present(q4)) sz2 = sz2+1
    if(present(q5)) sz2 = sz2+1
    if(present(q6)) sz2 = sz2+1

    sz1 = size(q, 1)

    num_loc_elem = nelem

    !---------------------------!
    ! coarsen / refine Solution !
    !---------------------------!

    !------------------------!
    ! (1) Get current levels !
    !------------------------!
    allocate(lvl_src(num_loc_elem))
    call p8esttonuma_get_element_lvl(lvl_src, num_loc_elem)

    !----------------------------------!
    ! (2) Mark elements for refinement !
    !----------------------------------!

    call mod_p4est_mark_elements(q, lvl_src)


    !----------------------------------------!
    ! (3) Coarsen, refine, and balance p4est !
    !----------------------------------------!
    call p8esttonuma_mark_neighbors_p4est(num_neigh_iter)
    call p8esttonuma_coarsen_refine_p4est(num_loc_elem)


    !--------------------!
    ! (4) Get new levels !
    !--------------------!
    allocate(lvl_dst(num_loc_elem))
    call p8esttonuma_get_element_lvl(lvl_dst, num_loc_elem)

    !-----------------------!
    ! (5) Transfer solution !
    !-----------------------!
    loc_npoin = num_loc_elem * npts
    allocate(qs(sz1, sz2, loc_npoin))

    call mod_p4est_transfer_q_3d(qs, q, 1, lvl_src, lvl_dst)
    if(present(q1)) call mod_p4est_transfer_q_3d(qs, q1, 2, lvl_src, lvl_dst)
    if(present(q2)) call mod_p4est_transfer_q_3d(qs, q2, 3, lvl_src, lvl_dst)
    if(present(q3)) call mod_p4est_transfer_q_3d(qs, q3, 4, lvl_src, lvl_dst)
    if(present(q4)) call mod_p4est_transfer_q_3d(qs, q4, 5, lvl_src, lvl_dst)
    if(present(q5)) call mod_p4est_transfer_q_3d(qs, q5, 6, lvl_src, lvl_dst)
    if(present(q6)) call mod_p4est_transfer_q_3d(qs, q6, 7, lvl_src, lvl_dst)

    deallocate(lvl_src)
    deallocate(lvl_dst)

    !------------------------------------!
    ! (6) Get the new parallel partition !
    !------------------------------------!
    call p8esttonuma_repartition(qid_src, qid_dst)
    ! if(irank == irank0) then
    !   print*,">>>",qid_src
    !   print*,">>>",qid_dst
    ! endif
    ! if(irank == irank0) then
    !   do rk = 0,nproc-1
    !     print*, rk
    !     print*, qid_src(rk), qid_src(rk+1)
    !     print*, qid_dst(rk), qid_dst(rk+1)
    !   enddo
    ! endif

    !-----------------------------------!
    ! (7) Parallel repartition solution !
    !-----------------------------------!

    !------------!
    ! Post recvs !
    !------------!
    recv_rank_start = binary_search(qid_src, qid_dst(irank  )    ) - 1
    recv_rank_end   = binary_search(qid_src, qid_dst(irank+1) - 1) - 1
    ! print*,"recv: ", irank, recv_rank_start, recv_rank_end

    allocate(recv_requests(recv_rank_start:recv_rank_end))
    allocate(recv_status(MPI_Status_size, recv_rank_start:recv_rank_end))

    ! determine the size of our comm buffer
    ! first determine our own overlap
    comm_start = max(qid_src(irank  ), qid_dst(irank  ))
    comm_end   = min(qid_src(irank+1), qid_dst(irank+1))
    comm_size  = max(0, comm_end - comm_start)

    ! total size - overlap size
    comm_size = qid_dst(irank+1) - qid_dst(irank) - comm_size

    allocate(q_recv(sz1, sz2, comm_size * npts))

    ! Actually post recvs. comm_oset: starting elements - 1 of recv array
    comm_oset = 0
    do rk = recv_rank_start,recv_rank_end
      comm_start = max(qid_dst(irank    ), qid_src(rk    ))
      comm_end   = min(qid_dst(irank + 1), qid_src(rk + 1))
      comm_size  = comm_end - comm_start

      if(rk .ne. irank) then
        ! Issue MPI_Irecv
        if(comm_size > 0) then
          call MPI_Irecv(q_recv(1, 1, comm_oset * npts + 1),                   &
            comm_size * npts * sz1 * sz2 , MPI_PRECISION, rk, tag,             &
            MPI_COMM_WORLD, recv_requests(rk), ierr)
          ! print*,"recv:", irank, rk, comm_oset, comm_start, comm_end
          comm_oset = comm_oset + comm_size
        else
          recv_requests(rk) = MPI_REQUEST_NULL
        end if
      else
        ! Set MPI_REQUEST_NULL for local wait
        recv_requests(rk) = MPI_REQUEST_NULL
      end if
    end do
    if (npts * comm_oset .ne. size(q_recv,3)) then
      print*, irank, "recv wrong", comm_oset, size(q_recv,3)
      stop "something went wrong with adapt repartition recv"
    endif

    !------------!
    ! Post sends !
    !------------!
    ! Determine the ranks I send to
    send_rank_start = binary_search(qid_dst, qid_src(irank    )    ) - 1
    send_rank_end   = binary_search(qid_dst, qid_src(irank + 1) - 1) - 1
    ! print*,"send: ", irank, send_rank_start, send_rank_end
    allocate(send_requests(send_rank_start:send_rank_end))
    allocate(send_status(MPI_Status_size, send_rank_start:send_rank_end))

    ! Actually post send. comm_oset: starting elements - 1 of src array
    comm_oset = 0
    do rk = send_rank_start,send_rank_end
      comm_start = max(qid_src(irank    ), qid_dst(rk    ))
      comm_end   = min(qid_src(irank + 1), qid_dst(rk + 1))
      comm_size  = comm_end - comm_start

      if(rk .ne. irank) then
        ! Issue MPI_Isend
        if(comm_size > 0) then
          call MPI_Isend(qs(1, 1, comm_oset * npts + 1),                       &
            comm_size * npts * sz1 * sz2 , MPI_PRECISION, rk, tag,             &
            MPI_COMM_WORLD, send_requests(rk), ierr)
          ! print*,"send:", irank, rk, comm_oset, comm_start, comm_end
          comm_oset = comm_oset + comm_size
        else
          send_requests(rk) = MPI_REQUEST_NULL
        endif
      else
        ! Set MPI_REQUEST_NULL for local wait
        send_requests(rk) = MPI_REQUEST_NULL
        comm_oset = comm_oset + comm_size
      end if

    end do
    if (npts * comm_oset .ne. size(qs,3)) then
      print*, irank, "send wrong", comm_oset, size(qs,3)
      stop "something went wrong with adapt repartition send"
    endif

    ! Wait on the recvs
    call MPI_Waitall(size(recv_requests,1), recv_requests, recv_status, ierr)

    !--------------------------------!
    ! (7) Set the new solution !
    !--------------------------------!
    num_loc_elem = qid_dst(irank + 1) - qid_dst(irank)
    loc_npoin = num_loc_elem * npts
    allocate(q(sz1, loc_npoin))
    call mod_p4est_adapt_set_q(q, qs, q_recv, 1, recv_rank_start,              &
      recv_rank_end, qid_src, qid_dst)
    if(present(q1)) then
      allocate(q1(sz1, loc_npoin))
      call mod_p4est_adapt_set_q(q1, qs, q_recv, 2, recv_rank_start,           &
        recv_rank_end, qid_src, qid_dst)
    endif
    if(present(q2)) then
      allocate(q2(sz1, loc_npoin))
      call mod_p4est_adapt_set_q(q2, qs, q_recv, 3, recv_rank_start,           &
        recv_rank_end, qid_src, qid_dst)
    endif
    if(present(q3)) then
      allocate(q3(sz1, loc_npoin))
      call mod_p4est_adapt_set_q(q3, qs, q_recv, 4, recv_rank_start,           &
        recv_rank_end, qid_src, qid_dst)
    endif
    if(present(q4)) then
      allocate(q4(sz1, loc_npoin))
      call mod_p4est_adapt_set_q(q4, qs, q_recv, 5, recv_rank_start,           &
        recv_rank_end, qid_src, qid_dst)
    endif
    if(present(q5)) then
      allocate(q5(sz1, loc_npoin))
      call mod_p4est_adapt_set_q(q5, qs, q_recv, 6, recv_rank_start,           &
        recv_rank_end, qid_src, qid_dst)
    endif
    if(present(q6)) then
      allocate(q6(sz1, loc_npoin))
      call mod_p4est_adapt_set_q(q6, qs, q_recv, 7, recv_rank_start,           &
        recv_rank_end, qid_src, qid_dst)
    endif

    ! Wait on sends before deallocating qs
    call MPI_Waitall(size(send_requests,1), send_requests, send_status, ierr)
    deallocate(qs)

    if(irank == irank0) print*, '----------- End mesh adapt   -------------'

    ! FIXME: Can we make these routines more efficient?
    call mod_p4est_create_grid(.false.)
    call mod_face_create_nc_list()

  end subroutine mod_p4est_adapt

  subroutine mod_p4est_mark_elements_random(q, lvl)

    use iso_c_binding, only: C_INT, C_INT8_T
    use mod_input, only: max_mesh_lvl
#ifdef __INTEL_COMPILER
    use ifport
#endif

    implicit none

    real, dimension(:, :), allocatable, intent(in):: q
    integer(C_INT), dimension(:), allocatable :: hadapt
    integer(C_INT8_T), dimension(:) :: lvl
    integer :: i
    integer :: a, h
#ifdef __PGI
    double precision :: rand
#endif

    allocate(hadapt(nelem))
    hadapt = 0
    call srand(irank)

    a = 1
    h = 0

    do i = 1, nelem
      if(mod(i, a) == 0) then
        a = floor(rand(0) * 13 + 8)
        h = floor(rand(0) * 3) - 1
      endif
      if(h > 0 .and. lvl(i) .ge. max_mesh_lvl) then
        hadapt(i) = 0
      else
        hadapt(i) = h
      endif
    end do

    if(init_refine) then
      do i = 1, nelem
        hadapt(i) = max(0, hadapt(i))
      enddo
    endif
    call p8esttonuma_mark_elements(hadapt)
    deallocate(hadapt)

  end subroutine mod_p4est_mark_elements_random

  subroutine mod_p4est_mark_elements_max_min(q, lvl)

    use iso_c_binding, only: C_INT, C_INT8_T
    use mod_input, only: max_mesh_lvl, amr_max_min_lim
    use mod_initial, only: nvar

    implicit none

    real, dimension(:, :), allocatable, intent(in):: q
    integer(C_INT), dimension(:), allocatable :: hadapt
    integer(C_INT8_T), dimension(:) :: lvl
    integer :: i, k
    real :: qmax, qmin

    allocate(hadapt(nelem))
    hadapt = 0

    do i = 1, nelem
      hadapt(i) = 0
      qmax = q(1, (i-1)*npts + 1)
      qmin = q(1, (i-1)*npts + 1)
      do k = 1,npts
        qmax = max(qmax, q(1, (i-1)*npts + k))
        qmin = min(qmin, q(1, (i-1)*npts + k))
      enddo
      if ((qmax-qmin) > amr_max_min_lim(1)/(2**lvl(i)) .and. &
          lvl(i) < max_mesh_lvl) then
        hadapt(i) = 1
      elseif ((qmax-qmin) < amr_max_min_lim(2)/(2**lvl(i)) ) then
        hadapt(i) = -1
      endif
    enddo

    if(init_refine) then
      do i = 1, nelem
        hadapt(i) = max(0, hadapt(i))
      enddo
    endif
    call p8esttonuma_mark_elements(hadapt)
    deallocate(hadapt)

  end subroutine mod_p4est_mark_elements_max_min

  subroutine mod_p4est_mark_threshold(q, lvl)

    use iso_c_binding, only: C_INT, C_INT8_T
    use mod_input, only: max_mesh_lvl, amr_threshold_lim, eqn_set, limit_threshold
    use mod_initial, only: nvar, q_ref

    implicit none

    real, dimension(:, :), allocatable, intent(in):: q
    integer(C_INT), dimension(:), allocatable :: hadapt
    integer(C_INT8_T), dimension(:) :: lvl
    integer :: i, k, j, kv, p
    real :: qval

    allocate(hadapt(nelem))
    hadapt = -1

    do kv = 1,size(amr_indicator_variables, 1)
      j = amr_indicator_variables(kv)
      if(j == 0) exit
      do i = 1, nelem
        do k = 1,npts
          qval = q(j, (i-1)*npts + k)
          if(j == 5) then
            if(eqn_set(1:5) == 'set2c') then
              p = (i-1)*npts + k
              qval = (q(5,p) + q_ref(5,p)) / (q(1,p) + q_ref(1,p)) &
                   - q_ref(5,p)/q_ref(1,p)
            else
              stop "mod_p4est: temp threshold only set up for set2c"
            endif
          endif
          if(abs(qval) > amr_threshold_lim(kv)) then
            if(lvl(i) < max_mesh_lvl) then
              hadapt(i) = 1
            else
              hadapt(i) = 0
            endif
          endif
        enddo
      enddo
    enddo

    if(init_refine) then
      do i = 1, nelem
        hadapt(i) = max(0, hadapt(i))
      enddo
    endif
    call p8esttonuma_mark_elements(hadapt)
    deallocate(hadapt)

  end subroutine mod_p4est_mark_threshold

  subroutine mod_p4est_vandermonde_Legendre(V, x)

    implicit none

    real, dimension(:,0:), intent(out) :: V
    real, dimension(:), intent(in) :: x
    integer :: Nrp, n, k

    Nrp = size(x, 1)
    n = size(V, 2) - 1

    ! edge cases
    V(:, 0) = 1
    if(n .ge. 1) then
      V(:, 1) = x
    endif

    ! recursion
    do k = 1, n-1
      V(:, k + 1) = ((2 * k + 1) * x * V(:, k) - k * V(:, k-1)) / (k + 1)
    enddo

    ! normalize
    do k = 0, n
      V(:, k) = V(:, k) * sqrt(real(2 * k + 1) / 2)
    enddo

  end subroutine mod_p4est_vandermonde_Legendre

  subroutine mod_p4est_build_Vinv(Vx_inv, x)

    use mod_legendre, only: legendre_gauss

    implicit none

    real, dimension(:, :), allocatable, intent(out) :: Vx_inv
    real, dimension(:), intent(in) :: x
    real, dimension(:), allocatable :: r_gl
    real, dimension(:), allocatable :: w_gl
    real, dimension(:), allocatable :: wx
    real, dimension(:, :), allocatable :: I_x2gl
    real, dimension(:, :), allocatable :: V_gl
    real, dimension(:, :), allocatable :: M
    integer :: Nrp, k

    ! Number of row points
    Nrp = size(x, 1)

    ! Get the Gauss points
    allocate(r_gl(Nrp))
    allocate(w_gl(Nrp))
    call legendre_gauss(Nrp, r_gl, w_gl)

    ! Build exact mass matrix
    allocate(M(Nrp, Nrp))
    M = 0
    do k = 1,Nrp
      M(k,k) = w_gl(k)
    enddo

    ! build interpolation matrix
    allocate(I_x2gl(Nrp, Nrp))
    allocate(wx(Nrp))
    call mod_p4est_barycentric_weights(x, wx)
    call mod_p4est_build_interpolation(I_x2gl, x, r_gl, wx)

    ! build Vandermonde
    allocate(V_gl(Nrp, Nrp))
    call mod_p4est_vandermonde_Legendre(V_gl, r_gl)

    ! Build the inverse Vandermonde
    allocate(Vx_inv(Nrp, Nrp))
    Vx_inv = matmul(transpose(V_gl), matmul(M, I_x2gl))

  end subroutine mod_p4est_build_Vinv

  subroutine mod_p4est_build_perfect_modal_coeff(b, m)
    implicit none

    real, dimension(:), allocatable, intent(out) :: b
    integer, intent(in) :: m
    real :: bN
    integer :: k

    allocate(b(m))
    bN = 0
    do k = 1, m-1
      bN = bN + 1 / real(k**(2*(m-1)))
    enddo
    bN = sqrt(bN)
    b(1) = 0
    do k = 1, m-1
      b(k+1) = 1 / (k ** (m-1) * bN)
    enddo

  end subroutine mod_p4est_build_perfect_modal_coeff

  subroutine mod_p4est_smoothness(s, q, p, ne, std, m, V_inv, b)
    use mod_input, only: amr_mark_modes_use_baseline_decay
    implicit none
    real, intent(inout) :: s
    real, dimension(:,:), intent(in) :: q
    real, dimension(0:), intent(out) :: p
    integer, intent(in) :: ne
    integer, intent(in) :: std
    integer, intent(in) :: m
    real, dimension(0:,0:), intent(in) :: V_inv
    real, dimension(0:), intent(in) :: b

    real :: qL2, pavg, xavg, numer, denom, tmp
    integer :: i, k, j
    real, dimension(0:m-1) :: qt

    do k = 1,size(amr_indicator_variables, 1)
      j = amr_indicator_variables(k)

      ! When this is zero we're done!
      if(j == 0) return

      do i = 0,m-1
        qt(i) = q(j, ne + i * std)
      enddo

      p = matmul(V_inv,qt)
      qL2 = sum(p**2)

      if(qL2 .le. amr_smoothness_qL2_limit(k)) cycle

      ! skyline pessimization
      if(amr_mark_modes_use_baseline_decay) then
        do i = 1,m-1
          p(i) = log10(sqrt(p(i)**2 + qL2 * b(i)**2))
        enddo
      else
        do i = 1,m-1
          p(i) = log10(abs(p(i)))
        enddo
      endif

      p(m-1) = max(p(m-1), p(m-2))

      pavg = p(m-1)
      xavg = log10(real(m-1))
      do i = m-2,1,-1
        p(i) = max(p(i+1), p(i))
        pavg = pavg + p(i)
        xavg = xavg + log10(real(i))
      enddo
      pavg = pavg / (m-1)
      xavg = xavg / (m-1)

      ! get the smoothness
      numer = 0
      denom = 0
      do i = 1,m-1
        tmp = log10(real(i)) - xavg
        numer = numer + (p(i) - pavg) * tmp
        denom = denom + tmp * tmp
      enddo
      s = min(s, -numer/denom)
    enddo

  end subroutine mod_p4est_smoothness

  ! This indicator is based on
  ! @article{klockner2011viscous,
  !   title={Viscous shock capturing in a time-explicit discontinuous Galerkin method},
  !   author={Kl{\"o}ckner, Andreas and Warburton, Tim and Hesthaven, Jan S},
  !   journal={Mathematical Modelling of Natural Phenomena},
  !   volume={6},
  !   number={3},
  !   pages={57--83},
  !   year={2011},
  !   publisher={EDP Sciences}
  ! }
  subroutine mod_p4est_mark_elements_modes(q, lvl)

    use iso_c_binding, only: C_INT, C_INT8_T
    use mod_input, only: max_mesh_lvl
    use mod_initial, only: nvar

    implicit none

    real, dimension(:, :), allocatable, intent(in):: q
    integer(C_INT), dimension(:), allocatable :: hadapt
    integer(C_INT8_T), dimension(:) :: lvl
    integer :: i, ix, iy, iz, k, ne, j, std
    real :: qmax, qmin
    real :: sx, sy, sz
    real :: qx(0:nglx-1), qy(0:ngly-1), qz(0:nglz-1)
    real :: qN


    ! Build the nodes to modes operators if necessary
    if(.not. allocated(Vx_inv)) then
      call mod_p4est_build_Vinv(Vx_inv, xglx)
      call mod_p4est_build_Vinv(Vy_inv, xgly)
      call mod_p4est_build_Vinv(Vz_inv, xglz)

      call mod_p4est_build_perfect_modal_coeff(bx, nglx)
      call mod_p4est_build_perfect_modal_coeff(by, ngly)
      call mod_p4est_build_perfect_modal_coeff(bz, nglz)
    endif

    allocate(hadapt(nelem))
    hadapt = 0

    do i = 1, nelem
      hadapt(i) = -1

      ! smoothness in x direction
      sx = nglx-1
      std = 1
      ix = 0
      do iz = 0,nglz-1
        do iy = 0,ngly-1
          ne = (i-1) * npts + ix + nglx * (iy + iz * ngly) + 1
          call mod_p4est_smoothness(sx, q, qx, ne, std, nglx, Vx_inv, bx)
        enddo
      enddo

      ! smoothness in y direction
      sy = ngly-1
      iy = 0
      std = nglx
      do iz = 0,nglz-1
        do ix = 0,nglx-1
          ne = (i-1) * npts + ix + nglx * (iy + iz * ngly) + 1
          call mod_p4est_smoothness(sy, q, qy, ne, std, ngly, Vy_inv, by)
        enddo
      enddo

      ! smoothness in z direction
      sz = nglz-1
      iz = 0
      std = nglx * ngly
      do iy = 0,ngly-1
        do ix = 0,nglx-1
          ne = (i-1) * npts + ix + nglx * (iy + iz * ngly) + 1
          call mod_p4est_smoothness(sz, q, qz, ne, std, nglz, Vz_inv, bz)
        enddo
      enddo

      if(min(sx,sy,sz) < amr_smoothness_limits(1) .and. (lvl(i) < max_mesh_lvl)) then
        hadapt(i) = 1
      elseif (min(sx,sy,sz) > amr_smoothness_limits(2)) then
        hadapt(i) = -1
      else
        hadapt(i) = 0
      endif
    enddo

    if(init_refine) then
      do i = 1, nelem
        hadapt(i) = max(0, hadapt(i))
      enddo
    endif
    call p8esttonuma_mark_elements(hadapt)
    deallocate(hadapt)

  end subroutine mod_p4est_mark_elements_modes


  subroutine mod_p4est_mark_elements_region(q, lvl)

    use iso_c_binding, only: C_INT, C_INT8_T
    use mod_input, only: max_mesh_lvl

    implicit none

    real, dimension(:, :), allocatable, intent(in):: q
    integer(C_INT8_T), dimension(:) :: lvl
    integer(C_INT), dimension(:), allocatable :: hadapt
    integer :: i

    allocate(hadapt(nelem))
    hadapt = 0

    do i = 1, nelem
      if (lvl(i) < max_mesh_lvl) then
        hadapt(i) = init_ref_crit(i)
      endif
    enddo

    if(init_refine) then
      do i = 1, nelem
        hadapt(i) = max(0, hadapt(i))
      enddo
    endif
    call p8esttonuma_mark_elements(hadapt)
    deallocate(hadapt)

  end subroutine mod_p4est_mark_elements_region

  subroutine mod_p4est_barycentric_weights(ra, wa)
    implicit none
    real, dimension(:), intent(in) :: ra
    real, dimension(:), intent(out) :: wa
    integer :: Nra, k, j

    Nra = size(ra, 1)

    do k = 1,Nra
      wa(k) = 1
      do j = 1,Nra
        if(j .ne. k) wa(k) = wa(k) * (ra(k) - ra(j))
      enddo
      wa(k) = 1 / wa(k)
    enddo

  end subroutine  mod_p4est_barycentric_weights

  subroutine mod_p4est_build_interpolation(I_a2b, ra, rb, wa)
    implicit none
    real, dimension(:, :), intent(out) :: I_a2b
    real, dimension(:), intent(in) :: ra, rb, wa
    integer :: Nra, Nrb, k, j, i
    real :: d

    Nra = size(ra,1)
    Nrb = size(rb,1)

    do k = 1,Nrb
      d = 0
      do j = 1,Nra
        if(approx_equal(rb(k), ra(j), sqrt(epsilon(real(1))), 0)) then
          do i = 1,Nra
            I_a2b(k, i) = 0.
          enddo
          I_a2b(k, j) = 1
          d  = 1
          exit
        end if
        d = d + wa(j) / (rb(k) - ra(j))
        I_a2b(k, j) = wa(j) / (rb(k) - ra(j))
      end do
      do j = 1,Nra
        I_a2b(k, j) = I_a2b(k, j) / d
      end do
    end do
  end subroutine  mod_p4est_build_interpolation

  subroutine mod_p4est_build_projection_1d(interp, project, ra)

    use mod_legendre, only: legendre_gauss

    implicit none

    real, dimension(:, :, :), allocatable, intent(out) :: interp
    real, dimension(:, :, :), allocatable, intent(out) :: project
    real, dimension(:), intent(in) :: ra
    real, dimension(:), allocatable :: ra_b, ra_t, wa, rb, wb, wb_gl
    integer :: Nrp, k
    real, dimension(:, :), allocatable :: tmp, I_a2b, I_b2a, M, MI

    ! Uncomment to check the projection
    ! real, dimension(:, :), allocatable :: chk
    ! real :: s
    ! integer :: j

    ! Number of points in stencil
    Nrp = size(ra,1)

    ! Create the barycentric weights
    allocate(wa(Nrp))
    call mod_p4est_barycentric_weights(ra, wa)

    ! Create the top and bottom grids
    allocate(ra_b(Nrp))
    allocate(ra_t(Nrp))
    do k = 1,Nrp
      ra_t(k) = (ra(k) + 1) / 2
      ra_b(k) = (ra(k) - 1) / 2
    enddo

    ! Create the interpolation matrices
    allocate(tmp(Nrp, Nrp))
    allocate(interp(Nrp, Nrp, 2))

    call mod_p4est_build_interpolation(tmp, ra, ra_b, wa)
    interp(:, :, 1) = tmp

    call mod_p4est_build_interpolation(tmp, ra, ra_t, wa)
    interp(:, :, 2) = tmp

    ! Build the mass matrix and its inverse
    allocate(rb(Nrp))
    allocate(wb_gl(Nrp))
    call legendre_gauss(Nrp, rb, wb_gl)

    allocate(wb(Nrp))
    call mod_p4est_barycentric_weights(rb, wb)

    allocate(I_a2b(Nrp, Nrp))
    call mod_p4est_build_interpolation(I_a2b, ra, rb, wa)

    allocate(I_b2a(Nrp, Nrp))
    call mod_p4est_build_interpolation(I_b2a, rb, ra, wb)

    allocate(M(Nrp, Nrp))
    allocate(MI(Nrp, Nrp))
    M = 0
    MI = 0

    do k = 1,Nrp
      M(k, k) = wb_gl(k)
      MI(k, k) = 1/wb_gl(k)
    end do

    M = matmul(transpose(I_a2b), matmul(M, I_a2b))
    MI = matmul(I_b2a, matmul(MI, transpose(I_b2a)))

    ! Create the projection matrices
    allocate(project(Nrp, Nrp, 2))

    project(:, :, 1) = matmul(MI, matmul(transpose(interp(:, :, 1)), M)) / 2
    project(:, :, 2) = matmul(MI, matmul(transpose(interp(:, :, 2)), M)) / 2

    ! Free up the memory we used
    deallocate(wa)
    deallocate(ra_b)
    deallocate(ra_t)
    deallocate(rb)
    deallocate(wb)
    deallocate(wb_gl)
    deallocate(I_a2b)
    deallocate(I_b2a)
    deallocate(M)
    deallocate(MI)

    ! Uncomment to check the projection
    ! allocate(chk(Nrp, Nrp))
    ! chk = 0
    ! do k = 1,Nrp
    !   chk(k,k) = 1
    ! end do
    ! chk = matmul(project(:, :, 2), interp(:, :, 2)) + &
    !       matmul(project(:,:,1), interp(:,:,1)) - &
    !       chk
    ! s = 0
    ! do k = 1,Nrp
    !   do j = 1,Nrp
    !     s = s + chk(k,j)**2
    !   end do
    ! end do
    ! print*,"If projection correct this should be close to zero: ", sqrt(s)
    ! deallocate(chk)

  end subroutine mod_p4est_build_projection_1d

  subroutine mod_p4est_transfer_q_3d(q_dst, q_src, n, lvl_src, lvl_dst)

    use iso_c_binding, only : C_INT8_T

    implicit none

    real, dimension(:, :, :), intent(inout) :: q_dst
    real, dimension(:, :), allocatable, intent(inout) :: q_src
    integer, intent(in) :: n
    integer(C_INT8_T), dimension(:), intent(in) :: lvl_src, lvl_dst
    integer :: num_elem, k_src, k_dst, o_src, o_dst
    integer :: ix, iy, iz
    integer :: jx, jy, jz
    integer :: cx, cy, cz
    integer :: nx, ny, nz, np
    integer :: f, nf
    real, dimension(:, :), allocatable :: qx, qxy


    nx=nglx
    ny=ngly
    nz=nglz
    np = nx * ny * nz
    nf = size(q_dst, 1)

    allocate(qx(nf, np))
    allocate(qxy(nf, np))

    if(.not.allocated(interp_x)) then
      call mod_p4est_build_projection_1d(interp_x, project_x, xglx)
      call mod_p4est_build_projection_1d(interp_y, project_y, xgly)
      call mod_p4est_build_projection_1d(interp_z, project_z, xglz)
    endif

    num_elem = size(lvl_dst, 1)

    k_dst = 1; o_dst = 0
    k_src = 1; o_src = 0

    do while (k_dst <= num_elem)
      ! Source and Destination are same level
      if (lvl_src(k_src) == lvl_dst(k_dst)) then
        do ix = 1,np
          do f = 1,nf
            q_dst(f, n, o_dst+ix) = q_src(f, o_src+ix)
          enddo
        enddo

        k_dst = k_dst + 1; o_dst = o_dst + np
        k_src = k_src + 1; o_src = o_src + np

      ! Destination is hanging
      elseif (lvl_src(k_src) < lvl_dst(k_dst)) then

        ! loop over children
        do cz = 1, 2
          do cy = 1, 2
            do cx = 1, 2

              qx = 0
              qxy = 0
              do ix = 1,np
                do f = 1,nf
                  q_dst(f, n, o_dst+ix) = 0
                enddo
              enddo

              ! interpolation in x: (I \otimes I \otimes A) * v
              do iz = 0,nz-1
                do iy = 0,ny-1
                  do ix = 0,nx-1
                    do jx = 0,nx-1
                      ! loop over fields
                      do f = 1,nf
                        qx(f, (iz * ny + iy) * nx + jx + 1) =                  &
                          qx(f, (iz * ny + iy) * nx + jx + 1) +                &
                          interp_x(jx + 1, ix + 1, cx) *                       &
                          q_src(f, o_src + (iz * ny + iy) * nx + ix + 1)
                        enddo
                    enddo
                  enddo
                enddo
              enddo

              ! interpolation in y: (I \otimes A \otimes I) * v
              do iz = 0,nz-1
                do iy = 0,ny-1
                  do jy = 0,ny-1
                    do ix = 0,nx-1
                      ! loop over fields
                      do f = 1,nf
                        qxy(f, (iz * ny + jy) * nx + ix + 1) =                 &
                          qxy(f, (iz * ny + jy) * nx + ix + 1) +               &
                          interp_y(jy + 1, iy + 1, cy) *                       &
                          qx(f, (iz * ny + iy) * nx + ix + 1)
                        enddo
                    enddo
                  enddo
                enddo
              enddo

              ! interpolation in z: (A \otimes I \otimes I) * v
              do iz = 0,nz-1
                do jz = 0,nz-1
                  do iy = 0,ny-1
                    do ix = 0,nx-1
                      ! loop over fields
                      do f = 1,nf
                        q_dst(f, n, o_dst + (jz * ny + iy) * nx + ix + 1) =    &
                          q_dst(f, n, o_dst + (jz * ny + iy) * nx + ix + 1) +  &
                          interp_z(jz + 1, iz + 1, cz) *                       &
                          qxy(f, (iz * ny + iy) * nx + ix + 1)
                        enddo
                    enddo
                  enddo
                enddo
              enddo

              k_dst = k_dst + 1; o_dst = o_dst + np
            enddo
          enddo
        enddo

        k_src = k_src + 1; o_src = o_src + np

      ! Source is hanging
      else
        ! Zero out the destination
        do ix = 1,np
          do f = 1,nf
            q_dst(f, n, o_dst+ix) = 0
          enddo
        enddo
        ! loop over children
        do cz = 1, 2
          do cy = 1, 2
            do cx = 1, 2

              qx = 0
              qxy = 0

              ! projection in x: (I \otimes I \otimes A) * v
              do iz = 0,nz-1
                do iy = 0,ny-1
                  do ix = 0,nx-1
                    do jx = 0,nx-1
                      ! loop over fields
                      do f = 1,nf
                        qx(f, (iz * ny + iy) * nx + jx + 1) =                  &
                          qx(f, (iz * ny + iy) * nx + jx + 1) +                &
                          project_x(jx + 1, ix + 1, cx) *                      &
                          q_src(f, o_src + (iz * ny + iy) * nx + ix + 1)
                      enddo
                    enddo
                  enddo
                enddo
              enddo

              ! projection in y: (I \otimes A \otimes I) * v
              do iz = 0,nz-1
                do iy = 0,ny-1
                  do jy = 0,ny-1
                    do ix = 0,nx-1
                      ! loop over fields
                      do f = 1,nf
                        qxy(f, (iz * ny + jy) * nx + ix + 1) =                 &
                          qxy(f, (iz * ny + jy) * nx + ix + 1) +               &
                          project_y(jy + 1, iy + 1, cy) *                      &
                          qx(f, (iz * ny + iy) * nx + ix + 1)
                      enddo
                    enddo
                  enddo
                enddo
              enddo

              ! projection in z: (A \otimes I \otimes I) * v
              do iz = 0,nz-1
                do jz = 0,nz-1
                  do iy = 0,ny-1
                    do ix = 0,nx-1
                      ! loop over fields
                      do f = 1,nf
                        q_dst(f, n, o_dst + (jz * ny + iy) * nx + ix + 1) =    &
                          q_dst(f, n, o_dst + (jz * ny + iy) * nx + ix + 1) +  &
                          project_z(jz + 1, iz + 1, cz) *                      &
                          qxy(f, (iz * ny + iy) * nx + ix + 1)
                      enddo
                    enddo
                  enddo
                enddo
              enddo

              k_src = k_src + 1; o_src = o_src + np
            enddo
          enddo
        end do
        k_dst = k_dst + 1; o_dst = o_dst + np
      endif
    end do

    deallocate(q_src)

  end subroutine mod_p4est_transfer_q_3d

  subroutine mod_p4est_adapt_set_q(q, qs, qr, n, recv_rank_start,       &
      recv_rank_end, qid_src, qid_dst)

    use iso_c_binding, only: C_INT64_T

    implicit none

    real, intent(out) :: q(:,0:)
    real, intent(in) :: qs(:,:,0:)
    real, intent(in) :: qr(:,:,0:)
    integer, intent(in) :: n
    integer, intent(in) :: recv_rank_start
    integer, intent(in) :: recv_rank_end
    integer(C_INT64_T), intent(in):: qid_src(0:), qid_dst(0:)
    integer(C_INT64_T) :: comm_start, comm_end, comm_size
    integer :: q_start, q_end, r_start, r_end, s_start, s_end

    q_start = 0;
    r_start = 0;

    ! lo ranks
    if(recv_rank_start < irank) then
      comm_start = qid_dst(irank)
      comm_end   = min(qid_dst(irank+1), qid_src(irank))
      comm_size = comm_end - comm_start

      if(comm_size > 0) then
        q_end = q_start + comm_size * npts
        r_end = r_start + comm_size * npts
        q(:, q_start : q_end - 1) = qr(:, n, r_start : r_end - 1)
        q_start = q_end
        r_start = r_end
      endif
    endif

    ! self
    if(recv_rank_start <= irank .and. irank <= recv_rank_end) then
      comm_start = max(qid_dst(irank  ), qid_src(irank  ))
      comm_end   = min(qid_dst(irank+1), qid_src(irank+1))
      comm_size = comm_end - comm_start

      if(comm_size > 0) then
        q_end = q_start + comm_size * npts
        s_start = (comm_start - qid_src(irank)) * npts
        s_end = s_start + comm_size * npts
        q(:, q_start : q_end - 1) = qs(:, n, s_start : s_end - 1)
        q_start = q_end
      endif
    endif

    ! hi
    if(irank < recv_rank_end) then
      comm_end   = qid_dst(irank + 1)
      comm_start = max(qid_dst(irank), qid_src(irank+1))
      comm_size = comm_end - comm_start

      if(comm_size > 0) then
        q_end = q_start + comm_size * npts
        r_end = r_start + comm_size * npts
        q(:, q_start : q_end - 1) = qr(:, n, r_start : r_end - 1)
        q_start = q_end
        r_start = r_end
      endif
    endif

    if(r_start .ne. size(qr,3)) stop "problem with repartition recv"
    if(q_start .ne. size(q, 2)) stop "problem with repartition total"

  end subroutine mod_p4est_adapt_set_q

  ! }}}

end module mod_p4est
