!----------------------------------------------------------------------!
!>@brief This module builds the Parallel (MPI) Communication
!>@author  Jim Kelly in 2010
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>@date F.X. Giraldo 1/2014 to include variable NGL orders and to use INTENT 
!>of passed arrays
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>@date S. Marras on 03/2014 to add ncol_l_max
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>@date S. Marras on 06/2014 to add nboun_poin_l
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>
!----------------------------------------------------------------------!
module mod_parallel
  
    use mod_basis, only: nglx, ngly, nglz
  
    use mod_global_grid, only: nelem_g, npoin_g, npoin_g_cg, nface_g, ncol_g, intma_g, &
        nelem_s, ele_col_g, flag_periodic_g, nz, nboun_poin_g
  
    use mod_input, only: geometry_type, nelx, nely, nelz
  
    public :: &
        ! Constructor Functions
        mod_parallel_reorder,&

        ! Variables
        nproc, &
        nelem_l, &
        nelemx_l, &
        nelemy_l, &
        nelemz_l, &
        npoin_l, &
        npoin_l_cg, &
        ncol_l,  &
        nbsido_l, &
        nbound_poin_l, &
        nface_l, &
        nface_boundary_l, &
        nelem_l_max, &
        nelemx_l_max, &
        nelemy_l_max, &
        nelemz_l_max, &
        npoin_l_max, &
        nface_l_max, &
        ncol_l_max, &
        nboun_max, &
        nboun_poin_l_max, &
        num_nbh_max, &
        face_type_l,&
        local_global_face_l, &
        local_global_elem_l, &
        global_proc, &
        local_global_poin_l_cg, &
        local_global_poin_l, &
        local_global_poin_periodic_l, &
        flag_periodic_l, &
        bsido_l, &
        nbh_proc_l, &
        num_nbh_l, &
        ipoin_proc_l, &
        num_nbh, &
       
        local_global_elem, &
        local_global_poin, &
        local_global_poin_cg, &
        local_global_poin_periodic, &
        local_global_face, &
        flag_periodic, &
        num_proc, &
        nbh_proc, &
        ipoin_proc, &
       
        num_send_recv, &
        num_send_recv_total, &
        nbh_send_recv, &
        nbh_send_recv_multi, &
        nbh_send_recv_half, &
        npoin, &
        npoin_cg, &
        nelem, &
        nelemx,&
        nelemy,&
        nelemz,&
        nface
  
    private
    !-----------------------------------------------------------------------
    !Arrays on the head node
    integer nproc, num_nbh, num_send_total, num_send_recv_total
    integer, dimension(:),     allocatable :: nelem_l, npoin_l, npoin_l_cg, nface_l, nface_boundary_l, nbsido_l
    integer, dimension(:),     allocatable :: nboun_poin_l
    integer, dimension(:),     allocatable :: nelemx_l, nelemy_l, nelemz_l
    integer, dimension(:),     allocatable :: ncol_l
    integer, dimension(:,:),   allocatable :: local_global_elem_l, global_local_elem_l, local_global_face_l,face_type_l
    integer, dimension(:),     allocatable :: global_proc
    integer, dimension(:,:),   allocatable :: local_global_poin_l_cg, local_global_poin_l, local_global_poin_periodic_l, global_local_poin_l
    integer, dimension(:,:),   allocatable :: flag_periodic_l
    integer, dimension(:,:),   allocatable :: nbh_proc_l
    integer, dimension(:),     allocatable :: num_nbh_l
    integer, dimension(:,:,:), allocatable :: bsido_l
    integer, dimension(:,:),   allocatable :: ipoin_proc_l
  
    !Arrays on individual processors
    integer, dimension(:),   allocatable :: local_global_elem, global_local_elem
    integer, dimension(:),   allocatable :: local_global_poin_cg, local_global_poin, local_global_poin_periodic,local_global_face
    integer, dimension(:),   allocatable :: flag_periodic
    integer, dimension(:),   allocatable :: nbh_proc
    integer, dimension(:),   allocatable :: ipoin_proc
    integer, dimension(:),   allocatable :: num_send_recv
    integer, dimension(:),   allocatable :: nbh_send_recv
    integer, dimension(:),   allocatable :: nbh_send_recv_multi
    integer, dimension(:),   allocatable :: nbh_send_recv_half
    integer npoin, npoin_cg, nelem, nelemx, nelemy, nelemz
    integer ncol_l_max
  !-----------------------------------------------------------------------
  
contains
  
    !-------------------------------------------
    !-  Reorder face_send
    !-------------------------------------------
    subroutine mod_parallel_reorder(face_send)
      use mod_input, only: space_method
      use mod_mpi_utilities, only: irank

      implicit none

      integer, dimension(:):: face_send
      integer:: ind, jnd, imult, inbh, ib, iface

      if (space_method == 'dg') then
        ind = 0
        jnd = 0
        do inbh = 1,num_nbh
          do ib = 1,num_send_recv(inbh)
            ind = ind + 1

            iface = nbh_send_recv(ind)
            do imult = 1,nbh_send_recv_multi(ind)
              jnd = jnd + 1
              face_send(jnd) = iface
            enddo
          end do
        end do
        if(jnd .ne. size(face_send)) then
          print*,"on rank = ", irank
          print*,"face_send mismatch"
          print*,"filled ", jnd, " of ", size(face_send)
          stop "** mod_parallel_reorder problem"
        endif
      endif
    end subroutine
  
end module mod_parallel




