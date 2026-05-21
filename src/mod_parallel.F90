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
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: July 10, 2023
!----------------------------------------------------------------------!

module mod_parallel

    implicit none

    private

    public :: parallel_CS, mod_parallel_reorder

    !-----------------------------------------------------------------------
    type :: parallel_CS

        !-------------------------------------------------------------------
        ! Arrays on the head node
        !-------------------------------------------------------------------
        integer :: nproc, num_nbh, num_send_total, num_send_recv_total

        integer, dimension(:),     allocatable :: nelem_l, npoin_l, npoin_l_cg, nface_l
        integer, dimension(:),     allocatable :: nboun_poin_l, nface_boundary_l, nbsido_l
        integer, dimension(:),     allocatable :: nelemx_l, nelemy_l, nelemz_l
        integer, dimension(:),     allocatable :: ncol_l
        integer, dimension(:,:),   allocatable :: local_global_elem_l, global_local_elem_l
        integer, dimension(:,:),   allocatable :: local_global_face_l, face_type_l
        integer, dimension(:),     allocatable :: global_proc
        integer, dimension(:,:),   allocatable :: local_global_poin_l_cg, local_global_poin_l
        integer, dimension(:,:),   allocatable :: local_global_poin_periodic_l, global_local_poin_l
        integer, dimension(:,:),   allocatable :: flag_periodic_l
        integer, dimension(:,:),   allocatable :: nbh_proc_l
        integer, dimension(:),     allocatable :: num_nbh_l
        integer, dimension(:,:,:), allocatable :: bsido_l
        integer, dimension(:,:),   allocatable :: ipoin_proc_l

        !-------------------------------------------------------------------
        ! Arrays on individual processors
        !-------------------------------------------------------------------
        integer, dimension(:), allocatable :: local_global_elem, global_local_elem
        integer, dimension(:), allocatable :: local_global_poin_cg, local_global_poin
        integer, dimension(:), allocatable :: local_global_poin_periodic, local_global_face
        integer, dimension(:), allocatable :: flag_periodic
        integer, dimension(:), allocatable :: nbh_proc
        integer, dimension(:), allocatable :: ipoin_proc
        integer, dimension(:), allocatable :: num_send_recv
        integer, dimension(:), allocatable :: nbh_send_recv
        integer, dimension(:), allocatable :: nbh_send_recv_multi
        integer, dimension(:), allocatable :: nbh_send_recv_half

        integer :: npoin, npoin_cg, nelem, nelemx, nelemy, nelemz
        integer :: ncol_l_max
        integer :: nelem_l_max, nelemx_l_max, nelemy_l_max, nelemz_l_max
        integer :: npoin_l_max, nface_l_max, nboun_max, nboun_poin_l_max
        integer :: num_nbh_max
        integer :: nface

    end type parallel_CS

contains

    !-------------------------------------------
    !-  Reorder face_send
    !-------------------------------------------
    subroutine mod_parallel_reorder(par, face_send)

        use mod_mpi_utilities, only: irank

        type(parallel_CS), intent(inout) :: par

        integer, dimension(:), intent(inout) :: face_send

        integer :: ind, jnd, imult, inbh, ib, iface

        ind = 0
        jnd = 0
        do inbh = 1, par%num_nbh
            do ib = 1, par%num_send_recv(inbh)
                ind = ind + 1

                iface = par%nbh_send_recv(ind)
                do imult = 1, par%nbh_send_recv_multi(ind)
                    jnd = jnd + 1
                    face_send(jnd) = iface
                end do
            end do
        end do

        if (jnd /= size(face_send)) then
            print *, "on rank = ", irank
            print *, "face_send mismatch"
            print *, "filled ", jnd, " of ", size(face_send)
            stop "** mod_parallel_reorder problem"
        end if

    end subroutine mod_parallel_reorder

end module mod_parallel
