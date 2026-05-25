!----------------------------------------------------------------------!
!>@brief This module builds the Reference Solution (Gradients and Pressure)
!>@author  Francis X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@ modified by Yao Gahounzo
!>      Computing PhD
!       Boise State University
!       Date: April 06, 2023
!       modified for exact integration ( added nq ad is_mlswe condtions)
!----------------------------------------------------------------------!
module mod_ref

    use mod_constants
    use mod_grid
    use mod_initial
    use mod_input
    use mod_parallel
    use mod_basis

    implicit none

    type :: mref
        real, dimension(:,:),   allocatable :: qb, grad_rho_ref, grad_theta_ref, grad_salinity_ref
        real, dimension(:,:),   allocatable :: grad_press_ref, grad_g0_ref
        real, dimension(:,:),   allocatable :: grad_qv_ref, grad_qc_ref, grad_qr_ref, grad_bathy
        real, dimension(:,:,:), allocatable :: grad_hB_ref, grad_phiA_ref, grad_rho_ref_layers
        real, dimension(:),     allocatable :: norm_inf_br, norm_inf_bu, norm_inf_bv, norm_inf_bw
        real, dimension(:),     allocatable :: norm_inf_bt, norm_inf_bene
        real, dimension(:,:),   allocatable :: norm_inf_brhs
        real, dimension(:,:),   allocatable :: d_rho_ref_dr, d_theta_ref_dr, d_g0_ref_dr
        real, dimension(:),     allocatable :: press_ref, g0, f0, h0, press, div_u_ref, f0_ref
        real, dimension(:),     allocatable :: g0_ref, dens_var
        real, dimension(:,:),   allocatable :: recv_data
        real, dimension(:,:,:), allocatable :: q_recv, q_send, q_recv_lap, q_send_lap
        real, dimension(:,:,:), allocatable :: q_send_bcl, q_recv_bcl, q_send_lap_bcl, q_recv_lap_bcl
        real, dimension(:),     allocatable :: recv_data_dg, send_data_dg, recv_data_dg_lap, send_data_dg_lap
        real, dimension(:),     allocatable :: recv_data_bcl, send_data_bcl, recv_data_lap_bcl, send_data_lap_bcl
        real, dimension(:,:,:), allocatable :: q_recv_quad, q_send_quad, lap_q_recv_df1, lap_q_send_df1
        real, dimension(:),     allocatable :: recv_data_dg_quad, send_data_dg_quad
        real, dimension(:),     allocatable :: lap_recv_data_dg_df1, lap_send_data_dg_df1
        real, dimension(:,:,:), allocatable :: q_send_csty, q_recv_csty
        real, dimension(:),     allocatable :: recv_data_csty, send_data_csty
        integer :: nmessage, nbtp_var, nboun_valid
        integer, allocatable :: face_pack_list(:)
    end type mref

    public :: mod_ref_create, mref

    private

contains

    subroutine mod_ref_create(G, inp, init, par, b, ref)

        implicit none

        type(grid), intent(in) :: G
        type(input), intent(in) :: inp
        type(initial), intent(in) :: init
        type(parallel_CS), intent(in) :: par
        type(basis), intent(in) :: b
        type(mref), intent(inout) :: ref

        integer i, j, ix, iz, ip, k, jj, inbh, ib, iface, imulti
        integer AllocateStatus

        real, dimension(:,:), allocatable :: f, dfdz, q_tempv
        real, dimension(:), allocatable :: q_temp
        real dpds(3), hmatrix(3,3), dpdr, dpdt, max_dp1, max_dp2

        !Size of Max DG Message
        ref%nmessage = 2*init%nvar + 4 !nvar+3 for inviscid dynamics
        if(inp%is_mlswe) ref%nmessage = 6 ! 2 for uv-momentum
        ref%nbtp_var = 3*inp%nlayers + 5 ! 3*nlayers for baroclinic variables + 5 for barotropic variables

        ! Compact list of valid MPI faces (face_type==2, imulti>0) in nbh_send_recv order.
        ! Computed once here so pack/unpack GPU kernels need no per-call pre-scans.
        if (allocated(ref%face_pack_list)) deallocate(ref%face_pack_list)
        allocate(ref%face_pack_list(G%nboun))
        ref%nboun_valid = 0
        jj = 1
        do inbh = 1, par%num_nbh
            do ib = 1, par%num_send_recv(inbh)
                iface  = par%nbh_send_recv(jj)
                imulti = par%nbh_send_recv_multi(jj)
                if (G%face_type(iface) == 2 .and. imulti > 0) then
                    ref%nboun_valid = ref%nboun_valid + 1
                    ref%face_pack_list(ref%nboun_valid) = iface
                end if
                jj = jj + 1
            end do
        end do

        if(allocated(ref%qb)) then
            deallocate(ref%qb, ref%press_ref, ref%press, ref%dens_var, ref%recv_data, &
                ref%norm_inf_brhs, ref%norm_inf_bene, &
                ref%norm_inf_br, ref%norm_inf_bu, ref%norm_inf_bv, ref%norm_inf_bw, ref%norm_inf_bt, &
                ref%grad_press_ref, ref%grad_rho_ref, ref%grad_theta_ref, ref%grad_salinity_ref, &
                ref%div_u_ref, ref%grad_bathy)
        endif
        if(allocated(ref%q_recv_quad)) then
            deallocate(ref%q_recv_quad, ref%q_send_quad)
        end if
        allocate( ref%qb(init%nvar,G%npoin), ref%press_ref(G%npoin), &
            ref%press(G%npoin), ref%dens_var(G%npoin), &
            ref%recv_data(init%nvar,G%nboun), &
            ref%norm_inf_brhs(init%nvar,G%nelem), ref%norm_inf_bene(G%nelem), &
            ref%norm_inf_br(G%nelem), ref%norm_inf_bu(G%nelem), &
            ref%norm_inf_bv(G%nelem), ref%norm_inf_bw(G%nelem), &
            ref%norm_inf_bt(G%nelem), ref%grad_press_ref(3,G%npoin), &
            ref%grad_rho_ref(3,G%npoin), ref%grad_theta_ref(3,G%npoin), &
            ref%grad_salinity_ref(3,G%npoin), ref%div_u_ref(G%npoin), &
            ref%grad_bathy(3,G%npoin), &
            ref%q_recv_quad(4,b%nq,G%nboun), ref%q_send_quad(4,b%nq,G%nboun), &
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 0**"

        if(allocated(ref%lap_q_recv_df1)) then
            deallocate(ref%lap_q_recv_df1, ref%lap_q_send_df1)
        endif
        allocate(ref%lap_q_recv_df1(5,b%ngl,G%nboun), ref%lap_q_send_df1(5,b%ngl,G%nboun), &
            stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 0**"

        if(allocated(ref%lap_recv_data_dg_df1)) then
            deallocate(ref%lap_recv_data_dg_df1, ref%lap_send_data_dg_df1)
        endif
        allocate(ref%lap_recv_data_dg_df1(5*b%ngl*G%nboun), &
            ref%lap_send_data_dg_df1(5*b%ngl*G%nboun), &
            stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 0**"

        if(allocated(ref%q_send)) then
            deallocate(ref%q_send, ref%q_recv, ref%q_send_lap, ref%q_recv_lap)
        endif
        allocate(ref%q_send(ref%nbtp_var,b%ngl,G%nboun), ref%q_recv(ref%nbtp_var,b%ngl,G%nboun), &
            ref%q_send_lap(10,b%ngl,G%nboun), ref%q_recv_lap(10,b%ngl,G%nboun), &
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 0**"

        ! if (inp%space_method(1:2) == 'dg') then
            if(allocated(ref%recv_data_dg)) then
                deallocate(ref%recv_data_dg, ref%send_data_dg, &
                    ref%recv_data_dg_lap, ref%send_data_dg_lap)
            endif
            allocate( ref%recv_data_dg(ref%nbtp_var*b%ngl*G%nboun), &
                ref%send_data_dg(ref%nbtp_var*b%ngl*G%nboun), &
                ref%recv_data_dg_lap(10*b%ngl*G%nboun), &
                ref%send_data_dg_lap(10*b%ngl*G%nboun), &
                stat=AllocateStatus )
            if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 1**"

            if(allocated(ref%recv_data_dg_quad)) then
                deallocate(ref%recv_data_dg_quad, ref%send_data_dg_quad)
            endif
            allocate(ref%recv_data_dg_quad(4*b%nq*G%nboun), &
                ref%send_data_dg_quad(4*b%nq*G%nboun), &
                stat=AllocateStatus )
            if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 1**"
        ! end if

        if(allocated(ref%q_send_bcl)) then
            deallocate(ref%q_send_bcl, ref%q_recv_bcl, ref%q_send_lap_bcl, ref%q_recv_lap_bcl)
        endif
        allocate(ref%q_send_bcl(3*inp%nlayers,b%ngl,G%nboun), &
            ref%q_recv_bcl(3*inp%nlayers,b%ngl,G%nboun), &
            ref%q_send_lap_bcl(5*inp%nlayers,b%ngl,G%nboun), &
            ref%q_recv_lap_bcl(5*inp%nlayers,b%ngl,G%nboun), &
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 0**"

        if(allocated(ref%recv_data_bcl)) then
            deallocate(ref%recv_data_bcl, ref%send_data_bcl, &
                ref%recv_data_lap_bcl, ref%send_data_lap_bcl)
        endif
        allocate( ref%recv_data_bcl(3*inp%nlayers*b%ngl*G%nboun), &
            ref%send_data_bcl(3*inp%nlayers*b%ngl*G%nboun), &
            ref%recv_data_lap_bcl(5*inp%nlayers*b%ngl*G%nboun), &
            ref%send_data_lap_bcl(5*inp%nlayers*b%ngl*G%nboun), &
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 1**"

        if(allocated(ref%q_send_csty)) then
            deallocate(ref%q_send_csty, ref%q_recv_csty)
        endif
        allocate(ref%q_send_csty(inp%nlayers+1,b%ngl,G%nboun), &
            ref%q_recv_csty(inp%nlayers+1,b%ngl,G%nboun), &
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 0**"

        if(allocated(ref%recv_data_csty)) then
            deallocate(ref%recv_data_csty, ref%send_data_csty)
        endif
        allocate( ref%recv_data_csty((inp%nlayers+1)*b%ngl*G%nboun), &
            ref%send_data_csty((inp%nlayers+1)*b%ngl*G%nboun), &
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 1**"

        !Initialize allocated arrays:
        ref%qb             = 0.0
        ref%grad_press_ref = 0.0
        ref%grad_rho_ref   = 0.0
        ref%grad_bathy     = 0.0

        !Construct NRBC Boundary Vector
        do i=1,G%npoin
            ref%qb(2,i) = init%q_ref(2,i)
            ref%qb(3,i) = init%q_ref(3,i)
            ref%qb(4,i) = init%q_ref(4,i)
        end do                    !i

    end subroutine mod_ref_create

end module mod_ref
