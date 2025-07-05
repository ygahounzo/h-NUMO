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

    use mod_constants, only: tol, gravity

    use mod_grid, only: nelem,  npoin, ncol, nz, node_column, nboun, npoin_cg

    use mod_initial, only: q_ref, nvar, nvart, kvector 

    use mod_input, only: space_method, nlayers, is_mlswe

    use mod_parallel, only: num_send_recv_total
    
    use mod_basis, only: ngl, nq

    public :: mod_ref_create, &
        recv_data, &
        grad_qv_ref, grad_qc_ref, grad_qr_ref, q_send, q_recv, &
        f0_ref, g0_ref, grad_g0_ref, d_g0_ref_dr, recv_data_dg, send_data_dg, nmessage, &
        grad_bathy, grad_hB_ref, grad_phiA_ref, grad_rho_ref_layers, &
        q_recv_quad, q_send_quad, recv_data_dg_quad, send_data_dg_quad, &
        lap_recv_data_dg_df1, lap_send_data_dg_df1, lap_q_recv_df1, lap_q_send_df1

    private

    !module variables and parameters
    real,    dimension(:,:), allocatable :: qb, grad_rho_ref, grad_theta_ref, grad_salinity_ref
    real,    dimension(:,:), allocatable :: grad_press_ref, grad_g0_ref
    real,    dimension(:,:), allocatable :: grad_qv_ref, grad_qc_ref, grad_qr_ref, grad_bathy
    real,    dimension(:,:,:), allocatable :: grad_hB_ref, grad_phiA_ref, grad_rho_ref_layers
    real,    dimension(:),   allocatable :: norm_inf_br, norm_inf_bu, norm_inf_bv, norm_inf_bw
    real,    dimension(:),   allocatable :: norm_inf_bt, norm_inf_bene
    real,    dimension(:,:), allocatable :: norm_inf_brhs
    real,    dimension(:,:), allocatable :: d_rho_ref_dr, d_theta_ref_dr, d_g0_ref_dr
    real,    dimension(:),   allocatable :: press_ref, g0, f0, h0, press, div_u_ref, f0_ref
    real,    dimension(:),   allocatable :: g0_ref, dens_var
    real,    dimension(:,:), allocatable :: recv_data
    real,    dimension(:,:,:), allocatable :: q_recv, q_send
    real,    dimension(:), allocatable :: recv_data_dg, send_data_dg
    real,    dimension(:,:,:), allocatable :: q_recv_quad, q_send_quad, lap_q_recv_df1
    real,    dimension(:,:,:), allocatable :: lap_q_send_df1
    real,    dimension(:), allocatable :: recv_data_dg_quad, send_data_dg_quad
    real,    dimension(:), allocatable :: lap_recv_data_dg_df1, lap_send_data_dg_df1
    integer :: nmessage

contains

    subroutine mod_ref_create()

        implicit none
        integer i, j, ix, iz, ip, k
        integer AllocateStatus

        real, dimension(:,:), allocatable :: f, dfdz, q_tempv
        real, dimension(:), allocatable :: q_temp
        real dpds(3), hmatrix(3,3), dpdr, dpdt, max_dp1, max_dp2

        !Size of Max DG Message
        nmessage=2*nvar+4 !nvar+3 for inviscid dynamics
        if(is_mlswe) nmessage = 6 ! 2 for uv-momentum

        if(allocated(qb)) then
            deallocate(qb, press_ref, press, dens_var, recv_data, norm_inf_brhs, norm_inf_bene, &
                norm_inf_br, norm_inf_bu, norm_inf_bv, norm_inf_bw, norm_inf_bt, &
                grad_press_ref, grad_rho_ref, grad_theta_ref, grad_salinity_ref,div_u_ref, &
                q_send, q_recv, grad_bathy)
        endif
        if(allocated(q_recv_quad)) then
            deallocate(q_recv_quad, q_send_quad)
        end if
        allocate( qb(nvar,npoin), press_ref(npoin),  press(npoin), dens_var(npoin), &
            recv_data(nvar,num_send_recv_total), &
            norm_inf_brhs(nvar,nelem), norm_inf_bene(nelem), &
            norm_inf_br(nelem), norm_inf_bu(nelem), norm_inf_bv(nelem), norm_inf_bw(nelem), &
            norm_inf_bt(nelem), grad_press_ref(3,npoin), grad_rho_ref(3,npoin), &
            grad_theta_ref(3,npoin), grad_salinity_ref(3,npoin), div_u_ref(npoin), &
            q_send(4,ngl,nboun),q_recv(4,ngl,nboun), grad_bathy(3,npoin), &
            q_recv_quad(4,nq,nboun), q_send_quad(4,nq,nboun),&
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 0**"

        if(allocated(lap_q_recv_df1)) then 
            deallocate(lap_q_recv_df1, lap_q_send_df1) 
        endif 
        allocate(lap_q_recv_df1(4,ngl,nboun), lap_q_send_df1(4,ngl,nboun), &
            stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 0**"

        if(allocated(lap_recv_data_dg_df1)) then
            deallocate(lap_recv_data_dg_df1, lap_send_data_dg_df1)
        endif
        allocate(lap_recv_data_dg_df1(4*ngl*nboun), lap_send_data_dg_df1(4*ngl*nboun), &
            stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 0**"

        if (space_method(1:2) == 'dg') then
            if(allocated(recv_data_dg)) then
                deallocate(recv_data_dg, send_data_dg)
            endif

            if(allocated(recv_data_dg_quad)) then
                deallocate(recv_data_dg_quad, send_data_dg_quad)
            endif

            allocate( recv_data_dg(4*ngl*nboun), &
                send_data_dg(4*ngl*nboun), &
                recv_data_dg_quad(4*nq*nboun), &
                send_data_dg_quad(4*nq*nboun), &
                stat=AllocateStatus )
            if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 1**"
        end if

        !Initialize allocated arrays:
        qb             = 0.0
        grad_press_ref = 0.0
        grad_rho_ref   = 0.0
        grad_bathy     = 0.0
        
        !Construct NRBC Boundary Vector
        do i=1,npoin
            qb(2,i)=q_ref(2,i)
            qb(3,i)=q_ref(3,i)
            qb(4,i)=q_ref(4,i)
        end do                    !i

    end subroutine mod_ref_create

end module mod_ref
