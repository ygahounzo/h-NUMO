!----------------------------------------------------------------------!
!>@brief This module builds the Reference Solution (Gradients and Pressure)
!>@author  Francis X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!----------------------------------------------------------------------!
module mod_ref

    use mod_constants, only: gamma, tol, p00, gravity

    use mod_grid, only: nelem,  npoin, ncol, nz, node_column, nboun, npoin_cg

    use mod_initial, only: q_ref, nvar, nvart, kvector, moist_coe, bathymetry, q_ref_layers

    use mod_input, only: si_dimension, eqn_set, delta, space_method, is_shallow, lsalinity, &
         lincompressible, locean, is_swe_layers, nlayers, is_mlswe

    use mod_parallel, only: num_send_recv_total
    
    use mod_basis, only: ngl, nq

    public :: &
        recv_data, &
        grad_qv_ref, grad_qc_ref, grad_qr_ref, q_send, q_recv, &
        f0_ref, g0_ref, grad_g0_ref, d_g0_ref_dr, recv_data_dg, send_data_dg, nmessage, &
        grad_bathy, grad_hB_ref, grad_phiA_ref, grad_rho_ref_layers, &
        q_recv_quad, q_send_quad, recv_data_dg_quad, send_data_dg_quad

    private
    !-----------------------------------------------------------------------

    !module variables and parameters
    real,    dimension(:,:), allocatable :: qb, grad_rho_ref, grad_theta_ref, grad_salinity_ref, grad_press_ref, grad_g0_ref
    real,    dimension(:,:), allocatable :: grad_qv_ref, grad_qc_ref, grad_qr_ref, grad_bathy
    real,    dimension(:,:,:), allocatable :: grad_hB_ref, grad_phiA_ref, grad_rho_ref_layers
    real,    dimension(:),   allocatable :: norm_inf_br, norm_inf_bu, norm_inf_bv, norm_inf_bw, norm_inf_bt, norm_inf_bene
    real,    dimension(:,:), allocatable :: norm_inf_brhs
    real,    dimension(:,:), allocatable :: d_rho_ref_dr, d_theta_ref_dr, d_g0_ref_dr
    real,    dimension(:),   allocatable :: press_ref, g0, f0, h0, press, div_u_ref, f0_ref, g0_ref, dens_var
    real,    dimension(:,:), allocatable :: recv_data
    real,    dimension(:,:,:,:), allocatable :: q_recv, q_send
    real,    dimension(:), allocatable :: recv_data_dg, send_data_dg
    real,    dimension(:,:,:), allocatable :: q_recv_quad, q_send_quad
    real,    dimension(:), allocatable :: recv_data_dg_quad, send_data_dg_quad
    integer :: nmessage
!-----------------------------------------------------------------------

contains

    !-----------------------------------------------------------------------
    subroutine mod_ref_create()

        implicit none
        integer i, j, ix, iz, ip, k
        integer AllocateStatus

        real, dimension(:,:), allocatable :: f, dfdz, q_tempv
        real, dimension(:), allocatable :: q_temp
        real dpds(3), hmatrix(3,3), dpdr, dpdt, max_dp1, max_dp2

        !Size of Max DG Message
        nmessage=2*nvar+4 !nvar+3 for inviscid dynamics and nvar for SIPG gradient + 1 for SIPG constant
        if(is_mlswe) nmessage = 6 ! 2 for uv-momentum and 2 for SIPG gradient + 1  for SIPG constant

        if(allocated(qb)) then
            deallocate(qb, press_ref, press, dens_var, recv_data, norm_inf_brhs, norm_inf_bene, &
                norm_inf_br, norm_inf_bu, norm_inf_bv, norm_inf_bw, norm_inf_bt, &
                grad_press_ref, grad_rho_ref, grad_theta_ref, grad_salinity_ref,div_u_ref, &
                q_send, q_recv, grad_bathy)
        endif
        if(allocated(q_recv_quad)) then
            deallocate(q_recv_quad, q_send_quad)
        end if
        allocate( qb(nvar,npoin), press_ref(npoin),  press(npoin), dens_var(npoin), recv_data(nvar,num_send_recv_total), &
            norm_inf_brhs(nvar,nelem), norm_inf_bene(nelem), &
            norm_inf_br(nelem), norm_inf_bu(nelem), norm_inf_bv(nelem), norm_inf_bw(nelem), norm_inf_bt(nelem), &
            grad_press_ref(3,npoin), grad_rho_ref(3,npoin), grad_theta_ref(3,npoin), grad_salinity_ref(3,npoin), &
            div_u_ref(npoin), q_send(nmessage,ngl,ngl,nboun),q_recv(nmessage,ngl,ngl,nboun), grad_bathy(3,npoin), &
            q_recv_quad(3,nq,nboun), q_send_quad(3,nq,nboun),&
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 0**"

        if (space_method(1:2) == 'dg') then
            if(allocated(recv_data_dg)) then
                deallocate(recv_data_dg, send_data_dg)
            endif

            if(allocated(recv_data_dg_quad)) then
                deallocate(recv_data_dg_quad, send_data_dg_quad)
            endif

            allocate( recv_data_dg(nmessage*ngl*ngl*nboun), &
                send_data_dg(nmessage*ngl*ngl*nboun), &
                recv_data_dg_quad(3*nq*nboun), &
                send_data_dg_quad(3*nq*nboun), &
                stat=AllocateStatus )
            if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ref 1**"
        end if

        !Initialize allocated arrays:
        qb             = 0.0
        press_ref      = p00
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
