!----------------------------------------------------------------------!
!>@brief This module contains Initial Conditions 
!>@author J.F. Kelly and F.X. Giraldo
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943-5216
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: April 03, 2023
!----------------------------------------------------------------------!
module mod_initial

    use mod_constants, only: gravity, earth_radius, omega

    use mod_grid, only:  npoin, coord, npoin_cg, nface, npoin_q

    use mod_basis, only: nq, npts, ngl

    use mod_input, only: time_initial, time_final, time_restart, time_scale, &
        nlayers, dt, dt_btp, is_mlswe, kstages
    
    use mod_initial_mlswe, only: bot_topo_derivatives, &
        wind_stress_coriolis, compute_reference_edge_variables, ssprk_coefficients, &
        compute_reference_edge_variables_df !, Tensor_product
    
    use mod_Tensorproduct, only: compute_gradient_quad, compute_gradient_df

    public :: &
        mod_initial_create, &
        create_kvector, &
        q_init, &
        hB_grad, phiA_grad, hA, &
        q_ref,&
        q_ref_layers, &
        coriolis_constant, kvector, shear_stress, bathymetry, &
        nvar, nvar_diag, nvart, ntracers, &
        nrhs_mxm, &
        height, &
        pi_values 

    public :: q_df_mlswe_init, pbprime, pbprime_df, &
        pbprime_face, one_over_pbprime, &
        one_over_pbprime_face, pbprime_edge, one_over_pbprime_edge, one_over_pbprime_df, & 
        qb_df_mlswe_init, alpha_mlswe, tau_wind, coriolis_quad, coriolis_df, & 
        coeff_pbpert_L,coeff_pbpert_R,coeff_pbub_LR, &
        coeff_mass_pbub_L,coeff_mass_pbub_R,coeff_mass_pbpert_LR, N_btp, zbot,zbot_face,zbot_df, &
        grad_zbot_quad, psih, dpsidx,dpsidy, indexq, wjac, fdt_bcl, fdt2_bcl, a_bcl, b_bcl, &
        qprime_df_init, one_over_pbprime_df_face, tau_wind_df, &
        ssprk_a, ssprk_beta, wjac_df,psih_df,dpsidx_df,dpsidy_df,index_df, &
        grad_zbot_df, pbprime_df_face

    public :: z_interface

    private

    !-----------------------------------------------------------------------
    real, dimension(:,:), allocatable :: q_init, q_exact, q_ref, kvector, q_sph, coord_sph
    real, dimension(:,:), allocatable :: pi_values, shear_stress, hA
    real, dimension(:,:,:), allocatable :: hB_grad, phiA_grad, q_ref_layers
    real, dimension(:), allocatable :: rho_layers, bathymetry
    real, dimension(:), allocatable:: height, coriolis_constant
    real, dimension(:,:,:), allocatable :: q_df_mlswe_init, qprime_df_init
    real, dimension(:), allocatable :: pbprime, pbprime_df, one_over_pbprime_df,one_over_pbprime 
    real, dimension(:,:,:), allocatable :: pbprime_face, one_over_pbprime_face
    real, dimension(:,:,:), allocatable :: one_over_pbprime_df_face, pbprime_df_face
    real, dimension(:,:), allocatable :: pbprime_edge, one_over_pbprime_edge
    real, dimension(:,:,:,:), allocatable :: qb_face_mlswe_init
    real, dimension(:,:), allocatable :: qb_df_mlswe_init, tau_wind, tau_wind_df
    real, dimension(:), allocatable :: coriolis_df,coriolis_quad
    real, dimension(:,:), allocatable :: coeff_pbpert_L,coeff_pbpert_R,coeff_pbub_LR
    real, dimension(:,:), allocatable :: coeff_mass_pbub_L, coeff_mass_pbub_R, coeff_mass_pbpert_LR
    real, dimension(:), allocatable :: zbot, zbot_df, fdt_bcl, fdt2_bcl, a_bcl, b_bcl
    real, dimension(:,:,:), allocatable :: zbot_face
    real, dimension(:,:), allocatable :: grad_zbot_quad, grad_zbot_df, z_interface

    real, dimension(:,:), allocatable :: psih, dpsidx,dpsidy, ssprk_a, psih_df,dpsidx_df,dpsidy_df
    integer, dimension(:,:), allocatable :: indexq, index_df
    real, dimension(:), allocatable :: wjac, ssprk_beta, wjac_df, alpha_mlswe

    integer :: nvar, nvart, nvar_diag
    integer :: nrhs_mxm, N_btp
  !-----------------------------------------------------------------------

    contains

    !-----------------------------------------------------------------------
    subroutine mod_initial_create()

        implicit none

        integer i, iperturbation, ip
        real time
        real x, y, z, xf, yf, zf, radius
        real lat, lon, press, temp, phis, ps
        real pb, tb, pi_f, pi_b, zradius

        !Define the number of prognostic variables
        nvar=5 !rho,u,v,w,theta
        nvart=nvar
        !Define the number of diagnostic variables
        nvar_diag=0

        !Store Number of RHS for MXM calls in CREATE_RHS_VOLUME
        nrhs_mxm=nvar + 1 !1 is for Pressure

        !Store Number of Tracers
        ntracers=nvar-5

        if(allocated(q_init)) deallocate(q_init,q_exact,q_ref,kvector,pi_values,height, &
                coriolis_constant, shear_stress)
        allocate( q_init(nvar,npoin), q_exact(nvar,npoin), q_ref(nvar,npoin), kvector(3,npoin), &
            pi_values(3,npoin), height(npoin), coriolis_constant(npoin), shear_stress(3,npoin))

        !hack to allocate layers stuff anyways
        allocate(rho_layers(1),bathymetry(npoin))

        if(is_mlswe) then
            if(allocated(q_df_mlswe_init)) deallocate(q_df_mlswe_init, pbprime, &
            pbprime_df, pbprime_face, one_over_pbprime, & 
            one_over_pbprime_face, pbprime_edge, one_over_pbprime_edge, one_over_pbprime_df, & 
            qb_df_mlswe_init, alpha_mlswe, tau_wind, coriolis_quad, &
            coriolis_df, coeff_pbpert_L, coeff_pbpert_R, coeff_pbub_LR, &
            coeff_mass_pbub_L, coeff_mass_pbub_R, coeff_mass_pbpert_LR, zbot, zbot_df, zbot_face, &
            grad_zbot_quad, qprime_df_init, one_over_pbprime_df_face, tau_wind_df,&
            ssprk_a,ssprk_beta, grad_zbot_df, &
            pbprime_df_face, z_interface)
            allocate(q_df_mlswe_init(3,npoin,nlayers), pbprime(npoin_q), pbprime_df(npoin), &
            pbprime_face(2,nq,nface), one_over_pbprime(npoin_q), &
            one_over_pbprime_face(2,nq,nface), pbprime_edge(nq,nface), &
            one_over_pbprime_edge(nq,nface), &
            one_over_pbprime_df(npoin), qb_df_mlswe_init(4,npoin), &
            alpha_mlswe(nlayers), tau_wind(2,npoin_q), coriolis_quad(npoin_q), coriolis_df(npoin), &
            coeff_mass_pbpert_LR(nq,nface), coeff_pbpert_L(nq,nface),coeff_pbpert_R(nq,nface), &
            coeff_pbub_LR(nq,nface), coeff_mass_pbub_L(nq,nface),coeff_mass_pbub_R(nq,nface), &
            zbot(npoin_q), zbot_df(npoin), zbot_face(2,nq,nface), grad_zbot_quad(2,npoin_q), &
            grad_zbot_df(2,npoin), psih(npts,npoin_q), dpsidx(npts,npoin_q), dpsidy(npts,npoin_q), &
            indexq(npts,npoin_q), wjac(npoin_q), fdt_bcl(npoin), fdt2_bcl(npoin), a_bcl(npoin), &
            b_bcl(npoin), qprime_df_init(3,npoin,nlayers), one_over_pbprime_df_face(2,ngl,nface), &
            tau_wind_df(2,npoin), ssprk_a(kstages,3), ssprk_beta(kstages), wjac_df(npoin), &
            psih_df(npts,npoin), dpsidx_df(npts,npoin),dpsidy_df(npts,npoin),index_df(npts,npoin), &
            pbprime_df_face(2,ngl,nface),z_interface(npoin,nlayers+1))

            q_df_mlswe_init = 0.0
            pbprime = 0.0
            pbprime_df = 0.0

        end if

        !Initialize q_* arrays to zero:
        q_init  = 0.0
        q_exact = 0.0
        q_ref   = 0.0

        time = 0.0
    
        if(is_mlswe) then

            call Tensor_product(wjac,psih,dpsidx,dpsidy,indexq, wjac_df,psih_df,dpsidx_df, &
                    dpsidy_df,index_df)

            call initial_conditions(q_df_mlswe_init, pbprime, pbprime_df, &
                pbprime_face, one_over_pbprime, one_over_pbprime_face, pbprime_edge, &
                one_over_pbprime_edge, one_over_pbprime_df, qb_df_mlswe_init, qprime_df_init, &
                alpha_mlswe, one_over_pbprime_df_face, pbprime_df_face, zbot_df,tau_wind_df, &
                z_interface)

            call compute_reference_edge_variables(coeff_pbpert_L,coeff_pbpert_R,coeff_pbub_LR, &
                coeff_mass_pbub_L, coeff_mass_pbub_R,coeff_mass_pbpert_LR, pbprime_face,alpha_mlswe)
                
            call bot_topo_derivatives(zbot,zbot_face,zbot_df)

            call compute_gradient_quad(grad_zbot_quad,zbot_df)
            call compute_gradient_df(grad_zbot_df,zbot_df)

            N_btp = ceiling(dt/dt_btp)
            dt_btp = dt/real(N_btp)

            call wind_stress_coriolis(tau_wind,coriolis_df,coriolis_quad, fdt_bcl, fdt2_bcl, &
                a_bcl, b_bcl, tau_wind_df)

            call ssprk_coefficients(ssprk_a,ssprk_beta)
        endif

        !Set-up Times
        time_initial=time_initial*time_scale
        time_final=time_final*time_scale
        time_restart=time_restart*time_scale

    end subroutine mod_initial_create

    !----------------------------------------------------------------------!
    !>@brief This subroutine constructs the KVECTOR on the sphere
    !> (called the Rvector in Notes and Papers)
    !----------------------------------------------------------------------!
    subroutine create_kvector(kvector,coord,npoin)

        implicit none

        !global arrays
        real kvector(3,npoin), coord(3,npoin)
        integer npoin

        !local
        integer ip
        real x, y, z, radius, xf, yf, zf

        do ip = 1,npoin
            x=coord(1,ip); y=coord(2,ip); z=coord(3,ip)
            radius=sqrt( dot_product(coord(:,ip),coord(:,ip)) )
            xf=x/radius; yf=y/radius; zf=z/radius
            kvector(1,ip) = xf
            kvector(2,ip) = yf
            kvector(3,ip) = zf
        end do

    end subroutine create_kvector

end module mod_initial
