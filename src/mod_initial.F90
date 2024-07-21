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
        geometry_type, icase, nlayers, dt, dt_btp, is_mlswe, kstages
    
    use mod_initial_mlswe, only: bot_topo_derivatives, &
        wind_stress_coriolis, compute_reference_edge_variables, Tensor_product, ssprk_coefficients
    
    use mod_Tensorproduct, only: compute_gradient_quad

    public :: &
        mod_initial_create, &
        create_kvector, &
        rhs_continuous, &
        q_init, &
        q_layers_init, &
        rho_layers, hB_grad, phiA_grad, hA, &
        q_exact, &
        q_ref,&
        q_ref_layers, &
        coriolis_constant, kvector, shear_stress, bathymetry, &
        nvar, nvar_diag, nvart, ntracers, &
        nrhs_mxm, &
        height, &
        pi_values,&
        mod_initial_create_height, &
        q_mlswe_init, qprime_mlswe_init, q_df_mlswe_init, pbprime, pbprime_df, q_mlswe_face_init, qprime_face_mlswe_init, pbprime_face, one_over_pbprime, &  ! added by Yao Gahounzo
        one_over_pbprime_face, pbprime_edge, one_over_pbprime_edge, dpprime_df_init, one_over_pbprime_df, &  ! added by Yao Gahounzo
        qb_mlswe_init, qb_face_mlswe_init, qb_df_mlswe_init, alpha_mlswe, layer_dz_eq, tau_wind, coriolis_quad, coriolis_df, & ! added by Yao Gahounzo
        coeff_pbpert_L,coeff_pbpert_R,coeff_pbub_LR, &
        coeff_mass_pbub_L,coeff_mass_pbub_R,coeff_mass_pbpert_LR, N_btp, zbot,zbot_face,zbot_df, grad_zbot_quad, &
        psih, dpsidx,dpsidy, indexq, wjac, fdt_btp, fdt2_btp, a_btp, b_btp, fdt_bcl, fdt2_bcl, a_bcl, b_bcl, a_bclp, b_bclp, qprime_df_init, one_over_pbprime_df_face, tau_wind_df, &
        ssprk_a, ssprk_beta, wjac_df,psih_df,dpsidx_df,dpsidy_df,index_df

    private

    !-----------------------------------------------------------------------
    real, dimension(:,:), allocatable :: q_init, q_exact, q_ref, kvector, q_sph, coord_sph, pi_values, shear_stress, hA
    real, dimension(:,:,:), allocatable :: q_layers_init, hB_grad, phiA_grad, q_ref_layers
    real, dimension(:), allocatable :: rho_layers, bathymetry
    real, dimension(:), allocatable:: height, coriolis_constant
    real, dimension(:,:,:), allocatable :: q_mlswe_init, q_df_mlswe_init, qprime_df_init
    real, dimension(:,:,:,:,:), allocatable :: q_mlswe_face_init
    real, dimension(:,:,:,:,:), allocatable :: qprime_face_mlswe_init
    real, dimension(:,:,:), allocatable :: qprime_mlswe_init
    real, dimension(:), allocatable :: pbprime, pbprime_df, one_over_pbprime_df,one_over_pbprime, layer_dz_eq, alpha_mlswe
    real, dimension(:,:,:), allocatable :: pbprime_face, one_over_pbprime_face, one_over_pbprime_df_face
    real, dimension(:,:), allocatable :: pbprime_edge, one_over_pbprime_edge
    real, dimension(:,:,:,:), allocatable :: qb_face_mlswe_init
    real, dimension(:,:), allocatable :: qb_mlswe_init, qb_df_mlswe_init, tau_wind, dpprime_df_init, tau_wind_df
    real, dimension(:), allocatable :: coriolis_df,coriolis_quad
    real, dimension(:,:), allocatable :: coeff_pbpert_L,coeff_pbpert_R,coeff_pbub_LR, coeff_mass_pbub_L, coeff_mass_pbub_R, coeff_mass_pbpert_LR
    real, dimension(:), allocatable :: zbot, zbot_df, fdt_btp, fdt2_btp, a_btp, b_btp, fdt_bcl, fdt2_bcl, a_bcl, b_bcl, a_bclp, b_bclp
    real, dimension(:,:,:), allocatable :: zbot_face
    real, dimension(:,:), allocatable :: grad_zbot_quad

    real, dimension(:,:), allocatable :: psih, dpsidx,dpsidy, ssprk_a, psih_df,dpsidx_df,dpsidy_df
    integer, dimension(:,:), allocatable :: indexq, index_df
    real, dimension(:), allocatable :: wjac, ssprk_beta, wjac_df

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

        if(allocated(q_init)) deallocate(q_init,q_exact,q_ref,kvector,pi_values,height,coriolis_constant, shear_stress)
        allocate( q_init(nvar,npoin),  &
            q_exact(nvar,npoin), &
            q_ref(nvar,npoin), &
            kvector(3,npoin), &
            pi_values(3,npoin), &
            height(npoin), &
            coriolis_constant(npoin), &
            shear_stress(3,npoin))

        !hack to allocate layers stuff anyways
        allocate(rho_layers(1),bathymetry(npoin))

        if(is_mlswe) then
            if(allocated(q_mlswe_init)) deallocate(q_mlswe_init, qprime_mlswe_init, q_df_mlswe_init, pbprime, pbprime_df, q_mlswe_face_init, qprime_face_mlswe_init, pbprime_face, one_over_pbprime, & 
            one_over_pbprime_face, pbprime_edge, one_over_pbprime_edge, dpprime_df_init, one_over_pbprime_df, & 
            qb_mlswe_init, qb_face_mlswe_init, qb_df_mlswe_init, alpha_mlswe, layer_dz_eq, tau_wind, coriolis_quad, coriolis_df, coeff_pbpert_L, coeff_pbpert_R, coeff_pbub_LR, &
            coeff_mass_pbub_L, coeff_mass_pbub_R, coeff_mass_pbpert_LR, zbot, zbot_df, zbot_face, grad_zbot_quad, qprime_df_init, one_over_pbprime_df_face, tau_wind_df,&
            ssprk_a,ssprk_beta)
            allocate(q_mlswe_init(3,npoin_q,nlayers), qprime_mlswe_init(3,npoin_q,nlayers), q_df_mlswe_init(3,npoin,nlayers), pbprime(npoin_q), pbprime_df(npoin), &
            q_mlswe_face_init(3,2,nq,nface,nlayers), qprime_face_mlswe_init(3,2,nq,nface,nlayers), pbprime_face(2,nq,nface), one_over_pbprime(npoin_q), &
            one_over_pbprime_face(2,nq,nface), pbprime_edge(nq,nface), one_over_pbprime_edge(nq,nface), dpprime_df_init(npoin,nlayers), &
            one_over_pbprime_df(npoin), qb_mlswe_init(4,npoin_q), qb_face_mlswe_init(4,2,nq,nface), qb_df_mlswe_init(4,npoin), &
            alpha_mlswe(nlayers), layer_dz_eq(nlayers), tau_wind(2,npoin_q), coriolis_quad(npoin_q), coriolis_df(npoin), &
            coeff_mass_pbpert_LR(nq,nface), coeff_pbpert_L(nq,nface),coeff_pbpert_R(nq,nface),coeff_pbub_LR(nq,nface), &
            coeff_mass_pbub_L(nq,nface),coeff_mass_pbub_R(nq,nface), &
            zbot(npoin_q), zbot_df(npoin), zbot_face(2,nq,nface), grad_zbot_quad(2,npoin_q), &
            psih(npoin_q,npts), dpsidx(npoin_q,npts), dpsidy(npoin_q,npts), indexq(npoin_q,npts), wjac(npoin_q), &
            fdt_btp(npoin), fdt2_btp(npoin), a_btp(npoin), b_btp(npoin), fdt_bcl(npoin), fdt2_bcl(npoin), a_bcl(npoin), &
            b_bcl(npoin), a_bclp(npoin), b_bclp(npoin), qprime_df_init(3,npoin,nlayers), one_over_pbprime_df_face(2,ngl,nface), tau_wind_df(2,npoin), &
            ssprk_a(kstages,3), ssprk_beta(kstages), wjac_df(npoin),psih_df(npoin,npts),dpsidx_df(npoin,npts),dpsidy_df(npoin,npts),index_df(npoin,npts))

            q_mlswe_init = 0.0
            qprime_mlswe_init = 0.0
            q_df_mlswe_init = 0.0
            pbprime = 0.0
            pbprime_df = 0.0
            q_mlswe_face_init = 0.
            qprime_face_mlswe_init = 0.0

        end if

        !Initialize q_* arrays to zero:
        q_init  = 0.0
        q_exact = 0.0
        q_ref   = 0.0

        time = 0.0
    
        if(is_mlswe) then

            call Tensor_product(wjac,psih,dpsidx,dpsidy,indexq, wjac_df,psih_df,dpsidx_df,dpsidy_df,index_df)

            call initial_conditions_mlswe(q_mlswe_init, qprime_mlswe_init, q_df_mlswe_init, pbprime, pbprime_df, q_mlswe_face_init, &
                qprime_face_mlswe_init, pbprime_face, one_over_pbprime, one_over_pbprime_face, pbprime_edge, one_over_pbprime_edge, &
                dpprime_df_init, one_over_pbprime_df, layer_dz_eq, qb_mlswe_init, qb_face_mlswe_init, qb_df_mlswe_init, qprime_df_init, &
                alpha_mlswe, one_over_pbprime_df_face,zbot_df, coord)

            call compute_reference_edge_variables(coeff_pbpert_L,coeff_pbpert_R,coeff_pbub_LR,coeff_mass_pbub_L, &
                coeff_mass_pbub_R,coeff_mass_pbpert_LR, pbprime_face,alpha_mlswe)

            call bot_topo_derivatives(zbot,zbot_face,zbot_df)

            call compute_gradient_quad(grad_zbot_quad,zbot_df)

            N_btp = ceiling(dt/dt_btp)
            dt_btp = dt/real(N_btp)

            call wind_stress_coriolis(tau_wind,coriolis_df,coriolis_quad,fdt_btp, fdt2_btp, a_btp, b_btp, fdt_bcl, fdt2_bcl, a_bcl, b_bcl,coord,icase, b_bclp, a_bclp, tau_wind_df)

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
    !----------------------------------------------------------------------!
    !>@brief This subroutine constructs height
    !----------------------------------------------------------------------!
    subroutine mod_initial_create_height(height,kvector,coord,npoin)

        implicit none

        !global arrays
        real height(npoin), kvector(3,npoin), coord(3,npoin)
        integer npoin

        !local
        integer ip
        real x, y, z, xf, yf, zf

        do ip = 1,npoin
            x=coord(1,ip); y=coord(2,ip); z=coord(3,ip)
            xf=kvector(1,ip); yf=kvector(2,ip); zf=kvector(3,ip);
            height(ip) = x*xf + y*yf + z*zf
        end do

        if (geometry_type(1:6) == 'sphere') then
            height(:) = height(:) - earth_radius
        end if

    end subroutine mod_initial_create_height

end module mod_initial
