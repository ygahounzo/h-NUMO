!----------------------------------------------------------------------!
!>@brief This module defines the Constants that were in PARAM.H
!>@author F.X. Giraldo and J.F. Kelly
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943-5216
!>
!> @date November 10 2014, S. Marras added sponge_lateralx_coe_east and sponge_lateralx_coe_west
!> @date November 10 2014, S. Marras added lforce_spongex, lforce_spongey, lforce_spongez
!>                         If not selected in input by the user, everything will behave as before, 
!>                         otherwise, the flag lsponge will become true even if the boundary
!>                         codes iboundary(1:6) are different from 6 in input. These are used in mod_bc.f90
!> @date December 8 2014,  S. Marras added dirichlet boundary conditions values
!> @date December 17 2014, S. Marras added a check on metis3D and periodic b.c.
!> @date March       2015, S. Marras added a lvisc_dynamics, lvisc_tracers, lnazarov_tracers
!> @date June 14, 2015, F.X. Giraldo added checks to the input data - space_method and explicit and si_dimension compatibilit
!> @date April 4, 2017, F.X. Giraldo added Incompressible NSE flags for LVISC_RECIPROCAL and others (forgot to update this line but PISO_TYPE was also added).
!----------------------------------------------------------------------!
module mod_input

   !    use mod_utilities, only : get_unit
   use mod_types, only : r8
 
   use mod_utilities, only : get_unit, lowercase, uppercase
 
   use iso_c_binding, only: C_INT32_T
 
   implicit none
 
   public :: namelist_input, &
        mod_input_create, &
        dt, dt0, dt1, dt2, &
        restoring_time, &
        lrestoring_sponge, &
        lpi_fraction, &
        time_initial, time_final, time_dynamic_amr, time_restart, time_scale, irestart_file_number,  &
        measurement_interval, &
        interval_start, &
        interval_end,   &
        dt_safety_factor, &
        lgravity, &
        lrestart_file, &
        icase, ti_method, ark2_type, si_method, si_dimension, &
        ti_method_btp, & ! Added by YaoG
        case9015_inlet_vel,&
        case9101_initial_time,&
        Iter_Type, Proj_Size, Focal_Dist_Adj, precon_order, &
        precon_type, precon_mode, PBNO_pNorm, precon_fname, vertical_precon, &
        Imag_Wt, Hmat_Size, &
        c_am2, lf2_filter, lf2_raw_filter, lf2_alpha, kstages, solver_type, lcg_proj, lf2_nu, &
        solver_tol, solver_tol_pressure, solver_tol_ALE, jfnk_tol, epsilon_jfnk, gamma1, delta, filter_mux, filter_muy, filter_muz, &
        ifilter, limit_solver_iter_press, limit_solver_iter_ALE,&
        filter_weight_type, filter_basis_type, &
        fname_root, out_type, lout_ascii, lout_asciimaya, format_vtk, nvtk_files, vtk_cell_type, &
        write_mesh, &
        visc, visc2, visc4, &
        visc_h, visc_v, &
        viscx, &
        viscy, &
        viscz, &
        diff_T, diff_Th, diff_Tv, &
        diff_S, diff_Sh, diff_Sv, &
        lsalinity, &
        which_eos, &
        locean,&
        lbalance_pressure, &
        SIPG_constant, &
        lvisc_reciprocal, & !FXG 4/5/17: New flag to divide by VISC to get proper convergence rates.
        lvisc_anisotropic, & !FXG 10/31/16: New flag to split Hor and Vert dirs.
        lvisc_dynamics, visc_dyn_flg,     &
        lvisc_tracers, visc_tracers_flg, &
        lnazarov_tracers, naza_tracers_flg, &
        vertical_viscosity, &
        alpha_thermal, alpha_salinity,&
        filter_tracers_flg, &
        lVMS, lshock_cpt, &
        lVMS_dyn, lshock_cpt_dyn, &
        tau_type, shock_type, &
        Pr, C1_in, C2_in, Cs_in, &
        Re_tau, &
        Re_bulk, &
        Mach, &
        lnazarov, &
        lapply_bcs_velocity, &
        lLAV, &
        llimit, llimit_below, llimit_above, limiter_type, limiter_qoi, &
        ltime_residual, time_residual_flg, &
        lLES, &
        lSMAG, &
        ladd_full_stress_tensor, &
        full_stress_flg, &
        ldamping_function, &
        damping_function_flg, &
        lforce_turb_momentum, &
        momentum_forcing_flg, &
        momentum_forcing_LES_flg, &
        luse_PSUP, &
        luse_min_element_length, &
        min_length_flg, max_length_flg, &
        moist_forcing_flg, moist_nocolumn_flg, rain_flg, buoyancy_flg, &
        lstabilize_mass_eqn, &
        mass_stabilize_flg, &
        nlaplacian, &
        nlaplacian_flg, &
        lanisotropic_laplacian, &
        ladapt_timestep, lprint_diagnostics, iprint_diagnostics, lcheck_all_conserved, lunit_test, unit_test_type, lcompute_barycenter, ljenkins_test, &
        ldynamics, lphysics, lsimple_physics, lkessler, lrain, lmoist_forcing, lmoist_column, lpassive, lpure_advection, ntracers_in, &
        kessler_production_constant, &
        sounding_file, lread_sound, lread_qvapor, lread_uvelo, lread_vvelo, &
        ctd_file_temp, ctd_file_salinity, &
        sounding_ncolumns, &
        uvelo_uniform, vvelo_uniform, wvelo_uniform, &
        lrotating_flow, &
        Mach1_in, Mach2_in, &
        bcast_type, &
        lout_rho, lout_uvelo, lout_vvelo, lout_wvelo, lout_velo, lout_dvelo, lout_theta, lout_sponge, &
        lout_thetae , lout_press , lout_temperature, &
        lout_previous_time_steps, &
        lout_rank, &
        lout_vorticity, &
        lout_xvorticity, &
        lout_yvorticity, &
        lout_zvorticity,&
        lout_rvorticity,&
        lout_tke, &
        lout_qvapor , lout_qcloud , lout_qrain, &
        lout_tau, lout_radius, lout_coords, &
        lout_spherical, &
        lwrite_time_averaging, &
        lout_statistics, &
        pert_type, &
        pert_type_tracers, &
        lleveque, &
        lburgers1d, &
        lburgers2d, &
        lsmolar_flow,&
        lslotted, &
        lout_spherical_shell, &
        zlevel_out, &
        lout_nc_3d, &
        p00_in, &
        nu_air_in, &
        lpert_theta, lpert_qv, lpert_qc, lpert_qr, lcylinder, &
        ldam, &
        lpert_velo, &
        ljet, &
        thetac_in, qvc_in, qcc_in, qrc_in, &
        thetac2_in, qvc2_in, qcc2_in, qrc2_in, &
        thetac3_in, qvc3_in, qcc3_in, qrc3_in, &
        thetac4_in, qvc4_in, qcc4_in, qrc4_in, &
        steepness_pert_in, &
        xradius_pert_in, yradius_pert_in, zradius_pert_in, &
        xradius2_pert_in, yradius2_pert_in, zradius2_pert_in, &
        xradius3_pert_in, yradius3_pert_in, zradius3_pert_in, &
        xradius4_pert_in, yradius4_pert_in, zradius4_pert_in, &
        xc_pert_in, yc_pert_in, zc_pert_in, &
        xc2_pert_in, yc2_pert_in, zc2_pert_in, &
        xc3_pert_in, yc3_pert_in, zc3_pert_in, &
        xc4_pert_in, yc4_pert_in, zc4_pert_in, &
        xmin_velo_in, xmax_velo_in, &
        ymin_velo_in, ymax_velo_in, &
        zmin_velo_in, zmax_velo_in, &
        space_method, &
        cgdg_method, &
        form_method, &
        dump_rhs, &
        n_corrections, &
        imass, &
        ad_mlswe, explt_coriolis, cd_mlswe, dp_cutoff1, dp_cutoff2, &
        dp_tau_bot, dp_tau_wind, dt_btp, method_visc,visc_mlswe,dpprime_visc_min, max_shear_dz, matlab_viz, &
        adjust_H_vertical_sum, is_mlswe_linear, botfr, mass_exact, bcl_flux, mlswe_bc_strong, dg_integ_exact, dump_data
 
   public :: nsmall, ti_alpha, ti_beta, dd_alpha
 
   public :: eqn_set, equations, lincompressible, piso_type, lpassive_advection, lhelmholtz, lALE, &
        lset_outlet_pressure
 
   public :: is_shallow, is_sphere, is_implicit, is_1d, is_imex, is_swe_layers, nonlinear_swe, is_mlswe
 
   public :: nelx, nely, nelz, nopx, nopy, nopz, xdims, ydims, ztop, zbottom, &
        nlayers, & !shallow water layers
        lread_input_dimensions, &
        delta_domain, & !used for turbulent channel flows
        nlon, nlat, &
        x_boundary, y_boundary, z_boundary, &
        x_periodic, y_periodic, z_periodic, &
        robin_bc_alpha, &
        lfree_slip_exact, &
        geometry_type, decomp_type, nproc_z, &
        LAM_flg, GCM_flg, &
        bc_tscale, bc_xscale, bc_yscale, bc_zscale, &
        ldirichlet_all, &
        ldirichlet_top, ldirichlet_bottom,  &
        ldirichlet_east, ldirichlet_west,   &
        ldirichlet_north, ldirichlet_south, &
        dirichlet_theta, &
        dirichlet_uvelo, &
        dirichlet_vvelo, &
        dirichlet_wvelo, &
        lviscous_boundary, &
        lviscous_ground, &
        lviscous_bottom, &
        lviscous_top, &
        lviscous_south, &
        lviscous_north, &
        lviscous_east, &
        lviscous_west, &
        lforce_spongex, &
        lforce_spongey, &
        lforce_spongez, &
        sponge_type, sponge_top_coe, sponge_lateralx_coe, sponge_lateralx_coe_east, sponge_lateralx_coe_west, sponge_lateraly_coe, &
        lsommerfeld, &
        lgrid_only, &
        lp4est, lp6est, lserial_grid_creation, lparallel_grid_creation, lio_grid_ascii, lwrite_grid_ascii, fname_initial, &
        restart_path, & ! Added by Yao G.
        nel_root_h, refinement_levels_h, nel_root_v, refinement_levels_v, &
        xstretch_coe, ystretch_coe, zstretch_coe, &
        lxstretch, lystretch, lzstretch, &
        xstretch_type, ystretch_type, zstretch_type, &
        lmountain,mount_type,mount_xc,mount_yc,mount_hm,mount_ac,mount_bc,&
        zgrid_type,hybrid_s,sleve_s1,sleve_s2, &
        lread_external_grid, read_external_grid_flg, is_non_conforming_flg, &
        interp_cg_flux_flg, &
        max_mesh_lvl, &
        nc_box_invert, &
        p4est_log_level, &
        xlim_min, xlim_max, ylim_min, ylim_max, zlim_min, zlim_max, &
        amr_indicator_variables, &
        amr_smoothness_limits, &
        amr_max_min_lim, &
        amr_threshold_lim, &
        amr_smoothness_qL2_limit, &
        amr_mark_max_min, &
        amr_mark_random, &
        amr_mark_threshold, &
        amr_mark_modes, &
        amr_mark_modes_use_baseline_decay, &
        amr_num_neigh_iter, &
        amr_mark_set2nc, &
        lread_bc,lgpu, numaocca_dir, threads_per_process, cpus_per_node, gpus_per_node, &
        platform, platformID, deviceID, platformWeight, &
        platform2, platformID2, deviceID2, platformWeight2, &
        luse_hybrid_cpu_gpu, Nelems, Nslices, NslicesV, &
        CudaCompilerFlags, OpenCLCompilerFlags, OpenMPCompilerFlags, &
        SerialCompilerFlags, PthreadsCompilerFlags, COICompilerFlags, &
        vectorization
 
   !!-- Shallow water inputs
   public :: ibathymetry, bathymetry_file, wave_file, &
        gravity_in, &
        xc_in, yc_in, rc_in, &
        xc2_in, yc2_in, rc2_in, &
        xc3_in, yc3_in, rc3_in, &
        xc4_in, yc4_in, rc4_in, &
        xc5_in, yc5_in, rc5_in, &
        xc6_in, yc6_in, rc6_in, &
        xc7_in, yc7_in, rc7_in, &
        beach_length_in, mesh_file, &
        xmin_topography_in, h0_in, free_surface_level_in, &
        lobstacle, &
        lobstacle_semicircle, &
        lobstacle_ellipsis, &
        lobstacle_distance_as_lambda_fraction, &
        lobstacle_distance_as_obstacle_length, &
        lobstacle_height_as_amplitude_fraction, &
        lobstacle_size_as_lambda_fraction, &
        lobstacle_crosssection_as_lambda_fraction, &
        Lobs_in, CXobs_in, Hobs_in, &
        Hobs2_in,  & !obstacle height
        Hobs3_in,  & !obstacle height
        Hobs4_in,  & !obstacle height
        Hobs5_in,  & !obstacle height
        ntimes_Lobs_in, &
        xdims_obstacle, ydims_obstacle, &
        xmin_obstacle, xmax_obstacle,   &
        ymin_obstacle, ymax_obstacle,   &
        x_obstacle_offset_in, &
        slope_in, &
        initial_shore_line_in, &
        dam_depth_in, dam_xlimit_in, &
        lwave, lread_external_bathy, lread_wave_from_file, wave_amplitude_in, &
        Lambda_in, lambda_fraction_in, obstacle_crosssection_lambda_fraction_in, wave_start_in, wave_type, &
        wave_min_crest_in, wave_max_crest_in, &
        carrier_scaling_in, carrier_x_offset_in, &
        wave_xcenter_in, wave_ycenter_in, &
        distance_lambda_fraction_in, &
        amplitude_height_fraction_in, &
        amplitude_in, amplitude2_in, amplitude3_in, &
        amplitude4_in, amplitude5_in, amplitude6_in, amplitude7_in, &
        hump_config, bathymetry_shift, limit_threshold, &
        amplitudes_in, nobstacles_in, radiusx_in, radiusy_in, xcc_in, ycc_in, &
        cms_coefficient, cms_coefficient2
 
              
   private
 
   !-----------------------------------------------------------------------
   ! Namelist Variables
   !-----------------------------------------------------------------------
   real(kind=r8) :: dt, dt0, dt1, dt2, dt_btp
   real(kind=r8) :: restoring_time = 1000000000
   logical       :: lrestoring_sponge = .false.
   real(kind=r8) :: time_initial, time_final, time_restart, lf2_filter, lf2_raw_filter, lf2_alpha, lf2_nu
   real(kind=r8) :: measurement_interval
   real(kind=r8) :: interval_start=0.0
   real(kind=r8) :: time_dynamic_amr=0.0
   real(kind=r8) :: interval_end=0.0
   real(kind=r8) :: dt_meas = 0.0
   real(kind=r8) :: gamma1, delta
   real(kind=r8) :: solver_tol = 1e-2
   real(kind=r8) :: jfnk_tol = 1e-4
   real(kind=r8) :: epsilon_jfnk = 1e-5
   real(kind=r8) :: filter_mux, filter_muy, filter_muz
   real(kind=r8) :: solver_tol_pressure = -1.0
   real(kind=r8) :: solver_tol_ALE = 1e-1
   real(kind=r8) :: case9015_inlet_vel = 1.0
   real(kind=r8) :: case9101_initial_time = 1.0
   integer       :: icase, kstages, ifilter
   integer       :: irestart_file_number = 0
   integer       :: limit_solver_iter_press = 500
   integer       :: limit_solver_iter_ALE = 100
   integer       :: nlaplacian
   integer       :: nlaplacian_flg = 0
   character     :: filter_basis_type*8, filter_weight_type*4
   character     :: fname_root*150, out_type*5, fname_initial*100
   character     :: restart_path*100
   logical       :: write_mesh=.false.
   character     :: ti_method*14, si_method*8, si_dimension*2
   character     :: ti_method_btp*14
   character(len=1) :: ark2_type = 'b' !default set to ARK2B
   character     :: sounding_file*124, ctd_file_temp*124, ctd_file_salinity*124
   logical       :: lprint_diagnostics
   integer       :: iprint_diagnostics = 1 !every how many iterations to print diagnostics
   logical       :: lcheck_all_conserved = .false.
   logical       :: ljenkins_test = .true.
   logical       :: lunit_test = .false.
 
   character     :: solver_type*9, real_string*9, unit_test_type*16
   logical       :: lcg_proj = .false.
   character     :: Iter_Type*5, precon_type*1, precon_mode*5, precon_fname*150
   character     :: space_method*3, cgdg_method*8
   character(len=17) :: form_method='strong' !strong form by default
   integer       :: precon_order, PBNO_pNorm, Hmat_Size, Proj_Size
   real(kind=r8) :: Focal_Dist_Adj, Imag_Wt
   logical       :: vertical_precon = .false.
 
 
   real(kind=r8) :: ad_mlswe = 0.0
   real(kind=r8) :: explt_coriolis = 0.0
   real(kind=r8) :: cd_mlswe = 0.0
   real(kind=r8) :: dp_cutoff1 = 0.0
   real(kind=r8) :: dp_cutoff2 = 0.0
   real(kind=r8) :: dp_tau_bot = 0.0
   real(kind=r8) :: dp_tau_wind = 0.0
   integer :: method_visc = 0
   integer :: botfr = 0
   integer :: adjust_H_vertical_sum = 0
   real(kind=r8) :: visc_mlswe = 0.0
   real(kind=r8) :: dpprime_visc_min = 0.0
   real(kind=r8) :: max_shear_dz = 0.0
   real(kind=r8) :: bcl_flux = 0.0
 
   !-----------------------------------------------------------------------
   ! Namelist Variables
   !-----------------------------------------------------------------------
   real(kind=r8), dimension(2) :: xdims, ydims
   real(kind=r8)     :: ztop
   real(kind=r8)     :: pi_trig
   real(kind=r8)     :: zbottom = 0.0
   real(kind=r8)     :: delta_domain = 1.0 !used for turbulent channel flows (default value is 1.0)
   real(kind=r8)     :: xstretch_coe, ystretch_coe, zstretch_coe
   integer           :: nelx, nely, nelz, nopx, nopy, nopz
   integer           :: nlon = 180 !default values
   integer           :: nlat = 45  !default values
   integer           :: nlayers = 1 !number of shallow water layers
   integer           :: nproc_z
   integer           :: zlevel_out = 1 !This integer is given by the user to decide what spherical level to write to a netcdf file. Value 1 is by default
   character(len=13) :: geometry_type
   character(len=7)  :: decomp_type
   character(len=24) :: sponge_type = 'no_sponge'
   character(len=8)  :: format_vtk  = 'BINARY' !can be ascii, binary
   character(len=8)  :: vtk_cell_type  = 'standard' !can be standard, lagrange
   integer           :: nvtk_files = 1 !number of files written in parallel (must be less than nproc)
 
   character(len=24) :: xstretch_type
   character(len=24) :: ystretch_type
   character(len=24) :: zstretch_type
   character(len=24) :: platform = 'OpenCL'
   character(len=24) :: platform2 = 'OpenMP'
   character(len=6)  :: vectorization = 'float4'
   character(len=150):: numaocca_dir = './numaocca/'
 
   character(len=512) :: OpenCLCompilerFlags = ' -cl-denorms-are-zero -cl-fast-relaxed-math -cl-finite-math-only -cl-mad-enable -cl-no-signed-zeros'
   character(len=512) :: CudaCompilerFlags = '--compiler-options -O2 --use_fast_math'
   character(len=512) :: OpenMPCompilerFlags = '-D__extern_always_inline=inline -O2'
   character(len=512) :: SerialCompilerFlags = '-D__extern_always_inline=inline -O2'
   character(len=512) :: PthreadsCompilerFlags = '-D__extern_always_inline=inline -O2'
   character(len=512) :: COICompilerFlags = '-D__extern_always_inline=inline -O2'
   !Vectorization OpenMP flags ===>>>  "-mavx"
   !Endevor OpenMP flags       ===>>>  "-xCOMMON-AVX512 -qopt-report=5 -qopt-report-file=stdout"
 
  
   !real :: dt_btp
 
   !-----------------------------------------------------------------------
   ! Default Namelist Values: Default Values
   !-----------------------------------------------------------------------
   real(kind=r8)         :: bc_tscale = 240.0
   real(kind=r8)         :: bc_xscale = -1000.0
   real(kind=r8)         :: bc_yscale = -1000.0
   real(kind=r8)         :: bc_zscale = -1000.0
   real(kind=r8)         :: sponge_top_coe = -1
   real(kind=r8)         :: sponge_lateralx_coe = -1
   real(kind=r8)         :: sponge_lateralx_coe_east = -1
   real(kind=r8)         :: sponge_lateralx_coe_west = -1
   real(kind=r8)         :: sponge_lateraly_coe = -1
 
   real(kind=r8)         :: visc = 0.0   !Default value if not given in input
   real(kind=r8)         :: visc2= 0.0   !Default value if not given in input
   real(kind=r8)         :: visc4= 0.0   !Default value if not given in input
   real(kind=r8)         :: viscx = 0.0
   real(kind=r8)         :: viscy = 0.0
   real(kind=r8)         :: viscz = 0.0
   real(kind=r8)         :: visc_h = 0.0 !Horizontal viscosity
   real(kind=r8)         :: visc_v = 0.0 !Vertical viscosity
   real(kind=r8)         :: diff_T = -1.0 !Thermal diffusivity
   real(kind=r8)         :: diff_Th = -1.0 !Horizontal thermal diffusivity
   real(kind=r8)         :: diff_Tv = -1.0 !Vertical thermal diffusivity
   real(kind=r8)         :: diff_S = -1.0 !Salinity diffusivity
   real(kind=r8)         :: diff_Sh = -1.0 !Horizontal salinity diffusivity
   real(kind=r8)         :: diff_Sv = -1.0 !Vertical salinity diffusivity
   logical               :: lsalinity = .false. !salinity flag
   character(len=32)     :: which_eos !which eos to use: 'linear','simplified' - quadratic in temperature
   logical               :: locean = .false. !ocean flag
 
   real(kind=r8)         :: alpha_thermal = 1.0 !Thermal expansion coefficient
   real(kind=r8)         :: alpha_salinity = 1.0 !Salinity contraction coefficient
   logical               :: lbalance_pressure = .true. !Balance pressure gradient with gravity due to constant density
 
 
   real(kind=r8)         :: SIPG_constant = 0
 
   integer, dimension(2)  :: x_boundary = 4
   integer, dimension(2)  :: y_boundary = 4
   integer, dimension(2)  :: z_boundary = 4
   integer                :: x_periodic = 0
   integer                :: y_periodic = 0
   integer                :: z_periodic = 0
   real(kind=r8)          :: robin_bc_alpha = 1.0 ! sets value for alpha in n \cdot \nabla q - \alpha q = \beta Robin BC
   logical                :: lfree_slip_exact = .false.
   integer                :: refinement_levels_h = 0
   integer                :: refinement_levels_v = 0
   integer                :: read_external_grid_flg = 0
   integer                :: full_stress_flg = 0
   integer                :: is_non_conforming_flg = 0
   logical                :: interp_cg_flux_flg = .false.
   integer                :: max_mesh_lvl = 0
   logical                :: nc_box_invert = .false.
   integer                :: p4est_log_level = 6
   real                   :: xlim_min = 0
   real                   :: xlim_max = 0
   real                   :: ylim_min = 0
   real                   :: ylim_max = 0
   real                   :: zlim_min = 0
   real                   :: zlim_max = 0
   integer, dimension(10) :: amr_indicator_variables = 0
   real, dimension(2)     :: amr_smoothness_limits = (/ 2.0, 3.0/)
   real, dimension(2)     :: amr_max_min_lim = (/ 1e-1, 1e-6/)
   real, dimension(10)    :: amr_threshold_lim = 1e16
   real, dimension(10)    :: amr_smoothness_qL2_limit = 1e-8
   logical                :: amr_mark_modes=.true.
   logical                :: amr_mark_modes_use_baseline_decay=.true.
   logical                :: amr_mark_max_min=.false.
   logical                :: amr_mark_random=.false.
   logical                :: amr_mark_threshold=.false.
   integer(C_INT32_T)     :: amr_num_neigh_iter = 3
   logical                :: amr_mark_set2nc = .false.
 
   integer               :: nel_root_h = 0
   integer               :: nel_root_v = 0
   integer               :: visc_dyn_flg     = 1   !This flag activates the laplacian computation for the dynamics equations (active by default)
   integer               :: visc_tracers_flg = 1   !This flag activates the laplacian computation for the tracers equations (active by default)
   integer               :: naza_tracers_flg = 0   !This flag activates the nazarov-LES model for the tracers equations (NOT active by default)
   real(kind=r8)         :: vertical_viscosity = 1.0 !1=On and 0=Off
   integer               :: filter_tracers_flg = 0 !This flag activates the filter computation for the tracers equations (INactive by default)
   integer               :: sounding_ncolumns = 5  !This indicates the number of columns of the sounding file
   !The default has 5 columns only and the order is the following:
   !    z(m)  theta(K)  qvapor(kg/kg)  uvelo(m/s)  vvelo(m/s)
   !Any additional variable should be added as column 6, 7, etc.
   !Ex.: to add a column with pressure we would have:
   !    z(m)  theta(K)  qvapor(kg/kg)  uvelo(m/s)  vvelo(m/s)  pressure(Pa)
   integer               :: ntracers_in = 3
   integer               :: platformID = 0
   integer               :: deviceID = 0
   integer               :: platformID2 = 1
   integer               :: deviceID2 = 0
   integer               :: Nelems = 1
   integer               :: Nslices = 1
   integer               :: NslicesV = -1
   integer               :: cpus_per_node = 1
   integer               :: gpus_per_node = 1
   integer               :: threads_per_process = 0
   integer               :: platformWeight = 1
   integer               :: platformWeight2 = 1
 
   integer               :: n_corrections = 1
   integer               :: imass = 1
 
   real(kind=r8)         :: mass_stabilize_flg = 0.0
   real(kind=r8)         :: min_length_flg
   real(kind=r8)         :: max_length_flg
   real(kind=r8)         :: moist_forcing_flg
   real(kind=r8)         :: moist_nocolumn_flg
   real(kind=r8)         :: momentum_forcing_flg = 0.0
   real(kind=r8)         :: momentum_forcing_LES_flg = 0.0
   real(kind=r8)         :: damping_function_flg = 0.0
   real(kind=r8)         :: time_residual_flg = 0.0
   real(kind=r8)         :: rain_flg
   real(kind=r8)         :: buoyancy_flg
   real(kind=r8)         :: dt_safety_factor = 0.5
   real(kind=r8)         :: thetac_in = 0.0
   real(kind=r8)         :: qvc_in = 0.0
   real(kind=r8)         :: qcc_in = 0.0
   real(kind=r8)         :: qrc_in = 0.0
   real(kind=r8)         :: thetac2_in = 0.0
   real(kind=r8)         :: qvc2_in = 0.0
   real(kind=r8)         :: qcc2_in = 0.0
   real(kind=r8)         :: qrc2_in = 0.0
   real(kind=r8)         :: thetac3_in = 0.0
   real(kind=r8)         :: qvc3_in = 0.0
   real(kind=r8)         :: qcc3_in = 0.0
   real(kind=r8)         :: qrc3_in = 0.0
   real(kind=r8)         :: thetac4_in = 0.0
   real(kind=r8)         :: qvc4_in = 0.0
   real(kind=r8)         :: qcc4_in = 0.0
   real(kind=r8)         :: qrc4_in = 0.0
   real(kind=r8)         :: xradius_pert_in = -1.0
   real(kind=r8)         :: yradius_pert_in = -1.0
   real(kind=r8)         :: zradius_pert_in = -1.0
   real(kind=r8)         :: xradius2_pert_in = -1.0
   real(kind=r8)         :: yradius2_pert_in = -1.0
   real(kind=r8)         :: zradius2_pert_in = -1.0
   real(kind=r8)         :: xradius3_pert_in = -1.0
   real(kind=r8)         :: yradius3_pert_in = -1.0
   real(kind=r8)         :: zradius3_pert_in = -1.0
   real(kind=r8)         :: xradius4_pert_in = -1.0
   real(kind=r8)         :: yradius4_pert_in = -1.0
   real(kind=r8)         :: zradius4_pert_in = -1.0
 
   real(kind=r8)         :: LAM_flg
   real(kind=r8)         :: GCM_flg
 
   integer               :: steepness_pert_in = 0 ! 0: sin-shaped rtb, >0: shape of rtb given by exp(-r**steepness_pert_in)
   ! for the default domain of 1000m xradius_pert_in=yradius_pert_in=zradius_pert_in should be chosen according to following table:
   ! steepness_pert_in  0   1     2     3     4     5     6
   ! radius_pert_in     250 125.1 139.1 143.6 144.8 144.8 144.3
   ! These numbers (for steepness_pert_in>0) were computed by Andreas Mueller (March 2014) with WolframAlpha using the following command:
   ! Table[FindRoot[NIntegrate[  r (  1/300 - 1/(300 + 0.5 0.5 (1 + Cos[pi r/250]))), {r, 0, 250}] -    NIntegrate[r (1/300 - 1/(300 + 0.5 E^-(r/x)^p)), {r, 0, 250}], {x, 125.}],{p,1,6}]
   ! this makes all the different bubbles have approximately the same total buoyancy
   real(kind=r8)         :: xc_pert_in = 1.0e6
   real(kind=r8)         :: yc_pert_in = 1.0e6
   real(kind=r8)         :: zc_pert_in = 1.0e6
   real(kind=r8)         :: xc2_pert_in = 1.0e6
   real(kind=r8)         :: yc2_pert_in = 1.0e6
   real(kind=r8)         :: zc2_pert_in = 1.0e6
   real(kind=r8)         :: xc3_pert_in = 1.0e6
   real(kind=r8)         :: yc3_pert_in = 1.0e6
   real(kind=r8)         :: zc3_pert_in = 1.0e6
   real(kind=r8)         :: xc4_pert_in = 1.0e6
   real(kind=r8)         :: yc4_pert_in = 1.0e6
   real(kind=r8)         :: zc4_pert_in = 1.0e6
   real(kind=r8)         :: mount_xc   = 0.0
   real(kind=r8)         :: mount_yc   = 0.0
   real(kind=r8)         :: mount_hm   = 0.0
   real(kind=r8)         :: mount_ac   = 1.0
   real(kind=r8)         :: mount_bc   = 1.0
   character(len=24)     :: mount_type = 'mountaintypehere_xxxxxx'
   character(len=12)     :: zgrid_type = 'sigma'
   real(kind=r8)         :: time_scale = 1.0_r8
   real(kind=r8)         :: uvelo_uniform = 0.0
   real(kind=r8)         :: vvelo_uniform = 0.0
   real(kind=r8)         :: wvelo_uniform = 0.0
   real(kind=r8)         :: Mach1_in    = 0.0
   real(kind=r8)         :: Mach2_in    = 0.1
   real(kind=r8)         :: hybrid_s    = 1000.0
   real(kind=r8)         :: sleve_s1    = 1000.0
   real(kind=r8)         :: sleve_s2    = 1000.0
   real(kind=r8)         :: xmin_velo_in = -6.0e8
   real(kind=r8)         :: xmax_velo_in = 6.0e8
   real(kind=r8)         :: ymin_velo_in = -6.0e8
   real(kind=r8)         :: ymax_velo_in = 6.0e8
   real(kind=r8)         :: zmin_velo_in = -6.0e8
   real(kind=r8)         :: zmax_velo_in = 6.0e8
   real(kind=r8)         :: p00_in = -999.999
   real(kind=r8)         :: nu_air_in = -1.0
   real(kind=r8)         :: Pr = 0.7            !Prandtl number for dry air. Default value
   real(kind=r8)         :: C1_in = 1.0         !Default C1 constant in Nazarov
   real(kind=r8)         :: C2_in = 0.5         !Default C2 constant in Nazarov
   real(kind=r8)         :: Cs_in = 0.14        !Default Cs constant in Smagorinsky
   real(kind=r8)         :: Re_tau = 1.0        !Friction Reynolds number
   real(kind=r8)         :: Re_bulk = 1.0       !Bulk Reynolds number (simply, Reynolds)
   real(kind=r8)         :: Mach = 0.1          !Mach number
 
   !
   ! The values of dirchlet_* are the actual (tota) physical values of the variable at the wall at hand
   !
   real(kind=r8), dimension(6) :: dirichlet_theta !dirichlet_theta contains a value of theta per face (6 faces in a box)
   !dirichlet_theta(1) -> bottom value
   !dirichlet_theta(2) -> top    value
   !dirichlet_theta(3) -> sourth value
   !dirichlet_theta(4) -> north  value
   !dirichlet_theta(5) -> west   value
   !dirichlet_theta(6) -> east   value
   real(kind=r8), dimension(6) :: dirichlet_uvelo !dirichlet_uvelo contains a value of uvelo per face (6 faces in a box)
   !dirichlet_uvelo(1) -> bottom value
   !dirichlet_uvelo(2) -> top    value
   !dirichlet_uvelo(3) -> sourth value
   !dirichlet_uvelo(4) -> north  value
   !dirichlet_uvelo(5) -> west   value
   !dirichlet_uvelo(6) -> east   value
   real(kind=r8), dimension(6) :: dirichlet_vvelo !dirichlet_vvelo contains a value of vvelo per face (6 faces in a box)
   !dirichlet_vvelo(1) -> bottom value
   !dirichlet_vvelo(2) -> top    value
   !dirichlet_vvelo(3) -> sourth value
   !dirichlet_vvelo(4) -> north  value
   !dirichlet_vvelo(5) -> west   value
   !dirichlet_vvelo(6) -> east   value
   real(kind=r8), dimension(6) :: dirichlet_wvelo !dirichlet_wvelo contains a value of wvelo per face (6 faces in a box)
   !dirichlet_wvelo(1) -> bottom value
   !dirichlet_wvelo(2) -> top    value
   !dirichlet_wvelo(3) -> sourth value
   !dirichlet_wvelo(4) -> north  value
   !dirichlet_wvelo(5) -> west   value
   !dirichlet_wvelo(6) -> east   value
 
   real(kind=r8)         :: c_am2 = 3.0/2.0 !AM2 constant from Durran-Blossey (2nd order)
 
   real(kind=r8)         :: kessler_production_constant = 0.1
 
   character(len=4)      :: bcast_type = 'mpi'
   character(len=12)     :: tau_type   = 'codina'
   character(len=12)     :: shock_type = 'codina'
   character(len=12)     :: eqn_set    = 'set2nc'
   character(len=12)     :: pert_type  = 'default'
   character(len=12)     :: pert_type_tracers  = 'default'
   character(len=12)     :: equations  = 'euler'
   character(len=12)     :: limiter_type = 'none'
   character(len=3)      :: piso_type = 'dsa'
 
   logical :: lread_input_dimensions = .false.
   logical :: lleveque       = .false.
   logical :: lburgers1d     = .false.
   logical :: lburgers2d     = .false.
   logical :: lsmolar_flow   = .false.
   logical :: lslotted       = .false.
   logical :: ldynamics      = .true.
   logical :: lpure_advection= .false.
   logical :: lphysics       = .false.
   logical :: lsimple_physics= .false.
   logical :: lkessler       = .false.
   logical :: lrain          = .true.
   logical :: lmoist_forcing = .true.
   logical :: lmoist_column  = .false.
   logical :: lpassive       = .false.
   logical :: lread_sound    = .false.
   logical :: lpert_theta    = .false.
   logical :: lpert_qv       = .false.
   logical :: lpert_qc       = .false.
   logical :: lpert_qr       = .false.
   logical :: lpert_velo     = .false.
   logical :: ljet           = .false.
   logical :: lcylinder      = .false.
   logical :: ldam           = .false.
   logical :: lread_qvapor   = .false.
   logical :: lread_uvelo    = .false.
   logical :: lread_vvelo    = .false.
   logical :: lrestart_file  = .false. !obsolete
   logical :: lmountain      = .false.
   logical :: lout_rho       = .false.
   logical :: lout_uvelo     = .false.
   logical :: lout_vvelo     = .false.
   logical :: lout_wvelo     = .false.
   logical :: lout_velo      = .true.
   logical :: lout_dvelo     = .false.
   logical :: lout_theta     = .true.
   logical :: lout_thetae    = .false.
   logical :: lout_press     = .true.
   logical :: lout_temperature=.false.
   logical :: lout_previous_time_steps = .false.
   logical :: lout_rank       = .false.
   logical :: lout_vorticity  = .false.
   logical :: lout_xvorticity = .false.
   logical :: lout_yvorticity = .false.
   logical :: lout_zvorticity = .false.
   logical :: lout_rvorticity = .false.
   logical :: lout_tke       = .false.
   logical :: lout_qvapor    = .false.
   logical :: lout_qcloud    = .false.
   logical :: lout_qrain     = .false.
   logical :: lout_sponge    = .false.
   logical :: lout_tau       = .false.
   logical :: lout_radius    = .false.
   logical :: lout_coords    = .false.
   logical :: lout_spherical = .false.
   logical :: lwrite_time_averaging = .false.
   logical :: lout_statistics     = .false.
   logical :: lsommerfeld    = .false.
   logical :: lxstretch      = .false.
   logical :: lystretch      = .false.
   logical :: lzstretch      = .false.
   logical :: lVMS           = .false.
   logical :: lVMS_dyn       = .false.
   logical :: lshock_cpt     = .false.
   logical :: lshock_cpt_dyn = .false.
   logical :: lapply_bcs_velocity = .false.
   logical :: lLAV           = .false.
   logical :: lnazarov       = .false.
   logical :: lnazarov_tracers=.false.
   logical :: llimit         = .false.
   logical :: llimit_below   = .false.
   logical :: llimit_above   = .false.
   logical :: lvisc_reciprocal = .false.
   logical :: lvisc_anisotropic = .false.
   logical :: lvisc_dynamics = .true.
   logical :: lvisc_tracers  = .true.
   logical :: ltime_residual = .false.
   logical :: lforce_turb_momentum = .false.
   logical :: lanisotropic_laplacian = .false.
   logical :: lviscous_boundary = .false.
   logical :: lviscous_ground   = .false.
   logical :: lviscous_bottom   = .false.
   logical :: lviscous_top      = .false.
   logical :: lviscous_south    = .false.
   logical :: lviscous_north    = .false.
   logical :: lviscous_east     = .false.
   logical :: lviscous_west     = .false.
   logical :: lstabilize_mass_eqn = .false.
   logical :: lLES           = .false.
   logical :: lSMAG          = .false.
   logical :: ladd_full_stress_tensor = .false.
   logical :: ldamping_function = .false.
   logical :: luse_PSUP      = .false.
   logical :: lgrid_only     = .false.
   logical :: lp4est         = .false.
   logical :: lp6est         = .false.
   logical :: lread_external_grid = .false.
   logical :: lset_outlet_pressure = .false.
   logical :: lread_bc       = .false.
   logical :: lserial_grid_creation = .true.
   logical :: lparallel_grid_creation = .false.
   logical :: lio_grid_ascii = .false.
   logical :: lwrite_grid_ascii = .false.
   logical :: lgravity       = .true.
   logical :: lout_ascii     = .false.
   logical :: lout_asciimaya = .false.
   logical :: lcompute_barycenter = .false.
   logical :: luse_min_element_length = .false.
   logical :: lrotating_flow = .false.
 
   logical :: ldirichlet_all    = .false.
   logical :: ldirichlet_top    = .false.
   logical :: ldirichlet_bottom = .false.
   logical :: ldirichlet_east   = .false.
   logical :: ldirichlet_west   = .false.
   logical :: ldirichlet_north  = .false.
   logical :: ldirichlet_south  = .false.
 
   logical :: lforce_spongex   = .false.
   logical :: lforce_spongey   = .false.
   logical :: lforce_spongez   = .false.
 
   logical :: ladapt_timestep = .false.
   logical :: lpi_fraction    = .false.
 
   logical :: lincompressible = .false.
   logical :: lpassive_advection = .false.
 
   !Input parameters to write a spherical shell (1 level) to a netcdf file:
   logical :: lout_spherical_shell = .false.
   logical :: lout_nc_3d = .true.
   logical :: dump_rhs = .false.
   logical :: lgpu = .false.
   logical :: luse_hybrid_cpu_gpu = .false.
 
   logical :: lhelmholtz = .false.
   logical :: lALE = .false.
   logical :: is_shallow = .false.
   integer :: nonlinear_swe = 1
   logical :: is_swe_layers = .false.
   logical :: is_sphere = .false.
   logical :: is_implicit = .false.
   logical :: is_imex = .false.
   logical :: is_1d = .false.  
   
   logical :: is_mlswe = .false.  !added by Yao Gahounzo
   logical :: matlab_viz = .false. !added by Yao Gahounzo
   logical :: is_mlswe_linear = .false. !added by Yao Gahounzo
   logical :: mass_exact = .false. !added by Yao Gahounzo
   logical :: mlswe_bc_strong = .false. !added by Yao Gahounzo
   logical :: dg_integ_exact = .true. ! added by Yao Gahounzo
   logical :: dump_data = .true. ! added by Yao Gahounzo

 
   !-----------------------------------------------------------------------
   ! Parameters
   !-----------------------------------------------------------------------
   character(len=*), parameter :: namelist_input='numo3d.in'
 
   !-----------------------------------------------------------------------
   ! Split-Explicit Parameters from Alex Reinecke
   !-----------------------------------------------------------------------
   real(kind=r8)      :: ti_alpha     = 0.0  !this is \theta from the writeup, we use ti_alpha=0.2
   real(kind=r8)      :: ti_beta      = 0.0  !this is \beta_d from the writeup, we use \beta_d = 0.2
   real(kind=r8)      :: dd_alpha     = 0.0 !this is the divergence damping coefficient for the explicit form of divergence damping.  We do not use this.
   integer            :: nsmall        = 6 !this is the number of small time steps per large step. 
 
   !-----------------------------------
   ! Shallow water variables
   !-----------------------------------
   real xc_in, yc_in, zc_in, rc_in
   real xc2_in, yc2_in, zc2_in, rc2_in
   real xc3_in, yc3_in, zc3_in, rc3_in
   real xc4_in, yc4_in, zc4_in, rc4_in
   real xc5_in, yc5_in, zc5_in, rc5_in
   real xc6_in, yc6_in, zc6_in, rc6_in
   real xc7_in, yc7_in, zc7_in, rc7_in
   real h0_in, beach_length_in, xmin_topography_in
   real xmin_in, xmax_in, ymin_in, ymax_in, zmin_in, zmax_in
   real x_obstacle_offset_in
   real xmax_obstacle, xmin_obstacle
   real ymax_obstacle, ymin_obstacle
   real xmax_water, xmin_water
   real ymax_water, ymin_water
 
   real xstretch_ref_in, zstretch_ref_in
 
   real, dimension(4) :: xdims_obstacle, ydims_obstacle
 
   real, allocatable :: amplitudes_in(:)
   real, allocatable :: radiusx_in(:), radiusy_in(:)
   real, allocatable :: xcc_in(:), ycc_in(:)
 
   real slope_in
   real initial_shore_line_in
   real time_scale_in
   real wave_amplitude_in
   real wave_min_crest_in
   real wave_max_crest_in
   real carrier_scaling_in
   real carrier_x_offset_in
   real wave_start_in
   real wave_xcenter_in
   real wave_ycenter_in
   real Lambda_in
   real distance_lambda_fraction_in
   real lambda_fraction_in !This value indicates a fraction of the wave_length
   real obstacle_crosssection_lambda_fraction_in
   real amplitude_height_fraction_in
   real amplitude_in, amplitude2_in, amplitude3_in
   real amplitude4_in, amplitude5_in, amplitude6_in, amplitude7_in
 
   real dam_depth_in
   real free_surface_level_in
   real, dimension(2) :: dam_xlimit_in, dam_ylimit_in
 
   real bathymetry_shift
   integer nobstacles_in
   integer hump_config
   integer ibathymetry
   integer nonlinear, balance_flag, output_flag, warp_grid
   character mesh_file*100, bathymetry_file*100, wave_file*100
 
   logical lout_tree
   logical lout_shoreline
   logical lwave
   logical lread_wave_from_file
   logical lread_external_bathy
   character wave_type*12
   real gravity_in
   real xdam_min, xdam_max
   real ydam_min, ydam_max
 
   logical luniform_flow
   real    flow_uvelocity_in
 
   real    ntimes_Lobs_in
   real    Lobs_in
   real    Hobs_in, Hobs2_in, Hobs3_in, Hobs4_in, Hobs5_in
   real    CXobs_in
 
   real :: limit_threshold = 1e-3
   integer :: limiter_qoi = 0
   real :: cms_coefficient = 0.0
   real :: cms_coefficient2 = 0.0
 
   logical llinear_pert
   logical lobstacle
   logical lobstacle_semicircle
   logical lobstacle_ellipsis
   logical lobstacle_size_as_lambda_fraction
   logical lobstacle_crosssection_as_lambda_fraction
   logical lobstacle_distance_as_lambda_fraction
   logical lobstacle_distance_as_obstacle_length
   logical lobstacle_height_as_amplitude_fraction
 
   !-----------------------------------
   ! End shallow water variables 
   !-----------------------------------
 
 contains
 
   !-----------------------------------------------------------------------
   subroutine mod_input_create(irank)
 
     implicit none
 
     integer irank, icheck
 
     integer :: funit
 
     !Namelist Input
 
     namelist /eqnset/ eqn_set, equations, lincompressible, piso_type, lpassive_advection, lhelmholtz, lALE, is_swe_layers, nonlinear_swe, is_mlswe
 
     namelist /input/ dt, dt_safety_factor, &
          restoring_time, &
          lrestoring_sponge, &
          lpi_fraction, &
          time_initial, time_final, time_dynamic_amr, time_restart, time_scale, irestart_file_number,&
          measurement_interval, &
          interval_start, &
          interval_end,   &
          lrestart_file, icase, &
          case9015_inlet_vel,&
          case9101_initial_time,&
          ti_method, ark2_type, si_method, si_dimension, &
          ti_method_btp, & ! Added by Yao G.
          Iter_Type, Proj_Size,Focal_Dist_Adj, precon_mode, precon_order, &
          precon_type, PBNO_pNorm, Imag_Wt, Hmat_Size, &
          c_am2, lf2_filter, lf2_raw_filter, lf2_alpha, lf2_nu, kstages, solver_type, lcg_proj, &
          solver_tol, solver_tol_pressure, solver_tol_ALE, jfnk_tol, epsilon_jfnk, gamma1, delta, &
          filter_mux, filter_muy, filter_muz, ifilter, &
          filter_weight_type, filter_basis_type, fname_root, out_type, lout_ascii, lout_asciimaya, format_vtk, nvtk_files, vtk_cell_type, &
          write_mesh, &
          fname_initial, limit_solver_iter_press, limit_solver_iter_ALE,&
          restart_path, & ! Added by Yao G.
          lgravity, &
          visc, visc2, visc4, &
          visc_h, visc_v, &
          viscx, &
          viscy, &
          viscz, &
          diff_T, diff_Th, diff_Tv, &
          diff_S, diff_Sh, diff_Sv, &
          lsalinity, &
          which_eos, &
          locean, &
          lbalance_pressure, &
          alpha_thermal, alpha_salinity,&
          SIPG_constant, &
          lvisc_anisotropic, &
          lvisc_dynamics, visc_dyn_flg,     &
          lvisc_tracers, visc_tracers_flg, &
          lnazarov_tracers, naza_tracers_flg, &
          vertical_viscosity, &
          filter_tracers_flg, &
          lVMS, lshock_cpt,   &
          lVMS_dyn, lshock_cpt_dyn,   &
          tau_type, shock_type, &
          Pr, C1_in, C2_in, Cs_in, &
          Re_tau, &
          Re_bulk, &
          Mach, &
          lLAV, &
          lnazarov, &
          llimit, llimit_below, llimit_above, limiter_type, limiter_qoi, &
          ltime_residual, &
          lLES, &
          lSMAG, &
          ladd_full_stress_tensor, &
          ldamping_function, &
          damping_function_flg, &
          luse_PSUP, &
          lforce_turb_momentum, &
          momentum_forcing_flg, &
          momentum_forcing_LES_flg, &
          luse_min_element_length, &
          min_length_flg, max_length_flg, &
          lstabilize_mass_eqn, &
          mass_stabilize_flg, &
          nlaplacian, &
          nlaplacian_flg, &
          lanisotropic_laplacian, &
          ladapt_timestep,lprint_diagnostics,iprint_diagnostics,lcheck_all_conserved,lunit_test,unit_test_type,lcompute_barycenter,ljenkins_test, &
          ldynamics, lphysics, lsimple_physics, lkessler, lrain, lmoist_forcing, lmoist_column, lpassive, lpure_advection, ntracers_in, &
          kessler_production_constant, &
          sounding_file, lread_sound, lread_qvapor, lread_uvelo, lread_vvelo, &
          ctd_file_temp, ctd_file_salinity, &
          sounding_ncolumns, &
          uvelo_uniform, vvelo_uniform, wvelo_uniform, &
          Mach1_in, Mach2_in, &
          lrotating_flow, &
          bcast_type, &
          lout_rho, lout_uvelo, lout_vvelo, lout_wvelo, lout_velo, lout_dvelo, lout_theta, lout_sponge, &
          lout_thetae, lout_press, lout_temperature, &
          lout_previous_time_steps, &
          lout_rank,  &
          lout_vorticity,  &
          lout_xvorticity, &
          lout_yvorticity, &
          lout_zvorticity, &
          lout_rvorticity, &
          lout_tke, &
          lout_qvapor, lout_qcloud, lout_qrain, &
          lout_tau, lout_radius, lout_coords, &
          lout_spherical, &
          lwrite_time_averaging, &
          lout_statistics, &
          p00_in, &
          nu_air_in, &
          pert_type, &
          pert_type_tracers, &
          lleveque, &
          lburgers1d, &
          lburgers2d, &
          lsmolar_flow, &
          lslotted, &
          lout_spherical_shell, &
          zlevel_out, &
          lout_nc_3d, &
          lpert_theta, lpert_qv, lpert_qc, lpert_qr, lcylinder, &
          ldam, &
          lpert_velo, &
          ljet, &
          thetac_in, qvc_in, qcc_in, qrc_in, &
          thetac2_in, qvc2_in, qcc2_in, qrc2_in, &
          thetac3_in, qvc3_in, qcc3_in, qrc3_in, &
          thetac4_in, qvc4_in, qcc4_in, qrc4_in, &
          steepness_pert_in, &
          steepness_pert_in, &
          xradius_pert_in, yradius_pert_in, zradius_pert_in, &
          xradius2_pert_in, yradius2_pert_in, zradius2_pert_in, &
          xradius3_pert_in, yradius3_pert_in, zradius3_pert_in, &
          xradius4_pert_in, yradius4_pert_in, zradius4_pert_in, &
          xc_pert_in, yc_pert_in, zc_pert_in, &
          xc2_pert_in, yc2_pert_in, zc2_pert_in, &
          xc3_pert_in, yc3_pert_in, zc3_pert_in, &
          xc4_pert_in, yc4_pert_in, zc4_pert_in, &
          xmin_velo_in, xmax_velo_in, &
          ymin_velo_in, ymax_velo_in, &
          zmin_velo_in, zmax_velo_in, &
          space_method,  cgdg_method, form_method, dump_rhs, &
          lgpu, numaocca_dir, Nelems, Nslices, NslicesV, vectorization, &
          platform, platformID, deviceID, platformWeight, platform2, platformID2, deviceID2, platformWeight2, &
          cpus_per_node, gpus_per_node, threads_per_process, luse_hybrid_cpu_gpu, &
                                 !!--Shallow
          ibathymetry, bathymetry_file, wave_file, &
          gravity_in, &
          xc_in, yc_in, rc_in, &
          xc2_in, yc2_in, rc2_in, &
          xc3_in, yc3_in, rc3_in, &
          xc4_in, yc4_in, rc4_in, &
          xc5_in, yc5_in, rc5_in, &
          xc6_in, yc6_in, rc6_in, &
          xc7_in, yc7_in, rc7_in, &
          beach_length_in, mesh_file, &
          xmin_topography_in, h0_in, free_surface_level_in, &
          lobstacle, &
          lobstacle_semicircle, &
          lobstacle_ellipsis, &
          lobstacle_distance_as_lambda_fraction, &
          lobstacle_distance_as_obstacle_length, &
          lobstacle_height_as_amplitude_fraction, &
          lobstacle_size_as_lambda_fraction, &
          lobstacle_crosssection_as_lambda_fraction, &
          Lobs_in, CXobs_in, Hobs_in, &
          Hobs2_in,  & !obstacle height
          Hobs3_in,  & !obstacle height
          Hobs4_in,  & !obstacle height
          Hobs5_in,  & !obstacle height
          ntimes_Lobs_in, &
          xdims_obstacle, ydims_obstacle, &
          xmin_obstacle, xmax_obstacle,   &
          ymin_obstacle, ymax_obstacle,   &
          x_obstacle_offset_in, &
          slope_in, &
          initial_shore_line_in, &
          dam_depth_in, dam_xlimit_in, &
          lwave, lread_external_bathy, lread_wave_from_file, wave_amplitude_in, &
          Lambda_in, lambda_fraction_in, obstacle_crosssection_lambda_fraction_in, wave_start_in, wave_type, &
          wave_min_crest_in, wave_max_crest_in, &
          carrier_scaling_in, carrier_x_offset_in, &
          wave_xcenter_in, wave_ycenter_in, &
          distance_lambda_fraction_in, &
          amplitude_height_fraction_in, &
          amplitude_in, amplitude2_in, amplitude3_in, &
          amplitude4_in, amplitude5_in, amplitude6_in, amplitude7_in, &
          hump_config, bathymetry_shift, limit_threshold, &
          amplitudes_in, nobstacles_in, radiusx_in, radiusy_in, xcc_in, ycc_in, &
          cms_coefficient, cms_coefficient2, n_corrections, imass, &
          ad_mlswe, explt_coriolis, cd_mlswe, dp_cutoff1, dp_cutoff2, dp_tau_bot, dp_tau_wind, dt_btp,method_visc,&
          visc_mlswe, dpprime_visc_min, max_shear_dz, matlab_viz, adjust_H_vertical_sum, is_mlswe_linear, botfr, &
          mass_exact, bcl_flux, mlswe_bc_strong, dg_integ_exact, dump_data
 
     namelist /gridnl/ nelx, nely, nelz, nopx, nopy, nopz, xdims, ydims, ztop, zbottom, &
          nlayers, &
          lread_input_dimensions, &
          delta_domain, & !used for turbulent channel flows
          nlon, nlat, &
          geometry_type, decomp_type, nproc_z, &
          x_boundary, y_boundary, z_boundary, &
          x_periodic, y_periodic, z_periodic, &
          robin_bc_alpha, &
          lset_outlet_pressure, &
          lfree_slip_exact, &
          bc_tscale, bc_xscale, bc_yscale, bc_zscale, &
          ldirichlet_all, &
          ldirichlet_top, ldirichlet_bottom,  &
          ldirichlet_east, ldirichlet_west,   &
          ldirichlet_north, ldirichlet_south, &
          dirichlet_theta, &
          dirichlet_uvelo, &
          dirichlet_vvelo, &
          dirichlet_wvelo, &
          lviscous_boundary, &
          lviscous_ground, &
          lviscous_top, &
          lviscous_south, &
          lviscous_north, &
          lviscous_east, &
          lviscous_west, &
          lforce_spongex, &
          lforce_spongey, &
          lforce_spongez, &
          sponge_type, sponge_top_coe, sponge_lateralx_coe,sponge_lateralx_coe_east, sponge_lateralx_coe_west, sponge_lateraly_coe, &
          lsommerfeld, &
          lgrid_only, &
          lp4est, lp6est, lio_grid_ascii, &
          lread_external_grid, &
          is_non_conforming_flg, &
          interp_cg_flux_flg, &
          max_mesh_lvl, &
          nc_box_invert, &
          p4est_log_level, &
          xlim_min, xlim_max, ylim_min, ylim_max, zlim_min, zlim_max, &
          amr_indicator_variables, &
          amr_smoothness_limits, &
          amr_max_min_lim, &
          amr_threshold_lim, &
          amr_smoothness_qL2_limit, &
          amr_mark_max_min, &
          amr_mark_random, &
          amr_mark_threshold, &
          amr_mark_modes, &
          amr_mark_modes_use_baseline_decay, &
          amr_num_neigh_iter, &
          amr_mark_set2nc, &
          lread_bc, &
          lserial_grid_creation, lparallel_grid_creation, lio_grid_ascii, lwrite_grid_ascii,&
          refinement_levels_h, refinement_levels_v, nel_root_v, nel_root_h,&
          xstretch_coe, ystretch_coe, zstretch_coe, &
          lxstretch, lystretch, lzstretch, &
          xstretch_type, ystretch_type, zstretch_type, &
          lmountain,mount_type,mount_xc,mount_yc,mount_hm,mount_ac,mount_bc,&
          zgrid_type,hybrid_s,sleve_s1,sleve_s2
 
     amr_indicator_variables(1) = 1
 
     !Read input namelist
     funit = get_unit()
     print*,"File:",namelist_input, funit
     open(funit,file=namelist_input)
     read(funit,input)
     close(funit)
 
     !Read eqnset namelist
     funit = get_unit()
     open(funit,file=namelist_input)
     read(funit,eqnset)
     close(funit)
 
     !Read gridnl namelist
     funit = get_unit()
     open(funit,file=namelist_input)
     read(funit,gridnl)
     close(funit)
 
     if(interp_cg_flux_flg .and. &
          .not. (form_method(1:17) == 'strong_no_product'  &
           .or.  form_method(1:4) == 'skew')) &
        stop "interp_cg_flux_flg requires form_method = 'strong_no_product'"
     geometry_type = lowercase(geometry_type)
     decomp_type = lowercase(decomp_type)
     ti_method = lowercase(ti_method)
     ti_method_btp = lowercase(ti_method_btp) !added by YaoG
     ark2_type = lowercase(ark2_type)
     si_method = lowercase(si_method)
     si_dimension = lowercase(si_dimension)
     Iter_Type = uppercase(Iter_Type)
     precon_mode = uppercase(precon_mode)
     precon_type = uppercase(precon_type)
     solver_type = lowercase(solver_type)
     filter_weight_type = lowercase(filter_weight_type)
     filter_basis_type = lowercase(filter_basis_type)
     out_type = lowercase(out_type)
     format_vtk = uppercase(format_vtk)
     vtk_cell_type = uppercase(vtk_cell_type)
     sponge_type = lowercase(sponge_type)
     pert_type   = lowercase(pert_type)
     pert_type_tracers = lowercase(pert_type_tracers)
 
     !Add DT to FNAME
     write(real_string,'(f9.4)') dt
     if (dt>=1000) then
        fname_root=trim(fname_root) // '_' // trim(real_string(1:9))
     else if (dt>=100) then
        fname_root=trim(fname_root) // '_' // trim(real_string(2:9))
     else if (dt>=10) then
        fname_root=trim(fname_root) // '_' // trim(real_string(3:9))
     else
        fname_root=trim(fname_root) // '_' // trim(real_string(4:9))
     end if
 
     if(lburgers1d) then
        pi_trig = 3.1415926535897932346_r8
        if(lpi_fraction) then
           dt = dt/pi_trig
           time_final = time_final/pi_trig
           time_restart = time_restart/pi_trig
        end if
     end if
 
     !More Viscosity Options
     if (visc>0 .and. nlaplacian==1) then
        visc2=visc
     end if
     if (visc>0 .and. nlaplacian==2) then
        visc4=visc
     end if
     if (visc2 > 0) then
        visc=visc2
        nlaplacian=1
     end if
     if (visc4 > 0) then
        visc=visc4
        nlaplacian=2
     end if
 
     if (.not.lvisc_anisotropic) then
        visc_h=visc
        visc_v=visc
     else if (lvisc_anisotropic) then
        visc=max(visc_h,visc_v)
     end if    
 
     !If full stress, then laplacian should not be added
     if(ladd_full_stress_tensor) then
        full_stress_flg = 1
        nlaplacian = 0
     end if
         
     !Default values for viscosity coefficient:
     if(nlaplacian > 0) then
        
        !If the anisotropic values of viscx,viscy,viscz
        !are > 0.0, then, force visc=1.0 so that it is
        !not accounted for in the computation and only
        !only viscx, viscy, and viscz are used:
 
        if(viscx > 0.0 .or. viscy > 0.0 .or. viscz > 0.0) then
           lanisotropic_laplacian=.true.
           visc = 1.0
        else
           !Set viscx,viscy,viscz to 1.0 so that the only visc coefficient
           !that is accounted for by default is visc...
           viscx = 1.0
           viscy = 1.0
           viscz = 1.0
 
           if(visc <= 0.0) then
              print*, ' STOP: YOU DEFINED nlaplacian= ', nlaplacian, ' but you did not gave a positive value to visc.'
              print*, ' Assign a value to visc in your input file and RELAUNCH the program.'
              print*, ' The program will STOP here (mod_input).'
              stop
           end if
           !in this case, visc is the only valueviscosity coefficient that
           !is used and it has the value given by the use in input.
           !the
        end if
 
 !       if (visc2 > 0 .or. visc4 > 0) then
 !          print*, ' STOP: YOU DEFINED nlaplacian= ', nlaplacian, ' but also defined VISC2 and VISC4.'
 !          print*, ' Choose either VISC2 and VISC4, or VISC and NAPLACIAN in the input file and RELAUNCH the program.'
 !          print*, ' The program will STOP here (mod_input).'
 !          stop
 !       end if
 
        !
        ! nlaplacian_flg is used by compute_strain_derivative_IBP
        ! to subtract dtau_ii/dx_i from the full stress if AD is also used
        nlaplacian_flg = 2       
     end if
     
     !Define VISC_T: if negative, then should be equal to VISC
     if (diff_T < 0) then
        diff_T=visc
     end if
 
     !Define VISC_S: if negative, then should be equal to VISC
     if (diff_S < 0) then
        diff_S=visc
     end if
 
     if (lincompressible) then
        if(.not.lvisc_anisotropic) then
           diff_Th=diff_T
           diff_Tv=diff_T
           diff_Sh=diff_S
           diff_Sv=diff_S
        else if(diff_Th<0.or.diff_Tv<0.or.diff_Sh<0.or.diff_Sv<0) then
           print*, " STOP: You asked for anisotropic thermal or salinity diffusivity, but did not specify values anisotropic."
           print*, ' The program will STOP here (mod_input).'
           stop
        end if
     end if
 
     !Switch off salinity when not doing incompressible
     if(.not.(lhelmholtz.and.piso_type=='fxg')) then
        lsalinity=.false.
        locean=.false.
     end if
 
     !some LAV options
     if (lLAV) then
        
        if (visc < epsilon(1.0)) then
           visc=epsilon(1.0) !touch viscosity to make it greater than 1 to form necessary arrays
        end if
        
        if (.not.ladd_full_stress_tensor) then
           if (nlaplacian == 0) then
              if (irank == 0) then
                 print*,' Error in MOD_INPUT; Incompatible input data.'
                 print*,' You selected lLAV=True but did not specificy NLAPLACIAN'
              end if
              stop
           end if
        end if
        
     end if
 
     !Add EQN_SET to FNAME
     if (eqn_set == 'set2nc' .or. eqn_set == 'set2c' .or. eqn_set == 'set3c') then
        fname_root=trim(fname_root) // '_' // trim(eqn_set)
     else
        if (irank == 0) then
           print*,' Error in MOD_INPUT; Incompatible input data.'
           print*,' eqn_set = ',eqn_set
        end if
        stop
     end if
 
     !Add SPACE_METHOD to FNAME
     if (space_method == 'cgc' .or. space_method == 'cgd' .or. space_method == 'dg') then
        fname_root=trim(fname_root) // '_' // trim(space_method)
     else
        if (irank == 0) then
           print*,' Error in MOD_INPUT; Incompatible input data.'
           print*,' space_method = ',space_method
        end if
        stop
     end if
 
     !Add TI_METHOD to FNAME
     if(is_mlswe) then
         if (ti_method_btp /= 'rk35' .and. ti_method_btp /= 'rk34' .and. ti_method_btp /= 'lsrk') then
            ti_method_btp = 'explicit'
         endif
        fname_root=trim(fname_root) // '_' // trim(ti_method) // '_' // trim(ti_method_btp)
     elseif (ti_method(1:2) == 'rk') then
        fname_root=trim(fname_root) // '_' // trim(ti_method)
        si_dimension = '0d' !to avoid issues with IMEX solver arrays being built
     else
        fname_root=trim(fname_root) // '_' // trim(ti_method) // '_' // trim(si_method) // '_' // trim(si_dimension)
     end if
 
     !Preconditioning File
     if (si_dimension == '3d' .and. delta >= 0) then
        precon_fname=trim(fname_root) // '.precon'
     end if
     if (precon_order == 0) precon_mode = 'SKIP'
     if (precon_mode == 'SKIP') precon_order = 0
 
     !Setup some Solver Defaults
     if (si_dimension == '1d') then
        !       solver_tol=1e-9
     else if (si_dimension == '3d') then
        !       solver_tol=1e-2
     end if
     if (si_method == 'schur'.and.(.not.lincompressible)) then
        gamma1=+1
     else if (si_method == 'no_schur') then
        gamma1=-1
     end if
 
     !Some options for grid-dependent schemes such as Nazarov
     min_length_flg = 0.0
     max_length_flg = 1.0
     if(luse_min_element_length .or. icase==71 .or. icase==7) then
        min_length_flg = 1.0
        max_length_flg = 0.0
     end if
     mass_stabilize_flg = 0.0
     if(lstabilize_mass_eqn) mass_stabilize_flg = 1.0
 
     !Some options for turbulence forcing:
     momentum_forcing_flg = 0.0
     momentum_forcing_LES_flg = 0.0
     if(lforce_turb_momentum) then
        if(lSMAG .or. lnazarov .or. lLES) then
           momentum_forcing_flg = 0.0
           momentum_forcing_LES_flg = 1.0
        else
           momentum_forcing_flg = 1.0
           momentum_forcing_LES_flg = 0.0
        end if
     end if
     
     time_residual_flg = 0.0
     if(ltime_residual) time_residual_flg = 1.0
     
     !Set LAM_flg and GCM_flg
     if (geometry_type == 'cube') then
        LAM_flg = 1.0
        GCM_flg = 0.0
     else if(geometry_type(1:6) == 'sphere') then
        LAM_flg = 0.0
        GCM_flg = 1.0
        is_sphere=.true.
     end if
 
     !Damping function for turbulence modeling (wall laws):
     damping_function_flg = 0.0
     if(ldamping_function) damping_function_flg = 1.0
 
     !Options for cloud microphysics like Kessler:
     buoyancy_flg = 1.0
     if(lpassive) buoyancy_flg = 0.0
     moist_nocolumn_flg = 1.0
     if(lmoist_column) moist_nocolumn_flg = 0.0
     if(.not. lkessler) lrain = .false.
     if(lkessler .or. lphysics) then
        lpassive = .false.
     end if
 
     ! lmoist_column always includes rain:
     if(.not. lrain) lmoist_column = .false.
     moist_forcing_flg = 0.0
     if(lmoist_forcing) moist_forcing_flg = 1.0
     rain_flg = 0.0
     if(lrain) rain_flg = 1.0
 
     !Setup some Solver checks
     if (si_dimension == '3d') then
        solver_type = 'iterative'
     else if (si_dimension == '1d') then
        !       solver_type = 'direct'
        !Iter_Type   = 'LU'
     end if
 
     !
     ! Perform some checks on the Initial Data
     !
     if (irank == 0) then
        if ( (si_dimension == '1d') .and. (decomp_type == 'metis3d' .or. &
             decomp_type == 'geom3' .or. &
             decomp_type == 'geom4') ) then
           print*,' Error in MOD_INPUT; Incompatible input data.'
           print*,' si_dimension = ',si_dimension
           print*,' decomp_type = ',decomp_type
           stop
        end if
 
        !Check Equation Sets and IMEX
        if (delta >= 0 .and. ti_method(1:2) /= 'rk' .and. eqn_set(1:6) == 'set3c') then
           print*,' Error in MOD_INPUT; Incompatible input data.'
           print*,' You asked for IMEX but IMEX only works for SET2NC and SET2C'
           stop
        end if
 
        !Check Equation Sets
        if (eqn_set(1:6) == 'set2nc' .or. eqn_set(1:5) == 'set2c' .or. eqn_set(1:5) == 'set3c') then
           !Everything OK
        else
           print*,' Error in MOD_INPUT; Incompatible input data.'
           print*,' This Equation Set is not supported: ',eqn_set
           stop
        end if
 
     end if
 
     !If periodicity is used, the GEOM1 should NOT be used for domain decomposition.
     !Only metis will work with periodic b.c.:
     if (irank == 0) then
        if(  x_boundary(1) == 3 .or.  x_boundary(2) == 3 .or. &
             y_boundary(1) == 3 .or.  y_boundary(2) == 3) then
 
           if(decomp_type(1:4) == 'geom') then
              print*,'!!!-----------------------------'
              print*,'!!! INCOMPATIBLE B.C. and Domain decomposition method!'
              print*,'    Periodic b.c. do NOT work with GEOM(1,2,...) decomposition method'
              print*,'    The decomposition method will be changed to metis2d'
              decomp_type = 'metis2d'
              print*,'    Now: decomp_type = ', decomp_type
              print*,'!!!-----------------------------'
           end if
 
           if(decomp_type(1:7) == 'metis3d') then
              print*,'!!!-----------------------------'
              print*,'!!! INCOMPATIBLE B.C. and Domain decomposition method!'
              print*,'    Periodic b.c. do NOT work with metis3d decomposition method'
              print*,'    The decomposition method will be changed to metis2d'
              decomp_type = 'metis2d'
              print*,'    Now: decomp_type = ', decomp_type
              print*,'!!!-----------------------------'
           end if
 
        end if
     end if
 
     !Run checks on the time averaging parameters:
     if(lwrite_time_averaging) then
        if(interval_start < 0.0 .or. interval_end < 0.0) then
           lwrite_time_averaging = .false.
 
           if(irank == 0) then
              print*, '---------------------------------------------------'
              print*,' INPUT WARNING:'
              print*,' lwrite_time_averaging = TRUE'
              print*,' however, interval_start or interval_start are < 0'
              print*,' CAREFUL: interval_start and interval_end are lengths and canNOT be negative'
              print*,' Please, correct your input'
              print*,' NO time Averaging will be performed'
              print*, '---------------------------------------------------'
           end if
        end if
 
        dt_meas = interval_end - interval_start
        if(dt > dt_meas) lwrite_time_averaging = .false.
 
     end if
 
     !
     !lviscous_bottom=lviscous_ground
     !just in case the user uses the keyword 'lviscous_bottom' instead of 'lviscous_ground'
     lviscous_bottom=lviscous_ground
 
     !Grid Generation and Graph Partitioning
     if (lp4est .or. lp6est) then
        lserial_grid_creation=.false.
        fname_root=trim(fname_root) // '_P4est'
     else
        lp4est=.false.
        lp6est=.false.
     end if
 
     if(lp6est) then
        if (geometry_type/='cube') then
           nelx=nel_root_h
           nely=nel_root_h
           nelz=nel_root_v
        endif
     end if
 
     !
     ! What viscosity are we using for dynamics and tracers?
     ! here you chose:
     !1) Dynamics:
     if(lnazarov) then
        visc_dyn_flg = 0 !FXG 10/7/16: This flag is no longer used by any piece of code
     end if
     !2) Tracers
     if(lnazarov_tracers) then
        visc_tracers_flg = 0
        naza_tracers_flg = 1
     end if
 
     !Andreas' flags for Marcin's test case:
     if(lpassive) then
        lmoist_forcing = .false.
        if(icase==6001) lmoist_column  = .true. !Steve Guimond's test case
     end if
 
     !Limiter flags:
     if(llimit) then
        llimit_below = .false.
        llimit_above = .false.
     end if
 
     if(llimit_below .or. llimit_above) llimit = .false.
 
     read_external_grid_flg = 0
     if(lread_external_grid) read_external_grid_flg = 1
 
 
     !Flags for APPLY_BOUNDARY_CONDITIONS
     if (lincompressible .or. cgdg_method == 'separate' .or. is_sphere .or. equations == 'shallow') then
        lapply_bcs_velocity=.true.
     end if
 
     !1D-IMEX Needs to APPLY BCs Strongly: otherwise, not stable
     if (si_dimension == '1d') then
        lapply_bcs_velocity=.true.
     end if
     
     if (x_boundary(1) == 3 .or. x_boundary(2) == 3 .or. x_periodic > 0) then
        x_periodic    = 1
        x_boundary(1) = 3
        x_boundary(2) = 3
     endif
 
     if (y_boundary(1) == 3 .or. y_boundary(2) == 3 .or. y_periodic > 0) then
        y_periodic    = 1
        y_boundary(1) = 3
        y_boundary(2) = 3
     endif
 
     if (z_boundary(1) == 3 .or. z_boundary(2) == 3 .or. z_periodic > 0) then
        z_periodic    = 1
        z_boundary(1) = 3
        z_boundary(2) = 3
     endif
 
     if (space_method == 'cgc' .and. (x_periodic+y_periodic+z_periodic > 0)) then
        if (irank == 0) then
           print*, '---------------------------------------------------'
           print*,' INPUT WARNING:'
           print*,' Space_Method and Periodic BCs conflict!'
           print*,' space_method = ',space_method
           print*,' x_periodic y_periodic z_periodic = ',x_periodic,y_periodic,z_periodic
           print*,' Please, correct your input and rerun'
           print*, '---------------------------------------------------'
        end if
        stop
     end if
 
     !Inundation Checks
     if (equations(1:7) == 'shallow') is_shallow=.true.
     if (.not. llimit) limit_threshold=-1e8 !large number will keep checks in the code from changing the water height.
 
     !check if Shallow Water tests are 1D
     if (is_shallow) then
        select case(icase)
        case(25,26,27)
           is_1d=.true.
        end select
     end if
 
     !Implicit Flag
     if (ti_method == 'rk') delta=-1.0
     is_implicit = (int(delta) >= 0 .and. ldynamics)
     is_imex = (is_implicit .and. ti_method(6:9) /= 'jfnk')
 
     !Checks on Space_Method and Equation_Sets for Implicit
     if (.not. lincompressible) then
 !!$       if (space_method /= 'cgc' .and. eqn_set == 'set2nc' .and. is_imex) then
 !!$          if (irank == 0) then
 !!$             print*, '---------------------------------------------------'
 !!$             print*,' INPUT WARNING:'
 !!$             print*,' Space_Method and Eqn_Set conflict for IMEX!'
 !!$             print*,' space_method = ',space_method
 !!$             print*,' eqn_set = ',eqn_set
 !!$             print*,' Please, correct your input and rerun'
 !!$             print*, '---------------------------------------------------'
 !!$          end if
 !!$          stop
 !!$       end if
 !!$       if (space_method == 'cgc' .and. cgdg_method== 'unified' .and. is_imex) then
 !!$          if (irank == 0) then
 !!$             print*, '---------------------------------------------------'
 !!$             print*,' INPUT WARNING:'
 !!$             print*,' CGC-Unified currently does not work with IMEX'
 !!$             print*, '---------------------------------------------------'
 !!$          end if
 !!$          stop
 !!$       end if
        if ( (space_method == 'cgd' .and. cgdg_method == 'unified') .or. space_method == 'dg') then
           if (is_implicit .and. si_method == 'schur') then
              if (irank == 0) then
                 print*, '---------------------------------------------------'
                 print*,' INPUT WARNING:'
                 print*,' CGD-Unified and DG do not currently work with Schur form!'
                 print*, '---------------------------------------------------'
              end if
              stop
           end if
        end if
     end if
 
 !!$    !SET2NC does not work for CGD/DG and IMEX (CGD explicit should work)
 !!$    if (.not. lincompressible) then
 !!$       if (space_method /= 'cgc' .and. eqn_set == 'set2nc') then
 !!$          if (si_dimension /= '0d') then
 !!$             if (irank == 0) then
 !!$                print*, '---------------------------------------------------'
 !!$                print*,' INPUT WARNING:'
 !!$                print*,' Space_Method, Equation Sets, and IMEX solver conflict!'
 !!$                print*,' lincompressible = ',lincompressible
 !!$                print*,' space_method = ',space_method
 !!$                print*,' eqn_set = ',eqn_set
 !!$                print*,' si_dimension = ',si_dimension
 !!$                print*,' Please, correct your input and rerun'
 !!$                print*, '---------------------------------------------------'
 !!$             end if
 !!$             stop
 !!$          end if
 !!$       end if
 !!$    end if
 
     !Checks on Viscosity for Incompressible
     if (lincompressible) then
        if (icase == 9001 .or. icase == 9002) lvisc_reciprocal=.true.
     end if
 
     !1D-IMEX and P6est
     if (si_dimension == '1d') then
        if (lp4est) then
           if(irank == 0) then
              print*, '---------------------------------------------------'
              print*,' INPUT WARNING:'
              print*,' 1D-IMEX cannot be used with P4est. Please switch flag to P6est.'
              print*,' Please, correct your input and rerun'
              print*, '---------------------------------------------------'
           end if
           stop
        end if
     end if
 
     !Viscosity Flags
     !!    if (lvisc_anisotropic) vertical_viscosity=0
 
     !Checks on Space_Method and IMEX Forms
     !    if (delta >= 0 .and. ti_method(1:2) /= 'rk' .and. si_method == 'schur') then
     !       if (space_method == 'dg' .or. (space_method == 'cgd' .and. cgdg_method /= 'separate') ) then
 !!!       icheck = 0
 !!!       if (space_method == 'dg') then
 !!!          icheck = 1
 !!!       else if (space_method == 'cgd') then
 !!!          if (eqn_set == 'set2c' .or. eqn_set == 'set3c') icheck = 1
 !!!       end if
 !!!       if (icheck == 1) then
     !
     !          if(irank == 0) then
     !             print*, '---------------------------------------------------'
     !             print*,' INPUT WARNING:'
     !             print*,' Space_Method and IMEX Forms conflict!'
     !             print*,' space_method = ',space_method
     !             print*,' si_method = ',si_method
     !             print*,' cgdg_method = ',cgdg_method
     !             print*,' Please, correct your input and rerun'
     !             print*, '---------------------------------------------------'
     !          end if
     !          stop
     !       end if
     !    end if
 
   end subroutine mod_input_create
 
 end module mod_input
 