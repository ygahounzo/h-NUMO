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
!
!>
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: April 03, 2024
!----------------------------------------------------------------------!
module mod_input

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
        lapply_bcs_velocity, &
        lLAV, &
        llimit, llimit_below, llimit_above, limiter_type, limiter_qoi, &
        ltime_residual, time_residual_flg, &
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
        ladapt_timestep, lprint_diagnostics, iprint_diagnostics, lcheck_all_conserved, lcompute_barycenter, ljenkins_test, &
        ldynamics, &
        ctd_file_temp, ctd_file_salinity, &
        lrotating_flow, &
        Mach1_in, Mach2_in, &
        bcast_type, &
        lout_spherical_shell, &
        zlevel_out, &
        lout_nc_3d, &
        p00_in, &
        lpert_theta, lpert_qv, lpert_qc, lpert_qr, lcylinder, &
        ldam, &
        lpert_velo, &
        ljet, &
        thetac_in, qvc_in, qcc_in, qrc_in, &
        thetac2_in, qvc2_in, qcc2_in, qrc2_in, &
        thetac3_in, qvc3_in, qcc3_in, qrc3_in, &
        thetac4_in, qvc4_in, qcc4_in, qrc4_in, &
        steepness_pert_in, &
        space_method, &
        dump_rhs, &
        n_corrections, &
        imass, &
        ad_mlswe, explt_coriolis, cd_mlswe, dp_cutoff1, dp_cutoff2, &
        dp_tau_bot, dp_tau_wind, dt_btp, method_visc,visc_mlswe,dpprime_visc_min, max_shear_dz, matlab_viz, &
        adjust_H_vertical_sum, is_mlswe_linear, botfr, mass_exact, bcl_flux, mlswe_bc_strong, dg_integ_exact, dump_data, flux_type
 
   public :: nsmall, ti_alpha, ti_beta, dd_alpha
 
   public :: eqn_set, equations, lincompressible, piso_type, lpassive_advection, lhelmholtz, lALE, &
        lset_outlet_pressure
 
   public :: is_shallow, is_sphere, is_implicit, is_1d, is_imex, is_swe_layers, nonlinear_swe, is_mlswe
 
   public :: nelx, nely, nelz, nopx, nopy, nopz, xdims, ydims, ztop, zbottom, &
        nlayers, & !shallow water layers
        delta_domain, & !used for turbulent channel flows
        nlon, nlat, &
        x_boundary, y_boundary, z_boundary, &
        x_periodic, y_periodic, z_periodic, &
        robin_bc_alpha, &
        lfree_slip_exact, &
        geometry_type, nproc_z, &
        bc_tscale, bc_xscale, bc_yscale, bc_zscale, &
        sponge_type, sponge_top_coe, sponge_lateralx_coe, sponge_lateralx_coe_east, sponge_lateralx_coe_west, sponge_lateraly_coe, &
        lsommerfeld, &
        lgrid_only, &
        lp4est, lp6est, lserial_grid_creation, lparallel_grid_creation, lio_grid_ascii, lwrite_grid_ascii, fname_initial, &
        restart_path, & ! Added by Yao G.
        nel_root_h, refinement_levels_h, nel_root_v, refinement_levels_v, &
        xstretch_coe, ystretch_coe, zstretch_coe, &
        lxstretch, lystretch, lzstretch, &
        xstretch_type, ystretch_type, zstretch_type, &
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
   character     :: ctd_file_temp*124, ctd_file_salinity*124
   logical       :: lprint_diagnostics
   integer       :: iprint_diagnostics = 1 !every how many iterations to print diagnostics
   logical       :: lcheck_all_conserved = .false.
   logical       :: ljenkins_test = .true.
 
   character     :: solver_type*9, real_string*9
   logical       :: lcg_proj = .false.
   character     :: Iter_Type*5, precon_type*1, precon_mode*5, precon_fname*150
   character     :: space_method*3
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
   real(kind=r8)         :: vertical_viscosity = 1.0 !1=On and 0=Off
   integer               :: filter_tracers_flg = 0 !This flag activates the filter computation for the tracers equations (INactive by default)
   
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
 
   integer               :: steepness_pert_in = 0 ! 0: sin-shaped rtb, >0: shape of rtb given by exp(-r**steepness_pert_in)
   character(len=12)     :: zgrid_type = 'sigma'
   real(kind=r8)         :: time_scale = 1.0_r8
   real(kind=r8)         :: Mach1_in    = 0.0
   real(kind=r8)         :: Mach2_in    = 0.1
   real(kind=r8)         :: hybrid_s    = 1000.0
   real(kind=r8)         :: sleve_s1    = 1000.0
   real(kind=r8)         :: sleve_s2    = 1000.0
   real(kind=r8)         :: p00_in = -999.999
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
   real(kind=r8), dimension(6) :: dirichlet_uvelo !dirichlet_uvelo contains a value of uvelo per face (6 faces in a box)
   real(kind=r8), dimension(6) :: dirichlet_vvelo !dirichlet_vvelo contains a value of vvelo per face (6 faces in a box)
   real(kind=r8), dimension(6) :: dirichlet_wvelo !dirichlet_wvelo contains a value of wvelo per face (6 faces in a box)

 
   real(kind=r8)         :: c_am2 = 3.0/2.0 !AM2 constant from Durran-Blossey (2nd order)
 
   real(kind=r8)         :: kessler_production_constant = 0.1
 
   character(len=4)      :: bcast_type = 'mpi'
   character(len=12)     :: tau_type   = 'codina'
   character(len=12)     :: shock_type = 'codina'
   character(len=12)     :: eqn_set    = 'set2nc'
   character(len=12)     :: equations  = 'euler'
   character(len=12)     :: limiter_type = 'none'
   character(len=3)      :: piso_type = 'dsa'
 
   logical :: ldynamics      = .true.
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
   logical :: lstabilize_mass_eqn = .false.
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
   character(len=12) :: flux_type = 'centered'

 
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
          llimit, llimit_below, llimit_above, limiter_type, limiter_qoi, &
          ltime_residual, &
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
          ladapt_timestep,lprint_diagnostics,iprint_diagnostics,lcheck_all_conserved,lcompute_barycenter,ljenkins_test, &
          ldynamics, &
          ctd_file_temp, ctd_file_salinity, &
          Mach1_in, Mach2_in, &
          lrotating_flow, &
          bcast_type, &
          p00_in, &
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
          space_method, dump_rhs, &
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
          mass_exact, bcl_flux, mlswe_bc_strong, dg_integ_exact, dump_data, flux_type
 
     namelist /gridnl/ nelx, nely, nelz, nopx, nopy, nopz, xdims, ydims, ztop, zbottom, &
          nlayers, &
          delta_domain, & !used for turbulent channel flows
          nlon, nlat, &
          geometry_type, nproc_z, &
          x_boundary, y_boundary, z_boundary, &
          x_periodic, y_periodic, z_periodic, &
          robin_bc_alpha, &
          lset_outlet_pressure, &
          lfree_slip_exact, &
          bc_tscale, bc_xscale, bc_yscale, bc_zscale, &
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
 
     geometry_type = lowercase(geometry_type)
     ti_method = lowercase(ti_method)
     ti_method_btp = lowercase(ti_method_btp) !added by YaoG
     si_method = lowercase(si_method)
     si_dimension = lowercase(si_dimension)
     filter_weight_type = lowercase(filter_weight_type)
     filter_basis_type = lowercase(filter_basis_type)
     out_type = lowercase(out_type)
     format_vtk = uppercase(format_vtk)
     vtk_cell_type = uppercase(vtk_cell_type)
 
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

      if (ti_method_btp /= 'rk35' .and. ti_method_btp /= 'rk34' .and. ti_method_btp /= 'lsrk') then
         ti_method_btp = 'explicit'
      endif
      fname_root=trim(fname_root) // '_' // trim(ti_method) // '_' // trim(ti_method_btp)

 
     ! Perform some checks on the Initial Data
     !
     if (irank == 0) then
 
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

 
     read_external_grid_flg = 0
     if(lread_external_grid) read_external_grid_flg = 1
 
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
 
   end subroutine mod_input_create
 
 end module mod_input
 