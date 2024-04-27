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
        time_initial, time_final, time_dynamic_amr, time_restart, time_scale, irestart_file_number,  &
        icase, ti_method, si_method, si_dimension, &
        ti_method_btp, & ! Added by YaoG
        filter_mux, filter_muy, filter_muz, &
        ifilter, kstages, &
        filter_weight_type, filter_basis_type, &
        fname_root, out_type, lout_ascii, lout_asciimaya, format_vtk, nvtk_files, vtk_cell_type, &
        write_mesh, &
        filter_tracers_flg, &
        ladapt_timestep, lprint_diagnostics, iprint_diagnostics, &
        bcast_type, &
        space_method, &
        imass, &
        ad_mlswe, cd_mlswe, dp_cutoff1, dp_cutoff2, &
        dp_tau_bot, dp_tau_wind, dt_btp, method_visc,visc_mlswe,dpprime_visc_min, max_shear_dz, matlab_viz, &
        adjust_H_vertical_sum, is_mlswe_linear, botfr, mass_exact, bcl_flux, mlswe_bc_strong, dg_integ_exact, dump_data, flux_type
 
   public :: eqn_set, is_mlswe
 
   public :: nelx, nely, nelz, nopx, nopy, nopz, xdims, ydims, ztop, zbottom, &
        nlayers, & !shallow water layers
        x_boundary, y_boundary, z_boundary, &
        x_periodic, y_periodic, z_periodic, &
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
        lread_external_grid, read_external_grid_flg, is_non_conforming_flg, &
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
   public :: gravity_in, &
        mesh_file, limit_threshold
         
   private
 
   !-----------------------------------------------------------------------
   ! Namelist Variables
   !-----------------------------------------------------------------------
   real(kind=r8) :: dt, dt0, dt1, dt2, dt_btp
   real(kind=r8) :: restoring_time = 1000000000
   logical       :: lrestoring_sponge = .false.
   real(kind=r8) :: time_initial, time_final, time_restart, lf2_filter, lf2_raw_filter, lf2_alpha, lf2_nu

   real(kind=r8) :: time_dynamic_amr=0.0
   real(kind=r8) :: interval_end=0.0
   real(kind=r8) :: filter_mux, filter_muy, filter_muz
   integer       :: icase, kstages, ifilter
   integer       :: irestart_file_number = 0
   character     :: filter_basis_type*8, filter_weight_type*4
   character     :: fname_root*150, out_type*5, fname_initial*100
   character     :: restart_path*100
   logical       :: write_mesh=.false.
   character     :: ti_method*14, si_method*8, si_dimension*2
   character     :: ti_method_btp*14
   logical       :: lprint_diagnostics
   integer       :: iprint_diagnostics = 1 !every how many iterations to print diagnostics
 
   character     :: space_method*3, real_string*9
 
 
   real(kind=r8) :: ad_mlswe = 0.0
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
 
   integer               :: imass = 1

   real(kind=r8)         :: time_scale = 1.0_r8

 
   character(len=4)      :: bcast_type = 'mpi'
   character(len=12)     :: eqn_set    = 'set2nc'
 
   logical :: lrestart_file  = .false. !obsolete
   logical :: lsommerfeld    = .false.
   logical :: lxstretch      = .false.
   logical :: lystretch      = .false.
   logical :: lzstretch      = .false.
   logical :: lgrid_only     = .false.
   logical :: lp4est         = .false.
   logical :: lp6est         = .false.
   logical :: lread_external_grid = .false.
   logical :: lread_bc       = .false.
   logical :: lserial_grid_creation = .true.
   logical :: lparallel_grid_creation = .false.
   logical :: lio_grid_ascii = .false.
   logical :: lwrite_grid_ascii = .false.
   logical :: lout_ascii     = .false.
   logical :: lout_asciimaya = .false.
   logical :: lcompute_barycenter = .false.
   logical :: luse_min_element_length = .false.
   logical :: lrotating_flow = .false.
 
   logical :: ladapt_timestep = .false.
 
   !Input parameters to write a spherical shell (1 level) to a netcdf file:
   logical :: lgpu = .false.
   logical :: luse_hybrid_cpu_gpu = .false.

   logical :: is_sphere = .false. 
   
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
 
   real time_scale_in
 
   real dam_depth_in
   real, dimension(2) :: dam_xlimit_in, dam_ylimit_in
 
   real bathymetry_shift
   integer hump_config
   integer ibathymetry
   integer nonlinear, balance_flag, output_flag, warp_grid
   character mesh_file*100, bathymetry_file*100
 
   logical lout_tree
   logical lout_shoreline
   logical lread_external_bathy
   real gravity_in
   real xdam_min, xdam_max
   real ydam_min, ydam_max
 
   real :: limit_threshold = 1e-3
 
   logical llinear_pert
 
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
 
     namelist /eqnset/ eqn_set, is_mlswe
 
     namelist /input/ dt, &
         restoring_time, &
         lrestoring_sponge, &
         time_initial, time_final, time_dynamic_amr, time_restart, time_scale, irestart_file_number,&
         lrestart_file, icase, &
         ti_method, si_method, si_dimension, &
         ti_method_btp, & ! Added by Yao G.
         kstages, &
         filter_mux, filter_muy, filter_muz, ifilter, &
         filter_weight_type, filter_basis_type, fname_root, out_type, lout_ascii, lout_asciimaya, format_vtk, nvtk_files, vtk_cell_type, &
         write_mesh, &
         fname_initial, &
         restart_path, & ! Added by Yao G.
         filter_tracers_flg, &
         ladapt_timestep,lprint_diagnostics,iprint_diagnostics, &
         bcast_type, &
         space_method, &
         lgpu, numaocca_dir, Nelems, Nslices, NslicesV, vectorization, &
         platform, platformID, deviceID, platformWeight, platform2, platformID2, deviceID2, platformWeight2, &
         cpus_per_node, gpus_per_node, threads_per_process, luse_hybrid_cpu_gpu, &
                              !!--Shallow
         ibathymetry, bathymetry_file, &
         gravity_in, limit_threshold, &
         mesh_file, &
         lread_external_bathy, &
         limit_threshold, &
         imass, &
         ad_mlswe, cd_mlswe, dp_cutoff1, dp_cutoff2, dp_tau_bot, dp_tau_wind, dt_btp,method_visc,&
         visc_mlswe, dpprime_visc_min, max_shear_dz, matlab_viz, adjust_H_vertical_sum, is_mlswe_linear, botfr, &
         mass_exact, bcl_flux, mlswe_bc_strong, dg_integ_exact, dump_data, flux_type
 
     namelist /gridnl/ nelx, nely, nelz, nopx, nopy, nopz, xdims, ydims, ztop, zbottom, &
          nlayers, &
          geometry_type, nproc_z, &
          x_boundary, y_boundary, z_boundary, &
          x_periodic, y_periodic, z_periodic, &
          bc_tscale, bc_xscale, bc_yscale, bc_zscale, &
          sponge_type, sponge_top_coe, sponge_lateralx_coe,sponge_lateralx_coe_east, sponge_lateralx_coe_west, sponge_lateraly_coe, &
          lsommerfeld, &
          lgrid_only, &
          lp4est, lp6est, lio_grid_ascii, &
          lread_external_grid, &
          is_non_conforming_flg, &
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
          lxstretch, lystretch, lzstretch
 
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
 
   end subroutine mod_input_create
 
 end module mod_input
 