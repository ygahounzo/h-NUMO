!----------------------------------------------------------------------!
!>@brief This module defines the Constants that were in PARAM.H
!>@author F.X. Giraldo and J.F. Kelly
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943-5216
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

   type input
      !-----------------------------------------------------------------------
      ! Namelist Variables
      !-----------------------------------------------------------------------
      real(kind=r8) :: dt, dt0, dt1, dt2, dt_btp
      real(kind=r8) :: restoring_time
      logical       :: lrestoring_sponge
      real(kind=r8) :: time_initial, time_final, time_restart
      real(kind=r8) :: lf2_filter, lf2_raw_filter, lf2_alpha, lf2_nu

      real(kind=r8) :: time_dynamic_amr
      real(kind=r8) :: interval_end
      real(kind=r8) :: filter_mux, filter_muy, filter_muz
      integer       :: kstages, kstages_bcl, ifilter
      integer       :: irestart_file_number
      character(len=8)   :: filter_basis_type
      character(len=4)   :: filter_weight_type
      character(len=150) :: fname_root
      character(len=5)   :: out_type
      character(len=100) :: fname_initial
      character(len=100) :: restart_path
      logical            :: write_mesh
      character(len=14)  :: ti_method_btp
      character(len=20)  :: test_case
      logical            :: lprint_diagnostics
      integer            :: iprint_diagnostics  ! every how many iterations to print diagnostics

      character(len=3)   :: space_method
      character(len=9)   :: real_string

      real(kind=r8) :: ad_mlswe
      real(kind=r8) :: cd_mlswe
      real(kind=r8) :: dp_tau_bot
      real(kind=r8) :: dp_tau_wind
      integer       :: method_visc
      integer       :: botfr
      integer       :: adjust_H_vertical_sum
      integer       :: adjust_bcl_mom_flux
      real(kind=r8) :: visc_mlswe
      real(kind=r8) :: max_shear_dz
      real(kind=r8) :: f0
      real(kind=r8) :: beta
      real(kind=r8) :: dry_cutoff

      !-----------------------------------------------------------------------
      ! Namelist Variables
      !-----------------------------------------------------------------------
      real(kind=r8), dimension(2) :: xdims, ydims
      real(kind=r8)     :: ztop
      real(kind=r8)     :: pi_trig
      real(kind=r8)     :: zbottom
      real(kind=r8)     :: delta_domain  ! used for turbulent channel flows (default value is 1.0)
      real(kind=r8)     :: xstretch_coe, ystretch_coe, zstretch_coe
      integer           :: nelx, nely, nopx, nopy
      integer           :: nopz   ! 0th order polynomial in the z direction
      integer           :: nelz   ! number of elements in the z direction
      integer           :: nlayers ! number of shallow water layers
      integer           :: nproc_z
      integer           :: zlevel_out  ! spherical level to write to netcdf (default 1)
      character(len=20) :: geometry_type
      integer           :: nvar_bcl, nvar_btp  ! number of baroclinic/barotropic prognostic variables, derived from geometry_type
      integer           :: ngrd_var            ! number of spatial gradient directions, derived from geometry_type
      integer           :: ngraduvw_var        ! width of the combined btp velocity-gradient array (4 cartesian, 9 sphere_hex)
      character(len=24) :: sponge_type
      character(len=8)  :: format_vtk      ! ascii or binary
      character(len=8)  :: vtk_cell_type   ! standard or lagrange
      integer           :: nvtk_files  ! number of files written in parallel (must be < nproc)

      character(len=24)  :: xstretch_type
      character(len=24)  :: ystretch_type
      character(len=24)  :: zstretch_type
      character(len=24)  :: platform
      character(len=24)  :: platform2
      character(len=6)   :: vectorization
      character(len=150) :: numaocca_dir

      character(len=512) :: OpenCLCompilerFlags
      character(len=512) :: CudaCompilerFlags
      character(len=512) :: OpenMPCompilerFlags
      character(len=512) :: SerialCompilerFlags
      character(len=512) :: PthreadsCompilerFlags
      character(len=512) :: COICompilerFlags

      !-----------------------------------------------------------------------
      ! Default Namelist Values
      !-----------------------------------------------------------------------
      real(kind=r8) :: bc_tscale
      real(kind=r8) :: bc_xscale
      real(kind=r8) :: bc_yscale
      real(kind=r8) :: bc_zscale
      real(kind=r8) :: sponge_top_coe
      real(kind=r8) :: sponge_lateralx_coe
      real(kind=r8) :: sponge_lateralx_coe_east
      real(kind=r8) :: sponge_lateralx_coe_west
      real(kind=r8) :: sponge_lateraly_coe

      real(kind=r8) :: SIPG_constant

      integer, dimension(2) :: x_boundary
      integer, dimension(2) :: y_boundary
      integer, dimension(2) :: z_boundary
      integer               :: x_periodic
      integer               :: y_periodic
      integer               :: z_periodic
      real(kind=r8)         :: robin_bc_alpha  ! alpha in: n.grad(q) - alpha*q = beta Robin BC
      logical               :: lfree_slip_exact
      integer               :: refinement_levels_h
      integer               :: read_external_grid_flg
      integer               :: full_stress_flg
      integer               :: is_non_conforming_flg
      logical               :: interp_cg_flux_flg
      integer               :: max_mesh_lvl
      logical               :: nc_box_invert
      integer               :: p4est_log_level
      real                  :: xlim_min, xlim_max
      real                  :: ylim_min, ylim_max
      real                  :: zlim_min, zlim_max
      integer, dimension(10) :: amr_indicator_variables
      real, dimension(2)     :: amr_smoothness_limits
      real, dimension(2)     :: amr_max_min_lim
      real, dimension(10)    :: amr_threshold_lim
      real, dimension(10)    :: amr_smoothness_qL2_limit
      logical                :: amr_mark_modes
      logical                :: amr_mark_modes_use_baseline_decay
      logical                :: amr_mark_max_min
      logical                :: amr_mark_random
      logical                :: amr_mark_threshold
      integer(C_INT32_T)     :: amr_num_neigh_iter
      logical                :: amr_mark_set2nc

      integer :: nel_root_h
      integer :: filter_tracers_flg  ! activates filter for tracer equations (inactive by default)

      integer :: platformID,  deviceID
      integer :: platformID2, deviceID2
      integer :: Nelems
      integer :: Nslices
      integer :: NslicesV
      integer :: cpus_per_node
      integer :: gpus_per_node
      integer :: threads_per_process
      integer :: platformWeight
      integer :: platformWeight2

      integer :: imass

      real(kind=r8) :: time_scale

      character(len=4)  :: bcast_type
      character(len=12) :: eqn_set

      logical :: lrestart_file  ! obsolete
      logical :: lsommerfeld
      logical :: lxstretch
      logical :: lystretch
      logical :: lzstretch
      logical :: lgrid_only
      logical :: lread_external_grid
      logical :: lread_bc
      logical :: lserial_grid_creation
      logical :: lparallel_grid_creation
      logical :: lwrite_grid_ascii
      logical :: lout_ascii
      logical :: lout_asciimaya
      logical :: lcompute_barycenter
      logical :: luse_min_element_length
      logical :: lrotating_flow

      logical :: ladapt_timestep

      logical :: lgpu
      logical :: luse_hybrid_cpu_gpu

      logical :: is_mlswe
      logical :: dg_integ_exact
      logical :: dump_data
      logical :: lcheck_conserved

      real :: time_scale_in
      real :: dam_depth_in
      real, dimension(2) :: dam_xlimit_in, dam_ylimit_in

      real    :: bathymetry_shift
      integer :: hump_config
      integer :: ibathymetry
      integer :: nonlinear, balance_flag, output_flag, warp_grid
      character(len=100) :: mesh_file, bathymetry_file

      character(len=20) :: bcl_time_method
      logical :: rk_bcl_FS

      logical :: lout_tree
      logical :: lout_shoreline
      logical :: lread_external_bathy
      real    :: gravity_in
      real    :: xdam_min, xdam_max
      real    :: ydam_min, ydam_max

      real :: limit_threshold

      logical :: llinear_pert
   
   end type input
 
   public :: mod_input_create, input
         
   private
 
   !-----------------------------------------------------------------------
   ! Namelist Variables
   !-----------------------------------------------------------------------
   
 
   !-----------------------------------
   ! End shallow water variables 
   !-----------------------------------
 
 contains
 
   !-----------------------------------------------------------------------
   subroutine mod_input_create(inp)
 
     implicit none

     type(input), intent(inout) :: inp
 
     integer icheck
 
     integer :: funit

     real(kind=r8) :: dt, dt0, dt1, dt2, dt_btp
      real(kind=r8) :: restoring_time = 1000000000
      logical       :: lrestoring_sponge = .false.
      real(kind=r8) :: time_initial, time_final, time_restart, lf2_filter, lf2_raw_filter, lf2_alpha, lf2_nu

      real(kind=r8) :: time_dynamic_amr=0.0
      real(kind=r8) :: interval_end=0.0
      real(kind=r8) :: filter_mux, filter_muy, filter_muz
      integer       :: kstages, kstages_bcl = 3, ifilter
      integer       :: irestart_file_number = 0
      character     :: filter_basis_type*8, filter_weight_type*4
      character     :: fname_root*150, out_type*5, fname_initial*100
      character     :: restart_path*100
      logical       :: write_mesh=.false.
      character     :: ti_method_btp*14, test_case*20
      logical       :: lprint_diagnostics
      integer       :: iprint_diagnostics = 1 !every how many iterations to print diagnostics
   
      character     :: space_method*3, real_string*9
   
   
      real(kind=r8) :: ad_mlswe = 0.0
      real(kind=r8) :: cd_mlswe = 0.0
      real(kind=r8) :: dp_tau_bot = 0.0
      real(kind=r8) :: dp_tau_wind = 0.0
      integer :: method_visc = 0
      integer :: botfr = 0
      integer :: adjust_H_vertical_sum = 2
      integer :: adjust_bcl_mom_flux = 1
      real(kind=r8) :: visc_mlswe = 0.0
      real(kind=r8) :: max_shear_dz = 0.0
      real(kind=r8) :: f0 = 0.0
      real(kind=r8) :: beta = 0.0
      real(kind=r8) :: dry_cutoff = 1.0e-10
   
      !-----------------------------------------------------------------------
      ! Namelist Variables
      !-----------------------------------------------------------------------
      real(kind=r8), dimension(2) :: xdims, ydims
      real(kind=r8)     :: ztop = 0.0
      real(kind=r8)     :: pi_trig
      real(kind=r8)     :: zbottom = 0.0
      real(kind=r8)     :: delta_domain = 1.0 !used for turbulent channel flows (default value is 1.0)
      real(kind=r8)     :: xstretch_coe, ystretch_coe, zstretch_coe
      integer           :: nelx, nely, nopx, nopy
      integer           :: nopz = 0.0 ! 0th order poynomial in the z direction
      integer           :: nelz = 1 !number of elements in the z direction
      integer           :: nlayers = 1 !number of shallow water layers
      integer           :: nproc_z
      integer           :: zlevel_out = 1 !This integer is given by the user to decide what spherical level to write to a netcdf file. Value 1 is by default
      character(len=20) :: geometry_type = 'cartesian'
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
      integer, dimension(2)  :: z_boundary = 0
      integer                :: x_periodic = 0
      integer                :: y_periodic = 0
      integer                :: z_periodic = 0
      real(kind=r8)          :: robin_bc_alpha = 1.0 ! sets value for alpha in n \cdot \nabla q - \alpha q = \beta Robin BC
      logical                :: lfree_slip_exact = .false.
      integer                :: refinement_levels_h = 0
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
      logical :: lread_external_grid = .false.
      logical :: lread_bc       = .false.
      logical :: lserial_grid_creation = .true.
      logical :: lparallel_grid_creation = .false.
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
      
      logical :: is_mlswe = .true.  
      logical :: dg_integ_exact = .true. 
      logical :: dump_data = .true. 
      logical :: lcheck_conserved = .false. 

   
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

      character(len=20) :: bcl_time_method = '2levels'
      logical :: rk_bcl_FS = .false.

     !Namelist Input
 
     namelist /eqnset/ eqn_set, is_mlswe
 
     namelist /input_nml/ dt, &
         restoring_time, &
         lrestoring_sponge, &
         time_initial, time_final, time_dynamic_amr, time_restart, time_scale, irestart_file_number,&
         lrestart_file, test_case, &
         ti_method_btp, &
         kstages, &
         kstages_bcl, &
         filter_mux, filter_muy, filter_muz, ifilter, &
         filter_weight_type, filter_basis_type, fname_root, out_type, lout_ascii, lout_asciimaya, format_vtk, nvtk_files, vtk_cell_type, &
         write_mesh, &
         fname_initial, &
         restart_path, & 
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
         ad_mlswe, cd_mlswe, dp_tau_bot, dp_tau_wind, dt_btp,method_visc,&
         visc_mlswe, max_shear_dz, adjust_H_vertical_sum, botfr, &
         dg_integ_exact, dump_data, lcheck_conserved, adjust_bcl_mom_flux, &
         f0, beta, dry_cutoff, &
         bcl_time_method, rk_bcl_FS
 
     namelist /gridnl/ nelx, nely, nelz, nopx, nopy, nopz, xdims, ydims, ztop, zbottom, &
         geometry_type, &
         nlayers, &
         nproc_z, &
         x_boundary, y_boundary, z_boundary, &
         x_periodic, y_periodic, z_periodic, &
         bc_tscale, bc_xscale, bc_yscale, bc_zscale, &
         sponge_type, sponge_top_coe, sponge_lateralx_coe,sponge_lateralx_coe_east, sponge_lateralx_coe_west, sponge_lateraly_coe, &
         lsommerfeld, &
         lgrid_only, &
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
         lserial_grid_creation, lparallel_grid_creation, lwrite_grid_ascii,&
         refinement_levels_h, nel_root_h,&
         xstretch_coe, ystretch_coe, zstretch_coe, &
         lxstretch, lystretch, lzstretch
 
      amr_indicator_variables(1) = 1
   
      !Read input namelist
      funit = get_unit()
      print*,"File:",namelist_input, funit
      open(funit,file=namelist_input)
      read(funit,input_nml)
      close(funit)
   
      !Read gridnl namelist
      funit = get_unit()
      open(funit,file=namelist_input)
      read(funit,gridnl)
      close(funit)
   
      ti_method_btp = lowercase(ti_method_btp) !added by YaoG
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
         fname_root=trim(fname_root) // '_' // trim(space_method)
   
      !Add TI_METHOD to FNAME
         fname_root=trim(fname_root) // '_' // trim(ti_method_btp)
   
      !Grid Generation and Graph Partitioning
         lserial_grid_creation=.false.

   
      read_external_grid_flg = 0
      if(lread_external_grid) read_external_grid_flg = 1
      ! sphere_ico builds its coarse icosahedral mesh as an ABAQUS .inp
      ! file (create_grid_sphere_ico_inp) and always feeds it to p4est
      ! through the read_external_grid_flg path, regardless of the
      ! user's lread_external_grid setting.
      if (trim(geometry_type) == 'sphere_ico') read_external_grid_flg = 1
   
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

      inp%dt                             = dt
      inp%dt0                            = dt0
      inp%dt1                            = dt1
      inp%dt2                            = dt2
      inp%dt_btp                         = dt_btp
      inp%restoring_time                 = restoring_time
      inp%lrestoring_sponge              = lrestoring_sponge
      inp%time_initial                   = time_initial
      inp%time_final                     = time_final
      inp%time_restart                   = time_restart
      inp%lf2_filter                     = lf2_filter
      inp%lf2_raw_filter                 = lf2_raw_filter
      inp%lf2_alpha                      = lf2_alpha
      inp%lf2_nu                         = lf2_nu
      inp%time_dynamic_amr               = time_dynamic_amr
      inp%interval_end                   = interval_end
      inp%filter_mux                     = filter_mux
      inp%filter_muy                     = filter_muy
      inp%filter_muz                     = filter_muz
      inp%kstages                        = kstages
      inp%kstages_bcl                    = kstages_bcl
      inp%ifilter                        = ifilter
      inp%irestart_file_number           = irestart_file_number
      inp%filter_basis_type              = filter_basis_type
      inp%filter_weight_type             = filter_weight_type
      inp%fname_root                     = fname_root
      inp%out_type                       = out_type
      inp%fname_initial                  = fname_initial
      inp%restart_path                   = restart_path
      inp%write_mesh                     = write_mesh
      inp%ti_method_btp                  = ti_method_btp
      inp%test_case                      = test_case
      inp%lprint_diagnostics             = lprint_diagnostics
      inp%iprint_diagnostics             = iprint_diagnostics
      inp%space_method                   = space_method
      inp%real_string                    = real_string
      inp%ad_mlswe                       = ad_mlswe
      inp%cd_mlswe                       = cd_mlswe
      inp%dp_tau_bot                     = dp_tau_bot
      inp%dp_tau_wind                    = dp_tau_wind
      inp%method_visc                    = method_visc
      inp%botfr                          = botfr
      inp%adjust_H_vertical_sum          = adjust_H_vertical_sum
      inp%adjust_bcl_mom_flux            = adjust_bcl_mom_flux
      inp%visc_mlswe                     = visc_mlswe
      inp%max_shear_dz                   = max_shear_dz
      inp%f0                             = f0
      inp%beta                           = beta
      inp%dry_cutoff                     = dry_cutoff
      inp%xdims                          = xdims
      inp%ydims                          = ydims
      inp%ztop                           = ztop
      inp%pi_trig                        = pi_trig
      inp%zbottom                        = zbottom
      inp%delta_domain                   = delta_domain
      inp%xstretch_coe                   = xstretch_coe
      inp%ystretch_coe                   = ystretch_coe
      inp%zstretch_coe                   = zstretch_coe
      inp%nelx                           = nelx
      inp%nely                           = nely
      inp%nopx                           = nopx
      inp%nopy                           = nopy
      inp%nopz                           = nopz
      inp%nelz                           = nelz
      inp%nlayers                        = nlayers
      inp%nproc_z                        = nproc_z
      inp%zlevel_out                     = zlevel_out
      inp%geometry_type                  = geometry_type
      if (trim(inp%geometry_type) == 'sphere_hex' .or. trim(inp%geometry_type) == 'sphere_ico') then
         inp%nvar_bcl                    = 4
         inp%nvar_btp                    = 5
         inp%ngrd_var                    = 3 ! 3 gradient directions for spherical coordinates
         inp%ngraduvw_var                = 9 ! du_dx,du_dy,dv_dx,dv_dy,gradw(3),graduvz(2)
      else
         inp%nvar_bcl                    = 3
         inp%nvar_btp                    = 4
         inp%ngraduvw_var                = 4 ! du_dx,du_dy,dv_dx,dv_dy
         inp%ngrd_var                    = 2 ! 2 gradient directions (d/dx, d/dy) for cartesian coordinates
      end if
      inp%sponge_type                    = sponge_type
      inp%format_vtk                     = format_vtk
      inp%vtk_cell_type                  = vtk_cell_type
      inp%nvtk_files                     = nvtk_files
      inp%xstretch_type                  = xstretch_type
      inp%ystretch_type                  = ystretch_type
      inp%zstretch_type                  = zstretch_type
      inp%platform                       = platform
      inp%platform2                      = platform2
      inp%vectorization                  = vectorization
      inp%numaocca_dir                   = numaocca_dir
      inp%OpenCLCompilerFlags            = OpenCLCompilerFlags
      inp%CudaCompilerFlags              = CudaCompilerFlags
      inp%OpenMPCompilerFlags            = OpenMPCompilerFlags
      inp%SerialCompilerFlags            = SerialCompilerFlags
      inp%PthreadsCompilerFlags          = PthreadsCompilerFlags
      inp%COICompilerFlags               = COICompilerFlags
      inp%bc_tscale                      = bc_tscale
      inp%bc_xscale                      = bc_xscale
      inp%bc_yscale                      = bc_yscale
      inp%bc_zscale                      = bc_zscale
      inp%sponge_top_coe                 = sponge_top_coe
      inp%sponge_lateralx_coe            = sponge_lateralx_coe
      inp%sponge_lateralx_coe_east       = sponge_lateralx_coe_east
      inp%sponge_lateralx_coe_west       = sponge_lateralx_coe_west
      inp%sponge_lateraly_coe            = sponge_lateraly_coe
      inp%SIPG_constant                  = SIPG_constant
      inp%x_boundary                     = x_boundary
      inp%y_boundary                     = y_boundary
      inp%z_boundary                     = z_boundary
      inp%x_periodic                     = x_periodic
      inp%y_periodic                     = y_periodic
      inp%z_periodic                     = z_periodic
      inp%robin_bc_alpha                 = robin_bc_alpha
      inp%lfree_slip_exact               = lfree_slip_exact
      inp%refinement_levels_h            = refinement_levels_h
      inp%read_external_grid_flg         = read_external_grid_flg
      inp%full_stress_flg                = full_stress_flg
      inp%is_non_conforming_flg          = is_non_conforming_flg
      inp%interp_cg_flux_flg             = interp_cg_flux_flg
      inp%max_mesh_lvl                   = max_mesh_lvl
      inp%nc_box_invert                  = nc_box_invert
      inp%p4est_log_level                = p4est_log_level
      inp%xlim_min                       = xlim_min
      inp%xlim_max                       = xlim_max
      inp%ylim_min                       = ylim_min
      inp%ylim_max                       = ylim_max
      inp%zlim_min                       = zlim_min
      inp%zlim_max                       = zlim_max
      inp%amr_indicator_variables        = amr_indicator_variables
      inp%amr_smoothness_limits          = amr_smoothness_limits
      inp%amr_max_min_lim                = amr_max_min_lim
      inp%amr_threshold_lim              = amr_threshold_lim
      inp%amr_smoothness_qL2_limit       = amr_smoothness_qL2_limit
      inp%amr_mark_modes                 = amr_mark_modes
      inp%amr_mark_modes_use_baseline_decay = amr_mark_modes_use_baseline_decay
      inp%amr_mark_max_min               = amr_mark_max_min
      inp%amr_mark_random                = amr_mark_random
      inp%amr_mark_threshold             = amr_mark_threshold
      inp%amr_num_neigh_iter             = amr_num_neigh_iter
      inp%amr_mark_set2nc                = amr_mark_set2nc
      inp%nel_root_h                     = nel_root_h
      inp%filter_tracers_flg             = filter_tracers_flg
      inp%platformID                     = platformID
      inp%deviceID                       = deviceID
      inp%platformID2                    = platformID2
      inp%deviceID2                      = deviceID2
      inp%Nelems                         = Nelems
      inp%Nslices                        = Nslices
      inp%NslicesV                       = NslicesV
      inp%cpus_per_node                  = cpus_per_node
      inp%gpus_per_node                  = gpus_per_node
      inp%threads_per_process            = threads_per_process
      inp%platformWeight                 = platformWeight
      inp%platformWeight2                = platformWeight2
      inp%imass                          = imass
      inp%time_scale                     = time_scale
      inp%bcast_type                     = bcast_type
      inp%eqn_set                        = eqn_set
      inp%lrestart_file                  = lrestart_file
      inp%lsommerfeld                    = lsommerfeld
      inp%lxstretch                      = lxstretch
      inp%lystretch                      = lystretch
      inp%lzstretch                      = lzstretch
      inp%lgrid_only                     = lgrid_only
      inp%lread_external_grid            = lread_external_grid
      inp%lread_bc                       = lread_bc
      inp%lserial_grid_creation          = lserial_grid_creation
      inp%lparallel_grid_creation        = lparallel_grid_creation
      inp%lwrite_grid_ascii              = lwrite_grid_ascii
      inp%lout_ascii                     = lout_ascii
      inp%lout_asciimaya                 = lout_asciimaya
      inp%lcompute_barycenter            = lcompute_barycenter
      inp%luse_min_element_length        = luse_min_element_length
      inp%lrotating_flow                 = lrotating_flow
      inp%ladapt_timestep                = ladapt_timestep
      inp%lgpu                           = lgpu
      inp%luse_hybrid_cpu_gpu            = luse_hybrid_cpu_gpu
      inp%is_mlswe                       = is_mlswe
      inp%dg_integ_exact                 = dg_integ_exact
      inp%dump_data                      = dump_data
      inp%lcheck_conserved               = lcheck_conserved
      inp%time_scale_in                  = time_scale_in
      inp%dam_depth_in                   = dam_depth_in
      inp%dam_xlimit_in                  = dam_xlimit_in
      inp%dam_ylimit_in                  = dam_ylimit_in
      inp%bathymetry_shift               = bathymetry_shift
      inp%hump_config                    = hump_config
      inp%ibathymetry                    = ibathymetry
      inp%nonlinear                      = nonlinear
      inp%balance_flag                   = balance_flag
      inp%output_flag                    = output_flag
      inp%warp_grid                      = warp_grid
      inp%mesh_file                      = mesh_file
      inp%bathymetry_file                = bathymetry_file
      inp%lout_tree                      = lout_tree
      inp%lout_shoreline                 = lout_shoreline
      inp%lread_external_bathy           = lread_external_bathy
      inp%gravity_in                     = gravity_in
      inp%xdam_min                       = xdam_min
      inp%xdam_max                       = xdam_max
      inp%ydam_min                       = ydam_min
      inp%ydam_max                       = ydam_max
      inp%limit_threshold                = limit_threshold
      inp%llinear_pert                   = llinear_pert
      inp%bcl_time_method                = bcl_time_method
      inp%rk_bcl_FS                      = rk_bcl_FS

   end subroutine mod_input_create
 
 end module mod_input
 
