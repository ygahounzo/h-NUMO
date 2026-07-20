!----------------------------------------------------------------------!
!> Write a parameter documentation file HNUMO_parameter_doc.all.
!> Called once at startup (rank 0 only) after all parameters are set.
!>
!> Format: NAME = value  <column 55>  ! [units] Description
!----------------------------------------------------------------------!
subroutine write_parameter_doc(gg, inp, b, numproc)

   use mod_global_grid, only: grid_global
   use mod_input,       only: input
   use mod_basis,       only: basis

   implicit none

   type(grid_global), intent(in) :: gg
   type(input),       intent(in) :: inp
   type(basis),       intent(in) :: b
   integer,           intent(in) :: numproc

   integer :: fu, nex, ney

   ! Column at which all '!' comments begin — change once here to re-align everything.
   integer, parameter :: CCOL = 55

   nex = inp%nelx * 2**inp%refinement_levels_h
   ney = inp%nely * 2**inp%refinement_levels_h

   open(newunit=fu, file='HNUMO_parameter_doc.all', status='replace', action='write')

   write(fu,'(a)') '! HNUMO parameter documentation'
   write(fu,'(a)') '! Auto-generated at run start.'
   write(fu,'(a)') '! Format:  NAME = value                              ! [units] Description'
   write(fu,'(a)') '!'

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Run Identity ==='
   write(fu,'(a)')
   write(fu,'(a,a,t55,a)')   'TEST_CASE = "',   trim(inp%test_case)//'"',      '! test case identifier'
   write(fu,'(a,a,t55,a)')   'EQN_SET = "',     trim(inp%eqn_set)//'",',       '! equation set (mlswe, ...)'
   write(fu,'(a,a,t55,a)')   'FNAME_ROOT = "',  trim(inp%fname_root)//'"',     '! output file root name'
   write(fu,'(a,i0,t55,a)')  'NPROC = ',        numproc,                       '! number of MPI ranks'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Grid ==='
   write(fu,'(a)')
   write(fu,'(a,a,t55,a)')   'GEOMETRY_TYPE = "', trim(inp%geometry_type)//'"', '! mesh geometry (sphere, cube, ...)'
   write(fu,'(a,i0,t55,a)')  'NLAYERS = ',         inp%nlayers,                 '! number of baroclinic layers'
   write(fu,'(a,i0,t55,a)')  'NELX_ROOT = ',       inp%nelx,                    '! root elements in x (before refinement)'
   write(fu,'(a,i0,t55,a)')  'NELY_ROOT = ',       inp%nely,                    '! root elements in y (before refinement)'
   write(fu,'(a,i0,t55,a)')  'REFINEMENT_LEVELS = ', inp%refinement_levels_h,   '! uniform h-refinement levels'
   write(fu,'(a,i0,t55,a)')  'NELX = ',            nex,                         '! effective elements in x after refinement'
   write(fu,'(a,i0,t55,a)')  'NELY = ',            ney,                         '! effective elements in y after refinement'
   write(fu,'(a,i0,t55,a)')  'NOPX = ',            inp%nopx,                    '! polynomial order in x'
   write(fu,'(a,i0,t55,a)')  'NOPY = ',            inp%nopy,                    '! polynomial order in y'
   write(fu,'(a,i0,t55,a)')  'NGLX = ',            b%ngl,                       '! quadrature points per edge (nopx+1)'
   write(fu,'(a,i0,t55,a)')  'NELEM_GLOBAL = ',    gg%nelem_g,                  '! total global elements'
   write(fu,'(a,i0,t55,a)')  'NPOIN_GLOBAL = ',    gg%npoin_g,                  '! total global DG DOFs'
   write(fu,'(a,i0,t55,a)')  'NBOUN_GLOBAL = ',    gg%nboun_g,                  '! total global MPI boundary faces'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Time Integration ==='
   write(fu,'(a)')
   write(fu,'(a,es14.6,t55,a)') 'TIME_INITIAL = ',    inp%time_initial,   '! [s] simulation start time'
   write(fu,'(a,es14.6,t55,a)') 'TIME_FINAL = ',      inp%time_final,     '! [s] simulation end time'
   write(fu,'(a,es14.6,t55,a)') 'TIME_SCALE = ',      inp%time_scale,     '! [s] time scale for diagnostic output'
   write(fu,'(a,es14.6,t55,a)') 'DT = ',              inp%dt,             '! [s] baroclinic time step'
   write(fu,'(a,es14.6,t55,a)') 'DT_BTP = ',          inp%dt_btp,         '! [s] barotropic time step'
   write(fu,'(a,a,t55,a)')      'TI_METHOD_BTP = "',  trim(inp%ti_method_btp)//'"',  '! barotropic time integration method'
   write(fu,'(a,a,t55,a)')      'BCL_TIME_METHOD = "',trim(inp%bcl_time_method)//'"', '! baroclinic time integration method'
   write(fu,'(a,i0,t55,a)')     'KSTAGES_BTP = ',     inp%kstages,        '! barotropic RK stages'
   write(fu,'(a,i0,t55,a)')     'KSTAGES_BCL = ',     inp%kstages_bcl,    '! baroclinic RK stages'
   write(fu,'(a,l1,t55,a)')     'RK_BCL_FS = ',            inp%rk_bcl_FS,            '! baroclinic free-surface coupling flag'
   write(fu,'(a,l1,t55,a)')     'IMPLICIT_CORIOLIS_SPH = ', inp%implicit_coriolis_sph, '! 3D implicit Coriolis for sphere (LSRK3 only)'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Physics ==='
   write(fu,'(a)')
   write(fu,'(a,es14.6,t55,a)') 'F0 = ',              inp%f0,             '! [s-1] reference Coriolis parameter'
   write(fu,'(a,es14.6,t55,a)') 'BETA = ',            inp%beta,           '! [m-1 s-1] beta-plane parameter'
   write(fu,'(a,es14.6,t55,a)') 'AD_MLSWE = ',        inp%ad_mlswe,       '! [m2 s-1] momentum diffusivity'
   write(fu,'(a,es14.6,t55,a)') 'CD_MLSWE = ',        inp%cd_mlswe,       '! [m2 s-1] thickness diffusivity'
   write(fu,'(a,es14.6,t55,a)') 'DP_TAU_BOT = ',      inp%dp_tau_bot,     '! bottom drag coefficient'
   write(fu,'(a,es14.6,t55,a)') 'DP_TAU_WIND = ',     inp%dp_tau_wind,    '! wind stress coefficient'
   write(fu,'(a,i0,t55,a)')     'BOTFR = ',           inp%botfr,          '! bottom friction type (0=none)'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Viscosity / Stabilization ==='
   write(fu,'(a)')
   write(fu,'(a,i0,t55,a)')     'METHOD_VISC = ',     inp%method_visc,    '! LDG viscosity (0=off, 1=LDG)'
   write(fu,'(a,es14.6,t55,a)') 'VISC_MLSWE = ',      inp%visc_mlswe,     '! [m2 s-1] explicit LDG viscosity coefficient'
   write(fu,'(a,es14.6,t55,a)') 'SIPG_CONSTANT = ',   inp%SIPG_constant,  '! SIP penalty multiplier (0=LDG, 1=Shahbazi)'
   write(fu,'(a,es14.6,t55,a)') 'C_APE = ',           inp%c_APE,          '! [m2] Chen (2025) APE stabilization coefficient'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Layer Thickness Limits ==='
   write(fu,'(a)')
   write(fu,'(a,es14.6,t55,a)') 'DRY_CUTOFF = ',      inp%dry_cutoff,     '! [m] minimum layer thickness (wetting/drying)'
   write(fu,'(a,es14.6,t55,a)') 'H_CUTOFF1 = ',       inp%h_cutoff1,      '! [m] lower health-check thickness threshold'
   write(fu,'(a,es14.6,t55,a)') 'H_CUTOFF2 = ',       inp%h_cutoff2,      '! [m] upper health-check thickness threshold'
   write(fu,'(a,i0,t55,a)')     'ADJUST_H_VERTICAL_SUM = ', inp%adjust_H_vertical_sum, &
                                                              '! enforce vertical sum conservation (0=off)'
   write(fu,'(a,i0,t55,a)')     'ADJUST_BCL_MOM_FLUX = ',   inp%adjust_bcl_mom_flux, &
                                                              '! adjust BCL momentum flux for consistency (0=off)'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Diagnostics and Output ==='
   write(fu,'(a)')
   write(fu,'(a,a,t55,a)')   'OUT_TYPE = "',          trim(inp%out_type)//'"',    '! output format (netcdf, vtk, ...)'
   write(fu,'(a,i0,t55,a)') 'NVTK_FILES = ',         inp%nvtk_files,             '! number of parallel VTK output files'
   write(fu,'(a,l1,t55,a)') 'LPRINT_DIAGNOSTICS = ', inp%lprint_diagnostics,     '! enable run-time diagnostics'
   write(fu,'(a,i0,t55,a)') 'IPRINT_DIAGNOSTICS = ', inp%iprint_diagnostics,     '! print diagnostics every N time steps'
   write(fu,'(a,l1,t55,a)') 'LCHECK_CONSERVED = ',   inp%lcheck_conserved,       '! check mass/energy conservation'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Restart ==='
   write(fu,'(a)')
   write(fu,'(a,es14.6,t55,a)') 'TIME_RESTART = ',         inp%time_restart,          '! [s] restart dump interval'
   write(fu,'(a,i0,t55,a)')     'IRESTART_FILE_NUMBER = ', inp%irestart_file_number,  '! restart file index to read (0=cold start)'
   write(fu,'(a)')

   close(fu)

   write(*,'(a)') 'HNUMO_parameter_doc.all written.'

end subroutine write_parameter_doc
