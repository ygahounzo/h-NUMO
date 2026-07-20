!----------------------------------------------------------------------!
!> Write a parameter documentation file HNUMO_parameter_doc.all.
!> Called once at startup (rank 0 only) after all parameters are set.
!>
!> Format: NAME = value  <column 55>  ! [units] Description
!----------------------------------------------------------------------!
subroutine write_parameter_doc(numproc)

   use mod_global_grid, only: npoin_g, nelem_g, nboun_g

   use mod_basis, only: nglx

   use mod_input, only: test_case, eqn_set, fname_root, space_method, &
        nlayers, nelx, nely, refinement_levels_h, nopx, nopy, &
        time_initial, time_final, time_scale, dt, dt_btp, ti_method_btp, kstages, &
        f0, beta, ad_mlswe, cd_mlswe, dp_tau_bot, dp_tau_wind, botfr, &
        method_visc, visc_mlswe, max_shear_dz, &
        adjust_H_vertical_sum, adjust_bcl_mom_flux, dg_integ_exact, &
        out_type, nvtk_files, format_vtk, lprint_diagnostics, iprint_diagnostics, &
        lcheck_conserved, dump_data, &
        time_restart, irestart_file_number

   implicit none

   integer, intent(in) :: numproc

   integer :: fu, nex, ney

   ! Column at which all '!' comments begin — change once here to re-align everything.
   integer, parameter :: CCOL = 55

   nex = nelx * 2**refinement_levels_h
   ney = nely * 2**refinement_levels_h

   open(newunit=fu, file='HNUMO_parameter_doc.all', status='replace', action='write')

   write(fu,'(a)') '! HNUMO parameter documentation'
   write(fu,'(a)') '! Auto-generated at run start.'
   write(fu,'(a)') '! Format:  NAME = value                              ! [units] Description'
   write(fu,'(a)') '!'

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Run Identity ==='
   write(fu,'(a)')
   write(fu,'(a,a,t55,a)')   'TEST_CASE = "',   trim(test_case)//'"',      '! test case identifier'
   write(fu,'(a,a,t55,a)')   'EQN_SET = "',     trim(eqn_set)//'"',        '! equation set (mlswe, ...)'
   write(fu,'(a,a,t55,a)')   'SPACE_METHOD = "',trim(space_method)//'"',   '! spatial discretization (cgc, cgd, dg)'
   write(fu,'(a,a,t55,a)')   'FNAME_ROOT = "',  trim(fname_root)//'"',     '! output file root name'
   write(fu,'(a,i0,t55,a)')  'NPROC = ',        numproc,                   '! number of MPI ranks'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Grid ==='
   write(fu,'(a)')
   write(fu,'(a,i0,t55,a)')  'NLAYERS = ',         nlayers,                   '! number of baroclinic layers'
   write(fu,'(a,i0,t55,a)')  'NELX_ROOT = ',       nelx,                      '! root elements in x (before refinement)'
   write(fu,'(a,i0,t55,a)')  'NELY_ROOT = ',       nely,                      '! root elements in y (before refinement)'
   write(fu,'(a,i0,t55,a)')  'REFINEMENT_LEVELS = ', refinement_levels_h,     '! uniform h-refinement levels'
   write(fu,'(a,i0,t55,a)')  'NELX = ',            nex,                       '! effective elements in x after refinement'
   write(fu,'(a,i0,t55,a)')  'NELY = ',            ney,                       '! effective elements in y after refinement'
   write(fu,'(a,i0,t55,a)')  'NOPX = ',            nopx,                      '! polynomial order in x'
   write(fu,'(a,i0,t55,a)')  'NOPY = ',            nopy,                      '! polynomial order in y'
   write(fu,'(a,i0,t55,a)')  'NGLX = ',            nglx,                      '! quadrature points per edge (nopx+1)'
   write(fu,'(a,i0,t55,a)')  'NELEM_GLOBAL = ',    nelem_g,                   '! total global elements'
   write(fu,'(a,i0,t55,a)')  'NPOIN_GLOBAL = ',    npoin_g,                   '! total global DG DOFs'
   write(fu,'(a,i0,t55,a)')  'NBOUN_GLOBAL = ',    nboun_g,                   '! total global MPI boundary faces'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Time Integration ==='
   write(fu,'(a)')
   write(fu,'(a,es14.6,t55,a)') 'TIME_INITIAL = ',    time_initial,   '! [s] simulation start time'
   write(fu,'(a,es14.6,t55,a)') 'TIME_FINAL = ',      time_final,     '! [s] simulation end time'
   write(fu,'(a,es14.6,t55,a)') 'TIME_SCALE = ',      time_scale,     '! [s] time scale for diagnostic output'
   write(fu,'(a,es14.6,t55,a)') 'DT = ',              dt,             '! [s] baroclinic time step'
   write(fu,'(a,es14.6,t55,a)') 'DT_BTP = ',          dt_btp,         '! [s] barotropic time step'
   write(fu,'(a,a,t55,a)')      'TI_METHOD_BTP = "',  trim(ti_method_btp)//'"',  '! barotropic time integration method'
   write(fu,'(a,i0,t55,a)')     'KSTAGES = ',         kstages,        '! barotropic RK stages'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Physics ==='
   write(fu,'(a)')
   write(fu,'(a,es14.6,t55,a)') 'F0 = ',              f0,             '! [s-1] reference Coriolis parameter'
   write(fu,'(a,es14.6,t55,a)') 'BETA = ',            beta,           '! [m-1 s-1] beta-plane parameter'
   write(fu,'(a,es14.6,t55,a)') 'AD_MLSWE = ',        ad_mlswe,       '! [m2 s-1] momentum diffusivity'
   write(fu,'(a,es14.6,t55,a)') 'CD_MLSWE = ',        cd_mlswe,       '! [m2 s-1] thickness diffusivity'
   write(fu,'(a,es14.6,t55,a)') 'DP_TAU_BOT = ',      dp_tau_bot,     '! bottom drag coefficient'
   write(fu,'(a,es14.6,t55,a)') 'DP_TAU_WIND = ',     dp_tau_wind,    '! wind stress coefficient'
   write(fu,'(a,i0,t55,a)')     'BOTFR = ',           botfr,          '! bottom friction type (0=none)'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Viscosity / Stabilization ==='
   write(fu,'(a)')
   write(fu,'(a,i0,t55,a)')     'METHOD_VISC = ',     method_visc,    '! LDG viscosity (0=off, 1=LDG)'
   write(fu,'(a,es14.6,t55,a)') 'VISC_MLSWE = ',      visc_mlswe,     '! [m2 s-1] explicit LDG viscosity coefficient'
   write(fu,'(a,es14.6,t55,a)') 'MAX_SHEAR_DZ = ',    max_shear_dz,   '! [s-1] maximum allowed vertical shear'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Layer Consistency ==='
   write(fu,'(a)')
   write(fu,'(a,i0,t55,a)')     'ADJUST_H_VERTICAL_SUM = ', adjust_H_vertical_sum, &
                                                              '! enforce vertical sum conservation (0=off)'
   write(fu,'(a,i0,t55,a)')     'ADJUST_BCL_MOM_FLUX = ',   adjust_bcl_mom_flux, &
                                                              '! adjust BCL momentum flux for consistency (0=off)'
   write(fu,'(a,l1,t55,a)')     'DG_INTEG_EXACT = ',        dg_integ_exact, &
                                                              '! use exact (over-integrated) DG quadrature'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Diagnostics and Output ==='
   write(fu,'(a)')
   write(fu,'(a,a,t55,a)')   'OUT_TYPE = "',          trim(out_type)//'"',        '! output format (netcdf, vtk, ...)'
   write(fu,'(a,i0,t55,a)') 'NVTK_FILES = ',         nvtk_files,                 '! number of parallel VTK output files'
   write(fu,'(a,a,t55,a)')   'FORMAT_VTK = "',        trim(format_vtk)//'"',      '! VTK format (ascii, binary)'
   write(fu,'(a,l1,t55,a)') 'LPRINT_DIAGNOSTICS = ', lprint_diagnostics,         '! enable run-time diagnostics'
   write(fu,'(a,i0,t55,a)') 'IPRINT_DIAGNOSTICS = ', iprint_diagnostics,         '! print diagnostics every N time steps'
   write(fu,'(a,l1,t55,a)') 'LCHECK_CONSERVED = ',   lcheck_conserved,           '! check mass/energy conservation'
   write(fu,'(a,l1,t55,a)') 'DUMP_DATA = ',          dump_data,                  '! write output data to disk'
   write(fu,'(a)')

   ! -----------------------------------------------------------------------
   write(fu,'(a)') '! === Restart ==='
   write(fu,'(a)')
   write(fu,'(a,es14.6,t55,a)') 'TIME_RESTART = ',         time_restart,          '! [s] restart dump interval'
   write(fu,'(a,i0,t55,a)')     'IRESTART_FILE_NUMBER = ', irestart_file_number,  '! restart file index to read (0=cold start)'
   write(fu,'(a)')

   close(fu)

   write(*,'(a)') 'HNUMO_parameter_doc.all written.'

end subroutine write_parameter_doc
