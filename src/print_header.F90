!----------------------------------------------------------------------!
!>@brief This subroutine prints out the input/output Headers: e.g, initial set-up, CPU time, etc.
!>@author  Francis X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@update Added GPU output by F.X. Giraldo on 6/6/2018
!>
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: April 03, 2024
!----------------------------------------------------------------------!
subroutine print_header(flag,numproc)
  
    use mod_global_grid, only: npoin_g, npoin_g_cg, nelem_g, nboun_g

    use mod_initial, only: nvar

    use mod_input, only: eqn_set, dt, time_initial, time_final, time_restart, &
        nelx, nely, nelz, nopx, nopy, nopz, ztop, &
        icase, space_method, ti_method, si_dimension, kstages, filter_mux, filter_muy, &
        filter_muz, ifilter, filter_weight_type, filter_basis_type, &
        geometry_type, fname_root, out_type, visc_mlswe, &
        lprint_diagnostics, lp4est, lp6est, time_scale, &
        refinement_levels_h, nel_root_h, refinement_levels_v, nel_root_v, flux_type

    implicit none

    real factor
    real kiter
    integer flag, numproc
    integer nex,ney,nez

    if (flag == 0) then
        print*,'-------------------Begin Simulation----------------------------'
    elseif (flag == 1) then
        print*,'----------------------End Simulation---------------------------'
    end if

    print*,'---------------------------------------------------------------'
    write(*,'("eqn_set = ",a)')eqn_set

    write(*,'("dt time_initial time_final time_restart time_scale = ",5(e12.4,1x))')dt,time_initial/time_scale,time_final/time_scale,time_restart/time_scale,time_scale
    write(*,'("nopx nopy nopz = ",3(i6,1x))')nopx,nopy,nopz
    !P4EST and/or P6est Output
    if (lp4est .or. lp6est) then
        write(*,'("--------P4est Parameters--------")')
        write(*,'("lp4est = ",l7)')lp4est
        write(*,'("lp6est = ",l7)')lp6est
        nex=nelx*2**refinement_levels_h
        ney=nely*2**refinement_levels_h
        if(lp4est) then
            nez=nelz*2**refinement_levels_h
        else
            nez=nelz*2**refinement_levels_v
        endif
        write(*,'("nel_root_x refinement_levels_x  nelx = ",3(i6,1x))')nelx,refinement_levels_h,nex
        write(*,'("nel_root_y refinement_levels_y  nely = ",3(i6,1x))')nely,refinement_levels_h,ney
        write(*,'("nel_root_z refinement_levels_z  nelz = ",3(i6,1x))')nelz,refinement_levels_v,nez
        write(*,'("--------P4est Parameters--------")')
    else
        write(*,'("nelx nely nelz = ",3(i6,1x))')nelx,nely,nelz
    end if
    write(*,'("ztop = ",1(e12.4,1x))')ztop
    write(*,'("icase  = ",1(i6,1x))')icase
    write(*,'("space_method = ",a)')space_method
    write(*,'("bcl_flux_type = ",a)')flux_type
    write(*,'("ti_method = ",a)')ti_method
    write(*,'("si_dimension = ",a)')si_dimension
    write(*,'("kstages = ",1(i6,1x))')kstages

    write(*,'("mux muy muz ifilter = ",3(e12.4,1x),i7)')filter_mux,filter_muy,filter_muz,ifilter
    write(*,'("filter_weight_type = ",a," filter_basis_type = ",a)')filter_weight_type,filter_basis_type
    write(*,'("geometry_type = ",a)')geometry_type
    write(*,'("fname_root = ",a)')fname_root
    write(*,'("out_type = ",a)')out_type

    write(*,'("viscosity = ",f6.3)')visc_mlswe

    write(*,'("nvar npoin npoin_cg nelem nboun = ",5(i9,1x))')nvar,npoin_g,npoin_g_cg,nelem_g,nboun_g
    write(*,'("lprint_diagnostics = ",l7)')lprint_diagnostics
    write(*,'("numproc = ",(i6,1x))')numproc
    print*,'---------------------------------------------------------------'
    print*
  
end subroutine print_header