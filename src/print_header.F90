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
    use mod_input, only: eqn_set, dt, dt_btp, time_initial, time_final, time_restart, &
        nelx, nely, nelz, nopx, nopy, nopz, ztop, &
        test_case, space_method, ti_method_btp, kstages, &
        fname_root, out_type, visc_mlswe, &
        lprint_diagnostics, time_scale, &
        refinement_levels_h, nel_root_h, nlayers

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

    write(*,'("dt dt_btp time_initial time_final time_restart time_scale = ",6(e12.4,1x))')dt, &
            dt_btp, time_initial/time_scale,time_final/time_scale,time_restart/time_scale,time_scale
    write(*,'("nopx nopy nopz = ",3(i6,1x))')nopx,nopy,nopz
    !P4EST and/or P6est Output

    write(*,'("--------P4est Parameters--------")')
    write(*,'("p4est = ",l7)') .true.
    nex=nelx*2**refinement_levels_h
    ney=nely*2**refinement_levels_h
    nez=nelz*2**refinement_levels_h

    write(*,'("nel_root_x refinement_levels_x  nelx = ",3(i6,1x))')nelx,refinement_levels_h,nex
    write(*,'("nel_root_y refinement_levels_y  nely = ",3(i6,1x))')nely,refinement_levels_h,ney
    write(*,'("--------P4est Parameters--------")')

    write(*,'("test_case  = ",a)')test_case
    write(*,'("space_method = ",a)')space_method
    write(*,'("ti_method_btp = ",a)')ti_method_btp
    write(*,'("kstages = ",1(i6,1x))')kstages
    write(*,'("fname_root = ",a)')fname_root
    write(*,'("out_type = ",a)')out_type

    write(*,'("viscosity = ",f6.3)')visc_mlswe

    write(*,'("nlayers npoin nelem nboun = ",5(i9,1x))') nlayers,npoin_g, &
        nelem_g,nboun_g
    write(*,'("lprint_diagnostics = ",l7)')lprint_diagnostics
    write(*,'("numproc = ",(i6,1x))')numproc
    print*,'---------------------------------------------------------------'
    print*
  
end subroutine print_header
