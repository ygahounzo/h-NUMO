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
subroutine print_header(gg, inp, flag, numproc)

    use mod_global_grid, only: grid_global
    use mod_input,       only: input

    implicit none

    type(grid_global), intent(in) :: gg
    type(input),       intent(in) :: inp
    integer, intent(in) :: flag, numproc

    real    :: factor, kiter
    integer :: nex, ney, nez

    if (flag == 0) then
        print *, '-------------------Begin Simulation----------------------------'
    elseif (flag == 1) then
        print *, '----------------------End Simulation---------------------------'
    end if

    print *, '---------------------------------------------------------------'
    write(*,'("eqn_set = ",a)') inp%eqn_set

    write(*,'("dt dt_btp time_initial time_final time_restart time_scale = ",6(e12.4,1x))') &
        inp%dt, inp%dt_btp, inp%time_initial/inp%time_scale, inp%time_final/inp%time_scale, &
        inp%time_restart/inp%time_scale, inp%time_scale
    write(*,'("nopx nopy nopz = ",3(i6,1x))') inp%nopx, inp%nopy, inp%nopz

    write(*,'("--------P4est Parameters--------")')
    write(*,'("p4est = ",l7)') .true.
    nex = inp%nelx * 2**inp%refinement_levels_h
    ney = inp%nely * 2**inp%refinement_levels_h
    nez = inp%nelz * 2**inp%refinement_levels_h

    write(*,'("nel_root_x refinement_levels_x  nelx = ",3(i6,1x))') inp%nelx, inp%refinement_levels_h, nex
    write(*,'("nel_root_y refinement_levels_y  nely = ",3(i6,1x))') inp%nely, inp%refinement_levels_h, ney
    write(*,'("--------P4est Parameters--------")')

    write(*,'("test_case  = ",a)') inp%test_case
    write(*,'("space_method = ",a)') inp%space_method
    write(*,'("ti_method_btp = ",a)') inp%ti_method_btp
    write(*,'("kstages = ",1(i6,1x))') inp%kstages
    write(*,'("fname_root = ",a)') inp%fname_root
    write(*,'("out_type = ",a)') inp%out_type

    if (inp%method_visc > 0) write(*,'("viscosity = ",f6.3)') inp%visc_mlswe

    write(*,'("nlayers npoin nelem nboun = ",5(i9,1x))') inp%nlayers, gg%npoin_g, &
        gg%nelem_g, gg%nboun_g
    write(*,'("lprint_diagnostics = ",l7)') inp%lprint_diagnostics
    write(*,'("numproc = ",(i6,1x))') numproc
    print *, '---------------------------------------------------------------'
    print *

end subroutine print_header
