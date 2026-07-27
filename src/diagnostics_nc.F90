subroutine diagnostics_nc(G, inp, b, init, gg, par, tsp, q, q_df, qb, itime, idone)

    use mod_grid,          only: grid
    use mod_input,         only: input
    use mod_basis,         only: basis
    use mod_global_grid,   only: grid_global
    use mod_parallel,      only: parallel_CS
    use mod_initial,       only: initial
    use mod_tensor,        only: tensor_CS
    use mod_gradient,      only: compute_vorticity_mlswe
    use mod_constants,     only: gravity
    use mod_mpi_utilities, only: irank, irank0
    use netcdf

    implicit none

    type(grid),        intent(in) :: G
    type(input),       intent(in) :: inp
    type(basis),       intent(in) :: b
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par
    type(initial),     intent(in) :: init
    type(tensor_CS),   intent(in) :: tsp

    real, intent(in)  :: q_df(inp%nvar_bcl, G%npoin, inp%nlayers)
    real, intent(in)  :: qb(inp%nvar_btp, G%npoin)
    integer, intent(in) :: itime, idone
    real, intent(out) :: q(5, G%npoin, inp%nlayers)

    character*5 :: tempchar
    character*4 :: num
    logical :: has_w
    integer :: i, j, k, iloop, ncid, dimids_2d(2)
    real :: q_gg(5, gg%npoin_g, inp%nlayers), ql(5, G%npoin), zbot_g(gg%npoin_g)
    real :: coord_dg_gathered(3, gg%npoin_g), qb_g(inp%nvar_btp, gg%npoin_g), q_g(5, gg%npoin_g)
    real, dimension(G%npoin,    inp%nlayers+1) :: mslwe_elevation
    real, dimension(gg%npoin_g, inp%nlayers+1) :: eta
    real, dimension(G%npoin)      :: rel_vort_l, abs_vort_l, pot_vort_l
    real, dimension(gg%npoin_g)   :: rel_vort_g1, abs_vort_g1, pot_vort_g1
    real :: rel_vort_gg(gg%npoin_g, inp%nlayers), abs_vort_gg(gg%npoin_g, inp%nlayers), pot_vort_gg(gg%npoin_g, inp%nlayers)
    character(len=100) :: fn

    character (len = *), parameter :: TIME_NAME      = "time"
    character (len = *), parameter :: NK_NAME        = "nlayers"
    character (len = *), parameter :: NPOIN_NAME     = "npoin"
    character (len = *), parameter :: DT_NAME        = "dt"
    character (len = *), parameter :: DT_BTP_NAME    = "dt_btp"
    character (len = *), parameter :: ZBOT_NAME      = "zbot"
    character (len = *), parameter :: X_NAME         = "x"
    character (len = *), parameter :: Y_NAME         = "y"
    character (len = *), parameter :: PB_NAME        = "pb"
    character (len = *), parameter :: PBUB_NAME      = "pbub"
    character (len = *), parameter :: PBVB_NAME      = "pbvb"
    character (len = *), parameter :: PBWB_NAME      = "pbwb"
    character (len = *), parameter :: H_NAME         = "h"
    character (len = *), parameter :: U_NAME         = "u"
    character (len = *), parameter :: V_NAME         = "v"
    character (len = *), parameter :: W_NAME         = "w"
    character (len = *), parameter :: ETA_NAME       = "eta"
    character (len = *), parameter :: INTERFACE_NAME = "zi"
    character (len = *), parameter :: RVORT_NAME     = "RelVorticity"
    character (len = *), parameter :: AVORT_NAME     = "AbsVorticity"
    character (len = *), parameter :: PVORT_NAME     = "PotVorticity"

    integer :: time_dimid, npoin_dimid, nlayers_dimid, parameter_dimid
    integer :: dt_varid, dt_btp_varid, zb_varid, pb_varid
    integer :: pbub_varid, pbvb_varid, pbwb_varid, h_varid, u_varid, v_varid, w_varid, e_varid
    integer :: zi_dimid, x_varid, y_varid
    integer :: rvort_varid, avort_varid, pvort_varid

    has_w = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

    do k = 1, inp%nlayers
        q(1,:,k) = (init%alpha_mlswe(k)/gravity)*q_df(1,:,k)
        q(2,:,k) = q_df(2,:,k) / q_df(1,:,k)
        q(3,:,k) = q_df(3,:,k) / q_df(1,:,k)
        if (has_w) then
            q(4,:,k) = q_df(4,:,k) / q_df(1,:,k)
        else
            q(4,:,k) = q_df(1,:,k)
        end if
    end do

    mslwe_elevation = 0.0
    mslwe_elevation(:, inp%nlayers+1) = init%zbot_df

    do k = inp%nlayers, 1, -1
        mslwe_elevation(:,k) = mslwe_elevation(:,k+1) + q(1,:,k)
    end do

    q(5,:,1) = mslwe_elevation(:,1)
    q(5,:,2:inp%nlayers) = mslwe_elevation(:,2:inp%nlayers)

    ! Absolute vorticity η=ζ+f and potential vorticity η/h per layer -- see
    ! mod_gradient.F90's compute_vorticity_mlswe for the formula (same routine
    ! used by the VTK path in write_output.F90). Gated by inp%lwrite_vorticity;
    ! the AbsVorticity/PotVorticity netCDF variables are only defined/written
    ! below when it is enabled.
    if (inp%lwrite_vorticity) then
        do k = 1, inp%nlayers
            call compute_vorticity_mlswe(G, b, tsp, init, has_w, q(2,:,k), q(3,:,k), q(4,:,k), q(1,:,k), &
                                          rel_vort_l, abs_vort_l, pot_vort_l)
            call gather_data(G, inp, gg, par, rel_vort_g1, rel_vort_l, 1)
            call gather_data(G, inp, gg, par, abs_vort_g1, abs_vort_l, 1)
            call gather_data(G, inp, gg, par, pot_vort_g1, pot_vort_l, 1)
            rel_vort_gg(:,k) = rel_vort_g1(:)
            abs_vort_gg(:,k) = abs_vort_g1(:)
            pot_vort_gg(:,k) = pot_vort_g1(:)
        end do
    end if

    ! Gather Data onto Head node
    do k = 1, inp%nlayers
        ql = q(:,:,k)
        call gather_data(G, inp, gg, par, q_g, ql, 5)
        q_gg(:,:,k) = q_g(:,:)
    enddo

    call gather_data(G, inp, gg, par, qb_g,               qb,             inp%nvar_btp)
    call gather_data(G, inp, gg, par, coord_dg_gathered,   G%coord,        3)
    call gather_data(G, inp, gg, par, zbot_g,              init%zbot_df,   1)

    eta(:, 1:inp%nlayers) = q_gg(5,:,1:inp%nlayers)
    eta(:, inp%nlayers+1) = zbot_g

    if (irank == irank0 .and. idone == 0) then

        ! Generate the name of the file to which the data will be written,
        ! and open the file.
        tempchar = 'mlswe'
        write(num,'(i4)') itime

        if(itime == 0) then
            num = '0000'
        else
            iloop = 3 - int(log10(real(itime)))
            do j = 1, iloop
                num(j:j) = '0'
            end do
        end if

        fn = tempchar//num//'.nc'
        ! Create NetCDF file
        call check(nf90_create(trim(fn), nf90_clobber+nf90_64bit_offset, ncid))

        ! Define dimensions
        call check(nf90_def_dim(ncid, TIME_NAME,      NF90_UNLIMITED,   time_dimid))
        call check(nf90_def_dim(ncid, NPOIN_NAME,     gg%npoin_g,       npoin_dimid))
        call check(nf90_def_dim(ncid, NK_NAME,        inp%nlayers,      nlayers_dimid))
        call check(nf90_def_dim(ncid, INTERFACE_NAME, inp%nlayers+1,    zi_dimid))

        call check(nf90_put_att(ncid, NF90_GLOBAL, "filename", trim(fn)))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "npoin", "Number of points in the mesh"))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "zi", "Number of interfaces"))

        ! Name each variable
        call check(nf90_def_var(ncid, DT_NAME,     NF90_DOUBLE, 1,          dt_varid))
        call check(nf90_put_att(ncid, dt_varid,     "name", "Baroclinic time step"))
        call check(nf90_put_att(ncid, dt_varid,     "units", "seconds"))
        call check(nf90_def_var(ncid, DT_BTP_NAME, NF90_DOUBLE, 1,          dt_btp_varid))
        call check(nf90_put_att(ncid, dt_btp_varid, "name", "Barotropic time step"))
        call check(nf90_put_att(ncid, dt_btp_varid, "units", "seconds"))
        call check(nf90_def_var(ncid, X_NAME, NF90_DOUBLE, npoin_dimid, x_varid))
        call check(nf90_put_att(ncid, x_varid, "name", "cartesian coordinates"))
        call check(nf90_put_att(ncid, x_varid, "axis", "X"))
        call check(nf90_def_var(ncid, Y_NAME, NF90_DOUBLE, npoin_dimid, y_varid))
        call check(nf90_put_att(ncid, y_varid, "name", "cartesian coordinates"))
        call check(nf90_put_att(ncid, y_varid, "axis", "Y"))
        call check(nf90_def_var(ncid, PB_NAME,   NF90_DOUBLE, npoin_dimid, pb_varid))
        call check(nf90_put_att(ncid, pb_varid,   "name", "Barotropic pressure pb"))
        call check(nf90_put_att(ncid, pb_varid,   "units", "N/m²"))
        call check(nf90_def_var(ncid, PBUB_NAME, NF90_DOUBLE, npoin_dimid, pbub_varid))
        call check(nf90_put_att(ncid, pbub_varid, "name", "Barotropic u-momentum"))
        call check(nf90_put_att(ncid, pbub_varid, "units", "kg·m/s"))
        call check(nf90_def_var(ncid, PBVB_NAME, NF90_DOUBLE, npoin_dimid, pbvb_varid))
        call check(nf90_put_att(ncid, pbvb_varid, "name", "Barotropic v-momentum"))
        call check(nf90_put_att(ncid, pbvb_varid, "units", "kg·m/s"))
        if (has_w) then
            call check(nf90_def_var(ncid, PBWB_NAME, NF90_DOUBLE, npoin_dimid, pbwb_varid))
            call check(nf90_put_att(ncid, pbwb_varid, "name", "Barotropic w-momentum"))
            call check(nf90_put_att(ncid, pbwb_varid, "units", "kg·m/s"))
        end if
        dimids_2d = (/npoin_dimid, nlayers_dimid/)
        call check(nf90_def_var(ncid, H_NAME, NF90_DOUBLE, dimids_2d, h_varid))
        call check(nf90_put_att(ncid, h_varid, "name", "Layer thickness"))
        call check(nf90_put_att(ncid, h_varid, "units", "m"))
        call check(nf90_def_var(ncid, U_NAME, NF90_DOUBLE, dimids_2d, u_varid))
        call check(nf90_put_att(ncid, u_varid, "name", "Baroclinic u-velocity"))
        call check(nf90_put_att(ncid, u_varid, "units", "m/s"))
        call check(nf90_def_var(ncid, V_NAME, NF90_DOUBLE, dimids_2d, v_varid))
        call check(nf90_put_att(ncid, v_varid, "name", "Baroclinic v-velocity"))
        call check(nf90_put_att(ncid, v_varid, "units", "m/s"))
        if (has_w) then
            call check(nf90_def_var(ncid, W_NAME, NF90_DOUBLE, dimids_2d, w_varid))
            call check(nf90_put_att(ncid, w_varid, "name", "Baroclinic w-velocity"))
            call check(nf90_put_att(ncid, w_varid, "units", "m/s"))
        end if
        dimids_2d = (/npoin_dimid, zi_dimid/)
        call check(nf90_def_var(ncid, ETA_NAME, NF90_DOUBLE, dimids_2d, e_varid))
        call check(nf90_put_att(ncid, e_varid, "name", "Interface Height Relative to Mean Sea Level"))
        call check(nf90_put_att(ncid, e_varid, "units", "m"))

        if (inp%lwrite_vorticity) then
            dimids_2d = (/npoin_dimid, nlayers_dimid/)
            call check(nf90_def_var(ncid, RVORT_NAME, NF90_DOUBLE, dimids_2d, rvort_varid))
            call check(nf90_put_att(ncid, rvort_varid, "name", "Relative vorticity"))
            call check(nf90_put_att(ncid, rvort_varid, "units", "s-1"))
            call check(nf90_def_var(ncid, AVORT_NAME, NF90_DOUBLE, dimids_2d, avort_varid))
            call check(nf90_put_att(ncid, avort_varid, "name", "Absolute vorticity"))
            call check(nf90_put_att(ncid, avort_varid, "units", "s-1"))
            call check(nf90_def_var(ncid, PVORT_NAME, NF90_DOUBLE, dimids_2d, pvort_varid))
            call check(nf90_put_att(ncid, pvort_varid, "name", "Potential vorticity"))
            call check(nf90_put_att(ncid, pvort_varid, "units", "m-1 s-1"))
        end if

        ! End define mode
        call check(nf90_enddef(ncid))

        ! Write data
        call check(nf90_put_var(ncid, dt_varid,     inp%dt))
        call check(nf90_put_var(ncid, dt_btp_varid, inp%dt_btp))
        call check(nf90_put_var(ncid, x_varid,    coord_dg_gathered(1,:)))
        call check(nf90_put_var(ncid, y_varid,    coord_dg_gathered(2,:)))
        call check(nf90_put_var(ncid, pb_varid,   qb_g(1,:)))
        call check(nf90_put_var(ncid, pbub_varid, qb_g(3,:)))
        call check(nf90_put_var(ncid, pbvb_varid, qb_g(4,:)))
        if (has_w) call check(nf90_put_var(ncid, pbwb_varid, qb_g(5,:)))
        call check(nf90_put_var(ncid, h_varid, q_gg(1,:,:)))
        call check(nf90_put_var(ncid, u_varid, q_gg(2,:,:)))
        call check(nf90_put_var(ncid, v_varid, q_gg(3,:,:)))
        if (has_w) call check(nf90_put_var(ncid, w_varid, q_gg(4,:,:)))
        call check(nf90_put_var(ncid, e_varid, eta(:,:)))
        if (inp%lwrite_vorticity) then
            call check(nf90_put_var(ncid, rvort_varid, rel_vort_gg(:,:)))
            call check(nf90_put_var(ncid, avort_varid, abs_vort_gg(:,:)))
            call check(nf90_put_var(ncid, pvort_varid, pot_vort_gg(:,:)))
        end if

        ! Close NetCDF file
        call check(nf90_close(ncid))

    end if

    contains

    !-----------------------!
    !>@brief Check status
    !-----------------------!
    subroutine check(status)
        integer, intent(in) :: status

        if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check

end subroutine diagnostics_nc
