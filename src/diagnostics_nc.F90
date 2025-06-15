subroutine diagnostics_nc(q,q_df,qb,itime,idone)

    use mod_input, only: nlayers, dt, dt_btp
    use mod_global_grid, only: coord_g, npoin_g
    use mod_grid, only: npoin, intma, coord, nelem
    use mod_mpi_utilities, only: irank, irank0
    use mod_initial, only: alpha_mlswe, zbot_df, one_over_pbprime_df
    use mod_constants, only: gravity
    use netcdf

    implicit none

    real, intent(in)    :: q_df(3,npoin,nlayers)
    real, intent(in)    :: qb(4,npoin)
    integer, intent(in) :: itime,idone
    real, intent(out)   :: q(5,npoin,nlayers)

    character*5 :: tempchar
    character*4 :: num
    integer     :: i,j,k,iloop, ncid, dimids_2d(2)
    real :: q_gg(5,npoin_g,nlayers), ql(5,npoin), zbot_g(npoin_g), coord_dg_gathered(3,npoin_g), qb_g(4,npoin_g), q_g(5,npoin_g)
    real, dimension(npoin,nlayers+1) :: mslwe_elevation
    real, dimension(npoin_g,nlayers+1) :: eta
    character(len=100) ::  fn

    character (len = *), parameter :: TIME_NAME="time"
    character (len = *), parameter :: NK_NAME="nlayers"
    character (len = *), parameter :: NPOIN_NAME="npoin"
    character (len = *), parameter :: DT_NAME="dt"
    character (len = *), parameter :: DT_BTP_NAME="dt_btp"
    character (len = *), parameter :: ZBOT_NAME="zbot"
    character (len = *), parameter :: X_NAME="x"
    character (len = *), parameter :: Y_NAME="y"
    character (len = *), parameter :: PB_NAME="pb"
    character (len = *), parameter :: PBUB_NAME="pbub"
    character (len = *), parameter :: PBVB_NAME="pbvb"
    character (len = *), parameter :: H_NAME="h"
    character (len = *), parameter :: U_NAME="u"
    character (len = *), parameter :: V_NAME="v"
    character (len = *), parameter :: ETA_NAME="eta"
    character (len = *), parameter :: INTERFACE_NAME="zi"

    integer :: time_dimid, npoin_dimid, nlayers_dimid, parameter_dimid
    integer :: dt_varid, dt_btp_varid, zb_varid, pb_varid
    integer :: pbub_varid, pbvb_varid, h_varid, u_varid, v_varid, e_varid
    integer :: zi_dimid, x_varid, y_varid

    do k = 1,nlayers
        q(1,:,k) = (alpha_mlswe(k)/gravity)*q_df(1,:,k)
        q(2,:,k) = q_df(2,:,k) / q_df(1,:,k)
        q(3,:,k) = q_df(3,:,k) / q_df(1,:,k)
        q(4,:,k) = q_df(1,:,k)
    end do

    mslwe_elevation = 0.0
    mslwe_elevation(:,nlayers+1) = zbot_df

    do k = nlayers,1,-1
        mslwe_elevation(:,k) = mslwe_elevation(:,k+1) + q(1,:,k)
    end do

    q(5,:,1) = qb(1,:)*one_over_pbprime_df(:) - 1.0
    q(5,:,2:nlayers) = mslwe_elevation(:,2:nlayers)

    ! Gather Data onto Head node
    do k = 1, nlayers
        ql = q(:,:,k)
        call gather_data(q_g,ql,5)
        q_gg(:,:,k) = q_g(:,:)
    enddo

    call gather_data(qb_g,qb,4)
    call gather_data(coord_dg_gathered,coord,3)
    call gather_data(zbot_g,zbot_df,1)

    eta(:,1:nlayers) = q_gg(5,:,1:nlayers)
    eta(:,nlayers+1) = zbot_g

    if (irank == irank0 .and. idone == 0) then

        ! Generate the name of the file to which the data will be written,
        ! and open the file.  
        tempchar = 'mlswe'
        write(num,'(i4)')itime

        if(itime == 0) then
            num = '0000'
        else
            iloop=3 - int(log10(real(itime)))
            do j=1,iloop
                num(j:j) = '0'
            end do
        end if

        fn = tempchar//num//'.nc'
        ! Create NetCDF file
        call check(nf90_create(trim(fn), nf90_clobber+nf90_64bit_offset, ncid))

        ! Define dimensions
        call check(nf90_def_dim(ncid, TIME_NAME, NF90_UNLIMITED, time_dimid))
        call check(nf90_def_dim(ncid, NPOIN_NAME, npoin_g, npoin_dimid))
        call check(nf90_def_dim(ncid, NK_NAME, nlayers, nlayers_dimid))
        call check(nf90_def_dim(ncid, INTERFACE_NAME, nlayers+1, zi_dimid))

        call check(nf90_put_att(ncid, NF90_GLOBAL, "filename", trim(fn)))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "npoin", "Number of points in the mesh"))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "zi", "Number of interfaces"))

        ! Name each variable
        call check(nf90_def_var(ncid, DT_NAME, NF90_DOUBLE, 1,dt_varid))
        call check(nf90_put_att(ncid, dt_varid, "name", "Baroclinic time step"))
        call check(nf90_put_att(ncid, dt_varid, "units", "seconds"))
        call check(nf90_def_var(ncid, DT_BTP_NAME, NF90_DOUBLE,1,dt_btp_varid))
        call check(nf90_put_att(ncid, dt_btp_varid, "name", "Barotropic time step"))
        call check(nf90_put_att(ncid, dt_btp_varid, "units", "seconds"))
        call check(nf90_def_var(ncid, X_NAME, NF90_DOUBLE,npoin_dimid,x_varid))
        call check(nf90_put_att(ncid, x_varid, "name", "cartesian coordinates"))
        call check(nf90_put_att(ncid, x_varid, "axis", "X"))
        call check(nf90_def_var(ncid, Y_NAME, NF90_DOUBLE,npoin_dimid,y_varid))
        call check(nf90_put_att(ncid, y_varid, "name", "cartesian coordinates"))
        call check(nf90_put_att(ncid, y_varid, "axis", "Y"))
        call check(nf90_def_var(ncid, PB_NAME, NF90_DOUBLE,npoin_dimid,pb_varid))
        call check(nf90_put_att(ncid, pb_varid, "name", "Barotropic pressure pb"))
        call check(nf90_put_att(ncid, pb_varid, "units", "N/m²"))
        call check(nf90_def_var(ncid, PBUB_NAME, NF90_DOUBLE,npoin_dimid,pbub_varid))
        call check(nf90_put_att(ncid, pbub_varid, "name", "Barotropic u-momentum"))
        call check(nf90_put_att(ncid, pbub_varid, "units", "kg·m/s"))
        call check(nf90_def_var(ncid, PBVB_NAME, NF90_DOUBLE,npoin_dimid,pbvb_varid))
        call check(nf90_put_att(ncid, pbvb_varid, "name", "Barotropic v-momentum"))
        call check(nf90_put_att(ncid, pbvb_varid, "units", "kg·m/s"))
        dimids_2d = (/npoin_dimid, nlayers_dimid/)
        call check(nf90_def_var(ncid, H_NAME, NF90_DOUBLE,dimids_2d,h_varid))
        call check(nf90_put_att(ncid, h_varid, "name", "Layer thickness"))
        call check(nf90_put_att(ncid, h_varid, "units", "m"))
        call check(nf90_def_var(ncid, U_NAME, NF90_DOUBLE,dimids_2d,u_varid))
        call check(nf90_put_att(ncid, u_varid, "name", "Baroclinic u-velocity"))
        call check(nf90_put_att(ncid, u_varid, "units", "m/s"))
        call check(nf90_def_var(ncid, V_NAME, NF90_DOUBLE,dimids_2d,v_varid))
        call check(nf90_put_att(ncid, v_varid, "name", "Baroclinic v-velocity"))
        call check(nf90_put_att(ncid, v_varid, "units", "m/s"))
        dimids_2d = (/npoin_dimid, zi_dimid/)
        call check(nf90_def_var(ncid, ETA_NAME, NF90_DOUBLE,dimids_2d,e_varid))
        call check(nf90_put_att(ncid, e_varid, "name", "Interface Height Relative to Mean Sea Level"))
        call check(nf90_put_att(ncid, e_varid, "units", "m"))

        ! End define mode
        call check(nf90_enddef(ncid))

        ! Write data
        call check(nf90_put_var(ncid, dt_varid, dt))
        call check(nf90_put_var(ncid, dt_btp_varid, dt_btp))
        call check(nf90_put_var(ncid, x_varid, coord_dg_gathered(1,:)))
        call check(nf90_put_var(ncid, y_varid, coord_dg_gathered(2,:)))
        call check(nf90_put_var(ncid, pb_varid, qb_g(1,:)))
        call check(nf90_put_var(ncid, pbub_varid, qb_g(3,:)))
        call check(nf90_put_var(ncid, pbvb_varid, qb_g(4,:)))
        call check(nf90_put_var(ncid, h_varid, q_gg(1,:,:)))
        call check(nf90_put_var(ncid, u_varid, q_gg(2,:,:)))
        call check(nf90_put_var(ncid, v_varid, q_gg(3,:,:)))
        call check(nf90_put_var(ncid, e_varid, eta(:,:)))

        ! Close NetCDF file
        call check(nf90_close(ncid))

    end if

    contains

    !-----------------------!
    !>@brief Check status
    !-----------------------!
    subroutine check(status)
        integer, intent ( in) :: status

        if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check
end subroutine diagnostics_nc



