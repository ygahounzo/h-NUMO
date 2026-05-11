subroutine diagnostics_btp_nc(qb,itime,idone)

   use mod_input, only: nlayers, dt, dt_btp, fname_root
   use mod_global_grid, only: coord_g, npoin_g
   use mod_grid, only: npoin, intma, coord, nelem
   use mod_mpi_utilities, only: irank, irank0
   use mod_initial, only: alpha_mlswe, zbot_df
   use mod_constants, only: gravity
   use netcdf

   implicit none

   real, intent(in)    :: qb(4,npoin)
   integer, intent(in) :: itime,idone

   character*3 :: tempchar
   character*4 :: num
   integer     :: i,j,k,iloop, ncid
   real :: zbot_g(npoin_g), coord_dg_gathered(3,npoin_g), qb_g(4,npoin_g)
   character(len=100) ::  fn

   character (len = *), parameter :: TIME_NAME="time"
   character (len = *), parameter :: NK_NAME="nlayers"
   character (len = *), parameter :: NPOIN_NAME="npoin"
   character (len = *), parameter :: DT_BTP_NAME="dt_btp"
   character (len = *), parameter :: ZBOT_NAME="zbot"
   character (len = *), parameter :: X_NAME="x"
   character (len = *), parameter :: Y_NAME="y"
   character (len = *), parameter :: H_NAME="h"
   character (len = *), parameter :: UB_NAME="ub"
   character (len = *), parameter :: VB_NAME="vb"
   character (len = *), parameter :: SSH_NAME="SSH"

   integer :: time_dimid, npoin_dimid, nlayers_dimid, parameter_dimid
   integer :: dt_varid, dt_btp_varid, zb_varid, ssh_varid
   integer :: ub_varid, vb_varid, h_varid
   integer :: zi_dimid, x_varid, y_varid
   real  :: q(4,npoin)

   q(1,:) = (alpha_mlswe(1)/gravity)*qb(1,:)
   q(2,:) = qb(3,:) / qb(1,:)
   q(3,:) = qb(4,:) / qb(1,:)
   q(4,:) = zbot_df + q(1,:)

   ! Gather Data onto Head node

   call gather_data(qb_g,q,4)
   call gather_data(coord_dg_gathered,coord,3)
   call gather_data(zbot_g,zbot_df,1)

   if (irank == irank0 .and. idone == 0) then

      ! Generate the name of the file to which the data will be written,
      ! and open the file.
      tempchar = trim(fname_root)
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

      call check(nf90_put_att(ncid, NF90_GLOBAL, "filename", trim(fn)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "npoin", "Number of points in the mesh"))

      ! Name each variable
      call check(nf90_def_var(ncid, DT_BTP_NAME, NF90_DOUBLE,1,dt_btp_varid))
      call check(nf90_put_att(ncid, dt_btp_varid, "name", "Barotropic time step"))
      call check(nf90_put_att(ncid, dt_btp_varid, "units", "seconds"))

      call check(nf90_def_var(ncid, X_NAME, NF90_DOUBLE,npoin_dimid,x_varid))
      call check(nf90_put_att(ncid, x_varid, "name", "cartesian coordinates"))
      call check(nf90_put_att(ncid, x_varid, "axis", "X"))

      call check(nf90_def_var(ncid, Y_NAME, NF90_DOUBLE,npoin_dimid,y_varid))
      call check(nf90_put_att(ncid, y_varid, "name", "cartesian coordinates"))
      call check(nf90_put_att(ncid, y_varid, "axis", "Y"))

      call check(nf90_def_var(ncid, H_NAME, NF90_DOUBLE,npoin_dimid,h_varid))
      call check(nf90_put_att(ncid, h_varid, "name", "Barotropic pressure pb"))
      call check(nf90_put_att(ncid, h_varid, "units", "m"))

      call check(nf90_def_var(ncid, UB_NAME, NF90_DOUBLE,npoin_dimid,ub_varid))
      call check(nf90_put_att(ncid, ub_varid, "name", "Barotropic u-velocity"))
      call check(nf90_put_att(ncid, ub_varid, "units", "m/s"))

      call check(nf90_def_var(ncid, VB_NAME, NF90_DOUBLE,npoin_dimid,vb_varid))
      call check(nf90_put_att(ncid, vb_varid, "name", "Barotropic v-velocitiy"))
      call check(nf90_put_att(ncid, vb_varid, "units", "m/s"))

      call check(nf90_def_var(ncid, SSH_NAME, NF90_DOUBLE,npoin_dimid,ssh_varid))
      call check(nf90_put_att(ncid, ssh_varid, "name", "SSH"))
      call check(nf90_put_att(ncid, ssh_varid, "units", "m"))

      call check(nf90_def_var(ncid, ZBOT_NAME, NF90_DOUBLE,npoin_dimid,zb_varid))
      call check(nf90_put_att(ncid, zb_varid, "name", "Total depdth at rest"))
      call check(nf90_put_att(ncid, zb_varid, "units", "m"))

      ! End define mode
      call check(nf90_enddef(ncid))

      ! Write data
      call check(nf90_put_var(ncid, dt_btp_varid, dt_btp))
      call check(nf90_put_var(ncid, x_varid, coord_dg_gathered(1,:)))
      call check(nf90_put_var(ncid, y_varid, coord_dg_gathered(2,:)))
      call check(nf90_put_var(ncid, h_varid, qb_g(1,:)))
      call check(nf90_put_var(ncid, ub_varid, qb_g(2,:)))
      call check(nf90_put_var(ncid, vb_varid, qb_g(3,:)))
      call check(nf90_put_var(ncid, ssh_varid, qb_g(4,:)))
      call check(nf90_put_var(ncid, zb_varid, zbot_g(:)))

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
end subroutine diagnostics_btp_nc
