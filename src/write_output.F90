!----------------------------------------------------------------------!
!>@brief Writes Output
!>@author  J.F. Kelly and F.X. Giraldo
!>Department of Applied Maths
!>Naval Postgraduate School
!>Monterey, CA 93943
!> @date March 12 2015, S. Marras added xml output from Michal's subroutine in the CG-only code
!----------------------------------------------------------------------!

subroutine write_output_mlswe(q_df,qb,fnp1,time,layer,qout)

    use mod_constants, only: gravity

      use mod_grid, only: npoin

      use mod_initial, only: nvar, alpha_mlswe, zbot_df

      use mod_input, only: fname_root, nlayers

      use mod_p4est, only: mod_p4est_dump_mesh

      use mod_mpi_utilities, only: irank0, irank

      use mod_parallel, only: num_send_recv_total

      use mpi

      implicit none

      real, intent(in)  :: q_df(3,npoin,nlayers), qb(4,npoin)
      integer, intent(in) :: layer
      real, intent(out) :: qout(nvar,npoin)
      real, dimension(npoin,nlayers+1) :: mslwe_elevation
      real time

      character :: fnp1*9, fnp*100, fnp2*4
      integer   :: ilevel, iloop, j, i, k

      integer   :: nproc, ierr

      ! Convert prognostic layer variables to physical output variables,
      ! matching the transformation done in diagnostics()
      qout(1,:) = (alpha_mlswe(layer)/gravity)*q_df(1,:,layer)
      qout(2,:) = q_df(2,:,layer) / q_df(1,:,layer)
      qout(3,:) = q_df(3,:,layer) / q_df(1,:,layer)
      qout(4,:) = q_df(1,:,layer)

      mslwe_elevation = 0.0
      mslwe_elevation(:,nlayers+1) = zbot_df
      do k = nlayers,1,-1
          mslwe_elevation(:,k) = mslwe_elevation(:,k+1) + (alpha_mlswe(k)/gravity)*q_df(1,:,k)
      end do
      qout(5,:) = mslwe_elevation(:,layer)

      fnp=trim(fname_root) // '_' // trim(fnp1) // '.vtk'
      call outvtk_g_binary_mlswe(qout,qb,fnp,time)

  end subroutine write_output_mlswe

  subroutine write_output_mlswe_global(qp,qb,qprime,fnp1,time,layer)

      use mod_mpi_utilities, only: irank0, irank
    
      use mod_parallel, only: num_send_recv_total

      use mod_input, only: fname_root
  
      use mod_grid, only: npoin
  
      use mpi
  
      implicit none
  
      real, intent(in)  :: qp(3,npoin), qb(3,npoin), qprime(3,npoin)
      integer, intent(in) :: layer
      real, dimension(:,:), allocatable  :: qq
      real time
  
      character :: fnp1*9, fnp*100, fnp2*4
      integer   :: ilevel, iloop, j, i
  
      integer   :: nproc, ierr

      allocate(qq(3,npoin))
      qq=qp

      fnp=trim(fname_root) // '_' // trim(fnp1) // '.vtk'
      call outvtk_g_binary_mlswe_global(qq,qb,qprime,fnp,time)
  
      deallocate(qq)
  
  end subroutine write_output_mlswe_global

!----------------------------!
!> @brief Timestamp
!----------------------------!
subroutine timestamp (stamp)

    implicit none

    character * (12) stamp
    character * ( 8 ) date
    character * ( 10 ) time
    call date_and_time ( date, time )
    time = time(1:4)
    stamp = date // time
end subroutine timestamp


