!----------------------------------------------------------------------!
!>@brief Writes Output
!>@author  J.F. Kelly and F.X. Giraldo
!>Department of Applied Maths
!>Naval Postgraduate School
!>Monterey, CA 93943
!> @date March 12 2015, S. Marras added xml output from Michal's subroutine in the CG-only code
!----------------------------------------------------------------------!

subroutine write_output_mlswe(qp,qb,fnp1,time,layer)

    use mod_constants, only: gravity
  
      use mod_convert_variables, only: mod_convert_variables_eqnset_to_set2nc
      
      use mod_grid, only: npoin
  
      use mod_initial, only: nvar, q_ref, q_exact, rho_layers, q_ref_layers
  
      use mod_input, only: icase, fname_root, nopz, nelz, out_type, lout_ascii, lout_asciimaya, space_method, nvtk_files, lout_vorticity, lp4est, lp6est, &
          zlevel_out, lout_spherical_shell, lout_nc_3d, lincompressible, locean, lsalinity, write_mesh, is_swe_layers, equations
  
      use mod_p4est, only: mod_p4est_dump_mesh
  
      use mod_interface, only: compute_curl, compute_vorticity
  
      use mod_mpi_utilities, only: irank0, irank
    
      use mod_parallel, only: num_send_recv_total
  
      use mpi
  
      implicit none
  
      real, intent(in)  :: qp(nvar,npoin), qb(4,npoin)
      integer, intent(in) :: layer
      real, dimension(:,:), allocatable  :: q_aux, qq, vorticity
      real time
  
      character :: fnp1*9, fnp*100, fnp2*4
      integer   :: ilevel, iloop, j, i
  
      integer   :: nproc, ierr
  
      !allocate(qq(nvar,npoin), q_aux(nvar,npoin))
      allocate(qq(nvar,npoin), vorticity(3,npoin),q_aux(nvar,npoin))
      qq=qp

      fnp=trim(fname_root) // '_' // trim(fnp1) // '.vtk'
      call outvtk_g_binary_mlswe(qq,qb,fnp,time)
  
      deallocate(qq,q_aux,vorticity)
  
  end subroutine write_output_mlswe

  subroutine write_output_mlswe_global(qp,qb,qprime,fnp1,time,layer)

      use mod_mpi_utilities, only: irank0, irank
    
      use mod_parallel, only: num_send_recv_total

      use mod_input, only: icase, fname_root, nopz, nelz, out_type, lout_ascii, lout_asciimaya, space_method, nvtk_files, lout_vorticity, lp4est, lp6est, &
          zlevel_out, lout_spherical_shell, lout_nc_3d, lincompressible, locean, lsalinity, write_mesh, is_swe_layers, equations
  
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


