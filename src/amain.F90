!-----------------------------------------------------------------!
!>@brief This code solves the 3D Euler Equations using CG and DG
!>@details It uses SET 2NC, 2C and 3C as described in Giraldo et al. SISC 2010 
!> using SE w/inexact integration and the following time-integrators: 
!> 1) explicit SSP RK2, RK3, RK34, RK35,
!> 2) IMEX 1D & 3D BDF2 (Schur and No Schur forms), 
!> 3) IMEX 1D & 3D LF2  (Schur and No Schur forms), and
!> 4) IMEX 1D & 3D ARK methods (Schur and No Schur Forms).
!> 5) IMEX 1D & 3D AM2, BDF23 methods (Schur and No Schur Forms).
!>The grids are Hexahedra but can be unstructured.
!>@author  F.X. Giraldo and J.F. Kelly 
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@date 1 August 2010 James F. Kelly
!> NUMA3dCG: A parallel implementation of NUMA for both spherical and 
!> Cartesian grids.  The spatial discretization is based on spectral elements
!> (or continuous Galerkin).  Domain decomposition utilizes either geometric
!> decomposition (e.g. NSEAM-style cube grid decomposition) or METIS-based
!> decomposition.
!>@date January 2014 F.X. Giraldo
!> Cleaned to remove redundant arrays, rename arrays to better names, and 
!> modified to handle Periodic Boundary conditions by
!>@date Dec 2014 Daniel S. Abdi
!> Modified to merge parallel CG/DG
!-----------------------------------------------------------------!
program numa3d

  use mod_time_loop_mlswe, only: time_loop, rhs_time, write_time

  use mod_types, only : r8

  use mod_input, only: lp4est, lp6est

  use mod_mpi_utilities

  use mod_p4est, only: mod_p4est_finalize

  implicit none

  !local variables
  real(kind=r8) :: time1, time2, cpu_time
  integer :: ierr, flag

  ! Initialize MPI
  call initialize_mpi_util(ierr)

  print*, "MPI initialized"

  ! Initialize Grid and other non-variable data
  call initialize_grid()

  ! Initialize field data
  call initialize_fields()

  !------Print Out Input Data
  if (irank == irank0) then
     flag=0
     call print_header(flag,numproc)
  end if
  
  rhs_time = 0
  !--------------------------------
  !   Perform Time-Integration
  !--------------------------------
  time1 = wtime()

  call time_loop()

  time2 = wtime()
  cpu_time=time2-time1

  !---------------------------
  !   Finalize
  !---------------------------
  if (irank == irank0) then
     write(*,'("Wall Clock Time (sec) : ",F12.4)') cpu_time
     write(*,'("RHS  CPU   Time (sec) : ",F12.4)') rhs_time

     call write_time()

  endif

  !-----Print Out Input Data
  if (irank == irank0) then
     flag=1
     call print_header(flag,numproc)
  end if

  if (lp4est .or. lp6est) call mod_p4est_finalize()

  !Finalize MPI
  call finalize_mpi_util(ierr)

end program numa3d

!---------------------------------------------------
!>@brief Re-initialize data after refinement
!---------------------------------------------------
subroutine initialize_fields()

  use mod_array_util, only: create_arrays, get_dimensions

  use mod_bc, only: mod_bc_create

  use mod_mpi_communicator, only: mod_mpi_communicator_create

  use mod_face, only: mod_face_create, mod_face_create_boundary, face_send

  use mod_filter, only: mod_filter_create_rhs

  use mod_grid, only: nboun

  use mod_initial, only: mod_initial_create, q_init, nvar

  use mod_input, only: &
       nopx, nopy, nopz, &
       space_method, &
       lp4est, lp6est

  use mod_metrics, only: mod_metrics_create_metrics, mod_metrics_create_mass

  use mod_mpi_utilities

  use mod_parallel, only: nproc, mod_parallel_reorder

  implicit none

  !local variables
  logical is_p4est

  integer :: ivar

  ! Is p4est
  is_p4est = lp4est .or. lp6est

  !-----Create FACE, METRICS, MASS Matrix-----!
  ! Create Face Data Structures on Processor
  call mod_face_create()
  if (irank == irank0) print *, "Faces Created"

  ! Create Metrics
  call mod_metrics_create_metrics()
  if (irank == irank0) print *, "Metrics Created"

  ! Create the Boundary of each Processor Element
  call mod_face_create_boundary(nboun)
  if (irank == irank0) print *, "Boundary Created"

   if (nproc > 1) then
      call mod_parallel_reorder(face_send)
   endif

  !----------------------------------------------------------------------------------------
  !    Create Initial Conditions, Boundary Conditions, reference fields, and Metric terms
  !----------------------------------------------------------------------------------------
  !Create Initial Conditions
  call mod_initial_create()
  if (irank == irank0) print *, "Initial Fields Created"

  ! Create Boundary Condtions
  call mod_bc_create() 
  if (irank == irank0) print *, "Boundary Conditions Created"

  ! Create Mass Matrix
  call mod_metrics_create_mass()
  if (irank == irank0) print *, "Mass Matrix Created"

  !--------------------------------------------------------------------------------
  !    Create Iterative Solver Arrays
  !--------------------------------------------------------------------------------

  !CPU only code

   !Create RHS Vectors for the Filter
   call mod_filter_create_rhs()
   if (irank == irank0) print *, "Filter RHS Created"

   ! Initialize local 2D and 3D arrays
   call get_dimensions()
   if (irank == irank0) print *, "GET_DIMENSIONS Created"

   call create_arrays()
   if (irank == irank0) print *, "CREATE_ARRAYS Created"

   !Create DG Communicator arrays
   if (space_method == 'dg') then
      call mod_mpi_communicator_create()
   end if

end subroutine initialize_fields

!---------------------------------------------------
!>@brief Create initial grid
!---------------------------------------------------
subroutine initialize_grid()

  use mod_basis, only: mod_basis_create

  use mod_constants, only: mod_constants_create

  use mod_grid, only: npoin, nelem, nboun, face, nface

  !use mod_initial, only: q_init, nvar

  use mod_input, only: mod_input_create, nopx, nopy, nopz, &
      space_method, lp4est, lp6est, icase

  use mod_mpi_utilities

  use mod_p4est, only: mod_p4est_create

  implicit none

  !local variables
  logical is_p4est

  ! Read NUMA3D.in
  call mod_input_create(irank)
  if (irank == irank0) then
     print *, "Input file Read"
  end if

  !Initialize Constants
  call mod_constants_create()
  if (irank == irank0) then
     print *, "Constants Created"
  end if

  !Create Basis Functions and Filters
  call mod_basis_create(nopx,nopy,nopz)
  if (irank == irank0) then
     print *, "LGL Basis Created"
  end if
  
  !--- Read case specific defaults
  call initial_grid_cube_shallow()

  !-------------------------------------------------------------
  !     Generate Global Grid and Domain Decomposition Graph
  !-------------------------------------------------------------
  is_p4est = lp4est .or. lp6est

   call mod_p4est_create()
   if (irank == irank0) print*,' P4EST Mesh Created'

end subroutine initialize_grid
