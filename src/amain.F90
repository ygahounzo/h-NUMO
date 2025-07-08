!=============================================================================================
!>@brief This code solves the multilayer shallow water equations (2D in the horizontal and Nl 
!> number of layers in the vertical) using DG
!> using exact integration and the following time-integrator: RK35
!> The grids are quadrilateral but can be unstructured.
!>@ author by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: July 02, 2023
!=============================================================================================

program numo3d

  use mod_time_loop, only: time_loop, rhs_time, write_time
  use mod_types, only : r8
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

  !   Perform Time-Integration
  time1 = wtime()

  call time_loop()

  time2 = wtime()
  cpu_time=time2-time1

  !   Finalize
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

  ! ----- call p4est
  call mod_p4est_finalize()

  !Finalize MPI
  call finalize_mpi_util(ierr)

end program numo3d

!---------------------------------------------------
!>@brief Re-initialize data after refinement
!---------------------------------------------------
subroutine initialize_fields()

  use mod_bc, only: mod_bc_create
  use mod_mpi_communicator, only: mod_mpi_communicator_create
  use mod_face, only: mod_face_create, mod_face_create_boundary, face_send
  use mod_grid, only: nboun
  use mod_initial, only: mod_initial_create, q_init, nvar
  use mod_input, only: nopx, nopy, nopz, space_method
  use mod_metrics, only: mod_metrics_create_metrics, mod_metrics_create_mass
  use mod_mpi_utilities
  use mod_parallel, only: nproc, mod_parallel_reorder
  use mod_variables, only: mod_allocate_mlswe
  use mod_ref, only: mod_ref_create

  implicit none

  !local variables
  logical is_p4est

  integer :: ivar

  ! Is p4est
  !is_p4est = .true.

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

  ! Create Initial Conditions, Boundary Conditions, reference fields, and Metric terms

  !Create Initial Conditions
  call mod_initial_create()
  if (irank == irank0) print *, "Initial Fields Created"

  ! Create Boundary Condtions
  call mod_bc_create() 
  if (irank == irank0) print *, "Boundary Conditions Created"

  ! Create Mass Matrix
  call mod_metrics_create_mass()
  if (irank == irank0) print *, "Mass Matrix Created"

  ! Create Reference Fields
  call mod_ref_create()
  if (irank == irank0) print *, "Reference Fields Created"

   !Create DG Communicator arrays
   if (space_method == 'dg') then
      call mod_mpi_communicator_create()
   end if

   ! allocate btp_variables

   call mod_allocate_mlswe()

end subroutine initialize_fields

!>@brief Create initial grid
subroutine initialize_grid()

  use mod_basis, only: mod_basis_create
  use mod_constants, only: mod_constants_create
  use mod_grid, only: npoin, nelem, nboun, face, nface
  use mod_input, only: mod_input_create, nopx, nopy, nopz, &
                        space_method
  use mod_mpi_utilities
  use mod_p4est, only: mod_p4est_create

  implicit none

  !local variables
  logical is_p4est

  ! Read NUMO3D.in
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
  call initial_grid_coord()

  ! Generate Global Grid and Domain Decomposition Graph

   call mod_p4est_create()
   if (irank == irank0) print*,' P4EST Mesh Created'

end subroutine initialize_grid
