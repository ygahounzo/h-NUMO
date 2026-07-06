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

  use mod_time_loop,       only: time_loop, rhs_time, write_time
  use mod_types,           only: r8
  use mod_mpi_utilities,   only: initialize_mpi_util, finalize_mpi_util, irank, irank0, numproc, wtime
  use mod_p4est,           only: mod_p4est_finalize
  use mod_grid,            only: grid
  use mod_input,           only: input
  use mod_basis,           only: basis
  use mod_global_grid,     only: grid_global
  use mod_parallel,        only: parallel_CS
  use mod_face,            only: face_CS
  use mod_variables,       only: btp_CS, bcl_CS
  use mod_initial,         only: initial
  use mod_tensor,          only: tensor_CS
  use mod_metrics,         only: metrics
  use mod_ref,             only: mref
  use mod_mpi_communicator, only: mpi_communicator

  implicit none

  type(grid)             :: G
  type(input)            :: inp
  type(basis)            :: b
  type(grid_global)      :: gg
  type(parallel_CS)      :: par
  type(face_CS)          :: mf
  type(btp_CS)           :: btp
  type(bcl_CS)           :: bcl
  type(initial)          :: init
  type(tensor_CS)        :: tsp
  type(metrics)          :: mt
  type(mref)             :: ref
  type(mpi_communicator) :: mpic

  real(kind=r8) :: time1, time2, cpu_time
  integer :: ierr, flag

  ! Initialize MPI
  call initialize_mpi_util(ierr)

  print*, "MPI initialized"

  ! Initialize Grid and other non-variable data
  call initialize_grid(G, inp, b, gg, par, mf, init)

  ! Initialize field data
  call initialize_fields(G, inp, b, par, mf, btp, bcl, init, tsp, mt, ref, mpic)

  !------Print Out Input Data
  if (irank == irank0) then
     flag=0
     call print_header(gg, inp, flag, numproc)
  end if

  rhs_time = 0

  !   Perform Time-Integration
  time1 = wtime()

  call time_loop(G, inp, b, gg, par, mf, btp, bcl, init, ref, mpic, tsp, mt)

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
     call print_header(gg, inp, flag, numproc)
  end if

  ! ----- call p4est
  call mod_p4est_finalize()

  !Finalize MPI
  call finalize_mpi_util(ierr)

end program numo3d

!---------------------------------------------------
!>@brief Re-initialize data after refinement
!---------------------------------------------------
subroutine initialize_fields(G, inp, b, par, mf, btp, bcl, init, tsp, mt, ref, mpic)

  use mod_bc,              only: mod_bc_create
  use mod_mpi_communicator, only: mod_mpi_communicator_create, mpi_communicator
  use mod_face,            only: mod_face_create, mod_face_create_boundary, face_CS
  use mod_grid,            only: grid
  use mod_initial,         only: mod_initial_create, initial
  use mod_input,           only: input
  use mod_basis,           only: basis
  use mod_metrics,         only: mod_metrics_create_metrics, mod_metrics_create_mass, metrics
  use mod_mpi_utilities,   only: irank, irank0
  use mod_parallel,        only: parallel_CS, mod_parallel_reorder
  use mod_variables,       only: mod_allocate_mlswe, btp_CS, bcl_CS
  use mod_ref,             only: mod_ref_create, mref
  use mod_tensor,          only: compute_tensor_product, tensor_CS

  implicit none

  type(grid),            intent(inout) :: G
  type(input),           intent(inout) :: inp
  type(basis),           intent(in)    :: b
  type(parallel_CS),     intent(inout) :: par
  type(face_CS),         intent(inout) :: mf
  type(btp_CS),          intent(inout) :: btp
  type(bcl_CS),          intent(inout) :: bcl
  type(initial),         intent(inout) :: init
  type(tensor_CS),       intent(inout) :: tsp
  type(metrics),         intent(inout) :: mt
  type(mref),            intent(inout) :: ref
  type(mpi_communicator), intent(inout) :: mpic

  !-----Create FACE, METRICS, MASS Matrix-----!
  ! Create Face Data Structures on Processor
  call mod_face_create(G, inp, b, mf)
  if (irank == irank0) print *, "Faces Created"

  ! Create Metrics
  call mod_metrics_create_metrics(G, b, inp, mt)
  if (irank == irank0) print *, "Metrics Created"

  ! Create the Boundary of each Processor Element
  call mod_face_create_boundary(G, mf, G%nboun)
  if (irank == irank0) print *, "Boundary Created"

   if (par%nproc > 1) then
      call mod_parallel_reorder(par, mf%face_send)
   endif

   ! Tensor products
   call compute_tensor_product(G, b, mt, tsp)
   if (irank == irank0) print *, "Tensor Products Computed"

  ! Create Initial Conditions, Boundary Conditions, reference fields, and Metric terms

   ! allocate btp_variables
   call mod_allocate_mlswe(inp, G, b, btp, bcl)

  !Create Initial Conditions
  call mod_initial_create(inp, b, G, mf, init, mt)
  if (irank == irank0) print *, "Initial Fields Created"

  ! Create Boundary Condtions
  call mod_bc_create(G, b, mf, par)
  if (irank == irank0) print *, "Boundary Conditions Created"

  ! Create Mass Matrix
  call mod_metrics_create_mass(G, b, par, mt)
  if (irank == irank0) print *, "Mass Matrix Created"

  ! Create Reference Fields
  call mod_ref_create(G, inp, init, par, b, ref)
  if (irank == irank0) print *, "Reference Fields Created"

   !Create DG Communicator arrays
   if (inp%space_method == 'dg') then
      call mod_mpi_communicator_create(par, mpic)
   end if

end subroutine initialize_fields

!>@brief Create initial grid
subroutine initialize_grid(G, inp, b, gg, par, mf, init)

  use mod_basis,       only: mod_basis_create, basis
  use mod_constants,   only: mod_constants_create
  use mod_global_grid, only: grid_global
  use mod_grid,        only: grid
  use mod_input,       only: mod_input_create, input
  use mod_mpi_utilities, only: irank, irank0, task_sync
  use mod_p4est,       only: mod_p4est_create
  use mod_parallel,    only: parallel_CS
  use mod_face,        only: face_CS
  use mod_initial,     only: initial

  implicit none

  type(grid),        intent(inout) :: G
  type(input),       intent(inout) :: inp
  type(basis),       intent(inout) :: b
  type(grid_global), intent(inout) :: gg
  type(parallel_CS), intent(inout) :: par
  type(face_CS),     intent(inout) :: mf
  type(initial),     intent(inout) :: init

  ! Read NUMO3D.in
  call mod_input_create(inp)
  if (irank == irank0) then
     print *, "Input file Read"
  end if

  !Initialize Constants
  call mod_constants_create()
  if (irank == irank0) then
     print *, "Constants Created"
  end if

  !Create Basis Functions and Filters
  call mod_basis_create(inp, b)
  if (irank == irank0) then
     print *, "LGL Basis Created"
  end if

  !--- Read case specific defaults
  call initial_grid_coord(inp)

  ! sphere_ico: generate the coarse icosahedral base mesh as an ABAQUS
  ! .inp file that p4est reads via the existing read_external_grid_flg
  ! path (p4est derives all tree adjacency itself). Rank 0 writes it;
  ! everyone else waits so the file is complete before p4est reads it.
  if (trim(inp%geometry_type) == 'sphere_ico') then
     if (irank == irank0) then
        ! nelx plays the same role here as it does for sphere_hex
        ! (elements per base panel edge) -- an arbitrary integer, not
        ! restricted to a power of 2 like refinement_levels_h is.
        call create_grid_sphere_ico_inp('EXTERNAL_MESH.inp', inp%nelx)
        print *, "Icosahedral base mesh (EXTERNAL_MESH.inp) Created"
     end if
     call task_sync()
  end if

  ! Generate Global Grid and Domain Decomposition Graph

   call mod_p4est_create(G, inp, b, gg, par, mf, init)
   if (irank == irank0) print*,' P4EST Mesh Created'

end subroutine initialize_grid
