!----------------------------------------------------------------------!
!>@brief This subroutine builds the RHS PreCOMMUNICATOR for DG
!>@author  Francis X. Giraldo on 9/2015
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@update September 9, 2016 by F.X. Giraldo to always use optimal communication.
!>In this case here we use non-blocking sends-receives for DG Overlapping
!>computation and communication
!>@modified by   Yao Gahounzo
!>     Computing PhD
!>     Boise State University
!----------------------------------------------------------------------!

subroutine btp_create_precommunicator(G, inp, b, mf, init, par, ref, mpic, q, qprime_df, nvarb)

   use mod_grid,             only: grid
   use mod_input,            only: input
   use mod_basis,            only: basis
   use mod_face,             only: face_CS
   use mod_initial,          only: initial
   use mod_parallel,         only: parallel_CS
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator

   implicit none

   type(grid),             intent(in)    :: G
   type(input),            intent(in)    :: inp
   type(basis),            intent(in)    :: b
   type(face_CS),          intent(in)    :: mf
   type(initial),          intent(in)    :: init
   type(parallel_CS),      intent(in)    :: par
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic

   !Global Arrays
   integer, intent(in) :: nvarb
   real, dimension(nvarb, G%npoin),            intent(inout) :: q
   real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in) :: qprime_df

   ! DG - Discontinuous communicator

   call pack_and_send_df_btp(G, inp, b, mf, init, par, ref, ref%send_data_dg, ref%recv_data_dg, q, qprime_df, nvarb, mpic%nreq, mpic%ireq, mpic%status)

end subroutine btp_create_precommunicator

subroutine btp_lap_create_precommunicator(G, b, mf, init, par, btp, ref, mpic, q, nvarb)

   use mod_grid,             only: grid
   use mod_basis,            only: basis
   use mod_face,             only: face_CS
   use mod_initial,          only: initial
   use mod_parallel,         only: parallel_CS
   use mod_variables,        only: btp_CS
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator

   implicit none

   type(grid),             intent(in)    :: G
   type(basis),            intent(in)    :: b
   type(face_CS),          intent(in)    :: mf
   type(initial),          intent(in)    :: init
   type(parallel_CS),      intent(in)    :: par
   type(btp_CS),           intent(in)    :: btp
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic

   !Global Arrays
   integer, intent(in) :: nvarb
   real, dimension(nvarb, G%npoin), intent(inout) :: q

   ! pack_and_send_df_btp_lap only communicates the cartesian-equivalent 4 components
   ! (du_dx,du_dy,dv_dx,dv_dy); the sphere-only slots 5:9 of btp_dpp_graduvw aren't
   ! yet consumed by the viscosity Laplacian, so they don't need to cross MPI boundaries.
   call pack_and_send_df_btp_lap(G, b, mf, init, par, ref,                 &
      q, btp%btp_dpp_graduvw(1:4,:), btp%pbprime_visc, &
      nvarb, mpic%nreq, mpic%ireq, mpic%status)

end subroutine btp_lap_create_precommunicator

subroutine bcl_create_precommunicator(G, inp, b, mf, par, ref, mpic, qprime_df)

   use mod_grid,             only: grid
   use mod_input,            only: input
   use mod_basis,            only: basis
   use mod_face,             only: face_CS
   use mod_parallel,         only: parallel_CS
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator

   implicit none

   type(grid),             intent(in)    :: G
   type(input),            intent(in)    :: inp
   type(basis),            intent(in)    :: b
   type(face_CS),          intent(in)    :: mf
   type(parallel_CS),      intent(in)    :: par
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic

   !Global Arrays
   real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in) :: qprime_df

   ! DG - Discontinuous communicator

   call pack_and_send_df_bcl(G, inp, b, mf, par, ref, ref%send_data_bcl, ref%recv_data_bcl, qprime_df, mpic%nreq, mpic%ireq, mpic%status)

end subroutine bcl_create_precommunicator

subroutine bcl_lap_create_precommunicator(G, inp, b, mf, par, ref, mpic, dpp_graduv, dpprime_visc)

   use mod_grid,             only: grid
   use mod_input,            only: input
   use mod_basis,            only: basis
   use mod_face,             only: face_CS
   use mod_parallel,         only: parallel_CS
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator

   implicit none

   type(grid),             intent(in)    :: G
   type(input),            intent(in)    :: inp
   type(basis),            intent(in)    :: b
   type(face_CS),          intent(in)    :: mf
   type(parallel_CS),      intent(in)    :: par
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic

   !Global Arrays
   real, dimension(4, G%npoin, inp%nlayers), intent(in) :: dpp_graduv
   real, dimension(G%npoin, inp%nlayers),    intent(in) :: dpprime_visc

   ! GPU pack + GPU-direct MPI non-blocking send/receive.
   call pack_and_send_df_bcl_lap(G, inp, b, mf, par, ref, ref%send_data_lap_bcl, ref%recv_data_lap_bcl, &
      dpp_graduv, dpprime_visc, mpic%nreq, mpic%ireq, mpic%status)

end subroutine bcl_lap_create_precommunicator

subroutine btp_create_postcommunicator(G, inp, b, mf, par, btp, init, ref, mpic, rhs, nvarb)

   use mod_grid,             only: grid
   use mod_input,            only: input
   use mod_basis,            only: basis
   use mod_face,             only: face_CS
   use mod_parallel,         only: parallel_CS
   use mod_variables,        only: btp_CS
   use mod_initial,          only: initial
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator

   implicit none

   type(grid),             intent(in)    :: G
   type(input),            intent(in)    :: inp
   type(basis),            intent(in)    :: b
   type(face_CS),          intent(in)    :: mf
   type(parallel_CS),      intent(in)    :: par
   type(btp_CS),           intent(inout) :: btp
   type(initial),          intent(in)    :: init
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic

   !Global Arrays
   integer, intent(in) :: nvarb
   real, dimension(3, G%npoin), intent(inout) :: rhs

   ! DG - Discontinuous communicator

   !To build inter-processor fluxes, All Procs Must Wait
   call mpi_waitall(mpic%nreq, mpic%ireq, mpic%status, mpic%ierr)
   call unpack_data_dg_general_df(G, b, par, ref%q_send, ref%q_recv, ref%send_data_dg, ref%recv_data_dg, ref%nbtp_var, ref%nboun_valid)
   call create_nbhs_face_df(G, inp, b, mf, par, btp, init, ref, rhs)

end subroutine btp_create_postcommunicator

subroutine create_rhs_lap_postcommunicator_df(G, b, mf, par, btp, ref, mpic, rhs, nvarb)

   use mod_grid,             only: grid
   use mod_basis,            only: basis
   use mod_face,             only: face_CS
   use mod_parallel,         only: parallel_CS
   use mod_variables,        only: btp_CS
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator

   implicit none

   type(grid),             intent(in)    :: G
   type(basis),            intent(in)    :: b
   type(face_CS),          intent(in)    :: mf
   type(parallel_CS),      intent(in)    :: par
   type(btp_CS),           intent(inout) :: btp
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic

   !Global Arrays
   integer, intent(in) :: nvarb
   real, dimension(nvarb, G%npoin), intent(inout) :: rhs

   ! DG - Discontinuous communicator

   !To build inter-processor fluxes, All Procs Must Wait
   call mpi_waitall(mpic%nreq, mpic%ireq, mpic%status, mpic%ierr)

   !Map Recv buffer to the boundary of the Receiver (unpack data)
   call unpack_data_dg_general_lap(G, b, par, ref, ref%q_send_lap, ref%q_recv_lap, ref%send_data_dg_lap, &
      ref%recv_data_dg_lap, nvarb)

   !Build Inviscid Fluxes On Element Boundary
   call create_nbhs_face_df_lap(G, b, mf, par, btp, rhs, ref, nvarb)

end subroutine create_rhs_lap_postcommunicator_df

subroutine bcl_create_postcommunicator(G, inp, b, mf, par, btp, init, ref, mpic, rhs)

   use mod_grid,             only: grid
   use mod_input,            only: input
   use mod_basis,            only: basis
   use mod_face,             only: face_CS
   use mod_parallel,         only: parallel_CS
   use mod_variables,        only: btp_CS
   use mod_initial,          only: initial
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator

   implicit none

   type(grid),             intent(in)    :: G
   type(input),            intent(in)    :: inp
   type(basis),            intent(in)    :: b
   type(face_CS),          intent(in)    :: mf
   type(parallel_CS),      intent(in)    :: par
   type(btp_CS),           intent(inout) :: btp
   type(initial),          intent(in)    :: init
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic

   !Global Arrays
   real, dimension(3, G%npoin, inp%nlayers), intent(inout) :: rhs

   ! DG - Discontinuous communicator

   !To build inter-processor fluxes, All Procs Must Wait
   call mpi_waitall(mpic%nreq, mpic%ireq, mpic%status, mpic%ierr)

   !Map Recv buffer to the boundary of the Receiver (unpack data)
   call unpack_data_dg_general_bcl(G, b, inp, par, ref%q_send_bcl, ref%q_recv_bcl, ref%send_data_bcl, ref%recv_data_bcl, ref%nboun_valid)

   !Build Inviscid Fluxes On Element Boundary
   call create_nbhs_face_bcl(G, inp, b, mf, par, btp, init, ref, rhs, ref%q_send_bcl, ref%q_recv_bcl)

end subroutine bcl_create_postcommunicator

subroutine bcl_create_postcommunicator_continuity(G, inp, b, mf, par, btp, ref, mpic, rhs)

   use mod_grid,             only: grid
   use mod_input,            only: input
   use mod_basis,            only: basis
   use mod_face,             only: face_CS
   use mod_parallel,         only: parallel_CS
   use mod_variables,        only: btp_CS
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator

   implicit none

   type(grid),             intent(in)    :: G
   type(input),            intent(in)    :: inp
   type(basis),            intent(in)    :: b
   type(face_CS),          intent(in)    :: mf
   type(parallel_CS),      intent(in)    :: par
   type(btp_CS),           intent(inout) :: btp
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic

   !Global Arrays
   real, dimension(G%npoin, inp%nlayers), intent(inout) :: rhs

   ! DG - Discontinuous communicator

   !To build inter-processor fluxes, All Procs Must Wait
   call mpi_waitall(mpic%nreq, mpic%ireq, mpic%status, mpic%ierr)

   !Map Recv buffer to the boundary of the Receiver (unpack data)
   call unpack_data_dg_general_bcl(G, b, inp, par, ref%q_send_bcl, ref%q_recv_bcl, ref%send_data_bcl, ref%recv_data_bcl, ref%nboun_valid)

   !Build Inviscid Fluxes On Element Boundary
   call create_nbhs_face_bcl_continuity(G, inp, b, mf, par, btp, rhs, ref%q_send_bcl, ref%q_recv_bcl, 0)

end subroutine bcl_create_postcommunicator_continuity

subroutine bcl_create_postcommunicator_momentum(G, inp, b, mf, par, btp, init, ref, mpic, rhs)

   use mod_grid,             only: grid
   use mod_input,            only: input
   use mod_basis,            only: basis
   use mod_face,             only: face_CS
   use mod_parallel,         only: parallel_CS
   use mod_variables,        only: btp_CS
   use mod_initial,          only: initial
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator

   implicit none

   type(grid),             intent(in)    :: G
   type(input),            intent(in)    :: inp
   type(basis),            intent(in)    :: b
   type(face_CS),          intent(in)    :: mf
   type(parallel_CS),      intent(in)    :: par
   type(btp_CS),           intent(inout) :: btp
   type(initial),          intent(in)    :: init
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic

   !Global Arrays
   real, dimension(2, G%npoin, inp%nlayers), intent(inout) :: rhs

   ! DG - Discontinuous communicator

   !To build inter-processor fluxes, All Procs Must Wait
   call mpi_waitall(mpic%nreq, mpic%ireq, mpic%status, mpic%ierr)

   !Map Recv buffer to the boundary of the Receiver (unpack data)
   call unpack_data_dg_general_bcl(G, b, inp, par, ref%q_send_bcl, ref%q_recv_bcl, ref%send_data_bcl, ref%recv_data_bcl, ref%nboun_valid)

   !Build Inviscid Fluxes On Element Boundary
   call create_nbhs_face_bcl_momentum(G, inp, b, mf, par, btp, init, rhs, ref%q_send_bcl, ref%q_recv_bcl, 0)

end subroutine bcl_create_postcommunicator_momentum

subroutine bcl_create_rhs_lap_postcommunicator_df(G, inp, b, mf, par, btp, ref, mpic, rhs)

   use mod_grid,             only: grid
   use mod_input,            only: input
   use mod_basis,            only: basis
   use mod_face,             only: face_CS
   use mod_parallel,         only: parallel_CS
   use mod_variables,        only: btp_CS
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator

   implicit none

   type(grid),             intent(in)    :: G
   type(input),            intent(in)    :: inp
   type(basis),            intent(in)    :: b
   type(face_CS),          intent(in)    :: mf
   type(parallel_CS),      intent(in)    :: par
   type(btp_CS),           intent(inout) :: btp
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic

   !Global Arrays
   real, dimension(2, G%npoin, inp%nlayers), intent(inout) :: rhs

   ! DG - Discontinuous communicator

   !To build inter-processor fluxes, All Procs Must Wait
   call mpi_waitall(mpic%nreq, mpic%ireq, mpic%status, mpic%ierr)

   ! GPU unpack from flat MPI buffer into q_send/q_recv.
   call unpack_data_dg_general_lap_bcl(G, b, par, ref, ref%q_send_lap_bcl, ref%q_recv_lap_bcl, ref%send_data_lap_bcl, &
      ref%recv_data_lap_bcl, inp%nlayers, ref%nboun_valid)

   ! GPU face-gang scatter into rhs on device.
   call create_nbhs_face_df_lap_bcl(G, b, mf, par, btp, ref, rhs, ref%q_send_lap_bcl, ref%q_recv_lap_bcl, inp%nlayers, 0)

end subroutine bcl_create_rhs_lap_postcommunicator_df

subroutine bcl_create_communicator(G, b, par, ref, mpic, q_face, nvarb, nlayers, nq)

   use mod_grid,             only: grid
   use mod_basis,            only: basis
   use mod_parallel,         only: parallel_CS
   use mod_ref,              only: mref
   use mod_mpi_communicator, only: mpi_communicator

   implicit none

   type(grid),             intent(in)    :: G
   type(basis),            intent(in)    :: b
   type(parallel_CS),      intent(in)    :: par
   type(mref),             intent(inout) :: ref
   type(mpi_communicator), intent(inout) :: mpic

   !Global Arrays
   integer, intent(in) :: nvarb, nlayers, nq
   real, dimension(nvarb, 2, nq, G%nface, nlayers), intent(inout) :: q_face

   !MPI Variables
   integer :: i, j, k, iv, e, ip
   real :: recv_data_dg_quad1(nvarb*nq*G%nboun*nlayers)
   real :: send_data_dg_quad1(nvarb*nq*G%nboun*nlayers)
   real :: q_recv_quad1(nvarb, nq, G%nboun, nlayers)
   real :: q_send_quad1(nvarb, nq, G%nboun, nlayers)

   recv_data_dg_quad1 = 0.0
   send_data_dg_quad1 = 0.0
   q_recv_quad1 = 0.0
   q_send_quad1 = 0.0

   ! DG - Discontinuous communicator

   !Load all the boundary data into a vector
   call pack_data_dg_quad_layer(send_data_dg_quad1, q_face, nvarb, nlayers, nq)

   !non-blocking sends-receives: message size=nmessage
   call send_bound_dg_general_quad_layer(G, b, par, send_data_dg_quad1, recv_data_dg_quad1, nvarb, &
      nlayers, nq, mpic%nreq, mpic%ireq, mpic%status)

   !To build inter-processor fluxes, All Procs Must Wait
   call mpi_waitall(mpic%nreq, mpic%ireq, mpic%status, mpic%ierr)

   !Map Recv buffer to the boundary of the Receiver (unpack data)
   call unpack_data_dg_general_quad_layer(G, par, q_send_quad1, q_recv_quad1, send_data_dg_quad1, &
      recv_data_dg_quad1, nvarb, nlayers, nq)

   !Build Inviscid Fluxes On Element Boundary
   call create_nbhs_face_quad_layer(G, b, par, q_face, q_send_quad1, q_recv_quad1, nvarb, nlayers, nq)

end subroutine bcl_create_communicator
