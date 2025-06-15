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

subroutine create_rhs_precommunicator_quad(q_face,nvarb)

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, nface

    use mod_initial, only: nvar

    use mod_ref, only: recv_data_dg_quad, send_data_dg_quad, nmessage

    use mod_basis, only: nq

    implicit none

    !Global Arrays
    real, dimension(nvarb,2,nq,nface), intent(in)  :: q_face
    integer, intent(in) :: nvarb

    !-------------------------------
    !     For DG
    !-------------------------------
    !Load all the boundary data into a vector
    call pack_data_dg_quad(send_data_dg_quad,q_face,nvarb)

    !non-blocking sends-receives: message size=nmessage
    call send_bound_dg_general_quad(send_data_dg_quad,recv_data_dg_quad,nvarb,nreq,ireq,status)

end subroutine create_rhs_precommunicator_quad

subroutine create_rhs_postcommunicator_quad(q_face,nvarb)

    use mod_basis, only: nq

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send_quad, q_recv_quad, recv_data_dg_quad, send_data_dg_quad, nmessage

    implicit none

    !Global Arrays
    real, dimension(nvarb,2,nq,nface), intent(inout) :: q_face
    integer, intent(in) :: nvarb

    integer :: multirate

    !MPI Variables
    integer :: i,j,k,iv,e,ip

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !To build inter-processor fluxes, All Procs Must Wait
    call mpi_waitall(nreq,ireq,status,ierr)

    !Map Recv buffer to the boundary of the Receiver (unpack data)
    call unpack_data_dg_general_quad(q_send_quad,q_recv_quad,send_data_dg_quad,recv_data_dg_quad,nvarb)

    !Build Inviscid Fluxes On Element Boundary - need to add multirate here
    call create_nbhs_face_quad(q_face,q_send_quad,q_recv_quad,nvarb,0)

end subroutine create_rhs_postcommunicator_quad

subroutine create_communicator_quad(q_face,nvarb)

    use mod_basis, only: nq

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send_quad, q_recv_quad, recv_data_dg_quad, send_data_dg_quad, nmessage

    implicit none

    !Global Arrays
    real, dimension(nvarb,2,nq,nface), intent(inout) :: q_face
    integer, intent(in) :: nvarb

    integer :: multirate

    !MPI Variables
    integer :: i,j,k,iv,e,ip
    real :: recv_data_dg_quad1(nvarb*nq*nboun)
    real :: send_data_dg_quad1(nvarb*nq*nboun)
    real :: q_recv_quad1(nvarb,nq,nboun), q_send_quad1(nvarb,nq,nboun)

    recv_data_dg_quad1 = 0.0
    send_data_dg_quad1 = 0.0
    q_recv_quad1 = 0.0
    q_send_quad1 = 0.0

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !Load all the boundary data into a vector
    call pack_data_dg_quad(send_data_dg_quad1,q_face,nvarb)

    !non-blocking sends-receives: message size=nmessage
    call send_bound_dg_general_quad(send_data_dg_quad1,recv_data_dg_quad1,nvarb,nreq,ireq,status)

    !To build inter-processor fluxes, All Procs Must Wait
    call mpi_waitall(nreq,ireq,status,ierr)

    !Map Recv buffer to the boundary of the Receiver (unpack data)
    call unpack_data_dg_general_quad(q_send_quad1,q_recv_quad1,send_data_dg_quad1,recv_data_dg_quad1,nvarb)

    !Build Inviscid Fluxes On Element Boundary - need to add multirate here
    call create_nbhs_face_quad(q_face,q_send_quad1,q_recv_quad1,nvarb,0)

end subroutine create_communicator_quad

subroutine create_communicator_df_1var(q0_df_face)

    use mod_basis, only: ngl

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send_quad, q_recv_quad, recv_data_dg_quad, send_data_dg_quad, nmessage

    implicit none

    !Global Arrays
    real, dimension(2,ngl,nface), intent(inout) :: q0_df_face

    integer :: multirate
    real, dimension(1,2,ngl,nface) :: q_df_face

    !MPI Variables
    integer :: i,j,k,iv,e,ip
    integer, parameter :: nvarb = 1
    real :: recv_data_dg_df1(nvarb*ngl*nboun)
    real :: send_data_dg_df1(nvarb*ngl*nboun)
    real :: q_recv_df1(nvarb,ngl,nboun), q_send_df1(nvarb,ngl,nboun)

    recv_data_dg_df1 = 0.0
    send_data_dg_df1 = 0.0
    q_recv_df1 = 0.0
    q_send_df1 = 0.0

    q_df_face(1,:,:,:) = q0_df_face(:,:,:)

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !Load all the boundary data into a vector
    call pack_data_dg_df(send_data_dg_df1,q_df_face,nvarb)

    !non-blocking sends-receives: message size=nmessage
    call send_bound_dg_general_df(send_data_dg_df1,recv_data_dg_df1,nvarb,nreq,ireq,status)

    !To build inter-processor fluxes, All Procs Must Wait
    call mpi_waitall(nreq,ireq,status,ierr)

    !Map Recv buffer to the boundary of the Receiver (unpack data)
    call unpack_data_dg_general_df(q_send_df1,q_recv_df1,send_data_dg_df1,recv_data_dg_df1,nvarb)

    !Build Inviscid Fluxes On Element Boundary - need to add multirate here
    call create_nbhs_face_df(q_df_face,q_send_df1,q_recv_df1,nvarb,0)

end subroutine create_communicator_df_1var

subroutine btp_create_precommunicator(q_df_face,nvarb)

    use mod_basis, only: ngl

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send, q_recv, recv_data_dg, send_data_dg

    implicit none

    !Global Arrays
    real, dimension(nvarb,2,ngl,nface), intent(inout) :: q_df_face
    integer, intent(in) :: nvarb

    integer :: multirate

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !Load all the boundary data into a vector
    call pack_data_dg_df(send_data_dg,q_df_face,nvarb)

    !non-blocking sends-receives: message size=nmessage
    call send_bound_dg_general_df(send_data_dg,recv_data_dg,nvarb,nreq,ireq,status)

end subroutine btp_create_precommunicator

subroutine btp_create_postcommunicator(q_df_face,nvarb)

    use mod_basis, only: ngl

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send, q_recv, recv_data_dg, send_data_dg

    implicit none

    !Global Arrays
    real, dimension(nvarb,2,ngl,nface), intent(inout) :: q_df_face
    integer, intent(in) :: nvarb

    integer :: multirate

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !To build inter-processor fluxes, All Procs Must Wait
    call mpi_waitall(nreq,ireq,status,ierr)

    !Map Recv buffer to the boundary of the Receiver (unpack data)
    call unpack_data_dg_general_df(q_send,q_recv,send_data_dg,recv_data_dg,nvarb)

    !Build Inviscid Fluxes On Element Boundary - need to add multirate here
    call create_nbhs_face_df(q_df_face,q_send,q_recv,nvarb,0)

end subroutine btp_create_postcommunicator


subroutine create_rhs_lap_precommunicator_df(q_df_face,nvarb)

    use mod_basis, only: ngl

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun

    use mod_initial, only: nvar

    use mod_ref, only: lap_recv_data_dg_df1, lap_send_data_dg_df1, lap_q_recv_df1, lap_q_send_df1, nmessage

    implicit none

    !Global Arrays
    real, dimension(nvarb,2,ngl,nface), intent(in) :: q_df_face
    integer, intent(in) :: nvarb

    integer :: multirate

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !Load all the boundary data into a vector
    call pack_data_dg_df(lap_send_data_dg_df1,q_df_face,nvarb)

    !non-blocking sends-receives: message size=nmessage
    call send_bound_dg_general_df(lap_send_data_dg_df1,lap_recv_data_dg_df1,nvarb,nreq,ireq,status)

end subroutine create_rhs_lap_precommunicator_df

subroutine create_rhs_lap_postcommunicator_df(q_df_face,nvarb)

    use mod_basis, only: ngl

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  nface,nboun

    use mod_ref, only: lap_recv_data_dg_df1, lap_send_data_dg_df1, lap_q_recv_df1, lap_q_send_df1, nmessage

    implicit none

    !Global Arrays
    real, dimension(nvarb,2,ngl,nface), intent(inout) :: q_df_face
    integer, intent(in) :: nvarb

    integer :: multirate

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !To build inter-processor fluxes, All Procs Must Wait
    call mpi_waitall(nreq,ireq,status,ierr)

    !Map Recv buffer to the boundary of the Receiver (unpack data)
    call unpack_data_dg_general_df(lap_q_send_df1,lap_q_recv_df1,lap_send_data_dg_df1,lap_recv_data_dg_df1,nvarb)

    !Build Inviscid Fluxes On Element Boundary - need to add multirate here
    call create_nbhs_face_df(q_df_face,lap_q_send_df1,lap_q_recv_df1,nvarb,0)

end subroutine create_rhs_lap_postcommunicator_df

subroutine create_communicator_quad_all(q_face,grad_uvdp_face,nvarb)

    use mod_basis, only: nq

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send_quad, q_recv_quad, recv_data_dg_quad, send_data_dg_quad, nmessage

    implicit none

    !Global Arrays
    real, dimension(nvarb,2,nq,nface), intent(inout) :: q_face
    real, dimension(nvarb,2,nq,nface), intent(inout) :: grad_uvdp_face
    integer, intent(in) :: nvarb


    !MPI Variables
    integer :: i,j,k,iv,e,ip
    real :: recv_data_dg_quad1(2*nvarb*nq*nboun)
    real :: send_data_dg_quad1(2*nvarb*nq*nboun)
    real :: q_recv_quad1(2*nvarb,nq,nboun), q_send_quad1(2*nvarb,nq,nboun)

    recv_data_dg_quad1 = 0.0
    send_data_dg_quad1 = 0.0
    q_recv_quad1 = 0.0
    q_send_quad1 = 0.0

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !Load all the boundary data into a vector
    call pack_data_dg_quad_all(send_data_dg_quad1,q_face,grad_uvdp_face,nvarb)

    !non-blocking sends-receives: message size=nmessage
    call send_bound_dg_general_quad(send_data_dg_quad1,recv_data_dg_quad1,2*nvarb,nreq,ireq,status)

    !To build inter-processor fluxes, All Procs Must Wait
    call mpi_waitall(nreq,ireq,status,ierr)

    !Map Recv buffer to the boundary of the Receiver (unpack data)
    call unpack_data_dg_general_quad(q_send_quad1,q_recv_quad1,send_data_dg_quad1,recv_data_dg_quad1,2*nvarb)

    !Build Inviscid Fluxes On Element Boundary - need to add multirate here
    call create_nbhs_face_quad_all(q_face,grad_uvdp_face,q_send_quad1,q_recv_quad1,nvarb,0)

end subroutine create_communicator_quad_all

subroutine bcl_create_communicator(q_face,nvarb,nlayers,nq)

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send_quad, q_recv_quad, recv_data_dg_quad, send_data_dg_quad, nmessage

    !use mod_input, only: nlayers

    implicit none

    !Global Arrays
    real, dimension(nvarb,2,nq,nface,nlayers), intent(inout) :: q_face
    integer, intent(in) :: nvarb, nlayers,nq

    integer :: multirate

    !MPI Variables
    integer :: i,j,k,iv,e,ip
    real :: recv_data_dg_quad1(nvarb*nq*nboun*nlayers)
    real :: send_data_dg_quad1(nvarb*nq*nboun*nlayers)
    real :: q_recv_quad1(nvarb,nq,nboun,nlayers), q_send_quad1(nvarb,nq,nboun,nlayers)

    recv_data_dg_quad1 = 0.0
    send_data_dg_quad1 = 0.0
    q_recv_quad1 = 0.0
    q_send_quad1 = 0.0

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !Load all the boundary data into a vector
    call pack_data_dg_quad_layer(send_data_dg_quad1,q_face,nvarb,nlayers,nq)

    !non-blocking sends-receives: message size=nmessage
    call send_bound_dg_general_quad_layer(send_data_dg_quad1,recv_data_dg_quad1,nvarb,nlayers,nq,nreq,ireq,status)

    !To build inter-processor fluxes, All Procs Must Wait
    call mpi_waitall(nreq,ireq,status,ierr)

    !Map Recv buffer to the boundary of the Receiver (unpack data)
    call unpack_data_dg_general_quad_layer(q_send_quad1,q_recv_quad1,send_data_dg_quad1,recv_data_dg_quad1,nvarb,nlayers,nq)

    !Build Inviscid Fluxes On Element Boundary - need to add multirate here
    call create_nbhs_face_quad_layer(q_face,q_send_quad1,q_recv_quad1,nvarb,nlayers,nq)

end subroutine bcl_create_communicator

subroutine create_communicator_quad_layer_all(qprime_face,q_face,nvarb,nlayers)

    use mod_basis, only: nq

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send_quad, q_recv_quad, recv_data_dg_quad, send_data_dg_quad, nmessage

    !use mod_input, only: nlayers

    implicit none

    !Global Arrays
    real, dimension(nvarb,2,nq,nface,nlayers), intent(inout) :: q_face, qprime_face
    integer, intent(in) :: nvarb, nlayers

    integer :: multirate

    !MPI Variables
    integer :: i,j,k,iv,e,ip, nvv
    real :: recv_data_dg_quad1(2*nvarb*nq*nboun*nlayers)
    real :: send_data_dg_quad1(2*nvarb*nq*nboun*nlayers)
    real :: q_recv_quad1(2*nvarb,nq,nboun,nlayers), q_send_quad1(2*nvarb,nq,nboun,nlayers)

    recv_data_dg_quad1 = 0.0
    send_data_dg_quad1 = 0.0
    q_recv_quad1 = 0.0
    q_send_quad1 = 0.0

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !Load all the boundary data into a vector
    call pack_data_dg_quad_layer_all(send_data_dg_quad1,q_face,qprime_face,nvarb,nlayers)
    
    !non-blocking sends-receives: message size=nmessage
    call send_bound_dg_general_quad_layer(send_data_dg_quad1,recv_data_dg_quad1,2*nvarb,nlayers,nq,nreq,ireq,status)

    !To build inter-processor fluxes, All Procs Must Wait
    call mpi_waitall(nreq,ireq,status,ierr)

    !Map Recv buffer to the boundary of the Receiver (unpack data)
    call unpack_data_dg_general_quad_layer(q_send_quad1,q_recv_quad1,send_data_dg_quad1,recv_data_dg_quad1,2*nvarb,nlayers,nq)

    !Build Inviscid Fluxes On Element Boundary - need to add multirate here
    call create_nbhs_face_quad_layer_all(q_face,qprime_face,q_send_quad1,q_recv_quad1,nvarb,nlayers)

end subroutine create_communicator_quad_layer_all

subroutine create_communicator_quad_1var(q_face)

    use mod_basis, only: nq

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send_quad, q_recv_quad, recv_data_dg_quad, send_data_dg_quad, nmessage

    implicit none

    !Global Arrays
    real, dimension(2,nq,nface), intent(inout) :: q_face

    integer :: multirate

    !MPI Variables
    integer :: i,j,k,iv,e,ip
    real :: recv_data_dg_quad1(nq*nboun)
    real :: send_data_dg_quad1(nq*nboun)
    real :: q_recv_quad1(nq,nboun), q_send_quad1(nq,nboun)

    recv_data_dg_quad1 = 0.0
    send_data_dg_quad1 = 0.0
    q_recv_quad1 = 0.0
    q_send_quad1 = 0.0

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !Load all the boundary data into a vector
    call pack_data_dg_quad_1v(send_data_dg_quad1,q_face)

    !non-blocking sends-receives: message size=nmessage
    call send_bound_dg_general_quad(send_data_dg_quad1,recv_data_dg_quad1,1,nreq,ireq,status)

    !To build inter-processor fluxes, All Procs Must Wait
    call mpi_waitall(nreq,ireq,status,ierr)

    !Map Recv buffer to the boundary of the Receiver (unpack data)
    call unpack_data_dg_general_quad_1v(q_send_quad1,q_recv_quad1,send_data_dg_quad1,recv_data_dg_quad1)

    !Build Inviscid Fluxes On Element Boundary - need to add multirate here
    call create_nbhs_face_quad_1v(q_face,q_send_quad1,q_recv_quad1,0)

end subroutine create_communicator_quad_1var


subroutine create_lap_postcommunicator_quad(rhs,nvarb)

    use mod_basis, only: nq

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send_quad, q_recv_quad, recv_data_dg_quad, send_data_dg_quad, nmessage

    implicit none

    !Global Arrays
    real, dimension(2,npoin), intent(inout) :: rhs
    integer, intent(in) :: nvarb

    integer :: multirate

    !MPI Variables
    integer :: i,j,k,iv,e,ip
    real :: recv_data_dg_quad1(nvarb*nq*nboun)
    real :: send_data_dg_quad1(nvarb*nq*nboun)
    real :: q_recv_quad1(nvarb,nq,nboun), q_send_quad1(nvarb,nq,nboun)

    recv_data_dg_quad1 = 0.0
    send_data_dg_quad1 = 0.0
    q_recv_quad1 = 0.0
    q_send_quad1 = 0.0

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !To build inter-processor fluxes, All Procs Must Wait
    call mpi_waitall(nreq,ireq,status,ierr)

    !Map Recv buffer to the boundary of the Receiver (unpack data)
    call unpack_data_dg_general_quad(q_send_quad1,q_recv_quad1,send_data_dg_quad1,recv_data_dg_quad1,nvarb)

    !Build Inviscid Fluxes On Element Boundary - need to add multirate here
    call create_nbhs_face_lap_quad_ip(rhs,q_send_quad1,q_recv_quad1,nvarb,0)

end subroutine create_lap_postcommunicator_quad

subroutine create_lap_precommunicator_quad(q_face,nvarb)

    use mod_basis, only: nq

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send_quad, q_recv_quad, recv_data_dg_quad, send_data_dg_quad, nmessage

    implicit none

    !Global Arrays
    real, dimension(nvarb,2,nq,nface), intent(in) :: q_face
    integer, intent(in) :: nvarb

    integer :: multirate

    !MPI Variables
    integer :: i,j,k,iv,e,ip
    real :: recv_data_dg_quad1(nvarb*nq*nboun)
    real :: send_data_dg_quad1(nvarb*nq*nboun)
    real :: q_recv_quad1(nvarb,nq,nboun), q_send_quad1(nvarb,nq,nboun)

    recv_data_dg_quad1 = 0.0
    send_data_dg_quad1 = 0.0
    q_recv_quad1 = 0.0
    q_send_quad1 = 0.0

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !Load all the boundary data into a vector
    call pack_data_dg_quad(send_data_dg_quad1,q_face,nvarb)

    !non-blocking sends-receives: message size=nmessage
    call send_bound_dg_general_quad(send_data_dg_quad1,recv_data_dg_quad1,nvarb,nreq,ireq,status)

end subroutine create_lap_precommunicator_quad

subroutine create_lap_precommunicator_quad_v1(grad_uvdp,nvarb)

    use mod_basis, only: nq

    use mod_mpi_communicator, only: ierr, ireq, nreq, status

    use mod_grid, only:  npoin, intma, nelem,nface,nboun, npoin_q

    use mod_initial, only: nvar

    use mod_metrics, only: massinv

    use mod_p4est, only: plist

    use mod_ref, only: q_send_quad, q_recv_quad, recv_data_dg_quad, send_data_dg_quad, nmessage

    implicit none

    !Global Arrays
    real, dimension(4,npoin_q), intent(in) :: grad_uvdp
    integer, intent(in) :: nvarb

    integer :: multirate

    !MPI Variables
    integer :: i,j,k,iv,e,ip
    real :: recv_data_dg_quad1(nvarb*nq*nboun)
    real :: send_data_dg_quad1(nvarb*nq*nboun)
    real :: q_recv_quad1(nvarb,nq,nboun), q_send_quad1(nvarb,nq,nboun)

    recv_data_dg_quad1 = 0.0
    send_data_dg_quad1 = 0.0
    q_recv_quad1 = 0.0
    q_send_quad1 = 0.0

    !-----------------------------------
    ! DG - Discontinuous communicator
    !-----------------------------------

    !Load all the boundary data into a vector
    call pack_data_dg_quad_lap(send_data_dg_quad1,grad_uvdp,nvarb)

    !non-blocking sends-receives: message size=nmessage
    call send_bound_dg_general_quad(send_data_dg_quad1,recv_data_dg_quad1,nvarb,nreq,ireq,status)


end subroutine create_lap_precommunicator_quad_v1

