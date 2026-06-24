
module mod_splitting

    ! ============================================================================================
    ! This module contains the routines for the barotropic-baroclinic splitting
    !   Author: Yao Gahounzo
    !   Computing PhD
    !   Boise State University
    !   Date: March 27, 2023
    ! It contains the following routines:
    ! - thickness: baroclinic substem for splitting system using two-level time integration
    ! - momentum: baroclinic substem for splitting system using two-level time integration
    ! These routines are based on Prof. Higdon 1D MLSWE code
    !
    ! ============================================================================================

    use mod_grid,             only: grid
    use mod_basis,            only: basis
    use mod_input,            only: input
    use mod_initial,          only: initial
    use mod_metrics,          only: metrics
    use mod_face,             only: face_CS
    use mod_tensor,           only: tensor_CS
    use mod_variables,        only: btp_CS, bcl_CS
    use mod_parallel,         only: parallel_CS
    use mod_ref,              only: mref
    use mod_mpi_communicator, only: mpi_communicator

    implicit none

    public :: create_rhs_bcl

contains

    subroutine create_rhs_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs, qprime_df, q_df)

        use mod_create_rhs_mlswe, only: bcl_rhs
        use mod_laplacian_quad,   only: bcl_create_laplacian

        implicit none

        type(grid),             intent(in)    :: G
        type(input),            intent(in)    :: inp
        type(basis),            intent(in)    :: b
        type(face_CS),          intent(in)    :: mf
        type(parallel_CS),      intent(in)    :: par
        type(btp_CS),           intent(inout) :: btp
        type(bcl_CS),           intent(inout) :: bcl
        type(initial),          intent(in)    :: init
        type(mref),             intent(inout) :: ref
        type(mpi_communicator), intent(inout) :: mpic
        type(tensor_CS),        intent(in)    :: tsp
        type(metrics),          intent(in)    :: mt

        real, dimension(3, G%npoin, inp%nlayers), intent(in)  :: qprime_df, q_df
        real, dimension(3, G%npoin, inp%nlayers), intent(out) :: rhs

        real, dimension(2, G%npoin, inp%nlayers) :: rhs_visc_bcl
        integer :: k

        rhs_visc_bcl = 0.0

        if (inp%method_visc > 0) call bcl_create_laplacian(G, inp, b, mf, par, btp, bcl, ref, mpic, mt, tsp, rhs_visc_bcl)

        call bcl_rhs(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs, qprime_df, q_df)

        do k = 1, inp%nlayers
            rhs(1,:,k) = mt%massinv(:)*rhs(1,:,k)
            rhs(2,:,k) = mt%massinv(:)*rhs(2,:,k) + rhs_visc_bcl(1,:,k)
            rhs(3,:,k) = mt%massinv(:)*rhs(3,:,k) + rhs_visc_bcl(2,:,k)
        end do

    end subroutine create_rhs_bcl

end module mod_splitting
