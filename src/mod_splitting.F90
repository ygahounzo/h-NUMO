
module mod_splitting

    ! ============================================================================================
    ! This module contains the routines for the barotropic-baroclinic splitting
    !   Author: Yao Gahounzo
    !   Computing PhD
    !   Boise State University
    !   Date: March 27, 2023
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
        use mod_laplacian_quad,   only: bcl_create_laplacian, compute_bcl_lap_z
        use mod_mpi_utilities,    only: irank, MPI_PRECISION
        use mod_constants,        only: gravity
        use mpi

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

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: qprime_df, q_df
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(out) :: rhs

        integer :: k, I, iv, nlayers_l, npoin_l, nvarb_l
        integer :: ierr_diag
        real    :: visc_term, lap_z_loc_max, lap_z_glb_max

        nlayers_l = inp%nlayers
        npoin_l   = G%npoin
        nvarb_l   = inp%nvar_bcl

        !$acc kernels present(bcl%rhs_visc_bcl)
        bcl%rhs_visc_bcl = 0.0
        !$acc end kernels

        ! ∇²z_k needed for APE (independent of nu_smag).
        if (inp%c_APE > 0.0) then
            !$acc update host(btp%ope2_ave_df, qprime_df)
            call compute_bcl_lap_z(G, inp, b, mf, par, ref, mpic, init, btp, bcl, mt, tsp, qprime_df, bcl%lap_z_df)
            ! lap_z_loc_max = maxval(abs(bcl%lap_z_df(:, 2:inp%nlayers)))
            ! call MPI_Reduce(lap_z_loc_max, lap_z_glb_max, 1, MPI_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr_diag)
            ! if (irank == 0) write(*,'("lap_z global max =",es12.4)') lap_z_glb_max
            !$acc update device(bcl%lap_z_df)
        end if

        call bcl_rhs(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs, qprime_df, q_df)
        ! rhs stays on device; apply mass scaling and viscous term on GPU.

        ! Momentum rows (u,v[,w]) get mass-matrix scaling plus viscosity.
        ! bcl%rhs_visc_bcl carries one row per momentum component (u,v[,w]),
        ! i.e. rows 1..nvar_bcl-1, matching rhs rows 2..nvar_bcl.
        !
        ! nu_smag is a strain-rate quantity recomputed fresh inside
        ! bcl_create_laplacian every call; rhs_visc_bcl = (visc+nu_smag)*massinv*L(q).
        if (inp%method_visc > 0) call bcl_create_laplacian(G, inp, b, mf, par, btp, bcl, ref, mpic, mt, tsp, bcl%rhs_visc_bcl)

        !$acc parallel loop gang collapse(2) &
        !$acc    present(rhs, bcl%rhs_visc_bcl, mt%massinv) private(iv, visc_term) &
        !$acc    firstprivate(nlayers_l, npoin_l, nvarb_l)
        do k = 1, nlayers_l
            do I = 1, npoin_l
                rhs(1,I,k) = mt%massinv(I)*rhs(1,I,k)
                !$acc loop seq
                do iv = 2, nvarb_l
                    visc_term = bcl%rhs_visc_bcl(iv-1,I,k)
                    rhs(iv,I,k) = mt%massinv(I)*rhs(iv,I,k) + visc_term
                end do
            end do
        end do
        !$acc end parallel loop

        !$acc update host(rhs)

    end subroutine create_rhs_bcl

end module mod_splitting
