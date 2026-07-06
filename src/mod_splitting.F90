
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

    public :: thickness, momentum, momentum_mass, create_rhs_bcl, rhs_momentum

contains

    subroutine thickness(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, qprime_df, q_df, qb_df)

        ! ========================================================================================
        ! This subroutine is used to predict or correct the layer thickness for the splitting
        ! system using two-level time integration
        ! The nodal points or degree of freedom of the layer thickness is stored in qprime_df(1,:,:)
        ! The quadrature points of the layer thickness dp is stored in q_df(1,:,:)
        !
        ! Enforce consistency between the layer masses and the barotropic mass.
        ! ========================================================================================

        use mod_create_rhs_mlswe, only: layer_mass_rhs

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

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(inout) :: qprime_df
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(inout) :: q_df
        real, intent(in) :: qb_df(inp%nvar_btp, G%npoin)

        real :: dp_advec(G%npoin, inp%nlayers), ope
        integer :: k, I

        call layer_mass_rhs(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, dp_advec, qprime_df)

        do k = 1, inp%nlayers
            q_df(1,:,k) = q_df(1,:,k) + inp%dt*dp_advec(:,k)
        end do

        do k = 1, inp%nlayers
            do I = 1, G%npoin
                ope = sum(q_df(1,I,:)) / init%pbprime_df(I)
                qprime_df(1,I,k) = q_df(1,I,k) / ope
            end do
        end do

    end subroutine thickness

    subroutine momentum(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, q_df, qprime_df, qb_df)

        ! ========================================================================================
        ! This subroutine is used to correct the layer momentum for the splitting system
        ! The nodal points or degree of freedom of the layer momentum is stored in q_df(2:3,:,:)
        ! ========================================================================================

        use mod_create_rhs_mlswe, only: rhs_layer_shear_stress
        use mod_layer_terms,      only: layer_mom_boundary_df, velocity_df, extract_velocity

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

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(inout) :: q_df
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)    :: qprime_df
        real, dimension(inp%nvar_btp, G%npoin),              intent(in)    :: qb_df

        ! NOTE: this two-level splitting scheme is not yet sphere-aware — the
        ! Coriolis semi-implicit rotation below is a 2D beta-plane formula and
        ! only ever touches q_df(2:3,...) (u,v); q_df(4,...) (w) is untouched.
        real, dimension(2, G%npoin,   inp%nlayers) :: q_df_temp, rhs_mom, rhs_stress
        real, dimension(inp%nvar_bcl-1, G%npoin, inp%nlayers) :: uv_df
        real, dimension(inp%nvar_bcl, G%npoin,   inp%nlayers) :: q_df3
        real, dimension(3, G%npoin_q, inp%nlayers) :: q
        real, dimension(G%npoin) :: tempu, tempv
        integer :: k, I

        q_df_temp = 0.0

        call rhs_momentum(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs_mom, qprime_df, q_df)

        do k = 1, inp%nlayers
            q_df_temp(1,:,k) = q_df(2,:,k) + inp%dt*rhs_mom(1,:,k)
            q_df_temp(2,:,k) = q_df(3,:,k) + inp%dt*rhs_mom(2,:,k)
        end do

        if(inp%ad_mlswe > 0.0) then
            do k = 1, inp%nlayers
                do I = 1, G%npoin
                    tempu(I) = q_df_temp(1,I,k) + init%fdt2_bcl(I)*q_df(3,I,k)
                    tempv(I) = q_df_temp(2,I,k) - init%fdt2_bcl(I)*q_df(2,I,k)
                    q_df3(2,I,k) = init%a_bcl(I)*tempu(I) + init%b_bcl(I)*tempv(I)
                    q_df3(3,I,k) = -init%b_bcl(I)*tempu(I) + init%a_bcl(I)*tempv(I)
                end do
                q_df3(1,:,k) = q_df(1,:,k)
            end do

            call velocity_df(G, inp, q_df3, qb_df)
            call rhs_layer_shear_stress(G, inp, b, init, tsp, rhs_stress, q_df3)

            do k = 1, inp%nlayers
                q_df_temp(1,:,k) = q_df_temp(1,:,k) + inp%dt*(mt%massinv(:)*rhs_stress(1,:,k))
                q_df_temp(2,:,k) = q_df_temp(2,:,k) + inp%dt*(mt%massinv(:)*rhs_stress(2,:,k))
            end do
        end if

        do k = 1, inp%nlayers
            tempu(:) = q_df_temp(1,:,k) + init%fdt2_bcl(:)*q_df(3,:,k)
            tempv(:) = q_df_temp(2,:,k) - init%fdt2_bcl(:)*q_df(2,:,k)
            q_df(2,:,k) = init%a_bcl(:)*tempu(:) + init%b_bcl(:)*tempv(:)
            q_df(3,:,k) = -init%b_bcl(:)*tempu(:) + init%a_bcl(:)*tempv(:)
        end do

        call layer_mom_boundary_df(G, inp, b, mf, init, q_df)

        call extract_velocity(G, inp, uv_df, q_df, qb_df)

        do k = 1, inp%nlayers
            q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
            q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
        end do

    end subroutine momentum

    subroutine momentum_mass(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, q_df, qprime_df, qb_df)

        ! ========================================================================================
        ! This subroutine is used to predict the layer mass and momentum for the splitting system
        ! The nodal points or degree of freedom of the layer mass is stored in q_df(1,:,:),
        ! momentum is stored in q_df(2:3,:,:)
        ! ========================================================================================

        use mod_create_rhs_mlswe, only: rhs_layer_shear_stress
        use mod_layer_terms,      only: layer_mom_boundary_df, velocity_df, &
                                        extract_qprime_df_face, extract_velocity

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

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(inout) :: q_df
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)    :: qprime_df
        real, dimension(inp%nvar_btp, G%npoin),              intent(in)    :: qb_df

        ! NOTE: this two-level splitting scheme is not yet sphere-aware — the
        ! Coriolis semi-implicit rotation below is a 2D beta-plane formula and
        ! only ever touches q_df(2:3,...) (u,v); q_df(4,...) (w) is untouched.
        real, dimension(2, G%npoin,   inp%nlayers) :: rhs_stress
        real, dimension(inp%nvar_bcl-1, G%npoin, inp%nlayers) :: uv_df
        real, dimension(inp%nvar_bcl, G%npoin,   inp%nlayers) :: q_df3, q_df_temp
        real, dimension(3, G%npoin_q, inp%nlayers) :: q
        real, dimension(G%npoin) :: tempu, tempv
        integer :: k, I

        call create_rhs_bcl(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, bcl%rhs_bcl, qprime_df, q_df)

        do k = 1, inp%nlayers
            q_df_temp(1,:,k) = q_df(1,:,k) + inp%dt*bcl%rhs_bcl(1,:,k)
            q_df_temp(2,:,k) = q_df(2,:,k) + inp%dt*bcl%rhs_bcl(2,:,k)
            q_df_temp(3,:,k) = q_df(3,:,k) + inp%dt*bcl%rhs_bcl(3,:,k)

            if(any(q_df(1,:,k) < 0.0)) then
                write(*,*) 'Negative mass in thickness at some points'
                stop
            end if
        end do

        if(inp%ad_mlswe > 0.0) then
            do k = 1, inp%nlayers
                do I = 1, G%npoin
                    tempu(I) = q_df_temp(2,I,k) + init%fdt2_bcl(I)*q_df(3,I,k)
                    tempv(I) = q_df_temp(3,I,k) - init%fdt2_bcl(I)*q_df(2,I,k)
                    q_df3(2,I,k) = init%a_bcl(I)*tempu(I) + init%b_bcl(I)*tempv(I)
                    q_df3(3,I,k) = -init%b_bcl(I)*tempu(I) + init%a_bcl(I)*tempv(I)
                end do
                q_df3(1,:,k) = q_df(1,:,k)
            end do

            call velocity_df(G, inp, q_df3, qb_df)
            call rhs_layer_shear_stress(G, inp, b, init, tsp, rhs_stress, q_df3)

            do k = 1, inp%nlayers
                q_df_temp(2,:,k) = q_df_temp(2,:,k) + inp%dt*(mt%massinv(:)*rhs_stress(1,:,k))
                q_df_temp(3,:,k) = q_df_temp(3,:,k) + inp%dt*(mt%massinv(:)*rhs_stress(2,:,k))
            end do
        end if

        do k = 1, inp%nlayers
            tempu(:) = q_df_temp(2,:,k) + init%fdt2_bcl(:)*q_df(3,:,k)
            tempv(:) = q_df_temp(3,:,k) - init%fdt2_bcl(:)*q_df(2,:,k)
            q_df(2,:,k) = init%a_bcl(:)*tempu(:) + init%b_bcl(:)*tempv(:)
            q_df(3,:,k) = -init%b_bcl(:)*tempu(:) + init%a_bcl(:)*tempv(:)
            q_df(1,:,k) = q_df_temp(1,:,k)
        end do

        call layer_mom_boundary_df(G, inp, b, mf, init, q_df)

        call extract_velocity(G, inp, uv_df, q_df, qb_df)

        do k = 1, inp%nlayers
            q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
            q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
        end do

    end subroutine momentum_mass

    subroutine rhs_momentum(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs_mom, qprime_df, q_df)

        use mod_create_rhs_mlswe, only: layer_momentum_rhs
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

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: qprime_df, q_df
        real, dimension(2, G%npoin, inp%nlayers), intent(out) :: rhs_mom

        integer :: k

        bcl%rhs_visc_bcl = 0.0

        if (inp%method_visc > 0) call bcl_create_laplacian(G, inp, b, mf, par, btp, bcl, ref, mpic, mt, tsp, bcl%rhs_visc_bcl)
        ! bcl_create_laplacian downloads rhs_visc_bcl to host internally; host copy is valid here.

        call layer_momentum_rhs(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs_mom, qprime_df, q_df)

        do k = 1, inp%nlayers
            rhs_mom(1,:,k) = mt%massinv(:)*rhs_mom(1,:,k) + bcl%rhs_visc_bcl(1,:,k)
            rhs_mom(2,:,k) = mt%massinv(:)*rhs_mom(2,:,k) + bcl%rhs_visc_bcl(2,:,k)
        end do

    end subroutine rhs_momentum

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

        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(in)  :: qprime_df, q_df
        real, dimension(inp%nvar_bcl, G%npoin, inp%nlayers), intent(out) :: rhs

        integer :: k, I, iv, nlayers_l, npoin_l, nvarb_l
        real    :: visc_term

        nlayers_l = inp%nlayers
        npoin_l   = G%npoin
        nvarb_l   = inp%nvar_bcl

        !$acc kernels present(bcl%rhs_visc_bcl)
        bcl%rhs_visc_bcl = 0.0
        !$acc end kernels

        if (inp%method_visc > 0) call bcl_create_laplacian(G, inp, b, mf, par, btp, bcl, ref, mpic, mt, tsp, bcl%rhs_visc_bcl)

        call bcl_rhs(G, inp, b, mf, par, btp, bcl, init, ref, mpic, tsp, mt, rhs, qprime_df, q_df)
        ! rhs stays on device; apply mass scaling and viscous term on GPU.

        ! Momentum rows (u,v[,w]) get mass-matrix scaling plus viscosity.
        ! bcl%rhs_visc_bcl carries one row per momentum component (u,v[,w]),
        ! i.e. rows 1..nvar_bcl-1, matching rhs rows 2..nvar_bcl.
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
