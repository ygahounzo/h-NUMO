! ================================================================================================
! This module contains the routines for the barotropic flux terms
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: March 27, 2023
! ================================================================================================

module mod_barotropic_terms

    implicit none

    public :: &
                btp_mom_boundary_df, &
                btp_bcl_coeffs_qdf

    contains

    subroutine btp_mom_boundary_df(G, b, mf, qb)

        use mod_basis,    only: basis
        use mod_grid,     only: grid
        use mod_face,     only: face_CS

        implicit none

        type(grid),    intent(in)    :: G
        type(basis),   intent(in)    :: b
        type(face_CS), intent(in)    :: mf

        real, intent(inout) :: qb(4,G%npoin)

        integer :: iface, il, jl, el, er, I, kl, n
        real :: nx, ny, unl

        !$acc parallel loop present(G%face, G%intma, mf%imapl, mf%normal_vector, qb) &
        !$acc    private(il, jl, kl, el, er, I, nx, ny, unl)
        do iface = 1, G%nface

            el = G%face(7,iface)
            er = G%face(8,iface)

            if (er == -4) then
                !$acc loop seq
                do n = 1, b%ngl
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    nx = mf%normal_vector(1,n,1,iface)
                    ny = mf%normal_vector(2,n,1,iface)
                    unl = qb(3,I)*nx + qb(4,I)*ny
                    !$acc atomic update
                    qb(3,I) = qb(3,I) - unl*nx
                    !$acc atomic update
                    qb(4,I) = qb(4,I) - unl*ny
                end do

            elseif (er == -2) then
                !$acc loop seq
                do n = 1, b%ngl
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    !$acc atomic write
                    qb(3,I) = 0.0
                    !$acc atomic write
                    qb(4,I) = 0.0
                end do
            end if
        end do
        !$acc end parallel loop

    end subroutine btp_mom_boundary_df

    subroutine btp_bcl_coeffs_qdf(G, inp, b, tsp, bcl, btp, qprime_df)

        ! Compute baroclinic coefficients in the advective barotropic momentum fluxes
        ! and in the barotropic pressure forcing.

        use mod_grid,      only: grid
        use mod_input,     only: input
        use mod_basis,     only: basis
        use mod_tensor,    only: tensor_CS
        use mod_variables, only: bcl_CS, btp_CS

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(tensor_CS), intent(in)    :: tsp
        type(bcl_CS),    intent(inout) :: bcl
        type(btp_CS),    intent(inout) :: btp

        real, dimension(3,G%npoin,inp%nlayers), intent(in) :: qprime_df

        integer :: Iq, ip, k, npoin_l, nlayers_l, npts_l, I
        real :: dhdx, dhdy, dpp
        real :: du_dx, du_dy, dv_dx, dv_dy
        real :: btg1, btg2, btg3, btg4, pbps

        npoin_l   = G%npoin
        nlayers_l = inp%nlayers
        npts_l    = b%npts

        ! Fused single kernel: gang over nodes, seq over layers and neighbours.
        ! Eliminates the graduv temporary and the sequential k-loop with its
        ! 3*nlayers kernel launches.  btp sums accumulated privately per node.
        !$acc parallel loop gang &
        !$acc    present(bcl%dpprime_visc, bcl%dpp_graduv, bcl%dpp_uvp, &
        !$acc            btp%btp_dpp_graduv, btp%pbprime_visc, qprime_df, &
        !$acc            tsp%index_df, tsp%dpsidx_df, tsp%dpsidy_df) &
        !$acc    firstprivate(npoin_l, nlayers_l, npts_l) &
        !$acc    private(dhdx, dhdy, dpp, du_dx, du_dy, dv_dx, dv_dy, &
        !$acc            btg1, btg2, btg3, btg4, pbps, I)
        do Iq = 1, npoin_l
            btg1 = 0.0; btg2 = 0.0; btg3 = 0.0; btg4 = 0.0; pbps = 0.0
            !$acc loop seq
            do k = 1, nlayers_l
                dpp = qprime_df(1,Iq,k)
                bcl%dpprime_visc(Iq,k) = dpp

                du_dx = 0.0; du_dy = 0.0; dv_dx = 0.0; dv_dy = 0.0
                !$acc loop seq
                do ip = 1, npts_l
                    I    = tsp%index_df(ip,Iq)
                    dhdx = tsp%dpsidx_df(ip,Iq)
                    dhdy = tsp%dpsidy_df(ip,Iq)
                    du_dx = du_dx + dhdx*qprime_df(2,I,k)
                    du_dy = du_dy + dhdy*qprime_df(2,I,k)
                    dv_dx = dv_dx + dhdx*qprime_df(3,I,k)
                    dv_dy = dv_dy + dhdy*qprime_df(3,I,k)
                end do

                bcl%dpp_graduv(1,Iq,k) = dpp*du_dx
                bcl%dpp_graduv(2,Iq,k) = dpp*du_dy
                bcl%dpp_graduv(3,Iq,k) = dpp*dv_dx
                bcl%dpp_graduv(4,Iq,k) = dpp*dv_dy

                bcl%dpp_uvp(1,Iq,k) = dpp*qprime_df(2,Iq,k)
                bcl%dpp_uvp(2,Iq,k) = dpp*qprime_df(3,Iq,k)

                btg1 = btg1 + bcl%dpp_graduv(1,Iq,k)
                btg2 = btg2 + bcl%dpp_graduv(2,Iq,k)
                btg3 = btg3 + bcl%dpp_graduv(3,Iq,k)
                btg4 = btg4 + bcl%dpp_graduv(4,Iq,k)
                pbps = pbps + dpp
            end do
            btp%btp_dpp_graduv(1,Iq) = btg1
            btp%btp_dpp_graduv(2,Iq) = btg2
            btp%btp_dpp_graduv(3,Iq) = btg3
            btp%btp_dpp_graduv(4,Iq) = btg4
            btp%pbprime_visc(Iq)     = pbps
        end do
        !$acc end parallel loop

    end subroutine btp_bcl_coeffs_qdf

end module mod_barotropic_terms
