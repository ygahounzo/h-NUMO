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

        !$acc parallel loop present(G, b, mf, qb)
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
                    unl     = qb(3,I)*nx + qb(4,I)*ny
                    qb(3,I) = qb(3,I) - unl*nx
                    qb(4,I) = qb(4,I) - unl*ny
                end do

            elseif (er == -2) then
                !$acc loop seq
                do n = 1, b%ngl
                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)
                    qb(3,I) = 0.0
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

        real, dimension(4,G%npoin) :: graduv
        integer :: k, I, Iq, ip, npoin_l, nlayers_l, npts_l
        real :: dhdx, dhdy

        npoin_l   = G%npoin
        nlayers_l = inp%nlayers
        npts_l    = b%npts

        !$acc data create(graduv) &
        !$acc      present(bcl%dpprime_visc, bcl%dpp_graduv, bcl%dpp_uvp, &
        !$acc              btp%btp_dpp_graduv, btp%pbprime_visc, qprime_df)

        !$acc kernels present(btp%btp_dpp_graduv, btp%pbprime_visc, bcl%dpprime_visc, qprime_df)
        btp%btp_dpp_graduv    = 0.0
        btp%pbprime_visc      = 0.0
        bcl%dpprime_visc(:,:) = qprime_df(1,:,:)
        !$acc end kernels

        ! Sequential over k: each GPU kernel completes before the next k starts.
        ! Accumulation into btp% arrays is correct — no concurrent k-iterations.
        do k = 1, nlayers_l

            ! Gather velocity gradient directly from qprime_df — no races on graduv.
            !$acc kernels present(graduv)
            graduv = 0.0
            !$acc end kernels

            !$acc parallel loop gang &
            !$acc    present(graduv, qprime_df, tsp%index_df, tsp%dpsidx_df, tsp%dpsidy_df) &
            !$acc    firstprivate(npoin_l, npts_l, k) private(I, dhdx, dhdy)
            do Iq = 1, npoin_l
                !$acc loop seq
                do ip = 1, npts_l
                    I    = tsp%index_df(ip,Iq)
                    dhdx = tsp%dpsidx_df(ip,Iq)
                    dhdy = tsp%dpsidy_df(ip,Iq)
                    graduv(1,Iq) = graduv(1,Iq) + dhdx*qprime_df(2,I,k)
                    graduv(2,Iq) = graduv(2,Iq) + dhdy*qprime_df(2,I,k)
                    graduv(3,Iq) = graduv(3,Iq) + dhdx*qprime_df(3,I,k)
                    graduv(4,Iq) = graduv(4,Iq) + dhdy*qprime_df(3,I,k)
                end do
            end do
            !$acc end parallel loop

            !$acc parallel loop gang &
            !$acc    present(bcl%dpp_graduv, bcl%dpprime_visc, bcl%dpp_uvp, &
            !$acc            btp%btp_dpp_graduv, btp%pbprime_visc, graduv, qprime_df) &
            !$acc    firstprivate(npoin_l, k)
            do I = 1, npoin_l
                bcl%dpp_graduv(1,I,k) = bcl%dpprime_visc(I,k)*graduv(1,I)
                bcl%dpp_graduv(2,I,k) = bcl%dpprime_visc(I,k)*graduv(2,I)
                bcl%dpp_graduv(3,I,k) = bcl%dpprime_visc(I,k)*graduv(3,I)
                bcl%dpp_graduv(4,I,k) = bcl%dpprime_visc(I,k)*graduv(4,I)

                bcl%dpp_uvp(1,I,k) = bcl%dpprime_visc(I,k)*qprime_df(2,I,k)
                bcl%dpp_uvp(2,I,k) = bcl%dpprime_visc(I,k)*qprime_df(3,I,k)

                btp%btp_dpp_graduv(1,I) = btp%btp_dpp_graduv(1,I) + bcl%dpp_graduv(1,I,k)
                btp%btp_dpp_graduv(2,I) = btp%btp_dpp_graduv(2,I) + bcl%dpp_graduv(2,I,k)
                btp%btp_dpp_graduv(3,I) = btp%btp_dpp_graduv(3,I) + bcl%dpp_graduv(3,I,k)
                btp%btp_dpp_graduv(4,I) = btp%btp_dpp_graduv(4,I) + bcl%dpp_graduv(4,I,k)
                btp%pbprime_visc(I)     = btp%pbprime_visc(I)     + bcl%dpprime_visc(I,k)
            end do
            !$acc end parallel loop

        end do

        !$acc end data

    end subroutine btp_bcl_coeffs_qdf

end module mod_barotropic_terms
