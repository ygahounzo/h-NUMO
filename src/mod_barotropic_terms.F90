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
                compute_gradient_uv, &
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

        integer :: iface, ilr, m, il, jl, el, er, ilocl, ilocr, I, kl, n
        real :: nx, ny, unl, upnl

        do iface = 1, G%nface

            el = G%face(7,iface)
            er = G%face(8,iface)

            if (er == -4) then
                do n = 1, b%ngl

                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)

                    nx = mf%normal_vector(1,n,1,iface)
                    ny = mf%normal_vector(2,n,1,iface)

                    unl    = qb(3,I)*nx + qb(4,I)*ny
                    qb(3,I) = qb(3,I) - unl*nx
                    qb(4,I) = qb(4,I) - unl*ny
                end do

            elseif (er == -2) then
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
        real, dimension(2,G%npoin) :: uv
        integer :: k, I

        btp%btp_dpp_graduv = 0.0
        btp%pbprime_visc   = 0.0

        bcl%dpprime_visc(:,:) = qprime_df(1,:,:)

        do k = 1, inp%nlayers

            uv = qprime_df(2:3,:,k)
            call compute_gradient_uv(G, b, tsp, graduv, uv)

            do I = 1, G%npoin

                bcl%dpp_graduv(1,I,k) = bcl%dpprime_visc(I,k)*graduv(1,I)
                bcl%dpp_graduv(2,I,k) = bcl%dpprime_visc(I,k)*graduv(2,I)
                bcl%dpp_graduv(3,I,k) = bcl%dpprime_visc(I,k)*graduv(3,I)
                bcl%dpp_graduv(4,I,k) = bcl%dpprime_visc(I,k)*graduv(4,I)

                bcl%dpp_uvp(1,I,k) = bcl%dpprime_visc(I,k)*qprime_df(2,I,k)
                bcl%dpp_uvp(2,I,k) = bcl%dpprime_visc(I,k)*qprime_df(3,I,k)

                btp%btp_dpp_graduv(:,I) = btp%btp_dpp_graduv(:,I) + bcl%dpp_graduv(:,I,k)
                btp%pbprime_visc(I)     = btp%pbprime_visc(I) + bcl%dpprime_visc(I,k)
            end do
        end do

    end subroutine btp_bcl_coeffs_qdf

    subroutine compute_gradient_uv(G, b, tsp, grad_uv, uv)

        use mod_grid,   only: grid
        use mod_basis,  only: basis
        use mod_tensor, only: tensor_CS

        implicit none

        type(grid),      intent(in)  :: G
        type(basis),     intent(in)  :: b
        type(tensor_CS), intent(in)  :: tsp

        real, dimension(2,G%npoin), intent(in)  :: uv
        real, dimension(4,G%npoin), intent(out) :: grad_uv

        integer :: Iq, I, ip
        real :: dhdx, dhdy

        grad_uv = 0.0

        do concurrent(Iq = 1:G%npoin)

            do ip = 1, b%npts

                I    = tsp%index_df(ip,Iq)
                dhdx = tsp%dpsidx_df(ip,Iq)
                dhdy = tsp%dpsidy_df(ip,Iq)

                grad_uv(1,Iq) = grad_uv(1,Iq) + dhdx*uv(1,I)
                grad_uv(2,Iq) = grad_uv(2,Iq) + dhdy*uv(1,I)
                grad_uv(3,Iq) = grad_uv(3,Iq) + dhdx*uv(2,I)
                grad_uv(4,Iq) = grad_uv(4,Iq) + dhdy*uv(2,I)

            end do
        end do

    end subroutine compute_gradient_uv

end module mod_barotropic_terms
