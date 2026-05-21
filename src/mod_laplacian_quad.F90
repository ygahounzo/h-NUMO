! =================================================================================================
! This module contains the routines for the barotropic and baroclinic horizontal viscosity terms
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: July 24, 2023
! =================================================================================================

module mod_laplacian_quad

    implicit none

    public :: bcl_create_laplacian, btp_create_laplacian

    contains

    subroutine btp_create_laplacian(G, inp, b, mf, par, btp, init, ref, mpic, tsp, rhs_lap, qb_df)

        use mod_grid,             only: grid
        use mod_input,            only: input
        use mod_basis,            only: basis
        use mod_face,             only: face_CS
        use mod_parallel,         only: parallel_CS
        use mod_variables,        only: btp_CS
        use mod_initial,          only: initial
        use mod_ref,              only: mref
        use mod_mpi_communicator, only: mpi_communicator
        use mod_tensor,           only: tensor_CS
        use mod_barotropic_terms, only: compute_gradient_uv

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
        type(tensor_CS),        intent(in)    :: tsp

        real, intent(out) :: rhs_lap(2,G%npoin)
        real, dimension(4,G%npoin), intent(in) :: qb_df

        real, dimension(2,G%npoin) :: Uk
        real, dimension(4,G%npoin) :: graduv

        Uk(1,:) = qb_df(3,:)/qb_df(1,:)
        Uk(2,:) = qb_df(4,:)/qb_df(1,:)

        call compute_gradient_uv(G, b, tsp, graduv, Uk)

        btp%graduvb_ave = btp%graduvb_ave + graduv

        call btp_lap_create_precommunicator(G, b, mf, init, par, btp, ref, mpic, graduv, 4)

        call btp_compute_laplacian(G, b, btp, tsp, rhs_lap, graduv)

        call create_rhs_laplacian_flux(G, b, mf, btp, rhs_lap, graduv)

        call create_rhs_lap_postcommunicator_df(G, b, mf, par, btp, ref, mpic, rhs_lap, 4)

        rhs_lap(1,:) = inp%visc_mlswe*rhs_lap(1,:)
        rhs_lap(2,:) = inp%visc_mlswe*rhs_lap(2,:)

    end subroutine btp_create_laplacian

    subroutine bcl_create_laplacian(G, inp, b, mf, par, btp, bcl, ref, mpic, mt, tsp, rhs_lap)

        use mod_grid,             only: grid
        use mod_input,            only: input
        use mod_basis,            only: basis
        use mod_face,             only: face_CS
        use mod_parallel,         only: parallel_CS
        use mod_variables,        only: btp_CS, bcl_CS
        use mod_ref,              only: mref
        use mod_mpi_communicator, only: mpi_communicator
        use mod_metrics,          only: metrics
        use mod_tensor,           only: tensor_CS

        implicit none

        type(grid),             intent(in)    :: G
        type(input),            intent(in)    :: inp
        type(basis),            intent(in)    :: b
        type(face_CS),          intent(in)    :: mf
        type(parallel_CS),      intent(in)    :: par
        type(btp_CS),           intent(in)    :: btp
        type(bcl_CS),           intent(in)    :: bcl
        type(mref),             intent(inout) :: ref
        type(mpi_communicator), intent(inout) :: mpic
        type(metrics),          intent(in)    :: mt
        type(tensor_CS),        intent(in)    :: tsp

        real, intent(out) :: rhs_lap(2,G%npoin,inp%nlayers)

        integer :: k

        call bcl_lap_create_precommunicator(G, inp, b, mf, par, ref, mpic, bcl%dpp_graduv, bcl%dpprime_visc)

        call bcl_compute_laplacian(G, inp, b, btp, bcl, tsp, rhs_lap)
        call bcl_create_rhs_laplacian_flux(G, inp, b, mf, btp, bcl, rhs_lap)

        call bcl_create_rhs_lap_postcommunicator_df(G, inp, b, mf, par, btp, ref, mpic, rhs_lap)

        do k = 1, inp%nlayers
            rhs_lap(1,:,k) = inp%visc_mlswe*mt%massinv(:)*rhs_lap(1,:,k)
            rhs_lap(2,:,k) = inp%visc_mlswe*mt%massinv(:)*rhs_lap(2,:,k)
        end do

    end subroutine bcl_create_laplacian

    subroutine btp_compute_laplacian(G, b, btp, tsp, lap_q, grad_dpuvp)

        use mod_grid,      only: grid
        use mod_basis,     only: basis
        use mod_variables, only: btp_CS
        use mod_tensor,    only: tensor_CS

        implicit none

        type(grid),      intent(in)  :: G
        type(basis),     intent(in)  :: b
        type(btp_CS),    intent(in)  :: btp
        type(tensor_CS), intent(in)  :: tsp

        real, intent(out) :: lap_q(2,G%npoin)
        real, dimension(4,G%npoin), intent(in) :: grad_dpuvp

        integer :: Iq, I, ip
        real :: wq, qq(4)

        lap_q = 0.0

        do concurrent(Iq = 1:G%npoin)

            wq = tsp%wjac_df(Iq)

            qq(1) = btp%pbprime_visc(Iq)*grad_dpuvp(1,Iq) + btp%btp_dpp_graduv(1,Iq)
            qq(2) = btp%pbprime_visc(Iq)*grad_dpuvp(2,Iq) + btp%btp_dpp_graduv(2,Iq)
            qq(3) = btp%pbprime_visc(Iq)*grad_dpuvp(3,Iq) + btp%btp_dpp_graduv(3,Iq)
            qq(4) = btp%pbprime_visc(Iq)*grad_dpuvp(4,Iq) + btp%btp_dpp_graduv(4,Iq)

            do ip = 1, b%npts
                I = tsp%index_df(ip,Iq)
                lap_q(1,I) = lap_q(1,I) - wq*(tsp%dpsidx_df(ip,Iq)*qq(1) + tsp%dpsidy_df(ip,Iq)*qq(2))
                lap_q(2,I) = lap_q(2,I) - wq*(tsp%dpsidx_df(ip,Iq)*qq(3) + tsp%dpsidy_df(ip,Iq)*qq(4))
            end do
        end do

    end subroutine btp_compute_laplacian

    subroutine bcl_compute_laplacian(G, inp, b, btp, bcl, tsp, lap_q)

        use mod_grid,      only: grid
        use mod_basis,     only: basis
        use mod_input,     only: input
        use mod_variables, only: btp_CS, bcl_CS
        use mod_tensor,    only: tensor_CS

        implicit none

        type(grid),      intent(in)  :: G
        type(basis),     intent(in)  :: b
        type(input),     intent(in)  :: inp
        type(btp_CS),    intent(in)  :: btp
        type(bcl_CS),    intent(in)  :: bcl
        type(tensor_CS), intent(in)  :: tsp

        real, intent(out) :: lap_q(2,G%npoin,inp%nlayers)

        integer :: Iq, I, ip, k
        real :: wq, qq(4)

        lap_q = 0.0

        do concurrent(k = 1:inp%nlayers, Iq = 1:G%npoin)

            wq = tsp%wjac_df(Iq)

            qq(1) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(1,Iq) + bcl%dpp_graduv(1,Iq,k)
            qq(2) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(2,Iq) + bcl%dpp_graduv(2,Iq,k)
            qq(3) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(3,Iq) + bcl%dpp_graduv(3,Iq,k)
            qq(4) = bcl%dpprime_visc(Iq,k)*btp%graduvb_ave(4,Iq) + bcl%dpp_graduv(4,Iq,k)

            do ip = 1, b%npts
                I = tsp%index_df(ip,Iq)
                lap_q(1,I,k) = lap_q(1,I,k) - wq*(tsp%dpsidx_df(ip,Iq)*qq(1) + tsp%dpsidy_df(ip,Iq)*qq(2))
                lap_q(2,I,k) = lap_q(2,I,k) - wq*(tsp%dpsidx_df(ip,Iq)*qq(3) + tsp%dpsidy_df(ip,Iq)*qq(4))
            end do
        end do

    end subroutine bcl_compute_laplacian

    subroutine create_rhs_laplacian_flux(G, b, mf, btp, rhs, gradq)

        use mod_grid,      only: grid
        use mod_basis,     only: basis
        use mod_face,      only: face_CS
        use mod_variables, only: btp_CS

        implicit none

        type(grid),    intent(in)    :: G
        type(basis),   intent(in)    :: b
        type(face_CS), intent(in)    :: mf
        type(btp_CS),  intent(inout) :: btp

        real, intent(inout) :: rhs(2,G%npoin)
        real, intent(in)    :: gradq(4,G%npoin)

        real, dimension(2) :: qu_mean, qv_mean
        real, dimension(4,2) :: flux_uv_visc_face
        real, dimension(2) :: qul, qur
        real, dimension(2) :: qvl, qvr
        real :: nx, ny, nz
        real :: wq, un
        integer :: iface, i, j, k, il, jl, kl, ir, jr, kr, el, er
        integer :: iel, ier, ilocl, ilocr, ip
        integer :: iquad, jquad, ivar
        real :: flux_qu, flux_qv, hi, mul, mur, c_jump, alpha, beta, iflux
        real, dimension(4,b%ngl) :: ql, qr
        real, dimension(5,b%ngl) :: btp_ql, btp_qr

        beta  = 0.5
        alpha = 1.0 - beta
        iflux = 0.0

        do concurrent(iface=1:G%nface, iquad=1:b%ngl)

            if (G%face_type(iface) == 2) cycle

            iel = G%face(7,iface)
            ier = G%face(8,iface)

            nx = mf%normal_vector(1,iquad,1,iface)
            ny = mf%normal_vector(2,iquad,1,iface)

            il = mf%imapl(1,iquad,1,iface)
            jl = mf%imapl(2,iquad,1,iface)
            kl = mf%imapl(3,iquad,1,iface)
            ip = G%intma(il,jl,kl,iel)

            ql(1,iquad) = gradq(1,ip)
            ql(2,iquad) = gradq(2,ip)
            ql(3,iquad) = gradq(3,ip)
            ql(4,iquad) = gradq(4,ip)

            btp_ql(1,iquad) = btp%btp_dpp_graduv(1,ip)
            btp_ql(2,iquad) = btp%btp_dpp_graduv(2,ip)
            btp_ql(3,iquad) = btp%btp_dpp_graduv(3,ip)
            btp_ql(4,iquad) = btp%btp_dpp_graduv(4,ip)
            btp_ql(5,iquad) = btp%pbprime_visc(ip)

            if (ier > 0) then

                ir = mf%imapr(1,iquad,1,iface)
                jr = mf%imapr(2,iquad,1,iface)
                kr = mf%imapr(3,iquad,1,iface)
                ip = G%intma(ir,jr,kr,ier)

                qr(1,iquad) = gradq(1,ip)
                qr(2,iquad) = gradq(2,ip)
                qr(3,iquad) = gradq(3,ip)
                qr(4,iquad) = gradq(4,ip)

                btp_qr(1,iquad) = btp%btp_dpp_graduv(1,ip)
                btp_qr(2,iquad) = btp%btp_dpp_graduv(2,ip)
                btp_qr(3,iquad) = btp%btp_dpp_graduv(3,ip)
                btp_qr(4,iquad) = btp%btp_dpp_graduv(4,ip)
                btp_qr(5,iquad) = btp%pbprime_visc(ip)

            else
                qr(1,iquad) = ql(1,iquad)
                qr(2,iquad) = ql(2,iquad)
                qr(3,iquad) = ql(3,iquad)
                qr(4,iquad) = ql(4,iquad)

                if (ier == -4) then

                    un = ql(1,iquad)*nx + ql(2,iquad)*ny
                    qr(1,iquad) = ql(1,iquad) - 2.0*un*nx
                    qr(2,iquad) = ql(2,iquad) - 2.0*un*ny

                    un = ql(3,iquad)*nx + ql(4,iquad)*ny
                    qr(3,iquad) = ql(3,iquad) - 2.0*un*nx
                    qr(4,iquad) = ql(4,iquad) - 2.0*un*ny

                    un = btp_ql(1,iquad)*nx + btp_ql(2,iquad)*ny
                    btp_qr(1,iquad) = btp_ql(1,iquad) - 2.0*un*nx
                    btp_qr(2,iquad) = btp_ql(2,iquad) - 2.0*un*ny

                    un = btp_ql(3,iquad)*nx + btp_ql(4,iquad)*ny
                    btp_qr(3,iquad) = btp_ql(3,iquad) - 2.0*un*nx
                    btp_qr(4,iquad) = btp_ql(4,iquad) - 2.0*un*ny
                end if
            end if

            btp%graduvb_face_ave(1,1,iquad,iface) = btp%graduvb_face_ave(1,1,iquad,iface) + ql(1,iquad)
            btp%graduvb_face_ave(2,1,iquad,iface) = btp%graduvb_face_ave(2,1,iquad,iface) + ql(2,iquad)
            btp%graduvb_face_ave(3,1,iquad,iface) = btp%graduvb_face_ave(3,1,iquad,iface) + ql(3,iquad)
            btp%graduvb_face_ave(4,1,iquad,iface) = btp%graduvb_face_ave(4,1,iquad,iface) + ql(4,iquad)

            btp%graduvb_face_ave(1,2,iquad,iface) = btp%graduvb_face_ave(1,2,iquad,iface) + qr(1,iquad)
            btp%graduvb_face_ave(2,2,iquad,iface) = btp%graduvb_face_ave(2,2,iquad,iface) + qr(2,iquad)
            btp%graduvb_face_ave(3,2,iquad,iface) = btp%graduvb_face_ave(3,2,iquad,iface) + qr(3,iquad)
            btp%graduvb_face_ave(4,2,iquad,iface) = btp%graduvb_face_ave(4,2,iquad,iface) + qr(4,iquad)

            do ivar = 1, 4
                flux_uv_visc_face(ivar,1) = btp_ql(5,iquad)*ql(ivar,iquad) + btp_ql(ivar,iquad)
                flux_uv_visc_face(ivar,2) = btp_qr(5,iquad)*qr(ivar,iquad) + btp_qr(ivar,iquad)
            end do

            qul(1) = flux_uv_visc_face(1,1)
            qul(2) = flux_uv_visc_face(2,1)
            qvl(1) = flux_uv_visc_face(3,1)
            qvl(2) = flux_uv_visc_face(4,1)

            qur(1) = flux_uv_visc_face(1,2)
            qur(2) = flux_uv_visc_face(2,2)
            qvr(1) = flux_uv_visc_face(3,2)
            qvr(2) = flux_uv_visc_face(4,2)

            qu_mean(1) = alpha*qul(1) + beta*qur(1)
            qu_mean(2) = alpha*qul(2) + beta*qur(2)
            qv_mean(1) = alpha*qvl(1) + beta*qvr(1)
            qv_mean(2) = alpha*qvl(2) + beta*qvr(2)

            wq = mf%jac_face(iquad,1,iface)

            flux_qu = (qu_mean(1) - iflux*qul(1))*nx + (qu_mean(2) - iflux*qul(2))*ny
            flux_qv = (qv_mean(1) - iflux*qvl(1))*nx + (qv_mean(2) - iflux*qvl(2))*ny

            do i = 1, b%ngl
                hi = b%psi(i,iquad)
                il = mf%imapl(1,i,1,iface)
                jl = mf%imapl(2,i,1,iface)
                kl = mf%imapl(3,i,1,iface)
                ip = G%intma(il,jl,kl,iel)

                rhs(1,ip) = rhs(1,ip) + wq*hi*flux_qu
                rhs(2,ip) = rhs(2,ip) + wq*hi*flux_qv
            end do

            if (ier > 0) then
                do i = 1, b%ngl
                    hi = b%psi(i,iquad)
                    ir = mf%imapr(1,i,1,iface)
                    jr = mf%imapr(2,i,1,iface)
                    kr = mf%imapr(3,i,1,iface)
                    ip = G%intma(ir,jr,kr,ier)

                    rhs(1,ip) = rhs(1,ip) - wq*hi*flux_qu
                    rhs(2,ip) = rhs(2,ip) - wq*hi*flux_qv
                end do
            end if
        end do !iface, iquad

    end subroutine create_rhs_laplacian_flux

    subroutine bcl_create_rhs_laplacian_flux(G, inp, b, mf, btp, bcl, rhs)

        use mod_grid,      only: grid
        use mod_basis,     only: basis
        use mod_input,     only: input
        use mod_face,      only: face_CS
        use mod_variables, only: btp_CS, bcl_CS

        implicit none

        type(grid),    intent(in)    :: G
        type(basis),   intent(in)    :: b
        type(input),   intent(in)    :: inp
        type(face_CS), intent(in)    :: mf
        type(btp_CS),  intent(in)    :: btp
        type(bcl_CS),  intent(in)    :: bcl

        real, intent(inout) :: rhs(2,G%npoin,inp%nlayers)

        real, dimension(2) :: qu_mean, qv_mean
        real, dimension(4,2) :: flux_uv_visc_face
        real, dimension(2) :: qul, qur
        real, dimension(2) :: qvl, qvr

        real :: nx, ny, nz
        real :: wq, un
        integer :: iface, i, j, il, jl, kl, ir, jr, kr, el, er
        integer :: iel, ier, ilocl, ilocr, ip
        integer :: iquad, jquad, ivar, k
        real :: flux_qu, flux_qv, hi, mul, mur, c_jump, alpha, beta, iflux
        real, dimension(5,b%ngl) :: ql, qr

        beta  = 0.5
        alpha = 1.0 - beta
        iflux = 0.0

        do concurrent(k = 1:inp%nlayers, iface=1:G%nface, iquad=1:b%ngl)

            if (G%face_type(iface) == 2) cycle

            iel = G%face(7,iface)
            ier = G%face(8,iface)

            nx = mf%normal_vector(1,iquad,1,iface)
            ny = mf%normal_vector(2,iquad,1,iface)

            il = mf%imapl(1,iquad,1,iface)
            jl = mf%imapl(2,iquad,1,iface)
            kl = mf%imapl(3,iquad,1,iface)
            ip = G%intma(il,jl,kl,iel)

            ql(1,iquad) = bcl%dpp_graduv(1,ip,k)
            ql(2,iquad) = bcl%dpp_graduv(2,ip,k)
            ql(3,iquad) = bcl%dpp_graduv(3,ip,k)
            ql(4,iquad) = bcl%dpp_graduv(4,ip,k)
            ql(5,iquad) = bcl%dpprime_visc(ip,k)

            if (ier > 0) then

                ir = mf%imapr(1,iquad,1,iface)
                jr = mf%imapr(2,iquad,1,iface)
                kr = mf%imapr(3,iquad,1,iface)
                ip = G%intma(ir,jr,kr,ier)

                qr(1,iquad) = bcl%dpp_graduv(1,ip,k)
                qr(2,iquad) = bcl%dpp_graduv(2,ip,k)
                qr(3,iquad) = bcl%dpp_graduv(3,ip,k)
                qr(4,iquad) = bcl%dpp_graduv(4,ip,k)
                qr(5,iquad) = bcl%dpprime_visc(ip,k)

            else
                qr(1,iquad) = ql(1,iquad)
                qr(2,iquad) = ql(2,iquad)
                qr(3,iquad) = ql(3,iquad)
                qr(4,iquad) = ql(4,iquad)
                qr(5,iquad) = ql(5,iquad)

                if (ier == -4) then

                    un = ql(1,iquad)*nx + ql(2,iquad)*ny
                    qr(1,iquad) = ql(1,iquad) - 2.0*un*nx
                    qr(2,iquad) = ql(2,iquad) - 2.0*un*ny

                    un = ql(3,iquad)*nx + ql(4,iquad)*ny
                    qr(3,iquad) = ql(3,iquad) - 2.0*un*nx
                    qr(4,iquad) = ql(4,iquad) - 2.0*un*ny
                end if
            end if

            do ivar = 1, 4
                flux_uv_visc_face(ivar,1) = ql(5,iquad)*btp%graduvb_face_ave(ivar,1,iquad,iface) + ql(ivar,iquad)
                flux_uv_visc_face(ivar,2) = qr(5,iquad)*btp%graduvb_face_ave(ivar,2,iquad,iface) + qr(ivar,iquad)
            end do

            qul(1) = flux_uv_visc_face(1,1)
            qul(2) = flux_uv_visc_face(2,1)
            qvl(1) = flux_uv_visc_face(3,1)
            qvl(2) = flux_uv_visc_face(4,1)

            qur(1) = flux_uv_visc_face(1,2)
            qur(2) = flux_uv_visc_face(2,2)
            qvr(1) = flux_uv_visc_face(3,2)
            qvr(2) = flux_uv_visc_face(4,2)

            qu_mean(1) = alpha*qul(1) + beta*qur(1)
            qu_mean(2) = alpha*qul(2) + beta*qur(2)
            qv_mean(1) = alpha*qvl(1) + beta*qvr(1)
            qv_mean(2) = alpha*qvl(2) + beta*qvr(2)

            wq = mf%jac_face(iquad,1,iface)

            flux_qu = (qu_mean(1) - iflux*qul(1))*nx + (qu_mean(2) - iflux*qul(2))*ny
            flux_qv = (qv_mean(1) - iflux*qvl(1))*nx + (qv_mean(2) - iflux*qvl(2))*ny

            do i = 1, b%ngl
                hi = b%psi(i,iquad)
                il = mf%imapl(1,i,1,iface)
                jl = mf%imapl(2,i,1,iface)
                kl = mf%imapl(3,i,1,iface)
                ip = G%intma(il,jl,kl,iel)

                rhs(1,ip,k) = rhs(1,ip,k) + wq*hi*flux_qu
                rhs(2,ip,k) = rhs(2,ip,k) + wq*hi*flux_qv
            end do !i

            if (ier > 0) then
                do i = 1, b%ngl
                    hi = b%psi(i,iquad)
                    ir = mf%imapr(1,i,1,iface)
                    jr = mf%imapr(2,i,1,iface)
                    kr = mf%imapr(3,i,1,iface)
                    ip = G%intma(ir,jr,kr,ier)

                    rhs(1,ip,k) = rhs(1,ip,k) - wq*hi*flux_qu
                    rhs(2,ip,k) = rhs(2,ip,k) - wq*hi*flux_qv
                end do !i
            end if !ier
        end do !k, iface, iquad

    end subroutine bcl_create_rhs_laplacian_flux

end module mod_laplacian_quad
