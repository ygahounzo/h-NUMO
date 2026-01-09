! =================================================================================================
! This module contains the routines for the barotropic and baroclinic horizontal viscosity terms 
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: July 24, 2023
! =================================================================================================

module mod_laplacian_quad

    use mod_grid, only: npoin, npoin_q, face, nface, intma_dg_quad, intma, mod_grid_get_face_ngl, &
                        mod_grid_get_face_nq
    use mod_face, only: imapl_q, imapr_q, normal_vector_q, jac_faceq, imapl, imapr, normal_vector, &
                        jac_face
    use mod_input, only: visc_mlswe, nlayers, xdims, ydims, nelx
    use mod_basis, only: ngl, nq, psiq, npts
    use mod_barotropic_terms, only: compute_gradient_uv, compute_gradient_uv_q
    use mod_metrics, only: massinv, jacq
    use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, psih_df, dpsidx_df,dpsidy_df, &
                            index_df, wjac_df
    use mod_variables, only: dpprime_visc_q

    implicit none

    public :: bcl_create_laplacian, btp_create_laplacian

    contains 

    ! Compute the Laplacian barotropic horizontal viscosity term using the LDG method and
    ! using nodal points
    subroutine btp_create_laplacian(rhs_lap,qb_df)

        ! Barotropic horizontal viscosity 

        use mod_variables, only: pbprime_visc, btp_dpp_graduv, graduvb_ave, graduvb_face_ave, &
                                    btp_graduv_dpp_face

        real, intent (out) :: rhs_lap(2,npoin)
        real, dimension(4,npoin), intent(in) :: qb_df

        real, dimension(2,npoin) :: Uk
        real, dimension(4, npoin) :: graduv

        Uk(1,:) = qb_df(3,:)/qb_df(1,:)
        Uk(2,:) = qb_df(4,:)/qb_df(1,:)

        ! Compute the auxilary LDG variable graduv
        call compute_gradient_uv(graduv, Uk)

        graduvb_ave = graduvb_ave + graduv

        ! Precommunication step 
        call btp_lap_create_precommunicator(graduv,4)

        ! Compute volume integral 
        call btp_compute_laplacian(rhs_lap,graduv)

        ! Compute interface integral 
        call create_rhs_laplacian_flux(rhs_lap,graduv)

        ! Postcommunication step 
        call create_rhs_lap_postcommunicator_df(rhs_lap,4)

        !RHS of viscosity 
        rhs_lap(1,:) = visc_mlswe*rhs_lap(1,:)
        rhs_lap(2,:) = visc_mlswe*rhs_lap(2,:)

    end subroutine btp_create_laplacian

    ! Compute the Laplacian baroclinic horizontal viscosity term using the LDG method and
    ! using nodal points
    subroutine bcl_create_laplacian(rhs_lap)

        use mod_variables, only: dpprime_visc, dpp_graduv, graduvb_ave, dpprime_visc

        real, intent (out) :: rhs_lap(2,npoin,nlayers)

        integer :: k

        call bcl_lap_create_precommunicator(dpp_graduv,dpprime_visc)

        call bcl_compute_laplacian(rhs_lap)
        call bcl_create_rhs_laplacian_flux(rhs_lap)

        call bcl_create_rhs_lap_postcommunicator_df(rhs_lap)

        do k = 1,nlayers
            rhs_lap(1,:,k) = visc_mlswe*massinv(:)*rhs_lap(1,:,k)
            rhs_lap(2,:,k) = visc_mlswe*massinv(:)*rhs_lap(2,:,k)
        end do

    end subroutine bcl_create_laplacian

    subroutine btp_compute_laplacian(lap_q,grad_dpuvp)

        use mod_variables, only: pbprime_visc, btp_dpp_graduv

        implicit none

        real, intent(out) :: lap_q(2,npoin)
        real, dimension(4,npoin), intent(in) :: grad_dpuvp

        integer :: Iq, I, ip
        real :: wq, qq(4)
        
        lap_q = 0.0

        ! Compute the gradient of q
        do concurrent(Iq = 1:npoin)

            wq = wjac_df(Iq)

            qq(1) = pbprime_visc(Iq)*grad_dpuvp(1,Iq) + btp_dpp_graduv(1,Iq)
            qq(2) = pbprime_visc(Iq)*grad_dpuvp(2,Iq) + btp_dpp_graduv(2,Iq)
            qq(3) = pbprime_visc(Iq)*grad_dpuvp(3,Iq) + btp_dpp_graduv(3,Iq)
            qq(4) = pbprime_visc(Iq)*grad_dpuvp(4,Iq) + btp_dpp_graduv(4,Iq)

            do ip = 1,npts

                I = index_df(ip,Iq)
                lap_q(1,I) = lap_q(1,I) - wq*(dpsidx_df(ip,Iq)*qq(1) + dpsidy_df(ip,Iq)*qq(2))
                lap_q(2,I) = lap_q(2,I) - wq*(dpsidx_df(ip,Iq)*qq(3) + dpsidy_df(ip,Iq)*qq(4))

            end do
        end do

    end subroutine btp_compute_laplacian

    subroutine bcl_compute_laplacian(lap_q)

        use mod_variables, only: dpprime_visc, graduvb_ave, dpp_graduv
        use mod_input, only: nlayers, dry_cutoff
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe

        implicit none

        real, intent(out) :: lap_q(2,npoin,nlayers)

        integer :: Iq, I, ip, k
        real :: wq, qq(4)

        lap_q = 0.0

        do concurrent(k = 1:nlayers, Iq = 1:npoin)
            ! Compute the gradient of q
            wq = wjac_df(Iq)

            qq(1) = dpprime_visc(Iq,k)*graduvb_ave(1,Iq) + dpp_graduv(1,Iq,k)
            qq(2) = dpprime_visc(Iq,k)*graduvb_ave(2,Iq) + dpp_graduv(2,Iq,k)
            qq(3) = dpprime_visc(Iq,k)*graduvb_ave(3,Iq) + dpp_graduv(3,Iq,k)
            qq(4) = dpprime_visc(Iq,k)*graduvb_ave(4,Iq) + dpp_graduv(4,Iq,k)
            
            do ip = 1,npts

                I = index_df(ip,Iq)
                lap_q(1,I,k) = lap_q(1,I,k) - wq*(dpsidx_df(ip,Iq)*qq(1) + dpsidy_df(ip,Iq)*qq(2))
                lap_q(2,I,k) = lap_q(2,I,k) - wq*(dpsidx_df(ip,Iq)*qq(3) + dpsidy_df(ip,Iq)*qq(4))

            end do
        end do

    end subroutine bcl_compute_laplacian

    subroutine create_rhs_laplacian_flux(rhs,gradq)

        use mod_basis, only: nq, psi
        use mod_variables, only: btp_dpp_graduv, graduvb_face_ave, pbprime_visc
        use mod_grid, only: face_type
    
        implicit none
    
        real, intent(inout) :: rhs(2,npoin)
        real, intent(in)    :: gradq(4,npoin)
    
        real, dimension(2) :: qu_mean, qv_mean
        real, dimension(4,2) :: flux_uv_visc_face
        real, dimension(2) :: qul,qur
        real, dimension(2) :: qvl,qvr
        real nx, ny, nz
        real wq, un
        integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er
        integer iel, ier, ilocl, ilocr, ip
        integer iquad, jquad, ivar
        real :: flux_qu, flux_qv, hi, mul, mur,c_jump, alpha, beta, iflux
        real, dimension(4,ngl) :: ql, qr
        real, dimension(5,ngl) :: btp_ql, btp_qr

        beta = 0.5
        alpha = 1.0 - beta

        iflux = 0.0
      
        !Construct FVM-type Operators
        do concurrent(iface=1:nface, iquad=1:ngl)

            !Skip boundary faces
            if (face_type(iface) == 2) cycle
            
            iel=face(7,iface)
            ier=face(8,iface)

            nx = normal_vector(1,iquad,1,iface)
            ny = normal_vector(2,iquad,1,iface)

            il = imapl(1,iquad,1,iface)
            jl = imapl(2,iquad,1,iface)
            kl = imapl(3,iquad,1,iface)
            ip = intma(il,jl,kl,iel)

            ! ql(:,iquad) = gradq(:,ip)
            ! btp_ql(1:4,iquad) = btp_dpp_graduv(1:4,ip)
            ! btp_ql(5,iquad) = pbprime_visc(ip)
            ql(1,iquad) = gradq(1,ip)
            ql(2,iquad) = gradq(2,ip)
            ql(3,iquad) = gradq(3,ip)
            ql(4,iquad) = gradq(4,ip)

            btp_ql(1,iquad) = btp_dpp_graduv(1,ip)
            btp_ql(2,iquad) = btp_dpp_graduv(2,ip)
            btp_ql(3,iquad) = btp_dpp_graduv(3,ip)
            btp_ql(4,iquad) = btp_dpp_graduv(4,ip)
            btp_ql(5,iquad) = pbprime_visc(ip)

            if (ier > 0) then

                ir = imapr(1,iquad,1,iface)
                jr = imapr(2,iquad,1,iface)
                kr = imapr(3,iquad,1,iface)
                ip = intma(ir,jr,kr,ier)

                ! qr(:,iquad) = gradq(:,ip)
                ! btp_qr(1:4,iquad) = btp_dpp_graduv(1:4,ip)
                ! btp_qr(5,iquad) = pbprime_visc(ip)
                qr(1,iquad) = gradq(1,ip)
                qr(2,iquad) = gradq(2,ip)
                qr(3,iquad) = gradq(3,ip)
                qr(4,iquad) = gradq(4,ip)

                btp_qr(1,iquad) = btp_dpp_graduv(1,ip)
                btp_qr(2,iquad) = btp_dpp_graduv(2,ip)
                btp_qr(3,iquad) = btp_dpp_graduv(3,ip)
                btp_qr(4,iquad) = btp_dpp_graduv(4,ip)
                btp_qr(5,iquad) = pbprime_visc(ip)

            else
                ! qr(:,iquad) = ql(:,iquad)
                ! btp_qr(:,iquad) = btp_ql(:,iquad)
                qr(1,iquad) = ql(1,iquad)
                qr(2,iquad) = ql(2,iquad)
                qr(3,iquad) = ql(3,iquad)
                qr(4,iquad) = ql(4,iquad)

                if(ier == -4) then

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

            ! graduvb_face_ave(:,1,iquad,iface) = graduvb_face_ave(:,1,iquad,iface) + ql(:,iquad)
            graduvb_face_ave(1,1,iquad,iface) = graduvb_face_ave(1,1,iquad,iface) + ql(1,iquad)
            graduvb_face_ave(2,1,iquad,iface) = graduvb_face_ave(2,1,iquad,iface) + ql(2,iquad)
            graduvb_face_ave(3,1,iquad,iface) = graduvb_face_ave(3,1,iquad,iface) + ql(3,iquad)
            graduvb_face_ave(4,1,iquad,iface) = graduvb_face_ave(4,1,iquad,iface) + ql(4,iquad)

            ! graduvb_face_ave(:,2,iquad,iface) = graduvb_face_ave(:,2,iquad,iface) + qr(:,iquad)
            graduvb_face_ave(1,2,iquad,iface) = graduvb_face_ave(1,2,iquad,iface) + qr(1,iquad)
            graduvb_face_ave(2,2,iquad,iface) = graduvb_face_ave(2,2,iquad,iface) + qr(2,iquad)
            graduvb_face_ave(3,2,iquad,iface) = graduvb_face_ave(3,2,iquad,iface) + qr(3,iquad)
            graduvb_face_ave(4,2,iquad,iface) = graduvb_face_ave(4,2,iquad,iface) + qr(4,iquad)

            do ivar = 1,4
                flux_uv_visc_face(ivar,1) = btp_ql(5,iquad)* &
                                            ql(ivar,iquad) + &
                                            btp_ql(ivar,iquad)

                flux_uv_visc_face(ivar,2) = btp_qr(5,iquad)* &
                                            qr(ivar,iquad) + &
                                            btp_qr(ivar,iquad)

            end do 

            qul(1) = flux_uv_visc_face(1,1)
            qul(2) = flux_uv_visc_face(2,1)
            qvl(1) = flux_uv_visc_face(3,1)
            qvl(2) = flux_uv_visc_face(4,1)

            qur(1) = flux_uv_visc_face(1,2)
            qur(2) = flux_uv_visc_face(2,2)
            qvr(1) = flux_uv_visc_face(3,2)
            qvr(2) = flux_uv_visc_face(4,2)

            ! The Flip-Flop flux of Cockburn & Shu 
            ! NOTE: beta=0.5 is the central flux
            qu_mean(1) = alpha*qul(1) + beta*qur(1)
            qu_mean(2) = alpha*qul(2) + beta*qur(2)
            qv_mean(1) = alpha*qvl(1) + beta*qvr(1)
            qv_mean(2) = alpha*qvl(2) + beta*qvr(2)

            wq=jac_face(iquad,1,iface)

            flux_qu = (qu_mean(1) - iflux*qul(1))*nx + (qu_mean(2) - iflux*qul(2))*ny
            flux_qv = (qv_mean(1) - iflux*qvl(1))*nx + (qv_mean(2) - iflux*qvl(2))*ny

            !  Do Gauss-Lobatto Integration
            do i=1,ngl

                hi = psi(i,iquad)

                il=imapl(1,i,1,iface)
                jl=imapl(2,i,1,iface)
                kl=imapl(3,i,1,iface)
                ip=intma(il,jl,kl,iel)
                
                !Update Flux
                rhs(1,ip) = rhs(1,ip) + wq*hi*flux_qu
                rhs(2,ip) = rhs(2,ip) + wq*hi*flux_qv
            enddo

            if (ier > 0) then
                do i=1,ngl

                    hi = psi(i,iquad)
                    ir=imapr(1,i,1,iface)
                    jr=imapr(2,i,1,iface)
                    kr=imapr(3,i,1,iface)

                    ip=intma(ir,jr,kr,ier)

                    !Update Flux
                    rhs(1,ip) = rhs(1,ip) - wq*hi*flux_qu
                    rhs(2,ip) = rhs(2,ip) - wq*hi*flux_qv

                end do
            end if
        end do !iface, iquad

    end subroutine create_rhs_laplacian_flux

    subroutine bcl_create_rhs_laplacian_flux(rhs)

        use mod_basis, only: nq, psi
        use mod_variables, only: graduv_dpp_face, graduvb_face_ave, dpp_graduv, dpprime_visc
        use mod_input, only: nlayers, dry_cutoff
        use mod_constants, only: gravity
        use mod_initial, only: alpha_mlswe
        use mod_grid, only: face_type
    
        implicit none
    
        real, intent(inout) :: rhs(2,npoin,nlayers)
        real, dimension(2) :: qu_mean, qv_mean
        real, dimension(4,2) :: flux_uv_visc_face
        real, dimension(2) :: qul,qur
        real, dimension(2) :: qvl,qvr
    
        real nx, ny, nz
        real wq, un
        integer iface, i, j, il, jl, kl, ir, jr, kr, el, er
        integer iel, ier, ilocl, ilocr, ip
        integer iquad, jquad, ivar, k
        real :: flux_qu, flux_qv, hi, mul, mur,c_jump, alpha, beta, iflux
        real, dimension(5,ngl) :: ql, qr

        beta = 0.5
        alpha = 1.0 - beta
        iflux = 0.0

        do concurrent(k = 1:nlayers, iface=1:nface, iquad=1:ngl)

            !Skip boundary faces
            if (face_type(iface) == 2) cycle

            iel=face(7,iface)
            ier=face(8,iface)

            nx = normal_vector(1,iquad,1,iface)
            ny = normal_vector(2,iquad,1,iface)

            il = imapl(1,iquad,1,iface)
            jl = imapl(2,iquad,1,iface)
            kl = imapl(3,iquad,1,iface)
            ip = intma(il,jl,kl,iel)

            ql(1,iquad) = dpp_graduv(1,ip,k)
            ql(2,iquad) = dpp_graduv(2,ip,k)
            ql(3,iquad) = dpp_graduv(3,ip,k)
            ql(4,iquad) = dpp_graduv(4,ip,k)
            ql(5,iquad) = dpprime_visc(ip,k)

            if (ier > 0) then

                ir = imapr(1,iquad,1,iface)
                jr = imapr(2,iquad,1,iface)
                kr = imapr(3,iquad,1,iface)
                ip = intma(ir,jr,kr,ier)

                qr(1,iquad) = dpp_graduv(1,ip,k)
                qr(2,iquad) = dpp_graduv(2,ip,k)
                qr(3,iquad) = dpp_graduv(3,ip,k)
                qr(4,iquad) = dpp_graduv(4,ip,k)
                qr(5,iquad) = dpprime_visc(ip,k)

            else
                qr(1,iquad) = ql(1,iquad)
                qr(2,iquad) = ql(2,iquad)
                qr(3,iquad) = ql(3,iquad)
                qr(4,iquad) = ql(4,iquad)
                qr(5,iquad) = ql(5,iquad)

                if(ier == -4) then 

                    un = ql(1,iquad)*nx + ql(2,iquad)*ny
                    qr(1,iquad) = ql(1,iquad) - 2.0*un*nx
                    qr(2,iquad) = ql(2,iquad) - 2.0*un*ny

                    un = ql(3,iquad)*nx + ql(4,iquad)*ny
                    qr(3,iquad) = ql(3,iquad) - 2.0*un*nx
                    qr(4,iquad) = ql(4,iquad) - 2.0*un*ny
                end if
            endif 

            do ivar = 1,4
                flux_uv_visc_face(ivar,1) = ql(5,iquad)* &
                                            graduvb_face_ave(ivar,1,iquad,iface) + &
                                            ql(ivar,iquad)

                flux_uv_visc_face(ivar,2) = qr(5,iquad)* &
                                            graduvb_face_ave(ivar,2,iquad,iface) + &
                                            qr(ivar,iquad)

            end do 

            qul(1) = flux_uv_visc_face(1,1)
            qul(2) = flux_uv_visc_face(2,1)
            qvl(1) = flux_uv_visc_face(3,1)
            qvl(2) = flux_uv_visc_face(4,1)

            qur(1) = flux_uv_visc_face(1,2)
            qur(2) = flux_uv_visc_face(2,2)
            qvr(1) = flux_uv_visc_face(3,2)
            qvr(2) = flux_uv_visc_face(4,2)

            ! The Flip-Flop flux of Cockburn & Shu 
            ! NOTE: beta=0.5 is the central flux
            qu_mean(1) = alpha*qul(1) + beta*qur(1)
            qu_mean(2) = alpha*qul(2) + beta*qur(2)
            qv_mean(1) = alpha*qvl(1) + beta*qvr(1)
            qv_mean(2) = alpha*qvl(2) + beta*qvr(2)

            wq = jac_face(iquad,1,iface)

            flux_qu = (qu_mean(1) - iflux*qul(1))*nx + (qu_mean(2) - iflux*qul(2))*ny
            flux_qv = (qv_mean(1) - iflux*qvl(1))*nx + (qv_mean(2) - iflux*qvl(2))*ny

            do i=1,ngl

                hi = psi(i,iquad)

                il=imapl(1,i,1,iface)
                jl=imapl(2,i,1,iface)
                kl=imapl(3,i,1,iface)
                ip=intma(il,jl,kl,iel)
                
                !Update Flux
                rhs(1,ip,k) = rhs(1,ip,k) + wq*hi*flux_qu
                rhs(2,ip,k) = rhs(2,ip,k) + wq*hi*flux_qv
            end do !i

            if (ier > 0) then
                do i=1,ngl

                    hi = psi(i,iquad)
                    ir=imapr(1,i,1,iface)
                    jr=imapr(2,i,1,iface)
                    kr=imapr(3,i,1,iface)
                    ip=intma(ir,jr,kr,ier)

                    !Update Flux
                    rhs(1,ip,k) = rhs(1,ip,k) - wq*hi*flux_qu
                    rhs(2,ip,k) = rhs(2,ip,k) - wq*hi*flux_qv
                end do !i
            end if ! ier
        end do !k, iface, iquad

    end subroutine bcl_create_rhs_laplacian_flux

end module mod_laplacian_quad
