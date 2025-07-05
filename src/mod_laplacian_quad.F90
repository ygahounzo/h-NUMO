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
    public :: bcl_create_laplacian_v2, btp_create_laplacian_v2

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
        real :: graduv_face(4,2,ngl,nface)
        integer :: iface, il, jl, kl, ir, jr, kr, iel, ier, iquad, Iq, c_jump,k, ivar
        real, dimension(4, npoin) :: graduv
        real :: un, nx, ny, ul, ur, vl, vr, dpp

        Uk(1,:) = qb_df(3,:)/qb_df(1,:)
        Uk(2,:) = qb_df(4,:)/qb_df(1,:)

        ! Compute the auxilary LDG variable graduv
        call compute_gradient_uv(graduv, Uk)

        graduvb_ave = graduvb_ave + graduv

        ! Get the auxilary LDG variable graduv at the elements faces
        do iface=1,nface

            !Store Left and Right Side Variables
            iel=face(7,iface)
            ier=face(8,iface)
    
            !----------------------------Left Element
            do iquad = 1,ngl

                il=imapl(1,iquad,1,iface)
                jl=imapl(2,iquad,1,iface)
                kl=imapl(3,iquad,1,iface)
                Iq = intma(il,jl,kl,iel)

                graduv_face(:,1,iquad,iface) = graduv(:,Iq) 

                if (ier > 0 ) then

                    ir=imapr(1,iquad,1,iface)
                    jr=imapr(2,iquad,1,iface)
                    kr=imapr(3,iquad,1,iface)
                    Iq=intma(ir,jr,kr,ier)

                    graduv_face(:,2,iquad,iface) = graduv(:,Iq)

                else
                    graduv_face(:,2,iquad,iface) = graduv_face(:,1,iquad,iface)

                    if(ier == -4) then 

                        nx = normal_vector(1,iquad,1,iface)
                        ny = normal_vector(2,iquad,1,iface)

                        un = graduv_face(1,1,iquad,iface)*nx + graduv_face(2,1,iquad,iface)*ny
                        graduv_face(1,2,iquad,iface) = graduv_face(1,1,iquad,iface) - 2.0*un*nx
                        graduv_face(2,2,iquad,iface) = graduv_face(2,1,iquad,iface) - 2.0*un*ny

                        un = graduv_face(3,1,iquad,iface)*nx + graduv_face(4,1,iquad,iface)*ny
                        graduv_face(3,2,iquad,iface) = graduv_face(3,1,iquad,iface) - 2.0*un*nx
                        graduv_face(4,2,iquad,iface) = graduv_face(4,1,iquad,iface) - 2.0*un*ny

                    end if 
                end if
            end do 
        end do

        ! Precommunication step 
        call create_rhs_lap_precommunicator_df(graduv_face,4)

        ! Compute volume integral 
        call btp_compute_laplacian(rhs_lap,graduv)

        ! Postcommunication step 
        call create_rhs_lap_postcommunicator_df(graduv_face,4)

        graduvb_face_ave = graduvb_face_ave + graduv_face

        ! Compute interface integral 
        call create_rhs_laplacian_flux(rhs_lap,graduv_face)

        !RHS of viscosity 
        rhs_lap(1,:) = visc_mlswe*massinv(:)*rhs_lap(1,:)
        rhs_lap(2,:) = visc_mlswe*massinv(:)*rhs_lap(2,:)

    end subroutine btp_create_laplacian

    ! Compute the Laplacian barotropic horizontal viscosity term using the LDG method and
    ! using quadrature points
    subroutine btp_create_laplacian_v2(rhs_lap,qprime_df,qb_df)

        ! Barotropic horizontal viscosity

        real, dimension(2,npoin), intent(out) :: rhs_lap
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(2,npoin) :: Uk
        real, dimension(4,npoin_q) :: flux_uv_visc
        real :: flux_uv_visc_face(4,2,nq,nface)
        integer :: iface, il, jl, kl, ir, jr, kr, iel, ier, iquad, Iq, c_jump,k
        real, dimension(2,2, npoin_q) :: graduv
        real :: rhs_temp(2,npoin), un, nx, ny, ul, ur, vl, vr, dpp

        rhs_lap = 0.0
        flux_uv_visc = 0.0

        do k = 1,nlayers

            Uk(1,:) = qprime_df(2,:,k) + qb_df(3,:)/qb_df(1,:)
            Uk(2,:) = qprime_df(3,:,k) + qb_df(4,:)/qb_df(1,:)

            ! Compute the auxilary LDG variable graduv
            call compute_gradient_uv_q(graduv, Uk)

            flux_uv_visc(1,:) = flux_uv_visc(1,:) + dpprime_visc_q(:,k)*graduv(1,1,:)
            flux_uv_visc(2,:) = flux_uv_visc(2,:) + dpprime_visc_q(:,k)*graduv(1,2,:)
            flux_uv_visc(3,:) = flux_uv_visc(3,:) + dpprime_visc_q(:,k)*graduv(2,1,:)
            flux_uv_visc(4,:) = flux_uv_visc(4,:) + dpprime_visc_q(:,k)*graduv(2,2,:)

        end do

        ! Get the auxilary LDG variable graduv at the elements faces
        do iface=1,nface

            !Store Left and Right Side Variables
            iel=face(7,iface)
            ier=face(8,iface)

            !----------------------------Left Element
            do iquad = 1,nq

                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                Iq = intma_dg_quad(il,jl,kl,iel)

                flux_uv_visc_face(1,1,iquad,iface) = flux_uv_visc(1,Iq)
                flux_uv_visc_face(2,1,iquad,iface) = flux_uv_visc(2,Iq)
                flux_uv_visc_face(3,1,iquad,iface) = flux_uv_visc(3,Iq)
                flux_uv_visc_face(4,1,iquad,iface) = flux_uv_visc(4,Iq)

                if (ier > 0 ) then

                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    Iq=intma_dg_quad(ir,jr,kr,ier)

                    flux_uv_visc_face(1,2,iquad,iface) = flux_uv_visc(1,Iq)
                    flux_uv_visc_face(2,2,iquad,iface) = flux_uv_visc(2,Iq)
                    flux_uv_visc_face(3,2,iquad,iface) = flux_uv_visc(3,Iq)
                    flux_uv_visc_face(4,2,iquad,iface) = flux_uv_visc(4,Iq)

                else

                    flux_uv_visc_face(1,2,iquad,iface) = flux_uv_visc_face(1,1,iquad,iface)
                    flux_uv_visc_face(2,2,iquad,iface) = flux_uv_visc_face(2,1,iquad,iface)
                    flux_uv_visc_face(3,2,iquad,iface) = flux_uv_visc_face(3,1,iquad,iface)
                    flux_uv_visc_face(4,2,iquad,iface) = flux_uv_visc_face(4,1,iquad,iface)

                    if(ier == -4) then
                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)

                        un = flux_uv_visc(1,Iq)*nx + flux_uv_visc(2,Iq)*ny
                        flux_uv_visc_face(1,2,iquad,iface) = flux_uv_visc(1,Iq) - 2.0*un*nx
                        flux_uv_visc_face(2,2,iquad,iface) = flux_uv_visc(2,Iq) - 2.0*un*ny

                        un = flux_uv_visc(3,Iq)*nx + flux_uv_visc(4,Iq)*ny
                        flux_uv_visc_face(3,2,iquad,iface) = flux_uv_visc(3,Iq) - 2.0*un*nx
                        flux_uv_visc_face(4,2,iquad,iface) = flux_uv_visc(4,Iq) - 2.0*un*ny

                    end if
                end if
            end do
        end do

        ! Communicate within processors
        call create_communicator_quad(flux_uv_visc_face,4)
        ! Compute volume integral 
        call compute_laplacian_quad(rhs_temp,flux_uv_visc)
        ! Compute interface integral 
        call create_rhs_laplacian_flux_quad(rhs_temp,flux_uv_visc_face)
        !RHS of the viscosity 
        rhs_lap(1,:) = visc_mlswe*massinv(:)*rhs_temp(1,:)
        rhs_lap(2,:) = visc_mlswe*massinv(:)*rhs_temp(2,:)

    end subroutine btp_create_laplacian_v2

    ! Compute the Laplacian baroclinic horizontal viscosity term using the LDG method and
    ! using nodal points
    subroutine bcl_create_laplacian(rhs_lap)

        use mod_variables, only: dpprime_visc, dpp_graduv, graduvb_ave

        real, intent (out) :: rhs_lap(2,npoin,nlayers)

        integer :: k
        real :: rhs_temp(2,npoin)

        rhs_lap = 0.0

        do k = 1,nlayers

            call bcl_compute_laplacian(rhs_temp,k)
            call bcl_create_rhs_laplacian_flux(rhs_temp, k)

            rhs_lap(1,:,k) = visc_mlswe*massinv(:)*rhs_temp(1,:)
            rhs_lap(2,:,k) = visc_mlswe*massinv(:)*rhs_temp(2,:)

        end do

    end subroutine bcl_create_laplacian

    ! Compute the Laplacian baroclinic horizontal viscosity term using the LDG method and
    ! using quadrature points
    subroutine bcl_create_laplacian_v2(rhs_lap,qprime_df)

        ! Baroclinic horizontal viscosity

        use mod_variables, only: uvb_ave_df

        implicit none

        real, intent (out) :: rhs_lap(2,npoin,nlayers)
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        real, dimension(2,npoin) :: Uk
        real, dimension(4,npoin_q,nlayers) :: flux_uv_visc
        real, dimension(2,2,nq,nface) :: Uk_face
        real :: flux_uv_visc_face(4,2,nq,nface,nlayers)
        integer :: iface, il, jl, kl, ir, jr, kr, iel, ier, iquad, Iq, c_jump,k
        real, dimension(2,2, npoin_q) :: graduv
        real :: rhs_temp(2,npoin), un(nlayers), nx, ny, ul, ur, vl, vr, dpp

        rhs_lap = 0.0

        do k = 1,nlayers
            Uk(:,:) = qprime_df(2:3,:,k) + uvb_ave_df(:,:)

            ! Compute the auxilary LDG variable graduv
            call compute_gradient_uv_q(graduv, Uk)

            flux_uv_visc(1,:,k) = dpprime_visc_q(:,k)*graduv(1,1,:)
            flux_uv_visc(2,:,k) = dpprime_visc_q(:,k)*graduv(1,2,:)
            flux_uv_visc(3,:,k) = dpprime_visc_q(:,k)*graduv(2,1,:)
            flux_uv_visc(4,:,k) = dpprime_visc_q(:,k)*graduv(2,2,:)

        end do

        ! Get the auxilary LDG variable graduv at the elements faces
        do iface=1,nface

            !Store Left and Right Side Variables
            iel=face(7,iface)
            ier=face(8,iface)

            ! Left Element
            do iquad = 1,nq

                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                Iq = intma_dg_quad(il,jl,kl,iel)

                flux_uv_visc_face(1,1,iquad,iface,:) = flux_uv_visc(1,Iq,:)
                flux_uv_visc_face(2,1,iquad,iface,:) = flux_uv_visc(2,Iq,:)
                flux_uv_visc_face(3,1,iquad,iface,:) = flux_uv_visc(3,Iq,:)
                flux_uv_visc_face(4,1,iquad,iface,:) = flux_uv_visc(4,Iq,:)

                if (ier > 0 ) then

                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    Iq=intma_dg_quad(ir,jr,kr,ier)

                    flux_uv_visc_face(1,2,iquad,iface,:) = flux_uv_visc(1,Iq,:)
                    flux_uv_visc_face(2,2,iquad,iface,:) = flux_uv_visc(2,Iq,:)
                    flux_uv_visc_face(3,2,iquad,iface,:) = flux_uv_visc(3,Iq,:)
                    flux_uv_visc_face(4,2,iquad,iface,:) = flux_uv_visc(4,Iq,:)

                else
                    flux_uv_visc_face(1,2,iquad,iface,:) = flux_uv_visc_face(1,1,iquad,iface,:)
                    flux_uv_visc_face(2,2,iquad,iface,:) = flux_uv_visc_face(2,1,iquad,iface,:)
                    flux_uv_visc_face(3,2,iquad,iface,:) = flux_uv_visc_face(3,1,iquad,iface,:)
                    flux_uv_visc_face(4,2,iquad,iface,:) = flux_uv_visc_face(4,1,iquad,iface,:)

                    if(ier == -4) then
                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)

                        un = flux_uv_visc(1,Iq,:)*nx + flux_uv_visc(2,Iq,:)*ny
                        flux_uv_visc_face(1,2,iquad,iface,:) = flux_uv_visc(1,Iq,:) - 2.0*un*nx
                        flux_uv_visc_face(2,2,iquad,iface,:) = flux_uv_visc(2,Iq,:) - 2.0*un*ny

                        un = flux_uv_visc(3,Iq,:)*nx + flux_uv_visc(4,Iq,:)*ny
                        flux_uv_visc_face(3,2,iquad,iface,:) = flux_uv_visc(3,Iq,:) - 2.0*un*nx
                        flux_uv_visc_face(4,2,iquad,iface,:) = flux_uv_visc(4,Iq,:) - 2.0*un*ny

                    end if
                end if
            end do
        end do

        ! Communication within processors
        call bcl_create_communicator(flux_uv_visc_face,4,nlayers,nq)

        ! Commpute bcl viscosity
        do k = 1,nlayers

            call compute_laplacian_quad(rhs_temp,flux_uv_visc(:,:,k))
            call create_rhs_laplacian_flux_quad(rhs_temp,flux_uv_visc_face(:,:,:,:,k))

            !Store Solution
            rhs_lap(1,:,k) = visc_mlswe*massinv(:)*rhs_temp(1,:)
            rhs_lap(2,:,k) = visc_mlswe*massinv(:)*rhs_temp(2,:)
        end do

    end subroutine bcl_create_laplacian_v2

    subroutine btp_compute_laplacian(lap_q,grad_dpuvp)

        use mod_variables, only: pbprime_visc, btp_dpp_graduv

        implicit none

        real, intent(out) :: lap_q(2,npoin)
        real, dimension(4,npoin), intent(in) :: grad_dpuvp

        integer :: Iq, I, ip
        real :: wq, qq(4)
        
        lap_q = 0.0

        ! Compute the gradient of q
        do Iq = 1,npoin

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

    subroutine bcl_compute_laplacian(lap_q,k)

        use mod_variables, only: dpprime_visc, graduvb_ave, dpp_graduv

        implicit none

        real, intent(out) :: lap_q(2,npoin)
        integer, intent(in) :: k

        integer :: Iq, I, ip
        real :: wq, qq(4)
        
        lap_q = 0.0

        ! Compute the gradient of q
        do Iq = 1,npoin

            wq = wjac_df(Iq)

            qq(1) = dpprime_visc(Iq,k)*graduvb_ave(1,Iq) + dpp_graduv(1,Iq,k)
            qq(2) = dpprime_visc(Iq,k)*graduvb_ave(2,Iq) + dpp_graduv(2,Iq,k)
            qq(3) = dpprime_visc(Iq,k)*graduvb_ave(3,Iq) + dpp_graduv(3,Iq,k)
            qq(4) = dpprime_visc(Iq,k)*graduvb_ave(4,Iq) + dpp_graduv(4,Iq,k)

            do ip = 1,npts

                I = index_df(ip,Iq)
                lap_q(1,I) = lap_q(1,I) - wq*(dpsidx_df(ip,Iq)*qq(1) + dpsidy_df(ip,Iq)*qq(2))
                lap_q(2,I) = lap_q(2,I) - wq*(dpsidx_df(ip,Iq)*qq(3) + dpsidy_df(ip,Iq)*qq(4))

            end do
        end do

    end subroutine bcl_compute_laplacian

    subroutine create_rhs_laplacian_flux(rhs,gradq_face)

        use mod_basis, only: nq, psi
        use mod_variables, only: btp_graduv_dpp_face
    
        implicit none
    
        real, intent(inout) :: rhs(2,npoin)
        real, intent(in)    :: gradq_face(4,2,ngl,nface)
    
        real, dimension(2) :: qu_mean, qv_mean
        real, dimension(4,2) :: flux_uv_visc_face
        real, dimension(2) :: qul,qur
        real, dimension(2) :: qvl,qvr
        real nx, ny, nz
        real wq, un
        integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er
        integer iel, ier, ilocl, ilocr, ip
        integer iquad, jquad, ivar
        real :: flux_qu, flux_qv, hi, mul, mur,c_jump, alpha, beta

        beta = 0.5
        alpha = 1.0 - beta
      
        !Construct FVM-type Operators
        do iface=1,nface
            
            iel=face(7,iface)
            ier=face(8,iface)
    
            do iquad = 1,ngl

                do ivar = 1,4
                    flux_uv_visc_face(ivar,1) = btp_graduv_dpp_face(5,1,iquad,iface)* &
                                                gradq_face(ivar,1,iquad,iface) + &
                                                btp_graduv_dpp_face(ivar,1,iquad,iface)
                                            
                    flux_uv_visc_face(ivar,2) = btp_graduv_dpp_face(5,2,iquad,iface)* &
                                                gradq_face(ivar,2,iquad,iface) + &
                                                btp_graduv_dpp_face(ivar,2,iquad,iface)

                end do 

                nx = normal_vector(1,iquad,1,iface)
                ny = normal_vector(2,iquad,1,iface)

                qul(:) = flux_uv_visc_face(1:2,1)
                qvl(:) = flux_uv_visc_face(3:4,1)
                qur(:) = flux_uv_visc_face(1:2,2)
                qvr(:) = flux_uv_visc_face(3:4,2)

                ! The Flip-Flop flux of Cockburn & Shu 
                ! NOTE: beta=0.5 is the central flux
                qu_mean(:) = alpha*qul(:) + beta*qur(:)
                qv_mean(:) = alpha*qvl(:) + beta*qvr(:)

                wq=jac_face(iquad,1,iface)

                flux_qu = (qu_mean(1) - qul(1)*nx) + (qu_mean(2) - qul(2)*ny)
                flux_qv = (qv_mean(1) - qvl(1)*nx) + (qv_mean(2) - qvl(2)*ny)

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

                    if (ier > 0) then

                        ir=imapr(1,i,1,iface)
                        jr=imapr(2,i,1,iface)
                        kr=imapr(3,i,1,iface)

                        ip=intma(ir,jr,kr,ier)

                        !Update Flux
                        rhs(1,ip) = rhs(1,ip) - wq*hi*flux_qu
                        rhs(2,ip) = rhs(2,ip) - wq*hi*flux_qv

                    end if !ier
                end do !i
            end do !iquad
        end do !iface
    
    end subroutine create_rhs_laplacian_flux

    subroutine bcl_create_rhs_laplacian_flux(rhs,k)

        use mod_basis, only: nq, psi
        use mod_variables, only: graduv_dpp_face, graduvb_face_ave
    
        implicit none
    
        real, intent(inout) :: rhs(2,npoin)
        integer, intent(in) :: k
        real, dimension(2) :: qu_mean, qv_mean
        real, dimension(4,2) :: flux_uv_visc_face
        real, dimension(2) :: qul,qur
        real, dimension(2) :: qvl,qvr
    
        real nx, ny, nz
        real wq, un
        integer iface, i, j, il, jl, kl, ir, jr, kr, el, er
        integer iel, ier, ilocl, ilocr, ip
        integer iquad, jquad, ivar
        real :: flux_qu, flux_qv, hi, mul, mur,c_jump, alpha, beta

        beta = 0.5
        alpha = 1.0 - beta
      
        !Construct FVM-type Operators
        do iface=1,nface

            iel=face(7,iface)
            ier=face(8,iface)
    
            !----------------------------Left Element
            do iquad = 1,ngl

                do ivar = 1,4
                    flux_uv_visc_face(ivar,1) = graduv_dpp_face(5,1,iquad,iface,k)* &
                                                graduvb_face_ave(ivar,1,iquad,iface) + &
                                                graduv_dpp_face(ivar,1,iquad,iface,k)
                                            
                    flux_uv_visc_face(ivar,2) = graduv_dpp_face(5,2,iquad,iface,k)* &
                                                graduvb_face_ave(ivar,2,iquad,iface) + &
                                                graduv_dpp_face(ivar,2,iquad,iface,k)

                end do 

                nx = normal_vector(1,iquad,1,iface)
                ny = normal_vector(2,iquad,1,iface)

                qul(:) = flux_uv_visc_face(1:2,1)
                qvl(:) = flux_uv_visc_face(3:4,1)
                qur(:) = flux_uv_visc_face(1:2,2)
                qvr(:) = flux_uv_visc_face(3:4,2)

                ! The Flip-Flop flux of Cockburn & Shu 
                ! NOTE: beta=0.5 is the central flux
                qu_mean(:) = alpha*qul(:) + beta*qur(:)
                qv_mean(:) = alpha*qvl(:) + beta*qvr(:)

                wq = jac_face(iquad,1,iface)

                flux_qu = (qu_mean(1) - qul(1)*nx) + (qu_mean(2) - qul(2)*ny)
                flux_qv = (qv_mean(1) - qvl(1)*nx) + (qv_mean(2) - qvl(2)*ny)

                do i=1,ngl

                    hi = psi(i,iquad)

                    il=imapl(1,i,1,iface)
                    jl=imapl(2,i,1,iface)
                    kl=imapl(3,i,1,iface)
                    ip=intma(il,jl,kl,iel)
                    
                    !Update Flux
                    rhs(1,ip) = rhs(1,ip) + wq*hi*flux_qu
                    rhs(2,ip) = rhs(2,ip) + wq*hi*flux_qv

                    if (ier > 0) then

                        ir=imapr(1,i,1,iface)
                        jr=imapr(2,i,1,iface)
                        kr=imapr(3,i,1,iface)
                        ip=intma(ir,jr,kr,ier)

                        !Update Flux
                        rhs(1,ip) = rhs(1,ip) - wq*hi*flux_qu
                        rhs(2,ip) = rhs(2,ip) - wq*hi*flux_qv
                    end if !ier
                end do !i
            end do !iquad
        end do !iface

    end subroutine bcl_create_rhs_laplacian_flux

    subroutine compute_laplacian_quad(lap_q,grad_dpuvp)

        implicit none

        real, intent(out) :: lap_q(2,npoin)
        real, dimension(4,npoin_q), intent(in) :: grad_dpuvp

        integer :: Iq, I, ip
        real :: wq, u_visc, v_visc

        lap_q = 0.0

        ! Compute the gradient of q
        do Iq = 1,npoin_q

            wq = wjac(Iq)

            do ip = 1,npts

                I = indexq(ip,Iq)
                u_visc = dpsidx(ip,Iq)*grad_dpuvp(1,Iq) + dpsidy(ip,Iq)*grad_dpuvp(2,Iq)
                v_visc = dpsidx(ip,Iq)*grad_dpuvp(3,Iq) + dpsidy(ip,Iq)*grad_dpuvp(4,Iq)

                lap_q(1,I) = lap_q(1,I) - wq*u_visc
                lap_q(2,I) = lap_q(2,I) - wq*v_visc

            end do
        end do

    end subroutine compute_laplacian_quad

    subroutine create_rhs_laplacian_flux_quad(rhs,gradq_face)

        implicit none

        real, intent(inout) :: rhs(2,npoin)
        real, intent(in)    :: gradq_face(4,2,nq,nface)

        real, dimension(2) :: qu_mean, qv_mean
        real, dimension(2) :: qul,qur
        real, dimension(2) :: qvl,qvr
        real nx, ny, nz
        real wq
        integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er
        integer iel, ier, ilocl, ilocr, ip
        integer iquad, jquad
        real :: flux_qu, flux_qv, hi, mul, mur,c_jump, alpha, beta

        beta = 0.5
        alpha = 1.0 - beta

        !Construct FVM-type Operators
        do iface=1,nface
            !Store Left and Right Side Variables
            ilocl=face(5,iface)
            iel=face(7,iface)
            ier=face(8,iface)

            !----------------------------Left Element
            do iquad = 1,nq

                qul(:) = gradq_face(1:2,1,iquad,iface)
                qvl(:) = gradq_face(3:4,1,iquad,iface)
                qur(:) = gradq_face(1:2,2,iquad,iface)
                qvr(:) = gradq_face(3:4,2,iquad,iface)

                ! The Flip-Flop flux of Cockburn & Shu 
                ! NOTE: beta=0.5 is the central flux
                qu_mean(:) = alpha*qul(:) + beta*qur(:)
                qv_mean(:) = alpha*qvl(:) + beta*qvr(:)

                wq=jac_faceq(iquad,1,iface)

                !Store normals
                nx = normal_vector_q(1,iquad,1,iface)
                ny = normal_vector_q(2,iquad,1,iface)

                flux_qu = (qu_mean(1) - qul(1)*nx) + (qu_mean(2) - qul(2)*ny)
                flux_qv = (qv_mean(1) - qvl(1)*nx) + (qv_mean(2) - qvl(2)*ny)

                !  Do Gauss-Lobatto Integration
                do i=1,ngl

                    hi = psiq(i,iquad)

                    il=imapl(1,i,1,iface)
                    jl=imapl(2,i,1,iface)
                    kl=imapl(3,i,1,iface)
                    ip=intma(il,jl,kl,iel)

                    !Update Flux
                    rhs(1,ip) = rhs(1,ip) + wq*hi*flux_qu
                    rhs(2,ip) = rhs(2,ip) + wq*hi*flux_qv

                    if (ier > 0) then

                        ir=imapr(1,i,1,iface)
                        jr=imapr(2,i,1,iface)
                        kr=imapr(3,i,1,iface)
                        ip=intma(ir,jr,kr,ier)

                        !Update Flux
                        rhs(1,ip) = rhs(1,ip) - wq*hi*flux_qu
                        rhs(2,ip) = rhs(2,ip) - wq*hi*flux_qv
                    end if !ier
                end do !i
            end do !iquad
        end do !iface

    end subroutine create_rhs_laplacian_flux_quad

end module mod_laplacian_quad
