! ===========================================================================================================================
! This module contains the routines for the barotropic and baroclinic horizontal viscosity terms 
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: July 24, 2023
! ==========================================================================================================================

module mod_laplacian_quad

    use mod_grid, only: npoin, npoin_q, face, nface, intma_dg_quad, intma, mod_grid_get_face_ngl, mod_grid_get_face_nq
    use mod_face, only: imapl_q, imapr_q, normal_vector_q, jac_faceq, imapl, imapr, normal_vector, jac_face
    use mod_input, only: visc_mlswe, nlayers, xdims, ydims, nelx
    use mod_basis, only: ngl, nq, psiq, npts
    use mod_barotropic_terms, only: compute_gradient_uv
    use mod_metrics, only: massinv, jacq
    use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, psih_df, dpsidx_df,dpsidy_df, index_df, wjac_df
    use mod_variables, only: dpprime_visc_q

    implicit none

    public :: bcl_create_laplacian, btp_create_laplacian

    contains 

    subroutine btp_create_laplacian(rhs_lap,qb_df)

        ! Barotropic horizontal viscosity 

        use mod_variables, only: pbprime_visc, btp_dpp_graduv, graduvb_ave, graduvb_face_ave, btp_graduv_dpp_face

        real, intent (out) :: rhs_lap(2,npoin)
        real, dimension(4,npoin), intent(in) :: qb_df

        real, dimension(2,npoin) :: Uk
        real :: graduv_face(4,2,ngl,nface)
        integer :: iface, il, jl, kl, ir, jr, kr, iel, ier, iquad, Iq, c_jump,k, ivar
        real, dimension(4, npoin) :: graduv
        real :: un, nx, ny, ul, ur, vl, vr, dpp

        !rhs_lap = 0.0

        Uk(1,:) = qb_df(3,:)/qb_df(1,:)
        Uk(2,:) = qb_df(4,:)/qb_df(1,:)

        call compute_gradient_uv(graduv, Uk)

        graduvb_ave = graduvb_ave + graduv

        do iface=1,nface

            !-------------------------------------
            !Store Left and Right Side Variables
            !-------------------------------------
            iel=face(7,iface)
            ier=face(8,iface)
    
            !----------------------------Left Element
            do iquad = 1,ngl
                !Get Pointers
                il=imapl(1,iquad,1,iface)
                jl=imapl(2,iquad,1,iface)
                kl=imapl(3,iquad,1,iface)
                Iq = intma(il,jl,kl,iel)

                !Variables
                graduv_face(:,1,iquad,iface) = graduv(:,Iq) 

                if (ier > 0 ) then

                    !Get Pointers
                    ir=imapr(1,iquad,1,iface)
                    jr=imapr(2,iquad,1,iface)
                    kr=imapr(3,iquad,1,iface)
                    Iq=intma(ir,jr,kr,ier)

                    !Variables
                    graduv_face(:,2,iquad,iface) = graduv(:,Iq)

                else
                    !default values

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
        call btp_compute_laplacian_IBP(rhs_lap,graduv)

        ! Postcommunication step 
        call create_rhs_lap_postcommunicator_df(graduv_face,4)

        graduvb_face_ave = graduvb_face_ave + graduv_face

        ! Compute interface integral 
        call create_rhs_laplacian_flux_SIPG(rhs_lap,graduv_face)

        !RHS of viscosity 
        rhs_lap(1,:) = visc_mlswe*massinv(:)*rhs_lap(1,:)
        rhs_lap(2,:) = visc_mlswe*massinv(:)*rhs_lap(2,:)

    end subroutine btp_create_laplacian

    subroutine bcl_create_laplacian(rhs_lap)

        use mod_variables, only: dpprime_visc, dpp_graduv, graduvb_ave

        real, intent (out) :: rhs_lap(2,npoin,nlayers)

        integer :: k
        real :: rhs_temp(2,npoin)

        rhs_lap = 0.0

        do k = 1,nlayers

            call bcl_compute_laplacian_IBP(rhs_temp,k)

            call bcl_create_rhs_laplacian_flux_SIPG(rhs_temp, k)

            rhs_lap(1,:,k) = visc_mlswe*massinv(:)*rhs_temp(1,:)
            rhs_lap(2,:,k) = visc_mlswe*massinv(:)*rhs_temp(2,:)

        end do

    end subroutine bcl_create_laplacian

    subroutine btp_compute_laplacian_IBP(lap_q,grad_dpuvp)

        use mod_variables, only: pbprime_visc, btp_dpp_graduv
        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(2,npoin)
        real, dimension(4,npoin), intent(in) :: grad_dpuvp

        integer :: Iq, I, ip
        real :: wq, qq(4)
        
        !initialize
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

    end subroutine btp_compute_laplacian_IBP

    subroutine bcl_compute_laplacian_IBP(lap_q,k)

        use mod_variables, only: dpprime_visc, graduvb_ave, dpp_graduv

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(2,npoin)
        integer, intent(in) :: k

        integer :: Iq, I, ip
        real :: wq, qq(4)
        
        !initialize
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

    end subroutine bcl_compute_laplacian_IBP

    subroutine create_rhs_laplacian_flux_SIPG(rhs,gradq_face)

        use mod_basis, only: nq, psi
        use mod_variables, only: btp_graduv_dpp_face
    
        implicit none
    
        !global arrays
        real, intent(inout) :: rhs(2,npoin)
        real, intent(in)    :: gradq_face(4,2,ngl,nface)
    
        !local arrays
        real, dimension(2) :: qu_mean, qv_mean
        real, dimension(4,2) :: flux_uv_visc_face
    
        !local variables
        real nx, ny, nz
        real wq, un
        integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er
        integer iel, ier, ilocl, ilocr, ip
        integer iquad, jquad, ivar
    
        real :: flux_qu, flux_qv, hi, mul, mur,c_jump
      
        !Construct FVM-type Operators
        do iface=1,nface
            
            iel=face(7,iface)
            ier=face(8,iface)
    
            
            do iquad = 1,ngl

                do ivar = 1,4
                    flux_uv_visc_face(ivar,1) = btp_graduv_dpp_face(5,1,iquad,iface)*gradq_face(ivar,1,iquad,iface) + &
                                                            btp_graduv_dpp_face(ivar,1,iquad,iface)
                                            
                    flux_uv_visc_face(ivar,2) = btp_graduv_dpp_face(5,2,iquad,iface)*gradq_face(ivar,2,iquad,iface) + &
                                                            btp_graduv_dpp_face(ivar,2,iquad,iface)

                end do 

                nx = normal_vector(1,iquad,1,iface)
                ny = normal_vector(2,iquad,1,iface)


                qu_mean(:) = 0.5*(flux_uv_visc_face(1:2,1) + flux_uv_visc_face(1:2,2))
                qv_mean(:) = 0.5*(flux_uv_visc_face(3:4,1) + flux_uv_visc_face(3:4,2))

                wq=jac_face(iquad,1,iface)

                flux_qu = nx*qu_mean(1) + ny*qu_mean(2)
                flux_qv = nx*qv_mean(1) + ny*qv_mean(2)

                !---------------------------------
                !  Do Gauss-Lobatto Integration
                !---------------------------------

                do i=1,ngl

                    hi = psi(i,iquad)

                    !Pointers
                    il=imapl(1,i,1,iface)
                    jl=imapl(2,i,1,iface)
                    kl=imapl(3,i,1,iface)
                    
                    ip=intma(il,jl,kl,iel)
                    
                    !Update Flux
                    rhs(1,ip) = rhs(1,ip) + wq*hi*flux_qu
                    rhs(2,ip) = rhs(2,ip) + wq*hi*flux_qv

                    if (ier > 0) then
                        !Pointers
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
    
    end subroutine create_rhs_laplacian_flux_SIPG

    subroutine bcl_create_rhs_laplacian_flux_SIPG(rhs,k)

        use mod_basis, only: nq, psi
        use mod_variables, only: graduv_dpp_face, graduvb_face_ave
    
        implicit none
    
        !global arrays
        real, intent(inout) :: rhs(2,npoin)
        integer, intent(in) :: k
    
        !local arrays
        real, dimension(2) :: qu_mean, qv_mean
        real, dimension(4,2) :: flux_uv_visc_face
    
        !local variables
        real nx, ny, nz
        real wq, un
        integer iface, i, j, il, jl, kl, ir, jr, kr, el, er
        integer iel, ier, ilocl, ilocr, ip
        integer iquad, jquad, ivar
    
        real :: flux_qu, flux_qv, hi, mul, mur,c_jump
      
        !Construct FVM-type Operators

        do iface=1,nface

            iel=face(7,iface)
            ier=face(8,iface)
    
            !----------------------------Left Element
            do iquad = 1,ngl

                do ivar = 1,4
                    flux_uv_visc_face(ivar,1) = graduv_dpp_face(5,1,iquad,iface,k)*graduvb_face_ave(ivar,1,iquad,iface) + &
                                                            graduv_dpp_face(ivar,1,iquad,iface,k)
                                            
                    flux_uv_visc_face(ivar,2) = graduv_dpp_face(5,2,iquad,iface,k)*graduvb_face_ave(ivar,2,iquad,iface) + &
                                                            graduv_dpp_face(ivar,2,iquad,iface,k)

                end do 

                nx = normal_vector(1,iquad,1,iface)
                ny = normal_vector(2,iquad,1,iface)

                qu_mean(:) = 0.5*(flux_uv_visc_face(1:2,1) + flux_uv_visc_face(1:2,2))
                qv_mean(:) = 0.5*(flux_uv_visc_face(3:4,1) + flux_uv_visc_face(3:4,2))

                wq = jac_face(iquad,1,iface)

                flux_qu = nx*qu_mean(1) + ny*qu_mean(2)
                flux_qv = nx*qv_mean(1) + ny*qv_mean(2)

                do i=1,ngl

                    hi = psi(i,iquad)

                    !Pointers
                    il=imapl(1,i,1,iface)
                    jl=imapl(2,i,1,iface)
                    kl=imapl(3,i,1,iface)
                    
                    ip=intma(il,jl,kl,iel)
                    
                    !Update Flux
                    rhs(1,ip) = rhs(1,ip) + wq*hi*flux_qu
                    rhs(2,ip) = rhs(2,ip) + wq*hi*flux_qv

                    if (ier > 0) then
                        !Pointers
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

    end subroutine bcl_create_rhs_laplacian_flux_SIPG

end module mod_laplacian_quad
