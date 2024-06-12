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
    use mod_input, only: visc_mlswe, nlayers, mass_exact, xdims, ydims, nelx
    use mod_basis, only: ngl, nq, psiq, npts
    use mod_barotropic_terms, only: evaluate_quprime2, compute_gradient_uv, compute_gradient_uv_v1
    use mod_metrics, only: massinv, jacq
    use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac, psih_df, dpsidx_df,dpsidy_df, index_df, wjac_df
    use mod_variables, only: dpprime_visc_q

    implicit none

    public :: bcl_create_laplacian, btp_create_laplacian_v1, &
                btp_create_laplacian, bcl_create_laplacian_v1, btp_create_laplacian_v3, bcl_create_laplacian_v2

    contains 

    subroutine btp_create_laplacian_v1(rhs_lap,qprime,qprime_face,qb,qb_face)

        ! Barotropic horizontal viscosity 

        real, intent (out) :: rhs_lap(2,npoin)
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face
        real, dimension(4,npoin_q), intent(in) :: qb
        real, dimension(4,2,nq,nface), intent(in) :: qb_face


        real, dimension(2,npoin_q) :: Uk
        real, dimension(4,npoin_q) :: flux_uv_visc
        real, dimension(2,2,nq,nface) :: Uk_face
        real :: flux_uv_visc_face(4,2,nq,nface)
        integer :: iface, il, jl, kl, ir, jr, kr, iel, ier, iquad, Iq, c_jump,k
        real, dimension(2,2, npoin_q) :: graduv
        real :: rhs_temp(2,npoin), un, nx, ny

        rhs_lap = 0.0
        flux_uv_visc = 0.0

        do k = 1,nlayers
            Uk(1,:) = qprime(2,:,k) + qb(3,:)/qb(1,:)
            Uk(2,:) = qprime(3,:,k) + qb(4,:)/qb(1,:)
            Uk_face(1,:,:,:) = qprime_face(2,:,:,:,k) + qb_face(3,:,:,:)/qb_face(1,:,:,:)
            Uk_face(2,:,:,:) = qprime_face(3,:,:,:,k) + qb_face(4,:,:,:)/qb_face(1,:,:,:)

            call evaluate_quprime2(graduv, Uk, Uk_face)


            flux_uv_visc(1,:) = flux_uv_visc(1,:) + dpprime_visc_q(:,k)*graduv(1,1,:)
            flux_uv_visc(2,:) = flux_uv_visc(2,:) + dpprime_visc_q(:,k)*graduv(1,2,:)

            flux_uv_visc(3,:) = flux_uv_visc(3,:) + dpprime_visc_q(:,k)*graduv(2,1,:)
            flux_uv_visc(4,:) = flux_uv_visc(4,:) + dpprime_visc_q(:,k)*graduv(2,2,:)

        end do
            
        do iface=1,nface
            !-------------------------------------
            !Store Left and Right Side Variables
            !-------------------------------------
            iel=face(7,iface)
            ier=face(8,iface)
    
            !----------------------------Left Element
            do iquad = 1,nq
                !Get Pointers
                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                Iq = intma_dg_quad(il,jl,kl,iel)

                !Variables
                flux_uv_visc_face(1,1,iquad,iface) = flux_uv_visc(1,Iq)
                flux_uv_visc_face(2,1,iquad,iface) = flux_uv_visc(2,Iq)

                flux_uv_visc_face(3,1,iquad,iface) = flux_uv_visc(3,Iq)
                flux_uv_visc_face(4,1,iquad,iface) = flux_uv_visc(4,Iq)

                if (ier > 0 ) then

                    !Get Pointers
                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    Iq=intma_dg_quad(ir,jr,kr,ier)

                    !Variables
                    flux_uv_visc_face(1,2,iquad,iface) = flux_uv_visc(1,Iq)
                    flux_uv_visc_face(2,2,iquad,iface) = flux_uv_visc(2,Iq)

                    flux_uv_visc_face(3,2,iquad,iface) = flux_uv_visc(3,Iq)
                    flux_uv_visc_face(4,2,iquad,iface) = flux_uv_visc(4,Iq)

                else
                    !default values

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
        ! Compute volume intergral 
        call compute_laplacian_IBP_set2nc_quad(rhs_temp,flux_uv_visc)
        ! Compute interface integral 
        call create_rhs_laplacian_flux_SIPG_quad(rhs_temp,flux_uv_visc_face)
        ! RHS of the viscosity term 
        rhs_lap(:,:) = visc_mlswe*rhs_temp(:,:)

        !Store Solution
        rhs_lap(1,:) = massinv(:)*rhs_lap(1,:)
        rhs_lap(2,:) = massinv(:)*rhs_lap(2,:)

    end subroutine btp_create_laplacian_v1

    subroutine btp_create_laplacian_v2(rhs_lap,qprime_df,qb_df, dpprime, qprime_face, qb_face)

        ! Barotropic horizontal viscosity 

        real, intent (out) :: rhs_lap(2,npoin)
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(npoin_q, nlayers), intent(in) :: dpprime
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face
        real, dimension(4,2,nq,nface), intent(in) :: qb_face

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

            call compute_gradient_uv(graduv, Uk)

            flux_uv_visc(1,:) = flux_uv_visc(1,:) + dpprime_visc_q(:,k)*graduv(1,1,:)
            flux_uv_visc(2,:) = flux_uv_visc(2,:) + dpprime_visc_q(:,k)*graduv(1,2,:)

            flux_uv_visc(3,:) = flux_uv_visc(3,:) + dpprime_visc_q(:,k)*graduv(2,1,:)
            flux_uv_visc(4,:) = flux_uv_visc(4,:) + dpprime_visc_q(:,k)*graduv(2,2,:)

        end do
            
        do iface=1,nface

            !-------------------------------------
            !Store Left and Right Side Variables
            !-------------------------------------
            iel=face(7,iface)
            ier=face(8,iface)
    
            !----------------------------Left Element
            do iquad = 1,nq
                !Get Pointers
                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                Iq = intma_dg_quad(il,jl,kl,iel)

                !Variables
                flux_uv_visc_face(1,1,iquad,iface) = flux_uv_visc(1,Iq)
                flux_uv_visc_face(2,1,iquad,iface) = flux_uv_visc(2,Iq)

                flux_uv_visc_face(3,1,iquad,iface) = flux_uv_visc(3,Iq)
                flux_uv_visc_face(4,1,iquad,iface) = flux_uv_visc(4,Iq)

                if (ier > 0 ) then

                    !Get Pointers
                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    Iq=intma_dg_quad(ir,jr,kr,ier)

                    !Variables
                    flux_uv_visc_face(1,2,iquad,iface) = flux_uv_visc(1,Iq)
                    flux_uv_visc_face(2,2,iquad,iface) = flux_uv_visc(2,Iq)

                    flux_uv_visc_face(3,2,iquad,iface) = flux_uv_visc(3,Iq)
                    flux_uv_visc_face(4,2,iquad,iface) = flux_uv_visc(4,Iq)

                else
                    !default values

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
        call compute_laplacian_IBP_set2nc_quad(rhs_temp,flux_uv_visc)
        ! Compute interface integral 
        call create_rhs_laplacian_flux_SIPG_quad(rhs_temp,flux_uv_visc_face)
        !RHS of the viscosity 
        rhs_lap(1,:) = visc_mlswe*massinv(:)*rhs_temp(1,:)
        rhs_lap(2,:) = visc_mlswe*massinv(:)*rhs_temp(2,:)

    end subroutine btp_create_laplacian_v2


    subroutine btp_create_laplacian(rhs_lap,qprime_df,qb_df)

        ! Barotropic horizontal viscosity 

        use mod_variables, only: dpprime_visc

        real, intent (out) :: rhs_lap(2,npoin)
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df
        real, dimension(4,npoin), intent(in) :: qb_df

        real, dimension(2,npoin) :: Uk
        real, dimension(4,npoin) :: flux_uv_visc
        real :: flux_uv_visc_face(4,2,ngl,nface)
        integer :: iface, il, jl, kl, ir, jr, kr, iel, ier, iquad, Iq, c_jump,k
        real, dimension(4, npoin) :: graduv
        real :: rhs_temp(2,npoin), un, nx, ny, ul, ur, vl, vr, dpp

        rhs_lap = 0.0
        flux_uv_visc = 0.0

        do k = 1,nlayers

            Uk(1,:) = qprime_df(2,:,k) + qb_df(3,:)/qb_df(1,:)
            Uk(2,:) = qprime_df(3,:,k) + qb_df(4,:)/qb_df(1,:)

            call compute_gradient_uv_v1(graduv, Uk)

            flux_uv_visc(1,:) = flux_uv_visc(1,:) + dpprime_visc(:,k)*graduv(1,:)
            flux_uv_visc(2,:) = flux_uv_visc(2,:) + dpprime_visc(:,k)*graduv(2,:)

            flux_uv_visc(3,:) = flux_uv_visc(3,:) + dpprime_visc(:,k)*graduv(3,:)
            flux_uv_visc(4,:) = flux_uv_visc(4,:) + dpprime_visc(:,k)*graduv(4,:)

        end do
            
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
                flux_uv_visc_face(1,1,iquad,iface) = flux_uv_visc(1,Iq)
                flux_uv_visc_face(2,1,iquad,iface) = flux_uv_visc(2,Iq)

                flux_uv_visc_face(3,1,iquad,iface) = flux_uv_visc(3,Iq)
                flux_uv_visc_face(4,1,iquad,iface) = flux_uv_visc(4,Iq)

                if (ier > 0 ) then

                    !Get Pointers
                    ir=imapr(1,iquad,1,iface)
                    jr=imapr(2,iquad,1,iface)
                    kr=imapr(3,iquad,1,iface)
                    Iq=intma(ir,jr,kr,ier)

                    !Variables
                    flux_uv_visc_face(1,2,iquad,iface) = flux_uv_visc(1,Iq)
                    flux_uv_visc_face(2,2,iquad,iface) = flux_uv_visc(2,Iq)

                    flux_uv_visc_face(3,2,iquad,iface) = flux_uv_visc(3,Iq)
                    flux_uv_visc_face(4,2,iquad,iface) = flux_uv_visc(4,Iq)

                else
                    !default values

                    flux_uv_visc_face(1,2,iquad,iface) = flux_uv_visc_face(1,1,iquad,iface)
                    flux_uv_visc_face(2,2,iquad,iface) = flux_uv_visc_face(2,1,iquad,iface)

                    flux_uv_visc_face(3,2,iquad,iface) = flux_uv_visc_face(3,1,iquad,iface)
                    flux_uv_visc_face(4,2,iquad,iface) = flux_uv_visc_face(4,1,iquad,iface)

                    if(ier == -4) then 
                        nx = normal_vector(1,iquad,1,iface)
                        ny = normal_vector(2,iquad,1,iface)

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
        call create_communicator_df(flux_uv_visc_face,4)
        ! Compute volume integral 
        call compute_laplacian_IBP_v2(rhs_temp,flux_uv_visc)
        ! Compute interface integral 
        call create_rhs_laplacian_flux_SIPG(rhs_temp,flux_uv_visc_face)
        !RHS of viscosity 
        rhs_lap(1,:) = visc_mlswe*massinv(:)*rhs_temp(1,:)
        rhs_lap(2,:) = visc_mlswe*massinv(:)*rhs_temp(2,:)

    end subroutine btp_create_laplacian

    subroutine btp_create_laplacian_v3(rhs_lap,qb_df)

        ! Barotropic horizontal viscosity 

        use mod_variables, only: pbprime_visc, btp_dpp_graduv, graduvb_ave, graduvb_face_ave, btp_graduv_dpp_face

        real, intent (out) :: rhs_lap(2,npoin)
        real, dimension(4,npoin), intent(in) :: qb_df

        real, dimension(2,npoin) :: Uk
        real :: graduv_face(4,2,ngl,nface)
        integer :: iface, il, jl, kl, ir, jr, kr, iel, ier, iquad, Iq, c_jump,k, ivar
        real, dimension(4, npoin) :: graduv
        real :: rhs_temp(2,npoin), un, nx, ny, ul, ur, vl, vr, dpp

        rhs_lap = 0.0

        Uk(1,:) = qb_df(3,:)/qb_df(1,:)
        Uk(2,:) = qb_df(4,:)/qb_df(1,:)

        call compute_gradient_uv_v1(graduv, Uk)

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
        call btp_compute_laplacian_IBP(rhs_temp,graduv)
        ! Postcommunication step 
        call create_rhs_lap_postcommunicator_df(graduv_face,4)

        graduvb_face_ave = graduvb_face_ave + graduv_face
        ! Compute interface integral 
        call create_rhs_laplacian_flux_SIPG_v1(rhs_temp,graduv_face)
        !RHS of viscosity 
        rhs_lap(1,:) = visc_mlswe*massinv(:)*rhs_temp(1,:)
        rhs_lap(2,:) = visc_mlswe*massinv(:)*rhs_temp(2,:)

    end subroutine btp_create_laplacian_v3

    subroutine bcl_create_laplacian(rhs_lap,qprime,qprime_face)

        ! Baroclinic horizontal viscosity 

        use mod_variables, only: uvb_ave,uvb_face_ave

        implicit none

        real, intent (out) :: rhs_lap(2,npoin,nlayers)
        real, dimension(3,npoin_q,nlayers), intent(in) :: qprime
        real, dimension(3,2,nq,nface,nlayers), intent(in) :: qprime_face
        !real, dimension(2,npoin_q), intent(in) :: uvb
        !real, dimension(2,2,nq,nface), intent(in) :: uvb_face


        real, dimension(2,npoin_q) :: Uk
        real, dimension(4,npoin_q,nlayers) :: flux_uv_visc
        real, dimension(2,2,nq,nface) :: Uk_face
        real :: flux_uv_visc_face(4,2,nq,nface,nlayers)
        integer :: iface, il, jl, kl, ir, jr, kr, iel, ier, iquad, Iq, c_jump,k
        real, dimension(2,2, npoin_q) :: graduv
        real :: rhs_temp(2,npoin), un(nlayers), nx, ny, ul, ur, vl, vr, dpp

        rhs_lap = 0.0

        do k = 1,nlayers
            Uk(:,:) = qprime(2:3,:,k) + uvb_ave(:,:)
            Uk_face(:,:,:,:) = qprime_face(2:3,:,:,:,k) + uvb_face_ave(:,:,:,:)

            call evaluate_quprime2(graduv, Uk, Uk_face)
            
            flux_uv_visc(1,:,k) = dpprime_visc_q(:,k)*graduv(1,1,:)
            flux_uv_visc(2,:,k) = dpprime_visc_q(:,k)*graduv(1,2,:)

            flux_uv_visc(3,:,k) = dpprime_visc_q(:,k)*graduv(2,1,:)
            flux_uv_visc(4,:,k) = dpprime_visc_q(:,k)*graduv(2,2,:)

        end do
            
        do iface=1,nface
            !-------------------------------------
            !Store Left and Right Side Variables
            !-------------------------------------
            iel=face(7,iface)
            ier=face(8,iface)
    
            !----------------------------Left Element
            do iquad = 1,nq
                !Get Pointers
                il=imapl_q(1,iquad,1,iface)
                jl=imapl_q(2,iquad,1,iface)
                kl=imapl_q(3,iquad,1,iface)
                Iq = intma_dg_quad(il,jl,kl,iel)

                !Variables
                flux_uv_visc_face(1,1,iquad,iface,:) = flux_uv_visc(1,Iq,:)
                flux_uv_visc_face(2,1,iquad,iface,:) = flux_uv_visc(2,Iq,:)

                flux_uv_visc_face(3,1,iquad,iface,:) = flux_uv_visc(3,Iq,:)
                flux_uv_visc_face(4,1,iquad,iface,:) = flux_uv_visc(4,Iq,:)

                if (ier > 0 ) then

                    !Get Pointers
                    ir=imapr_q(1,iquad,1,iface)
                    jr=imapr_q(2,iquad,1,iface)
                    kr=imapr_q(3,iquad,1,iface)
                    Iq=intma_dg_quad(ir,jr,kr,ier)

                    !Variables
                    flux_uv_visc_face(1,2,iquad,iface,:) = flux_uv_visc(1,Iq,:)
                    flux_uv_visc_face(2,2,iquad,iface,:) = flux_uv_visc(2,Iq,:)

                    flux_uv_visc_face(3,2,iquad,iface,:) = flux_uv_visc(3,Iq,:)
                    flux_uv_visc_face(4,2,iquad,iface,:) = flux_uv_visc(4,Iq,:)

                else
                    !default values
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

            call compute_laplacian_IBP_set2nc_quad(rhs_temp,flux_uv_visc(:,:,k))

            call create_rhs_laplacian_flux_SIPG_quad(rhs_temp,flux_uv_visc_face(:,:,:,:,k))

            !Store Solution
            rhs_lap(1,:,k) = visc_mlswe*massinv(:)*rhs_temp(1,:)
            rhs_lap(2,:,k) = visc_mlswe*massinv(:)*rhs_temp(2,:)

        end do

    end subroutine bcl_create_laplacian


    subroutine bcl_create_laplacian_v1(rhs_lap,qprime_df)

        ! Baroclinic horizontal viscosity 

        use mod_variables, only: uvb_ave_df, dpprime_visc

        implicit none

        real, intent (out) :: rhs_lap(2,npoin,nlayers)
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df


        real, dimension(2,npoin) :: Uk
        real, dimension(4,npoin,nlayers) :: flux_uv_visc
        real, dimension(2,2,nq,nface) :: Uk_face
        real :: flux_uv_visc_face(4,2,ngl,nface,nlayers)
        integer :: iface, il, jl, kl, ir, jr, kr, iel, ier, iquad, Iq, c_jump,k
        real, dimension(4, npoin) :: graduv
        real :: rhs_temp(2,npoin), un(nlayers), nx, ny, ul, ur, vl, vr, dpp

        rhs_lap = 0.0

        do k = 1,nlayers
            Uk(:,:) = qprime_df(2:3,:,k) + uvb_ave_df(:,:)
            
            call compute_gradient_uv_v1(graduv, Uk)
            
            flux_uv_visc(1,:,k) = dpprime_visc(:,k)*graduv(1,:)
            flux_uv_visc(2,:,k) = dpprime_visc(:,k)*graduv(2,:)

            flux_uv_visc(3,:,k) = dpprime_visc(:,k)*graduv(3,:)
            flux_uv_visc(4,:,k) = dpprime_visc(:,k)*graduv(4,:)

        end do
            
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
                flux_uv_visc_face(1,1,iquad,iface,:) = flux_uv_visc(1,Iq,:)
                flux_uv_visc_face(2,1,iquad,iface,:) = flux_uv_visc(2,Iq,:)

                flux_uv_visc_face(3,1,iquad,iface,:) = flux_uv_visc(3,Iq,:)
                flux_uv_visc_face(4,1,iquad,iface,:) = flux_uv_visc(4,Iq,:)

                if (ier > 0 ) then

                    !Get Pointers
                    ir=imapr(1,iquad,1,iface)
                    jr=imapr(2,iquad,1,iface)
                    kr=imapr(3,iquad,1,iface)
                    Iq=intma(ir,jr,kr,ier)

                    !Variables
                    flux_uv_visc_face(1,2,iquad,iface,:) = flux_uv_visc(1,Iq,:)
                    flux_uv_visc_face(2,2,iquad,iface,:) = flux_uv_visc(2,Iq,:)

                    flux_uv_visc_face(3,2,iquad,iface,:) = flux_uv_visc(3,Iq,:)
                    flux_uv_visc_face(4,2,iquad,iface,:) = flux_uv_visc(4,Iq,:)

                else
                    !default values
                    flux_uv_visc_face(1,2,iquad,iface,:) = flux_uv_visc_face(1,1,iquad,iface,:)
                    flux_uv_visc_face(2,2,iquad,iface,:) = flux_uv_visc_face(2,1,iquad,iface,:)

                    flux_uv_visc_face(3,2,iquad,iface,:) = flux_uv_visc_face(3,1,iquad,iface,:)
                    flux_uv_visc_face(4,2,iquad,iface,:) = flux_uv_visc_face(4,1,iquad,iface,:)

                    if(ier == -4) then 
                        nx = normal_vector(1,iquad,1,iface)
                        ny = normal_vector(2,iquad,1,iface)

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

        ! Communicate within processors 
        call bcl_create_communicator(flux_uv_visc_face,4,nlayers,ngl)

        ! Compute bcl viscosity 
        do k = 1,nlayers

            call compute_laplacian_IBP_v2(rhs_temp,flux_uv_visc(:,:,k))

            call create_rhs_laplacian_flux_SIPG(rhs_temp,flux_uv_visc_face(:,:,:,:,k))

            !Store Solution
            rhs_lap(1,:,k) = visc_mlswe*massinv(:)*rhs_temp(1,:)
            rhs_lap(2,:,k) = visc_mlswe*massinv(:)*rhs_temp(2,:)

        end do

    end subroutine bcl_create_laplacian_v1

    subroutine bcl_create_laplacian_v2(rhs_lap)

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

    end subroutine bcl_create_laplacian_v2

    subroutine create_rhs_laplacian_flux_SIPG_quad(rhs,gradq_face)
    
        implicit none
    
        !global arrays
        real, intent(inout) :: rhs(2,npoin)
        real, intent(in)    :: gradq_face(4,2,nq,nface)
    
        !local arrays
        real, dimension(2) :: qu_mean, qv_mean
        real, dimension(2) :: qul,qur
        real, dimension(2) :: qvl,qvr
    
        !local variables
        real nx, ny, nz
        real wq
        integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er
        integer iel, ier, ilocl, ilocr, ip
        integer iquad, jquad
    
        real :: flux_qu, flux_qv, hi, mul, mur,c_jump
      
        !Construct FVM-type Operators
        do iface=1,nface
            !-------------------------------------
            !Store Left and Right Side Variables
            !-------------------------------------
            ilocl=face(5,iface)
            iel=face(7,iface)
            ier=face(8,iface)
    
            !----------------------------Left Element
            do iquad = 1,nq

                qu_mean(:) = 0.5*(gradq_face(1:2,1,iquad,iface) + gradq_face(1:2,2,iquad,iface))
                qv_mean(:) = 0.5*(gradq_face(3:4,1,iquad,iface) + gradq_face(3:4,2,iquad,iface))

                wq=jac_faceq(iquad,1,iface)

                !Store normals
                nx = normal_vector_q(1,iquad,1,iface)
                ny = normal_vector_q(2,iquad,1,iface)

                flux_qu = nx*qu_mean(1) + ny*qu_mean(2)
                flux_qv = nx*qv_mean(1) + ny*qv_mean(2)

                !---------------------------------
                !  Do Gauss-Lobatto Integration
                !---------------------------------
                !--------------Left Side------------------!

                do i=1,ngl

                    hi = psiq(i,iquad)

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
    
    end subroutine create_rhs_laplacian_flux_SIPG_quad

    subroutine compute_laplacian_IBP_set2nc_quad(lap_q,grad_dpuvp)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(2,npoin)
        real, dimension(4,npoin_q), intent(in) :: grad_dpuvp

        integer :: Iq, I, ip
        real :: wq, u_visc, v_visc
        
        !initialize
        lap_q = 0.0

        ! Compute the gradient of q

        do Iq = 1,npoin_q

            wq = wjac(Iq)

            do ip = 1,npts

                I = indexq(Iq,ip)

                u_visc = dpsidx(Iq,ip)*grad_dpuvp(1,Iq) + dpsidy(Iq,ip)*grad_dpuvp(2,Iq)
                v_visc = dpsidx(Iq,ip)*grad_dpuvp(3,Iq) + dpsidy(Iq,ip)*grad_dpuvp(4,Iq)

                lap_q(1,I) = lap_q(1,I) - wq*u_visc
                lap_q(2,I) = lap_q(2,I) - wq*v_visc

            end do
        end do

    end subroutine compute_laplacian_IBP_set2nc_quad

    subroutine compute_laplacian_IBP_v2(lap_q,grad_dpuvp)

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(2,npoin)
        real, dimension(4,npoin), intent(in) :: grad_dpuvp

        integer :: Iq, I, ip
        real :: wq, u_visc, v_visc
        
        !initialize
        lap_q = 0.0

        ! Compute the gradient of q

        do Iq = 1,npoin

            wq = wjac_df(Iq)

            do ip = 1,npts

                I = index_df(Iq,ip)

                u_visc = dpsidx_df(Iq,ip)*grad_dpuvp(1,Iq) + dpsidy_df(Iq,ip)*grad_dpuvp(2,Iq)
                v_visc = dpsidx_df(Iq,ip)*grad_dpuvp(3,Iq) + dpsidy_df(Iq,ip)*grad_dpuvp(4,Iq)

                lap_q(1,I) = lap_q(1,I) - wq*u_visc
                lap_q(2,I) = lap_q(2,I) - wq*v_visc

            end do
        end do

    end subroutine compute_laplacian_IBP_v2

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

                I = index_df(Iq,ip)

                lap_q(1,I) = lap_q(1,I) - wq*(dpsidx_df(Iq,ip)*qq(1) + dpsidy_df(Iq,ip)*qq(2))
                lap_q(2,I) = lap_q(2,I) - wq*(dpsidx_df(Iq,ip)*qq(3) + dpsidy_df(Iq,ip)*qq(4))

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

                I = index_df(Iq,ip)

                lap_q(1,I) = lap_q(1,I) - wq*(dpsidx_df(Iq,ip)*qq(1) + dpsidy_df(Iq,ip)*qq(2))
                lap_q(2,I) = lap_q(2,I) - wq*(dpsidx_df(Iq,ip)*qq(3) + dpsidy_df(Iq,ip)*qq(4))

            end do
        end do

    end subroutine bcl_compute_laplacian_IBP

    subroutine compute_laplacian_IBP(lap_q,grad_dpuvp)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad
        use mod_basis, only: nq, psi, psix, psiy, psiz, dpsi, dpsix, dpsiy, dpsiz
        use mod_metrics, only: ksi_x,ksi_y,ksi_z, eta_x,eta_y,eta_z, zeta_x,zeta_y,zeta_z, jac

        implicit none

        !Global Arrays
        real, intent(out) :: lap_q(2,npoin)
        real, dimension(4,npoin), intent(in) :: grad_dpuvp

        integer :: Iq, I, ip, e, iquad, jquad, n, m
        real :: wq, u_visc, v_visc
        real :: h_e, h_n
        real :: e_x, e_y, n_x, n_y
        real :: dhdx, dhdy, hi
        
        !initialize
        lap_q = 0.0

        ! Compute the gradient of q


        do e = 1,nelem

            do jquad = 1,ngly
                do iquad = 1,nglx
                    
                    Iq = intma(iquad,jquad,1,e)

                    wq = jac(iquad,jquad,1,e)

                    e_x = ksi_x(iquad,jquad,1,e); e_y = ksi_y(iquad,jquad,1,e); 
                    n_x = eta_x(iquad,jquad,1,e); n_y = eta_y(iquad,jquad,1,e);

                    ip = 0
                    
                    do m = 1, ngly 
                        do n = 1, nglx 

                            I = intma(n,m,1,e)
                
                            ! Xi derivatives
                            h_e = dpsix(n, iquad) * psiy(m, jquad)

                            ! Eta derivatives
                            h_n = psix(n, iquad) * dpsiy(m, jquad)
                
                            ! Pressure terms
                            dhdx = h_e * e_x + h_n * n_x
                            dhdy = h_e * e_y + h_n * n_y

                            u_visc = dhdx*grad_dpuvp(1,Iq) + dhdy*grad_dpuvp(2,Iq)
                            v_visc = dhdx*grad_dpuvp(3,Iq) + dhdy*grad_dpuvp(4,Iq)

                            lap_q(1,I) = lap_q(1,I) - wq*u_visc
                            lap_q(2,I) = lap_q(2,I) - wq*v_visc

                        end do !n
                    end do !m

                end do !iquad
            end do !jquad

        end do !e

    end subroutine compute_laplacian_IBP


    subroutine create_rhs_laplacian_flux_SIPG(rhs,gradq_face)

        use mod_basis, only: nq, psi
    
        implicit none
    
        !global arrays
        real, intent(inout) :: rhs(2,npoin)
        real, intent(in)    :: gradq_face(4,2,ngl,nface)
    
        !local arrays
        real, dimension(2) :: qu_mean, qv_mean
        real, dimension(2) :: qul,qur
        real, dimension(2) :: qvl,qvr
    
        !local variables
        real nx, ny, nz
        real wq
        integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er
        integer iel, ier, ilocl, ilocr, ip
        integer iquad, jquad
    
        real :: flux_qu, flux_qv, hi, mul, mur,c_jump
      
        !Construct FVM-type Operators
        do iface=1,nface
            !-------------------------------------
            !Store Left and Right Side Variables
            !-------------------------------------
            ilocl=face(5,iface)
            iel=face(7,iface)
            ier=face(8,iface)
    
            !----------------------------Left Element
            do iquad = 1,ngl

                qu_mean(:) = 0.5*(gradq_face(1:2,1,iquad,iface) + gradq_face(1:2,2,iquad,iface))
                qv_mean(:) = 0.5*(gradq_face(3:4,1,iquad,iface) + gradq_face(3:4,2,iquad,iface))

                wq=jac_face(iquad,1,iface)

                !Store normals
                nx = normal_vector(1,iquad,1,iface)
                ny = normal_vector(2,iquad,1,iface)

                flux_qu = nx*qu_mean(1) + ny*qu_mean(2)
                flux_qv = nx*qv_mean(1) + ny*qv_mean(2)

                !---------------------------------
                !  Do Gauss-Lobatto Integration
                !---------------------------------
                !--------------Left Side------------------!

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

    subroutine create_rhs_laplacian_flux_SIPG_v1(rhs,gradq_face)

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
    
    end subroutine create_rhs_laplacian_flux_SIPG_v1

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