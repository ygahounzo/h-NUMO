! ===========================================================================================================================
! This module contains the routines for the baroclinic flux terms
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: March 27, 2023
! ==========================================================================================================================

module mod_layer_terms
    
    use mod_grid, only: npoin_q, nface, intma_dg_quad, face, npoin
    use mod_basis, only: nqx, nqy, nqz, nq
    use mod_input, only: nlayers
        
    implicit none

    public :: &
            shear_stress_system, layer_mom_boundary_df,  &
            evaluate_mom, velocity_df, evaluate_mom_face, &
            evaluate_bcl, evaluate_bcl_v1, &
            interpolate_qprime, evaluate_consistency_face, &
            extract_qprime_df_face, extract_dprime_df_face, interpolate_dpp

    contains

    subroutine interpolate_dpp(dprimeq,dprime_df)

        ! This routine interpolate layer mass from dofs to quad points

        use mod_basis, only: npts
        use mod_grid, only: npoin, npoin_q, intma, intma_dg_quad
        use mod_input, only: nlayers
        use mod_initial, only: psih, indexq

        implicit none
        real, dimension(npoin_q,nlayers), intent(inout) :: dprimeq
        real, dimension(npoin,nlayers), intent(inout) :: dprime_df

        integer :: k, Iq,I, ip
        real :: hi

        dprimeq = 0.0

        ! do k = 1, nlayers

        do Iq = 1,npoin_q
            do ip = 1,npts

                I = indexq(ip,Iq)
                hi = psih(ip,Iq)
                dprimeq(Iq,:) = dprimeq(Iq,:) + dprime_df(I,:)*hi

            end do
        end do

    end subroutine interpolate_dpp

    subroutine evaluate_consistency_face(mass_deficit_mass_face,dprime_df)

        ! This routines extracts the layer mass face values 

        use mod_grid, only:  npoin_q, npoin, intma, nface, face
        use mod_input, only: nlayers
        use mod_face, only: imapl, imapr
        use mod_variables, only: sum_layer_mass_flux_face, btp_mass_flux_face_ave
        use mod_initial, only: pbprime_face
        use mod_basis, only: ngl, psiq

        implicit none
        real, intent(out) :: mass_deficit_mass_face(2,2,nq,nface,nlayers)
        real, intent(in) :: dprime_df(npoin,nlayers)
    
        ! Local variables
        integer :: iface, iquad, k, il, jl, kl, el, er, ir, jr, kr, I, n
        real :: qprime_l, qprime_r, weights_face_l, weights_face_r, hi

        mass_deficit_mass_face = 0.0

        do k = 1,nlayers

            do iface = 1, nface

                !Store Left Side Variables

                el = face(7,iface)
                er = face(8,iface)

                do iquad = 1, nq
                    
                    qprime_l = 0.0; qprime_r = 0.0

                    do n = 1,ngl
                        hi = psiq(n,iquad)

                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)

                        I = intma(il,jl,kl,el)
                        qprime_l = qprime_l + hi*dprime_df(I,k)
                    enddo

                    if(er > 0) then
                        do n = 1,ngl
                            hi = psiq(n,iquad)
                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)
                           qprime_r = qprime_r + hi*dprime_df(I,k)
                        enddo
                    else
                        qprime_r = qprime_l
                    endif
  
                    weights_face_l = qprime_l / pbprime_face(1,iquad,iface)
                    weights_face_r = qprime_r / pbprime_face(2,iquad,iface)

                    mass_deficit_mass_face(1,1,iquad,iface,k) = weights_face_l*(btp_mass_flux_face_ave(1,iquad,iface) &
                                                                 - sum_layer_mass_flux_face(1,iquad,iface))
                    mass_deficit_mass_face(2,1,iquad,iface,k) = weights_face_l*(btp_mass_flux_face_ave(2,iquad,iface) &
                                                                 - sum_layer_mass_flux_face(2,iquad,iface))

                    mass_deficit_mass_face(1,2,iquad,iface,k) = weights_face_r*(btp_mass_flux_face_ave(1,iquad,iface) &
                                                                 - sum_layer_mass_flux_face(1,iquad,iface))
                    mass_deficit_mass_face(2,2,iquad,iface,k) = weights_face_r*(btp_mass_flux_face_ave(2,iquad,iface) &
                                                                 - sum_layer_mass_flux_face(2,iquad,iface))
                end do 
            end do
        end do

        call bcl_create_communicator(mass_deficit_mass_face,2,nlayers,nq)
    
    end subroutine evaluate_consistency_face

    subroutine velocity_df(q_df, qb_df)

        use mod_grid, only : npoin, intma, face, nface 
        use mod_basis, only : ngl
        use mod_input, only: nlayers
        use mod_face, only: imapl

        implicit none
        
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(4,npoin), intent(in) :: qb_df
        
        real, dimension(nlayers) :: a, b, c
        real, dimension(nlayers,2) :: r
        real, dimension(nlayers) :: weight
        real :: eps = 1.0e-12
        real :: ubar, vbar, mult, dpbar
        integer :: I, k, el, er, n , iface, il, jl, kl, Iq
        real :: uv_df(2,npoin,nlayers)

        uv_df = 0.0

        do k = 1,nlayers
            uv_df(1,:,k) = q_df(2,:,k) / q_df(1,:,k)
            uv_df(2,:,k) = q_df(3,:,k) / q_df(1,:,k)
        end do

        do I = 1, npoin
            ubar = 0.0
            vbar = 0.0
            
            do k = 1, nlayers
                ubar = ubar + uv_df(1,I,k) * q_df(1,I,k)
                vbar = vbar + uv_df(2,I,k) * q_df(1,I,k)
            end do

            if(qb_df(1,I) > 0.0) then
                
                ubar = ubar / qb_df(1,I)
                vbar = vbar / qb_df(1,I)

                do k = 1, nlayers

                    uv_df(1,I,k) = uv_df(1,I,k) - ubar + qb_df(3,I)/qb_df(1,I)
                    uv_df(2,I,k) = uv_df(2,I,k) - vbar + qb_df(4,I)/qb_df(1,I)
                end do

            else
                uv_df(:,I,:) = 0.0
            end if
        end do

        do k = 1,nlayers
            q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
            q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
        end do

    end subroutine velocity_df

    subroutine evaluate_bcl(qprime_df_face, q_df, qprime_df, qb_df)

        use mod_grid, only : npoin, intma, face, nface 
        use mod_basis, only : ngl
        use mod_input, only: nlayers
        use mod_initial, only: pbprime_df

        implicit none
        
        real, dimension(3,npoin,nlayers), intent(out) :: qprime_df
        real, dimension(3,2,ngl,nface,nlayers), intent(out) :: qprime_df_face
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        
        integer :: I, k, el, er, n , iface, il, jl, kl, Iq
        real :: uv_df(2,npoin,nlayers), one_plus_eta_temp(npoin)

        call extract_velocity(uv_df, q_df, qb_df)

        one_plus_eta_temp = 0.0

        do k = 1,nlayers
            q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
            q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)

            one_plus_eta_temp(:) = one_plus_eta_temp(:) + q_df(1,:,k)
        end do

        one_plus_eta_temp(:) = one_plus_eta_temp(:) / pbprime_df(:)

        call extract_velocity(uv_df, q_df, qb_df)

        do k = 1,nlayers
            qprime_df(1,:,k) = q_df(1,:,k) / one_plus_eta_temp(:)
            qprime_df(2,:,k) = uv_df(1,:,k) - qb_df(3,:)/qb_df(1,:)
            qprime_df(3,:,k) = uv_df(2,:,k) - qb_df(4,:)/qb_df(1,:)
        end do 

        call extract_qprime_df_face(qprime_df_face,qprime_df)

    end subroutine evaluate_bcl

    subroutine evaluate_bcl_v1(q_df, qprime_df, qb_df)

        use mod_grid, only : npoin, intma, face, nface 
        use mod_basis, only : ngl
        use mod_input, only: nlayers
        use mod_initial, only: pbprime_df

        implicit none
        
        real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(4,npoin), intent(in) :: qb_df
        
        integer :: I, k, el, er, n , iface, il, jl, kl, Iq
        real :: uv_df(2,npoin,nlayers), one_plus_eta_temp(npoin)

        call extract_velocity(uv_df, q_df, qb_df)

        do k = 1,nlayers
            q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
            q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
        end do

        call extract_velocity(uv_df, q_df, qb_df)

        do k = 1,nlayers
            qprime_df(2,:,k) = uv_df(1,:,k) - qb_df(3,:)/qb_df(1,:)
            qprime_df(3,:,k) = uv_df(2,:,k) - qb_df(4,:)/qb_df(1,:)
        end do 
        
    end subroutine evaluate_bcl_v1

    subroutine extract_velocity(uv_df, q_df, qb_df)

        use mod_grid, only : npoin, intma, face, nface 
        use mod_basis, only : ngl
        use mod_input, only: nlayers
        use mod_face, only: imapl

        implicit none
        
        real, dimension(3,npoin,nlayers), intent(in) :: q_df
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(2,npoin,nlayers), intent(out) :: uv_df
        
        real :: ubar, vbar
        integer :: I, k, el, er, n , iface, il, jl, kl, Iq

        uv_df = 0.0

        do k = 1,nlayers
            uv_df(1,:,k) = q_df(2,:,k) / q_df(1,:,k)
            uv_df(2,:,k) = q_df(3,:,k) / q_df(1,:,k)
        end do

        do I = 1, npoin
            ubar = 0.0
            vbar = 0.0
            
            do k = 1, nlayers
                ubar = ubar + uv_df(1,I,k) * q_df(1,I,k)
                vbar = vbar + uv_df(2,I,k) * q_df(1,I,k)
            end do

            if(qb_df(1,I) > 0.0) then
                
                ubar = ubar / qb_df(1,I)
                vbar = vbar / qb_df(1,I)

                do k = 1, nlayers

                    uv_df(1,I,k) = uv_df(1,I,k) - ubar + qb_df(3,I)/qb_df(1,I)
                    uv_df(2,I,k) = uv_df(2,I,k) - vbar + qb_df(4,I)/qb_df(1,I)
                end do

            else
                uv_df(:,I,:) = 0.0
            end if
        end do
        
    end subroutine extract_velocity

    subroutine evaluate_mom(q,q_df)

        ! Interpolate from dofs to quadrature points 

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad
        use mod_input, only: nlayers
        use mod_initial, only: psih, indexq

        implicit none

        real, intent(inout) :: q(3,npoin_q,nlayers)
        real, intent(in) :: q_df(3,npoin,nlayers)

        integer :: k, Iq, m, n, l,I, ip
        real :: hi

        q(2:3,:,:) = 0.0

        do Iq = 1,npoin_q
            do ip = 1,npts

                I = indexq(ip,Iq)
                hi = psih(ip,Iq)
                
                q(2,Iq,:) = q(2,Iq,:) + q_df(2,I,:)*hi
                q(3,Iq,:) = q(3,Iq,:) + q_df(3,I,:)*hi
            end do
        end do

    end subroutine evaluate_mom

    subroutine interpolate_qprime(qprime,qprime_face,qprime_df_face,qprime_df)

        ! Interpolate from dofs to quadrature points 

        use mod_basis, only: nq, npts, ngl
        use mod_grid, only:  npoin, npoin_q, intma, intma_dg_quad, face
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q, normal_vector_q, normal_vector, imapl, imapr
        use mod_initial, only: psih, indexq

        implicit none

        real, dimension(3,npoin_q,nlayers), intent(out) :: qprime
        real, dimension(3,2,nq,nface,nlayers), intent(out) :: qprime_face
        real, dimension(3,2,ngl,nface,nlayers), intent(out) :: qprime_df_face
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        integer :: k, Iq, I, ip, iface, iquad, el, er, il, jl, ir, jr, kl, kr, n
        real :: hi, nx, ny, un(nlayers)

        qprime = 0.0
        qprime_face = 0.0

        do Iq = 1,npoin_q
            do ip = 1,npts

                I = indexq(ip,Iq)
                hi = psih(ip,Iq)

                qprime(1,Iq,:) = qprime(1,Iq,:) + hi*qprime_df(1,I,:)
                qprime(2,Iq,:) = qprime(2,Iq,:) + hi*qprime_df(2,I,:)
                qprime(3,Iq,:) = qprime(3,Iq,:) + hi*qprime_df(3,I,:)
            end do
        end do
    
        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                il = imapl_q(1,iquad,1,iface)
                jl = imapl_q(2,iquad,1,iface)
                kl = imapl_q(3,iquad,1,iface)
                I = intma_dg_quad(il,jl,kl,el)

                qprime_face(1:3,1,iquad,iface,:) = qprime(1:3,I,:)

                if(er > 0) then

                    ir = imapr_q(1,iquad,1,iface)
                    jr = imapr_q(2,iquad,1,iface)
                    kr = imapr_q(3,iquad,1,iface)
                    I = intma_dg_quad(ir,jr,kr,er)

                    qprime_face(1:3,2,iquad,iface,:) = qprime(1:3,I,:)

                else 

                    qprime_face(1:3,2,iquad,iface,:) = qprime_face(1:3,1,iquad,iface,:)

                    if(er == -4) then
                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)
                        un = qprime(2,I,:)*nx + qprime(3,I,:)*ny
                        qprime_face(2,2,iquad,iface,:) = qprime(2,I,:) - 2.0*un*nx
                        qprime_face(3,2,iquad,iface,:) = qprime(3,I,:) - 2.0*un*ny
                    elseif(er == -2) then 
                        qprime_face(2:3,2,iquad,iface,:) = -qprime_face(2:3,1,iquad,iface,:)
                    end if
                end if

            end do

            do n = 1, ngl

                il = imapl(1,n,1,iface)
                jl = imapl(2,n,1,iface)
                kl = imapl(3,n,1,iface)
                I = intma(il,jl,kl,el)

                qprime_df_face(1:3,1,n,iface,:) = qprime_df(1:3,I,:)

                if(er > 0) then

                    ir = imapr(1,n,1,iface)
                    jr = imapr(2,n,1,iface)
                    kr = imapr(3,n,1,iface)
                    I = intma(ir,jr,kr,er)

                    qprime_df_face(1:3,2,n,iface,:) = qprime_df(1:3,I,:)

                else 

                    qprime_df_face(1:3,2,n,iface,:) = qprime_df_face(1:3,1,n,iface,:)

                    if(er == -4) then
                        nx = normal_vector(1,n,1,iface)
                        ny = normal_vector(2,n,1,iface)
                        un = qprime_df(2,I,:)*nx + qprime_df(3,I,:)*ny
                        qprime_df_face(2,2,n,iface,:) = qprime_df(2,I,:) - 2.0*un*nx
                        qprime_df_face(3,2,n,iface,:) = qprime_df(3,I,:) - 2.0*un*ny
                    elseif(er == -2) then
                        qprime_df_face(2:3,2,n,iface,:) = -qprime_df_face(2:3,1,n,iface,:)
                    end if
                end if
            end do
        end do

    end subroutine interpolate_qprime

    subroutine extract_qprime_df_face(qprime_df_face,qprime_df)

        ! Interpolate from dofs to quadrature points

        use mod_basis, only: nq, npts, ngl
        use mod_grid, only:  npoin, npoin_q, intma, intma_dg_quad, face
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q, normal_vector_q, normal_vector, imapl, imapr
        use mod_initial, only: psih, indexq

        implicit none

        real, dimension(3,2,ngl,nface,nlayers), intent(out) :: qprime_df_face
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        integer :: k, Iq, I, ip, iface, iquad, el, er, il, jl, ir, jr, kl, kr, n
        real :: hi, nx, ny, un(nlayers)

        qprime_df_face = 0.0

        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do n = 1, ngl

                il = imapl(1,n,1,iface)
                jl = imapl(2,n,1,iface)
                kl = imapl(3,n,1,iface)
                I = intma(il,jl,kl,el)

                qprime_df_face(1:3,1,n,iface,:) = qprime_df(1:3,I,:)

                if(er > 0) then

                    ir = imapr(1,n,1,iface)
                    jr = imapr(2,n,1,iface)
                    kr = imapr(3,n,1,iface)
                    I = intma(ir,jr,kr,er)

                    qprime_df_face(1:3,2,n,iface,:) = qprime_df(1:3,I,:)

                else

                    qprime_df_face(1:3,2,n,iface,:) = qprime_df_face(1:3,1,n,iface,:)

                    if(er == -4) then
                        nx = normal_vector(1,n,1,iface)
                        ny = normal_vector(2,n,1,iface)
                        un = qprime_df(2,I,:)*nx + qprime_df(3,I,:)*ny
                        qprime_df_face(2,2,n,iface,:) = qprime_df(2,I,:) - 2.0*un*nx
                        qprime_df_face(3,2,n,iface,:) = qprime_df(3,I,:) - 2.0*un*ny
                    elseif(er == -2) then
                        qprime_df_face(2:3,2,n,iface,:) = -qprime_df_face(2:3,1,n,iface,:)
                    end if
                end if
            end do
        end do

    end subroutine extract_qprime_df_face

    subroutine extract_dprime_df_face(dprime_df_face,dprime_df)

        ! Interpolate from dofs to quadrature points
        use mod_basis, only: nq, npts, ngl
        use mod_grid, only:  npoin, npoin_q, intma, intma_dg_quad, face
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q, normal_vector_q, normal_vector, imapl, imapr
        use mod_initial, only: psih, indexq

        implicit none

        real, dimension(2,ngl,nface,nlayers), intent(out) :: dprime_df_face
        real, dimension(npoin,nlayers), intent(in) :: dprime_df

        integer :: k, Iq, I, ip, iface, iquad, el, er, il, jl, ir, jr, kl, kr, n
        real :: hi, nx, ny, un(nlayers)

        dprime_df_face = 0.0

        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do n = 1, ngl

                il = imapl(1,n,1,iface)
                jl = imapl(2,n,1,iface)
                kl = imapl(3,n,1,iface)
                I = intma(il,jl,kl,el)

                dprime_df_face(1,n,iface,:) = dprime_df(I,:)

                if(er > 0) then

                    ir = imapr(1,n,1,iface)
                    jr = imapr(2,n,1,iface)
                    kr = imapr(3,n,1,iface)
                    I = intma(ir,jr,kr,er)

                    dprime_df_face(2,n,iface,:) = dprime_df(I,:)
                else
                    dprime_df_face(2,n,iface,:) = dprime_df_face(1,n,iface,:)
                end if
            end do
        end do

    end subroutine extract_dprime_df_face

    subroutine evaluate_mom_face(q_face, q)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz,nq
        use mod_grid, only:  npoin_q, intma_dg_quad, mod_grid_get_face_nq, nface,face
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q, normal_vector_q

        implicit none

        real, dimension(3, 2, nq, nface, nlayers), intent(inout) :: q_face
        real, dimension(3, npoin_q, nlayers), intent(in) :: q
    
        integer :: k, iface, iquad, ilocl, ilocr, el, er, il, jl, ir, jr, I, kl, kr, jquad
        real :: nx, ny, un(nlayers)
    
        q_face(2:3,:,:,:,:) = 0.0
    
        ! do k = 1, nlayers
        do iface = 1, nface

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1, nq

                il = imapl_q(1,iquad,1,iface)
                jl = imapl_q(2,iquad,1,iface)
                kl = imapl_q(3,iquad,1,iface)
                I = intma_dg_quad(il,jl,kl,el)

                q_face(2:3,1,iquad,iface,:) = q(2:3,I,:)

                if(er > 0) then

                    ir = imapr_q(1,iquad,1,iface)
                    jr = imapr_q(2,iquad,1,iface)
                    kr = imapr_q(3,iquad,1,iface)
                    I = intma_dg_quad(ir,jr,kr,er)

                    q_face(2:3,2,iquad,iface,:) = q(2:3,I,:)

                else 
                    q_face(2:3,2,iquad,iface,:) = q_face(2:3,1,iquad,iface,:)
                    if(er == -4) then

                        nx = normal_vector_q(1,iquad,1,iface)
                        ny = normal_vector_q(2,iquad,1,iface)

                        un = q(2,I,:)*nx + q(3,I,:)*ny
                        q_face(2,2,iquad,iface,:) = q(2,I,:) - 2.0*un*nx
                        q_face(3,2,iquad,iface,:) = q(3,I,:) - 2.0*un*ny
                    elseif(er == -2) then 
                        q_face(2:3,2,iquad,iface,:) = -q_face(2:3,1,iquad,iface,:)
                    end if
                end if
            end do
        end do
        ! end do
    
    end subroutine evaluate_mom_face

    subroutine layer_mom_boundary_df(q)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl,nq
        use mod_grid, only:  npoin, intma, nface, face,mod_grid_get_face_nq
        use mod_initial, only: pbprime_face
        use mod_face, only: imapl, imapr, normal_vector
        use mod_input, only: nlayers, mlswe_bc_strong

        implicit none

        real, intent(inout) :: q(2,npoin,nlayers)

        integer :: iface, n, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr,k
        real :: nx, ny, upnl(nlayers)

        do iface = 1, nface                 

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            if(er == -4) then
                do n = 1, ngl

                    il=imapl(1,n,1,iface)
                    jl=imapl(2,n,1,iface)
                    kl=imapl(3,n,1,iface)
                    I=intma(il,jl,kl,el)

                    nx = normal_vector(1,n,1,iface)
                    ny = normal_vector(2,n,1,iface)

                    upnl = q(1,I,:)*nx + q(2,I,:)*ny
                    q(1,I,:) = q(1,I,:) - upnl*nx
                    q(2,I,:) = q(2,I,:) - upnl*ny

                end do

            elseif(er == -2) then
                do n = 1, ngl

                    il=imapl(1,n,1,iface)
                    jl=imapl(2,n,1,iface)
                    kl=imapl(3,n,1,iface)
                    I=intma(il,jl,kl,el)

                    q(1,I,:) = 0.0
                    q(2,I,:) = 0.0

                end do
            end if
        end do
      
    end subroutine layer_mom_boundary_df 

    subroutine shear_stress_system(uv,q)

        use mod_constants, only : gravity
        use mod_initial, only : alpha_mlswe, coriolis_quad
        use mod_grid, only : nface, npoin, npoin_q
        use mod_input, only : dt, ad_mlswe, nlayers, max_shear_dz
  
        implicit none
        
        real, dimension(3,npoin_q,nlayers), intent(in) :: q
        real, dimension(2,npoin_q,nlayers), intent(out) :: uv
        real, dimension(npoin_q) :: coeff
        real, dimension(nlayers) :: a,b,c
        real, dimension(nlayers,2) :: r
        real :: mult
        integer :: Iq,k

        do Iq = 1,npoin_q
            coeff(Iq) = gravity*dt*max(sqrt(0.5*coriolis_quad(Iq)*ad_mlswe)/alpha_mlswe(1), &
                                  ad_mlswe/(alpha_mlswe(1) * max_shear_dz))
        end do
        
        do Iq = 1, npoin_q
            do k = 1, nlayers
                a(k) = -coeff(Iq)
                b(k) = q(1,Iq,k) + 2.0*coeff(Iq)
                c(k) = -coeff(Iq)
                r(k,1) = q(2,Iq,k)/q(1,Iq,k)
                r(k,2) = q(3,Iq,k)/q(1,Iq,k)
            end do
            
            b(1) = q(1,Iq,1) + coeff(Iq)
            b(nlayers) = q(1,Iq,nlayers) + coeff(Iq)
            a(1) = 0.0
            c(nlayers) = 0.0
            
            do k = 2, nlayers
                mult = a(k) / b(k-1)
                b(k) = b(k) - mult*c(k-1)
                r(k,1) = r(k,1) - mult*r(k-1,1)
                r(k,2) = r(k,2) - mult*r(k-1,2)
            end do
            
            r(nlayers,1) = r(nlayers,1) / b(nlayers)
            r(nlayers,2) = r(nlayers,2) / b(nlayers)
            uv(1,Iq,nlayers) = r(nlayers,1)
            uv(2,Iq,nlayers) = r(nlayers,2)
            
            do k = nlayers-1,1,-1
                r(k,1) = (r(k,1) - c(k)*r(k+1,1)) / b(k)
                r(k,2) = (r(k,2) - c(k)*r(k+1,2)) / b(k)
                uv(1,Iq,k) = r(k,1)
                uv(2,Iq,k) = r(k,2)
            end do
        end do
        
    end subroutine shear_stress_system

end module mod_layer_terms
