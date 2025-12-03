! ===============================================================================================
! This module contains the routines for the baroclinic flux terms
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: March 27, 2023
! ===============================================================================================

module mod_layer_terms
    
    use mod_grid, only: npoin_q, nface, intma_dg_quad, face, npoin
    use mod_basis, only: nqx, nqy, nqz, nq
    use mod_input, only: nlayers
        
    implicit none

    public :: layer_mom_boundary_df,  &
            velocity_df, &
            evaluate_consistency_face, &
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
        use mod_initial, only: pbprime_df
        use mod_basis, only: ngl, psiq

        implicit none
        real, intent(out) :: mass_deficit_mass_face(2,2,nq,nface,nlayers)
        real, intent(in) :: dprime_df(npoin,nlayers)
    
        ! Local variables
        integer :: iface, iquad, k, il, jl, kl, el, er, ir, jr, kr, I, n
        real :: qprime_l, qprime_r, weights_face_l, weights_face_r, hi
        real :: pbprime_l, pbprime_r

        mass_deficit_mass_face = 0.0

        do k = 1,nlayers

            do iface = 1, nface

                !Store Left Side Variables

                el = face(7,iface)
                er = face(8,iface)

                do iquad = 1, nq
                    
                    qprime_l = 0.0; qprime_r = 0.0
                    pbprime_l = 0.0; pbprime_r = 0.0

                    do n = 1,ngl
                        hi = psiq(n,iquad)

                        il = imapl(1,n,1,iface)
                        jl = imapl(2,n,1,iface)
                        kl = imapl(3,n,1,iface)

                        I = intma(il,jl,kl,el)
                        qprime_l = qprime_l + hi*dprime_df(I,k)
                        pbprime_l = pbprime_l + hi*pbprime_df(I)
                    enddo

                    if(er > 0) then
                        do n = 1,ngl
                            hi = psiq(n,iquad)
                            ir = imapr(1,n,1,iface)
                            jr = imapr(2,n,1,iface)
                            kr = imapr(3,n,1,iface)
                            I = intma(ir,jr,kr,er)
                            qprime_r = qprime_r + hi*dprime_df(I,k)
                            pbprime_r = pbprime_r + hi*pbprime_df(I)
                        enddo
                    else
                        qprime_r = qprime_l
                        pbprime_r = pbprime_l
                    endif
  
                    weights_face_l = qprime_l / pbprime_l
                    weights_face_r = qprime_r / pbprime_r

                    mass_deficit_mass_face(1,1,iquad,iface,k) = weights_face_l* &
                                                        (btp_mass_flux_face_ave(1,iquad,iface) &
                                                        - sum_layer_mass_flux_face(1,iquad,iface))
                    mass_deficit_mass_face(2,1,iquad,iface,k) = weights_face_l* &
                                                        (btp_mass_flux_face_ave(2,iquad,iface) &
                                                        - sum_layer_mass_flux_face(2,iquad,iface))

                    mass_deficit_mass_face(1,2,iquad,iface,k) = weights_face_r* &
                                                        (btp_mass_flux_face_ave(1,iquad,iface) &
                                                        - sum_layer_mass_flux_face(1,iquad,iface))
                    mass_deficit_mass_face(2,2,iquad,iface,k) = weights_face_r* &
                                                        (btp_mass_flux_face_ave(2,iquad,iface) &
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

    subroutine extract_qprime_df_face(qprime_df_face,qprime_df, q_df, qb_df)

        ! Interpolate from dofs to quadrature points

        use mod_basis, only: nq, npts, ngl
        use mod_grid, only:  npoin, npoin_q, intma, intma_dg_quad, face
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q, normal_vector_q, normal_vector, imapl, imapr
        use mod_initial, only: psih, indexq, pbprime_df

        implicit none

        real, dimension(3,2,ngl,nface,nlayers), intent(out) :: qprime_df_face
        real, dimension(3,npoin,nlayers), intent(out) :: qprime_df
        real, dimension(3,npoin,nlayers), intent(in) :: q_df
        real, dimension(4,npoin), intent(in) :: qb_df

        integer :: k, Iq, I, ip, iface, iquad, el, er, il, jl, ir, jr, kl, kr, n
        real :: hi, nx, ny, un(nlayers)
        real, dimension(npoin) :: one_plus_eta_temp
        real :: uv_df(2,npoin,nlayers)

        qprime_df = 0.0
        qprime_df_face = 0.0

        call extract_velocity(uv_df, q_df, qb_df)

        ! Prime variables at the dofs (nodal points) and quadrature points
        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)

        do k = 1,nlayers

            qprime_df(1,:,k) = q_df(1,:,k) / one_plus_eta_temp(:)
            qprime_df(2,:,k) = uv_df(1,:,k) - qb_df(3,:)/qb_df(1,:)
            qprime_df(3,:,k) = uv_df(2,:,k) - qb_df(4,:)/qb_df(1,:)

        end do

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

    subroutine layer_mom_boundary_df(q)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl,nq
        use mod_grid, only:  npoin, intma, nface, face,mod_grid_get_face_nq
        use mod_face, only: imapl, imapr, normal_vector
        use mod_input, only: nlayers

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

end module mod_layer_terms
