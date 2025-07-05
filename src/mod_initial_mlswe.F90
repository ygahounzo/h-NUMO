!----------------------------------------------------------------------!
!>@brief This module contains Initial Conditions settings
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: March 27, 2023
!----------------------------------------------------------------------!
module mod_initial_mlswe

    use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz
    use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad, mod_grid_get_face_nq
    use mod_basis, only: nq, psiq, psiqx, psiqy, psiqz, dpsiq, dpsiqx, dpsiqy, dpsiqz
    use mod_metrics, only: ksiq_x,ksiq_y,ksiq_z, etaq_x,etaq_y,etaq_z, zetaq_x,zetaq_y,zetaq_z

    public :: &
        bot_topo_derivatives, &
        compute_gradient_quad, &
        interpolate_from_dof_to_quad_uv_init, &
        interpolate_pbprime_init, wind_stress_coriolis, compute_reference_edge_variables, &
        map_deriv,interpolate_layer_from_quad_to_node_1d, ssprk_coefficients, &
        compute_reference_edge_variables_df !, Tensor_product

    private


    contains
    
    !> Compute the gradient of the bottom topography at quadrature points
    subroutine bot_topo_derivatives(zbot,zbot_face,zbot_df)

        use mod_grid, only : npoin, npoin_q, nface, intma_dg_quad, face, &
                            mod_grid_get_face_nq, coord, nelem, intma
        use mod_basis, only : nq, nglx, ngly, nqx, nqy
        use mod_input, only : nlayers, test_case, xdims, ydims
        use mod_face, only : imapl_q, imapr_q


        implicit none
        
        ! Declare output arguments
        real, dimension(npoin), intent(in) :: zbot_df
        real, dimension(npoin_q), intent(out) :: zbot
        real, dimension(2,nq,nface), intent(out) :: zbot_face
        
        ! Declare local variables
        integer :: il, jl, ir, jr,kl,kr, el, er, ilocl, ilocr, l, I1,I2, iface, &
                    iquad, jquad, nq_i, nq_j, plane_ij
        integer :: e,I, m,n,Iq
        real :: depth, x, y, r, xm, yl, hi

        
        ! Initialize output arrays

        do e = 1, nelem
            do jquad = 1, nqy
                do iquad = 1, nqx
                
                    Iq = intma_dg_quad(iquad, jquad, 1, e)

                    do m = 1, ngly
                        do n = 1, nglx
                            
                            I = intma(n, m, 1, e)
                            
                            hi = psiqx(n, iquad) * psiqy(m, jquad)
                            
                            zbot(Iq) = zbot(Iq) + zbot_df(I) * hi
                        end do 
                    end do 
                end do 
            end do 
        end do

        zbot_face = 0.0

        !open(10,file='zbot.txt',status='unknown')
        
        do iface=1,nface

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)
                
            do jquad = 1,1

                do iquad = 1, nq

                    il = imapl_q(1,iquad,jquad,iface)
                    jl = imapl_q(2,iquad,jquad,iface)
                    kl = imapl_q(3,iquad,jquad,iface)

                    I1 = intma_dg_quad(il,jl,kl,el)

                    zbot_face(1,iquad,iface) = zbot(I1)

                    if(er > 0) then

                        ir = imapr_q(1,iquad,jquad,iface)
                        jr = imapr_q(2,iquad,jquad,iface)
                        kr = imapr_q(3,iquad,jquad,iface)

                        I2 = intma_dg_quad(ir,jr,kr,er)

                        zbot_face(2,iquad,iface) = zbot(I2)

                    else

                        zbot_face(2,iquad,iface) = zbot_face(1,iquad,iface)

                    end if

                end do
            end do
        end do
        
    end subroutine bot_topo_derivatives

    !> Compute the gradient of the bottom topography at quadrature points
    subroutine interpolate_from_dof_to_quad_uv_init(q, q_df)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad
        use mod_basis, only: nq, psiq, psiqx, psiqy, psiqz
        use mod_input, only: nlayers

        implicit none
        real, dimension(5,npoin_q,nlayers), intent(inout) :: q
        real, dimension(5,npoin,nlayers), intent(in) :: q_df
        integer :: k, e, iquad, jquad, kquad, l, m, n, I, Iq
        real :: hi

        q(2:3,:,:) = 0.0
        
        do k = 1, nlayers
            do e = 1, nelem
                do kquad = 1, nqz
                    do jquad = 1, nqy
                        do iquad = 1, nqx
                        
                            Iq = intma_dg_quad(iquad, jquad, kquad, e)
                            
                            do l = 1, nglz
                                do m = 1, ngly
                                    do n = 1, nglx
                                        
                                        I = intma(n, m, l, e)
                                        
                                        hi = psiqx(n, iquad) * psiqy(m, jquad)!* psiqz(l, kquad)
                                        
                                        q(2,Iq,k) = q(2,Iq,k) + q_df(2,I,k) * hi
                                        q(3,Iq,k) = q(3,Iq,k) + q_df(3,I,k) * hi
                                        
                                    end do
                                end do
                            end do
                        end do
                        
                    end do
                end do
            end do
        end do

    end subroutine interpolate_from_dof_to_quad_uv_init

    !> Interpolate pressure perturbation from DOF to quadrature points
    subroutine interpolate_pbprime_init(pbprime,pbprime_face,pbprime_df_face,pbprime_df)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad, nface, face
        use mod_grid, only:  mod_grid_get_face_nq, mod_grid_get_face_ngl
        use mod_basis, only: nq, ngl
        use mod_input, only: nlayers
        use mod_face, only: imapl_q, imapr_q, imapl, imapr

        implicit none
        real, dimension(npoin_q), intent(out) :: pbprime
        real, dimension(2,nq,nface), intent(out) :: pbprime_face
        real, dimension(2,ngl,nface), intent(out) :: pbprime_df_face
        real, dimension(npoin), intent(in) :: pbprime_df
        integer :: k, e, iquad, jquad, kquad, l, m, n, I, Iq, iface, I1
        integer :: I2, nq_i, nq_j, plane_ij, ilocl, ilocr, el, er, il, jl, kl, ir, jr, kr
        integer :: ngl_i, ngl_j, ii, jj 
        real :: hi

        pbprime = 0.0
        pbprime_face = 0.0
        pbprime_df_face = 0.0
        
        do e = 1, nelem
            do kquad = 1, nqz
                do jquad = 1, nqy
                    do iquad = 1, nqx
                    
                        Iq = intma_dg_quad(iquad, jquad, kquad, e)
                        
                        do l = 1, nglz
                            do m = 1, ngly
                                do n = 1, nglx
                                    
                                    I = intma(n, m, l, e)
                                    
                                    hi = psiqx(n, iquad) * psiqy(m, jquad)!* psiqz(l, kquad)
                                    !print*, "hi = ", hi
                                    pbprime(Iq) = pbprime(Iq) + pbprime_df(I) * hi
                                    
                                end do
                            end do
                        end do
                    end do
                    
                end do
            end do
        end do

        do iface = 1, nface

            !Store Left Side Variables
            ilocl = face(5,iface)
            ilocr = face(6,iface)
            el = face(7,iface)
            er = face(8,iface)

            call mod_grid_get_face_nq(ilocl, nq_i, nq_j, plane_ij)
            call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

            do jquad = 1,nq_j
                do iquad = 1, nq_i

                    il = imapl_q(1,iquad,jquad,iface)
                    jl = imapl_q(2,iquad,jquad,iface)
                    kl = imapl_q(3,iquad,jquad,iface)
                    I1 = intma_dg_quad(il,jl,kl,el)

                    pbprime_face(1,iquad,iface) = pbprime(I1)
                    if(er > 0) then

                        ir = imapr_q(1,iquad,jquad,iface)
                        jr = imapr_q(2,iquad,jquad,iface)
                        kr = imapr_q(3,iquad,jquad,iface)
                        I2 = intma_dg_quad(ir,jr,kr,er)

                        pbprime_face(2,iquad,iface) = pbprime(I2)
                    else
                        pbprime_face(2,iquad,iface) = pbprime_face(1,iquad,iface)
                    end if
                end do
            end do

            do jj = 1,ngl_j
                do ii = 1, ngl_i

                    il = imapl(1,ii,jj,iface)
                    jl = imapl(2,ii,jj,iface)
                    kl = imapl(3,ii,jj,iface)
                    I1 = intma(il,jl,kl,el)

                    pbprime_df_face(1,ii,iface) = pbprime_df(I1)
                    if(er > 0) then

                        ir = imapr(1,ii,jj,iface)
                        jr = imapr(2,ii,jj,iface)
                        kr = imapr(3,ii,jj,iface)
                        I2 = intma(ir,jr,kr,er)

                        pbprime_df_face(2,ii,iface) = pbprime_df(I2)
                    else
                        pbprime_df_face(2,ii,iface) = pbprime_df_face(1,ii,iface)
                    endif
                enddo
            enddo
        end do

    end subroutine interpolate_pbprime_init

    !> Wind Stress and Coriolis Force
    subroutine wind_stress_coriolis(tau_wind,coriolis_df,coriolis_quad, fdt_bcl, fdt2_bcl, &
            a_bcl, b_bcl,tau_wind_df)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, npts
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad, coord
        use mod_constants, only: gravity, pi, tol, omega, earth_radius
        use mod_input, only: gravity_in, &
            nelx, nelz, &
            xdims, ydims, nlayers, dt, dt_btp, test_case, f0, beta
    
        implicit none

        real, dimension(2,npoin_q), intent(out) :: tau_wind
        real, dimension(npoin), intent(out) :: coriolis_df
        real, dimension(npoin_q), intent(out) :: coriolis_quad
        real, dimension(npoin), intent(out) :: fdt_bcl, fdt2_bcl, a_bcl, b_bcl

        real, dimension(2,npoin), intent(in) :: tau_wind_df

        integer :: k, e, iquad, jquad, kquad, l, m, n, I, Iq, ip
        real :: ym, Ly, y, hi, tau0, lat, omega1, sig, rho_air, w

        tau_wind = 0.0
        coriolis_df = 0.0
        coriolis_quad = 0.0

        gravity = 9.806
        
        Ly = ydims(2)
        ym = 0.5*Ly
        
        do I = 1, npoin
            y = coord(2,I)

            coriolis_df(I) = f0 + beta*(y - ym)
        end do


        do e = 1, nelem
            do kquad = 1, nqz
                do jquad = 1, nqy
                    do iquad = 1, nqx
                    
                        Iq = intma_dg_quad(iquad, jquad, kquad, e)
                        
                        do l = 1, nglz
                            do m = 1, ngly
                                do n = 1, nglx
                                    
                                    I = intma(n, m, l, e)
                                    
                                    hi = psiqx(n, iquad) * psiqy(m, jquad)
                                    
                                    coriolis_quad(Iq) = coriolis_quad(Iq) + coriolis_df(I) * hi
                                    
                                    tau_wind(1,Iq) = tau_wind(1,Iq) + tau_wind_df(1,I) * hi
                                    tau_wind(2,Iq) = tau_wind(2,Iq) + tau_wind_df(2,I) * hi
                                    
                                end do
                            end do
                        end do
                    end do
                    
                end do
            end do
        end do
        
        fdt_bcl = dt * coriolis_df
        fdt2_bcl = 0.5 * fdt_bcl
        a_bcl = 1.0 / (1.0 + fdt2_bcl**2)
        b_bcl = fdt2_bcl / (1.0 + fdt2_bcl**2)

    end subroutine wind_stress_coriolis

    !> Compute the wave speeds for linearized shallow equation at quadrature points
    subroutine compute_reference_edge_variables(coeff_pbpert_L,coeff_pbpert_R,coeff_pbub_LR, &
            coeff_mass_pbub_L, coeff_mass_pbub_R,coeff_mass_pbpert_LR, pbprime_face,alpha)

            use mod_input, only: nlayers
            use mod_basis, only: nq
            use mod_grid, only: nface
        
            implicit none

            real, dimension(2,nq,nface), intent(in) :: pbprime_face
            real, dimension(nlayers), intent(in) :: alpha
            real, dimension(nq,nface), intent(out) :: coeff_pbpert_L,coeff_pbpert_R,coeff_pbub_LR, &
                coeff_mass_pbub_L,coeff_mass_pbub_R,coeff_mass_pbpert_LR

            integer :: iface, iquad
            real :: c_minus, c_plus
        
            coeff_pbpert_L = 0.0
            coeff_pbpert_R = 0.0
            coeff_pbub_LR = 0.0
            coeff_mass_pbub_L = 0.0
            coeff_mass_pbub_R = 0.0
            coeff_mass_pbpert_LR = 0.0
        
            do iface = 1,nface
                do iquad = 1,nq
        
                    c_minus = sqrt(alpha(nlayers) * pbprime_face(2,iquad,iface))
                    c_plus  = sqrt(alpha(nlayers) * pbprime_face(1,iquad,iface))
                    if ((c_minus > 0.0) .or. (c_plus > 0.0)) then
                        coeff_pbpert_L(iquad,iface) = c_minus / (c_minus + c_plus)
                        coeff_pbpert_R(iquad,iface) = c_plus / (c_minus + c_plus)
                        coeff_pbub_LR(iquad,iface)  = 1.0 / (c_minus + c_plus)
                    end if
        
                    c_minus = sqrt(alpha(nlayers) * pbprime_face(2,iquad,iface))
                    c_plus  = sqrt(alpha(nlayers) * pbprime_face(1,iquad,iface))
                    if ((c_minus > 0.0) .or. (c_plus > 0.0)) then
                        coeff_mass_pbub_L(iquad,iface) = c_plus / (c_minus + c_plus)
                        coeff_mass_pbub_R(iquad,iface) = c_minus / (c_minus + c_plus)
                        coeff_mass_pbpert_LR(iquad,iface) = c_minus * c_plus / (c_minus + c_plus)
                    end if

                end do
            end do
        
    end subroutine compute_reference_edge_variables

    !> Compute the wave speeds for linearized shallow equation at dofs
    subroutine compute_reference_edge_variables_df(coeff_pbpert_L,coeff_pbpert_R,coeff_pbub_LR, &
            coeff_mass_pbub_L, coeff_mass_pbub_R,coeff_mass_pbpert_LR, pbprime_df_face,alpha)

        use mod_input, only: nlayers
        use mod_basis, only: ngl
        use mod_grid, only: nface

        implicit none

        real, dimension(2,ngl,nface), intent(in) :: pbprime_df_face
        real, dimension(nlayers), intent(in) :: alpha
        real, dimension(ngl,nface), intent(out) :: coeff_pbpert_L,coeff_pbpert_R,coeff_pbub_LR, &
            coeff_mass_pbub_L,coeff_mass_pbub_R,coeff_mass_pbpert_LR

        integer :: iface, n, il, jl, kl, ir, jr, kr, el, er, I1, I2
        real :: c_minus, c_plus, pl, pr

        coeff_pbpert_L = 0.0
        coeff_pbpert_R = 0.0
        coeff_pbub_LR = 0.0
        coeff_mass_pbub_L = 0.0
        coeff_mass_pbub_R = 0.0
        coeff_mass_pbpert_LR = 0.0

        do iface = 1,nface
            do n = 1,ngl

                c_minus = sqrt(alpha(nlayers) * pbprime_df_face(2,n,iface))
                c_plus  = sqrt(alpha(nlayers) * pbprime_df_face(1,n,iface))
                if ((c_minus > 0.0) .or. (c_plus > 0.0)) then
                    coeff_pbpert_L(n,iface) = c_plus / (c_minus + c_plus)
                    coeff_pbpert_R(n,iface) = c_minus / (c_minus + c_plus)
                    coeff_pbub_LR(n,iface)  = 1.0 / (c_minus + c_plus)
                end if

                c_minus = sqrt(alpha(nlayers) * pbprime_df_face(2,n,iface))
                c_plus  = sqrt(alpha(nlayers) * pbprime_df_face(1,n,iface))
                if ((c_minus > 0.0) .or. (c_plus > 0.0)) then
                    coeff_mass_pbub_L(n,iface) = c_plus / (c_minus + c_plus)
                    coeff_mass_pbub_R(n,iface) = c_minus / (c_minus + c_plus)
                    coeff_mass_pbpert_LR(n,iface) = (c_minus * c_plus) / (c_minus + c_plus)
                end if

            end do
        end do

    end subroutine compute_reference_edge_variables_df


!    subroutine Tensor_product(wjac,psih,dpsidx,dpsidy,indexq, wjac_df,psih_df,dpsidx_df, &
!               dpsidy_df,index_df)
!
!        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
!
!        use mod_basis, only: nglx, ngly, nglz, npts, dpsiqx, dpsiqy, dpsiqz, nqx, nqy, nqz, &
!            psiqx, psiqy, psiqz
!
!        use mod_basis, only: dpsix, dpsiy, dpsiz, psix, psiy, psiz
!
!        use mod_metrics, only: &
!            ksiq_x, ksiq_y, ksiq_z, &
!            etaq_x, etaq_y, etaq_z, &
!            zetaq_x, zetaq_y, zetaq_z, &
!            jacq, &
!            ksi_x, ksi_y, ksi_z, &
!            eta_x, eta_y, eta_z, &
!            jac
!
!        use mod_constants, only: gravity
!
!        !use mod_initial, only: grad_zbot_quad, tau_wind
!
!        implicit none 
!
!        real, dimension(npts,npoin_q), intent(out) :: psih, dpsidx,dpsidy
!        integer, dimension(npts,npoin_q) :: indexq
!        real, dimension(npoin_q), intent(out) :: wjac
!
!        real, dimension(npts,npoin), intent(out) :: psih_df, dpsidx_df,dpsidy_df
!        integer, dimension(npts,npoin) :: index_df
!        real, dimension(npoin), intent(out) :: wjac_df
!
!        integer :: e, jquad, iquad, Iq, ip, m, n, I
!        real :: h_e, h_n
!        real :: e_x, e_y, n_x, n_y
!      
!        wjac = 0.0
!        psih = 0.0
!        dpsidx = 0.0
!        dpsidy = 0.0
!        indexq = 0
!
!        wjac_df = 0.0
!        psih_df = 0.0
!        dpsidx_df = 0.0
!        dpsidy_df = 0.0
!        index_df = 0
!
!        do e = 1,nelem
!
!            do jquad = 1,nqy
!                do iquad = 1,nqx
!                    
!                    Iq = intma_dg_quad(iquad,jquad,1,e)
!
!                    wjac(Iq) = jacq(iquad,jquad,1,e)
!
!                    e_x = ksiq_x(iquad,jquad,1,e); e_y = ksiq_y(iquad,jquad,1,e); 
!                    n_x = etaq_x(iquad,jquad,1,e); n_y = etaq_y(iquad,jquad,1,e);
!
!                    ip = 0
!                    
!                    do m = 1, ngly 
!                        do n = 1, nglx 
!
!                            I = intma(n,m,1,e)
!                            ip = ip + 1
!                
!                            indexq(ip,Iq) = I
!                            psih(ip,Iq) = psiqx(n, iquad) * psiqy(m, jquad)
!                
!                            ! Xi derivatives
!                            h_e = dpsiqx(n, iquad) * psiqy(m, jquad)
!
!                            ! Eta derivatives
!                            h_n = psiqx(n, iquad) * dpsiqy(m, jquad)
!                
!                            ! Pressure terms
!                            dpsidx(ip,Iq) = h_e * e_x + h_n * n_x
!                            dpsidy(ip,Iq) = h_e * e_y + h_n * n_y
!
!                        end do !n
!                    end do !m
!
!                end do !iquad
!            end do !jquad
!
!            do jquad = 1,ngly
!                do iquad = 1,nglx
!                    
!                    Iq = intma(iquad,jquad,1,e)
!
!                    wjac_df(Iq) = jac(iquad,jquad,1,e)
!
!                    e_x = ksi_x(iquad,jquad,1,e); e_y = ksi_y(iquad,jquad,1,e); 
!                    n_x = eta_x(iquad,jquad,1,e); n_y = eta_y(iquad,jquad,1,e);
!
!                    ip = 0
!                    
!                    do m = 1, ngly 
!                        do n = 1, nglx 
!
!                            I = intma(n,m,1,e)
!                            ip = ip + 1
!                
!                            index_df(ip,Iq) = I
!                            psih_df(ip,Iq) = psix(n, iquad) * psiy(m, jquad)
!                
!                            ! Xi derivatives
!                            h_e = dpsix(n, iquad) * psiy(m, jquad)
!
!                            ! Eta derivatives
!                            h_n = psix(n, iquad) * dpsiy(m, jquad)
!                
!                            ! Pressure terms
!                            dpsidx_df(ip,Iq) = h_e * e_x + h_n * n_x
!                            dpsidy_df(ip,Iq) = h_e * e_y + h_n * n_y
!
!                        end do !n
!                    end do !m
!
!                end do !iquad
!            end do !jquad
!
!        end do !e
!
!    end subroutine Tensor_product

    subroutine ssprk_coefficients(ssprk_a,ssprk_beta)

        use mod_input, only: kstages, ti_method_btp

        real, dimension(kstages,3), intent(out) :: ssprk_a
        real, dimension(kstages), intent(out) :: ssprk_beta

        if(ti_method_btp == 'lsrk') then

            if(kstages == 5) then 

                ssprk_a(1,1) = 0.0
                ssprk_a(2,1) = -567301805773.0 / 1357537059087.0
                ssprk_a(3,1) = -2404267990393.0 / 2016746695238.0
                ssprk_a(4,1) = -3550918686646.0 / 2091501179385.0
                ssprk_a(5,1) = -1275806237668.0 / 842570457699.0

                ssprk_beta(1) = 1432997174477.0 / 9575080441755.0
                ssprk_beta(2) = 5161836677717.0 / 13612068292357.0
                ssprk_beta(3) = 1720146321549.0 / 2090206949498.0
                ssprk_beta(4) = 3134564353537.0 / 4481467310338.0
                ssprk_beta(5) = 2277821191437.0 / 14882151754819.0

                !ssprk_a(1,1) = 0.0
                !ssprk_a(2,1) = -4.344339134485095 !-567301805773.0 / 1357537059087.0
                !ssprk_a(3,1) = 0.0                !-2404267990393.0 / 2016746695238.0
                !ssprk_a(4,1) = 3.770024161386381  !-3550918686646.0 / 2091501179385.0
                !ssprk_a(5,1) = -0.046347284573284 !-1275806237668.0 / 842570457699.0

                !ssprk_beta(1) = 0.713497331193829 !1432997174477.0 / 9575080441755.0
                !ssprk_beta(2) = 0.133505249805329 !5161836677717.0 / 13612068292357.0
                !ssprk_beta(3) = 0.713497331193829 !1720146321549.0 / 2090206949498.0
                !ssprk_beta(4) = 0.149579395628565 !3134564353537.0 / 4481467310338.0
                !ssprk_beta(5) = 0.384471116121269 !2277821191437.0 / 14882151754819.0

            elseif(kstages == 14) then

                ssprk_a(1,1) = 0.0
                ssprk_a(2,1) = -0.7188012108672410
                ssprk_a(3,1) = -0.7785331173421570
                ssprk_a(4,1) = -0.0053282796654044
                ssprk_a(5,1) = -0.8552979934029281
                ssprk_a(6,1) = -3.9564138245774565
                ssprk_a(7,1) = -1.5780575380587385
                ssprk_a(8,1) = -2.0837094552574054
                ssprk_a(9,1) = -0.7483334182761610
                ssprk_a(10,1) = -0.7032861106563359
                ssprk_a(11,1) =  0.0013917096117681
                ssprk_a(12,1) = -0.0932075369637460
                ssprk_a(13,1) = -0.9514200470875948
                ssprk_a(14,1) = -7.1151571693922548

                ssprk_beta(1) = 0.0367762454319673d0;
                ssprk_beta(2) = 0.3136296607553959d0;
                ssprk_beta(3) = 0.1531848691869027d0;
                ssprk_beta(4) = 0.0030097086818182d0;
                ssprk_beta(5) = 0.3326293790646110d0;
                ssprk_beta(6) = 0.2440251405350864d0;
                ssprk_beta(7) = 0.3718879239592277d0;
                ssprk_beta(8) = 0.6204126221582444d0;
                ssprk_beta(9) = 0.1524043173028741d0;
                ssprk_beta(10) = 0.0760894927419266d0;
                ssprk_beta(11) = 0.0077604214040978d0;
                ssprk_beta(12) = 0.0024647284755382d0;
                ssprk_beta(13) = 0.0780348340049386d0;
                ssprk_beta(14) = 5.5059777270269628d0;
            end if 

        else 

            select case (kstages)
            case (1)               !RK1
                ssprk_a(1,1)=1.0;       ssprk_a(1,2)=0.0;  ssprk_a(1,3)=0.0; ssprk_beta(1)=1.0
            case (2)               !RK2
                ssprk_a(1,1)=1.0;       ssprk_a(1,2)=0.0;  ssprk_a(1,3)=0.0; ssprk_beta(1)=1.0
                ssprk_a(2,1)=1.0/2.0; ssprk_a(2,2)=1.0/2.0; ssprk_a(2,3)=0.0; ssprk_beta(2)=1.0/2.0
            case (3)              !RK3; SSP(3,3)
                ssprk_a(1,1)=1.0;       ssprk_a(1,2)=0.0;   ssprk_a(1,3)=0.0; ssprk_beta(1)=1.0
                ssprk_a(2,1)=3.0/4.0; ssprk_a(2,2)=1.0/4.0; ssprk_a(2,3)=0.0; ssprk_beta(2)=1.0/4.0
                ssprk_a(3,1)=1.0/3.0; ssprk_a(3,2)=2.0/3.0; ssprk_a(3,3)=0.0; ssprk_beta(3)=2.0/3.0
            case (4)              !RK3; SSP(4,3)
                ssprk_a(1,1)=1.0;       ssprk_a(1,2)=0.0;  ssprk_a(1,3)=0.0; ssprk_beta(1)=1.0/2.0
                ssprk_a(2,1)=0.0;       ssprk_a(2,2)=1.0;  ssprk_a(2,3)=0.0; ssprk_beta(2)=1.0/2.0
                ssprk_a(3,1)=2.0/3.0; ssprk_a(3,2)=1.0/3.0; ssprk_a(3,3)=0.0; ssprk_beta(3)=1.0/6.0
                ssprk_a(4,1)=0.0;       ssprk_a(4,2)=1.0;  ssprk_a(4,3)=0.0; ssprk_beta(4)=1.0/2.0
            case (5)              !RK3; SSP(5,3)
                ssprk_a(1,1)=1.0;       ssprk_a(1,2)=0.0; ssprk_a(1,3)=0.0
                ssprk_beta(1)=0.377268915331368
                ssprk_a(2,1)=0.0;       ssprk_a(2,2)=1.0; ssprk_a(2,3)=0.0
                ssprk_beta(2)=0.377268915331368
                ssprk_a(3,1)=0.355909775063326; ssprk_a(3,2)=0.644090224936674; ssprk_a(3,3)=0.0
                ssprk_beta(3)=0.242995220537396
                ssprk_a(4,1)=0.367933791638137; ssprk_a(4,2)=0.632066208361863; ssprk_a(4,3)=0.0
                ssprk_beta(4)=0.238458932846290
                ssprk_a(5,1)=0.0; ssprk_a(5,2)=0.762406163401431; ssprk_a(5,3)=0.237593836598569
                ssprk_beta(5)=0.287632146308408
            end select
        endif 

    end subroutine ssprk_coefficients

end module mod_initial_mlswe
