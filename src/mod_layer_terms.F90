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

        use mod_grid, only : npoin, intma, face, nface, nelem
        use mod_basis, only : ngl, npts
        use mod_input, only: nlayers
        use mod_face, only: imapl
        use mod_initial, only: alpha_mlswe
        use mod_constants, only: gravity

        implicit none
        
        real, dimension(3,npoin,nlayers), intent(inout) :: q_df
        real, dimension(4,npoin), intent(in) :: qb_df
        
        real, dimension(nlayers) :: a, b, c, u_ave, v_ave
        real, dimension(nlayers,2) :: r
        real, dimension(nlayers) :: weight
        real :: eps = 1.0e-3
        real :: ubar, vbar, mult, dpbar
        integer :: I, k, el, er, n , m, iface, il, jl, kl, Iq, e
        real :: uv_df(2,npoin,nlayers), pmax, pavg, uavg, vavg, pmin
        real :: dp_cutoff1, dp_cutoff2

        dp_cutoff1 = 0.10e4
        dp_cutoff2 = 0.20e4

        uv_df = 0.0

        ! do k = 1,nlayers
        !     uv_df(1,:,k) = q_df(2,:,k) / q_df(1,:,k)
        !     uv_df(2,:,k) = q_df(3,:,k) / q_df(1,:,k)
        ! end do

        do e = 1, nelem
            do k = 1, nlayers

                pmax = -1.0e10
                pavg = 0.0
                do m = 1, ngl
                    do n = 1, ngl
                        I = intma(n,m,1,e)
                        if(q_df(1,I,k) > pmax) pmax = q_df(1,I,k)
                        pavg = pavg + q_df(1,I,k)
                        uavg = uavg + q_df(2,I,k)
                        vavg = vavg + q_df(3,I,k)
                    end do
                end do
                pavg = pavg / npts
                uavg = uavg / npts
                vavg = vavg / npts
                weight(k) = (pmax - dp_cutoff1) / (dp_cutoff2 - dp_cutoff1)
                b(k) = 1.0
                r(k,1) = weight(k)*uavg / (pavg + eps)
                r(k,2) = weight(k)*vavg / (pavg + eps)
            end do

            if (nlayers >= 3) then
                do k = 2, nlayers-1
                    a(k) = -0.5 * (1.0 - weight(k))
                    c(k) = a(k)
                end do
            endif

            a(1) = 0.0
            c(1) = -(1.0 - weight(1))
            a(nlayers) = -(1.0 - weight(nlayers))
            c(nlayers) = 0.0

            do k = 2, nlayers
                mult = a(k) / b(k-1)
                b(k) = b(k) - mult * c(k-1)
                r(k,1) = r(k,1) - mult * r(k-1,1)
                r(k,2) = r(k,2) - mult * r(k-1,2)
            end do
            r(nlayers,1) = r(nlayers,1) / b(nlayers)
            r(nlayers,2) = r(nlayers,2) / b(nlayers)
            do k = nlayers-1, 1, -1
                r(k,1) = (r(k,1) - c(k)*r(k+1,1)) / b(k)
                r(k,2) = (r(k,2) - c(k)*r(k+1,2)) / b(k)
            end do

            do k = 1, nlayers
                u_ave(k) = r(k,1)
                v_ave(k) = r(k,2)
            end do

            do k = 1, nlayers
                pmin = 1.0e10
                do m = 1, ngl
                    do n = 1, ngl
                        I = intma(n,m,1,e)
                        if(q_df(1,I,k) < pmin) pmin = q_df(1,I,k)
                    end do
                end do
                weight(k) = (pmin - dp_cutoff1) / (dp_cutoff2 - dp_cutoff1)
                weight(k) = max(0.0, min(1.0, weight(k)))

                do m = 1, ngl
                    do n = 1, ngl
                        I = intma(n,m,1,e)
                        uv_df(1,I,k) = weight(k)*q_df(2,I,k) / (q_df(1,I,k) + eps) + (1.0 - weight(k))*u_ave(k)
                        uv_df(2,I,k) = weight(k)*q_df(3,I,k) / (q_df(1,I,k) + eps) + (1.0 - weight(k))*v_ave(k)
                    end do
                end do
            end do
        enddo

        do I = 1, npoin
            ubar = 0.0 ; vbar = 0.0
            
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
        use mod_input, only: nlayers, dry_cutoff
        use mod_face, only: imapl
        use mod_initial, only: alpha_mlswe
        use mod_constants, only: gravity

        implicit none
        
        real, dimension(3,npoin,nlayers), intent(in) :: q_df
        real, dimension(4,npoin), intent(in) :: qb_df
        real, dimension(2,npoin,nlayers), intent(out) :: uv_df
        
        real :: ubar, vbar
        integer :: I, k, el, er, n , iface, il, jl, kl, Iq

        uv_df = 0.0

        do k = 1,nlayers
            do I = 1, npoin
                ! uv_df(1,I,k) = q_df(2,I,k) / q_df(1,I,k)
                ! uv_df(2,I,k) = q_df(3,I,k) / q_df(1,I,k)

                if (q_df(1,I,k) > (gravity/alpha_mlswe(k)) * dry_cutoff) then
                    uv_df(1,I,k) = q_df(2,I,k) / q_df(1,I,k)
                    uv_df(2,I,k) = q_df(3,I,k) / q_df(1,I,k)
                end if
            end do
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

                    ! uv_df(1,I,k) = uv_df(1,I,k) - ubar + qb_df(3,I)/qb_df(1,I)
                    ! uv_df(2,I,k) = uv_df(2,I,k) - vbar + qb_df(4,I)/qb_df(1,I)

                    if (q_df(1,I,k) > (gravity/alpha_mlswe(k)) * dry_cutoff) then
                        uv_df(1,I,k) = uv_df(1,I,k) - ubar + qb_df(3,I)/qb_df(1,I)
                        uv_df(2,I,k) = uv_df(2,I,k) - vbar + qb_df(4,I)/qb_df(1,I)
                    end if
                end do

            else
                uv_df(:,I,:) = 0.0
            end if
        end do
        
    end subroutine extract_velocity

    subroutine extract_qprime_df_face(qprime_df, q_df, qb_df)

        ! Interpolate from dofs to quadrature points

        use mod_basis, only: nq, npts, ngl
        use mod_grid, only:  npoin, npoin_q, intma, intma_dg_quad, face
        use mod_input, only: nlayers, dry_cutoff
        use mod_face, only: imapl_q, imapr_q, normal_vector_q, normal_vector, imapl, imapr
        use mod_initial, only: psih, indexq, pbprime_df, alpha_mlswe
        use mod_constants, only: gravity

        implicit none

        real, dimension(3,npoin,nlayers), intent(out) :: qprime_df
        real, dimension(3,npoin,nlayers), intent(in) :: q_df
        real, dimension(4,npoin), intent(in) :: qb_df

        integer :: k, Iq, I, ip, iface, iquad, el, er, il, jl, ir, jr, kl, kr, n
        real :: hi, nx, ny, un(nlayers)
        real, dimension(npoin) :: one_plus_eta_temp
        real :: uv_df(2,npoin,nlayers)

        qprime_df = 0.0

        call extract_velocity(uv_df, q_df, qb_df)

        ! Prime variables at the dofs (nodal points) and quadrature points
        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)

        do k = 1,nlayers

            do I = 1, npoin
                qprime_df(1,I,k) = q_df(1,I,k) / one_plus_eta_temp(I)
                ! qprime_df(2,I,k) = uv_df(1,I,k) - qb_df(3,I)/qb_df(1,I)
                ! qprime_df(3,I,k) = uv_df(2,I,k) - qb_df(4,I)/qb_df(1,I)

                if (qprime_df(1,I,k) > (gravity/alpha_mlswe(k)) * dry_cutoff) then
                    qprime_df(2,I,k) = uv_df(1,I,k) - qb_df(3,I)/qb_df(1,I)
                    qprime_df(3,I,k) = uv_df(2,I,k) - qb_df(4,I)/qb_df(1,I)
                end if
            enddo

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

        real, intent(inout) :: q(3,npoin,nlayers)

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

                    upnl = q(2,I,:)*nx + q(3,I,:)*ny
                    q(2,I,:) = q(2,I,:) - upnl*nx
                    q(3,I,:) = q(3,I,:) - upnl*ny

                end do

            elseif(er == -2) then
                do n = 1, ngl

                    il=imapl(1,n,1,iface)
                    jl=imapl(2,n,1,iface)
                    kl=imapl(3,n,1,iface)
                    I=intma(il,jl,kl,el)

                    q(2,I,:) = 0.0
                    q(3,I,:) = 0.0

                end do
            end if
        end do
      
    end subroutine layer_mom_boundary_df 

    subroutine limit_layer_mass(q_df)
        
        ! This routine limits the layer mass to be positive

        use mod_input, only: nlayers
        use mod_initial, only: alpha_mlswe
        use mod_metrics, only: jac
        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, npts, wglx, wgly
        use mod_grid, only:  nelem, npoin, intma, coord
        use mod_initial, only: pbprime_df
        use mod_constants, only: gravity

        implicit none

        real, dimension(3,npoin,nlayers), intent(inout) :: q_df !, qprime_df

        integer :: i, j, k, e, ip, n, m
        real :: pmin, pavg, uavg, vavg, ds, contraction, dp_adjust, ope
        real :: u_adjust, v_adjust

        !First sweep upward from the bottom.

        do e = 1,nelem

            do k = nlayers, 2, -1
                pmin = 1.0e40
                pavg = 0.0
                uavg = 0.0
                vavg = 0.0
                do m = 1, ngly
                    do n = 1, nglx
                        ip = intma(n,m,1,e)
                        if (q_df(1,ip,k) < pmin) pmin = q_df(1,ip,k)

                        pavg = pavg + q_df(1,ip,k)
                        uavg = uavg + q_df(2,ip,k)
                        vavg = vavg + q_df(3,ip,k)
                    end do
                end do
                pavg = pavg / npts
                uavg = uavg / npts
                vavg = vavg / npts

                if (pavg < 0.0) then
                    write(*,*) 'Negative layer thickness in element  ', e, ',  layer  ', k
                    stop
                end if

                call get_contraction(contraction,e,k,pavg,pmin)

                do j = 1, ngly
                    do i = 1, nglx

                        ip = intma(i,j,1,e)

                        dp_adjust = (1.0-contraction)*q_df(1,ip,k)
                        q_df(1,ip,k) = q_df(1,ip,k) - dp_adjust
                        q_df(1,ip,k-1) = q_df(1,ip,k-1) + dp_adjust

                        ! u_adjust = (1.0-contraction)*q_df(2,ip,k)
                        ! q_df(2,ip,k) = q_df(2,ip,k) - u_adjust
                        ! q_df(2,ip,k-1) = q_df(2,ip,k-1) + u_adjust

                        ! v_adjust = (1.0-contraction)*q_df(3,ip,k)
                        ! q_df(3,ip,k) = q_df(3,ip,k) - v_adjust
                        ! q_df(3,ip,k-1) = q_df(3,ip,k-1) + v_adjust
                    end do
                end do
            end do
        enddo

        ! Then sweep downward from the surface

        do e = 1,nelem

            do k = 1, nlayers-1
                pmin = 1.0e40
                pavg = 0.0
                uavg = 0.0
                vavg = 0.0
                do m = 1, ngly
                    do n = 1, nglx
                        ip = intma(n,m,1,e)
                        if (q_df(1,ip,k) < pmin) pmin = q_df(1,ip,k)

                        pavg = pavg + q_df(1,ip,k)
                        uavg = uavg + q_df(2,ip,k)
                        vavg = vavg + q_df(3,ip,k)
                    end do
                end do
                pavg = pavg / npts
                uavg = uavg / npts
                vavg = vavg / npts

                if (pavg < 0.0) then
                    write(*,*) 'Negative layer thickness in element  ', e, ',  layer  ', k
                    stop
                end if

                call get_contraction(contraction,e,k,pavg,pmin)

                do j = 1, ngly
                    do i = 1, nglx

                        ip = intma(i,j,1,e)

                        dp_adjust = (1.0-contraction)*q_df(1,ip,k)
                        q_df(1,ip,k) = q_df(1,ip,k) - dp_adjust
                        q_df(1,ip,k+1) = q_df(1,ip,k+1) + dp_adjust

                        ! u_adjust = (1.0-contraction)*q_df(2,ip,k)
                        ! q_df(2,ip,k) = q_df(2,ip,k) - u_adjust
                        ! q_df(2,ip,k+1) = q_df(2,ip,k+1) + u_adjust

                        ! v_adjust = (1.0-contraction)*q_df(3,ip,k)
                        ! q_df(3,ip,k) = q_df(3,ip,k) - v_adjust
                        ! q_df(3,ip,k+1) = q_df(3,ip,k+1) + v_adjust
                    end do
                end do
            end do
        enddo

        ! Compute revised values of  dp',  based on the revised values of  dp.

        ! do ip = 1, npoin

        !     ope = sum(q_df(1,ip,:)) / pbprime_df(ip)

        !     do k = 1,nlayers

        !         qprime_df(1,ip,k) = q_df(1,ip,k) / ope
        !         ! qprime_df(2,ip,k) = uv_df(1,ip,k) - qb_df(3,ip)/qb_df(1,ip)
        !         ! qprime_df(3,ip,k) = uv_df(2,ip,k) - qb_df(4,ip)/qb_df(1,ip)
        !     enddo

        ! end do

    end subroutine limit_layer_mass

    subroutine get_contraction(contraction,i,k, dp_avg, dp_min)
        
        ! This routine computes the contraction factor for layer mass limiting
        use mod_input, only: dry_cutoff
        use mod_initial, only: alpha_mlswe
        use mod_constants, only: gravity

        implicit none

        integer, intent(in) :: i, k
        real, intent(in) :: dp_avg, dp_min
        real, intent(out) :: contraction

        real :: min_dp_frac, dp_contraction

        min_dp_frac = 0.
        dp_contraction = 0.0e3
        ! dp_contraction = (gravity/alpha_mlswe(k)) * dry_cutoff

        ! if (dp_avg < dp_contraction) then
        !     contraction = 0

        ! else if (dp_min >= min_dp_frac*dp_avg) then
        !     contraction = 1

        ! else
        !     contraction = (dp_avg - min_dp_frac*dp_avg) / (dp_avg - dp_min)  
        ! endif  

        if (dp_avg < dp_contraction) then
            contraction = 0

        else if (dp_min >= min_dp_frac*dp_avg) then
            contraction = 1

        else
            contraction = (dp_avg - min_dp_frac*dp_avg) / (dp_avg - dp_min)  
        endif  


    end subroutine get_contraction


    subroutine  limit_layer_momentum(q_df)

        use mod_input, only: nlayers
        use mod_initial, only: alpha_mlswe, zbot_df
        use mod_constants, only: gravity
        use mod_variables, only: z_interface_initial, z_init_flag_elem, frac_zdiff1, frac_zdiff2
        use mod_basis, only: nglx, ngly
        use mod_grid, only:  nelem, npoin, intma

        real, dimension(3,npoin,nlayers), intent(inout) :: q_df

        real :: z_init_max, z_now_max, zdiff_max, zdiff1, zdiff2, tmp, mom_limiter
        real, dimension(npoin,nlayers+1) :: z_interface_now
        integer :: i, j, k, m, n, e, ip, I1, I2

        frac_zdiff1  = 0.5  !  used to limit momentum at incrops on topog  
        frac_zdiff2 = 1.0   !  used to limit momentum at incrops on topog

        do e = 1,nelem

            ! Determine whether some limiting needs to be done on grid cell  e.

            if (minval(z_init_flag_elem(e,:)) == 0) then

                do j = 1, ngly
                    do i = 1, nglx
                        ip = intma(i,j,1,e)
                        z_interface_now(ip,nlayers+1) = zbot_df(ip)
                        do k = nlayers, 1, -1
                            z_interface_now(ip,k) = z_interface_now(ip,k+1) + (alpha_mlswe(k)/gravity)*q_df(1,ip,k)
                        end do
                    end do
                end do

                do k = 2, nlayers+1
                    if (z_init_flag_elem(e,k-1) == 0) then

                        z_init_max = -1.0e30
                        z_now_max = -1.0e30
                        do m = 1, ngly
                            do n = 1, nglx
                                ip = intma(n,m,1,e)
                                if (z_interface_initial(ip,k) > z_init_max) z_init_max = z_interface_initial(ip,k)
                                if (z_interface_now(ip,k) > z_now_max) z_now_max = z_interface_now(ip,k)
                            end do
                        end do

                        zdiff_max = abs(z_init_max - z_now_max)
                        
                        do j = 1, ngly
                            I1 = intma(1,j,1,e)
                            I2 = intma(nglx,j,1,e)
                            tmp = abs(z_interface_initial(I1,k-1) - z_interface_initial(I2,k-1))

                            zdiff1 = frac_zdiff1*tmp 
                            zdiff2 = frac_zdiff2*tmp
                            tmp = (zdiff_max - zdiff1) / (zdiff2 - zdiff1)
                            mom_limiter = max(0.0, min(1.0, tmp))

                            do i = 1, nglx
                                ip = intma(i,j,1,e)
                                q_df(2,ip,k) = mom_limiter*q_df(2,ip,k)
                                q_df(3,ip,k) = mom_limiter*q_df(3,ip,k)
                            end do

                            if (zdiff_max > zdiff2) z_init_flag_elem(e,k-1) = 1
                        end do

                        ! y-direction limiting
                        do i = 1, nglx
                            I1 = intma(i,1,1,e)
                            I2 = intma(i,ngly,1,e)
                            tmp = abs(z_interface_initial(I1,k-1) - z_interface_initial(I2,k-1))

                            zdiff1 = frac_zdiff1*tmp 
                            zdiff2 = frac_zdiff2*tmp
                            tmp = (zdiff_max - zdiff1) / (zdiff2 - zdiff1)
                            mom_limiter = max(0.0, min(1.0, tmp))

                            do j = 1, ngly
                                ip = intma(i,j,1,e)
                                q_df(2,ip,k) = mom_limiter*q_df(2,ip,k)
                                q_df(3,ip,k) = mom_limiter*q_df(3,ip,k)
                            end do
                            if (zdiff_max > zdiff2) z_init_flag_elem(e,k-1) = 1
                        enddo

                    end if
                end do
            end if
        end do
    end subroutine limit_layer_momentum

end module mod_layer_terms
