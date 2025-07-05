! ================================================================================================
! This module contains the routines for the barotropic flux terms
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: March 27, 2023
! ================================================================================================

module mod_barotropic_terms

    
    use mod_grid, only: npoin_q, nface, intma_dg_quad
    use mod_basis, only: nqx, nqy, nqz, nq
    use mod_input, only: nlayers
        
    implicit none

    public :: &
                btp_mom_boundary_df, &
                compute_gradient_uv, btp_extract_df, btp_extract_face, &
                btp_bcl_coeffs_qdf, compute_gradient_uv_q

    contains

    subroutine btp_extract_df(qb_df_face, qb_df)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl, nq, psiq
        use mod_grid, only:  npoin, intma, nface, face
        use mod_face, only: imapl, imapr, normal_vector

        implicit none

        real, intent(inout) :: qb_df_face(4,2,ngl,nface)
        real, intent(in) :: qb_df(4,npoin)

        integer :: iface, ilr, iquad,jquad, n, m, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr
        real :: un, nx, ny, hi
      
        ! Compute values of certain barotropic functions at quadrature points
        ! in each cell and endpoints of each cell, at barotropic
        ! time level m+1.
        ! The values at endpoints should be interpreted as one-sided limits.
      
        ! In the case of  pbub,  pbvb,  and  pbpert,  use degrees of freedom 
        ! and basis functions.  In the case of  ub  and  vb,  use division.
      
        qb_df_face = 0.0
      
        ! Evaluate  pbub  and  pbvb.
        do iface = 1, nface                  !i specifies the grid cell

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do n = 1,ngl

                ! Left
                il=imapl(1,n,1,iface)
                jl=imapl(2,n,1,iface)
                kl=imapl(3,n,1,iface)
                I=intma(il,jl,kl,el)

                qb_df_face(1:4,1,n,iface) = qb_df(1:4,I)

                if(er > 0) then

                    ! Right
                    ir=imapr(1,n,1,iface)
                    jr=imapr(2,n,1,iface)
                    kr=imapr(3,n,1,iface)
                    I=intma(ir,jr,kr,er)

                    qb_df_face(1:4,2,n,iface) = qb_df(1:4,I)
                else

                    qb_df_face(1:4,2,n,iface) = qb_df_face(1:4,1,n,iface)

                    if(er == -4) then

                        nx = normal_vector(1,n,1,iface)
                        ny = normal_vector(2,n,1,iface)

                        un = nx*qb_df_face(3,1,n,iface) + ny*qb_df_face(4,1,n,iface)

                        qb_df_face(3,2,n,iface) = qb_df_face(3,1,n,iface) - 2.0*un*nx
                        qb_df_face(4,2,n,iface) = qb_df_face(4,1,n,iface) - 2.0*un*ny

                    elseif(er == -2) then 
                        qb_df_face(3:4,2,n,iface) = -qb_df_face(3:4,1,n,iface)
                        
                    end if
                end if
            end do
        end do 
      
    end subroutine btp_extract_df 

    subroutine btp_extract_face(qb_face, qb_df_face)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl, nq, psiq
        use mod_grid, only:  npoin, intma, nface, face
        use mod_face, only: imapl, imapr, normal_vector_q

        implicit none

        real, intent(out) :: qb_face(4,2,nq,nface)
        real, intent(in) :: qb_df_face(4,2,ngl,nface)

        integer :: iface, ilr, iquad,jquad, n, m, il, jl, ir, jr, el, er, ilocl, ilocr, I,kl,kr
        real :: un, nx, ny, hi

        ! Compute values of certain barotropic functions at quadrature points
        ! in each cell and endpoints of each cell, at barotropic
        ! time level m+1.
        ! The values at endpoints should be interpreted as one-sided limits.
      
        ! In the case of  pbub,  pbvb,  and  pbpert,  use degrees of freedom 
        ! and basis functions.  In the case of  ub  and  vb,  use division.
      
        qb_face = 0.0
      
        ! Evaluate  pbub  and  pbvb.
        do iface = 1, nface                  !i specifies the grid cell

            !Store Left Side Variables
            el = face(7,iface)
            er = face(8,iface)

            do iquad = 1,nq
                ! Left
                do n = 1,ngl

                    hi = psiq(n,iquad)

                    il=imapl(1,n,1,iface)
                    jl=imapl(2,n,1,iface)
                    kl=imapl(3,n,1,iface)
                    I=intma(il,jl,kl,el)

                    qb_face(1:4,1,iquad,iface) = qb_face(1:4,1,iquad,iface) &
                                                + hi*qb_df_face(1:4,1,n,iface)
                end do 
                ! Right
                if(er > 0) then
                    do n = 1,ngl

                        hi = psiq(n,iquad)

                        ir=imapr(1,n,1,iface)
                        jr=imapr(2,n,1,iface)
                        kr=imapr(3,n,1,iface)
                        I=intma(ir,jr,kr,er)

                        qb_face(1:4,2,iquad,iface) = qb_face(1:4,2,iquad,iface) &
                                                    + hi*qb_df_face(1:4,2,n,iface)
                    end do 

                end if
            end do
        end do 
      
    end subroutine btp_extract_face

    subroutine btp_mom_boundary_df(qb)

        use mod_basis, only: nglx, ngly, nqx, nqy, nqz, ngl,nq
        use mod_grid, only:  npoin, intma, nface, face,mod_grid_get_face_nq
        use mod_initial, only: pbprime_face
        use mod_face, only: imapl, imapr, normal_vector

        implicit none

        real, intent(inout) :: qb(2,npoin)

        integer :: iface, ilr, m, il, jl, el, er, ilocl, ilocr, I, kl
        integer :: bcflag, n
        real :: nx, ny, unl, upnl

        do iface = 1, nface                  !i specifies the grid cell

            !Store Left Side Variables
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

                    unl = qb(1,I)*nx + qb(2,I)*ny
                    qb(1,I) = qb(1,I) - unl*nx
                    qb(2,I) = qb(2,I) - unl*ny
                end do

            elseif(er == -2) then
                do n = 1, ngl

                    il=imapl(1,n,1,iface)
                    jl=imapl(2,n,1,iface)
                    kl=imapl(3,n,1,iface)
                    I=intma(il,jl,kl,el)

                    qb(1,I) = 0.0
                    qb(2,I) = 0.0

                enddo
            end if
        end do
      
    end subroutine btp_mom_boundary_df
 
    subroutine btp_bcl_coeffs_qdf(qprime_df_face, qprime_df)

        ! Compute baroclinic coefficients in the advective barotropic momentum fluxes 
        ! and in the barotropic pressure forcing.

        use mod_grid, only: npoin_q, nface, npoin, face, intma
        use mod_input, only: nlayers
        use mod_basis, only: nqx, nqy, nqz, nq, ngl, npts, psiq
        use mod_initial, only: alpha_mlswe, indexq, psih
        use mod_Tensorproduct, only: compute_gradient_quad, interpolate_layer_from_quad_to_node_1d
        use mod_variables, only: Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, H_bcl_edge, &
                                    Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, &
                                    dpprime_visc,pbprime_visc,btp_dpp_graduv, btp_dpp_uvp, &
                                    dpp_uvp, dpp_graduv, graduv_dpp_face, btp_graduv_dpp_face 
        use mod_face, only: imapl_q, imapr_q, normal_vector_q, imapl, imapr, normal_vector

        implicit none

        real, dimension(3,2,ngl,nface,nlayers),intent(in)  :: qprime_df_face
        real, dimension(3,npoin,nlayers), intent(in) :: qprime_df

        real, dimension(nlayers+1) :: pprime_l, pprime_r
        real, dimension(nlayers+1) :: pprime
        real, dimension(npoin,nlayers+1) :: pprime_df
        real :: left_dp, right_dp
        real, dimension(4, npoin) :: graduv
        real :: qq(3), ql(3), qr(3)

        integer :: iface, il, jl, kl, ir, jr, kr, iel, ier, iquad, Iq, k, ivar, n, I, ip
        real :: un(nlayers), nx, ny, hi

        ! Compute Q_uu_dp, Q_uv_dp, Q_vv_dp at quadrature points
        Q_uu_dp = 0.0
        Q_uv_dp = 0.0
        Q_vv_dp = 0.0
        H_bcl = 0.0
        Q_uu_dp_edge = 0.0
        Q_uv_dp_edge = 0.0
        Q_vv_dp_edge = 0.0
        H_bcl_edge = 0.0

        btp_dpp_graduv = 0.0
        pbprime_visc = 0.0

        ! Compute Q_up_up_quad and Q_up_vp_quad at quadrature points.

        do Iq = 1,npoin_q
            pprime(1) = 0.0
            do k = 1, nlayers
                qq = 0.0
                do ip = 1,npts
                    I = indexq(ip,Iq)
                    hi = psih(ip,Iq)
                    qq(:) = qq(:) + hi*qprime_df(:,I,k)
                enddo
                Q_uu_dp(Iq) = Q_uu_dp(Iq) + qq(2)*(qq(2)*qq(1))
                Q_uv_dp(Iq) = Q_uv_dp(Iq) + qq(3)*(qq(2)*qq(1))
                Q_vv_dp(Iq) = Q_vv_dp(Iq) + qq(3)*(qq(3)*qq(1))

                ! Pressure term in barotropic momentum equation

                pprime(k+1) = pprime(k) + qq(1)
                H_bcl(Iq) = H_bcl(Iq) + 0.5*alpha_mlswe(k)*(pprime(k+1)**2 - pprime(k)**2)
            end do
        enddo

        pprime_df(:,1) = 0.0

        do k = 1,nlayers

            ! For viscosity terms

            call compute_gradient_uv(graduv, qprime_df(2:3,:,k))

            dpp_graduv(1,:,k) = dpprime_visc(:,k)*graduv(1,:)
            dpp_graduv(2,:,k) = dpprime_visc(:,k)*graduv(2,:)
            dpp_graduv(3,:,k) = dpprime_visc(:,k)*graduv(3,:)
            dpp_graduv(4,:,k) = dpprime_visc(:,k)*graduv(4,:)

            dpp_uvp(1,:,k) =  dpprime_visc(:,k)*qprime_df(2,:,k)
            dpp_uvp(2,:,k) =  dpprime_visc(:,k)*qprime_df(3,:,k)

            btp_dpp_graduv(:,:) = btp_dpp_graduv(:,:) + dpp_graduv(:,:,k)
            pbprime_visc(:) =  pbprime_visc(:) + dpprime_visc(:,k)

        end do

        do iface = 1, nface

            iel=face(7,iface)
            ier=face(8,iface)

            do iquad = 1,nq

                pprime_l(1) = 0.0
                pprime_r(1) = 0.0
                do k = 1, nlayers
                    ql = 0.0; qr = 0.0;
                    do n = 1, ngl
                        hi = psiq(n,iquad)
                        ql(:) = ql(:) + hi*qprime_df_face(:,1,n,iface,k)
                        qr(:) = qr(:) + hi*qprime_df_face(:,2,n,iface,k)
                    enddo

                    Q_uu_dp_edge(iquad,iface) = Q_uu_dp_edge(iquad,iface) &
                                                + 0.5*((ql(2)*ql(2)*ql(1)) + (qr(2)*qr(2)*qr(1)))
                    Q_uv_dp_edge(iquad,iface) = Q_uv_dp_edge(iquad,iface) &
                                                + 0.5*((ql(3)*ql(2)*ql(1)) + (qr(3)*qr(2)*qr(1)))
                    Q_vv_dp_edge(iquad,iface) = Q_vv_dp_edge(iquad,iface) &
                                                + 0.5*((ql(3)*ql(3)*ql(1)) + (qr(3)*qr(3)*qr(1)))

                    pprime_l(k+1) = pprime_l(k) + ql(1)
                    left_dp = 0.5*alpha_mlswe(k) *(pprime_l(k+1)**2 - pprime_l(k)**2)
                    pprime_r(k+1) = pprime_r(k) + qr(1)
                    right_dp = 0.5*alpha_mlswe(k) *(pprime_r(k+1)**2 - pprime_r(k)**2)

                    H_bcl_edge(iquad,iface) = H_bcl_edge(iquad,iface) + 0.5*(left_dp + right_dp)
                enddo
            end do

            do iquad = 1,ngl
                !Get Pointers
                il=imapl(1,iquad,1,iface)
                jl=imapl(2,iquad,1,iface)
                kl=imapl(3,iquad,1,iface)
                Iq = intma(il,jl,kl,iel)

                !Variables
                graduv_dpp_face(1,1,iquad,iface,:) = dpp_graduv(1,Iq,:)
                graduv_dpp_face(2,1,iquad,iface,:) = dpp_graduv(2,Iq,:)
                graduv_dpp_face(3,1,iquad,iface,:) = dpp_graduv(3,Iq,:)
                graduv_dpp_face(4,1,iquad,iface,:) = dpp_graduv(4,Iq,:)
                graduv_dpp_face(5,1,iquad,iface,:) = dpprime_visc(Iq,:)

                if (ier > 0 ) then

                    !Get Pointers
                    ir=imapr(1,iquad,1,iface)
                    jr=imapr(2,iquad,1,iface)
                    kr=imapr(3,iquad,1,iface)
                    Iq=intma(ir,jr,kr,ier)

                    !Variables
                    graduv_dpp_face(1,2,iquad,iface,:) = dpp_graduv(1,Iq,:)
                    graduv_dpp_face(2,2,iquad,iface,:) = dpp_graduv(2,Iq,:)
                    graduv_dpp_face(3,2,iquad,iface,:) = dpp_graduv(3,Iq,:)
                    graduv_dpp_face(4,2,iquad,iface,:) = dpp_graduv(4,Iq,:)
                    graduv_dpp_face(5,2,iquad,iface,:) = dpprime_visc(Iq,:)

                else
                    !default values
                    graduv_dpp_face(1,2,iquad,iface,:) = graduv_dpp_face(1,1,iquad,iface,:)
                    graduv_dpp_face(2,2,iquad,iface,:) = graduv_dpp_face(2,1,iquad,iface,:)
                    graduv_dpp_face(3,2,iquad,iface,:) = graduv_dpp_face(3,1,iquad,iface,:)
                    graduv_dpp_face(4,2,iquad,iface,:) = graduv_dpp_face(4,1,iquad,iface,:)
                    graduv_dpp_face(5,2,iquad,iface,:) = graduv_dpp_face(5,1,iquad,iface,:)

                    if(ier == -4) then
                        nx = normal_vector(1,iquad,1,iface)
                        ny = normal_vector(2,iquad,1,iface)

                        un = dpp_graduv(1,Iq,:)*nx + dpp_graduv(2,Iq,:)*ny
                        graduv_dpp_face(1,2,iquad,iface,:) = dpp_graduv(1,Iq,:) - 2.0*un*nx
                        graduv_dpp_face(2,2,iquad,iface,:) = dpp_graduv(2,Iq,:) - 2.0*un*ny

                        un = dpp_graduv(3,Iq,:)*nx + dpp_graduv(4,Iq,:)*ny
                        graduv_dpp_face(3,2,iquad,iface,:) = dpp_graduv(3,Iq,:) - 2.0*un*nx
                        graduv_dpp_face(4,2,iquad,iface,:) = dpp_graduv(4,Iq,:) - 2.0*un*ny

                    end if
                end if
            end do
        end do

        call bcl_create_communicator(graduv_dpp_face,5,nlayers,ngl)

        btp_graduv_dpp_face = 0.0

        do iface=1,nface
            do iquad = 1,ngl
                do k = 1,nlayers

                    btp_graduv_dpp_face(:,1,iquad,iface) = btp_graduv_dpp_face(:,1,iquad,iface) &
                                                            + graduv_dpp_face(:,1,iquad,iface,k)
                    btp_graduv_dpp_face(:,2,iquad,iface) = btp_graduv_dpp_face(:,2,iquad,iface) &
                                                            + graduv_dpp_face(:,2,iquad,iface,k)
                end do
            end do
        end do

    end subroutine btp_bcl_coeffs_qdf

    subroutine compute_gradient_uv(grad_uv,uv)

        use mod_basis, only:  npts
        use mod_grid, only:  npoin
        use mod_initial, only: dpsidx_df,dpsidy_df, index_df

        implicit none

        real, dimension(2,npoin), intent(in) :: uv
        real, dimension(4,npoin), intent(out) :: grad_uv
        integer :: Iq, I, ip
        real :: dhdx, dhdy
        
        grad_uv = 0.0
        
        !Construct Volume Integral Contributions
        do Iq = 1, npoin

            do ip = 1,npts

                I = index_df(ip,Iq)
                dhdx = dpsidx_df(ip,Iq)
                dhdy = dpsidy_df(ip,Iq)

                grad_uv(1,Iq) = grad_uv(1,Iq) + dhdx*uv(1,I)
                grad_uv(2,Iq) = grad_uv(2,Iq) + dhdy*uv(1,I)
                grad_uv(3,Iq) = grad_uv(3,Iq) + dhdx*uv(2,I)
                grad_uv(4,Iq) = grad_uv(4,Iq) + dhdy*uv(2,I)

            end do
        end do

    end subroutine compute_gradient_uv

    subroutine compute_gradient_uv_q(grad_uv,uv)

        use mod_basis, only: npts
        use mod_grid, only:  npoin_q, npoin
        use mod_initial, only: dpsidx,dpsidy, indexq

        implicit none

        real, dimension(2,npoin), intent(in) :: uv
        real, dimension(2,2,npoin_q), intent(out) :: grad_uv
        integer :: Iq, I, ip
        real :: dhdx, dhdy

        grad_uv = 0.0

        !Construct Volume Integral Contributions
        do Iq = 1, npoin_q

            do ip = 1,npts

                I = indexq(ip,Iq)
                dhdx = dpsidx(ip,Iq)
                dhdy = dpsidy(ip,Iq)

                grad_uv(1,1,Iq) = grad_uv(1,1,Iq) + dhdx*uv(1,I)
                grad_uv(1,2,Iq) = grad_uv(1,2,Iq) + dhdy*uv(1,I)
                grad_uv(2,1,Iq) = grad_uv(2,1,Iq) + dhdx*uv(2,I)
                grad_uv(2,2,Iq) = grad_uv(2,2,Iq) + dhdy*uv(2,I)

            end do
        end do

    end subroutine compute_gradient_uv_q

end module mod_barotropic_terms         
