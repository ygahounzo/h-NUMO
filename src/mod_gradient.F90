!----------------------------------------------------------------------!
!>@brief This module contains basic operators; derivatives (1st and 2nd),
!> curl, divergence and laplacian
!>@author  Francis X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@modified  Yao Gahounzo 04/03/2023
!>           Computing PhD
!>           Boise State University
!----------------------------------------------------------------------!
module mod_gradient

    use mod_basis

    public :: &
        compute_local_gradient, &
        compute_local_gradient_v3, &
        compute_local_gradient_quad_v3, &
        compute_vorticity_mlswe

    private

contains

    !----------------------------------------------------------------------!
    !>@brief Relative vorticity ζ, absolute vorticity η=ζ+f, and potential
    !> vorticity q=η/h for one MLSWE layer.  ζ is the local-vertical (radial,
    !> on the sphere) component
    !> of the curl of the layer velocity (u,v,w), computed via the same DG
    !> element-local strong-form derivative pattern as btp_bcl_coeffs_qdf
    !> (mod_barotropic_terms.F90).  One-shot diagnostic snapshot, so no LDG
    !> face correction or mass-matrix inversion is needed (cf. gradz in
    !> mod_create_rhs_mlswe.F90, same idea).
    !> w is unused (may hold anything, e.g. a placeholder) when has_w is
    !> .false. -- matches the Cartesian/no-w case, where kvector=(0,0,1) and
    !> ζ reduces to the standard planar ζ=∂v/∂x-∂u/∂y.
    !----------------------------------------------------------------------!
    subroutine compute_vorticity_mlswe(G, b, tsp, init, has_w, u, v, w, h, rel_vort, abs_vort, pot_vort)

        use mod_grid,    only: grid
        use mod_tensor,  only: tensor_CS
        use mod_initial, only: initial

        implicit none

        type(grid),      intent(in)  :: G
        type(basis),     intent(in)  :: b
        type(tensor_CS), intent(in)  :: tsp
        type(initial),   intent(in)  :: init
        logical,         intent(in)  :: has_w
        real,            intent(in)  :: u(G%npoin), v(G%npoin), w(G%npoin), h(G%npoin)
        real,            intent(out) :: rel_vort(G%npoin), abs_vort(G%npoin), pot_vort(G%npoin)

        integer :: Iq, ip, I
        real :: dhdx, dhdy, dhdz
        real :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz
        real :: zeta

        do Iq = 1, G%npoin
            du_dx = 0.0; du_dy = 0.0; du_dz = 0.0
            dv_dx = 0.0; dv_dy = 0.0; dv_dz = 0.0
            dw_dx = 0.0; dw_dy = 0.0; dw_dz = 0.0
            do ip = 1, b%npts
                I    = tsp%index_df(ip,Iq)
                dhdx = tsp%dpsidx_df(ip,Iq)
                dhdy = tsp%dpsidy_df(ip,Iq)
                dhdz = 0.0
                if (has_w) then
                    dhdx = dhdx + tsp%dpsidz_df_x(ip,Iq)
                    dhdy = dhdy + tsp%dpsidz_df_y(ip,Iq)
                    dhdz = tsp%dpsidz_df(ip,Iq) + tsp%dpsidz_df_z(ip,Iq)
                end if
                du_dx = du_dx + dhdx*u(I)
                du_dy = du_dy + dhdy*u(I)
                dv_dx = dv_dx + dhdx*v(I)
                dv_dy = dv_dy + dhdy*v(I)
                if (has_w) then
                    du_dz = du_dz + dhdz*u(I)
                    dv_dz = dv_dz + dhdz*v(I)
                    dw_dx = dw_dx + dhdx*w(I)
                    dw_dy = dw_dy + dhdy*w(I)
                    dw_dz = dw_dz + dhdz*w(I)
                end if
            end do

            ! ζ = k̂·(∇×u): reduces to the standard planar ζ=dv/dx-du/dy for
            ! Cartesian geometry, since kvector=(0,0,1) and has_w=.false. there.
            zeta = init%kvector(1,Iq)*(dw_dy - dv_dz) &
                 + init%kvector(2,Iq)*(du_dz - dw_dx) &
                 + init%kvector(3,Iq)*(dv_dx - du_dy)

            rel_vort(Iq) = zeta
            abs_vort(Iq) = zeta + init%coriolis_df(Iq)
            pot_vort(Iq) = abs_vort(Iq) / h(Iq)
        end do

    end subroutine compute_vorticity_mlswe

    !----------------------------------------------------------------------!
    !>@brief This subroutine computes the canonical derivative of the vector Q as one NxNDIMxN system.
    !>That is, it builds Q_e, Q_n, and Q_c.
    !>@author  Francis X. Giraldo on 12/2013
    !>                  Department of Applied Mathematics
    !>                  Naval Postgraduate School
    !>                  Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_local_gradient(q_e,q_n,q_c,q,b,ndim)

        implicit none

        type(basis), intent(in) :: b

        !global arrays
        real, intent(out) :: q_e(b%nglx,ndim,b%ngly,b%nglz)
        real, intent(out) :: q_n(b%nglx,ndim,b%ngly,b%nglz)
        real, intent(out) :: q_c(b%nglx,ndim,b%ngly,b%nglz)
        real, intent(in)  :: q(b%nglx,ndim,b%ngly,b%nglz)
        integer, intent(in) :: ndim

        !local arrays
        integer :: m1, m2, m3, k

        if(b%nglx > 1) then
            m1=b%nglx !nglx
            m2=b%nglx !nglx
            m3=ndim*b%ngly*b%nglz !nsize*ngly*nglz
            call mxm(b%dpsix_tr,m1,q,m2,q_e,m3)
        else
            q_e = 0
        endif

        if(b%ngly > 1) then
            m1=b%nglx*ndim !nglx*nsize
            m2=b%ngly !ngly
            m3=b%ngly !ngly
            do k=1,b%nglz !nglz
                call mxm(q(1,1,1,k),m1,b%dpsiy,m2,q_n(1,1,1,k),m3)
            enddo !k
        else
            q_n = 0
        endif

        if(b%nglz > 1) then
            m1=b%nglx*ndim*b%ngly !nglx*nsize*ngly
            m2=b%nglz !nglz
            m3=b%nglz !nglz
            call mxm(q,m1,b%dpsiz,m2,q_c,m3)
        else
            q_c = 0
        endif

    end subroutine compute_local_gradient

    !----------------------------------------------------------------------!
    !>@brief This subroutine computes the canonical derivative of the vector Q as (NDIM,N,N,N)
    !>That is, it builds Q_e, Q_n, and Q_c.
    !>@author  Francis X. Giraldo on 12/2013
    !>                  Department of Applied Mathematics
    !>                  Naval Postgraduate School
    !>                  Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine compute_local_gradient_v3(q_e,q_n,q_c,q,b,ndim)

        implicit none

        type(basis), intent(in) :: b

        !global arrays
        real, dimension(ndim,b%nglx,b%ngly,b%nglz), intent(out) :: q_e, q_n, q_c
        real, dimension(ndim,b%nglx,b%ngly,b%nglz), intent(in)  :: q
        integer, intent(in) :: ndim

        !local arrays
        integer :: m1, m2, m3, m, k
        real, dimension(b%nglx,b%ngly,b%nglz) :: qt, qt_e, qt_n, qt_c

        !Construct Local Derivatives
        do m=1,ndim
            !Store Input
            qt=q(m,:,:,:)

            !KSI Derivative
            if(b%nglx > 1) then
                m1=b%nglx !nglx
                m2=b%nglx !nglx
                m3=b%ngly*b%nglz !nsize*ngly*nglz
                call mxm(b%dpsix_tr,m1,qt,m2,qt_e,m3)
            else
                qt_e = 0
            endif

            !ETA Derivative
            if(b%ngly > 1) then
                m1=b%nglx !nglx
                m2=b%ngly !ngly
                m3=b%ngly !ngly
                do k=1,b%nglz
                    call mxm(qt(1,1,k),m1,b%dpsiy,m2,qt_n(1,1,k),m3)
                enddo !k
            else
                qt_n = 0
            endif

            !Zeta Derivative
            if(b%nglz > 1) then
                m1=b%nglx*b%ngly !nglx*nsize*ngly
                m2=b%nglz !nglz
                m3=b%nglz !nglz
                call mxm(qt,m1,b%dpsiz,m2,qt_c,m3)
            else
                qt_c = 0
            endif

            !Store Output
            q_e(m,:,:,:)=qt_e
            q_n(m,:,:,:)=qt_n
            q_c(m,:,:,:)=qt_c
        end do !m

    end subroutine compute_local_gradient_v3

    !----------------------------------------------------------------------!
    !>@brief This subroutine computes the canonical derivative of the vector Q as (Nx,Ny,Nz)
    !>That is, it builds Q_e, Q_n, and Q_c.
    !>@author  Yao Gahounzo 04/03/2023
    !>                  Computing PhD
    !>                  Boise State University
    !----------------------------------------------------------------------!
    subroutine compute_local_gradient_quad_v3(q_e,q_n,q_c,q,b)

        implicit none

        type(basis), intent(in) :: b

        !global arrays
        real, dimension(b%nqx,b%nqy,b%nqz), intent(out) :: q_e, q_n, q_c
        real, dimension(b%nglx,b%ngly,b%nglz), intent(in)  :: q
    
        !local arrays
        integer :: m1, m2, m3, m, k, l, kquad, jquad, iquad,n
        real :: hix, hiy, temp

        q_e = 0.0
        q_n = 0.0
        q_c = 0.0

        !Construct Local Derivatives

        !KSI Derivative

        !open(1,file='hix.dat',status='unknown')

        do kquad = 1,b%nqz
            do jquad = 1,b%nqy
                do iquad = 1,b%nqx

                    !temp = 0.0
                    
                    do l = 1,b%nglz
                        do m = 1,b%ngly
                            do n = 1,b%nglx

                                q_e(iquad,jquad,kquad) = q_e(iquad,jquad,kquad) &
                                                        + b%dpsiqx(n,iquad)*b%psiqy(m,jquad)*q(n,m,l)
                                q_n(iquad,jquad,kquad) = q_n(iquad,jquad,kquad) &
                                                        + b%psiqx(n,iquad)*b%dpsiqy(m,jquad)*q(n,m,l)

                                if(b%nglz > 1) then
                                    q_c(iquad,jquad,kquad) = q_c(iquad,jquad,kquad) &
                                            + b%psiqx(n,iquad)*b%psiqy(m,jquad)*b%dpsiz(l,kquad)*q(n,m,l)
                                endif
                                
                            enddo !n
                        enddo !m
                    enddo !l

                enddo !iquad
            enddo !jquad
        enddo !kquad

    end subroutine compute_local_gradient_quad_v3

 end module mod_gradient
