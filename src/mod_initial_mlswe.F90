!----------------------------------------------------------------------!
!>@brief This module contains Initial Conditions settings
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: March 27, 2023
!----------------------------------------------------------------------!
module mod_initial_mlswe

    use mod_input,   only: input

    use mod_basis,   only: basis
    use mod_grid,    only: grid
    use mod_metrics, only: metrics
    use mod_face,    only: face_CS
    use mod_tensor,  only: tensor_CS

    public :: &
        bot_topo_derivatives, &
        interpolate_pbprime_init, wind_stress_coriolis, &
        map_deriv, ssprk_coefficients, poslimiter, find_dry_elements, &
        check_layer_thickness, check_btp_thickness, &
        bcl_itime, trace_stage, trace_reset, trace_bcl_stage

    private

    ! Current BCL time-step index; set by the time loop before each BCL call
    ! so that check_btp_thickness can report "BTP failure at step N".
    integer, save :: trace_bcl_stage = 0
    integer, save :: trace_unit = 0, trace_sequence = 0
    logical, save :: trace_open = .false., trace_dumped = .false.
    character(len=2048), save :: trace_ring(128)
    integer, save :: trace_ring_count = 0
    real, allocatable, save :: trace_before(:,:,:)
    integer, save :: trace_before_step = -1, trace_before_stage = -1
    integer, save :: bcl_itime = 0


    contains
    
    !> Compute the gradient of the bottom topography at quadrature points
    subroutine bot_topo_derivatives(b, G, inp, mf, zbot,zbot_face,zbot_df)

        implicit none

        type(basis), intent(in) :: b
        type(grid), intent(in) :: G
        type(input), intent(in) :: inp
        type(face_CS), intent(in) :: mf

        real, dimension(G%npoin), intent(in) :: zbot_df
        real, dimension(G%npoin_q), intent(out) :: zbot
        real, dimension(2,b%nq,G%nface), intent(out) :: zbot_face
        
        ! Declare local variables
        integer :: il, jl, ir, jr,kl,kr, el, er, ilocl, ilocr, l, I1,I2, iface, &
                    iquad, jquad, nq_i, nq_j, plane_ij
        integer :: e,I, m,n,Iq
        real :: depth, x, y, r, xm, yl, hi

        
        ! Initialize output arrays

        do e = 1, G%nelem
            do jquad = 1, b%nqy
                do iquad = 1, b%nqx
                
                    Iq = G%intma_dg_quad(iquad, jquad, 1, e)

                    do m = 1, b%ngly
                        do n = 1, b%nglx
                            
                            I = G%intma(n, m, 1, e)
                            
                            hi = b%psiqx(n, iquad) * b%psiqy(m, jquad)
                            
                            zbot(Iq) = zbot(Iq) + zbot_df(I) * hi
                        end do 
                    end do 
                end do 
            end do 
        end do

        zbot_face = 0.0

        !open(10,file='zbot.txt',status='unknown')
        
        do iface=1,G%nface

            !Store Left Side Variables
            ilocl = G%face(5,iface)
            ilocr = G%face(6,iface)
            el = G%face(7,iface)
            er = G%face(8,iface)
                
            do jquad = 1,1

                do iquad = 1, b%nq

                    il = mf%imapl_q(1,iquad,jquad,iface)
                    jl = mf%imapl_q(2,iquad,jquad,iface)
                    kl = mf%imapl_q(3,iquad,jquad,iface)

                    I1 = G%intma_dg_quad(il,jl,kl,el)

                    zbot_face(1,iquad,iface) = zbot(I1)

                    if(er > 0) then

                        ir = mf%imapr_q(1,iquad,jquad,iface)
                        jr = mf%imapr_q(2,iquad,jquad,iface)
                        kr = mf%imapr_q(3,iquad,jquad,iface)

                        I2 = G%intma_dg_quad(ir,jr,kr,er)

                        zbot_face(2,iquad,iface) = zbot(I2)

                    else

                        zbot_face(2,iquad,iface) = zbot_face(1,iquad,iface)

                    end if

                end do
            end do
        end do
        
    end subroutine bot_topo_derivatives
   
    !> Interpolate pressure perturbation from DOF to quadrature points
    subroutine interpolate_pbprime_init(pbprime_df_face, pbprime_df, b, G, inp, mf)

        implicit none

        type(basis), intent(in) :: b
        type(grid), intent(in) :: G
        type(input), intent(in) :: inp
        type(face_CS), intent(in) :: mf

        real, dimension(G%npoin), intent(in) :: pbprime_df
        real, dimension(2,b%ngl,G%nface), intent(out) :: pbprime_df_face
        

        integer :: k, e, iquad, jquad, kquad, l, m, n, I, Iq, iface, I1
        integer :: I2, nq_i, nq_j, plane_ij, ilocl, ilocr, el, er, il, jl, kl, ir, jr, kr
        integer :: ngl_i, ngl_j, ii, jj 
        real :: hi

        pbprime_df_face = 0.0

        do iface = 1, G%nface

            !Store Left Side Variables
            ilocl = G%face(5,iface)
            ilocr = G%face(6,iface)
            el = G%face(7,iface)
            er = G%face(8,iface)

            do jj = 1,1
                do ii = 1, b%ngl

                    il = mf%imapl(1,ii,jj,iface)
                    jl = mf%imapl(2,ii,jj,iface)
                    kl = mf%imapl(3,ii,jj,iface)
                    I1 = G%intma(il,jl,kl,el)

                    pbprime_df_face(1,ii,iface) = pbprime_df(I1)
                    if(er > 0) then

                        ir = mf%imapr(1,ii,jj,iface)
                        jr = mf%imapr(2,ii,jj,iface)
                        kr = mf%imapr(3,ii,jj,iface)
                        I2 = G%intma(ir,jr,kr,er)

                        pbprime_df_face(2,ii,iface) = pbprime_df(I2)
                    else
                        pbprime_df_face(2,ii,iface) = pbprime_df_face(1,ii,iface)
                    endif
                enddo
            enddo
        end do

    end subroutine interpolate_pbprime_init

    !> Wind Stress and Coriolis Force
    subroutine wind_stress_coriolis(b, G, inp, tau_wind, coriolis_df, coriolis_quad, fdt_bcl, fdt2_bcl, &
            a_bcl, b_bcl, tau_wind_df, coriolis_3d_quad, kvector)

        use mod_constants, only: gravity, omega

        implicit none

        type(basis), intent(in) :: b
        type(grid),  intent(in) :: G
        type(input), intent(in) :: inp

        real, dimension(2,G%npoin),   intent(in)  :: tau_wind_df
        real, dimension(3,G%npoin),   intent(in)  :: kvector
        real, dimension(2,G%npoin_q), intent(out) :: tau_wind
        real, dimension(G%npoin),     intent(out) :: coriolis_df
        real, dimension(G%npoin_q),   intent(out) :: coriolis_quad
        real, dimension(G%npoin),     intent(out) :: fdt_bcl, fdt2_bcl, a_bcl, b_bcl
        real, dimension(3,G%npoin_q), intent(out) :: coriolis_3d_quad

        integer :: e, iquad, jquad, kquad, l, m, n, I, Iq
        real    :: ym, Ly, y, hi, x_k, y_k, z_k, f_k

        tau_wind         = 0.0
        coriolis_df      = 0.0
        coriolis_quad    = 0.0
        coriolis_3d_quad = 0.0

        Ly = inp%ydims(2)
        ym = 0.5*Ly

        ! f at DOF nodes: sphere uses 2*Omega*sin(lat) = 2*Omega*kvector(3)
        if (inp%geometry_type == 'sphere_hex' .or. inp%geometry_type == 'sphere_ico') then
            if (trim(inp%test_case) == 'gravity_wave') then
                ! Non-rotating test: f = 0 everywhere.
                do concurrent (I = 1:G%npoin)
                    coriolis_df(I) = 0.0
                end do
            else
                do concurrent (I = 1:G%npoin)
                    coriolis_df(I) = 2.0*omega*kvector(3,I)
                end do
            end if
        else
            do concurrent (I = 1:G%npoin)
                y = G%coord(2,I)
                coriolis_df(I) = inp%f0 + inp%beta*(y - ym)
            end do
        end if

        do concurrent (e = 1:G%nelem, kquad = 1:b%nqz, jquad = 1:b%nqy, iquad = 1:b%nqx)

            Iq = G%intma_dg_quad(iquad, jquad, kquad, e)

            do l = 1, b%nglz
                do m = 1, b%ngly
                    do n = 1, b%nglx

                        I  = G%intma(n, m, l, e)
                        hi = b%psiqx(n, iquad) * b%psiqy(m, jquad)

                        x_k = kvector(1,I)
                        y_k = kvector(2,I)
                        z_k = kvector(3,I)
                        f_k = coriolis_df(I)

                        coriolis_quad(Iq)      = coriolis_quad(Iq)      + f_k        * hi
                        coriolis_3d_quad(1,Iq) = coriolis_3d_quad(1,Iq) + (f_k*x_k) * hi
                        coriolis_3d_quad(2,Iq) = coriolis_3d_quad(2,Iq) + (f_k*y_k) * hi
                        coriolis_3d_quad(3,Iq) = coriolis_3d_quad(3,Iq) + (f_k*z_k) * hi

                        tau_wind(1,Iq) = tau_wind(1,Iq) + tau_wind_df(1,I) * hi
                        tau_wind(2,Iq) = tau_wind(2,Iq) + tau_wind_df(2,I) * hi

                    end do
                end do
            end do
        end do

        fdt_bcl  = inp%dt * coriolis_df
        fdt2_bcl = 0.5 * fdt_bcl
        a_bcl    = 1.0 / (1.0 + fdt2_bcl**2)
        b_bcl    = fdt2_bcl / (1.0 + fdt2_bcl**2)

    end subroutine wind_stress_coriolis

    subroutine ssprk_coefficients(inp, ssprk_a, ssprk_beta, ssprk_a_bcl, ssprk_beta_bcl)
        implicit none

        type(input), intent(in) :: inp
        real, dimension(inp%kstages,3),       intent(out) :: ssprk_a
        real, dimension(inp%kstages),         intent(out) :: ssprk_beta
        real, dimension(inp%kstages_bcl,3),   intent(out) :: ssprk_a_bcl
        real, dimension(inp%kstages_bcl),     intent(out) :: ssprk_beta_bcl

        if(inp%ti_method_btp == 'lsrk') then

            if(inp%kstages == 5) then 

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

                !init%ssprk_a(1,1) = 0.0
                !init%ssprk_a(2,1) = -4.344339134485095 !-567301805773.0 / 1357537059087.0
                !init%ssprk_a(3,1) = 0.0                !-2404267990393.0 / 2016746695238.0
                !init%ssprk_a(4,1) = 3.770024161386381  !-3550918686646.0 / 2091501179385.0
                !init%ssprk_a(5,1) = -0.046347284573284 !-1275806237668.0 / 842570457699.0

                !init%ssprk_beta(1) = 0.713497331193829 !1432997174477.0 / 9575080441755.0
                !init%ssprk_beta(2) = 0.133505249805329 !5161836677717.0 / 13612068292357.0
                !init%ssprk_beta(3) = 0.713497331193829 !1720146321549.0 / 2090206949498.0
                !init%ssprk_beta(4) = 0.149579395628565 !3134564353537.0 / 4481467310338.0
                !init%ssprk_beta(5) = 0.384471116121269 !2277821191437.0 / 14882151754819.0

            elseif(inp%kstages == 14) then

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

            select case (inp%kstages)
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
                ssprk_a(3,1)=0.355909775063327; ssprk_a(3,2)=1.0-ssprk_a(3,1); ssprk_a(3,3)=0.0
                ssprk_beta(3)=0.242995220537396
                ssprk_a(4,1)=0.367933791638137; ssprk_a(4,2)=1.0-ssprk_a(4,1); ssprk_a(4,3)=0.0
                ssprk_beta(4)=0.238458932846290
                ssprk_a(5,1)=0.0; ssprk_a(5,2)=0.762406163401431; ssprk_a(5,3)=1.0-ssprk_a(5,2)
                ssprk_beta(5)=0.287632146308408
            end select
        endif

        select case (inp%kstages_bcl)
        case (2)    ! SSP(2,2)
            ssprk_a_bcl(1,1)=1.0;       ssprk_a_bcl(1,2)=0.0;     ssprk_a_bcl(1,3)=0.0; ssprk_beta_bcl(1)=1.0
            ssprk_a_bcl(2,1)=1.0/2.0;   ssprk_a_bcl(2,2)=1.0/2.0; ssprk_a_bcl(2,3)=0.0; ssprk_beta_bcl(2)=1.0/2.0
        case (3)    ! SSP(3,3)
            ssprk_a_bcl(1,1)=1.0;       ssprk_a_bcl(1,2)=0.0;       ssprk_a_bcl(1,3)=0.0; ssprk_beta_bcl(1)=1.0
            ssprk_a_bcl(2,1)=3.0/4.0;   ssprk_a_bcl(2,2)=1.0/4.0;   ssprk_a_bcl(2,3)=0.0; ssprk_beta_bcl(2)=1.0/4.0
            ssprk_a_bcl(3,1)=1.0/3.0;   ssprk_a_bcl(3,2)=2.0/3.0;   ssprk_a_bcl(3,3)=0.0; ssprk_beta_bcl(3)=2.0/3.0
        end select

    end subroutine ssprk_coefficients

    subroutine poslimiter(b, G, inp, mt, q, alpha)

        use mod_constants, only: gravity

        implicit none

        type(basis),   intent(in)    :: b
        type(grid),    intent(in)    :: G
        type(input),   intent(in)    :: inp
        type(metrics), intent(in)    :: mt

        real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers), intent(inout) :: q
        real, dimension(inp%nlayers),                       intent(in)    :: alpha

        real    :: pmin, pavg, uavg, vavg, wavg, wsum, wjac, threshold, theta, denom
        real    :: d1, d2, d3, d4
        real    :: threshold_k(inp%nlayers)
        integer :: I, k, kc, n, m, e, ipass, nelem_l, nlayers_l, nglx_l, ngly_l, nglz_l, nqx_l, nqy_l
        integer :: nvar_l, npoin_l
        real    :: dry_cutoff_l
        logical :: has_w, has_nan_e

        nelem_l      = G%nelem
        nlayers_l    = inp%nlayers
        nglx_l       = b%nglx
        ngly_l       = b%ngly
        nglz_l       = b%nglz
        nqx_l        = b%nqx
        nqy_l        = b%nqy
        nvar_l       = inp%nvar_bcl
        npoin_l      = G%npoin
        dry_cutoff_l = inp%dry_cutoff
        has_w        = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

        do k = 1, nlayers_l
            threshold_k(k) = (gravity / alpha(k)) * dry_cutoff_l
        end do

        ! Zhang-Shu limiting of each layer toward its element mean, with the
        ! minimum taken over the nodes and the volume and face quadrature points
        ! (where the RHS evaluates the thickness). mt%jac carries the GLL weights,
        ! so the means below are the exact element means.
        !
        ! Pass 1 keeps the column: the change made to layer k (thickness and
        ! momentum) is subtracted from a neighbouring layer kc (the one below; the
        ! one above for the bottom layer), so sum_k q(:,I,k), and with it the
        ! barotropic pb and momentum, is unchanged at every node. The change has
        ! zero element mean, so every layer keeps its element mass and momentum.
        ! Pass 2 is the independent per-layer limiter, a fallback for the rare
        ! case where the compensation pushes a neighbour below its threshold.
        !
        ! DG elements own their nodes exclusively -- no cross-element races on q.
        !$acc parallel loop gang &
        !$acc    copyin(threshold_k) &
        !$acc    present(G%intma, b%psiqx, b%psiqy, mt%jac, q) &
        !$acc    firstprivate(nelem_l, nlayers_l, nglx_l, ngly_l, nglz_l, nqx_l, nqy_l, nvar_l, npoin_l, has_w) &
        !$acc    private(pmin, pavg, uavg, vavg, wavg, wsum, wjac, threshold, theta, denom, &
        !$acc            d1, d2, d3, d4, I, k, kc, n, m, ipass, has_nan_e)
        do e = 1, nelem_l
            !$acc loop seq
            do ipass = 1, 2
                !$acc loop seq
                do k = 1, nlayers_l
                    threshold = threshold_k(k)
                    kc = k + 1
                    if (k == nlayers_l) kc = k - 1

                    pavg      = 0.0
                    uavg      = 0.0
                    vavg      = 0.0
                    wavg      = 0.0
                    wsum      = 0.0
                    has_nan_e = .false.
                    !$acc loop seq
                    do m = 1, ngly_l
                        !$acc loop seq
                        do n = 1, nglx_l
                            I    = G%intma(n,m,1,e)
                            wjac = mt%jac(n,m,1,e)
                            ! IEEE: NaN /= NaN, so track NaN separately.
                            if (q(1,I,k) /= q(1,I,k)) has_nan_e = .true.
                            wsum = wsum + wjac
                            pavg = pavg + wjac * q(1,I,k)
                            uavg = uavg + wjac * q(2,I,k)
                            vavg = vavg + wjac * q(3,I,k)
                            if (has_w) wavg = wavg + wjac * q(4,I,k)
                        end do
                    end do
                    pavg = pavg / wsum
                    uavg = uavg / wsum
                    vavg = vavg / wsum
                    if (has_w) wavg = wavg / wsum

                    pmin = layer_min_e(q, nvar_l, npoin_l, nlayers_l, k, G%intma, nglx_l, ngly_l, nglz_l, &
                                       nelem_l, e, b%psiqx, b%psiqy, nqx_l, nqy_l)

                    if (pavg <= threshold .or. pavg /= pavg .or. has_nan_e) then
                        ! Entire element dry (or NaN) -- clamp to minimum, zero momentum.
                        !$acc loop seq
                        do m = 1, ngly_l
                            !$acc loop seq
                            do n = 1, nglx_l
                                I = G%intma(n,m,1,e)
                                if (ipass == 1 .and. nlayers_l > 1 .and. .not. has_nan_e) then
                                    q(1,I,kc) = q(1,I,kc) + (q(1,I,k) - threshold)
                                    q(2,I,kc) = q(2,I,kc) + q(2,I,k)
                                    q(3,I,kc) = q(3,I,kc) + q(3,I,k)
                                    if (has_w) q(4,I,kc) = q(4,I,kc) + q(4,I,k)
                                end if
                                q(1,I,k) = threshold
                                q(2,I,k) = 0.0
                                q(3,I,k) = 0.0
                                if (has_w) q(4,I,k) = 0.0
                            end do
                        end do
                    else if (pmin < threshold * (1.0 - 1.0e-10)) then
                        ! Zhang-Shu: scale toward element average so min(h) = threshold.
                        denom = pavg - pmin
                        theta = merge(min(1.0, (pavg - threshold) / denom), 0.0, denom /= 0.0)
                        !$acc loop seq
                        do m = 1, ngly_l
                            !$acc loop seq
                            do n = 1, nglx_l
                                I  = G%intma(n,m,1,e)
                                d1 = (theta - 1.0) * (q(1,I,k) - pavg)
                                d2 = (theta - 1.0) * (q(2,I,k) - uavg)
                                d3 = (theta - 1.0) * (q(3,I,k) - vavg)
                                d4 = 0.0
                                if (has_w) d4 = (theta - 1.0) * (q(4,I,k) - wavg)
                                q(1,I,k) = q(1,I,k) + d1
                                q(2,I,k) = q(2,I,k) + d2
                                q(3,I,k) = q(3,I,k) + d3
                                if (has_w) q(4,I,k) = q(4,I,k) + d4
                                if (ipass == 1 .and. nlayers_l > 1) then
                                    q(1,I,kc) = q(1,I,kc) - d1
                                    q(2,I,kc) = q(2,I,kc) - d2
                                    q(3,I,kc) = q(3,I,kc) - d3
                                    if (has_w) q(4,I,kc) = q(4,I,kc) - d4
                                end if
                            end do !n
                        end do !m
                    end if
                end do !k
            end do !ipass
        end do !e
        !$acc end parallel loop

    end subroutine poslimiter

    ! Minimum of layer k's thickness q(1,:,k) over element e: the GLL nodes, the
    ! nqx x nqy volume quadrature points, and the quadrature points on the four
    ! element edges (interpolated from the edge nodes; GLL includes the endpoints).
    function layer_min_e(q, nvar, npoin, nlay, k, intma, nglx, ngly, nglz, nelem, e, &
                         psiqx, psiqy, nqx, nqy) result(pmin)

        !$acc routine seq

        implicit none

        integer, intent(in) :: nvar, npoin, nlay, k, nglx, ngly, nglz, nelem, e, nqx, nqy
        real,    intent(in) :: q(nvar,npoin,nlay)
        integer, intent(in) :: intma(nglx,ngly,nglz,nelem)
        real,    intent(in) :: psiqx(nglx,nqx), psiqy(ngly,nqy)
        real :: pmin

        integer :: n, m, iq, jq
        real    :: s, sx

        pmin = 1.0e20
        do m = 1, ngly
            do n = 1, nglx
                pmin = min(pmin, q(1,intma(n,m,1,e),k))
            end do
        end do
        do jq = 1, nqy
            do iq = 1, nqx
                s = 0.0
                do m = 1, ngly
                    do n = 1, nglx
                        s = s + psiqx(n,iq) * psiqy(m,jq) * q(1,intma(n,m,1,e),k)
                    end do
                end do
                pmin = min(pmin, s)
            end do
        end do
        ! do jq = 1, nqy
        !     s = 0.0; sx = 0.0
        !     do m = 1, ngly
        !         s  = s  + psiqy(m,jq) * q(1,intma(1,m,1,e),k)
        !         sx = sx + psiqy(m,jq) * q(1,intma(nglx,m,1,e),k)
        !     end do
        !     pmin = min(pmin, s, sx)
        ! end do
        ! do iq = 1, nqx
        !     s = 0.0; sx = 0.0
        !     do n = 1, nglx
        !         s  = s  + psiqx(n,iq) * q(1,intma(n,1,1,e),k)
        !         sx = sx + psiqx(n,iq) * q(1,intma(n,ngly,1,e),k)
        !     end do
        !     pmin = min(pmin, s, sx)
        ! end do

    end function layer_min_e

    !------------------------------------------------------------------!
    ! Classify each element per layer:
    !   0 = fully wet, 1 = semi-dry (mixed), 2 = fully dry
    ! Call once per RK stage before the volume and surface RHS routines.
    !------------------------------------------------------------------!
    subroutine find_dry_elements(G, inp, b, tsp, q_df, alpha, dry_flg)

      use mod_constants, only: gravity

      implicit none

      type(grid),      intent(in) :: G
      type(input),     intent(in) :: inp
      type(basis),     intent(in) :: b
      type(tensor_CS), intent(in) :: tsp

      real,    dimension(3,G%npoin,inp%nlayers), intent(in)  :: q_df
      real,    dimension(inp%nlayers),            intent(in)  :: alpha
      integer, dimension(G%nelem,inp%nlayers),   intent(out) :: dry_flg

      integer :: ie, k, ip, I, Iq
      integer :: nelem_l, nlayers_l, npts_l
      real    :: threshold, dp_node, dry_cutoff_l
      logical :: any_wet, any_dry

      nelem_l      = G%nelem
      nlayers_l    = inp%nlayers
      npts_l       = b%npts
      dry_cutoff_l = real(inp%dry_cutoff)

      ! One gang per (element, layer) pair; inner ip loop is sequential.
      !$acc parallel loop gang collapse(2) &
      !$acc    present(tsp%indexq_e, tsp%indexq, q_df, alpha, dry_flg) &
      !$acc    firstprivate(nelem_l, nlayers_l, npts_l, dry_cutoff_l) &
      !$acc    private(any_wet, any_dry, dp_node, threshold, Iq, I)
      do ie = 1, nelem_l
        do k = 1, nlayers_l
          Iq        = tsp%indexq_e(1, ie)
          threshold = (gravity / alpha(k)) * dry_cutoff_l
          any_wet   = .false.
          any_dry   = .false.
          !$acc loop seq
          do ip = 1, npts_l
            I       = tsp%indexq(ip, Iq)
            dp_node = q_df(1, I, k)
            if (dp_node >= threshold) then
              any_wet = .true.
            else
              any_dry = .true.
            end if
          end do
          if (.not. any_dry) then
            dry_flg(ie, k) = 0
          else if (any_wet) then
            dry_flg(ie, k) = 1
          else
            dry_flg(ie, k) = 2
          end if
        end do
      end do
      !$acc end parallel loop

    end subroutine find_dry_elements

    subroutine check_layer_thickness(b, G, inp, q, label, stage)
        ! Scan q(1,:,:) element-by-element for negative or NaN dp.
        ! The first rank that finds a bad element prints the location and aborts.
        ! A device-side count of bad nodes runs first; q is only downloaded and
        ! scanned on the host when that count is non-zero.
        use mpi

        implicit none

        type(basis),  intent(in) :: b
        type(grid),   intent(in) :: G
        type(input),  intent(in) :: inp

        real, dimension(inp%nvar_bcl,G%npoin,inp%nlayers), intent(in) :: q
        character(len=*), intent(in) :: label
        integer,          intent(in) :: stage  ! pass 0 when there is no stage concept

        integer :: k, e, n, m, I, irank, ierr, nn
        integer :: nbad, npoin_l, nlayers_l
        real    :: dp, pmin, cx, cy, cz, radius, clat, clon
        real, parameter :: pi = acos(-1.0), rad2deg = 180.0/acos(-1.0)
        logical :: has_nan

        npoin_l   = G%npoin
        nlayers_l = inp%nlayers

        ! .not.(dp >= 0) is true for both negative and NaN dp.
        nbad = 0
        !$acc parallel loop gang vector collapse(2) present(q) &
        !$acc    firstprivate(npoin_l, nlayers_l) reduction(+:nbad)
        do k = 1, nlayers_l
            do I = 1, npoin_l
                if (.not. (q(1,I,k) >= 0.0)) nbad = nbad + 1
            end do
        end do
        !$acc end parallel loop

        if (nbad == 0) return

        !$acc update host(q(1:1,:,:))

        call mpi_comm_rank(mpi_comm_world, irank, ierr)

        do k = 1, inp%nlayers
            do e = 1, G%nelem
                pmin    = huge(1.0)
                has_nan = .false.
                cx = 0.0; cy = 0.0; cz = 0.0; nn = 0
                do m = 1, b%ngly
                    do n = 1, b%nglx
                        I  = G%intma(n,m,1,e)
                        dp = q(1,I,k)
                        if (dp /= dp) then
                            has_nan = .true.
                        else
                            pmin = min(pmin, dp)
                        end if
                        cx = cx + G%coord(1,I)
                        cy = cy + G%coord(2,I)
                        cz = cz + G%coord(3,I)
                        nn = nn + 1
                    end do
                end do
                if (.not. (pmin >= 0.0) .or. has_nan) then
                    ! Compute element centroid lat/lon from average 3D coordinates.
                    cx = cx / nn;  cy = cy / nn;  cz = cz / nn
                    radius = sqrt(cx**2 + cy**2 + cz**2)
                    clat   = asin(cz / radius) * rad2deg
                    clon   = atan2(cy, cx)     * rad2deg
                    if (clon < 0.0) clon = clon + 360.0

                    if (stage > 0) then
                        if (has_nan) then
                            write(*,'(A,A,": step ",I8,", stage ",I1,", layer ",I3,", element ",I6, &
                                     &", lat=",F7.2,", lon=",F7.2,", NaN in dp, rank ",I4)') &
                                'Fatal error ', trim(label), bcl_itime, stage, k, e, clat, clon, irank
                        else
                            write(*,'(A,A,": step ",I8,", stage ",I1,", layer ",I3,", element ",I6, &
                                     &", lat=",F7.2,", lon=",F7.2,", min dp=",ES12.4,", rank ",I4)') &
                                'Fatal error ', trim(label), bcl_itime, stage, k, e, clat, clon, pmin, irank
                        end if
                    else
                        if (has_nan) then
                            write(*,'(A,A,": step ",I8,", layer ",I3,", element ",I6, &
                                     &", lat=",F7.2,", lon=",F7.2,", NaN in dp, rank ",I4)') &
                                'Fatal error ', trim(label), bcl_itime, k, e, clat, clon, irank
                        else
                            write(*,'(A,A,": step ",I8,", layer ",I3,", element ",I6, &
                                     &", lat=",F7.2,", lon=",F7.2,", min dp=",ES12.4,", rank ",I4)') &
                                'Fatal error ', trim(label), bcl_itime, k, e, clat, clon, pmin, irank
                        end if
                    end if
                    call mpi_abort(mpi_comm_world, 1, ierr)
                end if
            end do
        end do

    end subroutine check_layer_thickness

    subroutine check_btp_thickness(b, G, inp, qb, label)
        ! Scan qb(1,:) element-by-element for negative or NaN barotropic dp.
        ! The first rank that finds a bad element prints the location and aborts.
        ! A device-side count of bad nodes runs first; qb is only downloaded and
        ! scanned on the host when that count is non-zero.
        use mpi

        implicit none

        type(basis),  intent(in) :: b
        type(grid),   intent(in) :: G
        type(input),  intent(in) :: inp

        real, dimension(inp%nvar_btp,G%npoin), intent(in) :: qb
        character(len=*), intent(in) :: label

        integer :: e, n, m, I, irank, ierr, nn
        integer :: nbad, npoin_l
        real    :: dp, pmin, cx, cy, cz, radius, clat, clon
        real, parameter :: rad2deg = 180.0/acos(-1.0)
        logical :: has_nan

        npoin_l = G%npoin

        ! .not.(dp >= 0) is true for both negative and NaN dp.
        ! Callers may have queued work on async queue 1 (barotropic sub-step).
        !$acc wait(1)
        nbad = 0
        !$acc parallel loop gang vector present(qb) &
        !$acc    firstprivate(npoin_l) reduction(+:nbad)
        do I = 1, npoin_l
            if (.not. (qb(1,I) >= 0.0)) nbad = nbad + 1
        end do
        !$acc end parallel loop

        if (nbad == 0) return

        !$acc update host(qb(1:1,:))

        call mpi_comm_rank(mpi_comm_world, irank, ierr)

        do e = 1, G%nelem
            pmin    = huge(1.0)
            has_nan = .false.
            cx = 0.0; cy = 0.0; cz = 0.0; nn = 0
            do m = 1, b%ngly
                do n = 1, b%nglx
                    I  = G%intma(n,m,1,e)
                    dp = qb(1,I)
                    if (dp /= dp) then
                        has_nan = .true.
                    else
                        pmin = min(pmin, dp)
                    end if
                    cx = cx + G%coord(1,I)
                    cy = cy + G%coord(2,I)
                    cz = cz + G%coord(3,I)
                    nn = nn + 1
                end do
            end do
            if (.not. (pmin >= 0.0) .or. has_nan) then
                ! Compute element centroid lat/lon from average 3D coordinates.
                cx = cx / nn;  cy = cy / nn;  cz = cz / nn
                radius = sqrt(cx**2 + cy**2 + cz**2)
                clat   = asin(cz / radius) * rad2deg
                clon   = atan2(cy, cx)     * rad2deg
                if (clon < 0.0) clon = clon + 360.0

                if (has_nan) then
                    write(*,'(A,A,": step ",I8,", element ",I6,", lat=",F7.2,", lon=",F7.2,", NaN in barotropic dp, rank ",I4)') &
                        'Fatal error ', trim(label), bcl_itime, e, clat, clon, irank
                else
                    write(*,'(A,A,": step ",I8,", element ",I6,", lat=",F7.2,", lon=",F7.2,", min barotropic dp=",ES12.4,", rank ",I4)') &
                        'Fatal error ', trim(label), bcl_itime, e, clat, clon, pmin, irank
                end if
                call mpi_abort(mpi_comm_world, 1, ierr)
            end if
        end do

    end subroutine check_btp_thickness

    function safe_div(n, d, altv) result(q)
        implicit none
        real, intent(in) :: n, d, altv
        real :: q
        if (exponent(n) - exponent(d) >= maxexponent(n) .or. d == 0.0) then
            q = altv
        else
            q = n / d
        end if
    end function safe_div

    ! Debug-only: copy per-element reductions, not the full state, after GPU scans.
    ! Layer quantities are geometric thickness [m]; BTP quantities are pressure [Pa].
    subroutine trace_reset()
        if (trace_open) close(trace_unit)
        trace_open = .false.
        trace_dumped = .false.
        trace_ring_count = 0
        trace_sequence = 0
        trace_before_step = -1
        trace_before_stage = -1
        if (allocated(trace_before)) deallocate(trace_before)
    end subroutine trace_reset

    subroutine trace_stage(b, G, inp, tsp, alpha, q, qb, phase, btp_only, substep, stage)
        use mpi
        use mod_constants, only: gravity
        implicit none
        type(basis), intent(in) :: b
        type(grid), intent(in) :: G
        type(input), intent(in) :: inp
        type(tensor_CS), intent(in) :: tsp
        real, intent(in) :: alpha(inp%nlayers)
        real, intent(in) :: q(inp%nvar_bcl,G%npoin,inp%nlayers), qb(inp%nvar_btp,G%npoin)
        character(len=*), intent(in) :: phase
        logical, intent(in) :: btp_only
        integer, intent(in) :: substep, stage
        ! Metrics: node min, volume min, face min, mean, integral, pressure
        ! residual, momentum residual, nonfinite count, node index, volume
        ! quadrature index, face sample index. Each element has its own output.
        real, allocatable :: stats(:,:,:)
        real :: h, v, w, vq, vs, weight, mass, area, residual, mom, denom, floor_k
        real :: xyz(3), radius, longitude, latitude, value, delta_mass
        real :: pi, hi
        integer :: e,k,j,n,m,ip,iqx,iqy,iq,side,I,loc,ne,nx,ny,nqx,nqy,nl,nv,nvb
        integer :: rank,ierr,metric,emin,changed,trigger_e,trigger_k,iu,c,slot
        integer :: nonfinite
        logical :: event, paired
        character(len=256) :: filename
        character(len=2048) :: line
        character(len=32), parameter :: names(8) = [character(len=32) :: &
            'node_min','volume_min','face_min','element_mean_min','integral', &
            'pressure_residual_max','momentum_residual_max','nonfinite_samples']
        if (.not. inp%ltrace_stages) return
        if (b%nglz /= 1 .or. b%nglx /= b%ngly) &
            error stop 'ltrace_stages requires equal-order horizontal quadrilateral elements'
        call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
        if (.not. trace_open) then
            write(filename,'("stage_trace_rank",i6.6,".csv")') rank
            open(newunit=trace_unit,file=trim(filename),status='replace',action='write')
            write(trace_unit,'(a)') 'sequence,step,bcl_stage,btp_substep,btp_stage,phase,rank,layer,' // &
                'metric,value,local_element,sample,longitude_deg,latitude_deg'
            trace_open = .true.
        end if
        ne=G%nelem; nx=b%nglx; ny=b%ngly; nqx=b%nqx; nqy=b%nqy
        if (ne == 0) return
        nl=merge(1,inp%nlayers,btp_only); nv=inp%nvar_bcl; nvb=inp%nvar_btp
        allocate(stats(11,ne,nl))
        pi=acos(-1.0)
        !$acc wait(1)
        !$acc parallel loop gang collapse(2) copyout(stats) &
        !$acc present(q,qb,G%intma,b%psiqx,b%psiqy,b%wglx,b%wgly,tsp%wjac_df,alpha) &
        !$acc firstprivate(ne,nx,ny,nqx,nqy,nl,nv,nvb,btp_only) &
        !$acc private(I,h,v,w,vq,vs,weight,mass,area,residual,mom,denom,nonfinite,loc,ip,side,iq,iqx,iqy,j,n,m)
        do k=1,nl
            do e=1,ne
                stats(:,e,k)=0.0
                stats(1:4,e,k)=huge(1.0)
                mass=0.0; area=0.0; nonfinite=0
                !$acc loop seq
                do m=1,ny
                    !$acc loop seq
                    do n=1,nx
                        I=G%intma(n,m,1,e)
                        if (btp_only) then
                            h=qb(1,I)
                        else
                            h=q(1,I,k)*alpha(k)/gravity
                        end if
                        if (.not.(abs(h)<=huge(h))) then
                            nonfinite=nonfinite+1
                        else if (h<stats(1,e,k)) then
                            stats(1,e,k)=h;stats(9,e,k)=real(n+(m-1)*nx)
                        end if
                        weight=tsp%wjac_df(I)
                        mass=mass+weight*h;area=area+weight
                        if (btp_only) then
                            !$acc loop seq
                            do j=2,nvb
                                if (.not.(abs(qb(j,I))<=huge(h))) nonfinite=nonfinite+1
                            end do
                        else
                            v=0.0
                            !$acc loop seq
                            do j=1,nl
                                v=v+q(1,I,j)
                            end do
                            denom=max(abs(qb(1,I)),tiny(1.0))
                            residual=abs(v-qb(1,I))/denom
                            stats(6,e,k)=max(stats(6,e,k),residual)
                            !$acc loop seq
                            do j=2,nv
                                v=0.0;w=0.0
                                !$acc loop seq
                                do ip=1,nl
                                    v=v+q(j,I,ip);w=w+abs(q(j,I,ip))
                                end do
                                mom=abs(v-qb(j+1,I))/max(abs(qb(j+1,I)),w,tiny(1.0))
                                stats(7,e,k)=max(stats(7,e,k),mom)
                                if (.not.(abs(v)<=huge(v))) nonfinite=nonfinite+1
                            end do
                            if (.not.(abs(residual)<=huge(residual))) nonfinite=nonfinite+1
                        end if
                    end do
                end do
                stats(4,e,k)=mass/area;stats(5,e,k)=mass
                !$acc loop seq
                do iqy=1,nqy
                    !$acc loop seq
                    do iqx=1,nqx
                        vq=0.0
                        !$acc loop seq
                        do m=1,ny
                            !$acc loop seq
                            do n=1,nx
                                I=G%intma(n,m,1,e)
                                if (btp_only) then
                                    h=qb(1,I)
                                else
                                    h=q(1,I,k)*alpha(k)/gravity
                                end if
                                vq=vq+b%psiqx(n,iqx)*b%psiqy(m,iqy)*h
                            end do
                        end do
                        if (.not.(abs(vq)<=huge(vq))) then
                            nonfinite=nonfinite+1
                        else if(vq<stats(2,e,k)) then
                            stats(2,e,k)=vq;stats(10,e,k)=real(iqx+(iqy-1)*nqx)
                        end if
                    end do
                end do
                ! Four edges, using the same 1-D interpolation tables as the volume rule.
                !$acc loop seq
                do side=1,4
                    !$acc loop seq
                    do iq=1,merge(nqy,nqx,side<=2)
                        vs=0.0
                        !$acc loop seq
                        do ip=1,merge(ny,nx,side<=2)
                            if(side<=2) then
                                n=merge(1,nx,side==1);m=ip;weight=b%psiqy(ip,iq)
                            else
                                n=ip;m=merge(1,ny,side==3);weight=b%psiqx(ip,iq)
                            end if
                            I=G%intma(n,m,1,e)
                            if(btp_only) then
                                h=qb(1,I)
                            else
                                h=q(1,I,k)*alpha(k)/gravity
                            end if
                            vs=vs+weight*h
                        end do
                        if (.not.(abs(vs)<=huge(vs))) then
                            nonfinite=nonfinite+1
                        else if(vs<stats(3,e,k)) then
                            stats(3,e,k)=vs;stats(11,e,k)=real(iq+(side-1)*max(nqx,nqy))
                        end if
                    end do
                end do
                stats(8,e,k)=real(nonfinite)
            end do
        end do
        !$acc end parallel loop
        event=.false.;trigger_e=0;trigger_k=0
        paired=.false.
        if(.not.btp_only) paired=allocated(trace_before) .and. &
            trace_before_step==bcl_itime .and. trace_before_stage==trace_bcl_stage
        do k=1,nl
            floor_k=inp%dry_cutoff
            if(btp_only) floor_k=0.0
            do metric=1,8
                if(metric<=4) then
                    emin=max(1,minloc(stats(metric,:,k),dim=1));value=stats(metric,emin,k)
                else if(metric==5.or.metric==8) then
                    emin=maxloc(stats(metric,:,k),dim=1);value=sum(stats(metric,:,k))
                else
                    emin=max(1,maxloc(stats(metric,:,k),dim=1));value=stats(metric,emin,k)
                end if
                emin=max(1,min(ne,emin))
                loc=0;xyz=0.0
                if(metric==1) then
                    loc=nint(stats(9,emin,k));n=mod(max(loc,1)-1,nx)+1;m=(max(loc,1)-1)/nx+1
                    xyz=G%coord(:,G%intma(n,m,1,emin))
                else if(metric==2) then
                    loc=nint(stats(10,emin,k));iqx=mod(max(loc,1)-1,nqx)+1;iqy=(max(loc,1)-1)/nqx+1
                    do m=1,ny
                        do n=1,nx
                            xyz=xyz+b%psiqx(n,iqx)*b%psiqy(m,iqy)*G%coord(:,G%intma(n,m,1,emin))
                        end do
                    end do
                else if(metric==3) then
                    loc=nint(stats(11,emin,k));side=(max(loc,1)-1)/max(nqx,nqy)+1
                    iq=mod(max(loc,1)-1,max(nqx,nqy))+1
                    do ip=1,merge(ny,nx,side<=2)
                        if(side<=2) then
                            n=merge(1,nx,side==1);m=ip;hi=b%psiqy(ip,iq)
                        else
                            n=ip;m=merge(1,ny,side==3);hi=b%psiqx(ip,iq)
                        end if
                        xyz=xyz+hi*G%coord(:,G%intma(n,m,1,emin))
                    end do
                else
                    do m=1,ny
                        do n=1,nx
                            xyz=xyz+G%coord(:,G%intma(n,m,1,emin))/real(nx*ny)
                        end do
                    end do
                end if
                radius=sqrt(sum(xyz**2));longitude=0.0;latitude=0.0
                if(radius>0.0) then
                    longitude=modulo(atan2(xyz(2),xyz(1))*180.0/pi,360.0)
                    latitude=asin(max(-1.0,min(1.0,xyz(3)/radius)))*180.0/pi
                end if
                ! BTP has no meaningful layer residuals; suppress those rows.
                if(btp_only.and.(metric==6.or.metric==7)) cycle
                write(line,'(5(i0,","),a,",",2(i0,","),a,",",es24.16e3,",",2(i0,","),es16.8,",",es16.8)') &
                    trace_sequence,bcl_itime,trace_bcl_stage,substep,stage,trim(phase),rank, &
                    merge(0,k,btp_only),trim(names(metric)),value,emin,loc,longitude,latitude
                call trace_line(line)
                if((metric<=3.and.value<=floor_k).or.(metric==8.and.value>0.0)) then
                    event=.true.;trigger_e=emin;trigger_k=k
                end if
                if(metric==6.and.index(phase,'after_velocity')>0.and.value>1.e-6) then
                    event=.true.;trigger_e=emin;trigger_k=k
                end if
            end do
            if(paired.and.index(phase,'after_limit')>0) then
                delta_mass=sum(stats(5,:,k))-sum(trace_before(1,:,k))
                changed=count(abs(stats(5,:,k)-trace_before(1,:,k)) > &
                    100.0*epsilon(1.0)*max(abs(trace_before(1,:,k)),tiny(1.0)))
                write(line,'(5(i0,","),a,",",2(i0,","),a,",",es24.16e3,",0,0,0,0")') &
                    trace_sequence,bcl_itime,trace_bcl_stage,substep,stage,trim(phase),rank,k,'limiter_mass_change',delta_mass
                call trace_line(line)
                write(line,'(5(i0,","),a,",",2(i0,","),a,",",es24.16e3,",0,0,0,0")') &
                    trace_sequence,bcl_itime,trace_bcl_stage,substep,stage,trim(phase),rank,k, &
                    'limiter_mean_changed_elements',real(changed)
                call trace_line(line)
                changed=count(trace_before(2,:,k)<=inp%dry_cutoff)
                write(line,'(5(i0,","),a,",",2(i0,","),a,",",es24.16e3,",0,0,0,0")') &
                    trace_sequence,bcl_itime,trace_bcl_stage,substep,stage,trim(phase),rank,k, &
                    'limiter_candidate_elements',real(changed)
                call trace_line(line)
                ! Candidate count is determined from pre-limiter nodal minima.

            end if
        end do
        if(.not.btp_only.and.index(phase,'before_limit')>0) then
            if(.not.allocated(trace_before)) allocate(trace_before(2,ne,inp%nlayers))
            trace_before(1,:,:)=stats(5,:,:)
            trace_before(2,:,:)=stats(1,:,:)
            trace_before_step=bcl_itime;trace_before_stage=trace_bcl_stage
        end if
        flush(trace_unit)
        if(event.and..not.trace_dumped) then
            write(filename,'("stage_trace_context_rank",i6.6,".csv")') rank
            open(newunit=iu,file=trim(filename),status='replace',action='write')
            do c=max(1,trace_ring_count-127),trace_ring_count
                slot=mod(c-1,128)+1;write(iu,'(a)') trim(trace_ring(slot))
            end do
            close(iu)
            ! Full state download only once, on this rank's first event.
            !$acc update host(q,qb)
            write(filename,'("stage_trace_event_rank",i6.6,".txt")') rank
            open(newunit=iu,file=trim(filename),status='replace',action='write')
            write(iu,*) 'step, BCL stage, BTP substep/stage, phase, BTP-only:', &
                bcl_itime,trace_bcl_stage,substep,stage,trim(phase),btp_only
            write(iu,*) 'local trigger element/layer:',trigger_e,trigger_k
            write(iu,*) 'BTP event q fields are frozen qprime, not layer prognostic state'
            do e=1,ne
                ! Include the trigger and face-neighbor elements.
                if(e/=trigger_e) then
                    changed=0
                    do j=1,G%nface
                        if ((G%face(7,j)==trigger_e.and.G%face(8,j)==e).or. &
                            (G%face(8,j)==trigger_e.and.G%face(7,j)==e)) changed=1
                    end do
                    if(changed==0) cycle
                end if
                do m=1,ny
                    do n=1,nx
                        I=G%intma(n,m,1,e)
                        write(iu,*) 'element,node,xyz,qb:',e,I,G%coord(:,I),qb(:,I)
                        do k=1,inp%nlayers
                            write(iu,*) 'layer,q:',k,q(:,I,k)
                        end do
                    end do
                end do
            end do
            close(iu);trace_dumped=.true.
        end if
        trace_sequence=trace_sequence+1
    end subroutine trace_stage

    subroutine trace_line(line)
        character(len=*), intent(in) :: line
        write(trace_unit,'(a)') trim(line)
        trace_ring_count=trace_ring_count+1
        trace_ring(mod(trace_ring_count-1,128)+1)=line
    end subroutine trace_line

end module mod_initial_mlswe
