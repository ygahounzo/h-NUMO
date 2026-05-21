!----------------------------------------------------------------------!
!>@brief This module contains Initial Conditions settings
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: March 27, 2023
!----------------------------------------------------------------------!
module mod_initial_mlswe

    use mod_basis,   only: basis
    use mod_grid,    only: grid
    use mod_metrics, only: metrics
    use mod_input,   only: input
    use mod_face,    only: face_CS

    public :: &
        bot_topo_derivatives, &
        interpolate_pbprime_init, wind_stress_coriolis, &
        map_deriv, ssprk_coefficients

    private


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
    subroutine wind_stress_coriolis(b, G, inp, tau_wind,coriolis_df,coriolis_quad, fdt_bcl, fdt2_bcl, &
            a_bcl, b_bcl,tau_wind_df)

        use mod_constants, only: gravity
    
        implicit none

        type(basis), intent(in) :: b
        type(grid), intent(in) :: G
        type(input), intent(in) :: inp
        
        real, dimension(2,G%npoin), intent(in) :: tau_wind_df
        real, dimension(2,G%npoin_q), intent(out) :: tau_wind
        real, dimension(G%npoin), intent(out) :: coriolis_df
        real, dimension(G%npoin_q), intent(out) :: coriolis_quad
        real, dimension(G%npoin), intent(out) :: fdt_bcl, fdt2_bcl, a_bcl, b_bcl

        integer :: k, e, iquad, jquad, kquad, l, m, n, I, Iq, ip
        real :: ym, Ly, y, hi, tau0, lat, sig, rho_air, w

        tau_wind = 0.0
        coriolis_df = 0.0
        coriolis_quad = 0.0

        gravity = 9.806
        
        Ly = inp%ydims(2)
        ym = 0.5*Ly

        do concurrent(I = 1:G%npoin)
            y = G%coord(2,I)
            coriolis_df(I) = inp%f0 + inp%beta*(y - ym)
        end do


        do concurrent (e = 1:G%nelem, kquad = 1:b%nqz, jquad = 1:b%nqy, iquad = 1:b%nqx)
                    
            Iq = G%intma_dg_quad(iquad, jquad, kquad, e)
            
            do l = 1, b%nglz
                do m = 1, b%ngly
                    do n = 1, b%nglx
                        
                        I = G%intma(n, m, l, e)
                        
                        hi = b%psiqx(n, iquad) * b%psiqy(m, jquad)
                        
                        coriolis_quad(Iq) = coriolis_quad(Iq) + coriolis_df(I) * hi
                        
                        tau_wind(1,Iq) = tau_wind(1,Iq) + tau_wind_df(1,I) * hi
                        tau_wind(2,Iq) = tau_wind(2,Iq) + tau_wind_df(2,I) * hi
                        
                    end do
                end do
            end do
        end do
        
        fdt_bcl = inp%dt * coriolis_df
        fdt2_bcl = 0.5 * fdt_bcl
        a_bcl = 1.0 / (1.0 + fdt2_bcl**2)
        b_bcl = fdt2_bcl / (1.0 + fdt2_bcl**2)

    end subroutine wind_stress_coriolis

    subroutine ssprk_coefficients(inp, ssprk_a,ssprk_beta)
        implicit none

        type(input), intent(in) :: inp
        real, dimension(inp%kstages,3), intent(out) :: ssprk_a
        real, dimension(inp%kstages), intent(out) :: ssprk_beta

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
