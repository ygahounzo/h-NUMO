! ===============================================================================================
! This module contains the routines for the baroclinic flux terms
!   Author: Yao Gahounzo
!   Computing PhD
!   Boise State University
!   Date: March 27, 2023
! ===============================================================================================

module mod_layer_terms

    use mod_grid,      only: grid
    use mod_basis,     only: basis
    use mod_input,     only: input
    use mod_initial,   only: initial
    use mod_face,      only: face_CS
    use mod_tensor,    only: tensor_CS
    use mod_variables, only: btp_CS, bcl_CS

    implicit none

    public :: layer_mom_boundary_df,    &
              velocity_df,              &
              evaluate_consistency_face, &
              extract_qprime_df_face,   &
              extract_dprime_df_face,   &
              interpolate_dpp

contains

    subroutine interpolate_dpp(G, inp, b, tsp, dprimeq, dprime_df)

        ! Interpolate layer mass from dofs to quad points

        implicit none

        type(grid),      intent(in)    :: G
        type(input),     intent(in)    :: inp
        type(basis),     intent(in)    :: b
        type(tensor_CS), intent(in)    :: tsp

        real, dimension(G%npoin_q, inp%nlayers), intent(inout) :: dprimeq
        real, dimension(G%npoin,   inp%nlayers), intent(inout) :: dprime_df

        integer :: Iq, I, ip
        real    :: hi

        dprimeq = 0.0

        do Iq = 1, G%npoin_q
            do ip = 1, b%npts
                I  = tsp%indexq(ip, Iq)
                hi = tsp%psih(ip, Iq)
                dprimeq(Iq,:) = dprimeq(Iq,:) + dprime_df(I,:) * hi
            end do
        end do

    end subroutine interpolate_dpp

    subroutine evaluate_consistency_face(G, inp, b, mf, btp, bcl, init, &
                                         mass_deficit_mass_face, dprime_df)

        ! Computes the mass deficit face values for consistency correction

        implicit none

        type(grid),    intent(in) :: G
        type(input),   intent(in) :: inp
        type(basis),   intent(in) :: b
        type(face_CS), intent(in) :: mf
        type(btp_CS),  intent(in) :: btp
        type(bcl_CS),  intent(in) :: bcl
        type(initial), intent(in) :: init

        real, intent(out) :: mass_deficit_mass_face(2, 2, b%nq, G%nface, inp%nlayers)
        real, intent(in)  :: dprime_df(G%npoin, inp%nlayers)

        integer :: iface, iquad, k, il, jl, kl, el, er, ir, jr, kr, I, n
        real    :: qprime_l, qprime_r, weights_face_l, weights_face_r, hi
        real    :: pbprime_l, pbprime_r

        mass_deficit_mass_face = 0.0

        do k = 1, inp%nlayers

            do iface = 1, G%nface

                el = G%face(7,iface)
                er = G%face(8,iface)

                do iquad = 1, b%nq

                    qprime_l = 0.0;  qprime_r = 0.0
                    pbprime_l = 0.0; pbprime_r = 0.0

                    do n = 1, b%ngl
                        hi = b%psiq(n,iquad)

                        il = mf%imapl(1,n,1,iface)
                        jl = mf%imapl(2,n,1,iface)
                        kl = mf%imapl(3,n,1,iface)

                        I = G%intma(il,jl,kl,el)
                        qprime_l  = qprime_l  + hi * dprime_df(I,k)
                        pbprime_l = pbprime_l + hi * init%pbprime_df(I)
                    enddo

                    if(er > 0) then
                        do n = 1, b%ngl
                            hi = b%psiq(n,iquad)
                            ir = mf%imapr(1,n,1,iface)
                            jr = mf%imapr(2,n,1,iface)
                            kr = mf%imapr(3,n,1,iface)
                            I = G%intma(ir,jr,kr,er)
                            qprime_r  = qprime_r  + hi * dprime_df(I,k)
                            pbprime_r = pbprime_r + hi * init%pbprime_df(I)
                        enddo
                    else
                        qprime_r  = qprime_l
                        pbprime_r = pbprime_l
                    endif

                    weights_face_l = qprime_l / pbprime_l
                    weights_face_r = qprime_r / pbprime_r

                    mass_deficit_mass_face(1,1,iquad,iface,k) = weights_face_l * &
                        (btp%btp_mass_flux_face_ave(1,iquad,iface) &
                         - bcl%sum_layer_mass_flux_face(1,iquad,iface))
                    mass_deficit_mass_face(2,1,iquad,iface,k) = weights_face_l * &
                        (btp%btp_mass_flux_face_ave(2,iquad,iface) &
                         - bcl%sum_layer_mass_flux_face(2,iquad,iface))

                    mass_deficit_mass_face(1,2,iquad,iface,k) = weights_face_r * &
                        (btp%btp_mass_flux_face_ave(1,iquad,iface) &
                         - bcl%sum_layer_mass_flux_face(1,iquad,iface))
                    mass_deficit_mass_face(2,2,iquad,iface,k) = weights_face_r * &
                        (btp%btp_mass_flux_face_ave(2,iquad,iface) &
                         - bcl%sum_layer_mass_flux_face(2,iquad,iface))
                end do
            end do
        end do

        call bcl_create_communicator(mass_deficit_mass_face, 2, inp%nlayers, b%nq)

    end subroutine evaluate_consistency_face

    subroutine velocity_df(G, inp, q_df, qb_df)

        implicit none

        type(grid),  intent(in)    :: G
        type(input), intent(in)    :: inp

        real, dimension(3, G%npoin, inp%nlayers), intent(inout) :: q_df
        real, dimension(4, G%npoin),              intent(in)    :: qb_df

        real    :: ubar, vbar
        integer :: I, k
        real    :: uv_df(2, G%npoin, inp%nlayers)

        uv_df = 0.0

        do k = 1, inp%nlayers
            uv_df(1,:,k) = q_df(2,:,k) / q_df(1,:,k)
            uv_df(2,:,k) = q_df(3,:,k) / q_df(1,:,k)
        end do

        do I = 1, G%npoin
            ubar = 0.0; vbar = 0.0

            do k = 1, inp%nlayers
                ubar = ubar + uv_df(1,I,k) * q_df(1,I,k)
                vbar = vbar + uv_df(2,I,k) * q_df(1,I,k)
            end do

            if(qb_df(1,I) > 0.0) then
                ubar = ubar / qb_df(1,I)
                vbar = vbar / qb_df(1,I)

                do k = 1, inp%nlayers
                    uv_df(1,I,k) = uv_df(1,I,k) - ubar + qb_df(3,I)/qb_df(1,I)
                    uv_df(2,I,k) = uv_df(2,I,k) - vbar + qb_df(4,I)/qb_df(1,I)
                end do
            else
                uv_df(:,I,:) = 0.0
            end if
        end do

        do k = 1, inp%nlayers
            q_df(2,:,k) = uv_df(1,:,k) * q_df(1,:,k)
            q_df(3,:,k) = uv_df(2,:,k) * q_df(1,:,k)
        end do

    end subroutine velocity_df

    subroutine extract_velocity(G, inp, uv_df, q_df, qb_df)

        implicit none

        type(grid),  intent(in)  :: G
        type(input), intent(in)  :: inp

        real, dimension(2, G%npoin, inp%nlayers), intent(out) :: uv_df
        real, dimension(3, G%npoin, inp%nlayers), intent(in)  :: q_df
        real, dimension(4, G%npoin),              intent(in)  :: qb_df

        real    :: ubar, vbar
        integer :: I, k

        uv_df = 0.0

        do k = 1, inp%nlayers
            uv_df(1,:,k) = q_df(2,:,k) / q_df(1,:,k)
            uv_df(2,:,k) = q_df(3,:,k) / q_df(1,:,k)
        end do

        do I = 1, G%npoin
            ubar = 0.0
            vbar = 0.0

            do k = 1, inp%nlayers
                ubar = ubar + (uv_df(1,I,k) * q_df(1,I,k))
                vbar = vbar + (uv_df(2,I,k) * q_df(1,I,k))
            end do

            if(qb_df(1,I) > 0.0) then
                ubar = ubar / qb_df(1,I)
                vbar = vbar / qb_df(1,I)

                do k = 1, inp%nlayers
                    uv_df(1,I,k) = uv_df(1,I,k) - (ubar - qb_df(3,I)/qb_df(1,I))
                    uv_df(2,I,k) = uv_df(2,I,k) - (vbar - qb_df(4,I)/qb_df(1,I))
                end do
            else
                uv_df(:,I,:) = 0.0
            end if
        end do

    end subroutine extract_velocity

    subroutine extract_qprime_df_face(G, inp, init, qprime_df, q_df, qb_df)

        implicit none

        type(grid),    intent(in) :: G
        type(input),   intent(in) :: inp
        type(initial), intent(in) :: init

        real, dimension(3, G%npoin, inp%nlayers), intent(out) :: qprime_df
        real, dimension(3, G%npoin, inp%nlayers), intent(in)  :: q_df
        real, dimension(4, G%npoin),              intent(in)  :: qb_df

        integer :: k, I
        real    :: ope
        real    :: uv_df(2, G%npoin, inp%nlayers)

        qprime_df = 0.0

        call extract_velocity(G, inp, uv_df, q_df, qb_df)

        do k = 1, inp%nlayers
            do I = 1, G%npoin
                ope = sum(q_df(1,I,:)) / init%pbprime_df(I)
                qprime_df(1,I,k) = q_df(1,I,k) / ope
                qprime_df(2,I,k) = uv_df(1,I,k) - qb_df(3,I)/qb_df(1,I)
                qprime_df(3,I,k) = uv_df(2,I,k) - qb_df(4,I)/qb_df(1,I)
            end do
        end do

    end subroutine extract_qprime_df_face

    subroutine extract_dprime_df_face(G, inp, b, mf, dprime_df_face, dprime_df)

        ! Extracts face values of layer mass at DG nodes

        implicit none

        type(grid),    intent(in) :: G
        type(input),   intent(in) :: inp
        type(basis),   intent(in) :: b
        type(face_CS), intent(in) :: mf

        real, dimension(2, b%ngl, G%nface, inp%nlayers), intent(out) :: dprime_df_face
        real, dimension(G%npoin, inp%nlayers),            intent(in)  :: dprime_df

        integer :: iface, el, er, il, jl, ir, jr, kl, kr, I, n

        dprime_df_face = 0.0

        do iface = 1, G%nface

            el = G%face(7,iface)
            er = G%face(8,iface)

            do n = 1, b%ngl

                il = mf%imapl(1,n,1,iface)
                jl = mf%imapl(2,n,1,iface)
                kl = mf%imapl(3,n,1,iface)
                I  = G%intma(il,jl,kl,el)

                dprime_df_face(1,n,iface,:) = dprime_df(I,:)

                if(er > 0) then
                    ir = mf%imapr(1,n,1,iface)
                    jr = mf%imapr(2,n,1,iface)
                    kr = mf%imapr(3,n,1,iface)
                    I  = G%intma(ir,jr,kr,er)
                    dprime_df_face(2,n,iface,:) = dprime_df(I,:)
                else
                    dprime_df_face(2,n,iface,:) = dprime_df_face(1,n,iface,:)
                end if
            end do
        end do

    end subroutine extract_dprime_df_face

    subroutine layer_mom_boundary_df(G, inp, b, mf, q)

        implicit none

        type(grid),    intent(in)    :: G
        type(input),   intent(in)    :: inp
        type(basis),   intent(in)    :: b
        type(face_CS), intent(in)    :: mf

        real, intent(inout) :: q(3, G%npoin, inp%nlayers)

        integer :: iface, n, il, jl, el, er, ilocl, ilocr, I, kl, k
        real    :: nx, ny, upnl(inp%nlayers)

        do iface = 1, G%nface

            ilocl = G%face(5,iface)
            ilocr = G%face(6,iface)
            el    = G%face(7,iface)
            er    = G%face(8,iface)

            if(er == -4) then
                do n = 1, b%ngl

                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)

                    nx = mf%normal_vector(1,n,1,iface)
                    ny = mf%normal_vector(2,n,1,iface)

                    upnl     = q(2,I,:)*nx + q(3,I,:)*ny
                    q(2,I,:) = q(2,I,:) - upnl*nx
                    q(3,I,:) = q(3,I,:) - upnl*ny

                end do

            elseif(er == -2) then
                do n = 1, b%ngl

                    il = mf%imapl(1,n,1,iface)
                    jl = mf%imapl(2,n,1,iface)
                    kl = mf%imapl(3,n,1,iface)
                    I  = G%intma(il,jl,kl,el)

                    q(2,I,:) = 0.0
                    q(3,I,:) = 0.0

                end do
            end if
        end do

    end subroutine layer_mom_boundary_df

end module mod_layer_terms
