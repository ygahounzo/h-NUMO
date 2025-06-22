subroutine Tensor_product(wjac,psih,dpsidx,dpsidy,indexq, wjac_df,psih_df,dpsidx_df, &
                        dpsidy_df,index_df)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma

        use mod_basis, only: nglx, ngly, nglz, npts, dpsiqx, dpsiqy, dpsiqz, nqx, nqy, nqz, &
            psiqx, psiqy, psiqz

        use mod_basis, only: dpsix, dpsiy, dpsiz, psix, psiy, psiz

        use mod_metrics, only: &
            ksiq_x, ksiq_y, ksiq_z, &
            etaq_x, etaq_y, etaq_z, &
            zetaq_x, zetaq_y, zetaq_z, &
            jacq, &
            ksi_x, ksi_y, ksi_z, &
            eta_x, eta_y, eta_z, &
            jac

        use mod_constants, only: gravity

        !use mod_initial, only: grad_zbot_quad, tau_wind

        implicit none

        real, dimension(npts,npoin_q), intent(out) :: psih, dpsidx,dpsidy
        integer, dimension(npts,npoin_q) :: indexq
        real, dimension(npoin_q), intent(out) :: wjac

        real, dimension(npts,npoin), intent(out) :: psih_df, dpsidx_df,dpsidy_df
        integer, dimension(npts,npoin) :: index_df
        real, dimension(npoin), intent(out) :: wjac_df

        integer :: e, jquad, iquad, Iq, ip, m, n, I
        real :: h_e, h_n
        real :: e_x, e_y, n_x, n_y

        wjac = 0.0
        psih = 0.0
        dpsidx = 0.0
        dpsidy = 0.0
        indexq = 0

        wjac_df = 0.0
        psih_df = 0.0
        dpsidx_df = 0.0
        dpsidy_df = 0.0
        index_df = 0

        do e = 1,nelem

            do jquad = 1,nqy
                do iquad = 1,nqx

                    Iq = intma_dg_quad(iquad,jquad,1,e)

                    wjac(Iq) = jacq(iquad,jquad,1,e)

                    e_x = ksiq_x(iquad,jquad,1,e); e_y = ksiq_y(iquad,jquad,1,e);
                    n_x = etaq_x(iquad,jquad,1,e); n_y = etaq_y(iquad,jquad,1,e);

                    ip = 0

                    do m = 1, ngly
                        do n = 1, nglx

                            I = intma(n,m,1,e)
                            ip = ip + 1

                            indexq(ip,Iq) = I
                            psih(ip,Iq) = psiqx(n, iquad) * psiqy(m, jquad)

                            ! Xi derivatives
                            h_e = dpsiqx(n, iquad) * psiqy(m, jquad)

                            ! Eta derivatives
                            h_n = psiqx(n, iquad) * dpsiqy(m, jquad)

                            ! Pressure terms
                            dpsidx(ip,Iq) = h_e * e_x + h_n * n_x
                            dpsidy(ip,Iq) = h_e * e_y + h_n * n_y

                        end do !n
                    end do !m

                end do !iquad
            end do !jquad

            do jquad = 1,ngly
                do iquad = 1,nglx

                    Iq = intma(iquad,jquad,1,e)

                    wjac_df(Iq) = jac(iquad,jquad,1,e)

                    e_x = ksi_x(iquad,jquad,1,e); e_y = ksi_y(iquad,jquad,1,e);
                    n_x = eta_x(iquad,jquad,1,e); n_y = eta_y(iquad,jquad,1,e);

                    ip = 0

                    do m = 1, ngly
                        do n = 1, nglx

                            I = intma(n,m,1,e)
                            ip = ip + 1

                            index_df(ip,Iq) = I
                            psih_df(ip,Iq) = psix(n, iquad) * psiy(m, jquad)

                            ! Xi derivatives
                            h_e = dpsix(n, iquad) * psiy(m, jquad)

                            ! Eta derivatives
                            h_n = psix(n, iquad) * dpsiy(m, jquad)

                            ! Pressure terms
                            dpsidx_df(ip,Iq) = h_e * e_x + h_n * n_x
                            dpsidy_df(ip,Iq) = h_e * e_y + h_n * n_y

                        end do !n
                    end do !m

                end do !iquad
            end do !jquad

        end do !e

end subroutine Tensor_product
