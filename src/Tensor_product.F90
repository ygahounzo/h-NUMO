subroutine Tensor_product(wjac,psih,dpsidx,dpsidy,indexq)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma

        use mod_basis, only: nglx, ngly, nglz, npts, dpsiqx, dpsiqy, dpsiqz, nqx, nqy, nqz, &
            psiqx, psiqy, psiqz

        use mod_metrics, only: &
            ksiq_x, ksiq_y, ksiq_z, &
            etaq_x, etaq_y, etaq_z, &
            zetaq_x, zetaq_y, zetaq_z, &
            jacq

        use mod_constants, only: gravity

        !use mod_initial, only: grad_zbot_quad, tau_wind

        implicit none 

        real, dimension(npoin_q,npts), intent(out) :: psih, dpsidx,dpsidy
        integer, dimension(npoin_q,npts) :: indexq
        real, dimension(npoin_q), intent(out) :: wjac

        integer :: e, jquad, iquad, Iq, ip, m, n, I
        real :: h_e, h_n
        real :: e_x, e_y, n_x, n_y
      
        wjac = 0.0
        psih = 0.0
        dpsidx = 0.0
        dpsidy = 0.0
        indexq = 0

        do e = 1,nelem

            do jquad = 1,nqy
                do iquad = 1,nqx
                    
                    Iq = intma_dg_quad(iquad,jquad,1,e)

                    wjac(Iq) = jacq(iquad,jquad,1,e)

                    e_x = ksiq_x(iquad,jquad,1,e); e_y = ksiq_y(iquad,jquad,1,e); 
                    n_x = etaq_x(iquad,jquad,1,e); n_y = etaq_y(iquad,jquad,1,e);

                    ip = 0;
                    
                    do m = 1, ngly 
                        do n = 1, nglx 

                            I = intma(n,m,1,e)
                            ip = ip + 1
                
                            indexq(Iq, ip) = I
                            psih(Iq, ip) = psiqx(n, iquad) * psiqy(m, jquad)
                
                            ! Xi derivatives
                            h_e = dpsiqx(n, iquad) * psiqy(m, jquad)

                            ! Eta derivatives
                            h_n = psiqx(n, iquad) * dpsiqy(m, jquad)
                
                            ! Pressure terms
                            dpsidx(Iq, ip) = h_e * e_x + h_n * n_x
                            dpsidy(Iq, ip) = h_e * e_y + h_n * n_y

                        end do !n
                    end do !m

                end do !iquad
            end do !jquad

        end do !e

    end subroutine Tensor_product