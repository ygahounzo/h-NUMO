!--------------------------------------------------------------------------
!@brief This module contains the interpolation from quad to nodal points.
!>@ author by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: July 02, 2023
!-------------------------------------------------------------------------
module mod_Tensorproduct

    use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz
    use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad, mod_grid_get_face_nq
    use mod_basis, only: nq, psiq, psiqx, psiqy, psiqz, dpsiq, dpsiqx, dpsiqy, dpsiqz
    use mod_metrics, only: ksiq_x,ksiq_y,ksiq_z, etaq_x,etaq_y,etaq_z, zetaq_x,zetaq_y,zetaq_z

    public :: compute_gradient_quad,interpolate_layer_from_quad_to_node_1d

    private

    contains


    ! subroutine compute_gradient_quad(grad_q_quad,q)

    !     use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, npts
    !     use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad
    !     use mod_basis, only: nq, psiq, psiqx, psiqy, psiqz, dpsiq, dpsiqx, dpsiqy, dpsiqz
    !     use mod_metrics, only: ksiq_x,ksiq_y,ksiq_z, etaq_x,etaq_y,etaq_z, zetaq_x,zetaq_y,zetaq_z

    !     use mod_initial, only: dpsidx,dpsidy, indexq

    !     implicit none

  
    !     real, dimension(npoin), intent(in) :: q
    !     real, dimension(2,npoin_q), intent(out) :: grad_q_quad
    !     integer :: e, iquad, jquad, kquad, l, m, n, Iq, I, ip
    !     real :: e_x, e_y, n_x, n_y, h_e, h_n,h_c, c_x, c_y, c_z, e_z, n_z, dhdx, dhdy
        
    !     grad_q_quad = 0.0d0
        
    !     !Construct Volume Integral Contributions

    !     do Iq = 1, npoin_q

    !         do ip = 1,npts

    !             I = indexq(Iq,ip)

    !             grad_q_quad(1,Iq) = grad_q_quad(1,Iq) + dpsidx(Iq,ip)*q(I)
    !             grad_q_quad(2,Iq) = grad_q_quad(2,Iq) + dpsidy(Iq,ip)*q(I)

    !         end do
    !     end do

    ! end subroutine compute_gradient_quad

    subroutine compute_gradient_quad(grad_q_quad,q)

        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz
        use mod_grid, only:  nelem, npoin, npoin_q, intma, intma_dg_quad
        use mod_basis, only: nq, psiq, psiqx, psiqy, psiqz, dpsiq, dpsiqx, dpsiqy, dpsiqz
        use mod_metrics, only: ksiq_x,ksiq_y,ksiq_z, etaq_x,etaq_y,etaq_z, zetaq_x,zetaq_y,zetaq_z

        implicit none

  
        real, dimension(npoin), intent(in) :: q
        real, dimension(2,npoin_q), intent(out) :: grad_q_quad
        integer :: e, iquad, jquad, kquad, l, m, n, Iq, I
        real :: e_x, e_y, n_x, n_y, h_e, h_n,h_c, c_x, c_y, c_z, e_z, n_z
        
        grad_q_quad = 0.0
        
        !Construct Volume Integral Contribution
        do e = 1, nelem
            
            !Loop Integration Points
            !do kquad = 1, nqz
                do jquad = 1, nqy
                    do iquad = 1, nqx
                        
        
                        e_x = ksiq_x(iquad,jquad,1,e); e_y = ksiq_y(iquad,jquad,1,e)
                        n_x = etaq_x(iquad,jquad,1,e); n_y = etaq_y(iquad,jquad,1,e)
                        
                        Iq = intma_dg_quad(iquad,jquad,1,e)
        
                        !Interpolate at Integration Points
                        !do l = 1, nglz
                            do m = 1, ngly
                                do n = 1, nglx
                                    
                                    I = intma(n,m,1,e)
                
                                    !Xi derivatives
                                    h_e = dpsiqx(n,iquad)*psiqy(m,jquad)
                                    !Eta derivatives
                                    h_n = psiqx(n,iquad)*dpsiqy(m,jquad)
                
                                    grad_q_quad(1,Iq) = grad_q_quad(1,Iq) + (h_e*e_x + h_n*n_x)*q(I)
                                    grad_q_quad(2,Iq) = grad_q_quad(2,Iq) + (h_e*e_y + h_n*n_y)*q(I)

                                end do
                            end do
                        !end do

                    end do !iquad
                end do !jquad
            !end do !kquad
        end do !e
        
    end subroutine compute_gradient_quad

    subroutine interpolate_layer_from_quad_to_node_1d(q_df,q)

        use mod_grid, only : npoin_q, npoin, nelem, intma_dg_quad, intma
        use mod_basis, only: npts, psiqx, psiqy
        use mod_metrics, only: massinv, jacq
        use mod_input, only: nlayers
        use mod_basis, only: nglx, ngly, nglz, nqx, nqy, nqz, psiqx, psiqy, psiqz, nq, psiq
        !use mod_initial, only: psih, dpsidx,dpsidy, indexq, wjac

        implicit none
        
        real, dimension(npoin), intent(out) :: q_df
        real, dimension(npoin_q), intent(in) :: q

        real :: wq, var_u, var_v, hi
        integer :: iquad, jquad, kquad, e, k, m, n, l, I, Iq, ip

        q_df = 0.0

        do e = 1, nelem
            do kquad = 1, nqz
                do jquad = 1, nqy
                    do iquad = 1, nqx
                    
                        Iq = intma_dg_quad(iquad, jquad, kquad, e)

                        wq = jacq(iquad,jquad,kquad,e)
                        var_u = q(Iq)
                        
                        do l = 1, nglz
                            do m = 1, ngly
                                do n = 1, nglx
                                    
                                    I = intma(n, m, l, e)
                                    
                                    hi = psiqx(n, iquad) * psiqy(m, jquad)
                                    q_df(I) = q_df(I) + wq*var_u*hi

                                end do
                            end do 
                        end do

                    end do 
                end do 
            end do
        end do

        q_df(:) = massinv(:)*q_df(:)
    
    end subroutine interpolate_layer_from_quad_to_node_1d

end module mod_Tensorproduct