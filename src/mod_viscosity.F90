!----------------------------------------------------------------------!
!>@brief This module contains the subroutines that create 
!> what is needed for the (hyper)-viscosity terms
!>@author Francis X. Giraldo
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@author Simone Marras
!> @date November 10 2014 simply changed kmatrix=0 to kmatrix=0.0d0
!>
!> @date November 11 2014, F. X. Giraldo: added mod_viscosity_create_Nazarov
!>
!----------------------------------------------------------------------!
module mod_viscosity
  
    use mod_basis, only: nglx, ngly, nglz

    use mod_grid, only: npoin, nelem, intma, npoin_q
  
    use mod_initial, only: nvar, kvector

    use mod_input, only: vertical_viscosity, visc, lvisc_anisotropic, &
         diff_T, diff_Th, diff_Tv, visc_h, visc_v, diff_S, diff_Sh, diff_Sv,&
         lsalinity, is_mlswe, visc_mlswe

    use mod_metrics, only: ksi_x,  ksi_y,  ksi_z, &
        eta_x,  eta_y,  eta_z, &
        zeta_x, zeta_y, zeta_z, jac

    use mod_parallel, only: num_send_recv_total

    public :: &
        mod_viscosity_create, &
        q_visc, &
        rhs_visc, &
        metrics_visc, &
        grad_q, &
        grad_q_mlswe, &   ! added by Yao Gahounzo
        visc_elem, &
        res_bq, &
        nvar_visc, nmetrics_visc, &
        recv_data_visc, &
        visc_pass, &
        q_visc_mlswe, &
        rhs_visc_mlswe
  
    private
    !-----------------------------------------------------------------------
  
    !module variables and parameters
    real, dimension(:,:),       allocatable :: q_visc, rhs_visc, recv_data_visc, q_visc_mlswe, rhs_visc_mlswe
    real, dimension(:,:,:,:,:), allocatable :: metrics_visc
    real, dimension(:,:,:),     allocatable :: grad_q, grad_q_mlswe
    real, dimension(:,:),       allocatable :: visc_elem
    real, dimension(:,:),       allocatable :: res_bq
    real, dimension(:,:),       allocatable :: visc_pass
    integer :: nvar_visc, nmetrics_visc
  !-----------------------------------------------------------------------
  
contains
  
    !-----------------------------------------------------------------------
    subroutine mod_viscosity_create()

        implicit none

        integer :: AllocateStatus
        integer :: e, i, j, k, l, ip
        real :: e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, wq, flag
        real :: a_e(3), a_n(3), a_c(3), g1, g2, g3, g4, g5, g6
        real, dimension(3,3) :: kmatrix

        !check to see how many viscous arrays we need per gridpoint.
        nvar_visc=nvar
        if(is_mlswe) then
            nvar_visc=2
        end if
        nmetrics_visc=6
        if (lvisc_anisotropic) then
            nvar_visc=2*nvar
            nmetrics_visc=12
        end if
        
        if (allocated(q_visc)) then
           deallocate(q_visc, rhs_visc, metrics_visc, grad_q, visc_elem, res_bq, recv_data_visc, visc_pass,grad_q_mlswe, q_visc_mlswe, rhs_visc_mlswe)
        endif
         
        allocate( q_visc(nvar_visc,npoin), rhs_visc(nvar_visc,npoin), &
            metrics_visc(nmetrics_visc,nglx,ngly,nglz,nelem), &
            !grad_q(3,nvar,npoin), &  ! commented by Yao G.
            grad_q(3,nvar_visc,npoin), & ! Yao G. added this line
            grad_q_mlswe(3,nvar_visc,npoin), &   ! added by Yao Gahounzo
            visc_elem(nvar_visc,npoin), &
            res_bq(nelem,nvar), &
            recv_data_visc(nvar_visc,num_send_recv_total), &
            visc_pass(nvar,2), &
            q_visc_mlswe(nvar_visc,npoin), &
            rhs_visc_mlswe(nvar_visc,npoin), &
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - MOD_VISCOSITY_CREATE **"

        !Initialize
        grad_q=0
        grad_q_mlswe = 0.0  ! added by Yao Gahounzo
        q_visc=0
        rhs_visc=0
        res_bq=0.0
        metrics_visc=0
        visc_elem=visc
        if(is_mlswe) then
            visc_elem = visc_mlswe
        end if
        visc_pass=visc
        q_visc_mlswe=0.0
        rhs_visc_mlswe=0.0

        !Fill visc_pass array for PISO
        visc_pass=visc
        visc_pass(5,:)=diff_T
        if(lsalinity) visc_pass(6,:)=diff_S
        if(lvisc_anisotropic) then
           visc_pass(2:4,1) = visc_h
           visc_pass(2:4,2) = visc_v
           visc_pass(5,1) = diff_Th
           visc_pass(5,2) = diff_Tv
           visc_elem(1:nvar,:)=visc_h
           visc_elem(nvar+1:nvar_visc,:)=visc_v

           if(lsalinity) then
              visc_pass(6,1) = diff_Sh
              visc_pass(6,2) = diff_Sv
           end if
        end if

        !Set Vertical_Viscosity Flag
        flag=vertical_viscosity !FXG:NEW_VISCOSITY

        !Form Viscosity Metric Terms
        do e=1,nelem
            do k=1,nglz
                do j=1,ngly
                    do i=1,nglx
                        wq=jac(i,j,k,e)

                        ip=intma(i,j,k,e)
                        e_x=ksi_x(i,j,k,e);  e_y=ksi_y(i,j,k,e);  e_z=ksi_z(i,j,k,e)
                        n_x=eta_x(i,j,k,e);  n_y=eta_y(i,j,k,e);  n_z=eta_z(i,j,k,e)
                        c_x=zeta_x(i,j,k,e); c_y=zeta_y(i,j,k,e); c_z=zeta_z(i,j,k,e)

                        !Kmatrix or Hmatrix (depending on if not separating or separating directions)
                        call mod_viscosity_create_hmatrix(kmatrix,kvector(:,ip),flag)
                        do l=1,3
                            a_e(l)=kmatrix(l,1)*e_x + kmatrix(l,2)*e_y + kmatrix(l,3)*e_z
                            a_n(l)=kmatrix(l,1)*n_x + kmatrix(l,2)*n_y + kmatrix(l,3)*n_z
                            a_c(l)=kmatrix(l,1)*c_x + kmatrix(l,2)*c_y + kmatrix(l,3)*c_z
                        end do

                        g1=a_e(1)*a_e(1) + a_e(2)*a_e(2) + a_e(3)*a_e(3)
                        g2=a_e(1)*a_n(1) + a_e(2)*a_n(2) + a_e(3)*a_n(3)
                        g3=a_e(1)*a_c(1) + a_e(2)*a_c(2) + a_e(3)*a_c(3)
                        g4=a_n(1)*a_n(1) + a_n(2)*a_n(2) + a_n(3)*a_n(3)
                        g5=a_n(1)*a_c(1) + a_n(2)*a_c(2) + a_n(3)*a_c(3)
                        g6=a_c(1)*a_c(1) + a_c(2)*a_c(2) + a_c(3)*a_c(3)

                        !Store Metrics for Diffusion Terms
                        metrics_visc(1,i,j,k,e)=wq*g1
                        metrics_visc(2,i,j,k,e)=wq*g2
                        metrics_visc(3,i,j,k,e)=wq*g3
                        metrics_visc(4,i,j,k,e)=wq*g4
                        metrics_visc(5,i,j,k,e)=wq*g5
                        metrics_visc(6,i,j,k,e)=wq*g6

                        if (lvisc_anisotropic) then
                            !Vmatrix (for vertical direction)
                            call mod_viscosity_create_vmatrix(kmatrix,kvector(:,ip),flag)
                            do l=1,3
                                a_e(l)=kmatrix(l,1)*e_x + kmatrix(l,2)*e_y + kmatrix(l,3)*e_z
                                a_n(l)=kmatrix(l,1)*n_x + kmatrix(l,2)*n_y + kmatrix(l,3)*n_z
                                a_c(l)=kmatrix(l,1)*c_x + kmatrix(l,2)*c_y + kmatrix(l,3)*c_z
                            end do

                            g1=a_e(1)*a_e(1) + a_e(2)*a_e(2) + a_e(3)*a_e(3)
                            g2=a_e(1)*a_n(1) + a_e(2)*a_n(2) + a_e(3)*a_n(3)
                            g3=a_e(1)*a_c(1) + a_e(2)*a_c(2) + a_e(3)*a_c(3)
                            g4=a_n(1)*a_n(1) + a_n(2)*a_n(2) + a_n(3)*a_n(3)
                            g5=a_n(1)*a_c(1) + a_n(2)*a_c(2) + a_n(3)*a_c(3)
                            g6=a_c(1)*a_c(1) + a_c(2)*a_c(2) + a_c(3)*a_c(3)

                            !Store Metrics for Diffusion Terms
                            metrics_visc(1+6,i,j,k,e)=wq*g1
                            metrics_visc(2+6,i,j,k,e)=wq*g2
                            metrics_visc(3+6,i,j,k,e)=wq*g3
                            metrics_visc(4+6,i,j,k,e)=wq*g4
                            metrics_visc(5+6,i,j,k,e)=wq*g5
                            metrics_visc(6+6,i,j,k,e)=wq*g6
                        end if !lvisc_anisotropic

                    end do !i
                end do !j
            end do !k
        end do !e

    end subroutine mod_viscosity_create


    !-----------------------------------------------------------------------
    subroutine mod_viscosity_create_hmatrix(kmatrix,kvector,flag)
    
        implicit none
    
        !global arrays
        real, intent(out)             :: kmatrix(3,3)
        real, intent(in)              :: kvector(3), flag

        !local arrays
        integer :: i,j

        kmatrix=0
        if (flag /= 0d0) then
            do i=1,3
                kmatrix(i,i)=1.0d0
            end do
        else if (flag == 0d0) then
            kmatrix(1,1)=1.0d0 - kvector(1)*kvector(1)
            kmatrix(1,2)=  - kvector(1)*kvector(2)
            kmatrix(1,3)=  - kvector(1)*kvector(3)
            kmatrix(2,1)=  kmatrix(1,2)
            kmatrix(2,2)=1.0d0 - kvector(2)*kvector(2)
            kmatrix(2,3)=  - kvector(2)*kvector(3)
            kmatrix(3,1)=  kmatrix(1,3)
            kmatrix(3,2)=  kmatrix(2,3)
            kmatrix(3,3)=1.0d0 - kvector(3)*kvector(3)
        end if

    end subroutine mod_viscosity_create_hmatrix

    !-----------------------------------------------------------------------
    subroutine mod_viscosity_create_vmatrix(kmatrix,kvector,flag)
    
        implicit none
    
        !global arrays
        real, intent(out)             :: kmatrix(3,3)
        real, intent(in)              :: kvector(3), flag

        !local arrays
        integer :: i,j

        kmatrix=0
        if (flag /= 0d0) then
            do i=1,3
                kmatrix(i,i)=0d0
            end do
        else if (flag == 0d0) then
            kmatrix(1,1)=kvector(1)*kvector(1)
            kmatrix(1,2)=kvector(1)*kvector(2)
            kmatrix(1,3)=kvector(1)*kvector(3)
            kmatrix(2,1)=kmatrix(1,2)
            kmatrix(2,2)=kvector(2)*kvector(2)
            kmatrix(2,3)=kvector(2)*kvector(3)
            kmatrix(3,1)=kmatrix(1,3)
            kmatrix(3,2)=kmatrix(2,3)
            kmatrix(3,3)=kvector(3)*kvector(3)
        end if

    end subroutine mod_viscosity_create_vmatrix

    !-----------------------------------------------------------------------
    subroutine mod_viscosity_create_kmatrix(kmatrix,kvector,flag)
    
        implicit none
    
        !global arrays
        real, intent(out)             :: kmatrix(3,3)
        real, intent(in)              :: kvector(3), flag

        !local arrays
        integer :: i,j

        kmatrix=0
        if (flag /= 0d0) then
            do i=1,3
                kmatrix(i,i)=1.0d0
            end do
        else if (flag == 0d0) then
            kmatrix(1,1)=1.0d0 - kvector(1)*kvector(1)
            kmatrix(1,2)=  - kvector(1)*kvector(2)
            kmatrix(1,3)=  - kvector(1)*kvector(3)
            kmatrix(2,1)=  kmatrix(1,2)
            kmatrix(2,2)=1.0d0 - kvector(2)*kvector(2)
            kmatrix(2,3)=  - kvector(2)*kvector(3)
            kmatrix(3,1)=  kmatrix(1,3)
            kmatrix(3,2)=  kmatrix(2,3)
            kmatrix(3,3)=1.0d0 - kvector(3)*kvector(3)
        end if

    end subroutine mod_viscosity_create_kmatrix

!-----------------------------------------------------------------------
!  subroutine mod_viscosity_create_v0()
!    
!    implicit none
!    
!    integer :: AllocateStatus
!    integer :: e, i, j, k
!    real :: e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z, wq, flag
!
!    allocate( q_visc(nvar,npoin), rhs_visc(nvar,npoin), &
!              metrics_visc(6,nglx,ngly,nglz,nelem), &
!              stat=AllocateStatus )
!    if (AllocateStatus /= 0) stop "** Not Enough Memory - MOD_VISCOSITY_CREATE **"
!
!    !Set Vertical_Viscosity Flag
!    flag=vertical_viscosity !FXG: BIC  (only valid for CUBE Geometry)
!
!    !Form Viscosity Metric Terms
!    do e=1,nelem
!       do k=1,nglz
!          do j=1,ngly
!             do i=1,nglx
!                wq=jac(i,j,k,e)
!                e_x=ksi_x(i,j,k,e);  e_y=ksi_y(i,j,k,e);  e_z=ksi_z(i,j,k,e)
!                n_x=eta_x(i,j,k,e);  n_y=eta_y(i,j,k,e);  n_z=eta_z(i,j,k,e)
!                c_x=zeta_x(i,j,k,e); c_y=zeta_y(i,j,k,e); c_z=zeta_z(i,j,k,e)
!
!                metrics_visc(1,i,j,k,e)=wq*( e_x*e_x + e_y*e_y + flag*e_z*e_z )
!                metrics_visc(2,i,j,k,e)=wq*( e_x*n_x + e_y*n_y + flag*e_z*n_z )
!                metrics_visc(3,i,j,k,e)=wq*( e_x*c_x + e_y*c_y + flag*e_z*c_z )
!                metrics_visc(4,i,j,k,e)=wq*( n_x*n_x + n_y*n_y + flag*n_z*n_z )
!                metrics_visc(5,i,j,k,e)=wq*( n_x*c_x + n_y*c_y + flag*n_z*c_z )
!                metrics_visc(6,i,j,k,e)=wq*( c_x*c_x + c_y*c_y + flag*c_z*c_z )
!             end do !i
!          end do !j
!       end do !k
!    end do !e
!
!  end subroutine mod_viscosity_create_v0



!-----------------------------------------------------------------------
!  subroutine mod_viscosity_create_kmatrix_v0(kmatrix,kvector,flag)
!    
!    implicit none
!    
!    !global arrays
!    real, intent(out)             :: kmatrix(3,3)
!    real, intent(in)              :: kvector(3), flag
!
!    !local arrays
!    integer :: i,j
!
!    kmatrix=0
!    if (flag /= 0d0) then
!       do i=1,3
!          kmatrix(i,i)=1.0d0
!       end do
!    else if (flag == 0d0) then
!       do i=1,3
!          do j=1,3
!             kmatrix(i,j)=-kvector(i)*kvector(j)
!          end do !j
!          kmatrix(i,i)=kmatrix(i,i) + 1.0d0
!       end do
!    end if
!
!  end subroutine mod_viscosity_create_kmatrix_v0

end module mod_viscosity
