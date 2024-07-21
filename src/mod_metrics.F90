!----------------------------------------------------------------------!
!>@brief This module builds the metrics: Ksi_X and Jacobian
!>@author  Francis X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>Upgraded to Parallel by James F. Kelly
!>15 April 2010 
!----------------------------------------------------------------------!
module mod_metrics

    use mod_basis, only: nglx, ngly, nglz, nqx,nqy,nqz, npts

    use mod_grid, only: npoin_cg, npoin, nelem, ncol, nz, nz_cg, &
        do_1dIMEX, npoin_massinv

    use mod_input, only: nelz, is_non_conforming_flg, space_method, is_mlswe

    use mod_parallel, only: num_nbh, num_send_recv_total

    use mod_mpi_utilities, only: irank

    public :: &
        mod_metrics_create_metrics, &
        mod_metrics_create_mass, &
        ksi_x, ksi_y, ksi_z, &
        eta_x, eta_y, eta_z, &
        zeta_x, zeta_y, zeta_z, &
        jac, xjac, mass, massinv, &
        jac_1d, mass_1d, massinv_1d, &
        ksiq_x,ksiq_y,ksiq_z, &
        etaq_x,etaq_y,etaq_z, &
        zetaq_x,zetaq_y,zetaq_z, jacq, xjacq

    private
    !-----------------------------------------------------------------------
    !module variables and parameters
    real, dimension(:,:,:,:), allocatable :: ksi_x,ksi_y,ksi_z, &
        eta_x,eta_y,eta_z, zeta_x,zeta_y,zeta_z, jac, xjac
    real, dimension(:,:,:,:), allocatable :: ksiq_x,ksiq_y,ksiq_z, &
        etaq_x,etaq_y,etaq_z, zetaq_x,zetaq_y,zetaq_z, jacq, xjacq
    real, dimension(:),       allocatable :: mass, massinv
    real, dimension(:,:),     allocatable :: mass_1d, massinv_1d
    real, dimension(:,:,:),   allocatable :: jac_1d
  !-----------------------------------------------------------------------

contains

    !-----------------------------------------------------------------------

    subroutine mod_metrics_create_metrics()

        implicit none

        integer :: AllocateStatus
        integer :: i, j

        if(allocated(ksi_x)) deallocate(ksi_x,ksi_y,ksi_z,eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,jac,xjac)

        allocate( ksi_x(nglx,ngly,nglz,nelem),ksi_y(nglx,ngly,nglz,nelem),ksi_z(nglx,ngly,nglz,nelem), &
            eta_x(nglx,ngly,nglz,nelem),eta_y(nglx,ngly,nglz,nelem),eta_z(nglx,ngly,nglz,nelem), &
            zeta_x(nglx,ngly,nglz,nelem), zeta_y(nglx,ngly,nglz,nelem), zeta_z(nglx,ngly,nglz,nelem), &
            jac(nglx,ngly,nglz,nelem), xjac(nglx,ngly,nglz,nelem),stat = AllocateStatus )

        if(is_mlswe) then ! quads data points added by Yao Gahounzo
            allocate( ksiq_x(nqx,nqy,nqz,nelem),ksiq_y(nqx,nqy,nqz,nelem),ksiq_z(nqx,nqy,nqz,nelem), & 
            etaq_x(nqx,nqy,nqz,nelem),etaq_y(nqx,nqy,nqz,nelem),etaq_z(nqx,nqy,nqz,nelem), &
            zetaq_x(nqx,nqy,nqz,nelem), zetaq_y(nqx,nqy,nqz,nelem), zetaq_z(nqx,nqy,nqz,nelem), &
            jacq(nqx,nqy,nqz,nelem), xjacq(nqx,nqy,nqz,nelem), &
            stat = AllocateStatus )
        end if  ! is_mlswe
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Metrics_Create_Metrics **"

        !Generate Volume Metric Terms
        call metrics(ksi_x,ksi_y,ksi_z,eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,jac,xjac)
        if(is_mlswe) then ! quads data points added by Yao Gahounzo
            call metrics_quad(ksiq_x,ksiq_y,ksiq_z,etaq_x,etaq_y,etaq_z,zetaq_x,zetaq_y,zetaq_z,jacq,xjacq)
        end if  ! is_mlswe

    end subroutine mod_metrics_create_metrics

    !-----------------------------------------------------------------------
    !>@brief Creates the (diag.) mass matrix and inverse mass matrix using a two step process:
    !> 1. Perform a local DSS (on processor)
    !> 2. Perform a global DSS (communicating across adjacent processors)
    !-----------------------------------------------------------------------
    subroutine mod_metrics_create_mass

        implicit none

        integer inbh, ib, ip, ip_g, ic, i, j, k, e
        integer :: AllocateStatus
        real, dimension(:), allocatable :: mass_recv
        real, dimension(:,:), allocatable :: mass_1d_continuous
        real, dimension(:), allocatable :: mass_cg

        if(allocated(mass)) deallocate(mass,massinv)
        allocate( mass(npoin), massinv(npoin_massinv),                         &
          mass_recv(num_send_recv_total), stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Metrics_Create_Mass **"

        !Compute Mass Matrix (on each element)
        call create_mass(mass,jac)

        !Compute Inverse Mass Matrix: This is what is used in the CODE

        !Do DSS on mass matrix
        call create_global_rhs(mass,mass_recv,1,0)
        massinv = 1.0/mass


        !Deallocate
        deallocate(mass_recv)

    end subroutine mod_metrics_create_mass

end module mod_metrics

