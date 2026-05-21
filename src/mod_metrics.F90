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

    use mod_basis
    use mod_grid
    use mod_input
    use mod_parallel

    type metrics
        ! module variables and parameters
        real, dimension(:,:,:,:), allocatable :: ksi_x,ksi_y,ksi_z, &
            eta_x,eta_y,eta_z, zeta_x,zeta_y,zeta_z, jac, xjac
        real, dimension(:,:,:,:), allocatable :: ksiq_x,ksiq_y,ksiq_z, &
            etaq_x,etaq_y,etaq_z, zetaq_x,zetaq_y,zetaq_z, jacq, xjacq
        real, dimension(:),       allocatable :: mass, massinv
    end type metrics    

    public :: &
        mod_metrics_create_metrics, &
        mod_metrics_create_mass, metrics

contains

    !-----------------------------------------------------------------------

    subroutine mod_metrics_create_metrics(G, b, inp, mt)

        implicit none

        type(grid), intent(in) :: G
        type(basis), intent(in) :: b
        type(input), intent(in) :: inp
        type(metrics), intent(inout) :: mt

        integer :: AllocateStatus
        integer :: i, j

        if(allocated(mt%ksi_x)) deallocate(mt%ksi_x,mt%ksi_y,mt%ksi_z,mt%eta_x,mt%eta_y,mt%eta_z,mt%zeta_x, &
                                mt%zeta_y,mt%zeta_z,mt%jac,mt%xjac)

        allocate( mt%ksi_x(b%nglx,b%ngly,b%nglz,G%nelem),mt%ksi_y(b%nglx,b%ngly,b%nglz,G%nelem), &
            mt%ksi_z(b%nglx,b%ngly,b%nglz,G%nelem), mt%eta_x(b%nglx,b%ngly,b%nglz,G%nelem),mt%eta_y(b%nglx,b%ngly,b%nglz,G%nelem), &
            mt%eta_z(b%nglx,b%ngly,b%nglz,G%nelem), mt%zeta_x(b%nglx,b%ngly,b%nglz,G%nelem), &
            mt%zeta_y(b%nglx,b%ngly,b%nglz,G%nelem), mt%zeta_z(b%nglx,b%ngly,b%nglz,G%nelem), &
            mt%jac(b%nglx,b%ngly,b%nglz,G%nelem), mt%xjac(b%nglx,b%ngly,b%nglz,G%nelem),stat = AllocateStatus )

        if(inp%is_mlswe) then ! quads data points added by Yao Gahounzo
            allocate( mt%ksiq_x(b%nqx,b%nqy,b%nqz,G%nelem),mt%ksiq_y(b%nqx,b%nqy,b%nqz,G%nelem), &
            mt%ksiq_z(b%nqx,b%nqy,b%nqz,G%nelem), mt%etaq_x(b%nqx,b%nqy,b%nqz,G%nelem),mt%etaq_y(b%nqx,b%nqy,b%nqz,G%nelem), &
            mt%etaq_z(b%nqx,b%nqy,b%nqz,G%nelem), mt%zetaq_x(b%nqx,b%nqy,b%nqz,G%nelem), mt%zetaq_y(b%nqx,b%nqy,b%nqz,G%nelem), &
            mt%zetaq_z(b%nqx,b%nqy,b%nqz,G%nelem), mt%jacq(b%nqx,b%nqy,b%nqz,G%nelem), mt%xjacq(b%nqx,b%nqy,b%nqz,G%nelem), &
            stat = AllocateStatus )
        end if  ! is_mlswe
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Metrics_Create_Metrics **"

        !Generate Volume Metric Terms
        call compute_metrics(mt%ksi_x,mt%ksi_y,mt%ksi_z,mt%eta_x,mt%eta_y,mt%eta_z,mt%zeta_x,mt%zeta_y,mt%zeta_z,mt%jac,mt%xjac, b, G)
        if(inp%is_mlswe) then ! quads data points added by Yao Gahounzo
            call metrics_quad(mt%ksiq_x,mt%ksiq_y,mt%ksiq_z,mt%etaq_x,mt%etaq_y,mt%etaq_z,mt%zetaq_x, &
            mt%zetaq_y,mt%zetaq_z,mt%jacq,mt%xjacq, b, G)
        end if  ! is_mlswe

    end subroutine mod_metrics_create_metrics

    !-----------------------------------------------------------------------
    !>@brief Creates the (diag.) mass matrix and inverse mass matrix using a two step process:
    !> 1. Perform a local DSS (on processor)
    !> 2. Perform a global DSS (communicating across adjacent processors)
    !-----------------------------------------------------------------------
    subroutine mod_metrics_create_mass(G, b, par, mt)

        implicit none

        type(grid),        intent(in)    :: G
        type(basis),       intent(in)    :: b
        type(parallel_CS), intent(in)    :: par
        type(metrics),     intent(inout) :: mt

        integer inbh, ib, ip, ip_g, ic, i, j, k, e
        integer :: AllocateStatus
        real, dimension(:), allocatable :: mass_recv

        if(allocated(mt%mass)) deallocate(mt%mass,mt%massinv)
        allocate( mt%mass(G%npoin), mt%massinv(G%npoin), &
          mass_recv(par%num_send_recv_total), stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Metrics_Create_Mass **"

        !Compute Mass Matrix (on each element)
        call create_mass(G, b, mt%mass, mt%jac)

        !Compute Inverse Mass Matrix: This is what is used in the CODE

        !Do DSS on mass matrix
        call create_global_rhs(mt%mass,mass_recv,1,0)
        mt%massinv = 1.0/mt%mass

        !Deallocate
        deallocate(mass_recv)

    end subroutine mod_metrics_create_mass

end module mod_metrics
