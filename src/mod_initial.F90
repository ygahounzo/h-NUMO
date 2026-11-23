!----------------------------------------------------------------------!
!>@brief This module contains Initial Conditions 
!>@author J.F. Kelly and F.X. Giraldo
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943-5216
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: April 03, 2023
!----------------------------------------------------------------------!
module mod_initial

    use mod_constants
    use mod_grid
    use mod_basis
    use mod_input
    
    use mod_initial_mlswe, only: bot_topo_derivatives, &
        wind_stress_coriolis, ssprk_coefficients

    use mod_tensor,  only: compute_gradient_quad, compute_gradient_df
    use mod_face,    only: face_CS
    use mod_metrics, only: metrics

    type initial
        ! module variables and parameters
        real, dimension(:,:), allocatable :: q_init, q_exact, q_ref, kvector, q_sph, coord_sph
        real, dimension(:,:), allocatable :: pi_values, shear_stress, hA
        real, dimension(:,:,:), allocatable :: hB_grad, phiA_grad, q_ref_layers
        real, dimension(:), allocatable :: rho_layers, bathymetry
        real, dimension(:), allocatable:: height, coriolis_constant
        real, dimension(:,:,:), allocatable :: q_df
        real, dimension(:), allocatable :: pbprime_df
        real, dimension(:,:,:), allocatable :: pbprime_df_face
        real, dimension(:,:,:,:), allocatable :: qb_face
        real, dimension(:,:), allocatable :: qb_df, tau_wind, tau_wind_df
        real, dimension(:), allocatable :: coriolis_df,coriolis_quad
        real, dimension(:,:), allocatable :: coriolis_3d_quad  ! (3,npoin_q): cf = f*r̂ at quad pts
        real, dimension(:), allocatable :: zbot, zbot_df, fdt_bcl, fdt2_bcl, a_bcl, b_bcl
        real, dimension(:,:,:), allocatable :: zbot_face
        real, dimension(:,:), allocatable :: grad_zbot_quad, grad_zbot_df, z_interface

        real, dimension(:,:), allocatable :: ssprk_a, ssprk_a_bcl
        real, dimension(:),   allocatable :: ssprk_beta, ssprk_beta_bcl, alpha_mlswe
        real, dimension(:),   allocatable :: ssprk_w  ! effective stage weights b_i of the BTP SSPRK step (sum 1)

        real, dimension(:,:), allocatable :: z_init_flag, z_interface_initial, z_init_flag_elem
        real, dimension(:),   allocatable :: zbot_df_init

        integer :: nvar, nvart, nvar_diag, ntracers
        integer :: nrhs_mxm, N_btp

    end type initial

    public :: mod_initial_create, initial

    private

    contains

    !-----------------------------------------------------------------------
    subroutine mod_initial_create(inp, b, G, mf, init, mt)

        implicit none

        ! global
        type(input),   intent(inout)        :: inp
        type(basis),   intent(in)           :: b
        type(grid),    intent(in)           :: G
        type(face_CS), intent(in)           :: mf
        type(metrics), intent(in), optional :: mt
        type(initial), intent(out)          :: init

        integer i, iperturbation, ip
        real time
        real x, y, z, xf, yf, zf, radius
        real lat, lon, press, temp, phis, ps
        real pb, tb, pi_f, pi_b, zradius
        integer :: nlayers, npoin, npoin_q, npts, nface, kstages, nq, k

        !Define the number of prognostic variables
        init%nvar = 5 !rho,u,v,w,theta
        init%nvart = init%nvar
        !Define the number of diagnostic variables
        init%nvar_diag = 0

        !Store Number of RHS for MXM calls in CREATE_RHS_VOLUME
        init%nrhs_mxm=init%nvar + 1 !1 is for Pressure

        !Store Number of Tracers
        init%ntracers=init%nvar-5

        npoin = G%npoin
        nlayers = inp%nlayers
        npoin_q = G%npoin_q
        kstages = inp%kstages
        nq = b%nq
        npts = b%npts
        nface = G%nface

        if(allocated(init%q_init)) deallocate(init%q_init,init%q_exact,init%q_ref,init%kvector,init%pi_values,init%height, &
                init%coriolis_constant, init%shear_stress)
        allocate( init%q_init(init%nvar,npoin), init%q_exact(init%nvar,npoin), init%q_ref(init%nvar,npoin), init%kvector(3,npoin), &
            init%pi_values(3,npoin), init%height(npoin), init%coriolis_constant(npoin), init%shear_stress(3,npoin))

        !Set-up kvector (must come after init%kvector is allocated)
        call create_kvector(init%kvector,G%coord,npoin,inp%geometry_type)

        !hack to allocate layers stuff anyways
        allocate(init%rho_layers(1),init%bathymetry(npoin))

        if(inp%is_mlswe) then
            if(allocated(init%q_df)) deallocate(init%q_df, init%pbprime_df, & 
            init%qb_df, init%alpha_mlswe, init%tau_wind, init%coriolis_quad, &
            init%coriolis_df, init%coriolis_3d_quad, init%zbot, init%zbot_df, init%zbot_face, &
            init%grad_zbot_quad, init%tau_wind_df,&
            init%ssprk_a,init%ssprk_beta, init%ssprk_w, init%ssprk_a_bcl, init%ssprk_beta_bcl, init%grad_zbot_df, &
            init%pbprime_df_face, init%z_interface, &
            init%z_init_flag, init%z_interface_initial,                     &
            init%z_init_flag_elem, init%zbot_df_init)
            allocate(init%q_df(inp%nvar_bcl,npoin,nlayers), init%pbprime_df(npoin), &
            init%qb_df(inp%nvar_btp,npoin), &
            init%alpha_mlswe(nlayers), init%tau_wind(2,npoin_q), init%coriolis_quad(npoin_q), init%coriolis_df(npoin), &
            init%coriolis_3d_quad(3,npoin_q), &
            init%zbot(npoin_q), init%zbot_df(npoin), init%zbot_face(2,nq,nface), init%grad_zbot_quad(inp%ngrd_var,npoin_q), &
            init%grad_zbot_df(inp%ngrd_var,npoin), &
            init%fdt_bcl(npoin), init%fdt2_bcl(npoin), init%a_bcl(npoin), &
            init%b_bcl(npoin), &
            init%tau_wind_df(2,npoin), init%ssprk_a(kstages,3), init%ssprk_beta(kstages), init%ssprk_w(kstages), &
            init%ssprk_a_bcl(inp%kstages_bcl,3), init%ssprk_beta_bcl(inp%kstages_bcl), &
            init%pbprime_df_face(2,b%ngl,nface),init%z_interface(npoin,nlayers+1), &
            init%z_init_flag(npoin,nlayers), init%z_interface_initial(npoin,nlayers+1),  &
            init%z_init_flag_elem(G%nelem,nlayers), init%zbot_df_init(npoin))

            init%q_df = 0.0
            init%pbprime_df = 0.0

        end if

        !Initialize q_* arrays to zero:
        init%q_init  = 0.0
        init%q_exact = 0.0
        init%q_ref   = 0.0

        time = 0.0
    
        if(inp%is_mlswe) then

            if (is_sphere_geometry(inp%geometry_type)) then
                call initial_conditions_sphere(init%q_df, init%pbprime_df, init%qb_df, &
                    init%alpha_mlswe, init%pbprime_df_face, init%zbot_df, init%tau_wind_df, init%kvector, &
                    nlayers, b, G, inp, mf)
                ! Build reference interface heights.
                ! z_interface_initial(:,k) = bottom of layer k-1 = top of layer k.
                init%z_interface_initial(:, nlayers+1) = init%zbot_df(:)
                do k = nlayers, 1, -1
                    init%z_interface_initial(:, k) = init%z_interface_initial(:, k+1) + &
                        (init%alpha_mlswe(k)/gravity) * init%q_df(1,:,k)
                end do
            else
                call initial_conditions(inp, G, b, mf, init)
            end if

            init%zbot_df_init = init%zbot_df
                
            call bot_topo_derivatives(b, G, inp, mf, init%zbot, init%zbot_face, init%zbot_df)

            if(present(mt)) then
                call compute_gradient_quad(G, inp, b, mt, init%grad_zbot_quad, init%zbot_df)
                call compute_gradient_df(G, inp, b, mt, init%grad_zbot_df, init%zbot_df)
            end if

            ! Collocated BTP: keep dt_btp at most 98% of its estimated stability limit
            if (inp%btp_collocation) call btp_colloc_dt_limit(inp, b, G, init)

            init%N_btp = ceiling(inp%dt/inp%dt_btp)
            inp%dt_btp = inp%dt/real(init%N_btp)

            ! No bcl-btp time splitting with a single layer: use
            ! dt_btp (sized for the barotropic CFL) as the master timestep
            ! directly, one barotropic step per call.
            if (inp%nlayers == 1) then
                init%N_btp = 1
                inp%dt     = inp%dt_btp
            end if

            call  wind_stress_coriolis(b, G, inp, init%tau_wind, init%coriolis_df, init%coriolis_quad, init%fdt_bcl, init%fdt2_bcl, &
                init%a_bcl, init%b_bcl, init%tau_wind_df, init%coriolis_3d_quad, init%kvector)
            call ssprk_coefficients(inp, init%ssprk_a, init%ssprk_beta, init%ssprk_w, init%ssprk_a_bcl, init%ssprk_beta_bcl)

        endif

        !Set-up Times
        inp%time_initial = inp%time_initial*inp%time_scale
        inp%time_final = inp%time_final*inp%time_scale
        inp%time_restart = inp%time_restart*inp%time_scale
        inp%write_restart_time = inp%write_restart_time*inp%time_scale

    end subroutine mod_initial_create

    subroutine btp_colloc_dt_limit(inp, b, G, init)

        ! Stability limit of the collocated barotropic solver (SSPRK(5,3), GLL nodes):
        ! nodal Courant number 1, dt_lim = min over elements of
        !     (smallest GLL node spacing) / max sqrt(alpha_L * p_b),
        ! with sqrt(alpha_L p_b) the barotropic gravity-wave speed (as in the Rusanov flux).
        ! Calibrated on the Cartesian bump (N = 4): estimate 1.74 s, measured limit 1.76 s.
        ! dt_btp is reduced to 98% of dt_lim if larger; a smaller dt_btp is kept.

        use mod_mpi_utilities, only: irank, irank0, MPI_PRECISION
        use mpi

        implicit none

        type(input),   intent(inout) :: inp
        type(basis),   intent(in)    :: b
        type(grid),    intent(in)    :: G
        type(initial), intent(in)    :: init

        integer :: ie, i, m, n, I00, IN0, I0N, ierr
        real    :: gap_x, gap_y, dx_len, dy_len, ds, c, dt_e, dt_loc, dt_lim

        ! Smallest GLL node gap as a fraction of the element edge
        gap_x = 1.0
        do i = 1, b%nglx - 1
            gap_x = min(gap_x, 0.5 * (b%xglx(i+1) - b%xglx(i)))
        end do
        gap_y = 1.0
        do i = 1, b%ngly - 1
            gap_y = min(gap_y, 0.5 * (b%xgly(i+1) - b%xgly(i)))
        end do

        dt_loc = huge(1.0)
        do ie = 1, G%nelem

            ! Element edge lengths from the corner nodes
            I00 = G%intma(1, 1, 1, ie)
            IN0 = G%intma(b%nglx, 1, 1, ie)
            I0N = G%intma(1, b%ngly, 1, ie)
            dx_len = sqrt(sum((G%coord(:,IN0) - G%coord(:,I00))**2))
            dy_len = sqrt(sum((G%coord(:,I0N) - G%coord(:,I00))**2))
            ds     = min(gap_x * dx_len, gap_y * dy_len)

            ! Largest barotropic gravity-wave speed in the element
            c = 0.0
            do m = 1, b%ngly
                do n = 1, b%nglx
                    c = max(c, sqrt(init%alpha_mlswe(inp%nlayers) * init%pbprime_df(G%intma(n, m, 1, ie))))
                end do
            end do

            if (c > 0.0) then
                dt_e   = ds / c
                dt_loc = min(dt_loc, dt_e)
            end if

        end do

        call mpi_allreduce(dt_loc, dt_lim, 1, MPI_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

        if (irank == irank0) then
            print '(a,es11.4,a,es11.4,a)', ' btp_collocation: estimated dt_btp limit ', dt_lim, &
                                           ' s (98%: ', 0.98*dt_lim, ' s)'
        end if

        if (inp%dt_btp > 0.98 * dt_lim) then
            if (irank == irank0) print '(a,es11.4,a,es11.4,a)', ' btp_collocation: dt_btp reduced from ', &
                                                                inp%dt_btp, ' to ', 0.98*dt_lim, ' s'
            inp%dt_btp = 0.98 * dt_lim
        end if

    end subroutine btp_colloc_dt_limit

    !----------------------------------------------------------------------!
    !>@brief This subroutine constructs the KVECTOR on the sphere
    !> (called the Rvector in Notes and Papers)
    !----------------------------------------------------------------------!
    subroutine create_kvector(kvector, coord, npoin, geometry_type)

        implicit none

        integer,          intent(in)  :: npoin
        real,             intent(out) :: kvector(3,npoin)
        real,             intent(in)  :: coord(3,npoin)
        character(len=*), intent(in)  :: geometry_type

        integer :: ip
        real    :: radius

        if (is_sphere_geometry(geometry_type)) then
            do ip = 1, npoin
                radius = sqrt(dot_product(coord(:,ip), coord(:,ip)))
                if (radius > 0.0) then
                    kvector(:,ip) = coord(:,ip) / radius
                else
                    kvector(1,ip) = 0.0; kvector(2,ip) = 0.0; kvector(3,ip) = 1.0
                end if
            end do
        else
            ! Cartesian: vertical direction is always z
            do ip = 1, npoin
                kvector(1,ip) = 0.0; kvector(2,ip) = 0.0; kvector(3,ip) = 1.0
            end do
        end if

    end subroutine create_kvector

end module mod_initial
