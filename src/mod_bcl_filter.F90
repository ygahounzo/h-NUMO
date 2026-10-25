!-----------------------------------------------------------------------
!>@brief Optional element-wise modal filter for the multilayer state.
!>
!> Off unless ifilter > 0. Every ifilter baroclinic steps the layer state
!> q_df and the barotropic state qb_df are filtered with the same
!> tensor-product nodal filter F = b%fx (x) b%fy, built in mod_basis by
!> filter_init from filter_weight_type ('erf' Boyd-Vandeven, 'exp', 'quad'),
!> filter_basis_type and the strengths filter_mux, filter_muy (blend with
!> the identity: 1 = full filter, 0 = no-op), as in NUMA.
!>
!> The filtered quantity is the departure from a reference state,
!>     q <- q_ref + F[ |J| (q - q_ref) ] / |J|,
!> with filter_ref = 'none' (default, as in NUMA): q_ref = 0, the full state
!> is filtered; 'initial': q_ref = the state at the start of the run, so fixed
!> features such as the bottom-topography imprint in the layer thicknesses are
!> not smoothed. Filtering |J| (q - q_ref) leaves each
!> element's GLL-quadrature integral unchanged (F keeps mode 0 and the higher
!> Legendre modes integrate to zero under GLL quadrature), so layer mass and
!> momentum are conserved element by element. F is linear and applied to
!> every layer and to the barotropic state alike, so sum_k dp_k = pb and the
!> layer/barotropic momentum consistency are preserved.
!>
!> filter_sensor = .true. applies the filter only where it is needed: per
!> element, s_e = max_k log10( E_top / E_all ) of the layer thicknesses, with
!> E_top the energy of the orthonormal Legendre modes with max(i,j) = N of
!> (q - q_ref) and E_all that of all modes of q, mean included (Persson &
!> Peraire 2006). The element
!> is blended toward its filtered value with weight eps_e = 0 below
!> filter_s0 - filter_kappa, 1 above filter_s0 + filter_kappa and a sine ramp
!> in between. eps_e is the same for every layer and the barotropic state, so
!> the conservation and consistency properties above are kept.
!-----------------------------------------------------------------------
module mod_bcl_filter

    use mod_grid,    only: grid
    use mod_basis,   only: basis
    use mod_input,   only: input
    use mod_metrics, only: metrics

    implicit none

    private
    public :: bcl_filter_init, bcl_filter_apply, bcl_filter_report

    real, allocatable :: jac_n(:)         ! |J| at each DG node (npoin)
    real, allocatable :: q_ref(:,:,:)     ! layer reference state (nvar_bcl,npoin,nlayers)
    real, allocatable :: qb_ref(:,:)      ! barotropic reference state (nvar_btp,npoin)
    real, allocatable :: q_wk(:,:,:)      ! filtered layer state (work)
    real, allocatable :: qb_wk(:,:)       ! filtered barotropic state (work)
    real, allocatable :: vinv(:,:)        ! GLL nodal -> orthonormal Legendre coefficients (1D)
    real, allocatable :: eps_e(:)         ! per-element filter weight from the sensor (nelem)
    integer :: n_active = 0               ! elements with eps_e > 0 in the last application (this rank)
    real    :: s_max = -huge(1.0)         ! largest sensor value in the last application (this rank)

contains

    subroutine bcl_filter_init(G, b, inp, mt, q_df, qb_df)

        implicit none

        type(grid),    intent(in) :: G
        type(basis),   intent(in) :: b
        type(input),   intent(in) :: inp
        type(metrics), intent(in) :: mt
        real, intent(in) :: q_df(inp%nvar_bcl, G%npoin, inp%nlayers)
        real, intent(in) :: qb_df(inp%nvar_btp, G%npoin)

        integer :: e, i, j, n
        real    :: p0, p1, p2, gam

        if (inp%ifilter <= 0) return

        if (b%nglx /= b%ngly) error stop "bcl_filter requires nopx == nopy"

        allocate(jac_n(G%npoin), &
                 q_ref(inp%nvar_bcl, G%npoin, inp%nlayers), qb_ref(inp%nvar_btp, G%npoin), &
                 q_wk(inp%nvar_bcl, G%npoin, inp%nlayers), qb_wk(inp%nvar_btp, G%npoin))

        do e = 1, G%nelem
            do j = 1, b%ngly
                do i = 1, b%nglx
                    jac_n(G%intma(i,j,1,e)) = abs(mt%xjac(i,j,1,e))
                end do
            end do
        end do

        if (inp%filter_ref == 'initial') then
            ! State at the start of the run (host copy is current here).
            q_ref  = q_df
            qb_ref = qb_df
        else
            q_ref  = 0.0
            qb_ref = 0.0
        end if

        ! Sensor transform: c_i = sum_n w_n phi_i(x_n) q_n / gamma_i with phi_i the
        ! orthonormal Legendre polynomials and gamma_i = sum_n w_n phi_i(x_n)^2.
        ! GLL quadrature is exact for phi_i phi_j except i = j = N, so this is the
        ! exact modal expansion of the nodal interpolant.
        allocate(vinv(b%nglx, b%nglx), eps_e(G%nelem))
        do i = 1, b%nglx
            gam = 0.0
            do n = 1, b%nglx
                p0 = 1.0; p1 = b%xglx(n)
                if (i == 1) then
                    p2 = p0
                else if (i == 2) then
                    p2 = p1
                else
                    do j = 2, i-1
                        p2 = ((2*j-1)*b%xglx(n)*p1 - (j-1)*p0)/j
                        p0 = p1; p1 = p2
                    end do
                end if
                vinv(i,n) = b%wglx(n) * p2 * sqrt((2*i-1)/2.0)
                gam = gam + b%wglx(n) * (p2 * sqrt((2*i-1)/2.0))**2
            end do
            vinv(i,:) = vinv(i,:) / gam
        end do
        eps_e = 1.0

        !$acc enter data copyin(b%fx, b%fy, jac_n, q_ref, qb_ref, vinv, eps_e) create(q_wk, qb_wk)

    end subroutine bcl_filter_init

    subroutine bcl_filter_apply(G, b, inp, pbprime_df, q_df, qb_df)

        implicit none

        type(grid),  intent(in)    :: G
        type(basis), intent(in)    :: b
        type(input), intent(in)    :: inp
        real, intent(in)    :: pbprime_df(G%npoin)
        real, intent(inout) :: q_df(inp%nvar_bcl, G%npoin, inp%nlayers)
        real, intent(inout) :: qb_df(inp%nvar_btp, G%npoin)

        integer :: e, i, j, l, m, k, iv, I1, Ilm, nelem_l, ngl_l, nlayers_l, nvb_l, nvt_l, npoin_l
        integer :: n_act
        real :: s, fij, cij, dij, etop, eall, se, smax, s0, kap, ep
        real, parameter :: pi = 3.14159265358979323846

        nelem_l   = G%nelem
        ngl_l     = b%nglx
        nlayers_l = inp%nlayers
        nvb_l     = inp%nvar_bcl
        nvt_l     = inp%nvar_btp
        npoin_l   = G%npoin

        ! Sensor: per-element weight eps_e from the top-mode energy of the
        ! filtered quantity (q - q_ref) relative to the total energy of the full
        ! thickness q. With filter_ref = 'initial' the fixed topography imprint
        ! is thus excluded from the top modes, while the normalisation (set by
        ! the mean thickness) and hence filter_s0 are the same for both options.
        if (inp%filter_sensor) then
            s0  = inp%filter_s0
            kap = inp%filter_kappa
            n_act = 0
            smax  = -huge(1.0)
            !$acc parallel loop gang vector present(G%intma, vinv, q_df, q_ref, eps_e) &
            !$acc   private(se, etop, eall, cij, dij, Ilm, i, j, l, m, k, ep)     &
            !$acc   firstprivate(nlayers_l, ngl_l, s0, kap)                      &
            !$acc   reduction(+:n_act) reduction(max:smax)
            do e = 1, nelem_l
                se = -huge(1.0)
                !$acc loop seq
                do k = 1, nlayers_l
                    etop = 0.0
                    eall = 0.0
                    !$acc loop seq
                    do j = 1, ngl_l
                        !$acc loop seq
                        do i = 1, ngl_l
                            cij = 0.0
                            dij = 0.0
                            !$acc loop seq
                            do m = 1, ngl_l
                                !$acc loop seq
                                do l = 1, ngl_l
                                    Ilm = G%intma(l,m,1,e)
                                    cij = cij + vinv(i,l)*vinv(j,m)*q_df(1,Ilm,k)
                                    dij = dij + vinv(i,l)*vinv(j,m)*(q_df(1,Ilm,k) - q_ref(1,Ilm,k))
                                end do
                            end do
                            eall = eall + cij*cij
                            if (max(i,j) == ngl_l) etop = etop + dij*dij
                        end do
                    end do
                    if (eall > 0.0 .and. etop > 0.0) se = max(se, log10(etop/eall))
                end do
                if (se <= s0 - kap) then
                    ep = 0.0
                else if (se >= s0 + kap) then
                    ep = 1.0
                else
                    ep = 0.5*(1.0 + sin(0.5*pi*(se - s0)/kap))
                end if
                eps_e(e) = ep
                if (ep > 0.0) n_act = n_act + 1
                smax = max(smax, se)
            end do
            !$acc end parallel loop
            n_active = n_act
            s_max    = smax
        end if

        ! Pass 1: filtered values into the work arrays (input left intact),
        ! blended with the unfiltered state by the element weight eps_e
        ! (eps_e = 1 everywhere without the sensor).
        !$acc parallel loop gang vector collapse(3)                                    &
        !$acc   present(G%intma, b%fx, b%fy, jac_n, q_ref, qb_ref, q_wk, qb_wk, q_df, qb_df, eps_e) &
        !$acc   private(I1, Ilm, s, fij, k, iv, l, m, ep)                                &
        !$acc   firstprivate(nlayers_l, nvb_l, nvt_l, ngl_l)
        do e = 1, nelem_l
            do j = 1, ngl_l
                do i = 1, ngl_l
                    I1 = G%intma(i,j,1,e)
                    ep = eps_e(e)
                    if (ep <= 0.0) then
                        !$acc loop seq
                        do k = 1, nlayers_l
                            !$acc loop seq
                            do iv = 1, nvb_l
                                q_wk(iv,I1,k) = q_df(iv,I1,k)
                            end do
                        end do
                        !$acc loop seq
                        do iv = 1, nvt_l
                            qb_wk(iv,I1) = qb_df(iv,I1)
                        end do
                        cycle
                    end if
                    !$acc loop seq
                    do k = 1, nlayers_l
                        !$acc loop seq
                        do iv = 1, nvb_l
                            s = 0.0
                            !$acc loop seq
                            do m = 1, ngl_l
                                !$acc loop seq
                                do l = 1, ngl_l
                                    Ilm = G%intma(l,m,1,e)
                                    fij = b%fx(i,l)*b%fy(j,m)
                                    s = s + fij*jac_n(Ilm)*(q_df(iv,Ilm,k) - q_ref(iv,Ilm,k))
                                end do
                            end do
                            q_wk(iv,I1,k) = q_df(iv,I1,k) + ep*(q_ref(iv,I1,k) + s/jac_n(I1) - q_df(iv,I1,k))
                        end do
                    end do
                    ! Barotropic: pb (1) and momenta (3:nvt); pb' (2) is rebuilt below.
                    !$acc loop seq
                    do iv = 1, nvt_l
                        if (iv == 2) cycle
                        s = 0.0
                        !$acc loop seq
                        do m = 1, ngl_l
                            !$acc loop seq
                            do l = 1, ngl_l
                                Ilm = G%intma(l,m,1,e)
                                fij = b%fx(i,l)*b%fy(j,m)
                                s = s + fij*jac_n(Ilm)*(qb_df(iv,Ilm) - qb_ref(iv,Ilm))
                            end do
                        end do
                        qb_wk(iv,I1) = qb_df(iv,I1) + ep*(qb_ref(iv,I1) + s/jac_n(I1) - qb_df(iv,I1))
                    end do
                end do
            end do
        end do
        !$acc end parallel loop

        ! Pass 2: copy back; pb' = pb - pbprime keeps the barotropic state consistent.
        !$acc parallel loop gang vector present(q_wk, qb_wk, q_df, qb_df, pbprime_df) &
        !$acc   private(k, iv) firstprivate(nlayers_l, nvb_l, nvt_l)
        do I1 = 1, npoin_l
            !$acc loop seq
            do k = 1, nlayers_l
                !$acc loop seq
                do iv = 1, nvb_l
                    q_df(iv,I1,k) = q_wk(iv,I1,k)
                end do
            end do
            qb_df(1,I1) = qb_wk(1,I1)
            qb_df(2,I1) = qb_wk(1,I1) - pbprime_df(I1)
            !$acc loop seq
            do iv = 3, nvt_l
                qb_df(iv,I1) = qb_wk(iv,I1)
            end do
        end do
        !$acc end parallel loop

    end subroutine bcl_filter_apply

    ! Rank 0 prints the global fraction of elements the sensor filtered in the
    ! last application, and the largest sensor value.
    subroutine bcl_filter_report(G, inp, itime)

        use mpi
        use mod_mpi_utilities, only: irank, irank0, MPI_PRECISION

        implicit none

        type(grid),  intent(in) :: G
        type(input), intent(in) :: inp
        integer,     intent(in) :: itime

        integer :: loc(2), glob(2), ierr
        real    :: smax_g

        if (inp%ifilter <= 0 .or. .not. inp%filter_sensor) return
        loc = [n_active, G%nelem]
        call mpi_allreduce(loc, glob, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(s_max, smax_g, 1, MPI_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
        if (irank == irank0) write(*,'(a,i8,a,i7,a,i7,a,f8.4,a,f7.3)') ' filter sensor: step ', itime, &
            '  filtered elements ', glob(1), ' /', glob(2), '  (', real(glob(1))/real(glob(2)), &
            ')  max s_e ', smax_g

    end subroutine bcl_filter_report

end module mod_bcl_filter
