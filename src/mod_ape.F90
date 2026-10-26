!-----------------------------------------------------------------------
!>@brief Optional available-potential-energy (APE) term of Chen (2025, QJRMS).
!>
!> Off unless c_APE > 0 (sphere, nlayers > 1). The Hamiltonian gains
!>     H_APE = sum_{k=2..nlayers} c_APE (rho_k - rho_{k-1})/2 g |grad z_k|^2
!> for the interior interfaces z_k = zbot + sum_{m>=k} h_m (the free surface,
!> k = 1, is excluded: Chen's nu_1 = 0). Its variational force on the
!> momentum dp_k u_k of layer k is
!>     F_k = g c_APE dp_k alpha_k sum_{j=2..k} (rho_j - rho_{j-1}) grad(lap z_j).
!>
!> Discretisation: nodal DG on the GLL points with the diagonal mass wjac_df,
!> BR1 (Bassi-Rebay) local DG with central fluxes:
!>     q   = grad z + M^-1 sum_faces psi (1/2)[z] n       (strong, lifted)
!>     L   = M^-1 ( -int grad psi . q + sum_faces psi {q}.n )   (weak)
!>     G_L = grad L + M^-1 sum_faces psi (1/2)[L] n
!> The two BR1 operators are negative adjoints in the wjac_df inner
!> product, so the energy that matches this force is
!>     E_APE = sum_j c_APE drho_j g/2 sum_I wjac_df(I) |q_j(I)|^2,
!> which ape_energy returns for energy_runtime.csv. Only the surface
!> gradient is used (dpsidx_df, dpsidy_df, dpsidz_df: grad xi and grad eta
!> are tangent to the sphere), so lap(const) = 0 exactly and F_k is tangent.
!>
!> ape_stage is called once per baroclinic stage, after btp_bcl_coeffs_qdf
!> and before the barotropic subcycle, from the stage state. It stores F_k at
!> the quadrature points (ape_force) and its column sum (ape_force_btp). The
!> baroclinic volume kernel adds ape_force(:,Iq,k) to layer k and the
!> barotropic volume kernel adds ape_force_btp(:,Iq), held fixed over the
!> subcycle like the other btp_bcl_coeffs_qdf terms, so sum_k (layer forcing)
!> = barotropic forcing exactly and the layer/barotropic momentum
!> consistency is unaffected.
!>
!> All kernels run on the device; MPI faces use their own GPU-direct
!> exchange (same face list and message layout as the bcl Laplacian).
!-----------------------------------------------------------------------
module mod_ape

    use mod_grid,     only: grid
    use mod_basis,    only: basis
    use mod_input,    only: input
    use mod_face,     only: face_CS
    use mod_parallel, only: parallel_CS
    use mod_ref,      only: mref
    use mod_tensor,   only: tensor_CS

    implicit none

    private
    public :: ape_init, ape_stage, ape_energy, ape_on, ape_force, ape_force_btp

    logical :: ape_on = .false.
    real, allocatable :: ape_force(:,:,:)     ! layer force at quad points (3,npoin_q,nlayers)
    real, allocatable :: ape_force_btp(:,:)   ! its column sum (3,npoin_q)

    integer :: nint = 0                       ! interior interfaces (nlayers-1)
    real, allocatable :: zf(:,:)              ! z_{j+1} at DG nodes (npoin,nint)
    real, allocatable :: qf(:,:,:)            ! BR1 gradient of z (3,npoin,nint)
    real, allocatable :: lf(:,:)              ! lap z (npoin,nint)
    real, allocatable :: gf(:,:,:)            ! BR1 gradient of lap z (3,npoin,nint)
    real, allocatable :: drho(:)              ! rho_{j+1} - rho_j (nint)
    real, allocatable :: sbuf(:), rbuf(:)     ! MPI face buffers
    integer, allocatable :: ireq(:)
    integer :: nreq = 0
    integer, parameter :: ape_tag = 77

contains

    subroutine ape_init(G, b, inp, ref, alpha_mlswe)

        implicit none

        type(grid),  intent(in) :: G
        type(basis), intent(in) :: b
        type(input), intent(in) :: inp
        type(mref),  intent(in) :: ref
        real,        intent(in) :: alpha_mlswe(inp%nlayers)

        integer :: j

        if (allocated(ape_force)) return

        ape_on = (inp%c_APE > 0.0 .and. inp%nlayers > 1)
        if (ape_on .and. inp%nvar_bcl /= 4) error stop "c_APE > 0 is implemented on the sphere only"

        if (.not. ape_on) then
            ! Placeholders so the rhs kernels can name them in present clauses.
            allocate(ape_force(3,1,1), ape_force_btp(3,1))
            ape_force = 0.0;  ape_force_btp = 0.0
            !$acc enter data copyin(ape_force, ape_force_btp)
            return
        end if

        nint = inp%nlayers - 1
        allocate(ape_force(3,G%npoin_q,inp%nlayers), ape_force_btp(3,G%npoin_q), &
                 zf(G%npoin,nint), qf(3,G%npoin,nint), lf(G%npoin,nint), gf(3,G%npoin,nint), &
                 drho(nint), sbuf(3*nint*b%ngl*max(ref%nboun_valid,1)), &
                 rbuf(3*nint*b%ngl*max(ref%nboun_valid,1)))
        do j = 1, nint
            drho(j) = 1.0/alpha_mlswe(j+1) - 1.0/alpha_mlswe(j)
        end do
        ape_force = 0.0;  ape_force_btp = 0.0

        !$acc enter data copyin(ape_force, ape_force_btp, drho) create(zf, qf, lf, gf, sbuf, rbuf)

    end subroutine ape_init

    ! Force for the current stage state q_df (on the device).
    subroutine ape_stage(G, b, inp, mf, par, ref, tsp, alpha_mlswe, zbot_df, q_df)

        use mod_constants, only: gravity

        implicit none

        type(grid),        intent(in) :: G
        type(basis),       intent(in) :: b
        type(input),       intent(in) :: inp
        type(face_CS),     intent(in) :: mf
        type(parallel_CS), intent(in) :: par
        type(mref),        intent(in) :: ref
        type(tensor_CS),   intent(in) :: tsp
        real, intent(in) :: alpha_mlswe(inp%nlayers), zbot_df(G%npoin)
        real, intent(in) :: q_df(inp%nvar_bcl, G%npoin, inp%nlayers)

        integer :: ie, iqx, iqy, Iq, n, m, I, k, nelem_l, nqx_l, nqy_l, nglx_l, ngly_l, nlayers_l
        real    :: hi, dp, ax, ay, az, gx, gy, gz, c_g, fx, fy, fz, sx, sy, sz

        if (.not. ape_on) return

        ! The barotropic volume kernels read ape_force_btp on async queue 1.
        !$acc wait

        call ape_interfaces(G, inp, alpha_mlswe, zbot_df, q_df)
        call br1_gradient(G, b, mf, par, ref, tsp, zf, qf)
        call br1_divergence(G, b, mf, par, ref, tsp, qf, lf)
        call br1_gradient(G, b, mf, par, ref, tsp, lf, gf)

        ! Force at the quad points, top to bottom: A accumulates
        ! sum_{j<k} drho_j G_L,j, so no per-thread layer arrays are needed.
        nelem_l = G%nelem;  nqx_l = b%nqx;  nqy_l = b%nqy
        nglx_l  = b%nglx;   ngly_l = b%ngly;  nlayers_l = inp%nlayers
        c_g     = gravity*inp%c_APE

        !$acc parallel loop gang vector collapse(3)                                    &
        !$acc   present(G%intma, G%intma_dg_quad, b%psiqx, b%psiqy, q_df, gf, drho,    &
        !$acc           alpha_mlswe, ape_force, ape_force_btp)                         &
        !$acc   firstprivate(nqx_l, nqy_l, nglx_l, ngly_l, nlayers_l, c_g)             &
        !$acc   private(Iq, n, m, I, k, hi, dp, ax, ay, az, gx, gy, gz, fx, fy, fz, sx, sy, sz)
        do ie = 1, nelem_l
        do iqy = 1, nqy_l
        do iqx = 1, nqx_l
            Iq = G%intma_dg_quad(iqx, iqy, 1, ie)
            ax = 0.0;  ay = 0.0;  az = 0.0
            sx = 0.0;  sy = 0.0;  sz = 0.0
            !$acc loop seq
            do k = 1, nlayers_l
                dp = 0.0;  gx = 0.0;  gy = 0.0;  gz = 0.0
                !$acc loop seq
                do m = 1, ngly_l
                    !$acc loop seq
                    do n = 1, nglx_l
                        I  = G%intma(n, m, 1, ie)
                        hi = b%psiqx(n, iqx)*b%psiqy(m, iqy)
                        dp = dp + hi*q_df(1,I,k)
                        if (k < nlayers_l) then
                            gx = gx + hi*gf(1,I,k)
                            gy = gy + hi*gf(2,I,k)
                            gz = gz + hi*gf(3,I,k)
                        end if
                    end do
                end do
                fx = c_g*alpha_mlswe(k)*dp*ax
                fy = c_g*alpha_mlswe(k)*dp*ay
                fz = c_g*alpha_mlswe(k)*dp*az
                ape_force(1,Iq,k) = fx
                ape_force(2,Iq,k) = fy
                ape_force(3,Iq,k) = fz
                sx = sx + fx;  sy = sy + fy;  sz = sz + fz
                if (k < nlayers_l) then
                    ax = ax + drho(k)*gx
                    ay = ay + drho(k)*gy
                    az = az + drho(k)*gz
                end if
            end do
            ape_force_btp(1,Iq) = sx
            ape_force_btp(2,Iq) = sy
            ape_force_btp(3,Iq) = sz
        end do
        end do
        end do
        !$acc end parallel loop

    end subroutine ape_stage

    ! Global E_APE of the state q_df (on the device); 0 when off.
    function ape_energy(G, b, inp, mf, par, ref, tsp, alpha_mlswe, zbot_df, q_df) result(e_global)

        use mpi
        use mod_mpi_utilities, only: MPI_PRECISION
        use mod_constants,     only: gravity

        implicit none

        type(grid),        intent(in) :: G
        type(basis),       intent(in) :: b
        type(input),       intent(in) :: inp
        type(face_CS),     intent(in) :: mf
        type(parallel_CS), intent(in) :: par
        type(mref),        intent(in) :: ref
        type(tensor_CS),   intent(in) :: tsp
        real, intent(in) :: alpha_mlswe(inp%nlayers), zbot_df(G%npoin)
        real, intent(in) :: q_df(inp%nvar_bcl, G%npoin, inp%nlayers)
        real :: e_global

        real    :: e_local
        integer :: I, j, npoin_l, nint_l, ierr

        e_global = 0.0
        if (.not. ape_on) return

        call ape_interfaces(G, inp, alpha_mlswe, zbot_df, q_df)
        call br1_gradient(G, b, mf, par, ref, tsp, zf, qf)

        npoin_l = G%npoin;  nint_l = nint
        e_local = 0.0
        !$acc parallel loop gang vector collapse(2) present(tsp%wjac_df, qf, drho) &
        !$acc   firstprivate(npoin_l, nint_l) reduction(+:e_local)
        do j = 1, nint_l
            do I = 1, npoin_l
                e_local = e_local + drho(j)*tsp%wjac_df(I)* &
                          (qf(1,I,j)**2 + qf(2,I,j)**2 + qf(3,I,j)**2)
            end do
        end do
        !$acc end parallel loop
        e_local = 0.5*gravity*inp%c_APE*e_local

        call mpi_allreduce(e_local, e_global, 1, MPI_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    end function ape_energy

    ! zf(:,j) = z_{j+1} = zbot + sum_{m>=j+1} (alpha_m/g) dp_m.
    subroutine ape_interfaces(G, inp, alpha_mlswe, zbot_df, q_df)

        use mod_constants, only: gravity

        implicit none

        type(grid),  intent(in) :: G
        type(input), intent(in) :: inp
        real, intent(in) :: alpha_mlswe(inp%nlayers), zbot_df(G%npoin)
        real, intent(in) :: q_df(inp%nvar_bcl, G%npoin, inp%nlayers)

        integer :: I, k, npoin_l, nlayers_l
        real    :: z

        npoin_l = G%npoin;  nlayers_l = inp%nlayers

        !$acc parallel loop gang vector present(zbot_df, q_df, alpha_mlswe, zf) &
        !$acc   firstprivate(npoin_l, nlayers_l) private(z, k)
        do I = 1, npoin_l
            z = zbot_df(I)
            !$acc loop seq
            do k = nlayers_l, 2, -1
                z = z + (alpha_mlswe(k)/gravity)*q_df(1,I,k)
                zf(I,k-1) = z
            end do
        end do
        !$acc end parallel loop

    end subroutine ape_interfaces

    ! v = BR1 gradient of s: strong element gradient plus the lifted jump
    ! (1/2)[s] n on every face.
    subroutine br1_gradient(G, b, mf, par, ref, tsp, s, v)

        implicit none

        type(grid),        intent(in)  :: G
        type(basis),       intent(in)  :: b
        type(face_CS),     intent(in)  :: mf
        type(parallel_CS), intent(in)  :: par
        type(mref),        intent(in)  :: ref
        type(tensor_CS),   intent(in)  :: tsp
        real, intent(in)  :: s(G%npoin, nint)
        real, intent(out) :: v(3, G%npoin, nint)

        integer :: I, Jn, ip, j, iface, n, IL, IR, kk, idx, npoin_l, npts_l, ngl_l, nface_l, nint_l, nbv_l
        real    :: gx, gy, gz, wq, d, nx, ny, nz

        npoin_l = G%npoin;  npts_l = b%npts;  ngl_l = b%ngl;  nface_l = G%nface
        nint_l  = nint;     nbv_l  = ref%nboun_valid

        call exchange_start(G, b, mf, par, ref, s, 1)

        !$acc parallel loop gang vector collapse(2)                                   &
        !$acc   present(tsp%index_df, tsp%dpsidx_df, tsp%dpsidy_df, tsp%dpsidz_df, s, v) &
        !$acc   firstprivate(npoin_l, npts_l, nint_l) private(gx, gy, gz, Jn)
        do j = 1, nint_l
            do I = 1, npoin_l
                gx = 0.0;  gy = 0.0;  gz = 0.0
                !$acc loop seq
                do ip = 1, npts_l
                    Jn = tsp%index_df(ip, I)
                    gx = gx + tsp%dpsidx_df(ip, I)*s(Jn, j)
                    gy = gy + tsp%dpsidy_df(ip, I)*s(Jn, j)
                    gz = gz + tsp%dpsidz_df(ip, I)*s(Jn, j)
                end do
                v(1,I,j) = gx;  v(2,I,j) = gy;  v(3,I,j) = gz
            end do
        end do
        !$acc end parallel loop

        ! Interior faces: the same jump term on both sides (n_R = -n_L and
        ! [s] changes sign with it).
        !$acc parallel loop gang vector collapse(2)                                   &
        !$acc   present(G%face, G%face_type, G%intma, mf%imapl, mf%imapr, mf%jac_face, mf%normal_vector, &
        !$acc           tsp%wjac_df, s, v)                                             &
        !$acc   firstprivate(nface_l, ngl_l, nint_l) private(IL, IR, wq, d, nx, ny, nz, j)
        do iface = 1, nface_l
            do n = 1, ngl_l
                if (G%face_type(iface) /= 2 .and. G%face(8,iface) > 0) then
                    IL = G%intma(mf%imapl(1,n,1,iface), mf%imapl(2,n,1,iface), mf%imapl(3,n,1,iface), G%face(7,iface))
                    IR = G%intma(mf%imapr(1,n,1,iface), mf%imapr(2,n,1,iface), mf%imapr(3,n,1,iface), G%face(8,iface))
                    wq = 0.5*mf%jac_face(n,1,iface)
                    nx = mf%normal_vector(1,n,1,iface)
                    ny = mf%normal_vector(2,n,1,iface)
                    nz = mf%normal_vector(3,n,1,iface)
                    !$acc loop seq
                    do j = 1, nint_l
                        d = wq*(s(IR,j) - s(IL,j))
                        !$acc atomic update
                        v(1,IL,j) = v(1,IL,j) + d*nx/tsp%wjac_df(IL)
                        !$acc atomic update
                        v(2,IL,j) = v(2,IL,j) + d*ny/tsp%wjac_df(IL)
                        !$acc atomic update
                        v(3,IL,j) = v(3,IL,j) + d*nz/tsp%wjac_df(IL)
                        !$acc atomic update
                        v(1,IR,j) = v(1,IR,j) + d*nx/tsp%wjac_df(IR)
                        !$acc atomic update
                        v(2,IR,j) = v(2,IR,j) + d*ny/tsp%wjac_df(IR)
                        !$acc atomic update
                        v(3,IR,j) = v(3,IR,j) + d*nz/tsp%wjac_df(IR)
                    end do
                end if
            end do
        end do
        !$acc end parallel loop

        call exchange_finish()
        if (nbv_l == 0) return

        ! MPI faces: the local element is side 1; neighbour values from rbuf.
        !$acc parallel loop gang vector collapse(2)                                   &
        !$acc   present(ref%face_pack_list, G%face, G%intma, mf%imapl, mf%jac_face, mf%normal_vector,  &
        !$acc           tsp%wjac_df, s, v, rbuf)                                       &
        !$acc   firstprivate(nbv_l, ngl_l, nint_l) private(iface, IL, wq, d, nx, ny, nz, j, idx)
        do kk = 1, nbv_l
            do n = 1, ngl_l
                iface = ref%face_pack_list(kk)
                IL = G%intma(mf%imapl(1,n,1,iface), mf%imapl(2,n,1,iface), mf%imapl(3,n,1,iface), G%face(7,iface))
                wq = 0.5*mf%jac_face(n,1,iface)/tsp%wjac_df(IL)
                nx = mf%normal_vector(1,n,1,iface)
                ny = mf%normal_vector(2,n,1,iface)
                nz = mf%normal_vector(3,n,1,iface)
                !$acc loop seq
                do j = 1, nint_l
                    idx = ((kk-1)*ngl_l + (n-1))*nint_l + j
                    d = wq*(rbuf(idx) - s(IL,j))
                    !$acc atomic update
                    v(1,IL,j) = v(1,IL,j) + d*nx
                    !$acc atomic update
                    v(2,IL,j) = v(2,IL,j) + d*ny
                    !$acc atomic update
                    v(3,IL,j) = v(3,IL,j) + d*nz
                end do
            end do
        end do
        !$acc end parallel loop

    end subroutine br1_gradient

    ! s = BR1 divergence of v: M^-1 ( -int grad psi . v + sum_faces psi {v}.n ).
    subroutine br1_divergence(G, b, mf, par, ref, tsp, v, s)

        implicit none

        type(grid),        intent(in)  :: G
        type(basis),       intent(in)  :: b
        type(face_CS),     intent(in)  :: mf
        type(parallel_CS), intent(in)  :: par
        type(mref),        intent(in)  :: ref
        type(tensor_CS),   intent(in)  :: tsp
        real, intent(in)  :: v(3, G%npoin, nint)
        real, intent(out) :: s(G%npoin, nint)

        integer :: I, Iq, ie, n, m, ipI, ipq, j, iface, IL, IR, kk, idx
        integer :: npoin_l, npts_l, ngl_l, nglx_l, ngly_l, nface_l, nint_l, nbv_l, nelem_l
        real    :: acc_s, wq, f

        npoin_l = G%npoin;  npts_l = b%npts;  ngl_l = b%ngl;  nface_l = G%nface
        nglx_l  = b%nglx;   ngly_l = b%ngly;  nelem_l = G%nelem
        nint_l  = nint;     nbv_l  = ref%nboun_valid

        call exchange_start(G, b, mf, par, ref, v, 3)

        ! Weak volume term, gathered per node (each DG node belongs to one
        ! element, so no atomics): node (n,m) of element ie is basis
        ! function ipI = n + (m-1)*nglx in the tensor tables.
        !$acc parallel loop gang vector collapse(4)                                   &
        !$acc   present(G%intma, tsp%index_df_elt, tsp%wjac_df, tsp%dpsidx_df,        &
        !$acc           tsp%dpsidy_df, tsp%dpsidz_df, v, s)                            &
        !$acc   firstprivate(nelem_l, nglx_l, ngly_l, npts_l, nint_l) private(I, ipI, Iq, acc_s)
        do j = 1, nint_l
        do ie = 1, nelem_l
        do m = 1, ngly_l
        do n = 1, nglx_l
            I   = G%intma(n, m, 1, ie)
            ipI = n + (m-1)*nglx_l
            acc_s = 0.0
            !$acc loop seq
            do ipq = 1, npts_l
                Iq = tsp%index_df_elt(ipq, ie)
                acc_s = acc_s + tsp%wjac_df(Iq)*(tsp%dpsidx_df(ipI,Iq)*v(1,Iq,j) + &
                        tsp%dpsidy_df(ipI,Iq)*v(2,Iq,j) + tsp%dpsidz_df(ipI,Iq)*v(3,Iq,j))
            end do
            s(I,j) = -acc_s
        end do
        end do
        end do
        end do
        !$acc end parallel loop

        ! Interior faces: central flux {v}.n_L, + on the left, - on the right.
        !$acc parallel loop gang vector collapse(2)                                   &
        !$acc   present(G%face, G%face_type, G%intma, mf%imapl, mf%imapr, mf%jac_face, mf%normal_vector, v, s) &
        !$acc   firstprivate(nface_l, ngl_l, nint_l) private(IL, IR, wq, f, j)
        do iface = 1, nface_l
            do n = 1, ngl_l
                if (G%face_type(iface) /= 2 .and. G%face(8,iface) > 0) then
                    IL = G%intma(mf%imapl(1,n,1,iface), mf%imapl(2,n,1,iface), mf%imapl(3,n,1,iface), G%face(7,iface))
                    IR = G%intma(mf%imapr(1,n,1,iface), mf%imapr(2,n,1,iface), mf%imapr(3,n,1,iface), G%face(8,iface))
                    wq = 0.5*mf%jac_face(n,1,iface)
                    !$acc loop seq
                    do j = 1, nint_l
                        f = wq*((v(1,IL,j) + v(1,IR,j))*mf%normal_vector(1,n,1,iface) + &
                                (v(2,IL,j) + v(2,IR,j))*mf%normal_vector(2,n,1,iface) + &
                                (v(3,IL,j) + v(3,IR,j))*mf%normal_vector(3,n,1,iface))
                        !$acc atomic update
                        s(IL,j) = s(IL,j) + f
                        !$acc atomic update
                        s(IR,j) = s(IR,j) - f
                    end do
                end if
            end do
        end do
        !$acc end parallel loop

        call exchange_finish()

        if (nbv_l > 0) then
            !$acc parallel loop gang vector collapse(2)                               &
            !$acc   present(ref%face_pack_list, G%face, G%intma, mf%imapl, mf%jac_face, mf%normal_vector, v, s, rbuf) &
            !$acc   firstprivate(nbv_l, ngl_l, nint_l) private(iface, IL, wq, f, j, idx)
            do kk = 1, nbv_l
                do n = 1, ngl_l
                    iface = ref%face_pack_list(kk)
                    IL = G%intma(mf%imapl(1,n,1,iface), mf%imapl(2,n,1,iface), mf%imapl(3,n,1,iface), G%face(7,iface))
                    wq = 0.5*mf%jac_face(n,1,iface)
                    !$acc loop seq
                    do j = 1, nint_l
                        idx = (((kk-1)*ngl_l + (n-1))*nint_l + (j-1))*3
                        f = wq*((v(1,IL,j) + rbuf(idx+1))*mf%normal_vector(1,n,1,iface) + &
                                (v(2,IL,j) + rbuf(idx+2))*mf%normal_vector(2,n,1,iface) + &
                                (v(3,IL,j) + rbuf(idx+3))*mf%normal_vector(3,n,1,iface))
                        !$acc atomic update
                        s(IL,j) = s(IL,j) + f
                    end do
                end do
            end do
            !$acc end parallel loop
        end if

        !$acc parallel loop gang vector collapse(2) present(tsp%wjac_df, s) firstprivate(npoin_l, nint_l)
        do j = 1, nint_l
            do I = 1, npoin_l
                s(I,j) = s(I,j)/tsp%wjac_df(I)
            end do
        end do
        !$acc end parallel loop

    end subroutine br1_divergence

    ! Pack a(1:nc, node, 1:nint) at the MPI-face nodes of the local element
    ! into sbuf and post the non-blocking GPU-direct sends/receives; the
    ! layout per face kk, node n, interface j, component ic is
    ! (((kk-1)*ngl + n-1)*nint + j-1)*nc + ic, as the face kernels read it.
    subroutine exchange_start(G, b, mf, par, ref, a, nc)

        use mpi
        use mod_mpi_utilities, only: MPI_PRECISION

        implicit none

        type(grid),        intent(in) :: G
        type(basis),       intent(in) :: b
        type(face_CS),     intent(in) :: mf
        type(parallel_CS), intent(in) :: par
        type(mref),        intent(in) :: ref
        integer,           intent(in) :: nc
        real,              intent(in) :: a(nc, G%npoin, nint)

        integer :: kk, n, j, ic, iface, IL, base, nbv_l, ngl_l, nint_l
        integer :: inbh, ib, i, jj, nqp, istart, iend, ierr

        nreq = 0
        nbv_l = ref%nboun_valid
        if (nbv_l == 0 .or. par%num_nbh == 0) return
        ngl_l = b%ngl;  nint_l = nint

        !$acc parallel loop gang vector collapse(2) present(ref%face_pack_list, G%face, G%intma, mf%imapl, a, sbuf) &
        !$acc   firstprivate(nbv_l, ngl_l, nint_l, nc) private(iface, IL, base, j, ic)
        do kk = 1, nbv_l
            do n = 1, ngl_l
                iface = ref%face_pack_list(kk)
                IL    = G%intma(mf%imapl(1,n,1,iface), mf%imapl(2,n,1,iface), mf%imapl(3,n,1,iface), G%face(7,iface))
                !$acc loop seq
                do j = 1, nint_l
                    base = (((kk-1)*ngl_l + (n-1))*nint_l + (j-1))*nc
                    !$acc loop seq
                    do ic = 1, nc
                        sbuf(base+ic) = a(ic,IL,j)
                    end do
                end do
            end do
        end do
        !$acc end parallel loop

        if (.not. allocated(ireq)) allocate(ireq(2*par%num_nbh))
        iend = 0
        jj   = 1
        !$acc host_data use_device(sbuf, rbuf)
        do inbh = 1, par%num_nbh
            nqp = 0
            do ib = 1, par%num_send_recv(inbh)
                do i = 1, par%nbh_send_recv_multi(jj)
                    nqp = nqp + ngl_l*nint_l*nc
                end do
                jj = jj + 1
            end do
            istart = iend + 1
            iend   = istart + nqp - 1
            if (nqp > 0) then
                call mpi_irecv(rbuf(istart:iend), nqp, MPI_PRECISION, par%nbh_proc(inbh)-1, ape_tag, &
                               mpi_comm_world, ireq(nreq+1), ierr)
                call mpi_isend(sbuf(istart:iend), nqp, MPI_PRECISION, par%nbh_proc(inbh)-1, ape_tag, &
                               mpi_comm_world, ireq(nreq+2), ierr)
            else
                ireq(nreq+1) = MPI_REQUEST_NULL
                ireq(nreq+2) = MPI_REQUEST_NULL
            end if
            nreq = nreq + 2
        end do
        !$acc end host_data

    end subroutine exchange_start

    subroutine exchange_finish()

        use mpi

        implicit none

        integer :: ierr

        if (nreq > 0) call mpi_waitall(nreq, ireq, MPI_STATUSES_IGNORE, ierr)
        nreq = 0

    end subroutine exchange_finish

end module mod_ape
