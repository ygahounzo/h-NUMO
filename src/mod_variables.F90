module mod_variables

   ! Pre-allocated variables for barotropic and baroclinic equation terms.
   ! Variables are grouped into two control-structure types:
   !   btp_CS  –  barotropic fields and time-averaged accumulators
   !   bcl_CS  –  baroclinic / BCL fields

   use mod_grid
   use mod_basis
   use mod_input

   implicit none
   private

   !  Barotropic control structure
   type, public :: btp_CS

      ! 1-D arrays
      real, dimension(:),       allocatable :: one_plus_eta, one_plus_eta_df
      real, dimension(:),       allocatable :: ope2_ave_df, one_plus_eta_out
      real, dimension(:),       allocatable :: pbprime_visc
      real, dimension(:),       allocatable :: pbprime_visc_scaled  ! A_H*pbprime_visc, DG-correct viscous flux (see btp_create_laplacian)
      real, dimension(:),       allocatable :: ope_ave, H_ave
      real, dimension(:),       allocatable :: ope2_ave

      ! 2-D volume arrays
      real, dimension(:,:),     allocatable :: Qu_ave  ! (nvel_btp,npoin_q): u-mom eqn flux row (x,y[,z])
      real, dimension(:,:),     allocatable :: Qv_ave  ! (nvel_btp,npoin_q): v-mom eqn flux row (x,y[,z])
      real, dimension(:,:),     allocatable :: Qw_ave  ! (nvel_btp,npoin_q): w-mom eqn flux row, sphere-only
      real, dimension(:,:),     allocatable :: tau_wind, tau_bot
      real, dimension(:,:),     allocatable :: btp_mass_flux_ave, uvb_ave, uvb_ave_df
      real, dimension(:,:),     allocatable :: btp_dpp_uvp, graduvb_ave
      real, dimension(:,:),     allocatable :: btp_dpp_graduvw  ! (ngraduvw_var,npoin): du_dx,du_dy,dv_dx,dv_dy [,graduvz(2)=du_dz,dv_dz, gradw(3)=dw_dx,dw_dy,dw_dz on sphere]
      real, dimension(:,:),     allocatable :: btp_dpp_graduvw_scaled  ! A_H*btp_dpp_graduvw, DG-correct viscous flux (see btp_create_laplacian)
      real, dimension(:,:),     allocatable :: tau_wind_ave, tau_bot_ave

      ! Face / edge arrays
      real, dimension(:,:),     allocatable :: one_plus_eta_edge, one_plus_eta_edge_2
      real, dimension(:,:),     allocatable :: one_plus_eta_edge_2_ave, H_face_ave

      real, dimension(:,:,:),   allocatable :: btp_mass_flux_face_ave
      real, dimension(:,:,:),   allocatable :: ope_face_ave, ope2_face_ave
      real, dimension(:,:,:),   allocatable :: Qu_face_ave, Qv_face_ave, Qw_face_ave

      real, dimension(:,:,:,:), allocatable :: uvb_face_ave
      real, dimension(:,:,:,:), allocatable :: btp_graduv_dpp_face, graduvb_face_ave

      ! RHS and state buffers
      real, dimension(:,:),     allocatable :: rhs_btp       ! (3,npoin)
      real, dimension(:,:),     allocatable :: rhs_btp_visc  ! (nvar_btp-2,npoin): u,v[,w]
      real, dimension(:),       allocatable :: nu_smag        ! (npoin): Smagorinsky eddy viscosity
      real, dimension(:,:),     allocatable :: qb_df         ! (nvar_btp,npoin)

      ! BTP time integrator work arrays (persistent GPU allocations)
      real, dimension(:,:),     allocatable :: qb0_df        ! (nvar_btp,npoin)
      real, dimension(:,:),     allocatable :: qb2_df        ! (nvar_btp,npoin)

      ! Precomputed BCL layer integrals at quad points — updated once per RK stage
      ! by btp_bcl_coeffs_qdf; fixed during BTP subcycling.
      real, dimension(:),       allocatable :: bcl_H         ! (npoin_q) baroclinic H_b
      real, dimension(:,:,:),   allocatable :: bcl_flux      ! (n,n,npoin_q) sum(pp*u_i*u_j), n=nvar_bcl-1
      real, dimension(:,:),     allocatable :: bcl_btp_flux  ! (nvar_bcl,npoin_q) bottom-layer (dp,u*dp,v*dp[,w*dp])
      real, dimension(:),       allocatable :: pbq           ! (npoin_q) background pressure

      ! Weighted flux/source coefficients of the BTP volume integral at each quad
      ! point, written and read within create_rhs_btp_volume_qdf_qp_gpu.
      ! F_vol(Iq,r,j): equation row r (mass, u, v, w); coefficient of j = 1 psi,
      ! 2 dpsi/dxi, 3 dpsi/deta (metric terms already applied). Quad-point index
      ! first so that neighbouring GPU threads access neighbouring addresses.
      real, dimension(:,:,:),   allocatable :: F_vol         ! (npoin_q,4,3)

      ! Per (element, node row m, quad column iqx) work array of
      ! create_rhs_btp_volume_qdf_qp_gpu (index iqx + (m-1)*nqx + (ie-1)*nqx*ngly):
      ! columns 1:nvar_btp hold qb_df summed along xi (first pass of the
      ! interpolation, kernels A0/A); after kernel A, columns 1:8 hold the eta
      ! sums of the test-function pass (kernels B1/B2).
      real, dimension(:,:),     allocatable :: qb_T          ! (nelem*ngly*nqx,max(nvar_btp,8))

      ! Weighted Riemann fluxes (mass, u, v, w) at each face quad point,
      ! index iquad + (iface-1)*nq: written by the flux kernel and scattered
      ! to rhs_btp by a second kernel in create_btp_fluxes_qdf_gpu.
      real, dimension(:,:),     allocatable :: flux_face_q   ! (nq*nface,4)

   end type btp_CS

   !  Baroclinic control structure
   type, public :: bcl_CS

      ! Viscosity / gradient arrays
      real, dimension(:,:,:),     allocatable :: dpp_uvp
      real, dimension(:,:,:),     allocatable :: dpp_graduvw  ! (ngraduvw_var,npoin,nlayers): du_dx,du_dy,dv_dx,dv_dy [,graduvz(2)=du_dz,dv_dz, gradw(3)=dw_dx,dw_dy,dw_dz on sphere]
      real, dimension(:,:,:),     allocatable :: dpp_graduvw_scaled  ! A_H*dpp_graduvw, DG-correct viscous flux (see bcl_create_laplacian)
      real, dimension(:,:,:,:,:), allocatable :: graduv_dpp_face
      real, dimension(:,:),     allocatable :: dpprime_visc, dpprime_visc_q
      real, dimension(:,:),     allocatable :: dpprime_visc_scaled  ! A_H*dpprime_visc, DG-correct viscous flux (see bcl_create_laplacian)

      ! BCL mass-flux accumulators
      real, dimension(:,:),   allocatable :: sum_layer_mass_flux
      real, dimension(:,:,:), allocatable :: sum_layer_mass_flux_face

      ! State vectors
      real, dimension(:,:,:), allocatable :: q_df       ! (nvar_bcl,npoin,nlayers)
      real, dimension(:,:,:), allocatable :: qprime_df  ! (nvar_bcl,npoin,nlayers)

      ! RHS buffers (persistent device allocations)
      real, dimension(:,:,:), allocatable :: rhs_bcl      ! (nvar_bcl,npoin,nlayers)
      real, dimension(:,:,:), allocatable :: rhs_visc_bcl ! (nvar_bcl-1,npoin,nlayers): u,v[,w]

      ! Time integrator work arrays (persistent GPU allocations, no per-call malloc)
      real, dimension(:,:,:), allocatable :: q0_df    ! (nvar_bcl,npoin,nlayers)
      real, dimension(:,:,:), allocatable :: q1_df    ! (nvar_bcl,npoin,nlayers) rk3 only
      real, dimension(:,:,:),   allocatable :: uv_df    ! (2,npoin,nlayers)
      real, dimension(:,:),     allocatable :: qbp_df   ! (nvar_btp,npoin)
      integer, dimension(:,:),  allocatable :: dry_flg  ! (nelem,nlayers): 0=wet,1=mixed,2=dry
      real,    dimension(:,:),  allocatable :: nu_smag  ! (npoin,nlayers): Smagorinsky eddy viscosity
      ! Weighted flux/source coefficients of the sphere BCL volume integral,
      ! written and read within create_rhs_dynamics_volume_bcl_sphere (GPU).
      ! F_bcl(Iq,r,j,k): equation row r (mass, u, v, w) of layer k; j = 1 psi,
      ! 2 dpsi/dxi, 3 dpsi/deta (metric terms already applied).
      real,    dimension(:,:,:,:), allocatable :: F_bcl ! (npoin_q,4,3,nlayers)

      ! Per-element layer work arrays of extract_velocity (GPU-thread index,
      ! the element, first): element-mean dp / momentum, element max / min dp
      ! (ev_min is overwritten by the pointwise blending weight), tridiagonal
      ! diagonal / upper diagonal, rhs overwritten by the solution.
      real, dimension(:,:),   allocatable :: ev_dp, ev_max, ev_min, ev_bc, ev_ctd ! (nelem,nlayers)
      real, dimension(:,:,:), allocatable :: ev_mdp, ev_r                         ! (nelem,3,nlayers)

      ! Per-(face, quad point) layer work arrays of Apply_bcl_fluxes_sphere
      ! (thread index ifq = iquad + (iface-1)*nq first):
      !   sf_z (ifq,side,k+1) face interface height; sf_pep/sf_pem, sf_zep/sf_zem
      !   (ifq,k+1) edge pressure / height, plus and minus side; sf_hl/sf_hr
      !   (ifq,k) clamped left/right thickness; sf_mdl/sf_mdr (ifq,c,k) left/
      !   right momentum*dp; sf_H (ifq,side,k) pressure term; sf_fdp, sf_fu,
      !   sf_fv, sf_fw (ifq,dir,k) upwind mass / momentum fluxes.
      real, dimension(:,:,:), allocatable :: sf_z, sf_mdl, sf_mdr, sf_H
      real, dimension(:,:,:), allocatable :: sf_fdp, sf_fu, sf_fv, sf_fw
      real, dimension(:,:),   allocatable :: sf_pep, sf_pem, sf_zep, sf_zem, sf_hl, sf_hr

   end type bcl_CS

   public :: mod_allocate_mlswe

contains

   subroutine mod_allocate_mlswe(inp, G, b, btp, bcl)

      type(input), intent(in) :: inp
      type(grid), intent(in) :: G
      type(basis), intent(in) :: b

      type(btp_CS), intent(inout) :: btp
      type(bcl_CS), intent(inout) :: bcl
      integer :: stat

      ! =========================================================
      !  btp_CS
      ! =========================================================

      if (allocated(btp%one_plus_eta)) then
         deallocate(                                                         &
            btp%one_plus_eta,       btp%one_plus_eta_df,                   &
            btp%ope2_ave_df,        btp%one_plus_eta_out,                  &
            btp%pbprime_visc,       btp%pbprime_visc_scaled,               &
            btp%ope_ave,            btp%H_ave,                             &
            btp%Qu_ave,             btp%Qv_ave,        btp%Qw_ave,        &
            btp%ope2_ave,                                                  &
            btp%tau_wind,           btp%tau_bot,                           &
            btp%btp_mass_flux_ave,  btp%uvb_ave,       btp%uvb_ave_df,    &
            btp%btp_dpp_graduvw,    btp%btp_dpp_graduvw_scaled,            &
            btp%btp_dpp_uvp,        btp%graduvb_ave,                       &
            btp%tau_wind_ave,       btp%tau_bot_ave,                       &
            btp%one_plus_eta_edge,  btp%one_plus_eta_edge_2,               &
            btp%one_plus_eta_edge_2_ave, btp%H_face_ave,                   &
            btp%btp_mass_flux_face_ave,                                    &
            btp%ope_face_ave,       btp%ope2_face_ave,                     &
            btp%Qu_face_ave,        btp%Qv_face_ave,   btp%Qw_face_ave,   &
            btp%uvb_face_ave,                                              &
            btp%btp_graduv_dpp_face, btp%graduvb_face_ave,                 &
            btp%rhs_btp,            btp%rhs_btp_visc,       btp%nu_smag, &
            btp%qb_df,  &
            btp%qb0_df,             btp%qb2_df,                         &
            btp%bcl_H,              btp%bcl_flux,  btp%bcl_btp_flux,    &
            btp%pbq)
      end if

      ! 1-D volume arrays
      allocate(                                                               &
         btp%one_plus_eta(G%npoin_q),                                         &
         btp%one_plus_eta_df(G%npoin),                                        &
         btp%ope2_ave_df(G%npoin),                                            &
         btp%one_plus_eta_out(G%npoin),                                       &
         btp%pbprime_visc(G%npoin),                                           &
         btp%pbprime_visc_scaled(G%npoin),                                    &
         btp%ope_ave(G%npoin_q),      btp%H_ave(G%npoin_q),                     &
         btp%ope2_ave(G%npoin_q),                                                &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (btp 1-D)"

      ! 2-D volume arrays
      allocate(                                                               &
         btp%Qu_ave(inp%nvar_btp-2,G%npoin_q), btp%Qv_ave(inp%nvar_btp-2,G%npoin_q), &
         btp%Qw_ave(inp%nvar_btp-2,G%npoin_q),                                &
         btp%tau_wind(2,G%npoin_q),                                           &
         btp%tau_bot(2,G%npoin_q),                                            &
         btp%btp_mass_flux_ave(inp%nvar_btp-2,G%npoin_q),                     &
         btp%uvb_ave(inp%nvar_btp-2,G%npoin_q),                               &
         btp%uvb_ave_df(inp%nvar_btp-2,G%npoin),                              &
         btp%btp_dpp_graduvw(inp%ngraduvw_var,G%npoin),                       &
         btp%btp_dpp_graduvw_scaled(inp%ngraduvw_var,G%npoin),                &
         btp%btp_dpp_uvp(2,G%npoin),                                          &
         btp%graduvb_ave(inp%ngraduvw_var,G%npoin),                          &
         btp%tau_wind_ave(2,G%npoin_q),                                       &
         btp%tau_bot_ave(inp%nvar_btp-2,G%npoin_q),                           &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (btp 2-D)"

      ! Face / edge arrays
      allocate(                                                               &
         btp%one_plus_eta_edge(b%nq,G%nface),                                   &
         btp%one_plus_eta_edge_2(b%nq,G%nface),                                 &
         btp%one_plus_eta_edge_2_ave(b%nq,G%nface),                             &
         btp%H_face_ave(b%nq,G%nface),                                          &
         btp%btp_mass_flux_face_ave(inp%nvar_btp-2,b%nq,G%nface),               &
         btp%ope_face_ave(2,b%nq,G%nface),                                      &
         btp%ope2_face_ave(2,b%nq,G%nface),                                     &
         btp%Qu_face_ave(inp%nvar_btp-2,b%nq,G%nface),                          &
         btp%Qv_face_ave(inp%nvar_btp-2,b%nq,G%nface),                          &
         btp%Qw_face_ave(inp%nvar_btp-2,b%nq,G%nface),                          &
         btp%uvb_face_ave(inp%nvar_btp-2,2,b%nq,G%nface),                       &
         btp%btp_graduv_dpp_face(5,2,b%ngl,G%nface),                            &
         btp%graduvb_face_ave(inp%ngraduvw_var,2,b%ngl,G%nface),               &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (btp face)"

      ! RHS and state buffers
      allocate(                                                               &
         btp%rhs_btp(inp%nvar_btp-1,G%npoin),                                 &
         btp%rhs_btp_visc(inp%nvar_btp-2,G%npoin),                           &
         btp%nu_smag(G%npoin),                                                &
         btp%qb_df(inp%nvar_btp,G%npoin),                                     &
         btp%qb0_df(inp%nvar_btp,G%npoin),                                    &
         btp%qb2_df(inp%nvar_btp,G%npoin),                                    &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (btp buffers)"

      ! Precomputed BCL layer integrals at quad points
      allocate(                                                               &
         btp%bcl_H(G%npoin_q),                                              &
         btp%bcl_flux(inp%nvar_bcl-1,inp%nvar_bcl-1,G%npoin_q),             &
         btp%bcl_btp_flux(inp%nvar_bcl,G%npoin_q),                          &
         btp%pbq(G%npoin_q),                                                &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (btp bcl precomp)"

      ! =========================================================
      !  bcl_CS
      ! =========================================================

      if (allocated(bcl%dpp_uvp)) then
         deallocate(                                                         &
            bcl%dpp_uvp,              bcl%dpp_graduvw,                    &
            bcl%dpp_graduvw_scaled,                                        &
            bcl%graduv_dpp_face,                                           &
            bcl%sum_layer_mass_flux,  bcl%sum_layer_mass_flux_face,        &
            bcl%q_df,                 bcl%qprime_df,                     &
            bcl%dpprime_visc,         bcl%dpprime_visc_q,                 &
            bcl%dpprime_visc_scaled,                                      &
            bcl%rhs_bcl,              bcl%rhs_visc_bcl,                  &
            bcl%q0_df,                bcl%q1_df,                         &
            bcl%uv_df,                bcl%qbp_df,                        &
            bcl%dry_flg,              bcl%nu_smag)
      end if

      allocate(                                                               &
         bcl%dpp_uvp(3,G%npoin,inp%nlayers),                                      &
         bcl%dpp_graduvw(inp%ngraduvw_var,G%npoin,inp%nlayers),                    &
         bcl%dpp_graduvw_scaled(inp%ngraduvw_var,G%npoin,inp%nlayers),             &
         bcl%graduv_dpp_face(5,2,b%ngl,G%nface,inp%nlayers),                        &
         bcl%sum_layer_mass_flux(2,G%npoin_q),                                &
         bcl%sum_layer_mass_flux_face(2,b%nq,G%nface),                          &
         bcl%q_df(inp%nvar_bcl,G%npoin,inp%nlayers),                              &
         bcl%qprime_df(inp%nvar_bcl,G%npoin,inp%nlayers),                         &
         bcl%dpprime_visc(G%npoin,inp%nlayers),                             &
         bcl%dpprime_visc_q(G%npoin_q,inp%nlayers),                          &
         bcl%dpprime_visc_scaled(G%npoin,inp%nlayers),                       &
         bcl%rhs_bcl(inp%nvar_bcl,G%npoin,inp%nlayers),                          &
         bcl%rhs_visc_bcl(inp%nvar_bcl-1,G%npoin,inp%nlayers),                  &
         bcl%q0_df(inp%nvar_bcl,G%npoin,inp%nlayers),                            &
         bcl%q1_df(inp%nvar_bcl,G%npoin,inp%nlayers),                            &
         bcl%uv_df(inp%nvar_bcl-1,G%npoin,inp%nlayers),                          &
         bcl%qbp_df(inp%nvar_btp,G%npoin),                                       &
         bcl%dry_flg(G%nelem,inp%nlayers),                                        &
         bcl%nu_smag(G%npoin,inp%nlayers),                                        &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (bcl)"

#ifdef _OPENACC
      ! GPU-only work arrays of the restructured kernels (not used by the CPU
      ! code paths, so not allocated in CPU builds).
      if (allocated(btp%F_vol)) then
         deallocate(btp%F_vol, btp%qb_T, btp%flux_face_q,                   &
            bcl%F_bcl,                                                    &
            bcl%ev_dp, bcl%ev_max, bcl%ev_min, bcl%ev_bc, bcl%ev_ctd,      &
            bcl%ev_mdp, bcl%ev_r,                                          &
            bcl%sf_z, bcl%sf_mdl, bcl%sf_mdr, bcl%sf_H,                    &
            bcl%sf_fdp, bcl%sf_fu, bcl%sf_fv, bcl%sf_fw,                   &
            bcl%sf_pep, bcl%sf_pem, bcl%sf_zep, bcl%sf_zem,                &
            bcl%sf_hl, bcl%sf_hr)
      end if
      allocate(                                                               &
         btp%F_vol(G%npoin_q,4,3),                                        &
         btp%qb_T(G%nelem*b%ngly*b%nqx,max(inp%nvar_btp,8)),                &
         btp%flux_face_q(b%nq*G%nface,4),                                   &
         bcl%F_bcl(G%npoin_q,4,3,inp%nlayers),                                   &
         bcl%ev_dp(G%nelem,inp%nlayers),   bcl%ev_max(G%nelem,inp%nlayers),      &
         bcl%ev_min(G%nelem,inp%nlayers),  bcl%ev_bc(G%nelem,inp%nlayers),       &
         bcl%ev_ctd(G%nelem,inp%nlayers),                                        &
         bcl%ev_mdp(G%nelem,3,inp%nlayers), bcl%ev_r(G%nelem,3,inp%nlayers),     &
         bcl%sf_z(b%nq*G%nface,2,inp%nlayers+1),                                 &
         bcl%sf_pep(b%nq*G%nface,inp%nlayers+1), bcl%sf_pem(b%nq*G%nface,inp%nlayers+1), &
         bcl%sf_zep(b%nq*G%nface,inp%nlayers+1), bcl%sf_zem(b%nq*G%nface,inp%nlayers+1), &
         bcl%sf_hl(b%nq*G%nface,inp%nlayers),   bcl%sf_hr(b%nq*G%nface,inp%nlayers),     &
         bcl%sf_mdl(b%nq*G%nface,3,inp%nlayers), bcl%sf_mdr(b%nq*G%nface,3,inp%nlayers), &
         bcl%sf_H(b%nq*G%nface,2,inp%nlayers),                                   &
         bcl%sf_fdp(b%nq*G%nface,3,inp%nlayers), bcl%sf_fu(b%nq*G%nface,3,inp%nlayers),  &
         bcl%sf_fv(b%nq*G%nface,3,inp%nlayers),  bcl%sf_fw(b%nq*G%nface,3,inp%nlayers),  &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (GPU work arrays)"
#endif
      bcl%nu_smag  = 0.0

   end subroutine mod_allocate_mlswe

end module mod_variables
