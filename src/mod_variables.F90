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
      real, dimension(:),       allocatable :: ope_ave, H_ave
      real, dimension(:),       allocatable :: Qu_ave, Qv_ave, Quv_ave, ope2_ave

      ! 2-D volume arrays
      real, dimension(:,:),     allocatable :: tau_wind, tau_bot
      real, dimension(:,:),     allocatable :: btp_mass_flux_ave, uvb_ave, uvb_ave_df
      real, dimension(:,:),     allocatable :: btp_dpp_graduv, btp_dpp_uvp, graduvb_ave
      real, dimension(:,:),     allocatable :: tau_wind_ave, tau_bot_ave

      ! Face / edge arrays
      real, dimension(:,:),     allocatable :: one_plus_eta_edge, one_plus_eta_edge_2
      real, dimension(:,:),     allocatable :: one_plus_eta_edge_2_ave, H_face_ave

      real, dimension(:,:,:),   allocatable :: btp_mass_flux_face_ave
      real, dimension(:,:,:),   allocatable :: ope_face_ave, ope2_face_ave
      real, dimension(:,:,:),   allocatable :: Qu_face_ave, Qv_face_ave, Quv_face_ave

      real, dimension(:,:,:,:), allocatable :: uvb_face_ave
      real, dimension(:,:,:,:), allocatable :: btp_graduv_dpp_face, graduvb_face_ave

      ! RHS and state buffers
      real, dimension(:,:),     allocatable :: rhs_btp       ! (3,npoin)
      real, dimension(:,:),     allocatable :: rhs_btp_visc  ! (2,npoin)
      real, dimension(:,:),     allocatable :: qb_df         ! (nvar_btp,npoin)

      ! BTP time integrator work arrays (persistent GPU allocations)
      real, dimension(:,:),     allocatable :: qb0_df        ! (nvar_btp,npoin)
      real, dimension(:,:),     allocatable :: qb2_df        ! (nvar_btp,npoin)

      ! Precomputed BCL layer integrals at quad points — updated once per RK stage
      ! by btp_bcl_coeffs_qdf; fixed during BTP subcycling.
      real, dimension(:),       allocatable :: bcl_H         ! (npoin_q) baroclinic H_b
      real, dimension(:),       allocatable :: bcl_uu        ! (npoin_q) sum(pp*u*u)
      real, dimension(:),       allocatable :: bcl_uv        ! (npoin_q) sum(pp*u*v)
      real, dimension(:),       allocatable :: bcl_vv        ! (npoin_q) sum(pp*v*v)
      real, dimension(:),       allocatable :: bcl_dpq       ! (npoin_q) bottom-layer dp
      real, dimension(:),       allocatable :: bcl_up_dpq    ! (npoin_q) bottom-layer u*dp
      real, dimension(:),       allocatable :: bcl_vp_dpq    ! (npoin_q) bottom-layer v*dp
      real, dimension(:),       allocatable :: pbq           ! (npoin_q) background pressure

   end type btp_CS

   !  Baroclinic control structure
   type, public :: bcl_CS

      ! Viscosity / gradient arrays
      real, dimension(:,:,:),     allocatable :: dpp_uvp, dpp_graduv
      real, dimension(:,:,:,:,:), allocatable :: graduv_dpp_face
      real, dimension(:,:),     allocatable :: dpprime_visc, dpprime_visc_q

      ! BCL mass-flux accumulators
      real, dimension(:,:),   allocatable :: sum_layer_mass_flux
      real, dimension(:,:,:), allocatable :: sum_layer_mass_flux_face

      ! State vectors
      real, dimension(:,:,:), allocatable :: q_df       ! (nvar_bcl,npoin,nlayers)
      real, dimension(:,:,:), allocatable :: qprime_df  ! (nvar_bcl,npoin,nlayers)

      ! RHS buffers (persistent device allocations)
      real, dimension(:,:,:), allocatable :: rhs_bcl      ! (nvar_bcl,npoin,nlayers)
      real, dimension(:,:,:), allocatable :: rhs_visc_bcl ! (2,npoin,nlayers)

      ! Time integrator work arrays (persistent GPU allocations, no per-call malloc)
      real, dimension(:,:,:), allocatable :: q0_df    ! (nvar_bcl,npoin,nlayers)
      real, dimension(:,:,:), allocatable :: q1_df    ! (nvar_bcl,npoin,nlayers) rk3 only
      real, dimension(:,:,:), allocatable :: uv_df    ! (2,npoin,nlayers)
      real, dimension(:,:),   allocatable :: qbp_df   ! (nvar_btp,npoin)

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
            btp%pbprime_visc,                                              &
            btp%ope_ave,            btp%H_ave,                             &
            btp%Qu_ave,             btp%Qv_ave,                            &
            btp%Quv_ave,            btp%ope2_ave,                          &
            btp%tau_wind,           btp%tau_bot,                           &
            btp%btp_mass_flux_ave,  btp%uvb_ave,       btp%uvb_ave_df,    &
            btp%btp_dpp_graduv,     btp%btp_dpp_uvp,   btp%graduvb_ave,   &
            btp%tau_wind_ave,       btp%tau_bot_ave,                       &
            btp%one_plus_eta_edge,  btp%one_plus_eta_edge_2,               &
            btp%one_plus_eta_edge_2_ave, btp%H_face_ave,                   &
            btp%btp_mass_flux_face_ave,                                    &
            btp%ope_face_ave,       btp%ope2_face_ave,                     &
            btp%Qu_face_ave,        btp%Qv_face_ave,   btp%Quv_face_ave,  &
            btp%uvb_face_ave,                                              &
            btp%btp_graduv_dpp_face, btp%graduvb_face_ave,                 &
            btp%rhs_btp,            btp%rhs_btp_visc,       btp%qb_df,  &
            btp%qb0_df,             btp%qb2_df,                         &
            btp%bcl_H,              btp%bcl_uu,    btp%bcl_uv,          &
            btp%bcl_vv,             btp%bcl_dpq,   btp%bcl_up_dpq,      &
            btp%bcl_vp_dpq,         btp%pbq)
      end if

      ! 1-D volume arrays
      allocate(                                                               &
         btp%one_plus_eta(G%npoin_q),                                         &
         btp%one_plus_eta_df(G%npoin),                                        &
         btp%ope2_ave_df(G%npoin),                                            &
         btp%one_plus_eta_out(G%npoin),                                       &
         btp%pbprime_visc(G%npoin),                                           &
         btp%ope_ave(G%npoin_q),      btp%H_ave(G%npoin_q),                     &
         btp%Qu_ave(G%npoin_q),       btp%Qv_ave(G%npoin_q),                    &
         btp%Quv_ave(G%npoin_q),      btp%ope2_ave(G%npoin_q),                  &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (btp 1-D)"

      ! 2-D volume arrays
      allocate(                                                               &
         btp%tau_wind(2,G%npoin_q),                                           &
         btp%tau_bot(2,G%npoin_q),                                            &
         btp%btp_mass_flux_ave(2,G%npoin_q),                                  &
         btp%uvb_ave(2,G%npoin_q),                                            &
         btp%uvb_ave_df(2,G%npoin),                                           &
         btp%btp_dpp_graduv(4,G%npoin),                                       &
         btp%btp_dpp_uvp(2,G%npoin),                                          &
         btp%graduvb_ave(4,G%npoin),                                          &
         btp%tau_wind_ave(2,G%npoin_q),                                       &
         btp%tau_bot_ave(2,G%npoin_q),                                        &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (btp 2-D)"

      ! Face / edge arrays
      allocate(                                                               &
         btp%one_plus_eta_edge(b%nq,G%nface),                                   &
         btp%one_plus_eta_edge_2(b%nq,G%nface),                                 &
         btp%one_plus_eta_edge_2_ave(b%nq,G%nface),                             &
         btp%H_face_ave(b%nq,G%nface),                                          &
         btp%btp_mass_flux_face_ave(2,b%nq,G%nface),                            &
         btp%ope_face_ave(2,b%nq,G%nface),                                      &
         btp%ope2_face_ave(2,b%nq,G%nface),                                     &
         btp%Qu_face_ave(2,b%nq,G%nface),                                       &
         btp%Qv_face_ave(2,b%nq,G%nface),                                       &
         btp%Quv_face_ave(2,b%nq,G%nface),                                      &
         btp%uvb_face_ave(2,2,b%nq,G%nface),                                    &
         btp%btp_graduv_dpp_face(5,2,b%ngl,G%nface),                            &
         btp%graduvb_face_ave(4,2,b%ngl,G%nface),                               &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (btp face)"

      ! RHS and state buffers
      allocate(                                                               &
         btp%rhs_btp(3,G%npoin),                                              &
         btp%rhs_btp_visc(2,G%npoin),                                         &
         btp%qb_df(inp%nvar_btp,G%npoin),                                     &
         btp%qb0_df(inp%nvar_btp,G%npoin),                                    &
         btp%qb2_df(inp%nvar_btp,G%npoin),                                    &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (btp buffers)"

      ! Precomputed BCL layer integrals at quad points
      allocate(                                                               &
         btp%bcl_H(G%npoin_q),      btp%bcl_uu(G%npoin_q),                    &
         btp%bcl_uv(G%npoin_q),    btp%bcl_vv(G%npoin_q),                    &
         btp%bcl_dpq(G%npoin_q),    btp%bcl_up_dpq(G%npoin_q),               &
         btp%bcl_vp_dpq(G%npoin_q), btp%pbq(G%npoin_q),                      &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (btp bcl precomp)"

      ! =========================================================
      !  bcl_CS
      ! =========================================================

      if (allocated(bcl%dpp_uvp)) then
         deallocate(                                                         &
            bcl%dpp_uvp,              bcl%dpp_graduv,                      &
            bcl%graduv_dpp_face,                                           &
            bcl%sum_layer_mass_flux,  bcl%sum_layer_mass_flux_face,        &
            bcl%q_df,                 bcl%qprime_df,                     &
            bcl%dpprime_visc,         bcl%dpprime_visc_q,                 &
            bcl%rhs_bcl,              bcl%rhs_visc_bcl,                  &
            bcl%q0_df,                bcl%q1_df,                         &
            bcl%uv_df,                bcl%qbp_df)
      end if

      allocate(                                                               &
         bcl%dpp_uvp(2,G%npoin,inp%nlayers),                                      &
         bcl%dpp_graduv(4,G%npoin,inp%nlayers),                                   &
         bcl%graduv_dpp_face(5,2,b%ngl,G%nface,inp%nlayers),                        &
         bcl%sum_layer_mass_flux(2,G%npoin_q),                                &
         bcl%sum_layer_mass_flux_face(2,b%nq,G%nface),                          &
         bcl%q_df(inp%nvar_bcl,G%npoin,inp%nlayers),                              &
         bcl%qprime_df(inp%nvar_bcl,G%npoin,inp%nlayers),                         &
         bcl%dpprime_visc(G%npoin,inp%nlayers),                             &
         bcl%dpprime_visc_q(G%npoin_q,inp%nlayers),                          &
         bcl%rhs_bcl(inp%nvar_bcl,G%npoin,inp%nlayers),                          &
         bcl%rhs_visc_bcl(2,G%npoin,inp%nlayers),                               &
         bcl%q0_df(inp%nvar_bcl,G%npoin,inp%nlayers),                            &
         bcl%q1_df(inp%nvar_bcl,G%npoin,inp%nlayers),                            &
         bcl%uv_df(2,G%npoin,inp%nlayers),                                       &
         bcl%qbp_df(inp%nvar_btp,G%npoin),                                       &
         stat=stat)
      if (stat /= 0) stop "** Not Enough Memory – mod_allocate_mlswe (bcl)"

   end subroutine mod_allocate_mlswe

end module mod_variables
