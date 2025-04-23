module mod_variables

    ! This routine contains pre-allocation variable for barotropic equations terms 

    use mod_grid, only: npoin_q, nface, npoin, face, intma_dg_quad
    use mod_basis, only: nq, ngl
    use mod_input, only: nlayers

    public :: Qu_face, Qv_face, one_plus_eta_face, flux_edge, &
                Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, H_bcl_edge, &
                Quu, Qvv, Quv, H, one_plus_eta, one_plus_eta_df, ope_ave_df, one_plus_eta_out, &
                Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave, &
                tau_bot, btp_mass_flux, H_face, one_plus_eta_edge_2, btp_mass_flux_ave, uvb_ave, &
                one_plus_eta_edge_2_ave, H_face_ave, tau_wind_ave, tau_bot_ave, &
                uvb_face_ave, btp_mass_flux_face_ave, ope_face_ave, Qu_face_ave, Qv_face_ave, Quv_face_ave, &
                pbprime_visc, one_plus_eta_edge, &
                mod_allocate_mlswe, uvb_ave_df, dpprime_visc,btp_dpp_graduv, btp_dpp_uvp, &
                dpp_uvp,dpp_graduv, graduv_dpp_face, btp_graduv_dpp_face, graduvb_face_ave, graduvb_ave, dpprime_visc_q, &
                Q_uu_dp_df, Q_uv_dp_df, Q_vv_dp_df, H_bcl_df, H_bcl_edge_df, &
                Q_uu_dp_edge_df, Q_uv_dp_edge_df, Q_vv_dp_edge_df, &
                tau_bot_ave_df, H_ave_df, Qu_ave_df, Quv_ave_df, Qv_ave_df, btp_mass_flux_ave_df

    public :: Quu_temp, Qvv_temp, Quv_temp, Q_uu_dp_temp, Q_uv_dp_temp, Q_vv_dp_temp, &
            Qu_ave_temp, Qv_ave_temp, Quv_ave_temp, tau_bot_ave_temp, &
            Q_uu_dp_edge_temp, Q_uv_dp_edge_temp, Q_vv_dp_edge_temp, tau_bot_temp

    public :: H_r,sum_layer_mass_flux, u_udp_temp, v_vdp_temp, p, z_elev, uvdp_temp, sum_layer_mass_flux_face, udp_left, &
        vdp_left, udp_right, vdp_right, u_vdp_temp, grad_z, flux_edge_bcl, u_edge, v_edge, udp_flux_edge, vdp_flux_edge, &
        H_r_face, flux_adjustment, flux_adjust_edge, tau_wind_int, tau_bot_int
        
    public :: one_plus_eta_edge_2_ave_df, H_face_ave_df, btp_mass_flux_face_ave_df, ope_face_ave_df, Qu_face_ave_df, Qv_face_ave_df
    public :: qbface_ave

    private 
    ! module variable and parameters 
    real, dimension(:,:,:), allocatable :: Qu_face, Qv_face, one_plus_eta_face, flux_edge
    real, dimension(:,:), allocatable :: Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, H_bcl_edge
    real, dimension(:), allocatable :: Q_uu_dp_df, Q_uv_dp_df, Q_vv_dp_df, H_bcl_df
    real, dimension(:,:), allocatable :: Q_uu_dp_edge_df, Q_uv_dp_edge_df, Q_vv_dp_edge_df, H_bcl_edge_df
    real, dimension(:), allocatable :: Quu, Qvv, Quv, H, one_plus_eta, one_plus_eta_df, ope_ave_df, one_plus_eta_out, pbprime_visc
    real, dimension(:), allocatable :: Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave
    real, dimension(:,:), allocatable :: tau_bot, btp_mass_flux, H_face, one_plus_eta_edge_2, btp_mass_flux_ave, uvb_ave, uvb_ave_df, dpprime_visc, dpprime_visc_q
    real, dimension(:), allocatable :: H_ave_df, Qu_ave_df, Quv_ave_df, Qv_ave_df
    real, dimension(:,:), allocatable :: btp_mass_flux_ave_df, tau_bot_ave_df

    real, dimension(:,:), allocatable :: one_plus_eta_edge_2_ave, H_face_ave, tau_wind_ave, tau_bot_ave, one_plus_eta_edge
    real, dimension(:,:,:,:), allocatable :: uvb_face_ave, btp_graduv_dpp_face, graduvb_face_ave
    real, dimension(:,:,:), allocatable :: btp_mass_flux_face_ave, ope_face_ave, Qu_face_ave, Qv_face_ave, Quv_face_ave
    real, dimension(:,:,:), allocatable :: btp_mass_flux_face_ave_df, ope_face_ave_df, Qu_face_ave_df, Qv_face_ave_df
    real, dimension(:,:), allocatable :: H_face_ave_df, one_plus_eta_edge_2_ave_df

    real, dimension(:), allocatable :: Quu_temp, Qvv_temp, Quv_temp, Q_uu_dp_temp, Q_uv_dp_temp, Q_vv_dp_temp
    real, dimension(:), allocatable :: Qu_ave_temp, Qv_ave_temp, Quv_ave_temp
    real, dimension(:,:,:), allocatable :: Q_uu_dp_edge_temp, Q_uv_dp_edge_temp, Q_vv_dp_edge_temp

    real, dimension(:,:,:), allocatable :: dpp_uvp, dpp_graduv
    real, dimension(:,:), allocatable :: btp_dpp_graduv, btp_dpp_uvp, graduvb_ave, tau_bot_temp, tau_bot_ave_temp
    real, dimension(:,:,:,:,:), allocatable :: graduv_dpp_face

    ! bcl variables 
    real, dimension(:,:), allocatable :: sum_layer_mass_flux, u_udp_temp, v_vdp_temp, H_r, p, z_elev
    real, dimension(:,:,:), allocatable :: uvdp_temp, sum_layer_mass_flux_face, udp_left, vdp_left, udp_right, vdp_right, &
                                            u_vdp_temp, grad_z, flux_adjustment, tau_wind_int, tau_bot_int
    real, dimension(:,:,:,:), allocatable :: flux_edge_bcl, u_edge, v_edge, udp_flux_edge, vdp_flux_edge, H_r_face, flux_adjust_edge
    real, dimension(:,:,:,:), allocatable :: qbface_ave

    contains

    subroutine mod_allocate_mlswe()

        implicit none
        integer AllocateStatus

        if(allocated(H_bcl)) then 
            deallocate(Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave, &
            Quu, Qvv, Quv, H, one_plus_eta, one_plus_eta_df, ope_ave_df, one_plus_eta_out, &
            Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, H_bcl_edge, &
            Qu_face, Qv_face, one_plus_eta_face, flux_edge, &
            tau_bot, btp_mass_flux, H_face, one_plus_eta_edge_2, btp_mass_flux_ave, uvb_ave, &
            one_plus_eta_edge_2_ave, H_face_ave, tau_wind_ave, tau_bot_ave, one_plus_eta_edge, &
            uvb_face_ave, btp_mass_flux_face_ave, ope_face_ave, Qu_face_ave, Qv_face_ave, Quv_face_ave, uvb_ave_df, &
            Q_uu_dp_df, Q_uv_dp_df, Q_vv_dp_df, H_bcl_df, H_bcl_edge_df, &
            Q_uu_dp_edge_df, Q_uv_dp_edge_df, Q_vv_dp_edge_df, &
            tau_bot_ave_df, H_ave_df, Qu_ave_df, Quv_ave_df, Qv_ave_df, btp_mass_flux_ave_df, &
            one_plus_eta_edge_2_ave_df, H_face_ave_df, btp_mass_flux_face_ave_df, &
            ope_face_ave_df, Qu_face_ave_df, Qv_face_ave_df, qbface_ave)
        endif 

        allocate(Q_uu_dp(npoin_q), Q_uv_dp(npoin_q), Q_vv_dp(npoin_q), H_bcl(npoin_q), &
            Q_uu_dp_edge(nq,nface), Q_uv_dp_edge(nq,nface), Q_vv_dp_edge(nq,nface), H_bcl_edge(nq,nface), &
            Qu_face(2,nq,nface), Qv_face(2,nq,nface), one_plus_eta_face(2,nq,nface), flux_edge(2,nq,nface), one_plus_eta_df(npoin), &
            tau_bot(2,npoin_q), btp_mass_flux(2,npoin_q), H_face(nq,nface), one_plus_eta_edge_2(nq,nface), &
            Quu(npoin_q), Qvv(npoin_q), Quv(npoin_q), H(npoin_q), one_plus_eta(npoin_q), &
            ope_ave(npoin_q), H_ave(npoin_q), Qu_ave(npoin_q), Qv_ave(npoin_q), Quv_ave(npoin_q), ope2_ave(npoin_q), &
            btp_mass_flux_ave(2,npoin_q), uvb_ave(2,npoin_q), ope_ave_df(npoin), uvb_face_ave(2,2,nq,nface), &
            btp_mass_flux_face_ave(2,nq,nface), ope_face_ave(2,nq,nface), H_face_ave(nq,nface), &
            Qu_face_ave(2,nq,nface), Qv_face_ave(2,nq,nface), Quv_face_ave(2,nq,nface), &
            one_plus_eta_out(npoin), tau_wind_ave(2,npoin_q), tau_bot_ave(2,npoin_q), one_plus_eta_edge(nq,nface), &
            one_plus_eta_edge_2_ave(nq,nface), uvb_ave_df(2,npoin), &
            Q_uu_dp_df(npoin), Q_uv_dp_df(npoin), Q_vv_dp_df(npoin), H_bcl_df(npoin), &
            H_bcl_edge_df(ngl,nface), Q_uu_dp_edge_df(ngl,nface), Q_uv_dp_edge_df(ngl,nface), Q_vv_dp_edge_df(ngl,nface), &
            tau_bot_ave_df(2,npoin), H_ave_df(npoin), Qu_ave_df(npoin), Quv_ave_df(npoin), Qv_ave_df(npoin), &
            btp_mass_flux_ave_df(2,npoin), &
            one_plus_eta_edge_2_ave_df(ngl,nface), H_face_ave_df(ngl,nface), btp_mass_flux_face_ave_df(2,ngl,nface), &
            ope_face_ave_df(2,ngl,nface), Qu_face_ave_df(2,ngl,nface), Qv_face_ave_df(2,ngl,nface), qbface_ave(4,2,ngl,nface), &
            stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - mod_variables"

        if(allocated(Q_uu_dp_temp)) then 
            deallocate(Quu_temp, Qvv_temp, Quv_temp, Q_uu_dp_temp, Q_uv_dp_temp, Q_vv_dp_temp, &
            Qu_ave_temp, Qv_ave_temp, Quv_ave_temp, tau_bot_ave_temp, &
            Q_uu_dp_edge_temp, Q_uv_dp_edge_temp, Q_vv_dp_edge_temp, tau_bot_temp)
        endif 

        allocate(Q_uu_dp_temp(npoin), Q_uv_dp_temp(npoin), Q_vv_dp_temp(npoin), &
            Q_uu_dp_edge_temp(2,ngl,nface), Q_uv_dp_edge_temp(2,ngl,nface), Q_vv_dp_edge_temp(2,ngl,nface), &
            tau_bot_temp(2,npoin), Quu_temp(npoin), Qvv_temp(npoin), Quv_temp(npoin), &
            Qu_ave_temp(npoin), Qv_ave_temp(npoin), Quv_ave_temp(npoin), tau_bot_ave_temp(2,npoin), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - mod_variables"

        if(allocated(dpprime_visc)) then 
            deallocate(dpprime_visc,pbprime_visc,btp_dpp_graduv, btp_dpp_uvp, dpp_uvp, dpp_graduv, graduv_dpp_face, btp_graduv_dpp_face, &
            graduvb_ave, graduvb_face_ave, dpprime_visc_q)
        endif 

        allocate(dpprime_visc(npoin,nlayers), pbprime_visc(npoin), btp_dpp_graduv(4,npoin),btp_dpp_uvp(2,npoin), dpp_uvp(2,npoin,nlayers),&
            dpp_graduv(4,npoin,nlayers), graduv_dpp_face(5,2,ngl,nface,nlayers),btp_graduv_dpp_face(5,2,ngl,nface), &
            graduvb_face_ave(4,2,ngl,nface), graduvb_ave(4,npoin), dpprime_visc_q(npoin_q,nlayers), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - mod_variables"



    ! ======== bcl variables =========

    if(allocated(H_r)) then 
        deallocate(H_r,sum_layer_mass_flux, u_udp_temp, v_vdp_temp, p, z_elev, uvdp_temp, sum_layer_mass_flux_face, udp_left, &
        vdp_left, udp_right, vdp_right, u_vdp_temp, grad_z, flux_edge_bcl, u_edge, v_edge, udp_flux_edge, vdp_flux_edge, H_r_face, &
        flux_adjustment, flux_adjust_edge, tau_wind_int, tau_bot_int)
    endif 

    allocate(H_r(npoin_q,nlayers),sum_layer_mass_flux(2,npoin_q), u_udp_temp(npoin_q, nlayers), v_vdp_temp(npoin_q, nlayers), &
        p(npoin_q,nlayers+1), z_elev(npoin,nlayers+1), uvdp_temp(2,npoin_q,nlayers), sum_layer_mass_flux_face(2,nq,nface), udp_left(nq,nface,nlayers), &
        vdp_left(nq,nface,nlayers), udp_right(nq,nface,nlayers), vdp_right(nq,nface,nlayers), u_vdp_temp(2, npoin_q, nlayers), grad_z(2,npoin_q,nlayers+1), &
        flux_edge_bcl(2,nq, nface, nlayers), u_edge(2,nq, nface, nlayers), v_edge(2,nq, nface, nlayers), udp_flux_edge(2, nq, nface, nlayers), &
        vdp_flux_edge(2, nq, nface, nlayers), H_r_face(2,nq,nface,nlayers), flux_adjustment(2,npoin_q, nlayers), &
        flux_adjust_edge(2,nq, nface, nlayers), tau_wind_int(2,npoin_q,nlayers), tau_bot_int(2,npoin_q,nlayers), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop "** Not Enough Memory - mod_variables" 
        
    end subroutine mod_allocate_mlswe

end module mod_variables
