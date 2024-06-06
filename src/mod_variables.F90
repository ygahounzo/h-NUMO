module mod_variables

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
                dpp_uvp,dpp_graduv, graduv_dpp_face, btp_graduv_dpp_face, graduvb_face_ave, graduvb_ave, dpprime_visc_q

    private 
    ! module variable and parameters 
    real, dimension(:,:,:), allocatable :: Qu_face, Qv_face, one_plus_eta_face, flux_edge
    real, dimension(:,:,:), allocatable :: Q_uu_dp_edge, Q_uv_dp_edge, Q_vv_dp_edge, H_bcl_edge
    real, dimension(:), allocatable :: Quu, Qvv, Quv, H, one_plus_eta, one_plus_eta_df, ope_ave_df, one_plus_eta_out, pbprime_visc
    real, dimension(:), allocatable :: Q_uu_dp, Q_uv_dp, Q_vv_dp, H_bcl, ope_ave, H_ave, Qu_ave, Qv_ave, Quv_ave, ope2_ave
    real, dimension(:,:), allocatable :: tau_bot, btp_mass_flux, H_face, one_plus_eta_edge_2, btp_mass_flux_ave, uvb_ave, uvb_ave_df, dpprime_visc, dpprime_visc_q

    real, dimension(:,:), allocatable :: one_plus_eta_edge_2_ave, H_face_ave, tau_wind_ave, tau_bot_ave, one_plus_eta_edge
    real, dimension(:,:,:,:), allocatable :: uvb_face_ave, btp_graduv_dpp_face, graduvb_face_ave
    real, dimension(:,:,:), allocatable :: btp_mass_flux_face_ave, ope_face_ave, Qu_face_ave, Qv_face_ave, Quv_face_ave

    real, dimension(:,:,:), allocatable :: dpp_uvp, dpp_graduv
    real, dimension(:,:), allocatable :: btp_dpp_graduv, btp_dpp_uvp, graduvb_ave
    real, dimension(:,:,:,:,:), allocatable :: graduv_dpp_face


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
            uvb_face_ave, btp_mass_flux_face_ave, ope_face_ave, Qu_face_ave, Qv_face_ave, Quv_face_ave, uvb_ave_df)
        endif 

        allocate(Q_uu_dp(npoin_q), Q_uv_dp(npoin_q), Q_vv_dp(npoin_q), H_bcl(npoin_q), &
            Q_uu_dp_edge(2,nq,nface), Q_uv_dp_edge(2,nq,nface), Q_vv_dp_edge(2,nq,nface), H_bcl_edge(2,nq,nface), &
            Qu_face(2,nq,nface), Qv_face(2,nq,nface), one_plus_eta_face(2,nq,nface), flux_edge(2,nq,nface), one_plus_eta_df(npoin), &
            tau_bot(2,npoin_q), btp_mass_flux(2,npoin_q), H_face(nq,nface), one_plus_eta_edge_2(nq,nface), &
            Quu(npoin_q), Qvv(npoin_q), Quv(npoin_q), H(npoin_q), one_plus_eta(npoin_q), &
            ope_ave(npoin_q), H_ave(npoin_q), Qu_ave(npoin_q), Qv_ave(npoin_q), Quv_ave(npoin_q), ope2_ave(npoin_q), &
            btp_mass_flux_ave(2,npoin_q), uvb_ave(2,npoin_q), ope_ave_df(npoin), uvb_face_ave(2,2,nq,nface), &
            btp_mass_flux_face_ave(2,nq,nface), ope_face_ave(2,nq,nface), H_face_ave(nq,nface), &
            Qu_face_ave(2,nq,nface), Qv_face_ave(2,nq,nface), Quv_face_ave(2,nq,nface), &
            one_plus_eta_out(npoin), tau_wind_ave(2,npoin_q), tau_bot_ave(2,npoin_q), one_plus_eta_edge(nq,nface), &
            one_plus_eta_edge_2_ave(nq,nface), uvb_ave_df(2,npoin), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - mod_variables"

        if(allocated(dpprime_visc)) then 
            deallocate(dpprime_visc,pbprime_visc,btp_dpp_graduv, btp_dpp_uvp, dpp_uvp, dpp_graduv, graduv_dpp_face, btp_graduv_dpp_face, &
            graduvb_ave, graduvb_face_ave, dpprime_visc_q)
        endif 

        allocate(dpprime_visc(npoin,nlayers), pbprime_visc(npoin), btp_dpp_graduv(4,npoin),btp_dpp_uvp(2,npoin), dpp_uvp(2,npoin,nlayers),&
            dpp_graduv(4,npoin,nlayers), graduv_dpp_face(5,2,ngl,nface,nlayers),btp_graduv_dpp_face(5,2,ngl,nface), &
            graduvb_face_ave(4,2,ngl,nface), graduvb_ave(4,npoin), dpprime_visc_q(npoin_q,nlayers), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - mod_variables"
        
    end subroutine mod_allocate_mlswe

end module mod_variables