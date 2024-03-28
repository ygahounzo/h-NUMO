subroutine ti_rk35_bcl(q_df,qb_df,qp_df_out)

	use mod_splitting_v2, only: thickness_v2, momentum_v2, momentum_mass_v2
	use mod_input, only: nlayers, dt, dt_btp, dpprime_visc_min, ti_method_btp
	use mod_grid, only: npoin, npoin_q, nface
	use mod_constants, only: gravity
	use mod_initial, only: alpha_mlswe, zbot_df, pbprime_df
	use mod_basis, only: nq
	use mod_rk_mlswe, only: ti_barotropic_ssprk_mlswe, ti_barotropic_rk_mlswe2, ti_barotropic_rk_mlswe3
    use mod_barotropic_terms, only: btp_evaluate_mom_dp, btp_evaluate_mom_dp_face
    use mod_layer_terms, only: evaluate_mom_mass, evaluate_mom_mass_face

	implicit none

	real, dimension(4,npoin), intent(inout) :: qb_df
	real, dimension(3,npoin,nlayers), intent(inout) :: q_df

	real, dimension(5,npoin,nlayers), intent(out) :: qp_df_out
	real, dimension(npoin) :: one_plus_eta

    real, dimension(3,npoin_q,nlayers) :: q
	real, dimension(3,2,nq,nface,nlayers) :: q_face
    real, dimension(4,npoin_q) :: qb
	real, dimension(4,2,nq,nface) :: qb_face
    real, dimension(3,npoin_q,nlayers) :: qprime
	real, dimension(3,2,nq,nface,nlayers) :: qprime_face
    real, dimension(npoin,nlayers) :: dpprime_df2
	real, dimension(3,npoin,nlayers) :: qprime_df

	real, dimension(nq,nface) :: one_plus_eta_edge_2_ave
	real, dimension(npoin_q) :: ope_ave, ope2_ave, H_ave, Qu_ave, Qv_ave, Quv_ave
	real, dimension(2,npoin_q) :: btp_mass_flux_ave, uvb_ave 
	real, dimension(npoin) :: ope_ave_df
	real, dimension(2,2,nq,nface) :: uvb_face_ave
	real, dimension(2,nq,nface) :: btp_mass_flux_face_ave, Qu_face_ave, Qv_face_ave, Quv_face_ave, ope_face_ave, ope2_face_ave
	real, dimension(nq,nface) :: H_face_ave
	real, dimension(2,npoin_q) :: tau_wind_ave, tau_bot_ave
	real, dimension(npoin,nlayers+1) :: mslwe_elevation

	real, dimension(3,npoin_q,nlayers) :: qprime_avg
	real, dimension(3,npoin_q,nlayers) :: qprime2
	real, dimension(3,2,nq,nface,nlayers) :: qprime_face2
	real, dimension(3,2,nq,nface,nlayers) :: qprime_face_avg
	real, dimension(3,npoin_q,nlayers) :: q2
	real, dimension(3,2,nq,nface,nlayers) :: q_face2, q_face1
	real, dimension(4,npoin) :: qbp_df
	real, dimension(3,npoin,nlayers) :: q_df2
	real, dimension(3,npoin,nlayers) :: qprime_df_avg, qprime_df2
	real, dimension(4,npoin_q) :: qbp
	real, dimension(4,2,nq,nface) :: qbp_face

	real, dimension(2,npoin) :: uvb_df_ave
    real, dimension(npoin) :: one_plus_eta_temp

	integer :: k, Iq, iquad, iface, ilr, flag_pred

	! Output variables
	real, dimension(2,nq, nface, nlayers) :: u_edge, v_edge
	real, dimension(2,nq,nface,nlayers) :: dprime_face_corr

	mslwe_elevation = 0.0

    ! Interpolate barotropic dofs to quadrature points
    call btp_evaluate_mom_dp(qb,qb_df)
    call btp_evaluate_mom_dp_face(qb_face, qb)

    ! Interpolate baroclinic dofs to quadrature points

    call evaluate_mom_mass(q,qprime,q_df,qb)
    call evaluate_mom_mass_face(q_face,qprime_face,q,qprime)

    ! Compue dpprime, uprime and vprime at the dofs 

    one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)

    do k = 1,nlayers
        
        qprime_df(1,:,k) = q_df(1,:,k) / one_plus_eta_temp(:)
        qprime_df(2,:,k) = q_df(2,:,k)/q_df(1,:,k) - qb_df(3,:)/qb_df(1,:)
        qprime_df(3,:,k) = q_df(3,:,k)/q_df(1,:,k) - qb_df(4,:)/qb_df(1,:)
    end do 

	! Prediction step

	flag_pred = 1

	!call create_communicator_quad_layer(qprime_face,3,nlayers)
	!call create_communicator_quad_layer(q_face,3,nlayers)

	call create_communicator_quad_layer_all(qprime_face,q_face,3,nlayers)

	qbp = qb
	qbp_face = qb_face
	qbp_df = qb_df

	!call ti_barotropic_ssprk_mlswe(one_plus_eta,one_plus_eta_edge_2_ave,uvb_ave, ope_ave, ope2_ave, &
	!	H_ave, Qu_ave, Qv_ave, Quv_ave, btp_mass_flux_ave, ope_ave_df, uvb_face_ave, &
	!	ope_face_ave, btp_mass_flux_face_ave, H_face_ave, &
	!	Qu_face_ave, Qv_face_ave, Quv_face_ave, tau_wind_ave, tau_bot_ave, qbp,qbp_face,&
	!	qbp_df,qprime,qprime_face, qprime_df, flag_pred, uvb_df_ave)

	call ti_barotropic_rk_mlswe3(one_plus_eta,one_plus_eta_edge_2_ave,uvb_ave, ope_ave, ope2_ave, &
		H_ave, Qu_ave, Qv_ave, Quv_ave, btp_mass_flux_ave, ope_ave_df, uvb_face_ave, &
		ope_face_ave, btp_mass_flux_face_ave, H_face_ave, &
		Qu_face_ave, Qv_face_ave, Quv_face_ave, tau_wind_ave, tau_bot_ave, qbp,qbp_face,&
		qbp_df,qprime,qprime_face, qprime_df, flag_pred, uvb_df_ave)

	qprime2 = qprime
	qprime_face2 = qprime_face
	q_df2 = q_df
	qprime_df2 = qprime_df
	q2 = q
	q_face2 = q_face

	call momentum_mass_v2(q2,qprime2,q_df2,q_face2,qprime_face2,qbp,qbp_face,ope_ave,one_plus_eta_edge_2_ave,uvb_ave,Qu_ave,&
        Qv_ave,Quv_ave,H_ave,uvb_face_ave,ope_face_ave,Qu_face_ave,Qv_face_ave,Quv_face_ave,H_face_ave,&
        qprime_df2,ope_ave_df,tau_bot_ave,tau_wind_ave,qprime_face2,flag_pred,&
        q,q_face,qbp_df,qprime, qprime_face, ope2_ave, btp_mass_flux_ave, btp_mass_flux_face_ave, uvb_df_ave)

	! Correction step

	! Communication of qprime_face2 values within the inter-processor boundary
	call create_communicator_quad_layer(qprime_face2,3,nlayers)

    ! Take the average between the predition and the previous time step values
	qprime_avg = 0.5*(qprime2 + qprime)
	qprime_face_avg = 0.5*(qprime_face2 + qprime_face)
	qprime_df_avg = 0.5*(qprime_df2 + qprime_df)
	
	flag_pred = 0

	!call ti_barotropic_ssprk_mlswe(one_plus_eta,one_plus_eta_edge_2_ave,uvb_ave, ope_ave, ope2_ave, &
	!	H_ave, Qu_ave, Qv_ave, Quv_ave, btp_mass_flux_ave, ope_ave_df, uvb_face_ave, &
	!	ope_face_ave, btp_mass_flux_face_ave, H_face_ave, &
	!	Qu_face_ave, Qv_face_ave, Quv_face_ave, tau_wind_ave, tau_bot_ave, qb,qb_face,&
	!	qb_df, qprime_avg,qprime_face_avg, qprime_df_avg, flag_pred, uvb_df_ave)

	call ti_barotropic_rk_mlswe3(one_plus_eta,one_plus_eta_edge_2_ave,uvb_ave, ope_ave, ope2_ave, &
		H_ave, Qu_ave, Qv_ave, Quv_ave, btp_mass_flux_ave, ope_ave_df, uvb_face_ave, &
		ope_face_ave, btp_mass_flux_face_ave, H_face_ave, &
		Qu_face_ave, Qv_face_ave, Quv_face_ave, tau_wind_ave, tau_bot_ave, qb,qb_face,&
		qb_df, qprime_avg,qprime_face_avg, qprime_df_avg, flag_pred, uvb_df_ave)

	q2 = q
	q_face2 = q_face

	qprime2 = qprime_avg
	qprime_face2 = qprime_face_avg

	!qprime_face_avg(1,:,:,:,:) = qprime_face(1,:,:,:,:) ! No correction for thickness

	call thickness_v2(q,qprime_avg, q_df, q_face, qprime_face_avg, u_edge, v_edge, uvb_ave, btp_mass_flux_ave, ope_ave, uvb_face_ave, ope_face_ave, &
		btp_mass_flux_face_ave, dpprime_df2, flag_pred, qb_df)

	dprime_face_corr = qprime_face_avg(1,:,:,:,:)

	! Communication of qprime_face_avg values within the processor boundary
	call create_communicator_quad_layer(qprime_face_avg(1,:,:,:,:),1,nlayers)

	qprime_df_avg(1,:,:) = 0.5*(qprime_df(1,:,:) + dpprime_df2(:,:))
	qprime_avg(1,:,:) = 0.5*(qprime(1,:,:) + qprime_avg(1,:,:))
	qprime_face_avg(1,:,:,:,:) = 0.5*(qprime_face(1,:,:,:,:) + qprime_face_avg(1,:,:,:,:))

	call momentum_v2(q,qprime_avg,q_df,q_face,qprime_face_avg,qb,qb_face,ope_ave,one_plus_eta_edge_2_ave,uvb_ave,u_edge,v_edge,Qu_ave,&
		Qv_ave,Quv_ave,H_ave,uvb_face_ave,ope_face_ave,Qu_face_ave,Qv_face_ave,Quv_face_ave,H_face_ave,&
		qprime_df_avg,ope_ave_df,tau_bot_ave,tau_wind_ave, qprime_face2,flag_pred,&
		q2,q_face2,qb_df,qprime2,qprime_face, ope2_ave,uvb_df_ave)

	do k = 1,nlayers
		qp_df_out(1,:,k) = (alpha_mlswe(k)/gravity)*q_df(1,:,k)
		qp_df_out(2,:,k) = q_df(2,:,k) / q_df(1,:,k)
		qp_df_out(3,:,k) = q_df(3,:,k) / q_df(1,:,k)
	end do

	mslwe_elevation(:,nlayers+1) = zbot_df

	do k = nlayers,1,-1
		mslwe_elevation(:,k) = mslwe_elevation(:,k+1) + qp_df_out(1,:,k)
	end do

	qp_df_out(4,:,1) = qb_df(3,:)
	qp_df_out(4,:,2) = qb_df(3,:)

	qp_df_out(5,:,1) = one_plus_eta(:)-1.0  !mslwe_elevation(:,1)
	!qp_df_out(5,:,1) = mslwe_elevation(:,1)
	qp_df_out(5,:,2:nlayers) = mslwe_elevation(:,2:nlayers)

end subroutine ti_rk35_bcl