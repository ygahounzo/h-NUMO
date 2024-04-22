subroutine ti_rk35_mlswe(q, q_df, q_face, qb, qb_face, qb_df, qprime, qprime_face,dpprime_df,qprime_df,qp_df_out)

	use mod_splitting, only: thickness, momentum, momentum_mass
	use mod_input, only: nlayers, ti_method_btp
	use mod_grid, only: npoin, npoin_q, nface
	use mod_constants, only: gravity
	use mod_initial, only: alpha_mlswe, zbot_df
	use mod_basis, only: nq
	use mod_rk_mlswe, only: ti_barotropic_rk_mlswe

	implicit none

	real, dimension(4,npoin_q), intent(inout) :: qb
	real, dimension(4,2,nq,nface), intent(inout) :: qb_face
	real, dimension(4,npoin), intent(inout) :: qb_df
	real, dimension(3,npoin_q,nlayers), intent(inout) :: qprime
	real, dimension(3,2,nq,nface,nlayers), intent(inout) :: qprime_face
	real, dimension(3,npoin_q,nlayers), intent(inout) :: q
	real, dimension(3,2,nq,nface,nlayers), intent(inout) :: q_face
	real, dimension(3,npoin,nlayers), intent(inout) :: q_df
	real, dimension(npoin,nlayers), intent(inout) :: dpprime_df
	real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df

	real, dimension(5,npoin,nlayers), intent(out) :: qp_df_out
	real, dimension(npoin) :: one_plus_eta

	real, dimension(nq,nface) :: one_plus_eta_edge_2_ave
	real, dimension(npoin_q) :: ope_ave, ope2_ave, H_ave, Qu_ave, Qv_ave, Quv_ave
	real, dimension(2,npoin_q) :: btp_mass_flux_ave, uvb_ave 
	real, dimension(npoin) :: ope_ave_df
	real, dimension(2,2,nq,nface) :: uvb_face_ave
	real, dimension(2,nq,nface) :: btp_mass_flux_face_ave, Qu_face_ave, Qv_face_ave, Quv_face_ave, ope_face_ave
	real, dimension(nq,nface) :: H_face_ave
	real, dimension(2,npoin_q) :: tau_wind_ave, tau_bot_ave
	real, dimension(npoin,nlayers+1) :: mslwe_elevation

	real, dimension(3,npoin_q,nlayers) :: qprime_avg, qprime2, qprime_corr, q2
	real, dimension(3,2,nq,nface,nlayers) :: qprime_face_corr, q_face2, qprime_face2, qprime_face_avg
	real, dimension(npoin,nlayers) :: dpprime_df2
	real, dimension(4,npoin) :: qbp_df
	real, dimension(3,npoin,nlayers) :: qprime_df_avg, qprime_df2, qprime_df_corr, q_df2
	real, dimension(4,npoin_q) :: qbp
	real, dimension(4,2,nq,nface) :: qbp_face
	real, dimension(nq,nface,nlayers) :: disp

	integer :: k, Iq, iquad, iface, ilr, flag_pred


	! Output variables
	real, dimension(2,nq, nface, nlayers) :: u_edge, v_edge
	real, dimension(2,nq,nface,nlayers) :: dprime_face_corr

	mslwe_elevation = 0.0

	! Prediction step

	flag_pred = 1

	!call create_communicator_quad_layer(qprime_face,3,nlayers)
	!call create_communicator_quad_layer(q_face,3,nlayers)

	call create_communicator_quad_layer_all(qprime_face,q_face,3,nlayers)

	qbp = qb
	qbp_face = qb_face
	qbp_df = qb_df


	call ti_barotropic_rk_mlswe(one_plus_eta,one_plus_eta_edge_2_ave,uvb_ave, ope_ave, ope2_ave, &
		H_ave, Qu_ave, Qv_ave, Quv_ave, btp_mass_flux_ave, ope_ave_df, uvb_face_ave, &
		ope_face_ave, btp_mass_flux_face_ave, H_face_ave, &
		Qu_face_ave, Qv_face_ave, Quv_face_ave, tau_wind_ave, tau_bot_ave, qbp,qbp_face,&
		qbp_df,qprime,qprime_face, qprime_df, flag_pred)


	qprime2 = qprime
	qprime_face2 = qprime_face
	q_df2 = q_df
	qprime_df2 = qprime_df
	q2 = q
	q_face2 = q_face


	call momentum_mass(q2,qprime2,q_df2,q_face2,qprime_face2,qbp,qbp_face,ope_ave,one_plus_eta_edge_2_ave,uvb_ave,Qu_ave,&
        Qv_ave,Quv_ave,H_ave,uvb_face_ave,ope_face_ave,Qu_face_ave,Qv_face_ave,Quv_face_ave,H_face_ave,&
        qprime_df2,ope_ave_df,tau_bot_ave,tau_wind_ave,qprime_face2,flag_pred,&
        q,q_face,qbp_df,qprime, qprime_face, ope2_ave, btp_mass_flux_ave, btp_mass_flux_face_ave)


	! Correction step

	! Communication of qprime_face2 values within the inter-processor boundary
	call create_communicator_quad_layer(qprime_face2,3,nlayers)

	qprime_avg = 0.5*(qprime2 + qprime)
	qprime_face_avg = 0.5*(qprime_face2 + qprime_face)

	qprime_df_avg = 0.5*(qprime_df2 + qprime_df)
	
	flag_pred = 0

	call ti_barotropic_rk_mlswe(one_plus_eta,one_plus_eta_edge_2_ave,uvb_ave, ope_ave, ope2_ave, &
		H_ave, Qu_ave, Qv_ave, Quv_ave, btp_mass_flux_ave, ope_ave_df, uvb_face_ave, &
		ope_face_ave, btp_mass_flux_face_ave, H_face_ave, &
		Qu_face_ave, Qv_face_ave, Quv_face_ave, tau_wind_ave, tau_bot_ave, qb,qb_face,&
		qb_df, qprime_avg,qprime_face_avg, qprime_df_avg, flag_pred)


	q2 = q
	q_face2 = q_face

	qprime2 = qprime_avg
	qprime_face2 = qprime_face_avg

	!qprime_face_avg(1,:,:,:,:) = qprime_face(1,:,:,:,:) ! No correction for thickness

	call thickness(q,qprime_avg, q_df, q_face, qprime_face_avg, u_edge, v_edge, uvb_ave, btp_mass_flux_ave, ope_ave, uvb_face_ave, ope_face_ave, &
		btp_mass_flux_face_ave, dpprime_df2, flag_pred, qb_df, disp)

	dprime_face_corr = qprime_face_avg(1,:,:,:,:)

	! Communication of qprime_face_avg values within the processor boundary
	call create_communicator_quad_layer(qprime_face_avg(1,:,:,:,:),1,nlayers)

	qprime_df_avg(1,:,:) = 0.5*(qprime_df(1,:,:) + dpprime_df2(:,:))

	qprime_corr(1,:,:) = 0.5*(qprime(1,:,:) + qprime_avg(1,:,:))
	qprime_face_corr(1,:,:,:,:) = 0.5*(qprime_face(1,:,:,:,:) + qprime_face_avg(1,:,:,:,:))
	qprime_corr(2:3,:,:) = qprime_avg(2:3,:,:)
	qprime_face_corr(2:3,:,:,:,:) = qprime_face_avg(2:3,:,:,:,:)
	! qprime_df = qprime_df_avg
	qprime_df_corr(1,:,:) = 0.5*(qprime_df(1,:,:) + dpprime_df2(:,:))
	qprime_df_corr(2:3,:,:) = qprime_df_avg(2:3,:,:)

	call momentum(q,qprime_corr,q_df,q_face,qprime_face_corr,qb,qb_face,ope_ave,one_plus_eta_edge_2_ave,uvb_ave,u_edge,v_edge,Qu_ave,&
		Qv_ave,Quv_ave,H_ave,uvb_face_ave,ope_face_ave,Qu_face_ave,Qv_face_ave,Quv_face_ave,H_face_ave,&
		qprime_df_corr,ope_ave_df,tau_bot_ave,tau_wind_ave, qprime_face2,flag_pred,&
		q2,q_face2,qb_df,qprime2,qprime_face, ope2_ave, disp)

	qprime(1,:,:) = qprime_avg(1,:,:)
	qprime(2:3,:,:) = qprime_corr(2:3,:,:)
	qprime_face(1,:,:,:,:) = dprime_face_corr
	qprime_face(2:3,:,:,:,:) = qprime_face_corr(2:3,:,:,:,:)
	
	qprime_df(1,:,:) = dpprime_df2
	qprime_df(2:3,:,:) = qprime_df_corr(2:3,:,:)

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

	qp_df_out(5,:,1) = one_plus_eta(:)-1.0!mslwe_elevation(:,1)
	qp_df_out(5,:,2:nlayers) = mslwe_elevation(:,2:nlayers)

end subroutine ti_rk35_mlswe