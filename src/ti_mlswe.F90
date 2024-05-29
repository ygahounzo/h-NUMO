! ===========================================================================================================================
! This module contains the routines for the predictor-corrector time method for MLSWE
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: April 27, 2023
! ==========================================================================================================================

subroutine ti_mlswe(q, q_df, q_face, qb, qb_face, qb_df, qprime, qprime_face,dpprime_df,qprime_df,qp_df_out)

	use mod_splitting, only: ti_barotropic, thickness, momentum
	use mod_input, only: nlayers, dt, dt_btp, dpprime_visc_min, ti_method_btp
	use mod_grid, only: npoin, npoin_q, nface
	use mod_constants, only: gravity
	use mod_initial, only: alpha_mlswe, zbot_df
	use mod_basis, only: nq
	use mod_variables, only: one_plus_eta_df

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

	real, dimension(npoin,nlayers+1) :: mslwe_elevation
	real, dimension(5,npoin,nlayers):: qp_df

	real, dimension(3,npoin_q,nlayers) :: qprime_avg
	real, dimension(3,npoin_q,nlayers) :: qprime1, qprime2
	real, dimension(3,npoin_q,nlayers) :: qprime_corr
	real, dimension(3,2,nq,nface,nlayers) :: qprime_face_corr
	real, dimension(3,2,nq,nface,nlayers) :: qprime_face1, qprime_face2
	real, dimension(3,2,nq,nface,nlayers) :: qprime_face_avg
	real, dimension(npoin,nlayers) :: dpprime_df1, dpprime_df2
	real, dimension(3,npoin_q,nlayers) :: q2
	real, dimension(3,2,nq,nface,nlayers) :: q_face2
	real, dimension(4,npoin) :: qbp_df
	real, dimension(3,npoin,nlayers) :: q_df1
	real, dimension(3,npoin,nlayers) :: qprime_df_avg, qprime_df2, qprime_df_corr
	real, dimension(4,npoin_q) :: qbp
	real, dimension(4,2,nq,nface) :: qbp_face
	real, dimension(npoin_q,nlayers) :: dpprime 

	integer :: k, Iq, iquad, iface, ilr, flag_pred

	! Output variables
	real, dimension(2,nq, nface, nlayers) :: u_edge, v_edge
	real, dimension(2,nq,nface,nlayers) :: dprime_face_corr

	mslwe_elevation = 0.0

	! Prediction step

	flag_pred = 1

	call create_communicator_quad_layer(qprime_face,3,nlayers)
	call create_communicator_quad_layer(q_face,3,nlayers)

	qbp = qb
	qbp_face = qb_face
	qbp_df = qb_df

	call ti_barotropic(qbp,qbp_face,qbp_df,qprime,qprime_face, qprime_df, flag_pred)

	q2 = q
	q_face2 = q_face
	qprime1 = qprime
	q_df1 = q_df
	qprime_face1 = qprime_face

	call thickness(q2,qprime1, q_df1, q_face2, qprime_face1, u_edge, v_edge, dpprime_df1, flag_pred, qbp_df)

	qprime2 = qprime
	qprime_face2 = qprime_face
	qprime_df2 = qprime_df

	call momentum(q2,qprime2,q_df1,q_face2,qprime_face2,qbp,qbp_face,u_edge,v_edge,&
		qprime_df2, qprime_face,flag_pred,q,q_face,qbp_df,qprime, qprime_face1)

	qprime2(1,:,:) = qprime1(1,:,:)
	qprime_face2(1,:,:,:,:) = qprime_face1(1,:,:,:,:)
	qprime_df2(1,:,:) = dpprime_df1(:,:)

	! Correction step

	! Communication of qprime_face2 values within the inter-processor boundary
	call create_communicator_quad_layer(qprime_face2,3,nlayers)

	qprime_avg = 0.5*(qprime2 + qprime)
	qprime_face_avg = 0.5*(qprime_face2 + qprime_face)

	! call create_communicator_quad_layer(qprime_face_avg,3,nlayers)

	qprime_df_avg = 0.5*(qprime_df2 + qprime_df)
	
	flag_pred = 0

	call ti_barotropic(qb,qb_face,qb_df, qprime_avg,qprime_face_avg, qprime_df_avg, flag_pred)


	q2 = q
	q_face2 = q_face

	qprime2 = qprime_avg
	qprime_face2 = qprime_face_avg

	qprime_face_avg(1,:,:,:,:) = qprime_face(1,:,:,:,:) ! No correction for thickness

	call thickness(q,qprime_avg, q_df, q_face, qprime_face_avg, u_edge, v_edge, dpprime_df2, flag_pred, qb_df)

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

	call momentum(q,qprime_corr,q_df,q_face,qprime_face_corr,qb,qb_face,u_edge,v_edge,&
		qprime_df_corr,qprime_face2,flag_pred,q2,q_face2,qb_df,qprime2,qprime_face)

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

	qp_df_out(5,:,1) = one_plus_eta_df(:)-1.0  !mslwe_elevation(:,1)
	qp_df_out(5,:,2:nlayers) = mslwe_elevation(:,2:nlayers)

end subroutine ti_mlswe