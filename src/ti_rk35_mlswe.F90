! ===========================================================================================================================
! This module contains the routines for the predictor-corrector (for baroclinic) and the RK35 time integration (for barotropic) methods
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: October 27, 2023
! ==========================================================================================================================

subroutine ti_rk35_mlswe(q, q_df, q_face, qb, qb_face, qb_df, qprime, qprime_face,dpprime_df,qprime_df,qp_df_out)

	! q: layer variable dp, u*dp, v*dp at quad points and their face values: q_face
	! q_df : layer variable dp, u*dp, v*dp at nodal (dof) point
	! qprime: value dp', u' and v' at quad points and their face values: qprime_face
	! qb : barotopic variable pb, pb_pert = pb'*eta, ub*pb, vb*pb at quad points and their face values: qb_face
	! qb_df: : barotopic variable pb, pb_pert = pb'*eta, ub*pb, vb*pb at nodal points 
	! qprime_df: value dp', u' and v' at nodal points
	! qp_df_out: output variable, thickness h_k, velocity u_k,v_k, free surface ssh

	use mod_splitting, only: thickness, momentum, momentum_mass
	use mod_input, only: nlayers, dpprime_visc_min
	use mod_grid, only: npoin, npoin_q, nface
	use mod_constants, only: gravity
	use mod_initial, only: alpha_mlswe, zbot_df
	use mod_basis, only: nq
	use mod_rk_mlswe, only: ti_barotropic_ssprk_mlswe
	use mod_barotropic_terms, only: btp_bcl_coeffs_qdf
	use mod_variables, only: one_plus_eta_df, dpprime_visc, dpprime_visc_q

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
	
	real, dimension(npoin,nlayers+1) :: mslwe_elevation

	real, dimension(3,npoin_q,nlayers) :: qprime_avg, qprime2, qprime_corr, q2
	real, dimension(3,2,nq,nface,nlayers) :: qprime_face_corr, q_face2, qprime_face2, qprime_face_avg
	real, dimension(npoin,nlayers) :: dpprime_df2
	real, dimension(4,npoin) :: qbp_df
	real, dimension(3,npoin,nlayers) :: qprime_df_avg, qprime_df2, qprime_df_corr, q_df2
	real, dimension(4,npoin_q) :: qbp
	real, dimension(4,2,nq,nface) :: qbp_face

	integer :: k


	! Output variables
	real, dimension(2,nq,nface,nlayers) :: dprime_face_corr

	mslwe_elevation = 0.0

	! ==================== Prediction step =================================

	!call create_communicator_quad_layer_all(qprime_face,q_face,3,nlayers)
	call bcl_create_communicator(qprime_face,3,nlayers,nq)

	qbp_df = qb_df

	dpprime_visc(:,:) = qprime_df(1,:,:)
	dpprime_visc_q(:,:) = qprime(1,:,:)

	call btp_bcl_coeffs_qdf(qprime,qprime_face, qprime_df)
	call ti_barotropic_ssprk_mlswe(qbp_df,qprime)

	qprime2 = qprime
	qprime_face2 = qprime_face
	q_df2 = q_df
	qprime_df2 = qprime_df

	call momentum_mass(qprime2,q_df2,qprime_face2,qprime_df2,qbp_df, qprime_face)

	! ==================== Correction step =================================

	! Communication of qprime_face2 values within the inter-processor boundary
	call bcl_create_communicator(qprime_face2,3,nlayers,nq)

	qprime_avg = 0.5*(qprime2 + qprime)
	qprime_face_avg = 0.5*(qprime_face2 + qprime_face)

	qprime_df_avg = 0.5*(qprime_df2 + qprime_df)

	dpprime_visc(:,:) = qprime_df_avg(1,:,:)
	dpprime_visc_q(:,:) = qprime_avg(1,:,:)

	call btp_bcl_coeffs_qdf(qprime_avg,qprime_face_avg, qprime_df_avg)
	call ti_barotropic_ssprk_mlswe(qb_df,qprime_avg)

	!qprime_face_avg(1,:,:,:,:) = qprime_face(1,:,:,:,:) ! No correction for thickness

	call thickness(qprime_avg, q_df, qprime_face_avg, dpprime_df2, qb_df)

	dprime_face_corr = qprime_face_avg(1,:,:,:,:)

	! Communication of qprime_face_avg values within the processor boundary
	call bcl_create_communicator(qprime_face_avg(1,:,:,:,:),1,nlayers,nq)

	qprime_df_avg(1,:,:) = 0.5*(qprime_df(1,:,:) + dpprime_df2(:,:))

	qprime_corr(1,:,:) = 0.5*(qprime(1,:,:) + qprime_avg(1,:,:))
	qprime_face_corr(1,:,:,:,:) = 0.5*(qprime_face(1,:,:,:,:) + qprime_face_avg(1,:,:,:,:))
	qprime_corr(2:3,:,:) = qprime_avg(2:3,:,:)
	qprime_face_corr(2:3,:,:,:,:) = qprime_face_avg(2:3,:,:,:,:)
	qprime_df_corr(1,:,:) = 0.5*(qprime_df(1,:,:) + dpprime_df2(:,:))
	qprime_df_corr(2:3,:,:) = qprime_df_avg(2:3,:,:)

	call momentum(qprime_corr,q_df,qprime_face_corr,qprime_df_corr,qb_df,qprime_face)

	qprime(1,:,:) = qprime_avg(1,:,:)
	qprime(2:3,:,:) = qprime_corr(2:3,:,:)
	qprime_face(1,:,:,:,:) = dprime_face_corr
	qprime_face(2:3,:,:,:,:) = qprime_face_corr(2:3,:,:,:,:)
	
	qprime_df(1,:,:) = dpprime_df2
	qprime_df(2:3,:,:) = qprime_df_corr(2:3,:,:)

	do k = 1,nlayers
		qp_df_out(1,:,k) = q_df(1,:,k) !(alpha_mlswe(k)/gravity)*q_df(1,:,k)
		!qp_df_out(2,:,k) = q_df(2,:,k) / q_df(1,:,k)
		!qp_df_out(3,:,k) = q_df(3,:,k) / q_df(1,:,k)
		qp_df_out(2,:,k) = qprime_df(2,:,k) + qb_df(3,:)/qb_df(1,:)
		qp_df_out(3,:,k) = qprime_df(3,:,k) + qb_df(4,:)/qb_df(1,:)
		qp_df_out(4,:,k) = q_df(1,:,k)
	end do

	mslwe_elevation(:,nlayers+1) = zbot_df

	do k = nlayers,1,-1
		mslwe_elevation(:,k) = mslwe_elevation(:,k+1) + qp_df_out(1,:,k)
	end do

	!qp_df_out(5,:,1) = one_plus_eta_df(:)-1.0!mslwe_elevation(:,1)
	qp_df_out(5,:,1) = mslwe_elevation(:,1)
	qp_df_out(5,:,2:nlayers) = mslwe_elevation(:,2:nlayers)

end subroutine ti_rk35_mlswe
