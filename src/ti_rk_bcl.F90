! ===========================================================================================================================
! This module contains the routines for the predictor-corrector (for baroclinic) and the RK35 time integration (for barotropic) methods
!   Author: Yao Gahounzo 
!   Computing PhD 
!   Boise State University
!   Date: October 27, 2023
! ==========================================================================================================================

subroutine ti_rk_bcl(q_df, qb_df, qprime_df)

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
	use mod_basis, only: nq, ngl
	use mod_rk_mlswe, only: ti_barotropic_ssprk_mlswe

	use mod_variables, only: one_plus_eta_df, dpprime_visc, dpprime_visc_q
	use mod_barotropic_terms, only: btp_bcl_coeffs_qdf
    use mod_layer_terms, only: interpolate_qprime, extract_qprime_df_face

	implicit none

	real, dimension(4,npoin), intent(inout) :: qb_df
	real, dimension(3,npoin,nlayers), intent(inout) :: q_df
	real, dimension(3,npoin,nlayers), intent(inout) :: qprime_df

	real, dimension(3,npoin_q,nlayers) :: qprime2, qprime
	real, dimension(3,2,nq,nface,nlayers) :: qprime_face2, qprime_face
	real, dimension(3,2,ngl,nface,nlayers) :: qprime_df_face2, qprime_df_face
	real, dimension(npoin,nlayers) :: dpprime_df2
	real, dimension(4,npoin) :: qbp_df
	real, dimension(3,npoin,nlayers) :: qprime_df2, q_df2

	integer :: k

	! ==================== Prediction step =================================

	call interpolate_qprime(qprime,qprime_face,qprime_df_face,qprime_df)
	call bcl_create_communicator(qprime_face,3,nlayers,nq)
        call extract_qprime_df_face(qprime_df_face,qprime_df)
	call bcl_create_communicator(qprime_df_face,3,nlayers,ngl)

	qbp_df = qb_df
	dpprime_visc(:,:) = qprime_df(1,:,:)

	!call btp_bcl_coeffs_qdf(qprime,qprime_face,qprime_df_face, qprime_df)
	call btp_bcl_coeffs_qdf(qprime_df_face, qprime_df)
	call ti_barotropic_ssprk_mlswe(qbp_df, qprime_df)

	qprime2 = qprime
	qprime_face2 = qprime_face
	q_df2 = q_df
	qprime_df2 = qprime_df
        qprime_df_face2 = qprime_df_face

	call momentum_mass(qprime2,q_df2,qprime_face2,qprime_df_face2,qprime_df2,qbp_df, qprime_face)

	! ==================== Correction step =================================

	! Communication of qprime_face2 values within the inter-processor boundary
	call bcl_create_communicator(qprime_face2,3,nlayers,nq)

	qprime2 = 0.5*(qprime2 + qprime)
	qprime_face2 = 0.5*(qprime_face2 + qprime_face)
	qprime_df2 = 0.5*(qprime_df2 + qprime_df)
        !call extract_qprime_df_face(qprime_df_face2,qprime_df2)
	call bcl_create_communicator(qprime_df_face2,3,nlayers,ngl)
        qprime_df_face2 = 0.5*(qprime_df_face + qprime_df_face2)

	dpprime_visc(:,:) = qprime_df2(1,:,:)

	!call btp_bcl_coeffs_qdf(qprime2,qprime_face2,qprime_df_face2, qprime_df2)
	call btp_bcl_coeffs_qdf(qprime_df_face2, qprime_df2)
	call ti_barotropic_ssprk_mlswe(qb_df,qprime_df2)

	call thickness(qprime_df2, q_df, qprime_face2, dpprime_df2, qb_df, qprime_df_face2)

	! Communication of qprime_face_avg values within the processor boundary
	call bcl_create_communicator(qprime_face2(1,:,:,:,:),1,nlayers,nq)
        call bcl_create_communicator(qprime_df_face2(1,:,:,:,:),1,nlayers,ngl)

        qprime_df2(1,:,:) = 0.5*(qprime_df(1,:,:) + dpprime_df2(:,:))
        qprime2(1,:,:) = 0.5*(qprime(1,:,:) + qprime2(1,:,:))
        qprime_face2(1,:,:,:,:) = 0.5*(qprime_face(1,:,:,:,:) + qprime_face2(1,:,:,:,:))
        qprime_df_face2(1,:,:,:,:) = 0.5*(qprime_df_face(1,:,:,:,:) + qprime_df_face2(1,:,:,:,:))
        
        call momentum(qprime2,q_df,qprime_face2,qprime_df2,qb_df,qprime_face, qprime_df_face2)

        qprime_df(1,:,:) = dpprime_df2
        qprime_df(2:3,:,:) = qprime_df2(2:3,:,:) 

end subroutine ti_rk_bcl
