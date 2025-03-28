!---------------------------------------------------------------------!
!> @brief This subroutine read output files to restart simulation
!> 
!> @author Yao Gahounzo on 01/2024
!---------------------------------------------------------------------!

subroutine read_mlswe(q_df,qb_df,fname)
  
    use mod_basis, only: nglx, ngly, nglz, CELL_CHILDREN, nsubcells

    use mod_grid, only: coord, npoin, nelem, intma
  
    use mod_initial, only: nvar

    use mod_input, only: format_vtk, nvtk_files,  fname_root, vtk_cell_type, nlayers

    use mod_parallel, only: npoin_l, npoin_l_max, nelem_l, nelem_l_max

    use mpi
  
    implicit none
  
    !global
    real :: q_df(3,npoin,nlayers), qb_df(3,npoin)
    character fname*100,fname_proc*200,proc_num*10
    character grid_file*100

    !local
    real, dimension(:,:), allocatable :: qb_grp
    real, dimension(:,:, :), allocatable :: q_grp

    real    :: x,y,z,r,lat,lon
    integer :: ncells, nchildren, nsize, nglm1, nglm13, npoly
    integer :: ivar, ifactor, ip, ic, ie, i, j, k, e, ndof, im, mp, me
    integer :: irank, ierr, nproc
    character :: form*10

    integer :: out_group, ngroups
    integer :: igroup, ig_rank, buffer
    integer :: npoin_grp, nelem_grp, ncells_grp, nproc_grp

    integer, dimension(:), allocatable  :: npoin_l_grp, nelem_l_grp
    integer :: npoin_l_max_grp, nelem_l_max_grp

    ngroups = 1

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    !create MPI groups
    call set_groups(out_group,ngroups,igroup,ig_rank)
    call mpi_comm_size(out_group,nproc_grp,ierr)

    !create cells
    ncells=nelem*nsubcells

    if(ig_rank==0) allocate(npoin_l_grp(nproc_grp),nelem_l_grp(nproc_grp))

    !count points and elements in the group
    call mpi_gather(npoin,1,mpi_integer,npoin_l_grp,1,mpi_integer,0,out_group,ierr)
    call mpi_gather(nelem,1,mpi_integer,nelem_l_grp,1,mpi_integer,0,out_group,ierr)

    !allocate memory on each group leader rank
    if(ig_rank==0) then       

        npoin_grp = sum(npoin_l_grp)
        nelem_grp = sum(nelem_l_grp)

        npoin_l_max_grp = maxval(npoin_l_grp)
        nelem_l_max_grp = maxval(nelem_l_grp)

        ncells_grp=nelem_grp*nsubcells
        allocate(q_grp(3,npoin_grp, nlayers))
        allocate(qb_grp(3,npoin_grp))

        call load_data_mlswe(fname, qb_grp, q_grp, npoin_grp, nlayers)

    end if

    ! !scatter data within groups
    call scatter_data_to(qb_grp, qb_df, 3, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)

    do k = 1, nlayers
        call scatter_data_to(q_grp(:,:,k), q_df(:,:,k), 3, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    end do
    
    if(ig_rank==0) then
       deallocate(q_grp, qb_grp) !not the best idea, but will work for now
       deallocate(npoin_l_grp,nelem_l_grp)
    end if

end subroutine read_mlswe

subroutine restart_mlswe(q_df,qb_df,q,qb,qprime,qprime_df,q_face,qprime_face,qb_face, qp_df_out, q_df_read, qb_df_read)

        use mod_grid, only : npoin_q, npoin, intma_dg_quad, intma,nface,face
        use mod_basis, only: npts
        use mod_input, only: nlayers
        use mod_initial, only: psih, indexq, pbprime_df, alpha_mlswe, pbprime, zbot_df
        use mod_face, only: imapl_q, imapr_q
        use mod_constants, only: gravity
        use mod_layer_terms, only: evaluate_mom, evaluate_mom_face, evaluate_dp, evaluate_dp_face

        real, dimension(3,npoin), intent(in) :: qb_df_read
        real, dimension(3,npoin,nlayers), intent(in) :: q_df_read

        real, dimension(4,npoin), intent(out) :: qb_df
        real, dimension(3,npoin,nlayers), intent(out) :: q_df
        real, dimension(4,npoin_q), intent(out) :: qb
        real, dimension(3,npoin_q,nlayers), intent(out) :: q

        real, dimension(3,npoin_q,nlayers), intent(out) :: qprime
        real, dimension(3,npoin,nlayers), intent(out) :: qprime_df
        real, dimension(3,2,nq,nface,nlayers), intent(out) :: q_face
        real, dimension(3,2,nq,nface,nlayers), intent(out) :: qprime_face
        real, dimension(4,2,nq,nface), intent(out) :: qb_face
        real, dimension(5,npoin,nlayers) :: qp_df_out

        real, dimension(npoin,nlayers+1) :: mslwe_elevation

        integer :: I, Iq, ip, iface, el, er, iquad, ir, jr, kr, il, jl, kl, n,k, ilocl, ilocr
        real :: hi
        real, dimension(npoin) :: one_plus_eta_temp

        q = 0.0
        qb = 0.0
        qprime = 0.0
        q_face = 0.0
        qprime_face = 0.0
        qb_face = 0.0

        ! Barotropic variables at the dofs (nodal points)

        qb_df(1,:) = qb_df_read(1,:)
        qb_df(3:4,:) = qb_df_read(2:3,:)
        qb_df(2,:) = qb_df(1,:) - pbprime_df(:)

        ! Interpolate to the quadrature points

        call btp_evaluate_mom_dp(qb,qb_df)
        call btp_evaluate_mom_dp_face(qb_face, qb)

        ! Baroclinic variables at the dofs (nodal points)

        !q_df(:,:,:) = q_df_read(:,:,:)

        do k = 1,nlayers

            q_df(1,:,k) = (gravity/alpha_mlswe(k))*q_df_read(1,:,k)
            q_df(2,:,k) = q_df_read(2,:,k)*q_df(1,:,k)
            q_df(3,:,k) = q_df_read(3,:,k)*q_df(1,:,k)

        end do

        ! Interpolate to the quadrature points

        call evaluate_dp(q,qprime,q_df, pbprime)
        call evaluate_dp_face(q_face, qprime_face,q, qprime)

        call evaluate_mom(q,q_df)
        call evaluate_mom_face(q_face, q)

        ! Prime variables at the dofs (nodal points) and quadrature points
        one_plus_eta_temp(:) = sum(q_df(1,:,:),dim=2) / pbprime_df(:)

        do k = 1,nlayers

            qprime_df(1,:,k) = q_df(1,:,k) / one_plus_eta_temp(:)

            qprime(2,:,k) = q(2,:,k)/q(1,:,k) - qb(3,:)/qb(1,:)
            qprime(3,:,k) = q(3,:,k)/q(1,:,k) - qb(4,:)/qb(1,:)
            qprime_df(2,:,k) = q_df(2,:,k)/q_df(1,:,k) - qb_df(3,:)/qb_df(1,:)
            qprime_df(3,:,k) = q_df(3,:,k)/q_df(1,:,k) - qb_df(4,:)/qb_df(1,:)

            qprime_face(2,:,:,:,k) = q_face(2,:,:,:,k)/q_face(1,:,:,:,k) - qb_face(3,:,:,:)/qb_face(1,:,:,:)
            qprime_face(3,:,:,:,k) = q_face(3,:,:,:,k)/q_face(1,:,:,:,k) - qb_face(4,:,:,:)/qb_face(1,:,:,:)
        end do

        ! Prepare output variables

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

        qp_df_out(5,:,1) = mslwe_elevation(:,1)
        qp_df_out(5,:,2) = mslwe_elevation(:,2)

end subroutine restart_mlswe

subroutine load_data_mlswe(name_data_file, qb_df_g, q_df_g, npoin_g, nlayers)

    character(len=*), intent(in) :: name_data_file

    real, dimension(3, npoin_g), intent(out) :: qb_df_g
    real, dimension(3, npoin_g, nlayers), intent(out) :: q_df_g

    real, dimension(2, npoin_g) :: coord_g
    real, dimension(npoin_g, nlayers+1) :: z_g

    integer :: count, nk, npoin_gg, i, j
    real :: dt, dt_btp

    open(unit=10, file=name_data_file, status='old', action='read')

    read(10, *) nk
    read(10, *) npoin_gg
    read(10, *) dt
    read(10, *) dt_btp

    read(10, *) (coord_g(:, i), i = 1, npoin_gg)
    
    do j = 1, 3
        read(10, *) (qb_df_g(j, i), i = 1, npoin_gg)
    end do

    read(10, *) (q_df_g(1, :, j), j = 1, nk)
    read(10, *) (q_df_g(2, :, j), j = 1, nk)
    read(10, *) (q_df_g(3, :, j), j = 1, nk)

    read(10, *) (z_g(:, j), j = 1, nk + 1)

    close(10)
end subroutine load_data_mlswe

!---------------------------------------------------------------------!
!> @brief This subroutine creates communicators for ngroups of 
!> processors (group_comm) and returns group number (igroup) and
!> rank within the group (ig_rank) 
!> 
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine set_groups(group_comm, ngroups, igroup, ig_rank)

    use mpi

    implicit none

    integer :: group_comm, ngroups

    integer :: irank, ierr, i
    integer :: igroup, ig_rank, group_size, nproc

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)
  
    !find desired group size
    group_size = ceiling(nproc/real(ngroups))

    if(group_size*(ngroups-1)==nproc) ngroups=ngroups-1

    !find to which group the process belongs and what is it's rank in the group
    igroup = int(irank/group_size)
    ig_rank = mod(irank,group_size)

    call mpi_comm_split (mpi_comm_world, igroup, irank, group_comm, ierr)

end subroutine set_groups
