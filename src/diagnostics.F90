subroutine diagnostics(q,q_df,qb,itime,idone)

    use mod_input, only: nlayers, dt, dt_btp
    use mod_global_grid, only: coord_g, npoin_g
    use mod_grid, only: npoin, intma, coord, nelem
    use mod_mpi_utilities, only: irank, irank0
    use mod_initial, only: alpha_mlswe, zbot_df, one_over_pbprime_df, z_interface
    use mod_constants, only: gravity

    implicit none

    real, intent(in)    :: q_df(3,npoin,nlayers)
    real, intent(in)    :: qb(4,npoin)
    integer, intent(in) :: itime,idone
    real, intent(out)   :: q(5,npoin,nlayers)

    character*5 :: tempchar
    character*4 :: num
    integer  :: AllocateStatus,i,j,k,iloop, I1
    real     :: q_gg(5,npoin_g,nlayers), ql(5,npoin), zbot_g(npoin_g)
    real :: coord_dg_gathered(3,npoin_g), qb_g(4,npoin_g), q_g(5,npoin_g)
    real, dimension(npoin,nlayers+1) :: mslwe_elevation

    do k = 1,nlayers
        q(1,:,k) = (alpha_mlswe(k)/gravity)*q_df(1,:,k)
        q(2,:,k) = q_df(2,:,k) / q_df(1,:,k)
        q(3,:,k) = q_df(3,:,k) / q_df(1,:,k)
        q(4,:,k) = q_df(1,:,k)
    end do

    mslwe_elevation = 0.0
    mslwe_elevation(:,nlayers+1) = zbot_df

    do k = nlayers,1,-1
        mslwe_elevation(:,k) = mslwe_elevation(:,k+1) + q(1,:,k)
    end do

    !do k = 1,nlayers+1
    !    mslwe_elevation(:,k) = z_interface(:,k)
    !enddo

    !q(5,:,1) = qb(2,:)*one_over_pbprime_df(:) 
    !q(5,:,1) = qb(1,:)*one_over_pbprime_df(:) - 1.0
    q(5,:,1) = mslwe_elevation(:,1)
    q(5,:,2:nlayers) = mslwe_elevation(:,2:nlayers)

    ! Gather Data onto Head node
    do k = 1, nlayers
        ql = q(:,:,k)
        call gather_data(q_g,ql,5)
        q_gg(:,:,k) = q_g(:,:)
    enddo

    call gather_data(qb_g,qb,4)
    call gather_data(coord_dg_gathered,coord,3)
    call gather_data(zbot_g,zbot_df,1)
    
    if (irank == irank0 .and. idone == 0) then

        ! Generate the name of the file to which the data will be written,
        ! and open the file.  
        tempchar = 'mlswe'
        write(num,'(i4)')itime

        if(itime == 0) then
            num = '0000'
        else
            iloop=3 - int(log10(real(itime)))
            do j=1,iloop
                num(j:j) = '0'
            end do
        end if
        open(unit=10, file = tempchar//num )

        ! Write model parameters and the degrees of freedom.

        write(10,'(i4)') nlayers
        write(10,'(i10)') npoin_g
        write(10,'(d23.16)') dt
        write(10,'(d23.16)') dt_btp
        write(10,'(d23.16)')  coord_dg_gathered(1:2,:)
        write(10,'(d23.16)')  qb_g(1,:)
        write(10,'(d23.16)')  qb_g(3,:)
        write(10,'(d23.16)')  qb_g(4,:)
        write(10,'(d23.16)')  q_gg(1,:,:)
        write(10,'(d23.16)')  q_gg(2,:,:)
        write(10,'(d23.16)')  q_gg(3,:,:)
        write(10,'(d23.16)')  q_gg(5,:,:)
        write(10,'(d23.16)')  zbot_g(:)

        close(10)
    end if

    return
end subroutine diagnostics
