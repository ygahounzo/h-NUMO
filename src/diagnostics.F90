subroutine diagnostics(q,qb,itime)

    use mod_input, only: nlayers, dt, dt_btp
    use mod_global_grid, only: coord_g, npoin_g
    use mod_grid, only: npoin, intma, coord, nelem
    use mod_mpi_utilities, only: irank, irank0
    use mod_initial, only: zbot_df

    implicit none

    real, intent(in) :: q(5,npoin,nlayers)
    real, intent(in) :: qb(4,npoin)
    integer, intent(in) :: itime

    character*7:: tempchar
    character*4 :: num
    integer AllocateStatus,i,j,k,iloop
    real :: q_gg(5,npoin_g,nlayers), ql(5,npoin), zbot_g(npoin_g), coord_dg_gathered(3,npoin_g), qb_g(4,npoin_g), q_g(5,npoin_g)

    ! Gather Data onto Head node
   
    do k = 1, nlayers
        ql = q(:,:,k)
        call gather_data(q_g,ql,5)
        q_gg(:,:,k) = q_g(:,:)
    enddo
    call gather_data(qb_g,qb,4)

    call gather_data(coord_dg_gathered,coord,3)

    call gather_data(zbot_g,zbot_df,1)

    if (irank == irank0) then

        ! Generate the name of the file to which the data will be written,
        ! and open the file.  
        tempchar = 'outborah'
        ! write(num,'(i7)') itime
        ! do i=1,7
        !     if (num(i:i).eq.' ') num(i:i) = '0'
        ! enddo
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
