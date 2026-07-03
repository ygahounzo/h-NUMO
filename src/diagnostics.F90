subroutine diagnostics(G, inp, gg, par, init, q, q_df, qb, itime, idone)

    use mod_grid,          only: grid
    use mod_input,         only: input
    use mod_global_grid,   only: grid_global
    use mod_parallel,      only: parallel_CS
    use mod_initial,       only: initial
    use mod_mpi_utilities, only: irank, irank0
    use mod_constants,     only: gravity

    implicit none

    type(grid),        intent(in) :: G
    type(input),       intent(in) :: inp
    type(grid_global), intent(in) :: gg
    type(parallel_CS), intent(in) :: par
    type(initial),     intent(in) :: init

    real, intent(in)    :: q_df(inp%nvar_bcl, G%npoin, inp%nlayers)
    real, intent(in)    :: qb(inp%nvar_btp, G%npoin)
    integer, intent(in) :: itime, idone
    real, intent(out)   :: q(5, G%npoin, inp%nlayers)

    character*5 :: tempchar
    character*4 :: num
    integer  :: i, j, k, iloop
    logical  :: has_w
    real     :: q_gg(5, gg%npoin_g, inp%nlayers), ql(5, G%npoin), zbot_g(gg%npoin_g)
    real :: coord_dg_gathered(3, gg%npoin_g), qb_g(4, gg%npoin_g), q_g(5, gg%npoin_g)
    real, dimension(G%npoin, inp%nlayers+1) :: mslwe_elevation

    has_w = (inp%nvar_bcl == 4) ! w (vertical momentum) is only carried on sphere_hex

    do k = 1, inp%nlayers
        q(1,:,k) = (init%alpha_mlswe(k)/gravity)*q_df(1,:,k)
        q(2,:,k) = q_df(2,:,k) / q_df(1,:,k)
        q(3,:,k) = q_df(3,:,k) / q_df(1,:,k)
        if (has_w) then
            q(4,:,k) = q_df(4,:,k) / q_df(1,:,k)
        else
            q(4,:,k) = q_df(1,:,k)
        end if
    end do

    mslwe_elevation = 0.0
    mslwe_elevation(:,inp%nlayers+1) = init%zbot_df

    do k = inp%nlayers, 1, -1
        mslwe_elevation(:,k) = mslwe_elevation(:,k+1) + q(1,:,k)
    end do

    q(5,:,1) = mslwe_elevation(:,1)
    q(5,:,2:inp%nlayers) = mslwe_elevation(:,2:inp%nlayers)

    ! Gather Data onto Head node
    do k = 1, inp%nlayers
        ql = q(:,:,k)
        call gather_data(G, inp, gg, par, q_g, ql, 5)
        q_gg(:,:,k) = q_g(:,:)
    enddo

    call gather_data(G, inp, gg, par, qb_g,             qb(1:4,:),    4)
    call gather_data(G, inp, gg, par, coord_dg_gathered, G%coord,     3)
    call gather_data(G, inp, gg, par, zbot_g,            init%zbot_df, 1)

    if (irank == irank0 .and. idone == 0) then

        tempchar = 'mlswe'
        write(num,'(i4)') itime

        if(itime == 0) then
            num = '0000'
        else
            iloop=3 - int(log10(real(itime)))
            do j=1,iloop
                num(j:j) = '0'
            end do
        end if
        open(unit=10, file = tempchar//num)

        write(10,'(i4)') inp%nlayers
        write(10,'(i10)') gg%npoin_g
        write(10,'(d23.16)') inp%dt
        write(10,'(d23.16)') inp%dt_btp
        write(10,'(d23.16)') coord_dg_gathered(1:2,:)
        write(10,'(d23.16)') qb_g(1,:)
        write(10,'(d23.16)') qb_g(3,:)
        write(10,'(d23.16)') qb_g(4,:)
        write(10,'(d23.16)') q_gg(1,:,:)
        write(10,'(d23.16)') q_gg(2,:,:)
        write(10,'(d23.16)') q_gg(3,:,:)
        write(10,'(d23.16)') q_gg(5,:,:)
        write(10,'(d23.16)') zbot_g(:)

        close(10)
    end if

    return
end subroutine diagnostics
