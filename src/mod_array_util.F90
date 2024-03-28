!----------------------------------------------------------------------!
!>@brief This module initializes arrays and copies them
!>@author  Sasa Gabersek on 10/2012
!>            Naval Research Laboratory
!>            Monterey, CA 93943-5216
!----------------------------------------------------------------------!
module mod_array_util

    use mod_basis, only: nglx, ngly, nglz

    use mod_grid, only: index2d, intma, intma_table, nelem, npoin

    use mod_initial, only: nvar

    use mod_types, only: r8

    implicit none

    private

    public :: &
        get_dimensions, &
        create_arrays, copy_arrays, delete_arrays, &
        qreg2, qreg3, nh_i, nz_i , rhs_continuous

    !> module variables and parameters
    real (kind=r8), dimension(:,:,:), allocatable :: qreg3
    real (kind=r8), dimension(:,:),   allocatable :: qreg2
    real (kind=r8), dimension(:,:),   allocatable :: rhs_continuous
    integer :: nh_i, nz_i
       
    logical :: is_mod_initialized = .false.
       
contains

    !----------------------------------------------------------------------!
    subroutine get_dimensions()

        implicit none

        save
   
        !> get the dimensions on this particular tile
        !> the original code, with 3d indexing
        !>    nx_i=maxval(index(1,:))-minval(index(1,:)) + 1
        !>    ny_i=maxval(index(2,:))-minval(index(2,:)) + 1
        !>    nz_i=maxval(index(3,:))-minval(index(3,:)) + 1

        !> the modified code, with 2d indexing (1 for horizontal direction, 1 for vertical)
        nh_i=maxval(index2d(1,:))-minval(index2d(1,:)) + 1
        nz_i=maxval(index2d(2,:))-minval(index2d(2,:)) + 1

    end subroutine get_dimensions


    !----------------------------------------------------------------------!
    !>
    !>@brief allocate 1D and 2D arrays with regular(global) indexing i,k (1:nh,nz)
    !> they are local arrays (on the processor)
    !>@date added RHS_CONTINUOUS by F.X. Giraldo to handle CGD
    !----------------------------------------------------------------------!
    subroutine create_arrays()

        implicit none

        !Local Variables
        integer :: AllocateStatus

        save

        if(is_mod_initialized) return

        call get_dimensions()

        !> get the dimensions on this particular tile
        !>    nh_i=maxval(index2d(1,:))-minval(index2d(1,:)) + 1
        !>    nz_i=maxval(index2d(2,:))-minval(index2d(2,:)) + 1

        !> allocate the arrays
        if(allocated(qreg2)) deallocate(qreg2, qreg3)
        allocate(qreg2(nvar,nh_i),qreg3(nvar,nh_i,nz_i), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - MOD_ARRAY_UTIL: CREATE_ARRAYS **"

        !> reset the initial values
        qreg2=0.0
        qreg3=0.0
        is_mod_initialized = .true.

    end subroutine create_arrays

    !----------------------------------------------------------------------!
    !>
    !>@brief copy array contents from element based arrays to global
    !>
    !----------------------------------------------------------------------!
    subroutine copy_arrays(q,direction)

        implicit none
        save

        integer :: i,j,k,ie
        integer :: ii,kk
        integer :: ip, ip1
        integer :: direction

        real (kind=r8) :: q(nvar,npoin)

        if(.not.is_mod_initialized) call create_arrays()

        !> direction=1: copy from element/tile based arrays onto the conventional indexing
        !> direction=-1: do the reverse

        do ie=1,nelem
            do i=1,nglx
                do j=1,ngly
                    do k=1,nglz
                        ip=intma_table(i,j,k,ie)
                        ip1=intma(i,j,k,ie)
                        ii=index2d(1,ip)
                        kk=index2d(2,ip)
                        if (direction.eq.1) then
                            qreg3(:,ii,kk)=q(:,ip1)
                        else
                            q(:,ip1)=qreg3(:,ii,kk)
                        endif
                    enddo
                enddo
            enddo
        enddo
      
        if (direction.eq.1) then
            qreg2(:,:)=qreg3(:,:,1)
        endif

    end subroutine copy_arrays

    !----------------------------------------------------------------------!
    subroutine delete_arrays()

        implicit none
        save

        is_mod_initialized = .false.

    end subroutine delete_arrays

end module mod_array_util

