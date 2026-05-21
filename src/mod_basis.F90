!----------------------------------------------------------------------
!>@brief This module builds the basis functions: derivatives and filters
!>@author  Francis X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@date 1/2014 F.X. Girlado
!> Modified to use variable order polynomials in all directions
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: April 03, 2023
!       modified for exact integration ( added nq ad is_mlswe condtions)
!----------------------------------------------------------------------
module mod_basis

    type basis
        ! module variables and parameters
        integer :: nop, ngl, ngl2, nglx, ngly, nglz, npts, nsubcells, FACE_LEN, FACE_CHILDREN
        integer :: CELL_CHILDREN, P4EST_FACES, P8EST_EDGES
        integer :: nqx, nqy, nqz, npts_quad, nq, nsubcells_quad, nq2
        logical :: is_2d
        real, dimension(:,:), allocatable :: psi, psix, psiy, psiz, &
            dpsi, dpsix, dpsiy, dpsiz, &
            dpsi_tr, dpsix_tr, dpsiy_tr, dpsiz_tr
        real, dimension(:),   allocatable :: xgl, xglx, xgly, xglz, &
            wgl, wglx, wgly, wglz
        real, dimension(:,:), allocatable :: psiq, psiqx, psiqy, psiqz, &
            dpsiq, dpsiqx, dpsiqy, dpsiqz, dpsiqx_tr
        real, dimension(:),   allocatable :: xnq, xnqx, xnqy, xnqz, wnq, wnqx, wnqy, wnqz
        integer :: nopx, nopy, nopz
    end type basis

    public :: &
        mod_basis_create, &
        basis
  
    private

  
    contains
  
    subroutine mod_basis_create(inp,b)

        use mod_legendre, only : legendre_gauss_lobatto, legendre_basis, lagrange_basis
        use mod_input,    only : input

        implicit none

        !global
        type(input), intent(in)    :: inp
        type(basis), intent(inout) :: b

        !global
        

        !local
        integer :: AllocateStatus, i, j
 
        b%nopx = inp%nopx
        b%nopy = inp%nopy
        b%nopz = inp%nopz
        b%nop = max(b%nopx,b%nopy,b%nopz)
        b%nglx = b%nopx+1
        b%ngly = b%nopy+1
        b%nglz = b%nopz+1
        b%ngl = max(b%nglx,b%ngly,b%nglz)
        b%ngl2 = b%ngl*b%ngl
        b%npts=b%nglx*b%ngly*b%nglz
        b%nsubcells=max(b%nglx-1,1)*max(b%ngly-1,1)*max(b%nglz-1,1)
        
        if(inp%dg_integ_exact) then 
            b%nqx = 2*b%nopx + 1
            b%nqy = 2*b%nopy + 1
            b%nqz = 2*b%nopz + 1
        else 
            b%nqx = 2*b%nopx - 1
            b%nqy = 2*b%nopy - 1
            b%nqz = 1
        end if

        if(inp%is_mlswe) b%nqz = 1

        b%nq = max(b%nqx,b%nqy,b%nqz)
        b%nq2 = b%nq*b%nq
        b%nsubcells_quad = max(b%nqx-1,1)*max(b%nqy-1,1)*max(b%nqz-1,1)

        b%npts_quad = b%nqx*b%nqy*b%nqz
        

        print*,"mod_basis: nglx, ngly, nglz",b%nglx,b%ngly,b%nglz
        if(b%nglx == 1 .or. b%ngly == 1 .or. b%nglz == 1) then
            b%is_2d = .true.
            b%FACE_CHILDREN = 2
            b%CELL_CHILDREN = 4
            b%P4EST_FACES   = 4
            b%P8EST_EDGES   = 0
        else
            b%is_2d = .false.
            b%FACE_CHILDREN = 4
            b%CELL_CHILDREN = 8
            b%P4EST_FACES   = 6
            b%P8EST_EDGES   = 12
        endif
        print*,"mod_basis: is_2d", b%is_2d
        
        if(inp%is_non_conforming_flg == 0) then
            b%FACE_LEN = 8
        else
            b%FACE_LEN = 7 + b%FACE_CHILDREN
        endif

        allocate( b%psi(b%ngl,b%ngl), b%psix(b%nglx,b%nglx), b%psiy(b%ngly,b%ngly), b%psiz(b%nglz,b%nglz), &
            b%dpsi(b%ngl,b%ngl), b%dpsix(b%nglx,b%nglx), b%dpsiy(b%ngly,b%ngly), b%dpsiz(b%nglz,b%nglz), &
            b%dpsi_tr(b%ngl,b%ngl), b%dpsix_tr(b%nglx,b%nglx), b%dpsiy_tr(b%ngly,b%ngly), b%dpsiz_tr(b%nglz,b%nglz), &
            b%xgl(b%ngl), b%xglx(b%nglx), b%xgly(b%ngly), b%xglz(b%nglz), &
            b%wgl(b%ngl), b%wglx(b%nglx), b%wgly(b%ngly), b%wglz(b%nglz), &
            stat=AllocateStatus)

        if(inp%is_mlswe) then
            allocate( b%psiq(b%ngl,b%nq), b%psiqx(b%nglx,b%nqx), b%psiqy(b%ngly,b%nqy), b%psiqz(b%nglz,b%nqz), &
                b%dpsiq(b%ngl,b%nq), b%dpsiqx(b%nglx,b%nqx), b%dpsiqy(b%ngly,b%nqy), b%dpsiqz(b%nglz,b%nqz), &
                b%xnq(b%nq), b%xnqx(b%nqx), b%xnqy(b%nqy), b%xnqz(b%nqz), &
                b%wnq(b%nq), b%wnqx(b%nqx), b%wnqy(b%nqy), b%wnqz(b%nqz), b%dpsiqx_tr(b%nqx,b%nglx), stat=AllocateStatus)

                print*, "Allocation Created"
        endif 
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Basis **"
    
        !Generate Gauss-Lobatto Points for NGL
        call legendre_gauss_lobatto(b%ngl,b%xgl,b%wgl)
        call legendre_gauss_lobatto(b%nglx,b%xglx,b%wglx)
        call legendre_gauss_lobatto(b%ngly,b%xgly,b%wgly)
        call legendre_gauss_lobatto(b%nglz,b%xglz,b%wglz)

        !Construct Legendre Cardinal Bases at Gauss-Lobatto Points
        call legendre_basis(b%ngl,b%xgl,b%psi,b%dpsi)
        call legendre_basis(b%nglx,b%xglx,b%psix,b%dpsix)
        call legendre_basis(b%ngly,b%xgly,b%psiy,b%dpsiy)
        call legendre_basis(b%nglz,b%xglz,b%psiz,b%dpsiz)

        if(inp%is_mlswe) then ! added by Yao Gahounzo

            call lagrange_basis(b%ngl,b%xgl,b%nq,b%xnq,b%wnq,b%psiq,b%dpsiq)
            call lagrange_basis(b%nglx,b%xglx,b%nqx,b%xnqx,b%wnqx,b%psiqx,b%dpsiqx)
            call lagrange_basis(b%ngly,b%xgly,b%nqy,b%xnqy,b%wnqy,b%psiqy,b%dpsiqy)
            call lagrange_basis(b%nglz,b%xglz,b%nqz,b%xnqz,b%wnqz,b%psiqz,b%dpsiqz)

            b%wnqz = 1.0
            b%wglz = 1.0

        end if

        !Construct Transpose of Differentiation Matrix to use with MXM routine
        b%dpsi_tr=transpose(b%dpsi)
        b%dpsix_tr=transpose(b%dpsix)
        b%dpsiy_tr=transpose(b%dpsiy)
        b%dpsiz_tr=transpose(b%dpsiz)
        
    end subroutine mod_basis_create
  
end module mod_basis
