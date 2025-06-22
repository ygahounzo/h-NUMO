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
  
    use mod_input, only: filter_weight_type, filter_basis_type, &
        filter_mux, filter_muy, filter_muz, &
        nopx, nopy, nopz, is_non_conforming_flg, is_mlswe, dg_integ_exact

    public :: &
        mod_basis_create, &
        nop, nopx, nopy, nopz, &
        ngl, ngl2, nglx, ngly, nglz, npts, nsubcells,&
        psi, psix, psiy, psiz, &
        dpsi, dpsix, dpsiy, dpsiz, &
        dpsi_tr, dpsix_tr, dpsiy_tr, dpsiz_tr, &
        xgl, xglx, xgly, xglz, &
        wgl, wglx, wgly, wglz, &
        f, fx, fy, fz, &
        f_tr, fx_tr, fy_tr, fz_tr, is_2d, FACE_LEN, FACE_CHILDREN, CELL_CHILDREN, &
        P4EST_FACES, P8EST_EDGES, &
        nqx, nqy, nqz, npts_quad, nq, nsubcells_quad, nq2, &
        wnq, wnqx, wnqy, wnqz, &
        psiq, psiqx, psiqy, psiqz, &
        dpsiq, dpsiqx, dpsiqy, dpsiqz, &
        xnq, xnqx, xnqy, xnqz
  
    private

    ! module variables and parameters
    integer :: nop, ngl, nglx, ngly, nglz, npts, nsubcells, FACE_LEN, FACE_CHILDREN
    integer :: CELL_CHILDREN, P4EST_FACES, P8EST_EDGES
    integer :: nqx, nqy, nqz, npts_quad, nq, nsubcells_quad, nq2
    logical :: is_2d
    real, dimension(:,:), allocatable :: psi, psix, psiy, psiz, &
        dpsi, dpsix, dpsiy, dpsiz, &
        dpsi_tr, dpsix_tr, dpsiy_tr, dpsiz_tr, &
        f, fx, fy, fz, &
        f_tr, fx_tr, fy_tr, fz_tr
    real, dimension(:),   allocatable :: xgl, xglx, xgly, xglz, &
        wgl, wglx, wgly, wglz
    real, dimension(:,:), allocatable :: psiq, psiqx, psiqy, psiqz, &
        dpsiq, dpsiqx, dpsiqy, dpsiqz, dpsiqx_tr
    real, dimension(:),   allocatable :: xnq, xnqx, xnqy, xnqz, wnq, wnqx, wnqy, wnqz

  
    contains
  
    subroutine mod_basis_create(nopx_loc,nopy_loc,nopz_loc)

        use mod_legendre, only : legendre_gauss_lobatto, legendre_basis, lagrange_basis

        implicit none

        !global
        integer, intent(in) :: nopx_loc, nopy_loc, nopz_loc

        !local
        integer :: AllocateStatus, i, j
 
        nopx=nopx_loc
        nopy=nopy_loc
        nopz=nopz_loc
        nop=max(nopx,nopy,nopz)
        nglx=nopx+1
        ngly=nopy+1
        nglz=nopz+1
        ngl=max(nglx,ngly,nglz)
        ngl2=ngl*ngl
        npts=nglx*ngly*nglz
        nsubcells=max(nglx-1,1)*max(ngly-1,1)*max(nglz-1,1)
        
        if(dg_integ_exact) then 
            nqx = 2*nopx + 1
            nqy = 2*nopy + 1
            nqz = 2*nopz + 1
        else 
            nqx = 2*nopx - 1
            nqy = 2*nopy - 1
            nqz = 2*nopz + 1
        end if

        if(is_mlswe) nqz = 1

        nq = max(nqx,nqy,nqz)
        nq2 = nq*nq
        nsubcells_quad = max(nqx-1,1)*max(nqy-1,1)*max(nqz-1,1);

        npts_quad = nqx*nqy*nqz
        

        print*,"mod_basis: nglx, ngly, nglz",nglx,ngly,nglz
        if(nglx == 1 .or. ngly == 1 .or. nglz == 1) then
            is_2d = .true.
            FACE_CHILDREN = 2
            CELL_CHILDREN = 4
            P4EST_FACES   = 4
            P8EST_EDGES   = 0
        else
            is_2d = .false.
            FACE_CHILDREN = 4
            CELL_CHILDREN = 8
            P4EST_FACES   = 6
            P8EST_EDGES   = 12
        endif
        print*,"mod_basis: is_2d", is_2d
        
        if(is_non_conforming_flg == 0) then
            FACE_LEN = 8
        else
            FACE_LEN = 7 + FACE_CHILDREN
        endif

        allocate( psi(ngl,ngl), psix(nglx,nglx), psiy(ngly,ngly), psiz(nglz,nglz), &
            dpsi(ngl,ngl), dpsix(nglx,nglx), dpsiy(ngly,ngly), dpsiz(nglz,nglz), &
            dpsi_tr(ngl,ngl), dpsix_tr(nglx,nglx), dpsiy_tr(ngly,ngly), dpsiz_tr(nglz,nglz), &
            xgl(ngl), xglx(nglx), xgly(ngly), xglz(nglz), &
            wgl(ngl), wglx(nglx), wgly(ngly), wglz(nglz), &
            f(ngl,ngl), fx(nglx,nglx), fy(ngly,ngly), fz(nglz,nglz), &
            f_tr(ngl,ngl), fx_tr(nglx,nglx), fy_tr(ngly,ngly), fz_tr(nglz,nglz),stat=AllocateStatus)

        if(is_mlswe) then ! added by Yao Gahounzo
            allocate( psiq(ngl,nq), psiqx(nglx,nqx), psiqy(ngly,nqy), psiqz(nglz,nqz), &
                dpsiq(ngl,nq), dpsiqx(nglx,nqx), dpsiqy(ngly,nqy), dpsiqz(nglz,nqz), &
                xnq(nq), xnqx(nqx), xnqy(nqy), xnqz(nqz), &
                wnq(nq), wnqx(nqx), wnqy(nqy), wnqz(nqz), dpsiqx_tr(nqx,nglx), stat=AllocateStatus)

                print*, "Allocation Created"
        endif 
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Basis **"
    
        !Generate Gauss-Lobatto Points for NGL
        call legendre_gauss_lobatto(ngl,xgl,wgl)
        call legendre_gauss_lobatto(nglx,xglx,wglx)
        call legendre_gauss_lobatto(ngly,xgly,wgly)
        call legendre_gauss_lobatto(nglz,xglz,wglz)

        !Construct Legendre Cardinal Bases at Gauss-Lobatto Points
        call legendre_basis(ngl,xgl,psi,dpsi)
        call legendre_basis(nglx,xglx,psix,dpsix)
        call legendre_basis(ngly,xgly,psiy,dpsiy)
        call legendre_basis(nglz,xglz,psiz,dpsiz)

        if(is_mlswe) then ! added by Yao Gahounzo

            call lagrange_basis(ngl,xgl,nq,xnq,wnq,psiq,dpsiq)
            call lagrange_basis(nglx,xglx,nqx,xnqx,wnqx,psiqx,dpsiqx)
            call lagrange_basis(ngly,xgly,nqy,xnqy,wnqy,psiqy,dpsiqy)
            call lagrange_basis(nglz,xglz,nqz,xnqz,wnqz,psiqz,dpsiqz)

            wnqz = 1.0
            wglz = 1.0

        end if

        !Construct Transpose of Differentiation Matrix to use with MXM routine
        dpsi_tr=transpose(dpsi)
        dpsix_tr=transpose(dpsix)
        dpsiy_tr=transpose(dpsiy)
        dpsiz_tr=transpose(dpsiz)

        ! Initialize Filter
        call filter_init(f,wgl,xgl,ngl,ngl-1,filter_mux,filter_weight_type,filter_basis_type)
        call filter_init(fx,wglx,xglx,nglx,nglx-1,filter_mux,filter_weight_type,filter_basis_type)
        call filter_init(fy,wgly,xgly,ngly,ngly-1,filter_muy,filter_weight_type,filter_basis_type)
        call filter_init(fz,wglz,xglz,nglz,nglz-1,filter_muz,filter_weight_type,filter_basis_type)

        !Store Filter Transpose for MXM
        f_tr=transpose(f)
        fx_tr=transpose(fx)
        fy_tr=transpose(fy)
        fz_tr=transpose(fz)
    
        
    end subroutine mod_basis_create
  
end module mod_basis
