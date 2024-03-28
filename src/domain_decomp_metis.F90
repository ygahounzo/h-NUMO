!----------------------------------------------------------------------!
!>@brief Create the generalized adjacency matrix for a 3D Structured Cartesian grid
!>@details This algorithm constructs xadj_gen in O(nelem) time, but is not 
!> general. The assumption on the grid is that it uses GEOM_TYPE=CUBE (with no Periodic BCs) or SPHERE_HEX
!>@author James F. Kellly 15 September 2010
!>
!>@date by: S. Marras, Aug 2013
!>
!>@date F.X. Giraldo January 2014 to handle Periodic BCS (in CREATE_GEN_ADJACENCY_GENERAL)
!>
!----------------------------------------------------------------------!
subroutine create_gen_adjacency_structured(xadj_gen,adjncy_gen,wgts_gen,nface)
  
    use mod_basis, only: nglx, ngly, nglz
  
    use mod_global_grid, only: intma_g, nelem_g
  
    use mod_input, only: nelx, nely, nelz

    implicit none

    !global arrays
    integer, intent(out) :: xadj_gen(nelem_g + 1)
    integer, intent(out) :: adjncy_gen(52*nelem_g)
    integer, intent(out) :: wgts_gen(52*nelem_g)
    integer, intent(in)  :: nface

    !local arrays
    integer :: ielez, ieley, ielex, ie, je
    integer :: jelez, jeley, jelex
    integer :: jx_start, jx_end, jy_start, jy_end, jz_start, jz_end
    integer :: k, ii, i, j, iface, jface
    integer :: ivert(8), jvert(8), facemap(6,4)
  
    if(nface>1) then
     
        facemap(1,1) = 4 !face 1 left to 4
        facemap(1,2) = 2 !face 1 right to 2
        facemap(1,3) = 6 !face 1 down to 6
        facemap(1,4) = 5 !face 1 up to 5
     
        facemap(2,1) = 1 !face 2 left to 1
        facemap(2,2) = 3 !face 2 right to 3
        facemap(2,3) = 6 !face 2 down to 6
        facemap(2,4) = 5 !face 2 up to 5
     
        facemap(3,1) = 2 !face 3 left to 2
        facemap(3,2) = 4 !face 3 right to 4
        facemap(3,3) = 6 !face 3 down to 6
        facemap(3,4) = 5 !face 3 up to 5
     
        facemap(4,1) = 3 !face 4 left to 3
        facemap(4,2) = 1 !face 4 right to 1
        facemap(4,3) = 6 !face 4 down to 6
        facemap(4,4) = 5 !face 4 up to 5
     
        facemap(5,1) = 4 !face 5 left to 4
        facemap(5,2) = 2 !face 5 right to 2
        facemap(5,3) = 1 !face 5 down to 1
        facemap(5,4) = 3 !face 5 up to 3
     
        facemap(6,1) = 4 !face 6 left to 4
        facemap(6,2) = 2 !face 6 right to 2
        facemap(6,3) = 3 !face 6 down to 3
        facemap(6,4) = 1 !face 6 up to 1
     
    end if
  
    ii = 1
    do ielez = 1, nelz
        do iface = 1, nface
            do ieley = 1,nely
                do ielex = 1,nelx
                    if(nface == 1)then
                        ie = ielex + nelx*(ieley -1) + nelx*nely*(ielez -1)
                    else
                        ie = ielex + nelx*(ieley -1) + (iface-1)*nelx*nely + (ielez - 1)*nface*nelx*nely
                    end if
                    xadj_gen(ie) = ii
              
                    ivert(1) = intma_g(1,1,1,ie)
                    ivert(2) = intma_g(1,1,nglz,ie)
                    ivert(3) = intma_g(1,ngly,1,ie)
                    ivert(4) = intma_g(1,ngly,nglz,ie)
                    ivert(5) = intma_g(nglx,1,1,ie)
                    ivert(6) = intma_g(nglx,1,nglz,ie)
                    ivert(7) = intma_g(nglx,ngly,1,ie)
                    ivert(8) = intma_g(nglx,ngly,nglz,ie)
              
                    ! loop through neighbors on same face as ie
                    jx_start = max(1,ielex - 1)
                    jx_end = min(nelx,ielex + 1)
                    jy_start = max(1,ieley - 1)
                    jy_end = min(nely,ieley + 1)
                    jz_start = max(1,ielez - 1)
                    jz_end = min(nelz,ielez + 1)
              
                    do jelez = jz_start, jz_end
                        do jeley = jy_start, jy_end
                            do jelex = jx_start, jx_end
                                if(nface == 1)then
                                    je = jelex + nelx*(jeley -1) + nelx*nely*(jelez -1)
                                else
                                    je = jelex + nelx*(jeley -1) + (iface-1)*nelx*nely + (jelez - 1)*nface*nelx*nely
                                end if
                                if (je /= ie) then
                                    jvert(1) = intma_g(1,1,1,je)
                                    jvert(2) = intma_g(1,1,nglz,je)
                                    jvert(3) = intma_g(1,ngly,1,je)
                                    jvert(4) = intma_g(1,ngly,nglz,je)
                                    jvert(5) = intma_g(nglx,1,1,je)
                                    jvert(6) = intma_g(nglx,1,nglz,je)
                                    jvert(7) = intma_g(nglx,ngly,1,je)
                                    jvert(8) = intma_g(nglx,ngly,nglz,je)
                          
                                    ! Find Intersection of ivert and jvert
                                    k = 0
                                    do i = 1,8
                                        do j = 1,8
                                            if ( ivert(i) == jvert(j) ) then
                                                k = k + 1
                                                exit
                                            end if
                                        end do
                                    end do
                          
                                    if ( k > 0 ) then
                                        adjncy_gen(ii) = je
                                        wgts_gen(ii) = k
                                        ii = ii + 1
                                    end if
                          
                                end if
                            end do
                        end do
                    end do
              
                    if(nface>1)then
                 
                        if(ielex == 1) then
                    
                            ! loop through neighbors on face left of ie
                            jface = facemap(iface,1)
                    
                            if(iface<5)then
                                jx_start = nelx
                                jx_end = nelx
                                jy_start = max(1,ieley - 1)
                                jy_end = min(nely,ieley + 1)
                            elseif(iface==5)then
                                jx_start = max(1,nely-ieley)
                                jx_end = min(nely,nely-ieley+2)
                                jy_start = nely
                                jy_end = nely
                            else
                                jx_start = max(1,ieley - 1)
                                jx_end = min(nely,ieley + 1)
                                jy_start = 1
                                jy_start = 1
                            endif
                    
                            jz_start = max(1,ielez - 1)
                            jz_end = min(nelz,ielez + 1)
                    
                            do jelez = jz_start, jz_end
                                do jeley = jy_start, jy_end
                                    do jelex = jx_start, jx_end
                                        je = jelex + nelx*(jeley -1) + (jface-1)*nelx*nely + (jelez - 1)*nface*nelx*nely
                                        if (je /= ie) then
                                            jvert(1) = intma_g(1,1,1,je)
                                            jvert(2) = intma_g(1,1,nglz,je)
                                            jvert(3) = intma_g(1,ngly,1,je)
                                            jvert(4) = intma_g(1,ngly,nglz,je)
                                            jvert(5) = intma_g(nglx,1,1,je)
                                            jvert(6) = intma_g(nglx,1,nglz,je)
                                            jvert(7) = intma_g(nglx,ngly,1,je)
                                            jvert(8) = intma_g(nglx,ngly,nglz,je)
                                
                                            ! Find Intersection of ivert and jvert
                                            k = 0
                                            do i = 1,8
                                                do j = 1,8
                                                    if ( ivert(i) == jvert(j) ) then
                                                        k = k + 1
                                                        exit
                                                    end if
                                                end do
                                            end do
                                
                                            if ( k > 0 ) then
                                                adjncy_gen(ii) = je
                                                wgts_gen(ii) = k
                                                ii = ii + 1
                                            end if
                                
                                        end if
                                    end do
                                end do
                            end do
                        endif
                 
                        if(ielex == nelx) then
                    
                            ! loop through neighbors on face right of ie
                            jface = facemap(iface,2)
                    
                            if(iface<5)then
                                jx_start = 1
                                jx_end = 1
                                jy_start = max(1,ieley - 1)
                                jy_end = min(nely,ieley + 1)
                            elseif(iface==5)then
                                jx_start = max(1,ieley - 1)
                                jx_end = min(nely,ieley + 1)
                                jy_start = nely
                                jy_end = nely
                            else
                                jx_start = max(1,nely-ieley)
                                jx_end = min(nely,nely-ieley+2)
                                jy_start = 1
                                jy_start = 1
                            endif
                    
                            jz_start = max(1,ielez - 1)
                            jz_end = min(nelz,ielez + 1)
                    
                            do jelez = jz_start, jz_end
                                do jeley = jy_start, jy_end
                                    do jelex = jx_start, jx_end
                                        je = jelex + nelx*(jeley -1) + (jface-1)*nelx*nely + (jelez - 1)*nface*nelx*nely
                                        if (je /= ie) then
                                            jvert(1) = intma_g(1,1,1,je)
                                            jvert(2) = intma_g(1,1,nglz,je)
                                            jvert(3) = intma_g(1,ngly,1,je)
                                            jvert(4) = intma_g(1,ngly,nglz,je)
                                            jvert(5) = intma_g(nglx,1,1,je)
                                            jvert(6) = intma_g(nglx,1,nglz,je)
                                            jvert(7) = intma_g(nglx,ngly,1,je)
                                            jvert(8) = intma_g(nglx,ngly,nglz,je)
                                
                                            ! Find Intersection of ivert and jvert
                                            k = 0
                                            do i = 1,8
                                                do j = 1,8
                                                    if ( ivert(i) == jvert(j) ) then
                                                        k = k + 1
                                                        exit
                                                    end if
                                                end do
                                            end do
                                
                                            if ( k > 0 ) then
                                                adjncy_gen(ii) = je
                                                wgts_gen(ii) = k
                                                ii = ii + 1
                                            end if
                                
                                        end if
                                    end do
                                end do
                            end do
                        endif
                 
                        if(ieley == 1) then
                    
                            ! loop through neighbors on face below ie
                            jface = facemap(iface,3)
                    
                            if(iface==1 .or. iface==5)then
                                jx_start = max(1,ielex - 1)
                                jx_end = min(nelx,ielex + 1)
                                jy_start = nely
                                jy_end = nely
                            elseif(iface==2)then
                                jx_start = nelx
                                jx_end = nelx
                                jy_start = max(1,nelx-ielex)
                                jy_end = min(nelx,nelx-ielex+2)
                            elseif(iface==3)then
                                jx_start = max(1,nelx-ielex)
                                jx_end = min(nelx,nelx-ielex+2)
                                jy_start = 1
                                jy_end = 1
                            elseif(iface==4)then
                                jx_start = 1
                                jx_end = 1
                                jy_start = max(1,ielex - 1)
                                jy_end = min(nelx,ielex + 1)
                            else
                                jx_start = max(1,nelx-ielex)
                                jx_end = min(nelx,nelx-ielex+2)
                                jy_start = 1
                                jy_start = 1
                            endif
                    
                            jz_start = max(1,ielez - 1)
                            jz_end = min(nelz,ielez + 1)
                    
                            do jelez = jz_start, jz_end
                                do jeley = jy_start, jy_end
                                    do jelex = jx_start, jx_end
                                        je = jelex + nelx*(jeley -1) + (jface-1)*nelx*nely + (jelez - 1)*nface*nelx*nely
                                        if (je /= ie) then
                                            jvert(1) = intma_g(1,1,1,je)
                                            jvert(2) = intma_g(1,1,nglz,je)
                                            jvert(3) = intma_g(1,ngly,1,je)
                                            jvert(4) = intma_g(1,ngly,nglz,je)
                                            jvert(5) = intma_g(nglx,1,1,je)
                                            jvert(6) = intma_g(nglx,1,nglz,je)
                                            jvert(7) = intma_g(nglx,ngly,1,je)
                                            jvert(8) = intma_g(nglx,ngly,nglz,je)
                                
                                            ! Find Intersection of ivert and jvert
                                            k = 0
                                            do i = 1,8
                                                do j = 1,8
                                                    if ( ivert(i) == jvert(j) ) then
                                                        k = k + 1
                                                        exit
                                                    end if
                                                end do
                                            end do
                                
                                            if ( k > 0 ) then
                                                adjncy_gen(ii) = je
                                                wgts_gen(ii) = k
                                                ii = ii + 1
                                            end if
                                
                                        end if
                                    end do
                                end do
                            end do
                        endif
                 
                        if(ieley == nely) then
                    
                            ! loop through neighbors on face above ie
                            jface = facemap(iface,4)
                    
                            if(iface==1 .or. iface==6)then
                                jx_start = max(1,ielex - 1)
                                jx_end = min(nelx,ielex + 1)
                                jy_start = 1
                                jy_end = 1
                            elseif(iface==2)then
                                jx_start = nelx
                                jx_end = nelx
                                jy_start = max(1,ielex - 1)
                                jy_end = min(nelx,ielex + 1)
                            elseif(iface==3)then
                                jx_start = max(1,nelx-ielex)
                                jx_end = min(nelx,nelx-ielex+2)
                                jy_start = nelx
                                jy_end = nelx
                            elseif(iface==4)then
                                jx_start = 1
                                jx_end = 1
                                jy_start = max(1,nelx-ielex)
                                jy_end = min(nelx,nelx-ielex+2)
                            else
                                jx_start = max(1,nelx-ielex)
                                jx_end = min(nelx,nelx-ielex+2)
                                jy_start = nelx
                                jy_start = nelx
                            endif
                    
                            jz_start = max(1,ielez - 1)
                            jz_end = min(nelz,ielez + 1)
                    
                            do jelez = jz_start, jz_end
                                do jeley = jy_start, jy_end
                                    do jelex = jx_start, jx_end
                                        je = jelex + nelx*(jeley -1) + (jface-1)*nelx*nely + (jelez - 1)*nface*nelx*nely
                                        if (je /= ie) then
                                            jvert(1) = intma_g(1,1,1,je)
                                            jvert(2) = intma_g(1,1,nglz,je)
                                            jvert(3) = intma_g(1,ngly,1,je)
                                            jvert(4) = intma_g(1,ngly,nglz,je)
                                            jvert(5) = intma_g(nglx,1,1,je)
                                            jvert(6) = intma_g(nglx,1,nglz,je)
                                            jvert(7) = intma_g(nglx,ngly,1,je)
                                            jvert(8) = intma_g(nglx,ngly,nglz,je)
                                
                                            ! Find Intersection of ivert and jvert
                                            k = 0
                                            do i = 1,8
                                                do j = 1,8
                                                    if ( ivert(i) == jvert(j) ) then
                                                        k = k + 1
                                                        exit
                                                    end if
                                                end do
                                            end do
                                
                                            if ( k > 0 ) then
                                                adjncy_gen(ii) = je
                                                wgts_gen(ii) = k
                                                ii = ii + 1
                                            end if
                                
                                        end if
                                    end do
                                end do
                            end do
                        endif
                    endif
                end do
            end do
        end do
    end do
  
    xadj_gen(nelem_g+ 1) = ii - 1
  
end subroutine create_gen_adjacency_structured

!----------------------------------------------------------------------!
!>@brief Create the generalized adjacency matrix for GENERAL GRIDS. 
!>@details It can use  any GEOM_TYPE but we use it for GEOM_TYPE=SPHERE_ICO and any general
!> unstructured grid. It is a general routine but not very efficient.
!> The adjacency matris is then used for the creation of NBH_PROC.
!> This routine counts neighbors as any element that shares faces, edges, or verticies
!> The (sparse) adjacency matrix adjncy_gen is stored in CSR format, with indices stored
!> in xadj_gen.  The weights of the graph (corresponding to the number of corner nodes shared)
!> is stored in wgts_gen.  The wgts can be
!> 4 = face shared, 2 = edge shared, 1 = vertex shared
!> The generalized adjacency graph does not require constructing face information
!>@author 21 July 2010
!>
!>@date F.X. Giraldo on January 2014
!> to handle PERIODIC BCs
!----------------------------------------------------------------------!
subroutine create_gen_adjacency_general(xadj_gen,adjncy_gen,wgts_gen)
  
    use mod_basis, only: nglx, ngly, nglz
  
    use mod_global_grid, only: nelem_g, intma_g
  
    implicit none

    ! xadj holds indexing info
    ! adjncy holds the adjacency info
    !global arrays
    integer, intent(out) :: xadj_gen(nelem_g + 1)
    integer, dimension(52*nelem_g), intent(out) :: adjncy_gen, wgts_gen

    !local arrays
    integer ie, je, i, j, k, ii, ipp
    integer ivert(8), jvert(8)
  
    !Initialize
    xadj_gen=0; adjncy_gen=0; wgts_gen=0

    !Begin
    ii = 1
    do ie = 1,nelem_g
        xadj_gen(ie) = ii
     
        ivert(1) = intma_g(1,1,1,ie)
        ivert(2) = intma_g(1,1,nglz,ie)
        ivert(3) = intma_g(1,ngly,1,ie)
        ivert(4) = intma_g(1,ngly,nglz,ie)
        ivert(5) = intma_g(nglx,1,1,ie)
        ivert(6) = intma_g(nglx,1,nglz,ie)
        ivert(7) = intma_g(nglx,ngly,1,ie)
        ivert(8) = intma_g(nglx,ngly,nglz,ie)
     
        do je = 1,nelem_g
            if (je /= ie) then
                jvert(1) = intma_g(1,1,1,je)
                jvert(2) = intma_g(1,1,nglz,je)
                jvert(3) = intma_g(1,ngly,1,je)
                jvert(4) = intma_g(1,ngly,nglz,je)
                jvert(5) = intma_g(nglx,1,1,je)
                jvert(6) = intma_g(nglx,1,nglz,je)
                jvert(7) = intma_g(nglx,ngly,1,je)
                jvert(8) = intma_g(nglx,ngly,nglz,je)
           
                ! Find Intersection of ivert and jvert
                k = 0
                ipp=0
                do i = 1,8
                    do j = 1,8
                        if ( ivert(i) == jvert(j) ) then
                            k = k + 1
                            exit
                        end if
                    end do
                end do
           
                if ( k > 0 ) then
                    adjncy_gen(ii) = je
                    wgts_gen(ii) = k
                    ii = ii + 1
                end if
           
            end if
        end do
    end do
    xadj_gen(nelem_g+1) = ii - 1
  
end subroutine create_gen_adjacency_general

!----------------------------------------------------------------------!
!>@brief Performs Domain Decomposition in the Horz. only (2D) using METIS
!>@author 16 March 2001
!----------------------------------------------------------------------!
subroutine domain_decomp_metis_adjacency2d(global_proc, nproc)

    use mod_adjncy, only: adjncy2d_face, xadj2d_face
  
    use mod_global_grid, only:  nelem_g, nelem_s, ele_col_g
  
    implicit none
  
    !global arrays
    integer, intent(in) :: nproc
    integer, intent(out) :: global_proc(nelem_g)

    !local arrays
    integer, allocatable :: options(:)
    integer, pointer :: vwgt=>null(), vsize=>null(), adjwgts=>null()
    real(8), pointer :: tpwgts=>null(), ubvec=>null()

    integer :: edgecut
    integer :: ie, je, iface, iproc, iprocr, ib_g, istart, iend, ie_g, ie_s
    integer, dimension(:), allocatable :: global_proc2d
    integer :: ierr,ncon
  
    !init
    global_proc = 0
    ncon=1
  
    !allocate local memory
    allocate (global_proc2d(nelem_s))
    global_proc2d = 0
  
    !initialize options
    allocate(options(0:40))
    call METIS_SetDefaultOptions(options)
    options(17)=1
  
    ! Call METIS
    if (nproc > 1 ) then

        !Metis 5.0 call
        print *, "Begin Calling METIS2D v5.0"
        call METIS_PartGraphRecursive(nelem_s,ncon,xadj2d_face,adjncy2d_face,vwgt,vsize,adjwgts,&
            nproc,tpwgts,ubvec,options,edgecut,global_proc2d)
        print *, "End Calling METIS2D"

        ! Check for a METIS failure
        call  check_global_proc(global_proc2d,nelem_s,nproc,ierr)
        if (ierr == 1) then
            print *, "METIS Failed!> Doing Ad-Hoc Domain Decomp."
            call adhoc_decomp(global_proc2d,nelem_s,nproc)
        end if

        !Now construct 3D global proc
        do ie_g = 1,nelem_g
            ie_s = ele_col_g(ie_g)
            global_proc(ie_g) = global_proc2d(ie_s)
        end do

        !Now construct 3D global proc
        do ie_g = 1,nelem_g
            ie_s = ele_col_g(ie_g)
            global_proc(ie_g) = global_proc2d(ie_s)
        end do
     
    else
        global_proc = 1
    end if
  
    !free
    deallocate(global_proc2d)
    deallocate(options)
end subroutine domain_decomp_metis_adjacency2d

!----------------------------------------------------------------------!
!>@brief Domain Decomposition using the adjacency graph (3D) via METIS
!> nboun_l_interp(iproc) = is the the number of boundaries that lie on adjacent processors
!----------------------------------------------------------------------!
subroutine domain_decomp_metis_adjacency3d( global_proc, nproc)
  
    use mod_adjncy, only: adjncy_gen, xadj_gen, adjncy_face, xadj_face, nadj_face, adjwgts_gen
  
    use mod_global_grid, only:  nelem_g, intma_g, npoin_g_cg
  
    implicit none
  
    !local arrays
    integer, allocatable :: options(:)
    integer, pointer :: vwgt=>null(), vsize=>null(), adjwgts=>null()
    real(8), pointer :: tpwgts=>null(), ubvec=>null()
  
    integer nproc
    integer global_proc(nelem_g)
    integer edgecut
    integer ie, je, iface, iproc, iprocr, ib_g, istart, iend, ie_g
    integer ierr, i, ncon
  
    !init
    ncon = 1
    global_proc = 0

    !initialize options
    allocate(options(0:40))
    call METIS_SetDefaultOptions(options)
    options(17)=1

    ! Call METIS
    if (nproc > 1) then

        !Metis 5.0 call
        print *, "Begin Calling METIS3D v5.0"
        call METIS_PartGraphRecursive(nelem_g,ncon,xadj_face,adjncy_face,vwgt,vsize,adjwgts,nproc,&
            tpwgts,ubvec,options,edgecut,global_proc)
        print *, "End Calling METIS3D"


        ! Check for a METIS failure
        call  check_global_proc(global_proc,nelem_g,nproc,ierr)
        ! If METIS Failed, perform an Ad-Hoc decomposition
        if (ierr == 1) then
            print *, "METIS Failed!> Doing Ad-Hoc Domain Decomp."
            call adhoc_decomp(global_proc,nelem_g,nproc)
        end if
     
    else
        global_proc = 1
    end if
    deallocate(options)
end subroutine domain_decomp_metis_adjacency3d

!----------------------------------------------------------------------!
!>@brief This subroutine creates coord_l, index_l, intma_l, 
!> bsido_l on the head node
!>
!>@date S. Marras June 2014 to add the saved_bdy_flg to handle Dirichlet boundary conditions
!----------------------------------------------------------------------!
subroutine domain_decomp_store_local(npoin_l,npoin_l_cg, nbsido_l,npoin_l_max, &
    nface_l,nface_l_max,&
    nelem_l,nelem_l_max,nelemx_l,nelemx_l_max,nelemy_l,nelemy_l_max,&
    nelemz_l,nelemz_l_max,nboun_max,nproc,local_global_elem_l,&
    local_global_poin_l,local_global_poin_l_cg, local_global_poin_periodic_l,local_global_face_l,global_proc)
  
    use mod_basis, only: nglx, ngly, nglz, FACE_LEN

    use mod_global_grid, only:  coord_g_cg, sigma_g, index_g, intma_g, &
        bsido_g, npoin_g, nelem_g, nboun_g, nz, npoin_r, &
        saved_coord, saved_sigma, saved_index, saved_intma, &
        saved_bsido, saved_face, iperiodic_g, &
        face_g, nface_g
  
    use mod_input, only: geometry_type, bcast_type
  
    use mod_parallel,only: face_type_l
  
    use mod_utilities, only : get_unit
  
    implicit none

    !Global Variables
    integer, intent(in) :: npoin_l_max,nboun_max,nproc,nface_l_max
    integer, intent(in) :: nelem_l_max, nelemx_l_max, nelemy_l_max, nelemz_l_max
    integer, dimension(nproc), intent(in) :: nelem_l, nelemx_l, nelemy_l, nelemz_l, npoin_l, nface_l, npoin_l_cg
    integer, dimension(npoin_l_max,nproc), intent(in) :: local_global_poin_l_cg, local_global_poin_l, local_global_poin_periodic_l
    integer, dimension(nface_l_max,nproc), intent(in) :: local_global_face_l
    integer, intent(in) :: local_global_elem_l(nelem_l_max,nproc), global_proc(nelem_g)
    integer, dimension(nproc), intent(out) :: nbsido_l

    !Local Variables
    real,    dimension(:,:),     allocatable :: coord_l
    real,    dimension(:),       allocatable :: sigma_l
    integer, dimension(:,:),     allocatable :: index_l, face_l
    integer, dimension(:,:,:,:), allocatable :: intma_l
    integer :: ip, ip_g, iproc, npoin, npoin_cg, ib, ib_g, ibb, ie_g, ipp_g, ie, iep_g, iface, iface_g
    integer :: ip1_g, ip2_g, ip3_g, ip4_g
    integer :: i, j,k, l,nelem,nelemx,nelemy,nelemz, nh,ierr,nface, itype
    integer :: iloop
    integer :: AllocateStatus
    integer, dimension (:,:), allocatable :: local_global_bsido, bsido_l
    integer, dimension(:),    allocatable :: global_local_elem, global_local_poin
    integer :: nbsido_l_max
    integer :: iunit
    character(len=14) :: fnp
  
    !Create Local Intma
    allocate(global_local_poin(npoin_g), global_local_elem(nelem_g), &
        local_global_bsido(nboun_g,nproc), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop "** Not Enough Memory - Domain_Decomp_Metis:Domain_Decomp_Store_Local 1**"
  
    !Initialize allocated memory
    global_local_poin  = 0
    global_local_elem  = 0
    local_global_bsido = 0
  
    ! First, count number of boundary points owned by each
    ! processor.
    nbsido_l = 0
    do ib = 1, nboun_g
        ie_g = bsido_g(5,ib)
        iproc = global_proc(ie_g)
        nbsido_l(iproc) = nbsido_l(iproc) + 1
        ibb = nbsido_l(iproc)
        local_global_bsido(ibb,iproc) = ib
    end do
  
    nbsido_l_max = maxval(nbsido_l)
    allocate(bsido_l(6,nbsido_l_max))
    !Initialize allocated memory
    bsido_l = 0

    ! Construct Local Coordinates Here
    !if (grid_gen > 0) then
    if (bcast_type=='mpi') then
        allocate(saved_coord(3,maxval(npoin_l_cg),nproc),&
            saved_sigma(maxval(npoin_l_cg),nproc),&
            saved_index(2,maxval(npoin_l_cg),nproc),&
            saved_intma(nglx,ngly,nglz,maxval(nelem_l),nproc),&
            saved_bsido(6,maxval(nbsido_l),nproc), &
            saved_face(FACE_LEN,maxval(nface_l),nproc),&
            stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Domain_Decomp_Metis:Domain_Decomp_Store_Local 2**"
        !   These are deallocated in mod_grid_create()
        !   after transferring data to all processors
     
        !initialize allocated variables:
        saved_coord = 0.0
        saved_index = 0
        saved_intma = 0
        saved_bsido = 0
        saved_face = 0

    endif
  
    ! Construct  Global/Local mappings
    global_local_elem = 0
    do iproc = 1,nproc
        do ie = 1,nelem_l(iproc)
            ie_g = local_global_elem_l(ie,iproc)
            global_local_elem(ie_g) =  ie
        end do
    end do
  
    do iproc = 1,nproc
        ! Create Name for File
        write(fnp,'(A5,I5.5,A4)') 'grid_',iproc,'.gri'

        npoin = npoin_l(iproc)
        npoin_cg = npoin_l_cg(iproc)
        nelem = nelem_l(iproc)
        nelemx= nelemx_l(iproc)
        nelemy= nelemy_l(iproc)
        nelemz= nelemz_l(iproc)
        nface = nface_l(iproc)

        allocate(coord_l(3,npoin_cg), sigma_l(npoin_cg), index_l(2,npoin_cg), &
            intma_l(nglx,ngly,nglz,nelem), &
            face_l(FACE_LEN,nface),&
            stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Domain_Decomp_Metis:Domain_Decomp_Store_Local 3**"
     
        !Initialize allocated variables
        coord_l = 0.0
        sigma_l = 0.0
        index_l = 0
        intma_l = 0
        face_l = 0
     
        !Consturct Local Grid
        do ip = 1,npoin_cg
            !Get Global Grid Point
            ip_g = local_global_poin_l_cg(ip,iproc)
            coord_l(:,ip) = coord_g_cg(:,ip_g)
            sigma_l(ip) = sigma_g(ip_g)
            global_local_poin(ip_g) = ip
        end do
     
        ! Assume no vertical decomposition for physics indexing!!!!!!!!!>
        if(geometry_type(1:4)=='cube')then
            nh = npoin_cg/nz
        else
            nh = npoin_cg/npoin_r
        endif
     
        ! Create a Local Intma
        do ie = 1,nelem_l(iproc)
            do k = 1,nglz
                do j = 1,ngly
                    do i = 1,nglx
                        ie_g = local_global_elem_l(ie,iproc)
                        ip_g = ( intma_g(i,j,k,ie_g) )
                        ip  = global_local_poin(ip_g)
                        intma_l(i,j,k,ie) = ip
                        index_l(2,ip) = index_g(2,ip_g)
                        index_l(1,ip) = ip - (index_l(2,ip)-1)*nh
                    end do
                end do
            end do
        end do
     
        ! Create BSIDO
        if (nbsido_l(iproc) > 0) then
        
            ! Now Populate Bsido_l
            do ib = 1,nbsido_l(iproc)
                ! Global BSIDO ID
                ib_g = local_global_bsido(ib,iproc)
                bsido_l(1:4,ib) = global_local_poin(bsido_g(1:4,ib_g))
                bsido_l(5,ib) = global_local_elem(bsido_g(5,ib_g))
                ! Assign Boundary Condition
                bsido_l(6,ib) = bsido_g(6,ib_g)
            end do
        
        end if


        do iface = 1, nface_l(iproc)
            iface_g = local_global_face_l(iface,iproc)
            face_l(1:6,iface) = face_g(1:6,iface_g)
            ! Type of Face
            itype = face_type_l(iface,iproc)
            ! Internal Face
            if (itype == 1) then
                face_l(7,iface) = global_local_elem(face_g(7,iface_g))
                face_l(8,iface) = global_local_elem(face_g(8,iface_g))
               ! Inter-Processor Face (element lies to the left)
            else if (itype == 2 ) then
                face_l(7,iface) = global_local_elem(face_g(7,iface_g))
                face_l(8,iface) = 0
               ! Inter-Processor Face (element lies to the right)
            else if (itype == 3) then
                face_type_l(iface,iproc) = 2
                face_l(5,iface) = face_g(6,iface_g)
                face_l(6,iface) = face_g(5,iface_g)
                face_l(7,iface) = global_local_elem(face_g(8,iface_g))
                face_l(8,iface) = 0
               ! Boundary Face
            else if (itype == 4) then
                face_l(7,iface) = global_local_elem(face_g(7,iface_g))
                face_l(8,iface) = face_g(8,iface_g)
            end if
        end do

        ! Communicate
        if (bcast_type=='file') then
            iunit=get_unit()
            open(iunit,file=fnp,form="unformatted")
            write(iunit) coord_l
            write(iunit) sigma_l
            write(iunit) index_l
            write(iunit) intma_l
            write(iunit) bsido_l(:,1:nbsido_l(iproc))
            write(iunit) face_l
            write(iunit) face_type_l(:,1:nface_l(iproc))
            close(iunit)
        elseif (bcast_type=='mpi') then
            do i = 1,3
                do j = 1,npoin_cg
                    saved_coord(i,j,iproc) = coord_l(i,j)
                enddo
            enddo
            do j = 1,npoin_cg
                saved_index(1,j,iproc) = index_l(1,j)
                saved_index(2,j,iproc) = index_l(2,j)
                saved_sigma(j,iproc)   = sigma_l(j)
            enddo
            do j = 1,nface_l(iproc)
                saved_face(:,j,iproc) = face_l(:,j)
            enddo
            do i = 1,nglx
                do j = 1,ngly
                    do k = 1,nglz
                        do l = 1,nelem_l(iproc)
                            saved_intma(i,j,k,l,iproc) = intma_l(i,j,k,l)
                        enddo
                    enddo
                enddo
            enddo
            do i = 1,6
                do j = 1,nbsido_l(iproc)
                    saved_bsido(i,j,iproc) = bsido_l(i,j)
                enddo
            enddo
        else
            print*,"Error in domain_decomp(): bcast_type is not valid."
            print*,"It should be either 'mpi' or 'file'."
            print*,"Instead it is '",bcast_type,"'."
            stop
        endif
     
        !Free local memory
        deallocate(coord_l,sigma_l,index_l,intma_l,face_l)
     
    end do !iproc
    !end if

    deallocate(bsido_l)
    deallocate(global_local_poin, global_local_elem)
  
end subroutine domain_decomp_store_local

!----------------------------------------------------------------------!
subroutine count_faces(nface_l,nface_boundary)
  
    use mod_adjncy, only: adjncy_face, xadj_face, nadj_face

    use mod_input, only: decomp_type
  
    use mod_global_grid, only: nelem_g, nboun_g, bsido_g
  
    use mod_parallel, only: global_proc, nproc
  
    implicit none
  
    !global arrays
    integer, intent(out) :: nface_l(nproc), nface_boundary(nproc)

    !local arrays
    integer :: nface_int(nproc)
    integer :: nface_iproc(nproc)
    integer :: nface_boun(nproc)
    integer :: ie, je, iproc, jproc, ied, istart, iend, ib, ibc, ipp
  
    ! Number of Faces = Internal Faces + Inter-Processor Faces + Physical Boundaries
    nface_l=0
    nface_boundary=0
    nface_int = 0
    nface_iproc = 0
    nface_boun = 0

    ! Count
    do ie = 1,nelem_g
        iproc = global_proc(ie)
        istart = xadj_face(ie)
        if (ie < nelem_g) then
            iend = xadj_face(ie + 1) - 1
        else
            iend = xadj_face(ie + 1)
        end if
        do ied = istart,iend
            je = adjncy_face(ied)
            jproc = global_proc(je)
            if (iproc == jproc) then
                nface_int(iproc) = nface_int(iproc) + 1
            else
                nface_iproc(iproc) = nface_iproc(iproc) + 1
            end if
        end do
    end do
  
    ! Get rid of redundant faces
    nface_int = nface_int/2
  
    !Now count Physical Boundaries
    do ib = 1,nboun_g
        ie = bsido_g(5,ib)
        ibc= bsido_g(6,ib)
        iproc = global_proc(ie)

        if (ibc /= 3 .or. decomp_type(1:4) == 'geom') then
            nface_boun(iproc) = nface_boun(iproc) + 1
        end if
     
    end do
  
    ! Total Number of Faces
    nface_l = nface_int + nface_iproc + nface_boun
    ! Total Number of Boundary Faces
    nface_boundary =  nface_iproc ! + nface_boun
  
  !Write Output
  
end subroutine count_faces

!----------------------------------------------------------------------!
subroutine count_elements(nelem_l)
  
    use mod_global_grid, only: nelem_g

    use mod_parallel, only: global_proc, nproc

    implicit none
  
    !global arrays
    integer, intent(out) :: nelem_l(nproc)

    !local arrays
    integer :: ie, iproc
  
    !Now Determine how many elements each processor owns
    nelem_l = 0
    do ie = 1,nelem_g
        iproc = global_proc(ie)
        nelem_l(iproc) = nelem_l(iproc) + 1
    end do
  
end subroutine count_elements

!----------------------------------------------------------------------!
subroutine count_elements_xyz(nelemx_l,nelemy_l,nelemz_l, nprocx, nprocy, nprocz)
  
    use mod_global_grid, only: nelem_g

    use mod_parallel, only: global_proc, nproc

    use mod_input, only: nelx, nely, nelz
  
    implicit none
  
    !global arrays
    integer, intent(out) :: nelemx_l(nproc),nelemy_l(nproc),nelemz_l(nproc)
    integer, intent(in)  :: nprocx, nprocy, nprocz

    !local arrays
    integer nelx_local, nely_local, nelz_local
    integer ie, iproc, iex,iey,iez
  
    nelx_local = nelx/nprocx
    nely_local = nely/nprocy
    nelz_local = nelz/nprocz
  
    !Now Determine how many elements each processor owns
    nelemx_l = 0
    nelemy_l = 0
    nelemz_l = 0
  
    ie = 0
    do ie = 1,nelem_g
        iproc = global_proc(ie)
     
        nelemx_l(iproc) = nelx_local
        nelemy_l(iproc) = nely_local
        nelemz_l(iproc) = nelz_local
    end do
  
end subroutine count_elements_xyz

!----------------------------------------------------------------------!
subroutine check_global_proc(global_proc,nelem,nproc,ierr)
  
    implicit none

    !global arrays
    integer, intent(in) :: nelem, nproc
    integer, intent(in) :: global_proc(nelem)
    integer, intent(out):: ierr

    !local arrays
    integer :: check(nproc)
    integer :: ie, iproc
  
    check = 0; ierr=0

    do ie = 1,nelem
        iproc = global_proc(ie)
        check(iproc) = 1
    end do
  
    do iproc = 1,nproc
        if (check(iproc) == 0) then
            ierr = 1
        end if
    end do
  
end subroutine check_global_proc

!----------------------------------------------------------------------!
subroutine adhoc_decomp(global_proc,nelem,nproc)
  
    implicit none

    !global arrays
    integer, intent(in) :: nelem, nproc
    integer, intent(out) :: global_proc(nelem)

    !local arrays
    integer :: i, i1, i2, nrem, ndiv, ndivp1
    integer :: j, iproc

    !initialize
    global_proc=0

    ndiv = nelem/nproc
    ndivp1 = ndiv + 1
  
    nrem = mod(nelem,nproc)
    !nrem2 = nproc - nrem
  
    !The first nrem procs get ndivp1 elements each
    do iproc = 1,nrem
        do i = 1,ndivp1
            j = ndivp1*(iproc - 1) + i
            global_proc(j) = iproc
        end do
    end do
  
    !The remaining processors get ndiv elements each
    do iproc = (nrem + 1),nproc
        do i = 1,ndiv
            j = ndivp1*nrem + ndiv*(iproc -nrem - 1) + i
            global_proc(j) = iproc
        end do
    end do
  
end subroutine adhoc_decomp


subroutine out_grid_info()
  
    use mod_parallel, only: nelem_l, npoin_l
  
end subroutine out_grid_info
