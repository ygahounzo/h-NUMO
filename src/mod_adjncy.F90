!----------------------------------------------------------------------!
!>@brief This module creates adjacency graphs in CSR format
!> for the global grid and for the processor space.
!>@details These adjacency graphs are later used in METIS-based
!> domain decomposition and in constructing inter-processor
!> communication
!>@author James F. Kelly 22 July 2010
!----------------------------------------------------------------------!
module mod_adjncy
  
    use mod_global_grid, only: nelem_g, npoin_g_cg, iperiodic_g
  
    use mod_input, only: geometry_type, decomp_type, icase
  
    public :: &
        mod_adjncy_create_adjncy_gen, &
        mod_adjncy_create_adjncy_face, &
        mod_adjncy_convert_adjncy2d_face, &
        mod_adjncy_create_adjncy_proc, &
        adjncy_face, adjncy2d_face, &
        xadj_face, xadj2d_face, &
        nadj_face, &
        nadj2d_face, &
        adjncy_gen, &
        xadj_gen, &
        adjwgts_gen, &
        nadj_gen, &
        adjncy_proc, &
        xadj_proc, &
        adjwgts_proc, &
        nadj_proc
  
    private
    !----------------------------------------------------------------------!
    integer nadj_face, nadj2d_face, nadj_gen, nadj_proc
    integer, dimension(:), allocatable :: adjncy_face, xadj_face
    integer, dimension(:), allocatable :: adjncy2d_face, xadj2d_face
    integer, dimension(:), allocatable :: adjncy_gen, xadj_gen, adjwgts_gen
    integer, dimension(:), allocatable :: adjncy_proc, xadj_proc, adjwgts_proc
  !-----------------------------------------------------------------------
  
contains
  
    !----------------------------------------------------------------------!
    !>@brief This routine builds the General Adjacency Matrix.
    !>@details Whereby for each element E, we check to see which other elements also
    !> use the 8 vertices of the element E. Therefore, we include face, edge and
    !> corner neighbors (for a structured hex grid, each element will have 26 neighbors)
    !>@author James F. Kelly 22 July 2010
    !----------------------------------------------------------------------!
    subroutine mod_adjncy_create_adjncy_gen()
    
        implicit none
    
        integer :: nface, e, i, istart, iend, AllocateStatus
    
        ! Maximum number of faces
        allocate(adjncy_gen(52*nelem_g), adjwgts_gen(52*nelem_g), &
            xadj_gen(nelem_g + 1), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Adjncy_Create_Adjncy_Gen **"
    
        ! Create a generalized Adjacency Matrix in CSR format
        if (geometry_type == 'cube' .or. geometry_type(1:10) == 'sphere_hex') then
            if (geometry_type == 'cube') then
                nface = 1
            else if (geometry_type(1:10) == 'sphere_hex') then
                nface = 6
            endif
            if (decomp_type(1:4) == 'geom') then
                call create_gen_adjacency_structured(xadj_gen,adjncy_gen,adjwgts_gen,nface)
            else if (decomp_type(1:5) == 'metis') then
                print*,' enter CREATE_GEN_ADJACENCY_GENERAL'
                call create_gen_adjacency_general(xadj_gen,adjncy_gen,adjwgts_gen)
                print*,' exit CREATE_GEN_ADJACENCY_GENERAL'
            end if
           ! else if (icase == 7) then
           !       call create_gen_adjacency_per(xadj_gen,adjncy_gen,adjwgts_gen,.false., .false., .false.)
        else
            call create_gen_adjacency_general(xadj_gen,adjncy_gen,adjwgts_gen)
        end if
        ! Get the number of "edges" in the generalized adjacency matrix
        ! "edges" here refers to all communications (face, edge, corners)
        nadj_gen = xadj_gen(nelem_g + 1)
        print *, "Number of Edges: ", nadj_gen

      !Write Out Data


    
    end subroutine mod_adjncy_create_adjncy_gen
  
    !----------------------------------------------------------------------!
    !>@brief This routine extracts the face neighbors only from the General
    !> AdjacencyGraph
    !>@author James F. Kelly 22 July 2010
    !>@date F.X. Giraldo January 2014
    !----------------------------------------------------------------------!
    subroutine mod_adjncy_create_adjncy_face()
    
        implicit none

        integer :: ied, ied1, ie_g, istart, iend, iface, e, i, AllocateStatus
    
        ! Count the number of edges (in other words, faces between elements)
        nadj_face = 0
        do ied = 1, nadj_gen
            if (adjwgts_gen(ied) == 4) then
                nadj_face = nadj_face + 1
            end if
        end do
    
        ! Allocate Memory
        allocate(adjncy_face(nadj_face), xadj_face(nelem_g + 1), stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_Ajncy_Create_Adjncy_Face **"
    
        ied1 = 1
        do ie_g = 1,nelem_g
            xadj_face(ie_g) = ied1
            istart =  xadj_gen(ie_g)
            if (ie_g < nelem_g) then
                iend = xadj_gen(ie_g + 1) - 1
            else
                iend = xadj_gen(ie_g + 1)
            end if
            do ied= istart,iend
                if (adjwgts_gen(ied) == 4) then
                    adjncy_face(ied1) = adjncy_gen(ied)
                    ied1 = ied1 + 1
                end if
            end do
        end do
        xadj_face(nelem_g + 1) = nadj_face
        print *, "Nadj_gen = ", nadj_gen
        print *, "Nadj_face= ", ied1, nadj_face
    
        !Write Out Data
        open(1,file='adjncy_face.out')
        do e=1,nelem_g
            istart=xadj_face(e)
            if (e < nelem_g) then
                iend=xadj_face(e+1) - 1
            else
                iend=xadj_face(e+1)
            end if
            do i=istart,iend
                write(1,*)e,adjncy_face(i)
            end do
        enddo
        close(1)

    end subroutine mod_adjncy_create_adjncy_face
  
    !----------------------------------------------------------------------!
    !>@brief Creates the Processor Adjacency Matrix from a METIS Decomposition
    !>@details type = cg ==> Processor Adj. Matrix for Continuous Galerkin (CG), which
    !>includes vertex, side, and face connectivity
    !>type = dg ==> Processor Adj. Matrix for Discontinous Galerkin (DG), which
    !>only includes face connectivity
    !----------------------------------------------------------------------!
    subroutine mod_adjncy_create_adjncy_proc(global_proc,nproc,type)
    
        implicit none

        !global
        integer, intent(in) :: nproc
        integer, intent(in) :: global_proc(nelem_g)
        character, intent(in) :: type*2

        !local
        integer :: ie_g, iproc, jproc, istart, iend, iface
        integer :: je_g, ied
        ! Allocate the full matrix, then compress later
        ! Probably a way to construct adj_proc in CSR format, however ....
        integer, dimension(:,:), allocatable :: A_proc
    
        !Temporary Full Matrix
        allocate(A_proc(nproc,nproc))
    
        A_proc = 0
    
        !DG Processor Adjacency Matrix uses the element adjacency matrix
        if (type == 'dg') then
            print *, "DG Type"
            do ie_g = 1,nelem_g
                iproc = global_proc(ie_g)
                istart =  xadj_face(ie_g)
                if (ie_g < nelem_g) then
                    iend = xadj_face(ie_g + 1) - 1
                else
                    iend = xadj_face(ie_g + 1)
                end if
          
                do iface = istart,iend
             
                    je_g = adjncy_face(iface)
                    jproc = global_proc(je_g)
             
                    if (jproc /= iproc) then
                        A_proc(iproc,jproc) = 1
                    end if
                end do
          
            end do
       
           !CG Processor Adjacency Matrix uses the "generalized" adjacency matrix
        else if (type == 'cg') then
            print *, "CG Type"
            do ie_g = 1,nelem_g
          
                iproc = global_proc(ie_g)
                istart =  xadj_gen(ie_g)
                if (ie_g < nelem_g) then
                    iend = xadj_gen(ie_g + 1) - 1
                else
                    iend = xadj_gen(ie_g + 1)
                end if
          
                do iface = istart,iend
             
                    je_g = adjncy_gen(iface)
                    jproc = global_proc(je_g)
             
                    if (jproc /= iproc) then
                        A_proc(iproc,jproc) = 1
                    end if
                end do
          
            end do
        end if
    
        ! Number of adjacencies is the number of non-zero entries
        nadj_proc = sum(A_proc)
    
        ! Now construct processor adjacency matrix in CSR format
        allocate(xadj_proc(nelem_g + 1), adjncy_proc(nadj_proc))
    
        xadj_proc(1) = 1
        ied = 0
    
        do ie_g = 1,nproc
            xadj_proc(ie_g + 1) = xadj_proc(ie_g) + sum(A_proc(ie_g,:))
            do je_g = 1,nproc
                if (A_proc(ie_g,je_g) > 0) then
                    ied = ied + 1
                    adjncy_proc(ied) = je_g
                end if
            end do
        end do
        xadj_proc(nproc + 1) = xadj_proc(nproc + 1) - 1
    
        ! Destroy Full Matrix and keep the CSR format
        deallocate(A_proc)
    
    end subroutine mod_adjncy_create_adjncy_proc
  
  
    !----------------------------------------------------------------------!
    !>@brief Convert a 3D adjacency matrix to a 2D adjaceny matrix (for DD in horz. only)
    !>@details his routine assumes that the vertical dimension is structured in ascending order (level 1, 2, 3, etc.)
    !>@author James F. Kelly: 16 March 2011
    !----------------------------------------------------------------------!
    subroutine  mod_adjncy_convert_adjncy2d_face()
    
        use mod_global_grid, only: nelem_s
    
        implicit none

        integer :: ie_s, je_s, istart, iend, iface, ii
    
        !STEP 1: Determine the number of edges in the horz. adjacency graph
        nadj2d_face = 0
        do ie_s = 1, nelem_s
            istart = xadj_face(ie_s)
            if (ie_s < nelem_g) then
                iend = xadj_face(ie_s + 1) - 1
            else
                iend = xadj_face(ie_s + 1)
            end if
            do iface = istart,iend
                je_s = adjncy_face(iface)
                if (je_s <= nelem_s) then
                    nadj2d_face = nadj2d_face + 1
                end if
            end do
        end do
    
        !STEP 2: Allocate space for 2D Adjacency graphs
        allocate(adjncy2d_face(nadj2d_face), xadj2d_face(nelem_s + 1))
    
        !STEP 3: Build the 2D Adjacency Matrix
        ii = 1
        do ie_s = 1, nelem_s
            istart = xadj_face(ie_s)
            xadj2d_face(ie_s) = ii
            if (ie_s < nelem_g) then
                iend = xadj_face(ie_s + 1) - 1
            else
                iend = xadj_face(ie_s + 1)
            end if
            do iface = istart,iend
                je_s = adjncy_face(iface)
                if (je_s <= nelem_s) then
                    adjncy2d_face(ii) = je_s
                    ii = ii + 1
                end if
            end do
        end do
    
        xadj2d_face(nelem_s + 1) = nadj2d_face
    
    end subroutine  mod_adjncy_convert_adjncy2d_face
  
end module mod_adjncy
