!---------------------------------------------------------------------!
!>@brief Writes out the graph for the processor space
!>@author James F. Kelly
!---------------------------------------------------------------------!
subroutine outgraph()
  
    use mod_adjncy, only: xadj_proc, adjncy_proc
  
    use mod_parallel, only: nproc
  
    implicit none

    !local arrays
    real, dimension(:,:), allocatable:: coord_proc
    integer :: ie_g, ngl3, iproc
    integer :: istart, iend, ied, jproc
  
    allocate(coord_proc(3,nproc))
    call compute_element_center(coord_proc)
  
    print *, "Coord Proc Constructed"
  
    open(1,file='coord_proc.out')
    do iproc = 1,nproc
        write(1,*) coord_proc(1,iproc), coord_proc(2,iproc), coord_proc(3,iproc)
    end do
    close(1)
  
    open(1,file='proc_graph.out')
    do iproc = 1,nproc
        istart = xadj_proc(iproc)
        if (iproc < nproc) then
            iend = xadj_proc(iproc + 1) - 1
        else
            iend = xadj_proc(iproc + 1)
        end if
        do ied = istart,iend
            jproc = adjncy_proc(ied)
            write(1,*) iproc, jproc, 1
        end do
    end do
    close(1)
  
end subroutine outgraph

!---------------------------------------------------------------------!
!>@brief Compute element center
!>@author James F. Kelly
!---------------------------------------------------------------------!
subroutine compute_element_center(coord_proc)
  
    use mod_basis, only: nglx, ngly, nglz

    use mod_global_grid, only: coord_g_cg, intma_g, nelem_g
  
    use mod_parallel, only: nproc,nelem_l, global_proc

    implicit none

    !global arrays
    real, intent(out) :: coord_proc(3,nproc)

    !local arrays
    integer :: i, j, k, ie_g, ngl3, iproc, ip_g
  
    coord_proc = 0
    ngl3 = nglx*ngly*nglz
    do ie_g = 1,nelem_g
        iproc = global_proc(ie_g)
        do k = 1,nglz
            do j = 1,ngly
                do i = 1,nglx
                    ip_g = intma_g(i,j,k,ie_g)
                    coord_proc(:,iproc) = coord_proc(:,iproc) + coord_g_cg(:,ip_g)/ngl3
                end do
            end do
        end do
    end do
  
    coord_proc(1,:) = coord_proc(1,:)/nelem_l
    coord_proc(2,:) = coord_proc(2,:)/nelem_l
    coord_proc(3,:) = coord_proc(3,:)/nelem_l
  
end subroutine compute_element_center
