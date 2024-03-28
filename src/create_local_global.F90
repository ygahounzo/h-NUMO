!----------------------------------------------------------------------!
!>@brief This subroutine creates local_global GridPoint-based mappings using 
!>global_proc which is generated via METIS
!>@author James F. Kelly 24 June 2010
!>@date F.X. Giraldo January 2014 to handle periodicity and to 
!>make INTENT IN and OUT explicit
!----------------------------------------------------------------------!
subroutine create_local_global_poin(local_global_poin_l,local_global_poin_l_cg, local_global_poin_periodic_l,local_global_elem_l,  &
    flag_periodic_l,npoin_l_max,npoin_l,npoin_l_cg, nelem_l,nelem_l_max,nproc)
  
    use mod_basis, only: nglx, ngly, nglz
  
    use mod_global_grid, only: npoin_g_cg, nelem_g, intma_g

    use mod_input, only: decomp_type, space_method
  
    implicit none
  
    !global arrays
    integer, intent(in) :: nproc, npoin_l_max, nelem_l_max, nelem_l(nproc)
    integer, dimension(npoin_l_max,nproc), intent(out) :: local_global_poin_l_cg, local_global_poin_l, local_global_poin_periodic_l, flag_periodic_l
    integer, intent(out) :: npoin_l(nproc), npoin_l_cg(nproc)
    integer, intent(in)  :: local_global_elem_l(nelem_l_max,nproc)

    !local arrays
    integer :: ie, i, j, k, npoin_cg, isnew, ip_g, ipp, ie_g, iproc, i2
    integer :: n_reshape, ngl_xyz, ngl_xy, ip, ipcg, i1, j1, k1, ie1, ii, AllocateStatus
    integer, dimension(:), allocatable :: intma_reshape, order_reshape

    !Initialize
    local_global_poin_l= 0
    local_global_poin_periodic_l= 0
    flag_periodic_l=0
    npoin_l = 0

    !Define some Constants: number of LGL Points per dimension
    ngl_xy  = nglx*ngly
    ngl_xyz = ngl_xy*nglz
  


    do iproc = 1,nproc
        ! Reshpae INTMA into a vector
        n_reshape = ngl_xyz*nelem_l(iproc)
        allocate(intma_reshape(n_reshape), order_reshape(n_reshape), &
            stat=AllocateStatus)
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Create_Local_Global_Poin 2 **"
     
        !Unravel INTMA_G into a vector
        do ie = 1,nelem_l(iproc)
            ie_g = local_global_elem_l(ie,iproc)
            do k = 1,nglz
                do j = 1,ngly
                    do i = 1,nglx
                        ii = i + (j - 1)*nglx + (k - 1)*ngl_xy + (ie - 1)*ngl_xyz
                        intma_reshape(ii) = intma_g(i,j,k,ie_g)
                    end do
                end do
            end do
        end do
     
        ! Sort INTMA_RESHAPE
        call quick_sort(intma_reshape,order_reshape,n_reshape)
     
        ip = 1
        ipcg = 1
        local_global_poin_l(ip,iproc) = intma_reshape(1)
        local_global_poin_l_cg(ipcg,iproc) = intma_reshape(1)
        local_global_poin_periodic_l(ip,iproc) = intma_reshape(1)
        flag_periodic_l(ip,iproc)=intma_reshape(1)
        do ii = 2,n_reshape
            if (space_method=='cgc') then
                if (intma_reshape(ii) /= intma_reshape(ii - 1)) then
                    ip = ip + 1
                    ipcg = ipcg + 1
                    local_global_poin_l(ip,iproc) = intma_reshape(ii)
                    local_global_poin_l_cg(ipcg,iproc) = intma_reshape(ii)
                    local_global_poin_periodic_l(ip,iproc) = intma_reshape(ii)
                    flag_periodic_l(ip,iproc)=intma_reshape(ii)
                end if
            else
                ip = ip + 1
                local_global_poin_l(ip,iproc) = intma_reshape(ii)
                if (intma_reshape(ii) /= intma_reshape(ii - 1)) then
                    ipcg = ipcg + 1
                    local_global_poin_l_cg(ipcg,iproc) = intma_reshape(ii)
                endif
                local_global_poin_periodic_l(ip,iproc) =  intma_reshape(ii)
                flag_periodic_l(ip,iproc)= intma_reshape(ii)
            end if
        end do
        npoin_l(iproc) = ip
        npoin_l_cg(iproc) = ipcg

        !Deallocate Arrays
        deallocate(intma_reshape, order_reshape)
     
    end do
  
end subroutine create_local_global_poin

!----------------------------------------------------------------------!
!>@brief This subroutine creates local_global element-based mappings usign global_proc which is
!> generated via METIS
!>@author James F. Kelly 24 June 2010
!----------------------------------------------------------------------!
subroutine create_local_global_elem(local_global_elem_l,nelem_l_max, &
    nproc,nelem_l,global_proc)
  
    use mod_global_grid, only: nelem_g
  
    implicit none

    !global arrays
    integer, intent(in)  :: nelem_l_max, nproc
    integer, intent(out) :: local_global_elem_l(nelem_l_max,nproc)
    integer, intent(in)  :: nelem_l(nproc), global_proc(nelem_g)

    !local arrays
    integer :: iproc, ie_l, ie_g
    integer :: iel(nproc)
  
    local_global_elem_l = 0
    iel = 0
  
    ! STEP 1: Create the LOCAL/GLOBAL MAPPING for ELEMENTS
    ! Loop over global elements
    do ie_g = 1,nelem_g
        ! Get the processor that element ie_g belongs to
        iproc = global_proc(ie_g)
     
        !Increment the local element number on processor iproc
        iel(iproc) = iel(iproc) + 1
        ie_l = iel(iproc)
        local_global_elem_l(ie_l,iproc) = ie_g
    end do
  
end subroutine create_local_global_elem

!-----------------------------------------------------------------!
!>@brief Creates Local/Global mappings for the face data strucutres
!>@details face_type = 1, 2, 3 or 4
!> 1 = internal face
!> 2 = inter-processor face, left element is on-processor, right-element is off-processor
!> 3 = inter-processor face, left element is off-processor, right-element is on-processor
!> 4 = physical boundary face
!>@author James F. Kelly 6 October 2010
!-----------------------------------------------------------------!
subroutine create_local_global_face(local_global_face_l,nface_l,face_type_l, &
    nface_l_max,global_proc,nproc)
  
    use mod_global_grid, only: intma_g, face_g, nface_g, nelem_g
  
    implicit none
  
    integer nproc, nface_l_max
    integer global_proc(nelem_g)
    integer nface_l(nproc)
    integer local_global_face_l(nface_l_max,nproc)
    integer face_type_l(nface_l_max,nproc)
    integer iproc, ier, iel, jproc
    integer iface, ifacel, ifacer
    integer il, jl, kl, ipl
    integer ir, jr, kr, ipr
    integer ib, ibb, itype
  
    nface_l = 0
  
    ! Loop over all the faces
    do iface = 1,nface_g
        iel = face_g(7,iface)
        ier = face_g(8,iface)
        iproc = global_proc(iel)
        if (ier > 0 ) then
            jproc = global_proc(ier)
        end if
     
        ! Internal Side
        if (ier > 0 .and. iproc == jproc) then
            nface_l(iproc) =  nface_l(iproc) + 1
            ifacel = nface_l(iproc)
            local_global_face_l(ifacel,iproc) = iface
            face_type_l(ifacel,iproc) = 1
        
        ! Inter-Processor Boundary
        else if (ier >0 .and. iproc /= jproc) then
            nface_l(iproc) =  nface_l(iproc) + 1
            ifacel = nface_l(iproc)
            nface_l(jproc) =  nface_l(jproc) + 1
            ifacer = nface_l(jproc)
            local_global_face_l(ifacel,iproc) = iface
            local_global_face_l(ifacer,jproc) = iface
            face_type_l(ifacel,iproc) = 2
            face_type_l(ifacer,jproc) = 3
        
        ! Physical Boundary
        else if (ier < 0) then
            nface_l(iproc) =  nface_l(iproc) + 1
            ifacel = nface_l(iproc)
            local_global_face_l(ifacel,iproc) = iface
            face_type_l(ifacel,iproc) = 4
        end if
    end do

end subroutine create_local_global_face
