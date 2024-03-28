!---------------------------------------------------------------------!
!> @brief This subroutine writes one VTU file per group of processors 
!> using DG storage. It assembles the data into master rank in each 
!> group. It uses xml vtk rather than legacy vtk.
!> 
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine outvtk_parallel_xml(q,q_ref,fname,vorticity,pressure)
  
    use mod_basis, only: nglx, ngly, nglz, CELL_CHILDREN, nsubcells
  
    use mod_constants, only: gravity, earth_radius, pi

    use mod_grid, only: coord, npoin, nelem, intma
  
    use mod_initial, only: nvar

    use mod_input, only: format_vtk, nvtk_files, lout_vorticity, lout_rank,    &
      fname_root, locean, lLAV, lout_tau, vtk_cell_type

    use mod_utilities, only: lowercase

    use mod_parallel, only: npoin_l, npoin_l_max, nelem_l, nelem_l_max

    use mod_viscosity, only: visc_elem

    use mpi
  
    implicit none
  
    !global
    real :: q(nvar,npoin), q_ref(nvar,npoin), vorticity(3,npoin), pressure(npoin)
    character fname*100,fname_proc*200,proc_num*10
    character grid_file*100

    !local
    integer, dimension(:,:), allocatable :: cells
    real, dimension(:,:), allocatable :: q_grp
    real, dimension(:,:), allocatable :: q_ref_grp
    real, dimension(:,:), allocatable :: visc_elem_grp
    real, dimension(:,:), allocatable :: coord_grp, vorticity_grp, rank_grp, rank, pressure_grp
    integer, dimension(:,:), allocatable :: cells_grp

    real    :: x,y,z,r,lat,lon
    integer :: ncells, nchildren, nsize, nglm1, nglm13, npoly
    integer :: ivar, ifactor, ip, ic, ie, i, j, k, e, ndof, im, mp, me
    integer :: irank, ierr, nproc
    character :: form*10

    integer :: out_group, ngroups
    integer :: igroup, ig_rank, buffer
    integer :: npoin_grp, nelem_grp, ncells_grp, nproc_grp

    integer, dimension(:), allocatable  :: npoin_l_grp, nelem_l_grp
    integer :: npoin_l_max_grp, nelem_l_max_grp


    ngroups = nvtk_files

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    if(nvtk_files>nproc) ngroups=nproc

    !create MPI groups
    call set_groups(out_group,ngroups,igroup,ig_rank)
    call mpi_comm_size(out_group,nproc_grp,ierr)

    !create cells
    ncells=nelem*nsubcells
    allocate(cells(CELL_CHILDREN,ncells))
    call create_cells(cells,ncells)

    !set-up directory to store data
    if(irank==0.and.ngroups>1) call system('mkdir -p ' // trim(fname_root))
    call mpi_barrier(mpi_comm_world,ierr)
    
    if(ig_rank==0) allocate(npoin_l_grp(nproc_grp),nelem_l_grp(nproc_grp))

    !count points and elements in the group
    call mpi_gather(npoin,1,mpi_integer,npoin_l_grp,1,mpi_integer,0,out_group,ierr)
    call mpi_gather(nelem,1,mpi_integer,nelem_l_grp,1,mpi_integer,0,out_group,ierr)

    !allocate memory on each group leader rank
    if(ig_rank==0) then       

        npoin_grp = sum(npoin_l_grp)
        nelem_grp = sum(nelem_l_grp)

        npoin_l_max_grp = maxval(npoin_l_grp)
        nelem_l_max_grp = maxval(nelem_l_grp)

        ncells_grp=nelem_grp*nsubcells
        allocate(q_grp(nvar,npoin_grp))
        allocate(q_ref_grp(nvar,npoin_grp))
        if (lLAV .and. lout_tau) allocate(visc_elem_grp(nvar,npoin_grp))
        allocate(coord_grp(3,npoin_grp))
        allocate(cells_grp(CELL_CHILDREN,ncells_grp))
        if (lout_vorticity) allocate(vorticity_grp(3,npoin_grp))
        allocate(pressure_grp(1,npoin_grp))
        if (lout_rank) allocate(rank_grp(1, npoin_grp))
        coord_grp=0

    end if
    if (lout_rank) then
        allocate(rank(1,npoin))
        rank = irank
    endif

    ! !gather data within groups
    call gather_data_from(q_grp,     q,     nvar, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    call gather_data_from(q_ref_grp,     q_ref,     nvar, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    if (lLAV .and. lout_tau) call gather_data_from(visc_elem_grp, visc_elem, nvar, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    call gather_data_from(coord_grp, coord, 3,    npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    call gather_data_from(pressure_grp, pressure, 1,    npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
  
    if(lout_vorticity) call gather_data_from(vorticity_grp, vorticity, 3,    npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    if(lout_rank) call gather_data_from(rank_grp, rank, 1,    npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    
    call gather_connectivity_from(cells_grp,cells,ncells,ncells_grp,npoin_l_grp, npoin_l_max_grp, nproc_grp,&
        nelem_l_grp, nelem_l_max_grp, out_group)

    if(ig_rank==0) then
        call get_proc_num(proc_num,igroup)
        if(ngroups>1) then
           fname_proc = trim(fname)//'_p'// trim(proc_num)
           fname_proc = trim(fname_root)//'/'//trim(fname_proc) !add directory name
        else
           fname_proc = trim(fname)
        end if

        ! JEK: 10/18/18
        !
        ! To get the VTK to work optimally, we should really interpolate to an
        ! equally spaced grid since this is what the VTK LAGRANGE HEXAHEDRON
        ! elemement uses.
        !
        ! (Still looks better I think that the old way, using VTK_HEXAHEDRON
        ! elements since you do not see the DOF mesh)
        if(format_vtk=='ascii'.or.format_vtk=='ASCII') then
          if(vtk_cell_type == 'LAGRANGE') then
            call write_local_xml_ascii_lagrange(fname_proc, q_grp, npoin_grp,  &
              ncells_grp/nsubcells, nvar, coord_grp, vorticity_grp, rank_grp,  &
              pressure_grp)
          else
            call write_local_xml_ascii(fname_proc,q_grp,npoin_grp,ncells_grp,nvar,coord_grp,cells_grp, vorticity_grp, rank_grp, pressure_grp)
          endif
         else
          if(vtk_cell_type == 'LAGRANGE') then
            call write_local_xml_binary_lagrange(fname_proc, q_grp,            &
              visc_elem_grp, npoin_grp, ncells_grp/nsubcells, nvar, coord_grp, &
              vorticity_grp, rank_grp, pressure_grp)
          else
            call write_local_xml_binary(fname_proc,q_grp,q_ref_grp,visc_elem_grp,npoin_grp,ncells_grp,nvar,coord_grp,cells_grp, vorticity_grp, rank_grp, pressure_grp)
          endif
        end if

        deallocate(q_grp,coord_grp,cells_grp, pressure_grp) !not the best idea, but will work for now
        deallocate(q_ref_grp)
        if(lLAV .and. lout_tau) deallocate(visc_elem_grp)
        if(lout_vorticity) deallocate(vorticity_grp)
        if(lout_rank) deallocate(rank_grp)
    end if
    
    if(irank==0.and.ngroups>1) call write_parallel_vtu(fname,ngroups,format_vtk)

    if(ig_rank==0) deallocate(npoin_l_grp,nelem_l_grp)

end subroutine outvtk_parallel_xml

!---------------------------------------------------------------------!
!> @brief This subroutine writes one VTU file per group of processors 
!> using DG storage. It assembles the data into master rank in each 
!> group. It uses xml vtk rather than legacy vtk and prints out the grid
!> 
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine outvtk_parallel_xml_grid(fname)
  
    use mod_basis, only: nglx, ngly, nglz, nsubcells, CELL_CHILDREN
  
    use mod_constants, only: gravity, earth_radius, pi

    use mod_grid, only: coord, npoin, nelem, intma
  
    use mod_initial, only: nvar

    use mod_input, only: format_vtk, nvtk_files, vtk_cell_type

    use mod_parallel, only: npoin_l, npoin_l_max, nelem_l, nelem_l_max

    use mpi
  
    implicit none
  
    !global
    character fname*100,fname_proc*200,proc_num*10
    character grid_file*100

    !local
    integer, dimension(:,:), allocatable :: cells
    real, dimension(:,:), allocatable :: coord_grp
    integer, dimension(:,:), allocatable :: cells_grp

    real    :: x,y,z,r,lat,lon
    integer :: ncells, nsize, nglm1, nglm13, npoly
    integer :: ivar, ifactor, ip, ic, ie, i, j, k, e, ndof, im, mp, me
    integer :: irank, ierr, nproc
    character :: form*10

    integer :: out_group, ngroups
    integer :: igroup, ig_rank, buffer
    integer :: npoin_grp, nelem_grp, ncells_grp, nproc_grp

    integer, dimension(:), allocatable  :: npoin_l_grp, nelem_l_grp
    integer :: npoin_l_max_grp, nelem_l_max_grp

    ngroups = nvtk_files

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    !create MPI groups
    call set_groups(out_group,ngroups,igroup,ig_rank)
    call mpi_comm_size(out_group,nproc_grp,ierr)

    !create cells
    ncells=nelem*nsubcells
    allocate(cells(CELL_CHILDREN,ncells))
    call create_cells(cells,ncells)

    if(ig_rank==0) allocate(npoin_l_grp(nproc_grp),nelem_l_grp(nproc_grp))

    !count points and elements in the group
    call mpi_gather(npoin,1,mpi_integer,npoin_l_grp,1,mpi_integer,0,out_group,ierr)
    call mpi_gather(nelem,1,mpi_integer,nelem_l_grp,1,mpi_integer,0,out_group,ierr)

    !allocate memory on each group leader rank
    if(ig_rank==0) then

        npoin_grp = sum(npoin_l_grp)
        nelem_grp = sum(nelem_l_grp)

        npoin_l_max_grp = maxval(npoin_l_grp)
        nelem_l_max_grp = maxval(nelem_l_grp)

        ncells_grp=nelem_grp*nsubcells
        allocate(coord_grp(3,npoin_grp))
        allocate(cells_grp(CELL_CHILDREN,ncells_grp))

        coord_grp=0

    end if

    ! !gather data within groups
    call gather_data_from(coord_grp, coord, 3,    npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
  
    call gather_connectivity_from(cells_grp,cells,ncells,ncells_grp,npoin_l_grp, npoin_l_max_grp, nproc_grp,&
        nelem_l_grp, nelem_l_max_grp, out_group)
    
    if(ig_rank==0) then
        call get_proc_num(proc_num,igroup)
        fname_proc = trim(fname)//'_p'// trim(proc_num)

        if(vtk_cell_type == 'LAGRANGE') stop 'outvtk_parallel_xml_grid with LAGRANGE not implemented'
        if(format_vtk=='ascii'.or.format_vtk=='ASCII') then
            call write_local_xml_ascii_grid(fname_proc,npoin_grp,ncells_grp,coord_grp,cells_grp)
        else
            call write_local_xml_binary_grid(fname_proc,npoin_grp,ncells_grp,coord_grp,cells_grp)
        end if

        deallocate(coord_grp,cells_grp) !not the best idea, but will work for now
    end if

    if(irank==0) call write_parallel_vtu(fname,ngroups,format_vtk)

    if(ig_rank==0) deallocate(npoin_l_grp,nelem_l_grp)

end subroutine outvtk_parallel_xml_grid

subroutine outvtk_parallel_xml_mesh(fname)

  use mod_basis, only: nglx, ngly, nglz, nsubcells, CELL_CHILDREN

  use mod_constants, only: gravity, earth_radius, pi

  use mod_grid, only: coord, npoin, nelem, intma

  use mod_initial, only: nvar

  use mod_input, only: format_vtk, nvtk_files, fname_root, vtk_cell_type

  use mod_parallel, only: npoin_l, npoin_l_max, nelem_l, nelem_l_max

  use mpi

  implicit none

  !global
  character fname*100,fname_proc*200,proc_num*10

  !local
  integer, dimension(:,:), allocatable :: cells
  real, dimension(:,:), allocatable :: coord_grp
  integer, dimension(:,:), allocatable :: cells_grp

  real    :: x,y,z,r,lat,lon
  integer :: ncells, nsize, nglm1, nglm13, npoly
  integer :: ivar, ifactor, ip, ic, ie, i, j, k, e, ndof, im, mp, me
  integer :: irank, ierr, nproc
  character :: form*10

  integer :: out_group, ngroups
  integer :: igroup, ig_rank, buffer
  integer :: npoin_grp, nelem_grp, ncells_grp, nproc_grp

  integer, dimension(:), allocatable  :: npoin_l_grp, nelem_l_grp
  integer :: npoin_l_max_grp, nelem_l_max_grp
  integer :: nsubcells_back

  nsubcells_back = nsubcells
  nsubcells = 1

  ngroups = nvtk_files

  call mpi_comm_rank(mpi_comm_world,irank,ierr)
  call mpi_comm_size(mpi_comm_world,nproc,ierr)

  !create MPI groups
  call set_groups(out_group,ngroups,igroup,ig_rank)
  call mpi_comm_size(out_group,nproc_grp,ierr)

  !create cells
  ncells=nelem*nsubcells
  allocate(cells(CELL_CHILDREN,ncells))
  call create_mesh_cells(cells,ncells)

  if(ig_rank==0) allocate(npoin_l_grp(nproc_grp),nelem_l_grp(nproc_grp))

  !count points and elements in the group
  call mpi_gather(npoin,1,mpi_integer,npoin_l_grp,1,mpi_integer,0,out_group,ierr)
  call mpi_gather(nelem,1,mpi_integer,nelem_l_grp,1,mpi_integer,0,out_group,ierr)

  !allocate memory on each group leader rank
  if(ig_rank==0) then

    npoin_grp = sum(npoin_l_grp)
    nelem_grp = sum(nelem_l_grp)

    npoin_l_max_grp = maxval(npoin_l_grp)
    nelem_l_max_grp = maxval(nelem_l_grp)

    ncells_grp=nelem_grp*nsubcells
    allocate(coord_grp(3,npoin_grp))
    allocate(cells_grp(CELL_CHILDREN,ncells_grp))

    coord_grp=0

  end if

  ! !gather data within groups
  call gather_data_from(coord_grp, coord, 3,    npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)

  call gather_connectivity_from(cells_grp,cells,ncells,ncells_grp,npoin_l_grp, npoin_l_max_grp, nproc_grp,&
    nelem_l_grp, nelem_l_max_grp, out_group)

  if(ig_rank==0) then
    call get_proc_num(proc_num,igroup)
    if(ngroups>1) then
      fname_proc = trim(fname)//'_p'// trim(proc_num)
      fname_proc = trim(fname_root)//'/'//trim(fname_proc) !add directory name
    else
      fname_proc = trim(fname)
    end if

    if(vtk_cell_type == 'LAGRANGE') stop 'outvtk_parallel_xml_mesh with LAGRANGE not implemented'

    if(format_vtk=='ascii'.or.format_vtk=='ASCII') then
      call write_local_xml_ascii_grid(fname_proc,npoin_grp,ncells_grp,coord_grp,cells_grp)
    else
      call write_local_xml_binary_grid(fname_proc,npoin_grp,ncells_grp,coord_grp,cells_grp)
    end if

    deallocate(coord_grp,cells_grp) !not the best idea, but will work for now
  end if

  if(irank==0.and.ngroups>1) call write_parallel_vtu(fname,ngroups,format_vtk)

  if(ig_rank==0) deallocate(npoin_l_grp,nelem_l_grp)

  nsubcells = nsubcells_back

end subroutine outvtk_parallel_xml_mesh

!---------------------------------------------------------------------!
!> @brief This subroutine creates communicators for ngroups of 
!> processors (group_comm) and returns group number (igroup) and
!> rank within the group (ig_rank) 
!> 
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine set_groups(group_comm, ngroups, igroup, ig_rank)

    use mpi

    implicit none

    integer :: group_comm, ngroups

    integer :: irank, ierr, i
    integer :: igroup, ig_rank, group_size, nproc

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)
  
    !find desired group size
    group_size = ceiling(nproc/real(ngroups))

    if(group_size*(ngroups-1)==nproc) ngroups=ngroups-1

    !find to which group the process belongs and what is it's rank in the group
    igroup = int(irank/group_size)
    ig_rank = mod(irank,group_size)

    call mpi_comm_split (mpi_comm_world, igroup, irank, group_comm, ierr)

end subroutine set_groups

!---------------------------------------------------------------------!
!> @brief This subroutine creates cells array for VTK using INTMA 
!> 
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine create_cells(cells,ncells)

    use mod_basis, only: nglx, ngly, nglz, CELL_CHILDREN

    use mod_grid, only:  npoin, nelem, intma
  
    implicit none

    integer                      :: ncells
    integer, dimension(CELL_CHILDREN,ncells) :: cells

    integer :: ic, i,j,k,e, ii, jj, kk

    ic=0
  
    do e=1,nelem
        !construct cells in each element
        do k=1,max(nglz-1,1)
            do j=1,max(ngly-1,1)
                do i=1,max(nglx-1,1)
                    ii=min(i+1,nglx)
                    jj=min(j+1,ngly)
                    kk=min(k+1,nglz)
                    ic=ic+1
                    
                    if(nglx == 1) then
                        cells(1,ic)= intma( i, j, k,e) -1
                        cells(2,ic)= intma( i,jj, k,e) -1
                        cells(3,ic)= intma( i,jj,kk,e) -1
                        cells(4,ic)= intma( i, j,kk,e) -1
                    else if(ngly == 1) then
                        cells(1,ic)= intma( i, j, k,e) -1
                        cells(2,ic)= intma(ii, j, k,e) -1
                        cells(3,ic)= intma(ii, j,kk,e) -1
                        cells(4,ic)= intma( i, j,kk,e) -1
                    else if(nglz == 1) then
                        cells(1,ic)= intma( i, j, k,e) -1
                        cells(2,ic)= intma(ii, j, k,e) -1
                        cells(3,ic)= intma(ii,jj, k,e) -1
                        cells(4,ic)= intma( i,jj, k,e) -1
                    else
                        cells(1,ic)= intma( i, j, k,e) -1
                        cells(2,ic)= intma(ii, j, k,e) -1
                        cells(3,ic)= intma(ii,jj, k,e) -1
                        cells(4,ic)= intma( i,jj, k,e) -1
                        cells(5,ic)= intma( i, j,kk,e) -1
                        cells(6,ic)= intma(ii, j,kk,e) -1
                        cells(7,ic)= intma(ii,jj,kk,e) -1
                        cells(8,ic)= intma( i,jj,kk,e) -1
                    endif
                    
                end do !i
            end do !j
        end do
    end do !ie

  
end subroutine create_cells

subroutine create_mesh_cells(cells,ncells)

  use mod_basis, only: nglx, ngly, nglz, CELL_CHILDREN

  use mod_grid, only:  npoin, nelem, intma

  implicit none

  integer                      :: ncells
  integer, dimension(CELL_CHILDREN,ncells) :: cells

  integer :: ic, i,j,k,e, ii, jj, kk

  ic=0

  do e=1,nelem
    !construct cells in each element
    ii=nglx
    jj=ngly
    kk=nglz
    ic=ic+1

    if(nglx == 1) then
      cells(1,ic)= intma( 1, 1, 1,e) - 1
      cells(2,ic)= intma( 1,jj, 1,e) - 1
      cells(3,ic)= intma( 1,jj,kk,e) - 1
      cells(4,ic)= intma( 1, 1,kk,e) - 1
    else if(ngly == 1) then
      cells(1,ic)= intma( 1, 1, 1,e) - 1
      cells(2,ic)= intma(ii, 1, 1,e) - 1
      cells(3,ic)= intma(ii, 1,kk,e) - 1
      cells(4,ic)= intma( 1, 1,kk,e) - 1
    else if(nglz == 1) then
      cells(1,ic)= intma( 1, 1, 1,e) - 1
      cells(2,ic)= intma(ii, 1, 1,e) - 1
      cells(3,ic)= intma(ii,jj, 1,e) - 1
      cells(4,ic)= intma( 1,jj, 1,e) - 1
    else
      cells(1,ic)= intma( 1, 1, 1,e) - 1
      cells(2,ic)= intma(ii, 1, 1,e) - 1
      cells(3,ic)= intma(ii,jj, 1,e) - 1
      cells(4,ic)= intma( 1,jj, 1,e) - 1
      cells(5,ic)= intma( 1, 1,kk,e) - 1
      cells(6,ic)= intma(ii, 1,kk,e) - 1
      cells(7,ic)= intma(ii,jj,kk,e) - 1
      cells(8,ic)= intma( 1,jj,kk,e) - 1
    endif
  end do !ie

end subroutine create_mesh_cells


!---------------------------------------------------------------------!
!> @brief This subroutine creates processor file label 
!> 
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine get_proc_num(proc_num,irank)

    implicit none

    character :: proc_num*10
    integer :: irank

    if(irank<10) then
        write(proc_num,'(A4,I1)')'0000',irank
    else if(irank<100) then
        write(proc_num,'(A3,I2)')'000',irank
    else if(irank<1000) then
        write(proc_num,'(A2,I3)')'00',irank
    else if(irank<10000) then
        write(proc_num,'(A1,I3)')'0',irank
    else
        write(proc_num,'(I5)')irank
    end if
  
end subroutine get_proc_num

!---------------------------------------------------------------------!
!> @brief This subroutine writes a VTK file for each processor 
!> using DG storage rather than assembling the data into global
!> storage. It uses xml vtk rather than legacy vtk.
!> 
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine write_local_xml_ascii(fname,q,npoin,ncells,nvar,coord,cells,vorticity,rank,pressure)

    use mod_basis, only: CELL_CHILDREN
    
    use mod_input, only: lout_vorticity, icase, lout_rank, lsalinity, lincompressible
  
    implicit none

    !global arrays
    character                    :: fname*200
    real, dimension(nvar,npoin)  :: q
    real, dimension(3,npoin)     :: coord, vorticity
    real, dimension(1,npoin)     :: rank
    real, dimension(npoin)     :: pressure
    integer, dimension(CELL_CHILDREN,ncells) :: cells
    integer                      :: npoin, ncells, nvar

    !local arrays
    integer :: i
    real    :: xfactor, yfactor, zfactor

    xfactor=1.0; yfactor=1.0; zfactor=1.0
    if (icase == 1) then
        xfactor=10.0
        yfactor=10.0
    end if

    fname=trim(fname)//'.vtu'

    !Open VTK file
    open(1,file=fname)
  
    write(1,'(a)')'<?xml version="1.0"?>'
    write(1,'(a)')'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(1,'(a)')'<UnstructuredGrid>'
    write(1,'(a,I20,a,I20,a)')'<Piece NumberOfPoints="',npoin,'" NumberOfCells="',ncells,'">'
  
    !Write out point coordinates
    write(1,'(a)')'      <Points>'
    write(1,'(a)')'        <DataArray type="Float32" Name="position" NumberOfComponents="3" format="ascii">'
    do i=1,npoin
        write(1,'(10X,3E24.16)')coord(1,i)/xfactor, coord(2,i)/yfactor, coord(3,i)/zfactor
    end do
    write(1,'(a)')'        </DataArray>'
    write(1,'(a)')'      </Points>'
  
    !Write out connectivity
    write(1,'(a)')'      <Cells>'
    write(1,'(a)')'        <DataArray type="Int32" Name="connectivity" format="ascii">'
    do i=1,ncells
        write(1,'(10X,I8)') cells(:,i)
    end do
    write(1,'(a)')'        </DataArray>'
    write(1,'(a)')'        <DataArray type="Int32" Name="offsets" format="ascii">'
    do i=1,ncells
        write(1,'(10X,I20)')i*CELL_CHILDREN
    end do
    write(1,'(a)')'        </DataArray>'
    write(1,'(a)')'        <DataArray type="UInt8" Name="types" format="ascii">'
    do i=1,ncells
        write(1,'(10X,1I2)')12
    end do
    write(1,'(a)')'        </DataArray>'
    write(1,'(a)')'      </Cells>'

    !write point data
    if(lout_rank) then
       write(1,'(a)')'      <PointData Scalars="rank, density, theta" Vectors="velocity">'
    else if(lsalinity) then
       write(1,'(a)')'      <PointData Scalars="density, theta, salinity, pressure" Vectors="velocity">'
    else
       write(1,'(a)')'      <PointData Scalars="density, theta" Vectors="velocity">'
    endif

    if(lout_rank) then
        write(1,'(a)')'        <DataArray type="Float32" Name="rank" format="ascii">'
        do i=1,npoin
            write(1,'(10X,E24.16)')rank(1,i)
        end do
        write(1,'(a)')'        </DataArray>'
    endif

    write(1,'(a)')'        <DataArray type="Float32" Name="density" format="ascii">'
    do i=1,npoin
        write(1,'(10X,E24.16)')q(1,i)
    end do
    write(1,'(a)')'        </DataArray>'

    write(1,'(a)')'        <DataArray type="Float32" Name="theta" format="ascii">'
    do i=1,npoin
        write(1,'(10X,E24.16)')q(5,i)
    end do
    write(1,'(a)')'        </DataArray>'

    if(lsalinity) then
       write(1,'(a)')'        <DataArray type="Float32" Name="salinity" format="ascii">'
       do i=1,npoin
          write(1,'(10X,E24.16)')q(6,i)
       end do
       write(1,'(a)')'        </DataArray>'

       write(1,'(a)')'        <DataArray type="Float32" Name="pressure" format="ascii">'
       do i=1,npoin
          write(1,'(10X,E24.16)')pressure(i)
       end do
       write(1,'(a)')'        </DataArray>'
    end if

    write(1,'(a)')'        <DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="ascii">'
    do i=1,npoin
        write(1,'(10X,3E24.16)')q(2:4,i)
    end do

    if(lout_vorticity) then
        write(1,'(a)')'        <DataArray type="Float32" Name="vorticity" NumberOfComponents="3" format="ascii">'
        do i=1,npoin
            write(1,'(10X,3E24.16)')vorticity(1:3,i)
        end do
    end if
  
    write(1,'(a)')'        </DataArray>'
    write(1,'(a)')'       </PointData>'
    write(1,'(a)')'       <CellData>'
    write(1,'(a)')'       </CellData>'
  
    write(1,'(a)')'</Piece>'
    write(1,'(a)')'</UnstructuredGrid>'
    write(1,'(a)')'</VTKFile>'

    close(1)

end subroutine write_local_xml_ascii

subroutine lagrange_cells_3d(Nq, npts, cells)
  implicit none

  integer, intent(in) :: Nq, npts
  integer, intent(out) :: cells(npts)

  integer :: i, j, k, n, m, c

  cells(1) = 0 + 0 * Nq + 0 * Nq * Nq
  cells(2) = (Nq-1) + 0 * Nq + 0 * Nq * Nq
  cells(3) = (Nq-1) + (Nq-1) * Nq + 0 * Nq * Nq
  cells(4) = 0 + (Nq-1) * Nq + 0 * Nq * Nq
  cells(5) = 0 + 0 * Nq + (Nq-1) * Nq * Nq
  cells(6) = (Nq-1) + 0 * Nq + (Nq-1) * Nq * Nq
  cells(7) = (Nq-1) + (Nq-1) * Nq + (Nq-1) * Nq * Nq
  cells(8) = 0 + (Nq-1) * Nq + (Nq-1) * Nq * Nq

  ! edges
  c = 9
  do n = 1,Nq-2
    i = n
    j = 0
    k = 0
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo
  do n = 1,Nq-2
    i = Nq-1
    j = n
    k = 0
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo
  do n = 1,Nq-2
    i = n
    j = Nq-1
    k = 0
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo
  do n = 1,Nq-2
    i =  0
    j = n
    k = 0
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo

  do n = 1,Nq-2
    i = n
    j = 0
    k = Nq-1
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo
  do n = 1,Nq-2
    i = Nq-1
    j = n
    k = Nq-1
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo
  do n = 1,Nq-2
    i = n
    j = Nq-1
    k = Nq-1
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo
  do n = 1,Nq-2
    i =  0
    j = n
    k = Nq-1
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo

  do n = 1,Nq-2
    i = 0
    j = 0
    k = n
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo
  do n = 1,Nq-2
    i = Nq-1
    j = 0
    k = n
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo
  do n = 1,Nq-2
    i = 0
    j = Nq-1
    k = n
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo
  do n = 1,Nq-2
    i =  Nq-1
    j = Nq-1
    k = n
    cells(c) = i + j * Nq + k * Nq * Nq
    c = c + 1
  enddo

  ! faces
  do m = 1,Nq-2
    do n = 1,Nq-2
      i =  0
      j = n
      k = m
      cells(c) = i + j * Nq + k * Nq * Nq
      c = c + 1
    enddo
  enddo
  do m = 1,Nq-2
    do n = 1,Nq-2
      i =  Nq-1
      j = n
      k = m
      cells(c) = i + j * Nq + k * Nq * Nq
      c = c + 1
    enddo
  enddo
  do m = 1,Nq-2
    do n = 1,Nq-2
      i =  n
      j = 0
      k = m
      cells(c) = i + j * Nq + k * Nq * Nq
      c = c + 1
    enddo
  enddo
  do m = 1,Nq-2
    do n = 1,Nq-2
      i =  n
      j = Nq-1
      k = m
      cells(c) = i + j * Nq + k * Nq * Nq
      c = c + 1
    enddo
  enddo
  do m = 1,Nq-2
    do n = 1,Nq-2
      i =  n
      j = m
      k = 0
      cells(c) = i + j * Nq + k * Nq * Nq
      c = c + 1
    enddo
  enddo
  do m = 1,Nq-2
    do n = 1,Nq-2
      i =  n
      j = m
      k = Nq-1
      cells(c) = i + j * Nq + k * Nq * Nq
      c = c + 1
    enddo
  enddo

  ! volume
  do k = 1,Nq-2
    do j = 1,Nq-2
      do i = 1,Nq-2
        cells(c) = i + j * Nq + k * Nq * Nq
        c = c + 1
      enddo
    enddo
  enddo



end subroutine lagrange_cells_3d

subroutine write_local_xml_ascii_lagrange(fname,q,npoin,ncells,nvar,coord,vorticity,rank,pressure)

  use mod_basis, only: CELL_CHILDREN, npts, ngl, nglx, ngly, nglz

  use mod_input, only: lout_vorticity, icase, lout_rank, lsalinity

  implicit none

  !global arrays
  character                    :: fname*200
  real, dimension(nvar,npoin)  :: q
  real, dimension(3,npoin)     :: coord, vorticity
  real, dimension(1,npoin)     :: rank
  real, dimension(npoin)     :: pressure
  integer                      :: npoin, ncells, nvar

  !local arrays
  integer :: i, j
  real    :: xfactor, yfactor, zfactor
  integer :: cells(npts)

  xfactor=1.0; yfactor=1.0; zfactor=1.0
  if (icase == 1) then
    xfactor=10.0
    yfactor=10.0
  end if

  fname=trim(fname)//'.vtu'

  if (ngl .ne. nglx .and. ngl .ne. ngly .and. ngl .ne. nglz) &
    stop 'nglx == ngly == nlgz required for write_local_xml_ascii_lagrange'

  call lagrange_cells_3d(ngl, npts, cells)

  !Open VTK file
  open(1,file=fname)

  write(1,'(a)')'<?xml version="1.0"?>'
  write(1,'(a)')'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(1,'(a)')'<UnstructuredGrid>'
  write(1,'(a,I20,a,I20,a)')'<Piece NumberOfPoints="',npoin,'" NumberOfCells="',ncells,'">'

  !Write out point coordinates
  write(1,'(a)')'      <Points>'
  write(1,'(a)')'        <DataArray type="Float32" Name="position" NumberOfComponents="3" format="ascii">'
  do i=1,npoin
    write(1,'(10X,3E24.16)')coord(1,i)/xfactor, coord(2,i)/yfactor, coord(3,i)/zfactor
  end do
  write(1,'(a)')'        </DataArray>'
  write(1,'(a)')'      </Points>'

  !Write out connectivity
  write(1,'(a)')'      <Cells>'
  write(1,'(a)')'        <DataArray type="Int32" Name="connectivity" format="ascii">'
  do i=1,ncells
    write(1,'(a)', advance="no") '        '
    do j = 1,npts
      write(1,'(I8)', advance="no") cells(j) + (i-1) * npts
      if (modulo(j, 8) == 0 .and. j .ne. npts) then
        write(1,'(a)') ''
        write(1,'(a)', advance="no") '        '
      endif
    end do
    write(1,'(a)') ''
  end do
  write(1,'(a)')'        </DataArray>'
  write(1,'(a)')'        <DataArray type="Int32" Name="offsets" format="ascii">'
  do i=1,ncells
    write(1,'(10X,I20)')i*npts
  end do
  write(1,'(a)')'        </DataArray>'
  write(1,'(a)')'        <DataArray type="UInt8" Name="types" format="ascii">'
  do i=1,ncells
    write(1,'(10X,1I2)')72
  end do
  write(1,'(a)')'        </DataArray>'
  write(1,'(a)')'      </Cells>'

  !write point data
  if(lout_rank) then
    write(1,'(a)')'      <PointData Scalars="rank, density, theta" Vectors="velocity">'
  else if(lsalinity) then
    write(1,'(a)')'      <PointData Scalars="density, theta, salinity, pressure" Vectors="velocity">'
  else
    write(1,'(a)')'      <PointData Scalars="density, theta" Vectors="velocity">'
  endif

  if(lout_rank) then
    write(1,'(a)')'        <DataArray type="Float32" Name="rank" format="ascii">'
    do i=1,npoin
      write(1,'(10X,E24.16)')rank(1,i)
    end do
    write(1,'(a)')'        </DataArray>'
  endif

  write(1,'(a)')'        <DataArray type="Float32" Name="density" format="ascii">'
  do i=1,npoin
    write(1,'(10X,E24.16)')q(1,i)
  end do
  write(1,'(a)')'        </DataArray>'

  write(1,'(a)')'        <DataArray type="Float32" Name="theta" format="ascii">'
  do i=1,npoin
    write(1,'(10X,E24.16)')q(5,i)
  end do
  write(1,'(a)')'        </DataArray>'

  if(lsalinity) then
    write(1,'(a)')'        <DataArray type="Float32" Name="salinity" format="ascii">'
    do i=1,npoin
      write(1,'(10X,E24.16)')q(6,i)
    end do
    write(1,'(a)')'        </DataArray>'

    write(1,'(a)')'        <DataArray type="Float32" Name="pressure" format="ascii">'
    do i=1,npoin
      write(1,'(10X,E24.16)')pressure(i)
    end do
    write(1,'(a)')'        </DataArray>'
  end if

  write(1,'(a)')'        <DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="ascii">'
  do i=1,npoin
    write(1,'(10X,3E24.16)')q(2:4,i)
  end do

  if(lout_vorticity) then
    write(1,'(a)')'        <DataArray type="Float32" Name="vorticity" NumberOfComponents="3" format="ascii">'
    do i=1,npoin
      write(1,'(10X,3E24.16)')vorticity(1:3,i)
    end do
  end if

  write(1,'(a)')'        </DataArray>'
  write(1,'(a)')'       </PointData>'
  write(1,'(a)')'       <CellData>'
  write(1,'(a)')'       </CellData>'

  write(1,'(a)')'</Piece>'
  write(1,'(a)')'</UnstructuredGrid>'
  write(1,'(a)')'</VTKFile>'

  close(1)

end subroutine write_local_xml_ascii_lagrange

subroutine write_local_xml_ascii_grid(fname,npoin,ncells,coord,cells)

    use mod_input, only: icase
    
    use mod_basis, only: CELL_CHILDREN

    implicit none

    character                    :: fname*200
    real, dimension(3,npoin)     :: coord
    integer, dimension(CELL_CHILDREN,ncells) :: cells
    integer                      :: npoin, ncells, nvar
    integer :: i
    real :: xfactor, yfactor, zfactor

    xfactor=1.0; yfactor=1.0; zfactor=1.0
    if (icase == 1) then
        xfactor=10.0
        yfactor=10.0
    end if

    fname=trim(fname)//'.vtu'

    !Open VTK file
    open(1,file=fname)
  
    write(1,'(a)')'<?xml version="1.0"?>'
    write(1,'(a)')'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(1,'(a)')'<UnstructuredGrid>'
    write(1,'(a,I20,a,I20,a)')'<Piece NumberOfPoints="',npoin,'" NumberOfCells="',ncells,'">'
  
    !Write out point coordinates
    write(1,'(a)')'      <Points>'
    write(1,'(a)')'        <DataArray type="Float32" Name="position" NumberOfComponents="3" format="ascii">'
    do i=1,npoin
        write(1,'(10X,3E24.16)')coord(1,i)/xfactor, coord(2,i)/yfactor, coord(3,i)/zfactor
    end do
    write(1,'(a)')'        </DataArray>'
    write(1,'(a)')'      </Points>'
  
    !Write out connectivity
    write(1,'(a)')'      <Cells>'
    write(1,'(a)')'        <DataArray type="Int32" Name="connectivity" format="ascii">'
    do i=1,ncells
        write(1,'(10X,I8)') cells(:,i)
    end do
    write(1,'(a)')'        </DataArray>'
    write(1,'(a)')'        <DataArray type="Int32" Name="offsets" format="ascii">'
    do i=1,ncells
        write(1,'(10X,I20)')i*CELL_CHILDREN
    end do
    write(1,'(a)')'        </DataArray>'
    write(1,'(a)')'        <DataArray type="UInt8" Name="types" format="ascii">'
    do i=1,ncells
        write(1,'(10X,1I2)')12
    end do
    write(1,'(a)')'        </DataArray>'
    write(1,'(a)')'      </Cells>'

    write(1,'(a)')'</Piece>'
    write(1,'(a)')'</UnstructuredGrid>'
    write(1,'(a)')'</VTKFile>'

    close(1)

end subroutine write_local_xml_ascii_grid

!---------------------------------------------------------------------!
!> @brief This subroutine reads one VTU file per group of processors 
!> using DG storage. It scatters the data from master rank to processors 
!> in the group. It uses xml vtk rather than legacy vtk.
!> 
!> @author Michal A. Kopera on 07/2017
!---------------------------------------------------------------------!
subroutine outvtk_parallel_xml_read(q,q_ref,fname)
  
    use mod_basis, only: nglx, ngly, nglz, CELL_CHILDREN, nsubcells

    use mod_grid, only: coord, npoin, nelem, intma
  
    use mod_initial, only: nvar

    use mod_input, only: format_vtk, nvtk_files, lout_vorticity, lout_rank, fname_root, vtk_cell_type, is_swe_layers

    use mod_parallel, only: npoin_l, npoin_l_max, nelem_l, nelem_l_max

    use mpi
  
    implicit none
  
    !global
    real :: q(nvar,npoin), q_ref(nvar,npoin)
    character fname*100,fname_proc*200,proc_num*10
    character grid_file*100

    !local
    real, dimension(:,:), allocatable :: q_grp, q_ref_grp

    real    :: x,y,z,r,lat,lon
    integer :: ncells, nchildren, nsize, nglm1, nglm13, npoly
    integer :: ivar, ifactor, ip, ic, ie, i, j, k, e, ndof, im, mp, me
    integer :: irank, ierr, nproc
    character :: form*10

    integer :: out_group, ngroups
    integer :: igroup, ig_rank, buffer
    integer :: npoin_grp, nelem_grp, ncells_grp, nproc_grp

    integer, dimension(:), allocatable  :: npoin_l_grp, nelem_l_grp
    integer :: npoin_l_max_grp, nelem_l_max_grp

    ngroups = nvtk_files

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    !create MPI groups
    call set_groups(out_group,ngroups,igroup,ig_rank)
    call mpi_comm_size(out_group,nproc_grp,ierr)

    !create cells
    ncells=nelem*nsubcells

    if(ig_rank==0) allocate(npoin_l_grp(nproc_grp),nelem_l_grp(nproc_grp))

    !count points and elements in the group
    call mpi_gather(npoin,1,mpi_integer,npoin_l_grp,1,mpi_integer,0,out_group,ierr)
    call mpi_gather(nelem,1,mpi_integer,nelem_l_grp,1,mpi_integer,0,out_group,ierr)

    !allocate memory on each group leader rank
    if(ig_rank==0) then       


       
        npoin_grp = sum(npoin_l_grp)
        nelem_grp = sum(nelem_l_grp)

        npoin_l_max_grp = maxval(npoin_l_grp)
        nelem_l_max_grp = maxval(nelem_l_grp)

        ncells_grp=nelem_grp*nsubcells
        allocate(q_grp(nvar,npoin_grp))
        allocate(q_ref_grp(nvar,npoin_grp))

        call get_proc_num(proc_num,igroup)
        if(ngroups>1) then
           fname_proc = trim(fname)//'_p'// trim(proc_num)
           fname_proc = trim(fname_root)//'/'//trim(fname_proc)
        else
           fname_proc = trim(fname)
        end if

        ! JEK: 10/18/18
        !
        ! This routine may work fine. It doesn't seem to depend on the cells so
        ! probably does just causing code to bomb here so someone who knows more
        ! can check it.
        !
        ! NOTE: If we switch to equally spaced interpolation for LANGRANGE cells
        !       this definitely won't work....
        if(vtk_cell_type == 'LAGRANGE') &
          stop 'outvtk_parallel_xml_read with LAGRANGE not implemented (it may work, not sure... JEK)'

        if(format_vtk=='ascii'.or.format_vtk=='ASCII') then
            call read_local_xml_ascii(fname_proc,q_grp,npoin_grp,ncells_grp,nvar)
        else
            call read_local_xml_binary(fname_proc,q_grp,q_ref_grp,npoin_grp,ncells_grp,nvar)
        end if

    end if

    ! !scatter data within groups
    call scatter_data_to(q_grp,     q,     nvar, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)

    if(is_swe_layers) then
       ! !scatter data within groups
       call scatter_data_to(q_ref_grp,     q_ref,     nvar, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    end if
    
    if(ig_rank==0) then
       deallocate(q_grp, q_ref_grp) !not the best idea, but will work for now
       deallocate(npoin_l_grp,nelem_l_grp)
    end if

end subroutine outvtk_parallel_xml_read

!---------------------------------------------------------------------!
!> @brief This subroutine reads one VTU file per group of processors 
!> using DG storage. It scatters the data from master rank to processors 
!> in the group. It uses xml vtk rather than legacy vtk.
!> 
!> @author Michal A. Kopera on 07/2017
!
!> Modified by Yao Gahounzo on 10/2023
!---------------------------------------------------------------------!
subroutine outvtk_parallel_xml_read_mlswe(q,qb,qprime,fname)
  
    use mod_basis, only: nglx, ngly, nglz, CELL_CHILDREN, nsubcells

    use mod_grid, only: coord, npoin, nelem, intma
  
    use mod_initial, only: nvar

    use mod_input, only: format_vtk, nvtk_files, lout_vorticity, lout_rank, fname_root, vtk_cell_type, is_swe_layers

    use mod_parallel, only: npoin_l, npoin_l_max, nelem_l, nelem_l_max

    use mpi
  
    implicit none
  
    !global
    real :: q(3,npoin), qb(3,npoin), qprime(3,npoin)
    character fname*100,fname_proc*200,proc_num*10
    character grid_file*100

    !local
    real, dimension(:,:), allocatable :: q_grp, qb_grp, qprime_grp

    real    :: x,y,z,r,lat,lon
    integer :: ncells, nchildren, nsize, nglm1, nglm13, npoly
    integer :: ivar, ifactor, ip, ic, ie, i, j, k, e, ndof, im, mp, me
    integer :: irank, ierr, nproc
    character :: form*10

    integer :: out_group, ngroups
    integer :: igroup, ig_rank, buffer
    integer :: npoin_grp, nelem_grp, ncells_grp, nproc_grp

    integer, dimension(:), allocatable  :: npoin_l_grp, nelem_l_grp
    integer :: npoin_l_max_grp, nelem_l_max_grp

    ngroups = nvtk_files

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    !create MPI groups
    call set_groups(out_group,ngroups,igroup,ig_rank)
    call mpi_comm_size(out_group,nproc_grp,ierr)

    !create cells
    ncells=nelem*nsubcells

    if(ig_rank==0) allocate(npoin_l_grp(nproc_grp),nelem_l_grp(nproc_grp))

    !count points and elements in the group
    call mpi_gather(npoin,1,mpi_integer,npoin_l_grp,1,mpi_integer,0,out_group,ierr)
    call mpi_gather(nelem,1,mpi_integer,nelem_l_grp,1,mpi_integer,0,out_group,ierr)

    !allocate memory on each group leader rank
    if(ig_rank==0) then       


       
        npoin_grp = sum(npoin_l_grp)
        nelem_grp = sum(nelem_l_grp)

        npoin_l_max_grp = maxval(npoin_l_grp)
        nelem_l_max_grp = maxval(nelem_l_grp)

        ncells_grp=nelem_grp*nsubcells
        allocate(q_grp(3,npoin_grp))
        allocate(qb_grp(3,npoin_grp))
        allocate(qprime_grp(3,npoin_grp))

        call get_proc_num(proc_num,igroup)
        if(ngroups>1) then
           fname_proc = trim(fname)//'_p'// trim(proc_num)
           fname_proc = trim(fname_root)//'/'//trim(fname_proc)
        else
           fname_proc = trim(fname)
        end if

        ! JEK: 10/18/18
        !
        ! This routine may work fine. It doesn't seem to depend on the cells so
        ! probably does just causing code to bomb here so someone who knows more
        ! can check it.
        !
        ! NOTE: If we switch to equally spaced interpolation for LANGRANGE cells
        !       this definitely won't work....
        if(vtk_cell_type == 'LAGRANGE') &
          stop 'outvtk_parallel_xml_read with LAGRANGE not implemented (it may work, not sure... JEK)'


        call read_local_xml_ascii_mlswe(fname_proc,q_grp,qb_grp,qprime_grp,npoin_grp,ncells_grp,3)
        ! call read_local_xml_binary_mlswe(fname,q_grp,qb_grp,qprime_grp,npoin,ncells,3)

    end if

    ! !scatter data within groups
    call scatter_data_to(q_grp,q, 3, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    call scatter_data_to(qprime_grp,qprime, 3, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    call scatter_data_to(qb_grp,qb, 3, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    
    if(ig_rank==0) then
       deallocate(q_grp, qb_grp, qprime_grp) !not the best idea, but will work for now
       deallocate(npoin_l_grp,nelem_l_grp)
    end if

end subroutine outvtk_parallel_xml_read_mlswe


subroutine read_mlswe(q_df,qb_df,fname)
  
    use mod_basis, only: nglx, ngly, nglz, CELL_CHILDREN, nsubcells

    use mod_grid, only: coord, npoin, nelem, intma
  
    use mod_initial, only: nvar

    use mod_input, only: format_vtk, nvtk_files, lout_vorticity, lout_rank, fname_root, vtk_cell_type, is_swe_layers, nlayers

    use mod_parallel, only: npoin_l, npoin_l_max, nelem_l, nelem_l_max

    use mpi
  
    implicit none
  
    !global
    real :: q_df(3,npoin,nlayers), qb_df(3,npoin)
    character fname*100,fname_proc*200,proc_num*10
    character grid_file*100

    !local
    real, dimension(:,:), allocatable :: qb_grp
    real, dimension(:,:, :), allocatable :: q_grp

    real    :: x,y,z,r,lat,lon
    integer :: ncells, nchildren, nsize, nglm1, nglm13, npoly
    integer :: ivar, ifactor, ip, ic, ie, i, j, k, e, ndof, im, mp, me
    integer :: irank, ierr, nproc
    character :: form*10

    integer :: out_group, ngroups
    integer :: igroup, ig_rank, buffer
    integer :: npoin_grp, nelem_grp, ncells_grp, nproc_grp

    integer, dimension(:), allocatable  :: npoin_l_grp, nelem_l_grp
    integer :: npoin_l_max_grp, nelem_l_max_grp

    ngroups = 1

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    !create MPI groups
    call set_groups(out_group,ngroups,igroup,ig_rank)
    call mpi_comm_size(out_group,nproc_grp,ierr)

    !create cells
    ncells=nelem*nsubcells

    if(ig_rank==0) allocate(npoin_l_grp(nproc_grp),nelem_l_grp(nproc_grp))

    !count points and elements in the group
    call mpi_gather(npoin,1,mpi_integer,npoin_l_grp,1,mpi_integer,0,out_group,ierr)
    call mpi_gather(nelem,1,mpi_integer,nelem_l_grp,1,mpi_integer,0,out_group,ierr)

    !allocate memory on each group leader rank
    if(ig_rank==0) then       

        npoin_grp = sum(npoin_l_grp)
        nelem_grp = sum(nelem_l_grp)

        npoin_l_max_grp = maxval(npoin_l_grp)
        nelem_l_max_grp = maxval(nelem_l_grp)

        ncells_grp=nelem_grp*nsubcells
        allocate(q_grp(3,npoin_grp, nlayers))
        allocate(qb_grp(3,npoin_grp))

        call load_data_mlswe(fname, qb_grp, q_grp, npoin_grp, nlayers)

    end if

    ! !scatter data within groups
    call scatter_data_to(qb_grp, qb_df, 3, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)

    do k = 1, nlayers
        call scatter_data_to(q_grp(:,:,k), q_df(:,:,k), 3, npoin, npoin_grp, npoin_l_grp, npoin_l_max_grp, nproc_grp, out_group)
    end do
    
    if(ig_rank==0) then
       deallocate(q_grp, qb_grp) !not the best idea, but will work for now
       deallocate(npoin_l_grp,nelem_l_grp)
    end if

end subroutine read_mlswe

subroutine load_data_mlswe(name_data_file, qb_df_g, q_df_g, npoin_g, nlayers)

    character(len=*), intent(in) :: name_data_file

    real, dimension(3, npoin_g), intent(out) :: qb_df_g
    real, dimension(3, npoin_g, nlayers), intent(out) :: q_df_g

    real, dimension(2, npoin_g) :: coord_g
    real, dimension(npoin_g, nlayers+1) :: z_g

    integer :: count, nk, npoin_gg, i, j
    real :: dt, dt_btp

    open(unit=10, file=name_data_file, status='old', action='read')

    read(10, *) nk
    read(10, *) npoin_gg
    read(10, *) dt
    read(10, *) dt_btp

    read(10, *) (coord_g(:, i), i = 1, npoin_gg)
    
    do j = 1, 3
        read(10, *) (qb_df_g(j, i), i = 1, npoin_gg)
    end do

    read(10, *) (q_df_g(1, :, j), j = 1, nk)
    read(10, *) (q_df_g(2, :, j), j = 1, nk)
    read(10, *) (q_df_g(3, :, j), j = 1, nk)

    read(10, *) (z_g(:, j), j = 1, nk + 1)

    close(10)
end subroutine load_data_mlswe

!---------------------------------------------------------------------!
!> @brief This subroutine reads a VTK file for each processor 
!> using DG storage rather than assembling the data into global
!> storage. It uses xml vtk rather than legacy vtk.
!> 
!> @author Michal A. Kopera on 11/2016
!---------------------------------------------------------------------!
subroutine read_local_xml_ascii(fname,q,npoin,ncells,nvar)

  use mod_input, only: lsalinity
  
  implicit none

  character                    :: fname*100, dump*200
  real, dimension(nvar,npoin)  :: q
  real, dimension(3,npoin)     :: coord, vorticity
  real, dimension(npoin)     :: pressure
  integer                      :: npoin, ncells, nvar
  

  integer :: i,j


!  print*,"read_local",trim(fname)

  fname=trim(fname)//'.vtu'

  print*,"READING:",fname
  print*,"READING:",trim(fname)



  !Open VTK file
  open(1,file=fname,status='OLD',action='read')
  
  read(1,'(a)')dump
  read(1,'(a)')dump
  read(1,'(a)')dump
  read(1,'(a)')dump
  !read(1,'(a,I10,a,I10,a)')dump,npoin_read,dump,ncells_read,dump

!  write(*,'(a,I10,a,I10,a)')dump,npoin_read,dump,ncells_read,dump
!  print*,npoin, fname
  
!   if(npoin_read.ne.npoin) then
!      print*,"[READ_LOCAL_XML_ASCII] READ ERROR!"
!      stop
!   end if


  !Write out point coordinates
  read(1,'(a)')dump
  read(1,'(a)')dump
  do i=1,npoin
     read(1,'(a)')dump
  end do
  read(1,'(a)')dump
  read(1,'(a)')dump!'      </Points>'
  
  !Write out connectivity
  read(1,'(a)')dump!'      <Cells>'
  read(1,'(a)')dump!'        <DataArray type="Int32" Name="connectivity" format="ascii">'
  do i=1,ncells
     read(1,'(a)')dump! cells(:,i)
  end do
  read(1,'(a)')dump!'        </DataArray>'
  read(1,'(a)')dump!'        <DataArray type="Int32" Name="offsets" format="ascii">'
  do i=1,ncells
     read(1,'(a)')dump!i*8
  end do
  read(1,'(a)')dump!'        </DataArray>'
  read(1,'(a)')dump!'        <DataArray type="UInt8" Name="types" format="ascii">'
  do i=1,ncells
     read(1,'(a)')dump!12
  end do
  read(1,'(a)')dump!'        </DataArray>'
  read(1,'(a)')dump!'      </Cells>'

  !write point data
  read(1,'(a)')dump!'      <PointData Scalars="density, theta" Vectors="velocity">'

  read(1,'(a)')dump!'        <DataArray type="Float32" Name="density" format="ascii">'

  !write(*,'(a)')dump
  do i=1,npoin
     read(1,'(10X,E24.16)')q(1,i)
  end do
  read(1,*)!'        </DataArray>'

  read(1,*)!'        <DataArray type="Float32" Name="theta" format="ascii">'
  do i=1,npoin
     read(1,'(10X,E24.16)')q(5,i)
  end do
  read(1,*)!'        </DataArray>'

  if(lsalinity) then
     read(1,*)!'        <DataArray type="Float32" Name="theta" format="ascii">'
     do i=1,npoin
        read(1,'(10X,E24.16)')q(6,i)
     end do
     read(1,*)!'        </DataArray>'

     read(1,*)!'        <DataArray type="Float32" Name="theta" format="ascii">'
     do i=1,npoin
        read(1,'(10X,E24.16)')pressure(i)
     end do
     read(1,*)!'        </DataArray>'
  end if

  read(1,*)!'        <DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="ascii">' 
  do i=1,npoin
     read(1,'(10X,3E24.16)')q(2:4,i)
  end do

  close(1)

end subroutine read_local_xml_ascii

!---------------------------------------------------------------------!
!> @brief This subroutine reads a VTK file for each processor 
!> using DG storage rather than assembling the data into global
!> storage. It uses xml vtk rather than legacy vtk.
!> 
!> @author Michal A. Kopera on 11/2016
!---------------------------------------------------------------------!
subroutine read_local_xml_ascii_mlswe(fname,q,qb,qprime,npoin,ncells,nvar)

  use mod_input, only: lsalinity, restart_path
  
  implicit none

  character                    :: fname*100, dump*200
  real, dimension(3,npoin)  :: q, qb, qprime
  real, dimension(3,npoin)     :: coord
  integer                      :: npoin, ncells, nvar
  

  integer :: i,j
  real :: dq(2)


!  print*,"read_local",trim(fname)

  fname=trim(fname)//'.vtk'
  ! fname='mlswe_100.0000_set2c_dg_btp_bcl_explicit_P4est_l001_0002.vtk'

  print*,"READING:",fname
  ! print*,"READING:",trim(fname)



  !Open VTK file
  ! open(1,file=fname,status='OLD',action='read')
  open(1,file=fname,status='OLD',action='read')
  
  read(1,'(a)')dump
  read(1,'(a)')dump
  read(1,'(a)')dump
  read(1,'(a)')dump
  !read(1,'(a,I10,a,I10,a)')dump,npoin_read,dump,ncells_read,dump

!  write(*,'(a,I10,a,I10,a)')dump,npoin_read,dump,ncells_read,dump
!  print*,npoin, fname
  
!   if(npoin_read.ne.npoin) then
!      print*,"[READ_LOCAL_XML_ASCII] READ ERROR!"
!      stop
!   end if


  !Write out point coordinates
  read(1,'(a)')dump
  ! print*, dump
  ! read(1,'(a)')dump
  do i=1,npoin
     read(1,'(a)')dump
  end do
  read(1,'(a)')dump
  ! print*, dump
  ! read(1,'(a)')dump!'      </Points>'
  ! print*, dump
  
  !Write out connectivity
  ! read(1,'(a)')dump!'      <Cells>'
  ! print*, dump
  ! read(1,'(a)')dump!'        <DataArray type="Int32" Name="connectivity" format="ascii">'
  do i=1,ncells*5
     read(1,'(a)')dump! cells(:,i)
    !  print*, dump, ncells
  end do
  read(1,'(a)')dump!'        </DataArray>'
  ! print*, dump
  ! read(1,'(a)')dump!'        <DataArray type="Int32" Name="offsets" format="ascii">'
  do i=1,ncells
     read(1,'(a)')dump!i*8
  end do

  !write point data
  read(1,'(a)')dump!'      <PointData Scalars="dp, dprime, pb" Vectors="velocity, veloPrime, veloBaro">'
  ! print*, dump

  read(1,'(a)')dump!'        <DataArray type="Float32" Name="dp" format="ascii">'
  read(1,'(a)')dump!' 
  ! print*, dump
  do i=1,npoin
     read(1,'(E41.33)')q(1,i)  !'(10X,E24.16)'
    !  print*, q(1,i)
  end do
  !read(1,*)!'        </DataArray>'

  read(1,'(a)')dump!'        <DataArray type="Float32" Name="dpprime" format="ascii">'
  read(1,'(a)')dump!' 
  ! print*, dump
  do i=1,npoin
     read(1,'(E41.33)')qprime(1,i)
  end do
  ! read(1,*)!'        </DataArray>'

  read(1,'(a)')dump!'        <DataArray type="Float32" Name="pbpert" format="ascii">'
  ! print*, dump
  read(1,'(a)')dump!' 
  ! print*, dump
  do i=1,npoin
     read(1,'(E41.33)')qb(1,i)
    !  print*, qb(1,i)
  end do
  ! read(1,*)!'        </DataArray>'

  read(1,'(a)')dump!' !'        <DataArray type="Float32" Name="velocity" NumberOfComponents="2" format="ascii">' 
  print*, dump
  ! read(1,'(a)')dump!' 
  ! print*, dump
  do i=1,npoin
     read(1,'(2E41.33)') q(2:3,i)
    !  print*, dq, npoin
    !  print*, q(2:3,i)
  end do
  ! read(1,*)!'       </DataArray>'

  read(1,'(a)')dump!' !'        <DataArray type="Float32" Name="veloPrime" NumberOfComponents="2" format="ascii">' 
  print*, dump
  do i=1,npoin
     read(1,'(2E41.33)')qprime(2:3,i)
  end do
  ! read(1,*)!'        </DataArray>'

  read(1,'(a)')dump!' !'        <DataArray type="Float32" Name="veloBaro" NumberOfComponents="2" format="ascii">' 
  print*, dump
  do i=1,npoin
     read(1,'(2E41.33)')qb(2:3,i)
  end do

  close(1)

end subroutine read_local_xml_ascii_mlswe


!---------------------------------------------------------------------!
!> @brief This subroutine writes a binary VTK file for each processor 
!> using DG storage rather than assembling the data into global
!> storage. It uses xml vtk rather than legacy vtk.
!> 
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine write_local_xml_binary(fname,q,q_ref,visc_elem,npoin,ncells,nvar,coord,cells,vorticity,rank,pressure)

    use mod_basis, only: CELL_CHILDREN, is_2d
    
    use mod_input, only: lout_vorticity, lout_rank, lsalinity, lout_tau, lLAV, is_swe_layers
    
    implicit none

    character                    :: fname*200
    real, dimension(nvar,npoin)  :: q
    real, dimension(nvar,npoin)  :: q_ref
    real, dimension(nvar,npoin)  :: visc_elem
    real, dimension(3,npoin)     :: coord, vorticity
    real, dimension(1,npoin)     :: rank
    real, dimension(npoin)     :: pressure
    integer, dimension(CELL_CHILDREN,ncells) :: cells
    integer                      :: npoin, ncells, nvar
  
    integer                      :: fid
    character(len=500)           :: buf_string
    character                    :: lf, off_str*32
    integer                      :: i, j, offset, int, etype
    real                         :: float
  
    integer                      :: nbytes_scal, nbytes_vec, nbytes_coord, nbytes_conn, nbytes_offset, nbytes_etype

    nbytes_scal   =         npoin * sizeof(float)
    nbytes_vec    = 3  *    npoin * sizeof(float)
    nbytes_coord  = 3  *    npoin * sizeof(float)
    nbytes_conn   = CELL_CHILDREN  *   ncells * sizeof(  int)
    nbytes_offset =        ncells * sizeof(  int)
    nbytes_etype  =        ncells * sizeof(  int)
    
    if(is_2d) then
        etype=9
    else
        etype=12
    endif
    
    fid=111

    call write_xml_header_binary(fname,fid,npoin,ncells,nvar,CELL_CHILDREN)
    
    !WRITE DATA NOW
    if(lout_rank) then
        write(fid) nbytes_scal    , (rank(1,i),i=1,npoin)         ! rank
    end if
    write(fid) nbytes_scal    , (q(1,i),i=1,npoin)              ! rho

    if(is_swe_layers) then
       write(fid) nbytes_scal    , (q_ref(1,i),i=1,npoin)       ! rho_ref
    end if
    
    write(fid) nbytes_scal    , (q(5,i),i=1,npoin)              ! theta
    if(lLAV .and. lout_tau) then
       write(fid) nbytes_scal    , (visc_elem(5,i),i=1,npoin)    ! SGS
    end if
    if(lsalinity) then
       write(fid) nbytes_scal , (q(6,i),i=1,npoin)              ! salinity
       write(fid) nbytes_scal , (pressure(i),i=1,npoin)         ! pressure
    end if
    write(fid) nbytes_vec     , ((q(j,i),j=2,4),i=1,npoin)      ! velocity

    if(lout_vorticity) then
        write(fid) nbytes_vec     , ((vorticity(j,i),j=1,3),i=1,npoin)      ! vorticity
    end if
    
    write(fid) nbytes_coord   , ((coord(j,i),j=1,3),i=1,npoin)  ! coordinates
    write(fid) nbytes_conn    , ((cells(j,i),j=1,CELL_CHILDREN),i=1,ncells) ! cell connectivity
    write(fid) nbytes_offset  , (i*CELL_CHILDREN,i=1,ncells)    ! cell offsets
    write(fid) nbytes_etype   , (etype,i=1,ncells)              ! cell type

    !END WRITING DATA

    call write_xml_footer_binary(fid)

end subroutine write_local_xml_binary

subroutine write_local_xml_binary_lagrange(fname, q, visc_elem, npoin, ncells, &
    nvar, coord, vorticity, rank, pressure)

  use mod_basis, only: CELL_CHILDREN, is_2d, npts, ngl, nglx, ngly, nglz

  use mod_input, only: lout_vorticity, lout_rank, lsalinity, lout_tau, lLAV

  implicit none

  character                    :: fname*200
  real, dimension(nvar,npoin)  :: q
  real, dimension(nvar,npoin)  :: visc_elem
  real, dimension(3,npoin)     :: coord, vorticity
  real, dimension(1,npoin)     :: rank
  real, dimension(npoin)     :: pressure
  integer                      :: npoin, ncells, nvar

  integer                      :: fid
  character(len=500)           :: buf_string
  character                    :: lf, off_str*32
  integer                      :: i, j, offset, int, etype
  real                         :: float

  integer                      :: nbytes_scal, nbytes_vec, nbytes_coord, nbytes_conn, nbytes_offset, nbytes_etype

  integer :: cells(npts)

  nbytes_scal   =         npoin * sizeof(float)
  nbytes_vec    = 3  *    npoin * sizeof(float)
  nbytes_coord  = 3  *    npoin * sizeof(float)
  nbytes_conn   = npts *  ncells * sizeof(  int)
  nbytes_offset =         ncells * sizeof(  int)
  nbytes_etype  =         ncells * sizeof(  int)

  if(is_2d) then
    stop 'write_local_xml_binary_lagrange not implemented in 2d'
  else
    etype=72
  endif

  if (ngl .ne. nglx .and. ngl .ne. ngly .and. ngl .ne. nglz) &
    stop 'nglx == ngly == nlgz required for write_local_xml_ascii_lagrange'

  call lagrange_cells_3d(ngl, npts, cells)

  fid=111

  call write_xml_header_binary(fname,fid,npoin,ncells,nvar, npts)

  !WRITE DATA NOW
  if(lout_rank) then
    write(fid) nbytes_scal    , (rank(1,i),i=1,npoin)         ! rank
  end if
  write(fid) nbytes_scal    , (q(1,i),i=1,npoin)              ! rho
  write(fid) nbytes_scal    , (q(5,i),i=1,npoin)              ! theta
  if(lLAV .and. lout_tau) then
    write(fid) nbytes_scal    , (visc_elem(5,i),i=1,npoin)    ! SGS
  end if
  if(lsalinity) then
    write(fid) nbytes_scal , (q(6,i),i=1,npoin)              ! salinity
    write(fid) nbytes_scal , (q(1,i),i=1,npoin)         ! pressure
  end if
  write(fid) nbytes_vec     , ((q(j,i),j=2,4),i=1,npoin)      ! velocity

  if(lout_vorticity) then
    write(fid) nbytes_vec     , ((vorticity(j,i),j=1,3),i=1,npoin)      ! vorticity
    write(fid) nbytes_vec     , ((vorticity(j,i),j=1,3),i=1,npoin)      ! vorticity
  end if

  write(fid) nbytes_coord   , ((coord(j,i),j=1,3),i=1,npoin)  ! coordinates
  write(fid) nbytes_conn    , ((cells(j) + (i-1)*npts,j=1,npts),i=1,ncells) ! cell connectivity
  write(fid) nbytes_offset  , (npts*i,i=1,ncells)    ! cell offsets
  write(fid) nbytes_etype   , (etype,i=1,ncells)              ! cell type

  !END WRITING DATA

  call write_xml_footer_binary(fid)

end subroutine write_local_xml_binary_lagrange

subroutine write_local_xml_binary_grid(fname,npoin,ncells,coord,cells)

    use mod_basis, only: CELL_CHILDREN, is_2d
    
    implicit none

    character                    :: fname*200
    real, dimension(3,npoin)     :: coord
    integer, dimension(CELL_CHILDREN,ncells) :: cells
    integer                      :: npoin, ncells, nvar
  
    integer                      :: fid
    character(len=500)           :: buf_string
    character                    :: lf, off_str*32
    integer                      :: i, j, offset, int, etype
    real                         :: float
  
    integer                      :: nbytes_scal, nbytes_vec, nbytes_coord, nbytes_conn, nbytes_offset, nbytes_etype

    nbytes_scal   =         npoin * sizeof(float)
    nbytes_vec    = 3  *    npoin * sizeof(float)
    nbytes_coord  = 3  *    npoin * sizeof(float)
    nbytes_conn   = CELL_CHILDREN  *   ncells * sizeof(  int)
    nbytes_offset =        ncells * sizeof(  int)
    nbytes_etype  =        ncells * sizeof(  int)
  
    if(is_2d) then
        etype=9
    else
        etype=12
    endif
    
    call write_xml_header_binary(fname,fid,npoin,ncells,nvar,CELL_CHILDREN)

    !WRITE DATA NOW
    write(fid) nbytes_coord   , ((coord(j,i),j=1,3),i=1,npoin)  ! coordinates
    write(fid) nbytes_conn    , ((cells(j,i),j=1,CELL_CHILDREN),i=1,ncells) ! cell connectivity
    write(fid) nbytes_offset  , (i*CELL_CHILDREN,i=1,ncells)                ! cell offsets
    write(fid) nbytes_etype   , (etype,i=1,ncells)              ! cell type

    !END WRITING DATA

    call write_xml_footer_binary(fid)

end subroutine write_local_xml_binary_grid

!---------------------------------------------------------------------!
!> @brief This subroutine writes a VTK header for xml file
!> @author Michal A. Kopera on 02/2015
!---------------------------------------------------------------------!
subroutine write_xml_header_binary(fname,fid,npoin,ncells,nvar,cell_children)

    use mod_basis, only: is_2d
    
    use mod_input, only: lout_vorticity, lout_rank, lsalinity, lout_tau, lLAV, lincompressible, is_swe_layers

    implicit none

    character                    :: fname*200
    integer                      :: fid
    integer                      :: npoin, ncells, nvar, cell_children
  
    character(len=500)           :: buf_string
    character(len=10)            :: density_str
    character                    :: lf, off_str*32
    integer                      :: i, j, offset, int, etype
    real                         :: float
  
    integer                      :: nbytes_scal, nbytes_vec, nbytes_coord, nbytes_conn, nbytes_offset, nbytes_etype

    nbytes_scal   =         npoin * sizeof(float)
    nbytes_vec    = 3  *    npoin * sizeof(float)
    nbytes_coord  = 3  *    npoin * sizeof(float)
    nbytes_conn   = cell_children  *   ncells * sizeof(  int)
    nbytes_offset =        ncells * sizeof(  int)
    nbytes_etype  =        ncells * sizeof(  int)
  

    if(is_2d) then
        etype=9
    else
        etype=12
    endif

    lf = char(10) !end of line character

    fname=trim(fname)//'.vtu'

    !Open VTK file
#ifdef __INTEL_COMPILER
  open(fid,file=trim(fname), form='binary') !here I20 can also set little_endian if needed
#else
    open(fid,file=trim(fname), access='stream') !here I20 can also set little_endian if needed
#endif
  
    !Write file header
    buf_string = '<?xml version="1.0"?>'//lf
    write(fid)trim(buf_string)

    buf_string = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
    write(fid)trim(buf_string)

    buf_string = '  <UnstructuredGrid>'//lf
    write(fid)trim(buf_string)

    write(buf_string,'(a,I20,a,I20,a)')'    <Piece NumberOfPoints="',npoin,'" NumberOfCells="',ncells,'">'
    buf_string=trim(buf_string)//lf
    write(fid)trim(buf_string)

    !write point data header
    buf_string = '      <PointData>'//lf
    write(fid)trim(buf_string)

    offset=0
    if(lout_rank) then
        write(off_str,'(I20)') offset
        buf_string = '        <DataArray type="Float64" Name="rank" format="appended" offset="'//trim(off_str)//'"/>'//lf
        write(fid)trim(buf_string)
        offset=offset + sizeof(int) + nbytes_scal
    end if
    
    write(off_str,'(I20)') offset
    density_str = 'density'
    if(is_swe_layers) then
       density_str = "height"
    elseif(lincompressible) then
       density_str='press'
    end if
       
    buf_string = '        <DataArray type="Float64" Name="'//trim(density_str)//'" format="appended" offset="'//trim(off_str)//'"/>'//lf
    write(fid)trim(buf_string)

    !here put in the header for reference height and bathymetry also
    if(is_swe_layers) then
       offset=offset + sizeof(int) + nbytes_scal
       write(off_str,'(I20)') offset
       density_str = "h_ref"
       buf_string = '        <DataArray type="Float64" Name="'//trim(density_str)//'" format="appended" offset="'//trim(off_str)//'"/>'//lf
       write(fid)trim(buf_string)

    end if

    offset=offset + sizeof(int) + nbytes_scal
    write(off_str,'(I20)') offset
    buf_string = '        <DataArray type="Float64" Name="theta" format="appended" offset="'//trim(off_str)//'"/>'//lf
    write(fid)trim(buf_string)

    
    if(lLAV .and. lout_tau) then
       offset=offset +  sizeof(int) + nbytes_scal
       write(off_str,'(I20)') offset
       buf_string = '        <DataArray type="Float64" Name="SGS_mu" format="appended" offset="'//trim(off_str)//'"/>'//lf
       write(fid)trim(buf_string)
    end if
    
    if(lsalinity) then
       offset=offset + sizeof(int) + nbytes_scal
       write(off_str,'(I20)') offset
       buf_string = '        <DataArray type="Float64" Name="salinity" format="appended" offset="'//trim(off_str)//'"/>'//lf
       write(fid)trim(buf_string)

       offset=offset + sizeof(int) + nbytes_scal
       write(off_str,'(I20)') offset
       buf_string = '        <DataArray type="Float64" Name="pressure" format="appended" offset="'//trim(off_str)//'"/>'//lf
       write(fid)trim(buf_string)
    end if
    
    offset=offset +  sizeof(int) + nbytes_scal
    write(off_str,'(I20)') offset
    buf_string = '        <DataArray type="Float64" Name="velocity" NumberOfComponents="3" format="appended" offset="'//trim(off_str)//'"/>'//lf
    write(fid)trim(buf_string)
    
    if(lout_vorticity) then
       offset=offset +  sizeof(int) + nbytes_vec
       write(off_str,'(I20)') offset
       buf_string = '        <DataArray type="Float64" Name="vorticity" NumberOfComponents="3" format="appended" offset="'//trim(off_str)//'"/>'//lf
       write(fid)trim(buf_string)
    end if

    buf_string = '       </PointData>'//lf
    write(fid)trim(buf_string)

    buf_string = '       <CellData> </CellData>'//lf
    write(fid)trim(buf_string)
  
    !Write out point coordinates header
    buf_string = '      <Points>'//lf
    write(fid)trim(buf_string)

    offset = offset + sizeof(int) + nbytes_vec
    write(off_str,'(I20)') offset
    buf_string = '      <DataArray type="Float64" Name="position" NumberOfComponents="3" format="appended" offset="'//trim(off_str)//'"/>'//lf
    write(fid)trim(buf_string)

    buf_string = '      </Points>'//lf
    write(fid)trim(buf_string)

    !Write out connectivity header
    buf_string = '      <Cells>'//lf
    write(fid)trim(buf_string)

    offset=offset + sizeof(int) + nbytes_coord
    write(off_str,'(I20)') offset
    buf_string = '        <DataArray type="Int32" Name="connectivity" format="appended" offset="'//trim(off_str)//'"/>'//lf
    write(fid)trim(buf_string)

    offset=offset + sizeof(int) + nbytes_conn
    write(off_str,'(I20)') offset
    buf_string = '        <DataArray type="Int32" Name="offsets" format="appended" offset="'//trim(off_str)//'"/>'//lf
    write(fid)trim(buf_string)

    offset=offset + sizeof(int) + nbytes_offset
    write(off_str,'(I20)') offset
    buf_string = '        <DataArray type="Int32" Name="types" format="appended" offset="'//trim(off_str)//'"/>'//lf
    write(fid)trim(buf_string)

    buf_string = '      </Cells>'//lf
    write(fid)trim(buf_string)
  
    buf_string = '</Piece>'//lf
    write(fid)trim(buf_string)

    buf_string = '</UnstructuredGrid>'//lf
    write(fid)trim(buf_string)

    buf_string = '<AppendedData encoding="raw">'//lf
    write(fid)trim(buf_string)

    buf_string = '_'
    write(fid)trim(buf_string)

end subroutine write_xml_header_binary

!---------------------------------------------------------------------!
!> @brief This subroutine writes a VTK footer for xml file
!> @author Michal A. Kopera on 02/2015
!---------------------------------------------------------------------!
subroutine write_xml_footer_binary(fid)
  
    implicit none

    integer :: fid

    character(len=500)           :: buf_string
    character                    :: lf

    lf=char(10)
  
    buf_string = lf//'</AppendedData>'//lf
    write(fid)trim(buf_string)

    buf_string = '</VTKFile>'//lf
    write(fid)trim(buf_string)

    close(fid)

end subroutine write_xml_footer_binary

!---------------------------------------------------------------------!
!> @brief This subroutine reads a VTK header for xml file 
!> actually skips all the data there
!> @author Michal A. Kopera on 07/2017
!---------------------------------------------------------------------!
subroutine read_xml_header_binary(fname,fid)

    implicit none

    character                    :: fname*200, char12*12, char1
    integer                      :: fid
  
    character(len=500)           :: buf_string

    logical                      :: found_data=.false.

    fname=trim(fname)//'.vtu'

    !Open VTK file
#ifdef __INTEL_COMPILER
  open(fid,file=trim(fname), form='binary') !here I20 can also set little_endian if needed
#else
    open(fid,file=trim(fname), access='stream') !here I20 can also set little_endian if needed
#endif

    !reset variables
    char1 = ''
    found_data = .false.
    
    do while(.not.found_data)
       do while (char1.ne.'<')
          read(fid) char1
          !print*,char1
       end do
       read(fid) char1
       if(char1=='A') then !Looking for ApendedData - a bit rusty but works for now
          do while (char1.ne.'>')
             read(fid) char1
          end do
          !read newline
          read(fid) char1
         !print*,char1
          !read leading underscore
          read(fid) char1
          !print*,char1
          found_data=.true.          
       end if
    end do
  

end subroutine read_xml_header_binary

subroutine read_xml_header_binary_mlswe(fname,fid)

    implicit none

    character                    :: fname*100, char12*12, char1
    integer                      :: fid
  
    character(len=500)           :: buf_string

    logical                      :: found_data=.false.

    fname=trim(fname)//'.vtk'

    !Open VTK file
#ifdef __INTEL_COMPILER
  open(fid,file=trim(fname), form='binary') !here I20 can also set little_endian if needed
#else
    open(fid,file=trim(fname), access='stream') !here I20 can also set little_endian if needed
#endif

    !reset variables
    char1 = ''
    found_data = .false.
    
    do while(.not.found_data)
       do while (char1.ne.'<')
          read(fid) char1
          !print*,char1
       end do
       read(fid) char1
       if(char1=='A') then !Looking for ApendedData - a bit rusty but works for now
          do while (char1.ne.'>')
             read(fid) char1
          end do
          !read newline
          read(fid) char1
         !print*,char1
          !read leading underscore
          read(fid) char1
          !print*,char1
          found_data=.true.          
       end if
    end do
  

end subroutine read_xml_header_binary_mlswe


!---------------------------------------------------------------------!
!> @brief This subroutine reads a binary VTK file for each processor 
!> using DG storage rather than assembling the data into global
!> storage. It uses xml vtk rather than legacy vtk.
!> 
!> @author Michal A. Kopera on 07/2017
!---------------------------------------------------------------------!
subroutine read_local_xml_binary(fname,q,q_ref,npoin,ncells,nvar)

    use mod_basis, only: CELL_CHILDREN, is_2d
    
    use mod_input, only: lout_vorticity, lout_rank, lsalinity, is_swe_layers

    implicit none

    character                    :: fname*200, dump*200
    real, dimension(nvar,npoin)  :: q
    real, dimension(nvar,npoin)  :: q_ref
    real, dimension(npoin)       :: pressure
!    integer, dimension(CELL_CHILDREN,ncells) :: cells
    integer                      :: npoin, ncells, nvar
  
    integer                      :: fid
    character(len=500)           :: buf_string
    character                    :: lf, off_str*32
    integer                      :: i, j, offset, int, etype, dump_int
    real                         :: float
  
    integer                      :: nbytes_scal, nbytes_vec, nbytes_coord, nbytes_conn, nbytes_offset, nbytes_etype

    nbytes_scal   =         npoin * sizeof(float)
    nbytes_vec    = 3  *    npoin * sizeof(float)
    nbytes_coord  = 3  *    npoin * sizeof(float)
    nbytes_conn   = CELL_CHILDREN  *   ncells * sizeof(  int)
    nbytes_offset =        ncells * sizeof(  int)
    nbytes_etype  =        ncells * sizeof(  int)
  
  
    fid=111

    print*,"scalars",nbytes_scal,nbytes_vec,nbytes_coord,nbytes_offset,nbytes_etype

    print*,"Reading header",fname
    call read_xml_header_binary(fname,fid)
  
    !READ DATA NOW
    if(lout_rank) then
        read(fid) nbytes_scal    , (dump_int,i=1,npoin)         ! rank - don't need to store this
    end if
    read(fid) nbytes_scal    , (q(1,i),i=1,npoin)              ! rho

    print*,"maxvalin read",maxval(q(1,:)),npoin

    
    if(is_swe_layers) then
       read(fid) nbytes_scal, (q_ref(1,i),i=1,npoin)           ! rho_ref
    end if

    print*,"maxvalin read ref",maxval(q_ref(1,:)),npoin
    
    read(fid) nbytes_scal    , (q(5,i),i=1,npoin)              ! theta
    if(lsalinity) then
       read(fid) nbytes_scal , (q(6,i),i=1,npoin)              ! salinity
       read(fid) nbytes_scal , (pressure(i),i=1,npoin)              ! pressure
    end if
    read(fid) nbytes_vec     , ((q(j,i),j=2,4),i=1,npoin)      ! velocity

    !I don't need the rest of the data

    ! if(lout_vorticity) then
!         read(fid) nbytes_vec     , ((vorticity(j,i),j=1,3),i=1,npoin)      ! vorticity
!     end if

!     read(fid) nbytes_coord   , ((coord(j,i),j=1,3),i=1,npoin)  ! coordinates
!     read(fid) nbytes_conn    , ((cells(j,i),j=1,CELL_CHILDREN),i=1,ncells) ! cell connectivity
!     read(fid) nbytes_offset  , (i*CELL_CHILDREN,i=1,ncells)    ! cell offsets
!     read(fid) nbytes_etype   , (etype,i=1,ncells)              ! cell type

    !END READING DATA


    !call write_xml_footer_binary(fid)


    print*,"Closing ",fid,fname
    close(fid)

end subroutine read_local_xml_binary


subroutine read_local_xml_binary_mlswe(fname,q,qb,qprime,npoin,ncells,nvar)

    use mod_basis, only: CELL_CHILDREN, is_2d
    
    use mod_input, only: lout_vorticity, lout_rank, lsalinity, is_swe_layers

    implicit none

    character                    :: fname*100, dump*100
    real, dimension(3,npoin)  :: q
    real, dimension(3,npoin)  :: qb
    real, dimension(3,npoin)  :: qprime
    integer                      :: npoin, ncells, nvar
  
    integer                      :: fid
    character(len=500)           :: buf_string
    character                    :: lf, off_str*32
    integer                      :: i, j, offset, int, etype, dump_int
    real                         :: float
  
    integer                      :: nbytes_scal, nbytes_vec, nbytes_coord, nbytes_conn, nbytes_offset, nbytes_etype

    nbytes_scal   =         npoin * sizeof(float)
    nbytes_vec    = 2  *    npoin * sizeof(float)
    nbytes_coord  = 3  *    npoin * sizeof(float)
    nbytes_conn   = CELL_CHILDREN  *   ncells * sizeof(  int)
    nbytes_offset =        ncells * sizeof(  int)
    nbytes_etype  =        ncells * sizeof(  int)
  
  
    fid=111

    print*,"scalars",nbytes_scal,nbytes_vec,nbytes_coord,nbytes_offset,nbytes_etype

    print*,"Reading header",fname
    call read_xml_header_binary_mlswe(fname,fid)
  
    !READ DATA NOW
    if(lout_rank) then
        read(fid) nbytes_scal    , (dump_int,i=1,npoin)         ! rank - don't need to store this
    end if
    read(fid) nbytes_scal    , (q(1,i),i=1,npoin)              ! rho

    print*,"maxvalin read",maxval(q(1,:)),npoin

    
    read(fid) nbytes_scal, (qprime(1,i),i=1,npoin)           ! rho_ref

    print*,"maxvalin read ref",maxval(qprime(1,:)),npoin

    read(fid) nbytes_scal, (qb(1,i),i=1,npoin)           ! rho_ref

    print*,"maxvalin read ref",maxval(qb(1,:)),npoin
    
    read(fid) nbytes_vec     , ((q(j,i),j=2,3),i=1,npoin)      ! velocity
    read(fid) nbytes_vec     , ((qprime(j,i),j=2,3),i=1,npoin)      ! velocity prime
    read(fid) nbytes_vec     , ((qb(j,i),j=2,3),i=1,npoin)      ! velocity barotropic


    print*,"Closing ",fid,fname
    close(fid)

end subroutine read_local_xml_binary_mlswe



!---------------------------------------------------------------------!
!> @brief This subroutine writes a parallel VTU (.pvtu) file which puts 
!> together all VTU files written separately in parallel.
!> 
!> @author Michal A. Kopera on 01/2015
!---------------------------------------------------------------------!
subroutine write_parallel_vtu(fname,nproc,format_vtk)

    use mod_input, only: lout_vorticity, lout_rank, fname_root, lsalinity, lLAV, lout_tau

    implicit none

    character                    :: fname*100, format_vtk*8
    integer                      :: nproc

    character :: fname_proc*200, proc_num*10, type*7, form*8
    integer :: i

    !define data type parameters
    if (format_vtk=='ascii'.or.format_vtk=='ASCII') then
        type='Float32'
        form='ascii'
    else
        type='Float64'
        form='appended'
    end if
    
    fname_proc=trim(fname)//'.pvtu'

    !Open VTK file
    open(1,file=fname_proc)
  
    write(1,'(a)')'<?xml version="1.0"?>'
    write(1,'(a)')'<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(1,'(a)')'<PUnstructuredGrid GhostLevel="0">'

    write(1,'(a)')'      <PPoints>'
    write(1,'(a)')'        <PDataArray type="'//trim(type)//'" Name="position" NumberOfComponents="3" format="'//trim(form)//'"/>'
    write(1,'(a)')'      </PPoints>'

    if(lout_rank) then
        write(1,'(a)')'      <PPointData Scalars="rank, density, theta" Vectors="velocity">'
        write(1,'(a)')'        <PDataArray type="'//trim(type)//'" Name="rank" format="'//trim(form)//'"/>'
    else if (lsalinity) then
       write(1,'(a)')'      <PPointData Scalars="density, theta, salinity, pressure" Vectors="velocity">'
    else if (lLAV .and. lout_tau) then
       write(1,'(a)')'      <PPointData Scalars="density, theta, SGS_mu" Vectors="velocity">'
    else
        write(1,'(a)')'      <PPointData Scalars="density, theta" Vectors="velocity">'
    endif
    write(1,'(a)')'        <PDataArray type="'//trim(type)//'" Name="density" format="'//trim(form)//'"/>'
    write(1,'(a)')'        <PDataArray type="'//trim(type)//'" Name="theta" format="'//trim(form)//'"/>'
    if(lsalinity) then
       write(1,'(a)')'        <PDataArray type="'//trim(type)//'" Name="salinity" format="'//trim(form)//'"/>'
       write(1,'(a)')'        <PDataArray type="'//trim(type)//'" Name="pressure" format="'//trim(form)//'"/>'
    end if
    
    if(lLAV .and. lout_tau) then
       write(1,'(a)')'        <PDataArray type="'//trim(type)//'" Name="SGS_mu" format="'//trim(form)//'"/>'
    end if
    
    write(1,'(a)')'        <PDataArray type="'//trim(type)//'" Name="velocity" format="'//trim(form)//'" NumberOfComponents="3"/>'

    if(lout_vorticity)   write(1,'(a)')'        <PDataArray type="'//trim(type)//'" Name="vorticity" format="'//trim(form)//'" NumberOfComponents="3"/>'

    write(1,'(a)')'      </PPointData>'
    write(1,'(a)')'      <PCellData>'
    write(1,'(a)')'      </PCellData>'

    do i=1,nproc
        call get_proc_num(proc_num,i-1)
        fname_proc = trim(fname_root)//'/'//trim(fname)//'_p'// trim(proc_num)//'.vtu'
        write(1,'(a,a,a)')'      <Piece Source="',trim(fname_proc),'"/>'
    end do

    write(1,'(a)')'</PUnstructuredGrid>'
    write(1,'(a)')'</VTKFile>'

    close(1)

end subroutine write_parallel_vtu


!=====================================================================
! UNUSED
!=====================================================================


!---------------------------------------------------------------------!
!> @brief This subroutine writes one VTU file per group of processors 
!> using DG storage. It uses double buffer to overlap assembly of
!> the data and writing to disk.
!> Still a prototype and not tested well
!> 
!> @author Michal A. Kopera on 02/2015
!---------------------------------------------------------------------!
subroutine outvtk_parallel_db_xml(q,fname)
  
    use mod_basis, only: nglx, ngly, nglz, nsubcells, CELL_CHILDREN, is_2d
  
    use mod_constants, only: gravity, earth_radius, pi

    use mod_grid, only: coord, npoin, nelem, intma
  
    use mod_initial, only: nvar

    use mod_input, only: format_vtk, nvtk_files

    use mod_parallel, only: npoin_l, npoin_l_max, nelem_l, nelem_l_max

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION
  
    implicit none
  
    !global
    real :: q(nvar,npoin)
    character fname*200,fname_proc*200,proc_num*10
 
    !local
    integer, dimension(:,:), allocatable :: cells
    real, dimension(:,:), allocatable :: rscal_buff
    real, dimension(:,:,:), allocatable :: rvec_buff
    integer, dimension(:,:,:), allocatable :: ivec_buff

    real    :: x,y,z,r,lat,lon
    integer :: ncells, nsize, nglm1, nglm13, npoly
    integer :: ivar, ifactor, ip, ic, ie, i, j, k, e, ndof, im, mp, me
    integer :: irank, ierr, nproc
    character :: form*10

    integer :: out_group, ngroups, etype, offset
    integer :: igroup, ig_rank, buffer
    integer :: npoin_grp, nelem_grp, ncells_grp, nproc_grp

    integer, dimension(:), allocatable  :: npoin_l_grp, nelem_l_grp
    integer :: npoin_l_max_grp, nelem_l_max_grp
    integer :: request, tag, status(MPI_STATUS_SIZE)

    integer :: nbytes_scal, nbytes_vec, nbytes_coord, nbytes_conn, nbytes_offset, nbytes_etype

    integer :: fid, ibuf_write, ibuf_recv, npoin_recv, npoin_write

    real :: float
    integer :: int

    ngroups = nvtk_files
    if(is_2d) then
        etype=9
    else
        etype=12
    endif

    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    !create MPI groups
    call set_groups(out_group,ngroups,igroup,ig_rank)
    call mpi_comm_size(out_group,nproc_grp,ierr)

    !create cells
    ncells=nelem*nsubcells
    allocate(cells(CELL_CHILDREN,ncells))
    call create_cells(cells,ncells)

    if(ig_rank==0) allocate(npoin_l_grp(nproc_grp),nelem_l_grp(nproc_grp))

    !count points and elements in the group
    call mpi_gather(npoin,1,mpi_integer,npoin_l_grp,1,mpi_integer,0,out_group,ierr)
    call mpi_gather(nelem,1,mpi_integer,nelem_l_grp,1,mpi_integer,0,out_group,ierr)

    !allocate memory on each group leader rank
    if(ig_rank==0) then

        npoin_grp = sum(npoin_l_grp)
        nelem_grp = sum(nelem_l_grp)

        npoin_l_max_grp = maxval(npoin_l_grp)
        nelem_l_max_grp = maxval(nelem_l_grp)

        ncells_grp=nelem_grp*nsubcells

        nbytes_scal   =         npoin_grp * sizeof(float)
        nbytes_vec    = 3  *    npoin_grp * sizeof(float)
        nbytes_coord  = 3  *    npoin_grp * sizeof(float)
        nbytes_conn   = CELL_CHILDREN  *   ncells_grp * sizeof(  int)
        nbytes_offset =        ncells_grp * sizeof(  int)
        nbytes_etype  =        ncells_grp * sizeof(  int)


        allocate(rscal_buff(npoin_l_max_grp,2)) !allocate double buffer
        allocate(rvec_buff(3,npoin_l_max_grp,2)) !allocate double buffer
        allocate(ivec_buff(CELL_CHILDREN,nelem_l_max_grp*nsubcells,2)) !allocate double buffer

        call get_proc_num(proc_num,igroup)
        fname_proc = trim(fname)//'_p'// trim(proc_num)

        !open file and write file header
        call write_xml_header_binary(fname_proc,fid,npoin_grp,ncells_grp,nvar,CELL_CHILDREN)

    end if

    request=mpi_request_null
    tag=111

    !! BEGIN WRITING DENSITY

    !initialize first buffer with rank 0 data
    if (ig_rank==0) then
        rscal_buff(:,1)=q(1,:)
        write(fid) nbytes_scal
    end if

    !go over all processes in the group, getting data while writing the other buffer
    do i=1,nproc_grp-1

        ibuf_write = mod(i+1,2)+1
        ibuf_recv = mod(i,2)+1

        ! process i sends data to group leader
        if (ig_rank==i) call mpi_isend(q(1,:),npoin,MPI_PRECISION,0,tag,out_group,request,ierr)
     
        if(ig_rank==0) then
            ! group leader receives to ibuf_recv

            npoin_recv  = npoin_l_grp(i+1)
            npoin_write = npoin_l_grp(i)

            call mpi_irecv(rscal_buff(:,ibuf_recv),npoin_recv,MPI_PRECISION,i,tag,out_group,request,ierr)

            !but in the meantime writes previously stored data in ibuf_write
            write(fid) (rscal_buff(j,ibuf_write), j=1,npoin_write)
        end if

        !wait till communication done before sending/receiving another batch
        call mpi_wait(request,status,ierr)

    end do

    !final write of last received buffer
    if(ig_rank==0) write(fid) (rscal_buff(j,ibuf_recv), j=1,npoin_recv)

    !! END WRITING DENSITY


    !! BEGIN WRITING THETA

    !initialize first buffer with rank 0 data
    if (ig_rank==0) then
        rscal_buff(:,1)=q(5,:)
        write(fid) nbytes_scal
    end if

    !go over all processes in the group, getting data while writing the other buffer
    do i=1,nproc_grp-1

        ibuf_write = mod(i+1,2)+1
        ibuf_recv = mod(i,2)+1

        ! process i sends data to group leader
        if (ig_rank==i) call mpi_isend(q(5,:),npoin,MPI_PRECISION,0,tag,out_group,request,ierr)
     
        if(ig_rank==0) then
            ! group leader receives to ibuf_recv

            npoin_recv  = npoin_l_grp(i+1)
            npoin_write = npoin_l_grp(i)

            call mpi_irecv(rscal_buff(:,ibuf_recv),npoin_recv,MPI_PRECISION,i,tag,out_group,request,ierr)

            !but in the meantime writes previously stored data in ibuf_write
            write(fid) (rscal_buff(j,ibuf_write), j=1,npoin_write)
        end if

        !wait till communication done before sending/receiving another batch
        call mpi_wait(request,status,ierr)

    end do

    !final write of last received buffer
    if(ig_rank==0) write(fid) (rscal_buff(j,ibuf_recv), j=1,npoin_recv)

    !! END WRITING THETA

    !! BEGIN WRITING VELOCITY

    !initialize first buffer with rank 0 data
    if (ig_rank==0) then
        rvec_buff(:,:,1)=q(2:4,:)
        write(fid) nbytes_vec
    end if

    !go over all processes in the group, getting data while writing the other buffer
    do i=1,nproc_grp-1

        ibuf_write = mod(i+1,2)+1
        ibuf_recv = mod(i,2)+1

        ! process i sends data to group leader
        if (ig_rank==i) call mpi_isend(q(2:4,:),3*npoin,MPI_PRECISION,0,tag,out_group,request,ierr)
     
        if(ig_rank==0) then
            ! group leader receives to ibuf_recv

            npoin_recv  = npoin_l_grp(i+1)
            npoin_write = npoin_l_grp(i)

            call mpi_irecv(rvec_buff(:,:,ibuf_recv),3*npoin_recv,MPI_PRECISION,i,tag,out_group,request,ierr)

            !but in the meantime writes previously stored data in ibuf_write
            write(fid) ((rvec_buff(k,j,ibuf_write), k=1,3), j=1,npoin_write)
        end if

        !wait till communication done before sending/receiving another batch
        call mpi_wait(request,status,ierr)

    end do

    !final write of last received buffer
    if(ig_rank==0) write(fid) ((rvec_buff(k,j,ibuf_recv), k=1,3), j=1,npoin_recv)

    !! END WRITING VELOCITY


    !! BEGIN WRITING COORD

    !initialize first buffer with rank 0 data
    if (ig_rank==0) then
        rvec_buff(:,:,1)=coord
        write(fid) nbytes_coord
    end if

    !go over all processes in the group, getting data while writing the other buffer
    do i=1,nproc_grp-1

        ibuf_write = mod(i+1,2)+1
        ibuf_recv = mod(i,2)+1

        ! process i sends data to group leader
        if (ig_rank==i) call mpi_isend(coord,3*npoin,MPI_PRECISION,0,tag,out_group,request,ierr)
     
        if(ig_rank==0) then
            ! group leader receives to ibuf_recv

            npoin_recv  = npoin_l_grp(i+1)
            npoin_write = npoin_l_grp(i)

            call mpi_irecv(rvec_buff(:,:,ibuf_recv),3*npoin_recv,MPI_PRECISION,i,tag,out_group,request,ierr)

            !but in the meantime writes previously stored data in ibuf_write
            write(fid) ((rvec_buff(k,j,ibuf_write), k=1,3), j=1,npoin_write)
        end if

        !wait till communication done before sending/receiving another batch
        call mpi_wait(request,status,ierr)

    end do

    !final write of last received buffer
    if(ig_rank==0) write(fid) ((rvec_buff(k,j,ibuf_recv), k=1,3), j=1,npoin_recv)

    !! END WRITING COORD

    !! BEGIN WRITING CELLS

    !initialize first buffer with rank 0 data
    if (ig_rank==0) then
        ivec_buff(:,:,1)=cells
        write(fid) nbytes_conn
        offset=0
    end if

    !go over all processes in the group, getting data while writing the other buffer
    do i=1,nproc_grp-1

        ibuf_write = mod(i+1,2)+1
        ibuf_recv = mod(i,2)+1

        ! process i sends data to group leader
        if (ig_rank==i) call mpi_isend(cells,CELL_CHILDREN*ncells,mpi_integer,0,tag,out_group,request,ierr)
     
        if(ig_rank==0) then
            ! group leader receives to ibuf_recv

            npoin_recv  = nelem_l_grp(i+1) *nsubcells
            npoin_write = nelem_l_grp(i)   *nsubcells

            call mpi_irecv(ivec_buff(:,:,ibuf_recv),CELL_CHILDREN*npoin_recv,mpi_integer,i,tag,out_group,request,ierr)

            !but in the meantime writes previously stored data in ibuf_write
            write(fid) ((ivec_buff(k,j,ibuf_write) + offset, k=1,CELL_CHILDREN), j=1,npoin_write)
            !make sure cells are offset correctly
            offset = offset+npoin_l_grp(i)
        end if

        !wait till communication done before sending/receiving another batch
        call mpi_wait(request,status,ierr)

    end do

    !final write of last received buffer
    if(ig_rank==0) write(fid) ((ivec_buff(k,j,ibuf_recv) + offset, k=1,CELL_CHILDREN), j=1,npoin_recv)

    !! END WRITING CELLS

  
    if(ig_rank==0) then
        !! WRITE OFFSETS AND TYPES NORMALLY
        write(fid) nbytes_offset  , (i*CELL_CHILDREN,i=1,ncells_grp)                ! cell offsets
        write(fid) nbytes_etype   , (etype,i=1,ncells_grp)              ! cell type
     
        call write_xml_footer_binary(fid)
        deallocate(rscal_buff,rvec_buff,ivec_buff)
        deallocate(npoin_l_grp,nelem_l_grp)
    end if

    if(irank==0) call write_parallel_vtu(fname,ngroups,format_vtk)

end subroutine outvtk_parallel_db_xml
