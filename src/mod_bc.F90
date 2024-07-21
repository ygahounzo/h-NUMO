!---------------------------------------------------------------------!
!>@brief This module contains Boundary Condition information (excluding iperiodic,
!>  which is contained in mod_grid)
!>@author James F. Kelly 29 November 2009
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943-5216
!>
!>@date 29 August 2011 F.X. Giraldo
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943-5216
!----------------------------------------------------------------------!
module mod_bc

    use mod_basis, only: ngl, nglx, ngly, nglz

    use mod_global_grid, only: xperiodic, yperiodic, zperiodic, iperiodic_g

    use mod_grid, only: npoin, nbsido, coord, nelem

    use mod_input, only: icase, sponge_type

    use mod_types, only : r8

    public :: &
        mod_bc_create, &
        mod_bc_create_iperiodic, &
        read_bc, &
        create_bc_list,&
        normals, &
        bc_list, bc_count, & !for identifying boundary points not on boundary element faces (unstructured meshes)
        vc_el_type
    

    private
    !-----------------------------------------------------------------------
    real,    dimension(:,:),    allocatable :: normals
    integer, dimension(:),      allocatable :: ipoin_bound, ip_bound_bc
    integer, dimension(:,:),    allocatable :: bc_list
    integer, dimension(4) :: bc_count
    integer, dimension(:),      allocatable :: vc_el_type
    integer npoin_bound
    integer npoin_errmask, num_bound_bc
    logical ldirichlet
  !-----------------------------------------------------------------------  

contains

    !-----------------------------------------------------------------------
    subroutine mod_bc_create()

        implicit none

        integer AllocateStatus
        integer ierr
        
        ! Allocate Memory for BC data structures

        allocate(normals(3,npoin), ipoin_bound(npoin), stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_BC **"

        !Construct PMATRIX for NFBCs
        call create_nfbc_vector(normals,ipoin_bound,npoin_bound,ldirichlet)

        call create_bc_list()
        
    end subroutine mod_bc_create

    !-----------------------------------------------------------------------!
    !@brief This subroutine reads boundary condition patches. First, read bc.inp,
    !       which specifies the number of boundary condition files to read,
    !       followed by the file names themselves and their associated boundary
    !       condition type. Each boundary condition file has two header lines. The
    !       third line specifies the number of points in each coordinate
    !       direction.
    !
    !       :bc.inp - format:
    !       4
    !       "bc1.dat" 4
    !       "bc2.dat" 4
    !       "bc3.dat" 3
    !       "bc4.dat" 3
    !
    !       :boundary condition file - format. Ex. "bc1.dat":
    !       junk header
    !       junk header
    !       2 2
    !       0.0 0.0 0.0
    !       0.0 1.0 0.0
    !       0.0 0.0 1.0
    !       0.0 1.0 1.0
    !
    !@author Nathan A. Wukie on 03/2015
    !         University of Cincinnati
    !         Department of Aerospace Engineering
    !         Cincinnati, OH 45220
    !-----------------------------------------------------------------------!
    subroutine read_bc(face,nface,bsido,nbsido)
        use mod_types,    only: r8

        use mod_basis, only: FACE_LEN

        use mod_grid,     only: coord_cg

        implicit none

        integer                       :: nface,nbsido,npoin_cg,iboun,AllocateStatus
        integer                       :: nfiles,nptsi,nptsj,npts,bc_type,num_match
        integer                       :: ifiles,ipoint_f,ipoint_bc,ip,iface
        integer, dimension(FACE_LEN,nface)   :: face
        integer, dimension(6,nbsido)  :: bsido
        character(len=100)            :: bc_file
        real(kind=r8)                 :: tol,face_x,face_y,face_z
        logical                       :: lmatch
        real(kind=r8), dimension(:), allocatable :: bc_x,bc_y,bc_z

        ! set boundary condition match tolerance
        tol = 1.e-5

        ! open boundary condition input file
        open(1,file="bc.inp")
        read(1,*) nfiles  !number of boundary condition files to read

        do ifiles = 1,nfiles

            ! get new boundary condition file name and associated type
            read(1,*) bc_file, bc_type

            ! open new boundary condition file
            open(2,file=bc_file)

            ! discard first two lines. they are just extra information from the
            ! grid generator
            read(2,*)
            read(2,*)

            ! get number of points in each direction, (i,j) - structured
            read(2,*) nptsi,nptsj
            npts = nptsi*nptsj

            allocate(bc_x(npts),bc_y(npts),bc_z(npts),stat=AllocateStatus)
            if (AllocateStatus /= 0) stop "Memory error - read_bc"
         
            ! read boundary condition patch coordinates
            do ipoint_bc = 1,npts
                read(2,*) bc_x(ipoint_bc), bc_y(ipoint_bc), bc_z(ipoint_bc)
            end do

            ! loop through grid boundary faces
            do iboun = 1,nbsido

                ! loop through the points of the current face
                num_match = 0
                do ipoint_f = 1,4
                    ip = bsido(ipoint_f,iboun)
                    face_x = coord_cg(1,ip)
                    face_y = coord_cg(2,ip)
                    face_z = coord_cg(3,ip)

                    do ipoint_bc = 1,npts
                        ! check if current face point is included in the boundary
                        ! condition patch
                        lmatch = (abs(face_x - bc_x(ipoint_bc)) < tol) .and. &
                            (abs(face_y - bc_y(ipoint_bc)) < tol) .and. &
                            (abs(face_z - bc_z(ipoint_bc)) < tol)

                        if (lmatch) then
                            num_match = num_match + 1
                            exit
                        end if

                    end do !ipoint_bc
                end do ! ipoint_f

                ! if all points are included in the current boundary condition
                ! patch, then set boundary condition type
                if (num_match == 4) then
                    bsido(6,iboun) = bc_type
                end if

            end do !iboun


            do iface = 1,nface

                ! loop through the points of the current face
                num_match = 0
                do ipoint_f = 1,4
                    ip = face(ipoint_f,iface)
                    face_x = coord_cg(1,ip)
                    face_y = coord_cg(2,ip)
                    face_z = coord_cg(3,ip)

                    do ipoint_bc = 1,npts
                        ! check if current face point is included in the boundary
                        ! condition patch
                        lmatch = (abs(face_x - bc_x(ipoint_bc)) < tol) .and. &
                            (abs(face_y - bc_y(ipoint_bc)) < tol) .and. &
                            (abs(face_z - bc_z(ipoint_bc)) < tol)

                        if (lmatch) then
                            num_match = num_match + 1
                            exit
                        end if

                    end do !ipoint_bc
                end do !ipoint_f

                ! if all points are included in the current boundary condition
                ! patch, then set boundary condition type
                if (num_match == 4) then
                    face(8,iface) = -bc_type
                end if

            end do !iface

            deallocate(bc_x,bc_y,bc_z)

        end do !ifiles


    end subroutine read_bc


    subroutine create_bc_list()

      use mod_grid, only: npoin, mod_grid_get_face_ngl, nface, face, intma, coord
      
      use mod_face, only: imapl

      use mod_parallel, only : num_send_recv_total

      use mpi

      implicit none

      !Global Arrays    
!      integer, dimension(:), allocatable :: bc_list
!      integer :: bc_count
      
     
      !Local Arrays
      integer :: iface, ier, ilocl, iel, i, j, k, il, jl, kl, ip
      integer :: ngl_i, ngl_j,plane_ij, irank, ierr
      real, dimension (4,num_send_recv_total) :: bc_marker_nbh
      real, dimension(4,npoin) :: bc_marker  !one marker for velocity, one for temperature, one for salinity, one for mesh velocity
      integer, dimension(4) :: bcflag !same here

      call mpi_comm_rank(mpi_comm_world,irank,ierr)

      bc_marker = 0.0
      
      !------------------------------
      !loop over boundary points
      !------------------------------
      do iface=1, nface
         ier=face(8,iface)
         if (ier >= 0) cycle
        
         ilocl=face(5,iface)
         iel=face(7,iface)
         call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

         !----  boundary condition for U, T, S
         bcflag(1) = mod(ier,10)
         bcflag(2) = ier/10
         bcflag(2) = mod(bcflag(2),10)
         bcflag(3) = ier/100
         bcflag(3) = mod(bcflag(3),100)
         bcflag(4) = mod(ier,10)

         do j=1,ngl_j
            do i=1,ngl_i
               il=imapl(1,i,j,iface)
               jl=imapl(2,i,j,iface)
               kl=imapl(3,i,j,iface)
               ip=intma(il,jl,kl,iel)
               
               do k=1,4 !loop over different variables (U,T,S) and mesh velocity

                  if (bcflag(k) == -3) then !non-reflecting
                     !do nothing
                  else if (bcflag(k) == -4) then !no-flux (velocity-neumann, pressure-dirichlet)
                     !don't know how to classify this BC
                  else if (bcflag(k) == -6) then !outflow
                     !same as above
                  else if (bcflag(k) == -5) then !inflow (velocity-dirichlet, pressure-neumann)
                     bc_marker(k,ip) = 1
                  end if

                  if(k==4.and.bcflag(k)<0) then
                     bc_marker(k,ip) = 1 !mark all boundary points
                  end if
               end do
            end do
         end do
      end do


      call create_global_rhs(bc_marker,bc_marker_nbh,4,0)

      !remove points that lie on boundary faces to keep only boundary vertices (that lie on a non-boundary face of an element)
      do iface=1, nface
         ier=face(8,iface)
         if (ier >= 0) cycle
        
         ilocl=face(5,iface)
         iel=face(7,iface)
         call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

         !----  boundary condition for U, T, S
         bcflag(1) = mod(ier,10)

         !temperature
         bcflag(2) = ier/10
         bcflag(2) = mod(bcflag(2),10)

         !salinity
         bcflag(3) = ier/100
         bcflag(3) = mod(bcflag(3),100)

         bcflag(4) = mod(ier,10)                

         do j=1,ngl_j
            do i=1,ngl_i
               il=imapl(1,i,j,iface)
               jl=imapl(2,i,j,iface)
               kl=imapl(3,i,j,iface)
               ip=intma(il,jl,kl,iel)

               do k=1,4
                 
                  if (bcflag(k) == -3) then !non-reflecting
                     !do nothing
                  else if (bcflag(k) == -4) then !no-flux (velocity-neumann, pressure-dirichlet)
                     !don't know how to classify this BC
                  else if (bcflag(k) == -6) then !outflow
                     !same as above
                  else if (bcflag(k) == -5) then !inflow (velocity-dirichlet, pressure-neumann)
                     bc_marker(k,ip) = 0
                  end if

                  if(k==4.and.bcflag(k)<0) bc_marker(k,ip)=0

               end do
            end do
         end do
      end do

      bc_count=0
      !count how many markers we need
      do k=1,4
         do i=1,npoin
            if(bc_marker(k,i).gt.0) bc_count(k)=bc_count(k)+1
         end do
      end do

      if(allocated(bc_list)) deallocate(bc_list)
      allocate(bc_list(4,maxval(bc_count)))
      bc_list=0
      do k=1,4
         ip=0
         do i=1,npoin
            if(bc_marker(k,i).gt.0) then
               ip=ip+1
               bc_list(k,ip)=i
            end if
         end do
      end do

!      bc_count=0

!       do k=1,4
!          print*,"bc_list",k,irank,"::",bc_count(k),"->",bc_list(k,:)
!          do i=1,bc_count(k)
!             print*,"bc point",irank,bc_list(k,i),coord(:,bc_list(k,i))
!          end do
!       end do
    end subroutine create_bc_list


end module mod_bc
