!---------------------------------------------------------------------!
!>@brief This subroutine writes a GKS file for FXG's GKS.F plotting package.
!>@author  F.X. Giraldo on 8/97 
!>@date James F. Kelly Modified for 3D
!>@date 11/2009 F.X. Giraldo Rewritten to spit out a 2D slice 
!>@date 8/2010 James F. Kelly Parallelization
!---------------------------------------------------------------------!
subroutine outmesh_2d_g(q,fname,time)
  
    use mod_basis, only: ngl, nglx, ngly, nglz
  
    use mod_constants, only: earth_radius

    use mod_face, only: imapl
  
    use mod_global_grid, only: coord_g_cg, intma_g, npoin_g_cg
  
    use mod_grid, only: coord_cg, intma, bsido, npoin_cg, nelem, nbsido, face, nface
  
    use mod_initial, only: nvar
  
    use mod_input, only: icase, geometry_type
  
    use mod_parallel, only: nproc
  
    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none
  
    !global
    real q(nvar,npoin_cg), time
    character fname*72
  
    !local
    real yt, y, ymax, ymin, ycheck, rfactor
    real xmax, xmin
    integer face_pointer(nface), ip_pointer(npoin_cg), ip_pointer_back(npoin_cg)
    integer intma2d(ngl,ngl,nelem)
    integer ivar, ifactor, nop, ip, ie, i, j, k, mi, mj, mk, ix, iy, iz
    integer npoin2d, nelem2d, nboun2d, iface, iel, ilocl, il, jl, kl, iflag
    character text*72
  
    !MPI Variables
    integer irank, ierr
    !Global Data
    real, dimension(:,:), allocatable :: q_g, q_s
    integer, dimension(:), allocatable :: ip_pointer_g
    integer npoin2d_g, npoin2dmax
    integer npoin2d_l(nproc)
    integer nelem2d_g, ngl_i, ngl_j

    !Initialize
    intma2d=0

    !Get IRANK (should use MPI_UTILITY for this
    call mpi_comm_rank(mpi_comm_world,irank,ierr)

    !Extract 2D Slice
    if (irank == 0) then
        xmax=maxval(coord_g_cg(1,:))
        xmin=minval(coord_g_cg(1,:))
        ymax=maxval(coord_g_cg(2,:))
        ymin=minval(coord_g_cg(2,:))
    end if

    call mpi_bcast(xmax,1,MPI_PRECISION,0,mpi_comm_world,ierr)
    call mpi_bcast(ymax,1,MPI_PRECISION,0,mpi_comm_world,ierr)
    call mpi_bcast(xmin,1,MPI_PRECISION,0,mpi_comm_world,ierr)
    call mpi_bcast(ymin,1,MPI_PRECISION,0,mpi_comm_world,ierr)

    ip_pointer=0
    ip_pointer_back=0
    face_pointer=0
  
    !Determine Which variable to Print
    ivar=5
    ifactor=1
    yt=0.5*(ymax + ymin)
    ix=1; iy=2; iz=3
    select case (icase)
        case (0) !Passive Advection
            ivar=1
            ifactor=1
        case (1) !Inertia-Gravity Wave
            ifactor=10
        case (2) !2D Rising Thermal Bubble (in x-z)
            yt=0.5*(ymax + ymin)
        case (20) !2D Large Rising Thermal Bubble (in x-z)
            yt=0.5*(ymax + ymin)
        case (200) !2D Rising Thermal Bubble with Tracer (in x-z)
            yt=0.5*(ymax + ymin)
        case (21) !3D Rising Thermal Bubble
            yt=0.5*(ymax + ymin)
        case (22) !3D Rising Thermal Bubble:Catalano
            yt=0.5*(ymax + ymin)
        case (23) !2D Rising Thermal Bubble (in y-z)
            yt=0.5*( xmax + xmin )
            ix=2
            iy=1
            iz=3
        case (4) !Linear Hydrostatic Ridge (in x-z)
            ivar=4
            yt=0.5*( ymax + ymin )
        case (40) !Linear Hydrostatic Ridge (in x-z)
            ivar=4
            yt=0.5*( ymax + ymin )
        case (41) !Linear Hydrostatic Mountain
            ivar=4
            yt=0.5*( ymax + ymin )
        case (42) !Linear Hydrostatic Ridge (in y-z)
            ivar=4
            yt=0.5*( xmax + xmin )
            ix=2
            iy=1
            iz=3
        case (43) !Linear Hydrostatic Mountain
            ivar=4
            yt=0.5*( ymax + ymin )
        case (5) !3D Rising Bubble on Sphere
            yt=0
            ix=1
            iy=2
            iz=3
        case (8) !Gravity-Waves on the sphere
            yt = 0
            ix = 1
            iy = 2
            iz = 3
    end select
  
    !Normalize Coordinates if on Sphere
    rfactor=1
    if (geometry_type(1:6) == 'sphere') rfactor=earth_radius
  
    !Extract Faces on the Slice
    npoin2d=0
    nelem2d=0
    do iface=1,nface
        iflag=0
        do j=1,4
            y=coord_cg(iy,face(j,iface))
            ycheck=(y-yt)
            if (abs(y-yt) > 1e-5) then
                iflag=1
                exit
            end if
        end do !j
     
        !Store Face
        if (iflag == 0) then
            nelem2d=nelem2d + 1
            face_pointer(nelem2d)=iface
            iel=face(7,iface)
            ilocl=face(5,iface)

            !Check the Dimensions of the Face
            if (ilocl == 1 .or. ilocl == 2) then !Zeta=+/-1
                ngl_i=nglx
                ngl_j=ngly
            elseif (ilocl == 3 .or. ilocl == 4) then !Eta=+/-1
                ngl_i=nglx
                ngl_j=nglz
            elseif (ilocl == 5 .or. ilocl == 6) then !Ksi=+/-1
                ngl_i=ngly
                ngl_j=nglz
            end if

            !Loop through Face points and Form Pointer
            do j=1,ngl_j
                do i=1,ngl_i
                    il=imapl(1,i,j,iface)
                    jl=imapl(2,i,j,iface)
                    kl=imapl(3,i,j,iface)
                    ip=intma(il,jl,kl,iel)
                    intma2d(i,j,nelem2d)=ip
                    if (ip_pointer_back(ip)==0) then
                        npoin2d=npoin2d + 1
                        ip_pointer(npoin2d)=ip
                        ip_pointer_back(ip)=npoin2d
                    end if
              
                end do !i
            end do !j
        end if
    end do !i
  
    !>  print *, irank, npoin2d, nelem2d
    ! Build a Slice
    allocate(q_s(nvar,npoin2d))
    do i = 1,npoin2d
        ip=ip_pointer(i)
        q_s(:,i) = q(:,ip)
    end do
  
    !Get the maximum number of slice points on each processors
    call mpi_reduce(npoin2d,npoin2dmax,1,mpi_integer,mpi_max,0, &
        mpi_comm_world,ierr)
    !Get the total number of points
    call mpi_reduce(npoin2d,npoin2d_g,1,mpi_integer,mpi_sum,0, &
        mpi_comm_world,ierr)
    !Get the total number of elements
    call mpi_reduce(nelem2d,nelem2d_g,1,mpi_integer,mpi_sum,0, &
        mpi_comm_world,ierr)
  
    if (irank == 0) then
        allocate(q_g(nvar,npoin2d_g))
        allocate(ip_pointer_g(npoin2d_g))
    end if
  
    ! Gather Data onto Head node
    call gather_slice_data(q_g,q_s,nvar,ip_pointer,npoin2d,npoin2d_l,npoin2dmax,npoin2d_g,ip_pointer_g)
  
    !Constants
    nboun2d=0
  
    !Open GKS file
    if (irank == 0) then
        open(1,file=fname)
     
        text=" Time "
        write(1,'(a)')text
        write(1,'(e16.8)')time
     
        text=" NPOIN  NELEM NBOUN NOPX NOPY NOPZ "
        write(1,'(a)')text
        write(1,'(4(i7,1x))')npoin2d_g,nelem2d_g,nboun2d,nglx-1,ngly-1,nglz-1
     
        text=" COORDINATES "
        write(1,'(a)')text
        do i=1,npoin2d_g
            ip=ip_pointer_g(i)
            write(1,'(i7,1x,2(e16.8,1x))')i,coord_g_cg(ix,ip)/ifactor/rfactor, coord_g_cg(iz,ip)/rfactor
        end do !ip
     
        text=" ELEMENT CONNECTIVITY "
        write(1,'(a)')text
        do ie=1,nelem2d
            do j=1,ngl
                mj=(j-1)*ngl
                do i=1,ngl
                    if (intma2d(i,j,ie) > 0) then
                        mi=mj + i
                        ip=ip_pointer_back(intma2d(i,j,ie))
                        write(1,'(3(i7,1x))')ie,mi,ip
                    end if
                end do !i
            end do !j
        end do !ip
     
        text=" BOUNDARY SIDES "
        write(1,'(a)')text
        do i=1,nboun2d
            write(1,*)i, (bsido(j,i), j=1,4)
        end do
     
        text=" DEPENDENT VARIABLES: Q "
        write(1,'(a)')text
        do i=1,npoin2d_g
            write(1,*)i,q_g(1,i),q_g(ix+1,i),q_g(iz+1,i),q_g(5,i)
        end do
     
        close(1)
     
    end if
  
end subroutine outmesh_2d_g

!----------------------------------------------------!
!>@brief Gather slice data
!>@author F.X. Giraldo on 8/97 
!----------------------------------------------------!
subroutine gather_slice_data(q_g,q,nvar,ip_pointer, &
    npoin2d,npoin2d_l,npoin2dmax,npoin2d_g,ip_pointer_g)
  
    use mod_global_grid, only: npoin_g_cg, nelem_g

    use mod_grid, only: npoin_cg

    use mod_parallel, only: local_global_poin_l, npoin_l, nproc

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION
  
    integer npoin2d, npoin2dmax, npoin2d_g
    integer npoin2d_l(nproc)
    integer ip_pointer(npoin2d)
    integer ip_pointer_g(npoin2d_g)
    real q(nvar,npoin2d)
    real q_g(nvar,npoin2d_g)
    integer displs1(nproc), displs2(nproc)
    integer ierr, irank, iproc, ip, ip_g
    integer ipp, ii, jj, isnew, nvar
    real, dimension(:,:,:), allocatable :: q_l
    integer, dimension(:,:), allocatable :: ip_pointer_l
  
    call mpi_comm_rank(mpi_comm_world,irank,ierr)
  
    !Gather number of slice points on head node
    call mpi_gather(npoin2d,1,mpi_integer,npoin2d_l,1,mpi_integer, &
        0, mpi_comm_world,ierr)
  
    if (irank == 0) then
        allocate (q_l(nvar,npoin2dmax,nproc))
        allocate(ip_pointer_l(npoin2dmax,nproc))
    end if
  
    do iproc = 1,nproc
        displs1(iproc) = nvar*(iproc - 1)*npoin2dmax
        displs2(iproc) = (iproc - 1)*npoin2dmax
    end do
    !print *, "DISPLS computed"
  
    ! Gather all the data onto the head node
    call mpi_gatherv(q,nvar*npoin2d,MPI_PRECISION, &
        q_l,nvar*npoin2d_l, displs1,MPI_PRECISION, 0, mpi_comm_world,ierr)
  
    call mpi_gatherv(ip_pointer,npoin2d,mpi_integer, &
        ip_pointer_l,npoin2d_l, displs2,mpi_integer, 0, mpi_comm_world,ierr)
  
    !print *, "Data Gathered"
    if (irank ==0) then
        ii = 0
        do iproc = 1,nproc
            do ipp = 1,npoin2d_l(iproc)
           
                ip = ip_pointer_l(ipp,iproc)
                ip_g = local_global_poin_l(ip,iproc)
                isnew = 1
                do jj = 1,ii
                    if (ip_g == ip_pointer_g(jj)) then
                        isnew = 0
                    end if
                end do
                if (isnew == 1) then
                    ii = ii + 1
                    ip_pointer_g(ii) = ip_g
                    q_g(:,ii) = q_l(:,ipp,iproc)
                end if
           
            end do
        end do
        npoin2d_g = ii
        deallocate(q_l)
    end if
  
end subroutine gather_slice_data


