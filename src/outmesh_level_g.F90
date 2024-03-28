!---------------------------------------------------------------------!
!>@brief This subroutine extracts a vertical level, and writes variables to a file
!>@author James F. Kelly 6 September 2010
!>@date 22 April 2011
!> itype = 1 ==> perturbations, itype = 2 ==> total fields
!---------------------------------------------------------------------!
subroutine outmesh_level_g(q,fname,time,ipoin_r,itype)
  
    use mod_basis, only: ngl, nglx, ngly, nglz
  
    use mod_constants, only: pi, p00, kappa

    use mod_global_grid, only: coord_g_cg, intma_g, npoin_g_cg, intma_s, npoin_s, intma_lev_g, nelem_s
  
    use mod_grid, only: coord_cg, intma, bsido, npoin_cg, nelem, nbsido
  
    use mod_initial, only: nvar, q_ref

    use mod_input, only: time_scale

    use mod_parallel, only: nproc, local_global_poin

    use mpi

    use mod_mpi_utilities, only: MPI_PRECISION

    implicit none
  
    !global
    real q(nvar,npoin_cg), time
    character fname*72
    integer ipoin_r
    real rlevel
    integer itype
  
    real, dimension(:,:), allocatable :: coord_s, q_s, q_g
    real, dimension(:,:,:), allocatable :: q_l
    integer, dimension(:), allocatable :: level_global, ip_pointer, displs1, displs2, npoin2d_l
    integer, dimension(:,:), allocatable :: ip_pointer_l
    

    !local
    integer irank, ierr, iproc, ie, mj,nboun2d, mi, nop
    integer ips, i, j, ip_s, ie_s, ip_g, ip, ipp
    real x, y, z, r, r2, pressure_total, press
    real olon, olat, uc, vc, wc, us, vs, eps, pio4
    real rho, theta, temp
    integer npoin2dmax, npoin2d
    character text*72
  
    call mpi_comm_rank(mpi_comm_world,irank,ierr)

    !Construct Level-Based Coords on the Head Node
    call mpi_bcast(npoin_s,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nelem_s,1,mpi_integer,0,mpi_comm_world,ierr)
    allocate(level_global(npoin_s))
  
    if (irank == 0) then
        allocate(coord_s(3,npoin_s))

        do ie_s = 1,nelem_s
            do j = 1,ngly
                do i = 1,nglx
                    ip_s = intma_s(i,j,ie_s)      !Get Grid Point on the sphere
                    ip_g = intma_lev_g(i,j,ie_s,ipoin_r)    !Get Global Grid Point
                    coord_s(:,ip_s) = coord_g_cg(:,ip_g)       !Populate Data on the sphere
                    level_global(ip_s) = ip_g               !Map Level --> Global Grid Points
                end do
            end do
        end do
    end if

    !Allocate Space for the Data on each processor
    allocate(q_s(nvar,npoin_s),ip_pointer(npoin_s))

    !call mpi_barrier(mpi_comm_world,ierr)
    call mpi_bcast(level_global,npoin_s,mpi_integer,0,mpi_comm_world,ierr)
  
    !Extract Data on the specified vertical level
    if (itype == 1) then
        ipp = 0
        do ip = 1,npoin_cg
            ip_g = local_global_poin(ip)
            do ip_s = 1,npoin_s
                if (level_global(ip_s) == ip_g) then
                    ipp = ipp + 1
                    q_s(:,ipp) = q(:,ip)
                    ip_pointer(ipp) = ip_s
                end if
            end do
        end do
    end if

    if (itype == 2) then
        ipp = 0
        do ip = 1,npoin_cg
            ip_g = local_global_poin(ip)
            do ip_s = 1,npoin_s
                if (level_global(ip_s) == ip_g) then
                    ipp = ipp + 1
                    q_s(1,ipp)   = q(1,ip) + q_ref(1,ip)
                    q_s(2:4,ipp) = q(2:4,ip)
                    q_s(5,ipp)   = q(5,ip) + q_ref(5,ip)
                    ip_pointer(ipp) = ip_s
                end if
            end do
        end do
    end if

    npoin2d = ipp

    ! Get the maximum number of level points
    call mpi_reduce(npoin2d,npoin2dmax,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  
    ! Allocate Memory on the Head Node
    if (irank == 0) then
        allocate(q_l(nvar,npoin2dmax,nproc))
        allocate(q_g(nvar,npoin_s))
        allocate(npoin2d_l(nproc))
        allocate(ip_pointer_l(npoin2dmax,nproc))
    end if

    !Gather number of data points onto the head node
    call mpi_gather(npoin2d,1,mpi_integer,npoin2d_l,1,mpi_integer, &
        0, mpi_comm_world,ierr)
  
    allocate(displs1(nproc),displs2(nproc))
    ! Consturct Pointers for MPI_GATHERV
    do iproc = 1,nproc
        displs1(iproc) = (iproc - 1)*npoin2dmax*nvar
        displs2(iproc) = (iproc - 1)*npoin2dmax
    end do

    ! Gather Sliced Data onto the head node
    call mpi_gatherv(q_s,nvar*npoin2d,MPI_PRECISION, &
        q_l,nvar*npoin2d_l, displs1,MPI_PRECISION, 0, mpi_comm_world,ierr)
    call mpi_gatherv(ip_pointer,npoin2d,mpi_integer, &
        ip_pointer_l,npoin2d_l, displs2,mpi_integer, 0, mpi_comm_world,ierr)
 
    !OK, Now Construct the data in the correct order
    if (irank == 0) then
        do iproc = 1,nproc
            do ipp = 1,npoin2d_l(iproc)
                ip_s = ip_pointer_l(ipp,iproc)   !Get Local Grid Point
                q_g(:,ip_s) = q_l(:,ipp,iproc)
            end do
        end do
  
    ! Convert to NSEAM-type variables (pressure, zonal vel., meridional vel., and temp)
    !call nseam_out(q_g,coord_s,npoin_s)

    end if

    ! Head Node Writes Out Data
    if (irank == 0) then
        open(1,file=fname)

        text=" Time "
        write(1,'(a)')text
        write(1,'(e16.8)')time/time_scale

        text=" NPOIN  NELEM NOP "
        write(1,'(a)')text
        write(1,'(4(i7,1x))')npoin_s,nelem_s,nglx-1

        text=" COORDINATES "
        write(1,'(a)')text
        do i=1,npoin_s
            r=sqrt( dot_product(coord_s(:,i),coord_s(:,i)) )
            write(1,'(i7,1x,2(e16.8,1x))')i,coord_s(1,i)/r,coord_s(2,i)/r, coord_s(3,i)/r
        end do

        text=" ELEMENT CONNECTIVITY "
        write(1,'(a)')text
        do ie=1,nelem_s
            do j=1,ngly
                mj=(j-1)*nglx
                do i=1,nglx
                    mi=mj + i
                    ip= intma_s(i,j,ie)
                    write(1,'(3(i7,1x))')ie,mi,ip
                end do !i
            end do !j
        end do !ie

        eps=1.0e-6
        pio4=pi/4
        text=" DEPENDENT VARIABLES: Q "
        write(1,'(a)')text
        do i=1,npoin_s
            x=coord_s(1,i)
            y=coord_s(2,i)
            z=coord_s(3,i)
            r=sqrt( x*x + y*y + z*z )
            x=x/r
            y=y/r
            z=z/r
            olon=atan2(y,x+eps)
            olat=asin(z)
        
            !Compute Spherical Velocities
            rho=q_g(1,i)
            uc=q_g(2,i)
            vc=q_g(3,i)
            wc=q_g(4,i)
            theta=q_g(5,i)
            if (itype == 2) then
                press=pressure_total(rho,theta)
                temp=theta*(press/p00)**kappa
            end if
            us=vc*cos(olon) - uc*sin(olon)
            vs=-uc*cos(olon)*sin(olat) - vc*sin(olon)*sin(olat) + wc*cos(olat)
            if (itype == 1) then
                write(1,*)i,rho,us,vs,theta
            else if (itype == 2) then
                write(1,*)i,press,us,vs,temp
            end if
        end do !i
        close(1)
    end if

end subroutine outmesh_level_g

!----------------------------------------------------!
!>@brief Calculates total pressure
!----------------------------------------------------!
function pressure_total(rho,theta)
  
    use mod_constants, only: rgas, p00, cp, cv

    implicit none

    !global arrays
    real c, rho, theta, pressure_total

    pressure_total=p00*(rho*rgas*theta/p00)**(cp/cv)

end function pressure_total
