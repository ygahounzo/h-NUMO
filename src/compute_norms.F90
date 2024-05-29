!----------------------------------------------------------------------!
!>@brief This subroutine computes the L1, L2, and L infinity Norms
!>@author  Francis X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: July 02, 2023
!----------------------------------------------------------------------!
subroutine compute_norms(l1_norm,l2_norm,l8_norm,time_norm,q0,inorm,time)
  
    use mod_basis, only: nglx, ngly, nglz

    use mod_constants, only: nnorm, tol
  
    use mod_convert_variables, only: mod_convert_variables_eqnset_to_set2nc

    use mod_grid, only: intma, npoin, nelem
  
    use mod_initial, only : q_exact, q_ref, nvar

    use mod_input, only : icase

    use mod_metrics, only: jac
  
    use mod_parallel, only :ipoin_proc

    use mpi
  
    implicit none
  
    !Global Arrays
    real, intent(in)    :: q0(nvar,npoin), time
    real, intent(inout) :: l1_norm(nvar,nnorm), l2_norm(nvar,nnorm), l8_norm(nvar,nnorm), time_norm(nnorm)
    integer, intent(in) :: inorm
  
    !Local Arrays
    real, dimension(:,:), allocatable :: qn, qe
    real, dimension(:),   allocatable :: l2_norm_top, l2_norm_bot,l1_norm_top, l1_norm_bot,l8_norm_top, l8_norm_bot
    real, dimension(:),   allocatable :: l2_norm_top_g, l2_norm_bot_g,l1_norm_top_g, l1_norm_bot_g,l8_norm_top_g, l8_norm_bot_g
    real q_max, qe_max, qee
    real wq, dq(nvar)
    integer AllocateStatus
    integer ip, i, j, k, l, e, ierr, irank
  
    call mpi_comm_rank(mpi_comm_world,irank,ierr)
  
    allocate( qn(nvar,npoin), qe(nvar,npoin), l2_norm_top(nvar), l2_norm_bot(nvar), &
        l1_norm_top(nvar), l1_norm_bot(nvar), l8_norm_top(nvar), l8_norm_bot(nvar), stat=AllocateStatus )
    if (AllocateStatus /= 0) stop "** Not Enough Memory - COMPUTE_NORMS Local **"
    if (irank == 0) then
        allocate( l2_norm_top_g(nvar), l2_norm_bot_g(nvar),l1_norm_top_g(nvar), &
            l1_norm_bot_g(nvar),l8_norm_top_g(nvar), l8_norm_bot_g(nvar) , stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - COMPUTE_NORMS Global **"
    end if
  
    !Initialize
    l1_norm_top=0; l1_norm_bot=0
    l2_norm_top=0; l2_norm_bot=0
    l8_norm_top=0; l8_norm_bot=0

    !Compute Total Primitive Variables
    call mod_convert_variables_eqnset_to_set2nc(qn,q0,q_ref)
    call mod_convert_variables_eqnset_to_set2nc(qe,q_exact,q_ref)

    !Compute Solution Conservation Variables
    if (icase == 7) then
       qn=q0      + q_ref
       qe=q_exact + q_ref
        do i=1,npoin
            qn(2:5,i)=qn(1,i)*qn(2:5,i)
            qe(2:5,i)=qe(1,i)*qe(2:5,i)
        end do
    end if

    !loop thru the elements
    do e=1,nelem
        do k=1,nglz
            do j=1,ngly
                do i=1,nglx
                    ip=intma(i,j,k,e)
                    wq=jac(i,j,k,e)
                    do l=1,nvar
                        dq(l)=qn(l,ip) - qe(l,ip)
                        l1_norm_top(l)=l1_norm_top(l) + wq*abs(dq(l))
                        l2_norm_top(l)=l2_norm_top(l) + wq*dq(l)*dq(l)
                        l8_norm_top(l)=max( l8_norm_top(l), abs(dq(l)) )

                        !check for zero value
                        qee=abs(qe(l,ip))
                        if (qee < tol) qee=1.0
                        l1_norm_bot(l)=l1_norm_bot(l) + wq*qee
                        l2_norm_bot(l)=l2_norm_bot(l) + wq*qee*qee
                        l8_norm_bot(l)=max( l8_norm_bot(l), qee)
                    end do !l
                end do !i
            end do !j
        end do !k
    end do !e

    !Communicate Across all Processors
    call mpi_reduce(l1_norm_top,l1_norm_top_g,nvar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call mpi_reduce(l1_norm_bot,l1_norm_bot_g,nvar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call mpi_reduce(l2_norm_top,l2_norm_top_g,nvar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call mpi_reduce(l2_norm_bot,l2_norm_bot_g,nvar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call mpi_reduce(l8_norm_top,l8_norm_top_g,nvar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call mpi_reduce(l8_norm_bot,l8_norm_bot_g,nvar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)

    time_norm(inorm)=time
    if (irank == 0) then

        l1_norm(:,inorm) = l1_norm_top_g(:)/l1_norm_bot_g(:)
        l2_norm(:,inorm) = sqrt(l2_norm_top_g(:)/l2_norm_bot_g(:))
        l8_norm(:,inorm) = l8_norm_top_g(:)/l8_norm_bot_g(:)

        !Print out Norms
        print*,'---------------------------------------------------------------'
        print*,'---------ERROR NORMS-------------------------------------------'
        print*,'---------------------------------------------------------------'
        write(*,'("Rho:    L1 L2 L8 = ",3(e13.5,1x))')l1_norm(1,inorm),l2_norm(1,inorm),l8_norm(1,inorm)
        write(*,'("U:      L1 L2 L8 = ",3(e13.5,1x))')l1_norm(2,inorm),l2_norm(2,inorm),l8_norm(2,inorm)
        write(*,'("V:      L1 L2 L8 = ",3(e13.5,1x))')l1_norm(3,inorm),l2_norm(3,inorm),l8_norm(3,inorm)
        write(*,'("W:      L1 L2 L8 = ",3(e13.5,1x))')l1_norm(4,inorm),l2_norm(4,inorm),l8_norm(4,inorm)
        write(*,'("Theta:  L1 L2 L8 = ",3(e13.5,1x))')l1_norm(5,inorm),l2_norm(5,inorm),l8_norm(5,inorm)
        print*,'---------------------------------------------------------------'
        print*,'---------------------------------------------------------------'
    end if
  
end subroutine compute_norms


!----------------------------------------------------------------------!
!>@brief This subroutine computes the L1, L2, and L infinity Norms
!>@author  Francis X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>@date 6 September 2010 James F. Kelly
!> Parallelization
!>@date 2 May 2016 Michal A. Kopera
!> Simplified for a single vector
!----------------------------------------------------------------------!
subroutine compute_norms_vector(l1_norm,l2_norm,l8_norm,q0,qe)
  
  use mod_basis, only: nglx, ngly, nglz

  use mod_constants, only: nnorm, tol
  
  use mod_grid, only: intma, npoin, nelem

  use mod_input, only : icase 

  use mod_metrics, only: jac
  
  use mod_parallel, only :ipoin_proc
  
  use mpi
  
  implicit none
  
  !Global Arrays
  real, intent(in)    :: q0(npoin), qe(npoin)
  real, intent(inout) :: l1_norm, l2_norm, l8_norm
  
  !Local Arrays
  real:: l2_norm_top, l2_norm_bot,l1_norm_top, l1_norm_bot,l8_norm_top, l8_norm_bot
  real:: l2_norm_top_g, l2_norm_bot_g,l1_norm_top_g, l1_norm_bot_g,l8_norm_top_g, l8_norm_bot_g
  real q_max, qe_max, qee
  real wq, dq
  integer ip, i, j, k, l, e, ierr, irank
  
  call mpi_comm_rank(mpi_comm_world,irank,ierr)
  
  
  !Initialize
  l1_norm_top=0; l1_norm_bot=0
  l2_norm_top=0; l2_norm_bot=0
  l8_norm_top=0; l8_norm_bot=0


  !loop thru the elements
  do e=1,nelem
     do k=1,nglz
        do j=1,ngly
           do i=1,nglx
              ip=intma(i,j,k,e)
              wq=jac(i,j,k,e)

              dq=q0(ip) - qe(ip)
              l1_norm_top=l1_norm_top + wq*abs(dq)
              l2_norm_top=l2_norm_top + wq*dq*dq
              l8_norm_top=max( l8_norm_top, abs(dq) )
              
              !check for zero value
              qee=abs(qe(ip))
              if (qee < tol) qee=1.0
              l1_norm_bot=l1_norm_bot + wq*qee
              l2_norm_bot=l2_norm_bot + wq*qee*qee
              l8_norm_bot=max( l8_norm_bot, qee)
           
           end do !i
        end do !j
     end do !k
  end do !e   

  !Communicate Across all Processors
  call mpi_reduce(l1_norm_top,l1_norm_top_g,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(l1_norm_bot,l1_norm_bot_g,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(l2_norm_top,l2_norm_top_g,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(l2_norm_bot,l2_norm_bot_g,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(l8_norm_top,l8_norm_top_g,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(l8_norm_bot,l8_norm_bot_g,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)

  if (irank == 0) then 

     l1_norm = l1_norm_top_g/l1_norm_bot_g
     l2_norm = sqrt(l2_norm_top_g/l2_norm_bot_g)
     l8_norm = l8_norm_top_g/l8_norm_bot_g

  end if
  
end subroutine compute_norms_vector

