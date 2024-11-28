!----------------------------------------------------------------------!
!>@brief This subroutine prints out Diagnostics
!>@author  Francis X. Giraldo on 11/2009
!>           Department of Applied Mathematics
!>           Naval Postgraduate School
!>           Monterey, CA 93943-5216
!>
!>@ modified by Yao Gahounzo 
!>      Computing PhD 
!       Boise State University
!       Date: April 03, 2024
!----------------------------------------------------------------------!

subroutine print_diagnostics_mlswe(q_mlswe,qb,time,itime,dt,idone,&
   mass_conserv0_g,cfl,cflu,ntime,fnp11, unit0)

   use mpi

   use mod_constants, only: nnorm

   use mod_global_grid, only: npoin_g

   use mod_grid, only: npoin, ncol

   use mod_initial, only: nvar

   use mod_input, only: lprint_diagnostics, si_dimension, ti_method, &
      icase, fname_root, time_scale, nlayers, dt_btp, lcheck_conserved

   use mod_mpi_utilities, only: MPI_PRECISION

   implicit none

   !global arrays
   real, intent(in)  :: q_mlswe(5,npoin,nlayers), qb(4,npoin)
   real, intent(in)  :: time, dt, mass_conserv0_g(nlayers)
   integer, intent(in) :: itime, idone, ntime, unit0
   character, intent(in) :: fnp11*18

   !local arrays
   real, dimension(:,:), allocatable :: q
   real, dimension(nvar,nlayers) :: qmax_layers, qmin_layers
   real, dimension(4,nlayers) :: cfl_vector_layers
   real :: qmax(nvar),   qmin(nvar), qbmax(4), qbmin(4)
   real :: qmax_g(nvar), qmin_g(nvar), qbmax_g(4), qbmin_g(4)
   real :: cfl_vector(5), cfl_vector_g(5), cfl, cflu
   real :: min_dx_vec(2),min_dx_vec_g(2)
   real :: xm1(nlayers)
   integer :: ierr, irank, i, j, ncol_g, m, ll
   real :: mass_conserv_l, mass_conserv(nlayers), mass_conserv_g
   character :: fname_conservation*72,fname_qmax_qmin*72
   character(len=3), dimension(5) :: fileds

   allocate(q(nvar,npoin))

   !Get Processor ID number
   call mpi_comm_rank(mpi_comm_world,irank,ierr)


   !possibly start layers loop here 

   do ll=1,nlayers

      q = q_mlswe(:,:,ll)

      !Compute Mass

      if(lcheck_conserved) then 
      
         call compute_conserved(mass_conserv_l,q(1,:))
         mass_conserv_g = 0.0

         call mpi_reduce(mass_conserv_l,mass_conserv_g,1,MPI_PRECISION,mpi_sum,0,mpi_comm_world,ierr)

         mass_conserv(ll) = mass_conserv_g
      end if
      
      !Calculate max on processor
       
      do i=1,nvar
         qmax(i)=maxval(q(i,:))
         qmin(i)=minval(q(i,:))
      end do !i
       
       
      !Get Global Max
      call mpi_reduce(qmax,qmax_g,nvar,MPI_PRECISION,&
         mpi_max,0,mpi_comm_world,ierr)
      
      !Get Global Min
      call mpi_reduce(qmin,qmin_g,nvar,MPI_PRECISION,&
         mpi_min,0,mpi_comm_world,ierr)
      
      qmax_layers(:,ll) = qmax_g
      qmin_layers(:,ll) = qmin_g

   end do
   

   if (irank == 0 .and. lcheck_conserved) then 
      
      do ll = 1,nlayers
         xm1(ll) = abs(mass_conserv(ll) - mass_conserv0_g(ll))/mass_conserv0_g(ll)
      end do 

      ! Save mass conservation to a file
      if (idone == 0) write(unit0, fnp11) itime, mass_conserv

   end if

   do i = 1,4
      qbmax(i)=maxval(qb(i,:))
      qbmin(i)=minval(qb(i,:))
   end do

   !Get Global Max
   call mpi_reduce(qbmax,qbmax_g,4,MPI_PRECISION,&
      mpi_max,0,mpi_comm_world,ierr)
   
   !Get Global Min
   call mpi_reduce(qbmin,qbmin_g,4,MPI_PRECISION,&
      mpi_min,0,mpi_comm_world,ierr)


   !Each Proc. computes a CFL
   call courant_mlswe(cfl_vector,q_mlswe,qb,dt,dt_btp,nlayers,min_dx_vec)
      
   call mpi_reduce(cfl_vector,cfl_vector_g,5,MPI_PRECISION,&
      mpi_max,0,mpi_comm_world,ierr)
   cfl=max(cfl_vector_g(1),cfl_vector_g(2))
   cflu=max(cfl_vector_g(3),cfl_vector_g(4))

   call mpi_reduce(min_dx_vec,min_dx_vec_g,2,MPI_PRECISION,&
      mpi_min,0,mpi_comm_world,ierr)
   
   !cfl_vector_layers(:,ll) = cfl_vector_g

   if (irank == 0 .and. idone == 0) then
      print*,'==============================================================='
      write(*,'("itime time dt dt_btp = ",i8,1x,2(es13.5,1x),2(es13.5,1x))')itime,time/time_scale, dt, dt_btp
      write(*,'("CFL_H = ",e11.4," CFL_B = ",e11.4," CFL = ",f9.4)')cfl_vector_g(1),cfl_vector_g(2), cfl_vector_g(5)
      write(*,'("CFLU_H = ",e11.4," CFLU_B = ",e11.4)')cfl_vector_g(3),cfl_vector_g(4)
      write(*,'("dx_min = ",e11.4," dy_min = ",e11.4)')min_dx_vec_g(1),min_dx_vec_g(2)
      print*,'---------------------------------------------------------------'
      do ll = 1,nlayers
         write(*,'("Layer = ",i8)')ll
         write(*,'("Mass Loss   = ",1(e16.8,1x))') xm1(ll)
         do i=1,nvar
            write(*,'("Q: i    Max/Min = ",i3,1x,2(e24.12,1x))')i,qmax_layers(i,ll), qmin_layers(i,ll)
         end do !i
         print*,'---------------------------------------------------------------'
      end do
      print*,'---------------------------------------------------------------'
      write(*,*)'Barotropic'
      do i=1,4
         write(*,'("Qb: i    Max/Min = ",i3,1x,2(e24.12,1x))')i,qbmax_g(i), qbmin_g(i)
      end do !i
        
      print*,'==============================================================='
   else if (irank == 0 .and. idone == 1) then
      print*,'---------------------------------------------------------------'
      write(*,'(" **Simulation Finished**")')
      write(*,'("itime time dt dt_btp = ",i8,1x,2(es13.5,1x),2(es13.5,1x))')itime,time/time_scale,dt, dt_btp
      write(*,'("CFL_H = ",e11.4," CFL_B = ",e11.4," CFL = ",f9.4)')cfl_vector_g(1),cfl_vector_g(2), cfl_vector_g(5)
      write(*,'("CFLU_H = ",e11.4," CFLU_B = ",e11.4)')cfl_vector_g(3),cfl_vector_g(4)
      print*,'---------------------------------------------------------------'

      fileds(1) = "h"
      fileds(2) = "u"
      fileds(3) = "v"
      fileds(4) = "dp"
      fileds(5) = "ssh"
      
      open(unit=100, file = 'mlswe_FIN.txt')
      do ll = 1,nlayers
         write(*,'("Layer = ",i8)')ll
         write(*,'("Mass Loss  = ",1(e16.8,1x))') xm1(ll)
         write(100,'("Layer = ",i8)')ll
         write(100,'("Mass Loss  = ",1(e16.8,1x))') xm1(ll)
         do i=1,nvar
            write(*,'("Q: i    Max/Min = ",i3,1x,2(e24.12,1x))')i,qmax_layers(i,ll), qmin_layers(i,ll)

            if(i /= 4) then
               write(100,'("Fields:   Max/Min = ",(A),1x,2(e24.12,1x))')fileds(i),qmax_layers(i,ll), qmin_layers(i,ll)
            end if 
         end do !i
         print*,'---------------------------------------------------------------'
      end do
      close(100)

   end if !irank=0
        
   deallocate(q)
  
end subroutine print_diagnostics_mlswe
