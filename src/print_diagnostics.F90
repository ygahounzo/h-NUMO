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

subroutine print_diagnostics_mlswe(G, inp, b, tsp, init, q_mlswe, qb, time, itime, dt, idone, &
   mass_conserv0_g, ntime, fnp11, unit0)

   use mpi
   use mod_mpi_utilities, only: MPI_PRECISION
   use mod_grid,          only: grid
   use mod_input,         only: input
   use mod_basis,         only: basis
   use mod_tensor,     only: tensor_CS
   use mod_initial,       only: initial

   implicit none

   type(grid),      intent(in) :: G
   type(input),     intent(in) :: inp
   type(basis),     intent(in) :: b
   type(tensor_CS), intent(in) :: tsp
   type(initial),   intent(in) :: init

   real, intent(in) :: q_mlswe(5, G%npoin, inp%nlayers), qb(4, G%npoin)
   real, intent(in) :: time, dt, mass_conserv0_g(inp%nlayers)
   integer, intent(in) :: itime, idone, ntime, unit0
   character, intent(in) :: fnp11*18

   real, dimension(:,:), allocatable :: q
   real, dimension(init%nvar, inp%nlayers) :: qmax_layers, qmin_layers
   real :: qmax(init%nvar),   qmin(init%nvar), qbmax(4), qbmin(4)
   real :: qmax_g(init%nvar), qmin_g(init%nvar), qbmax_g(4), qbmin_g(4)
   real :: cfl_vector(2), cfl_vector_g(2)
   real :: min_dx_vec(2), min_dx_vec_g(2)
   real :: xm1(inp%nlayers)
   integer :: ierr, irank, i, j, ncol_g, m, ll
   real :: mass_conserv_l, mass_conserv(inp%nlayers), mass_conserv_g
   character(len=3), dimension(5) :: fileds

   allocate(q(init%nvar, G%npoin))

   call mpi_comm_rank(mpi_comm_world, irank, ierr)

   do ll = 1, inp%nlayers

      q = q_mlswe(:,:,ll)

      if(inp%lcheck_conserved) then

         call compute_conserved(G, b, tsp, mass_conserv_l, q(1,:))
         mass_conserv_g = 0.0

         call mpi_reduce(mass_conserv_l, mass_conserv_g, 1, MPI_PRECISION, mpi_sum, &
                         0, mpi_comm_world, ierr)

         mass_conserv(ll) = mass_conserv_g
      end if

      do i = 1, init%nvar
         qmax(i) = maxval(q(i,:))
         qmin(i) = minval(q(i,:))
      end do

      call mpi_reduce(qmax, qmax_g, init%nvar, MPI_PRECISION, &
         mpi_max, 0, mpi_comm_world, ierr)

      call mpi_reduce(qmin, qmin_g, init%nvar, MPI_PRECISION, &
         mpi_min, 0, mpi_comm_world, ierr)

      qmax_layers(:,ll) = qmax_g
      qmin_layers(:,ll) = qmin_g

   end do

   if (irank == 0 .and. inp%lcheck_conserved) then

      do ll = 1, inp%nlayers
         xm1(ll) = abs(mass_conserv(ll) - mass_conserv0_g(ll)) / mass_conserv0_g(ll)
      end do

      if (idone == 0) write(unit0, fnp11) itime, mass_conserv

   end if

   do i = 1, 4
      qbmax(i) = maxval(qb(i,:))
      qbmin(i) = minval(qb(i,:))
   end do

   call mpi_reduce(qbmax, qbmax_g, 4, MPI_PRECISION, &
      mpi_max, 0, mpi_comm_world, ierr)

   call mpi_reduce(qbmin, qbmin_g, 4, MPI_PRECISION, &
      mpi_min, 0, mpi_comm_world, ierr)

   call courant_mlswe(G, b, cfl_vector, q_mlswe, qb, dt, inp%dt_btp, inp%nlayers, min_dx_vec)

   call mpi_reduce(cfl_vector, cfl_vector_g, 2, MPI_PRECISION, &
      mpi_max, 0, mpi_comm_world, ierr)

   call mpi_reduce(min_dx_vec, min_dx_vec_g, 2, MPI_PRECISION, &
      mpi_min, 0, mpi_comm_world, ierr)

   if (irank == 0 .and. idone == 0) then
      print *, '======================================================================='
      write(*,'("itime time dt dt_btp = ",i8,1x,2(es13.5,1x),2(es13.5,1x))') itime, &
               time/inp%time_scale, dt, inp%dt_btp
      write(*,'("CFL_B = ",e11.4," CFL = ",e11.4)') cfl_vector_g(1), cfl_vector_g(2)
      write(*,'("dx_min = ",e11.4," dy_min = ",e11.4)') min_dx_vec_g(1), min_dx_vec_g(2)
      print *, '-----------------------------------------------------------------------'
      do ll = 1, inp%nlayers
         write(*,'("Layer = ",i8)') ll
         write(*,'("Mass Loss   = ",1(e22.8,1x))') xm1(ll)
         do i = 1, init%nvar
            write(*,'("Q: i    Max/Min = ",i3,1x,2(e24.12,1x))') i, qmax_layers(i,ll), &
                     qmin_layers(i,ll)
         end do
         print *, '---------------------------------------------------------------------'
      end do
      print *, '------------------------------------------------------------------------'
      write(*,*) 'Barotropic'
      do i = 1, 4
         write(*,'("Qb: i    Max/Min = ",i3,1x,2(e24.12,1x))') i, qbmax_g(i), qbmin_g(i)
      end do

      print *, '========================================================================'
   else if (irank == 0 .and. idone == 1) then
      print *, '------------------------------------------------------------------------'
      write(*,'(" **Simulation Finished**")')
      write(*,'("itime time dt dt_btp = ",i8,1x,2(es13.5,1x),2(es13.5,1x))') itime, &
               time/inp%time_scale, dt, inp%dt_btp
      write(*,'("CFL_B = ",e11.4," CFL = ",e11.4)') cfl_vector_g(1), cfl_vector_g(2)
      print *, '------------------------------------------------------------------------'

      fileds(1) = "h"
      fileds(2) = "u"
      fileds(3) = "v"
      fileds(4) = "dp"
      fileds(5) = "ssh"

      open(unit=100, file = 'mlswe_FIN.txt')
      do ll = 1, inp%nlayers
         write(*,'("Layer = ",i8)') ll
         write(*,'("Mass Loss  = ",1(e16.8,1x))') xm1(ll)
         write(100,'("Layer = ",i8)') ll
         write(100,'("Mass Loss  = ",1(e16.8,1x))') xm1(ll)
         do i = 1, init%nvar
            write(*,'("Q: i    Max/Min = ",i3,1x,2(e24.12,1x))') i, qmax_layers(i,ll), &
                  qmin_layers(i,ll)
            if(i /= 4) then
               write(100,'("Fields:   Max/Min = ",(A),1x,2(e24.12,1x))') fileds(i), &
                  qmax_layers(i,ll), qmin_layers(i,ll)
            end if
         end do
         print *, '---------------------------------------------------------------------'
      end do
      close(100)

   end if

   deallocate(q)

end subroutine print_diagnostics_mlswe
