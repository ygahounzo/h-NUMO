!----------------------------------------------------------------------!
!>@brief This module contains routines to convert variables 
!> from SET2NC to Set2c and 3C and vice versa.
!>@author F.X. Giraldo on 1/2014
!> Department of Applied Mathematics
!> Naval Postgraduate School
!> Monterey, CA 93943-5216
!----------------------------------------------------------------------!
module mod_convert_variables

  use mod_constants, only: gravity, gamma, rgas, cp, cv, p00

  use mod_grid, only: npoin, coord

  use mod_initial, only: nvar, height, rho_layers

  use mod_input, only: eqn_set, icase, llimit, limit_threshold, is_shallow, equations, lsalinity, &
       is_swe_layers

  public :: &
       mod_convert_variables_set2nc_to_eqnset, &
       mod_convert_variables_eqnset_to_set2nc, &
       mod_convert_variables_eqnset_to_primitive_perturbation, &
       mod_convert_variables_primitive_perturbation_to_eqnset

contains

  !----------------------------------------------------------------------!
  !>@brief This subroutine converts variables from SET2NC to EQNSET: e.g., used in INITIAL_CONDITIONS
  !>@author  Francis X. Giraldo on 1/2014
  !>           Department of Applied Mathematics
  !>           Naval Postgraduate School
  !>           Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine mod_convert_variables_set2nc_to_eqnset(qout,qout_ref,qout_exact,q,q_ref,q_exact,pi_values)

    implicit none

    !global arrays
    real, intent(in) :: q(nvar,npoin), q_ref(nvar,npoin), q_exact(nvar,npoin), pi_values(3,npoin)
    real, intent(out):: qout(nvar,npoin), qout_ref(nvar,npoin), qout_exact(nvar,npoin)

    !Local arrays
    real :: rho_k, u_k, v_k, w_k, theta_k, e_k, p_k, t_k, pi_k
    real :: rho_ref, u_ref, v_ref, w_ref, theta_ref, e_ref, p_ref, t_ref, pi_ref
    real :: rho_exact, u_exact, v_exact, w_exact, theta_exact, e_exact, p_exact, t_exact, pi_exact
    real :: radius, theta0
    real :: c
    integer :: i, m

    c = cv/rgas

    if (eqn_set(1:6) == 'set2nc') then
       !Do Nothing
    else if (eqn_set(1:5) == 'set2c') then
       do i=1,npoin
          rho_k=q(1,i)
          rho_ref=q_ref(1,i)
          rho_exact=q_exact(1,i)

          !Build Conservation, Moisture and Tracers Variables
          qout(1,i)=rho_k
          qout_ref(1,i)=rho_ref
          qout_exact(1,i)=rho_exact
          do m=2,nvar
             qout(m,i)    =rho_k*q(m,i)
             qout_ref(m,i)=rho_ref*q_ref(m,i)
             qout_exact(m,i)=rho_exact*q_exact(m,i)
          end do
       end do
    else if (eqn_set(1:5) == 'set3c') then
       do i=1,npoin

          !Total variables
          rho_k   = q(1,i) !rho
          u_k     = q(2,i) !uvelo
          v_k     = q(3,i) !vvelo
          w_k     = q(4,i) !wvelo
          theta_k = q(5,i) !theta

          p_k     = p00*(rho_k*rgas*theta_k/p00)**gamma
          t_k     = theta_k*(p_k/p00)**(rgas/cp)
          !pi_k    = pi_values(1,i) !this should be pi_k in the icase of initial_conditions
          !t_k     = theta_k*pi_k   !Tempe

          radius = height(i)
          e_k     = cv*t_k + 0.5*(u_k*u_k + v_k*v_k + w_k*w_k) + gravity*radius !E/rho

          !Reference variables:
          rho_ref   = q_ref(1,i)
          u_ref     = q_ref(2,i)
          v_ref     = q_ref(3,i)
          w_ref     = q_ref(4,i)
          theta_ref = q_ref(5,i)

          p_ref     = p00*(rho_ref*rgas*theta_ref/p00)**gamma
          t_ref     = theta_ref*(p_ref/p00)**(rgas/cp)
          !pi_ref    = pi_values(2,i)
          !t_ref     = theta_ref*pi_ref

          e_ref     = cv*t_ref + 0.5*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref) + gravity*radius

          !-----------------------------------------

          !-----------------------------------------
          rho_exact  = q_exact(1,i)
          u_exact    = q_exact(2,i)
          v_exact    = q_exact(3,i)
          w_exact    = q_exact(4,i)
          theta_exact= q_exact(5,i)

          p_exact     = p00*(rho_exact*rgas*theta_exact/p00)**gamma
          t_exact     = theta_exact*(p_exact/p00)**(rgas/cp)
          !pi_exact   = pi_values(3,i)
          !t_exact    = theta_exact*pi_exact

          e_exact    = cv*t_exact + 0.5*(u_exact*u_exact + v_exact*v_exact + w_exact*w_exact) + gravity*radius

          !Build Conservation, Moisture and Tracers Variables
          qout(1,i)    =rho_k
          qout_ref(1,i)=rho_ref
          qout_exact(1,i)=rho_exact
          do m=2,4
             qout(m,i)    =rho_k*q(m,i)
             qout_ref(m,i)=rho_ref*q_ref(m,i)
             qout_exact(m,i)=rho_exact*q_exact(m,i)
          end do
          qout(5,i)=rho_k*e_k
          qout_ref(5,i)=rho_ref*e_ref
          qout_exact(5,i)=rho_exact*e_exact
          do m=6,nvar
             qout(m,i)    =rho_k*q(m,i)
             qout_ref(m,i)=rho_ref*q_ref(m,i)
             qout_exact(m,i)=rho_exact*q_exact(m,i)
          end do
       end do !i=1,npoin
    end if !EQN_SET

  end subroutine mod_convert_variables_set2nc_to_eqnset

  !----------------------------------------------------------------------!
  !>@brief This subroutine converts variables from EQNSET to SET2NC: e.g., used in PRINT_DIAGNOSTICS
  !>@author  Francis X. Giraldo on 1/2014
  !>           Department of Applied Mathematics
  !>           Naval Postgraduate School
  !>           Monterey, CA 93943-5216
  !>----------------------------------------------------------------------!
  subroutine mod_convert_variables_eqnset_to_set2nc(qout,q,q_ref)

    implicit none

    !global arrays
    real, intent(in) :: q(nvar,npoin), q_ref(nvar,npoin)
    real, intent(out):: qout(nvar,npoin)

    !Local arrays
    real :: rho_k, u_k, v_k, w_k, theta_k, e_k, p_k, t_k, pi_k
    real :: rho_ref, u_ref, v_ref, w_ref, theta_ref, e_ref, p_ref, t_ref, pi_ref
    real :: radius, theta0, rfactor
    integer :: i, m

    rfactor=1
    if (is_shallow) rfactor=gravity

    !hack for now, should be current_layer
    if(is_swe_layers) rfactor=rho_layers(1)

    if (eqn_set(1:6) == 'set2nc') then
       qout(1,:)=q(1,:)
       do m=2,4
          qout(m,:)=q(m,:) - q_ref(m,:)
       end do
       qout(5,:)=q(5,:)

       do m=6,nvar
          qout(m,:)=q(m,:) - q_ref(m,:)
       end do

       if(lsalinity) qout(6,:) = q(6,:)

    else if (eqn_set(1:5) == 'set2c') then
       do i=1,npoin

          rho_k   = q(1,i) + q_ref(1,i)
          rho_ref = q_ref(1,i)

          !Variables
!          rho_k   = q(1,i) + q_ref(1,i)
!          u_k     = q(2,i)/rho_k
!          v_k     = q(3,i)/rho_k
!          w_k     = q(4,i)/rho_k
!          theta_k =( q(5,i) + q_ref(5,i))/rho_k
          !-----------------------------------------
!          rho_ref   = q_ref(1,i)
!          u_ref     = q_ref(2,i)/rho_ref
!          v_ref     = q_ref(3,i)/rho_ref
!          w_ref     = q_ref(4,i)/rho_ref
!          theta_ref = q_ref(5,i)/rho_ref

          !Perturbation Variables
          if (is_shallow .and. llimit .and. rho_k <= limit_threshold) then
             rho_k=limit_threshold
             u_k=0
             v_k=0
             w_k=0
             theta_k=0
          else
             u_k     = q(2,i)/rho_k
             v_k     = q(3,i)/rho_k
             w_k     = q(4,i)/rho_k
             theta_k =( q(5,i) + q_ref(5,i))/rho_k
          end if

          !Background Variables
          if (is_shallow .and. llimit .and. abs(rho_ref) < 1e-16) then
             rho_ref=0
             u_ref=0
             v_ref=0
             w_ref=0
             theta_ref=0
          else
             u_ref     = q_ref(2,i)/rho_ref
             v_ref     = q_ref(3,i)/rho_ref
             w_ref     = q_ref(4,i)/rho_ref
             theta_ref = q_ref(5,i)/rho_ref
          end if

          !Store Set2NC variables
          qout(1,i) = (rho_k - rho_ref)/rfactor
          qout(2,i) = u_k - u_ref
          qout(3,i) = v_k - v_ref
          qout(4,i) = w_k - w_ref
          qout(5,i) = theta_k - theta_ref

          !Moisture and Tracers
          do m=6,nvar
             qout(m,i)=q(m,i)/rho_k - q_ref(m,i)/rho_ref
          end do
       end do !i
    else if (eqn_set(1:5) == 'set3c') then
       do i=1,npoin

          !Store Variables
          rho_k = q(1,i) + q_ref(1,i)
          u_k   = q(2,i)/rho_k
          v_k   = q(3,i)/rho_k
          w_k   = q(4,i)/rho_k
          e_k   =( q(5,i) + q_ref(5,i))/rho_k
          !e_k   = q(5,i)/rho_k
          !-----------------------------------------
          rho_ref =q_ref(1,i)
          if(abs(rho_ref) <= 1.0e-16) then
             u_ref = 0.0
             v_ref = 0.0
             w_ref = 0.0
             e_ref = 0.0
             theta_ref = 0.0
          else
             u_ref   =q_ref(2,i)/rho_ref
             v_ref   =q_ref(3,i)/rho_ref
             w_ref   =q_ref(4,i)/rho_ref
             e_ref   =q_ref(5,i)/rho_ref
          end if

          radius = height(i)

          !Compute Pressure, PI, and Temperature
          p_k=(gamma-1)*rho_k*( e_k - 0.5*(u_k*u_k + v_k*v_k + w_k*w_k) - gravity*radius )

          if(icase == 5000 .or. icase == 5001) then
             theta_k   = p_k/(rho_k*rgas)
             p_ref     = 0.0
             theta_ref = 0.0
          else
             if(abs(rho_ref) <= 1.0e-16) then
                theta_k = p00/(rho_k*rgas)*( p_k/p00 )**(1.0/gamma)
                theta_ref = 0.0
                p_ref     = 0.0
             else
                theta_k=p00/(rho_k*rgas)*( p_k/p00 )**(1.0/gamma)
                p_ref=(gamma-1)*rho_ref*( e_ref - 0.5*(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref) - gravity*radius )
                theta_ref=p00/(rho_ref*rgas)*( p_ref/p00 )**(1.0/gamma)
             end if
          end if

          !Store Set 2NC variables
          qout(1,i) = rho_k - rho_ref
          qout(2,i) = u_k - u_ref
          qout(3,i) = v_k - v_ref
          qout(4,i) = w_k - w_ref
          qout(5,i) = theta_k - theta_ref

          !Moisture and Tracers
          if(abs(rho_ref) <= 1.0e-16) then
             do m=6,nvar
                qout(m,:)=q(m,:)/rho_k
             end do
          else
             do m=6,nvar
                qout(m,:)=q(m,:)/rho_k - q_ref(m,:)/rho_ref
             end do
          end if

       end do
    end if !EQN_SET

  end subroutine mod_convert_variables_eqnset_to_set2nc

  !----------------------------------------------------------------------!
  !>@brief This subroutine converts variables from EQNSET to Primitive Perturbation variables: e.g., used in CREATE_RHS_LAPLACIAN and MOD_LIMITER
  !>@author  Francis X. Giraldo on 11/2016
  !>           Department of Applied Mathematics
  !>           Naval Postgraduate School
  !>           Monterey, CA 93943-5216
  !>----------------------------------------------------------------------!
  subroutine mod_convert_variables_eqnset_to_primitive_perturbation(qout,q,q_ref)

    implicit none

    !global arrays
    real, intent(in) :: q(nvar,npoin), q_ref(nvar,npoin)
    real, intent(out):: qout(nvar,npoin)

    !Local arrays
    real :: r_k, rho_k, u_k, v_k, w_k, theta_k
    real :: rho_ref, u_ref, v_ref, w_ref, theta_ref
    integer :: i, m

    if (eqn_set(1:6) == 'set2nc') then
       qout(1,:)=q(1,:)
       do m=2,4
          qout(m,:)=q(m,:) - q_ref(m,:)
       end do
       qout(5,:)=q(5,:)
       do m=6,nvar
          qout(m,:)=q(m,:) - 0*q_ref(m,:) !already subtracted in Initial_Condition
       end do
    else if (eqn_set(1:5) == 'set2c' .or. eqn_set(1:5) == 'set3c') then
       do i=1,npoin
          !Variables
          rho_ref   = q_ref(1,i)
!!          if ( abs(rho_ref) < 1e-16) rho_ref=1e-16 !FXG: SWE 2/13/17
          u_ref     = q_ref(2,i)/rho_ref
          v_ref     = q_ref(3,i)/rho_ref
          w_ref     = q_ref(4,i)/rho_ref
          theta_ref = q_ref(5,i)/rho_ref
          !-----------------------------------------
          rho_k   = q(1,i) !perturbation
          r_k     = rho_k + rho_ref !FXG: SWE 1/28/17
          if (is_shallow .and. llimit .and. r_k <= limit_threshold) then
             r_k=limit_threshold
             u_k=0
             v_k=0
             w_k=0
             theta_k=0
             u_ref=0
             v_ref=0
             w_ref=0
             theta_ref=0
          else
             u_k     = ( q(2,i) + q_ref(2,i) )/r_k
             v_k     = ( q(3,i) + q_ref(3,i) )/r_k
             w_k     = ( q(4,i) + q_ref(4,i) )/r_k
             theta_k = ( q(5,i) + q_ref(5,i) )/r_k
          end if

          !Store Set 2NC variables
          qout(1,i) = rho_k !perturbation
!          qout(1,i) = r_k 
          qout(2,i) = u_k - u_ref
          qout(3,i) = v_k - v_ref
          qout(4,i) = w_k - w_ref
          qout(5,i) = theta_k - theta_ref

          !Moisture and Tracers
          do m=6,nvar
             qout(m,i)=q(m,i)/r_k - 0*q_ref(m,i)/rho_ref !already subtracted in Initial_Condition
          end do
       end do !i
    end if !EQN_SET

  end subroutine mod_convert_variables_eqnset_to_primitive_perturbation

  !----------------------------------------------------------------------!
  !>@brief This subroutine converts variables from Primitive Perturbation variables to EQNSET: e.g., used in CREATE_RHS_LAPLACIAN and MOD_LIMITER
  !>@author  Francis X. Giraldo on 11/2016
  !>           Department of Applied Mathematics
  !>           Naval Postgraduate School
  !>           Monterey, CA 93943-5216
  !>----------------------------------------------------------------------!
  subroutine mod_convert_variables_primitive_perturbation_to_eqnset(qout,q,q_ref)

    implicit none

    !global arrays
    real, intent(in) :: q(nvar,npoin), q_ref(nvar,npoin)
    real, intent(out):: qout(nvar,npoin)

    !Local arrays
    real :: r_k, rho_k, u_k, v_k, w_k, theta_k
    real :: rho_ref, u_ref, v_ref, w_ref, theta_ref
    integer :: i, m

    if (eqn_set(1:6) == 'set2nc') then
       qout(1,:)=q(1,:)
       do m=2,4
          qout(m,:)=q(m,:) + q_ref(m,:)
       end do
       qout(5,:)=q(5,:)
       do m=6,nvar
          qout(m,:)=q(m,:) + 0*q_ref(m,:) !should be a perturbation
       end do
    else if (eqn_set(1:5) == 'set2c' .or. eqn_set(1:5) == 'set3c') then
       do i=1,npoin
          !Variables
          rho_ref   = q_ref(1,i)
          if ( abs(rho_ref) < 1e-16) rho_ref=1e-16 !FXG: SWE 2/13/17
          u_ref     = q_ref(2,i)/rho_ref
          v_ref     = q_ref(3,i)/rho_ref
          w_ref     = q_ref(4,i)/rho_ref
          theta_ref = q_ref(5,i)/rho_ref
          !-----------------------------------------
          rho_k   = q(1,i) !perturbation
          r_k     = rho_k + rho_ref !FXG: SWE 1/28/17
          u_k     = q(2,i) + u_ref
          v_k     = q(3,i) + v_ref
          w_k     = q(4,i) + w_ref
          theta_k = q(5,i) + theta_ref

          !Store Set 2NC variables
          qout(1,i) = rho_k
          qout(2,i) = r_k*u_k - 0*q_ref(2,i) !can take off 0 if moved to Q_REF
          qout(3,i) = r_k*v_k - 0*q_ref(3,i)
          qout(4,i) = r_k*w_k - 0*q_ref(4,i)
          qout(5,i) = r_k*theta_k - q_ref(5,i)

          !Moisture and Tracers
          do m=6,nvar
             qout(m,i)=r_k*q(m,i) + 0*q_ref(m,i) !should be perturbation
          end do
       end do !i
    end if !EQN_SET

  end subroutine mod_convert_variables_primitive_perturbation_to_eqnset


end module mod_convert_variables
