! Optional internal-wave LLF penalty for the split layer fluxes.
! This is an experimental resting-column speed estimate, not an entropy-
! stability or arbitrary-shear hyperbolicity guarantee for the split scheme.
module mod_bcl_dissipation
    implicit none
    private
    public :: internal_wave_speed, internal_wave_penalty
contains
    real function internal_wave_speed(n, alpha, state, gravity) result(speed)
        !$acc routine seq
        integer, intent(in) :: n
        real, intent(in) :: alpha(n), state(4,n), gravity
        real :: h, rho, total_h, total_mass, pairs, previous_rho, density_step
        real :: h1, h2, reduced, discriminant
        integer :: k

        speed = 0.0
        if (n < 2) return
        if (n == 2) then
            h1 = max(0.0, alpha(1)*state(1,1)/gravity)
            h2 = max(0.0, alpha(2)*state(1,2)/gravity)
            total_h = h1+h2
            if (total_h <= 0.0) return
            reduced = max(0.0, (alpha(1)-alpha(2))/alpha(1))
            discriminant = max(0.0, total_h**2-4.0*reduced*h1*h2)
            speed = sqrt(2.0*gravity*reduced*h1*h2/(total_h+sqrt(discriminant)))
            return
        end if

        ! Bound c_internal,max^2 by trace(S)-R(x), where
        ! S_ij=g*sqrt(h_i*h_j)*min(rho_i,rho_j)/sqrt(rho_i*rho_j)
        ! is symmetric positive semidefinite and x_i=sqrt(rho_i*h_i).
        ! The result is g*sum_{i<j} h_i*h_j*(rho_j-rho_i)/sum rho_i*h_i.
        ! Accumulate density differences without subtracting large sums.
        total_h = 0.0
        total_mass = 0.0
        pairs = 0.0
        density_step = 0.0
        previous_rho = 1.0/alpha(1)
        do k = 1, n
            rho = 1.0/alpha(k)
            h = max(0.0, alpha(k)*state(1,k)/gravity)
            density_step = density_step + (rho-previous_rho)*total_h
            pairs = pairs + h*density_step
            total_h = total_h+h
            total_mass = total_mass+rho*h
            previous_rho = rho
        end do
        if (total_mass > 0.0) speed = sqrt(gravity*max(0.0,pairs)/total_mass)
    end function internal_wave_speed

    subroutine internal_wave_penalty(n, alpha, left, right, factor, dry_cutoff, gravity, penalty)
        !$acc routine seq
        integer, intent(in) :: n
        real, intent(in) :: alpha(n), left(4,n), right(4,n), factor, dry_cutoff, gravity
        real, intent(out) :: penalty(4,n)
        real :: diss, hmin
        integer :: k, iv

        penalty = 0.0
        if (factor <= 0.0 .or. n < 2) return
        ! Keep the existing wet/dry treatment: disable the added penalty for
        ! a face if either column has a layer at/below the dry threshold.
        do k = 1, n
            hmin = alpha(k)*min(left(1,k),right(1,k))/gravity
            if (hmin <= dry_cutoff*(1.0+64.0*epsilon(hmin))) return
        end do
        diss = 0.5*factor*max(internal_wave_speed(n,alpha,left,gravity), &
                              internal_wave_speed(n,alpha,right,gravity))
        do k = 1, n
            do iv = 1, 4
                penalty(iv,k) = diss*(left(iv,k)-right(iv,k))
            end do
        end do
    end subroutine internal_wave_penalty
end module mod_bcl_dissipation
