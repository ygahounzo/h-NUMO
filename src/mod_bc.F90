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

    use mod_input, only: icase, decomp_type, sponge_type

    use mod_types, only : r8

    public :: &
        mod_bc_create, &
        mod_bc_create_iperiodic, &
        read_bc, &
        create_bc_list,&
        aa, &
        bb, &
        errmask, &
        lsponge, &
        npoin_errmask, &
        normals, &
        bc_list, bc_count, & !for identifying boundary points not on boundary element faces (unstructured meshes)
        vc_el_type, &
        mod_bc_init_restoring, mod_bc_apply_restoring, &
        robin_bc_beta 
    

    private
    !-----------------------------------------------------------------------
    real,    dimension(:),      allocatable :: aa, bb, errmask
    real,    dimension(:),      allocatable :: robin_bc_beta !values for beta in n\nabla q = \alpha q + \beta
    real,    dimension(:,:),    allocatable :: normals
    integer, dimension(:),      allocatable :: ipoin_bound, ip_bound_bc
    integer, dimension(:,:),    allocatable :: bc_list
    integer, dimension(4) :: bc_count
    integer, dimension(:),      allocatable :: vc_el_type
    integer npoin_bound
    integer npoin_errmask, num_bound_bc
    logical lsponge, ldirichlet
  !-----------------------------------------------------------------------  

contains

    !-----------------------------------------------------------------------
    subroutine mod_bc_create()

        implicit none

        integer AllocateStatus
        integer ierr
!
!         ! Allocate Memory for BC data structures
!         if(allocated(aa)) then
!             deallocate(aa, bb, errmask)
!         endif
!         allocate(aa(npoin), bb(npoin), errmask(npoin), stat=AllocateStatus )
!         if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_BC **"
!
!         !Construct Sponge Layer for NRBCs
!         bb = 0.0
!         aa = 1.0
!         errmask = 0.0
        
        ! Allocate Memory for BC data structures
        if(allocated(aa)) then
           deallocate(normals, aa, bb, errmask, ipoin_bound, robin_bc_beta)
        endif
        allocate(normals(3,npoin), aa(npoin), bb(npoin), &
           errmask(npoin), ipoin_bound(npoin), robin_bc_beta(npoin), stat=AllocateStatus )
        if (AllocateStatus /= 0) stop "** Not Enough Memory - Mod_BC **"

        !Initialize sponge arrays to zero:
        bb = 0.0
        aa = 1.0
        errmask = 0.0

        !Construct PMATRIX for NFBCs
!        call create_pmatrix(pmatrix,normals,ipoin_bound,npoin_bound,ldirichlet)
        call create_nfbc_vector(normals,ipoin_bound,npoin_bound,ldirichlet)

        call create_sponge(aa,bb,errmask,npoin_errmask,lsponge)

        call create_bc_list()
        
    end subroutine mod_bc_create
  
    !----------------------------------------------------------------------!
    !@brief This subroutine computes the Sponge (NRBC/ABC) Layer.
    !@author  Francis X. Giraldo on 11/2009
    !           Department of Applied Mathematics
    !           Naval Postgraduate School
    !           Monterey, CA 93943-5216
    !----------------------------------------------------------------------!
    subroutine create_sponge(alpha,beta,errmask,npoin_errmask,lsponge)

        use mod_input, only: geometry_type

        use mod_global_grid, only: xmin, xmax, ymin, ymax, zmax, iboundary

        use mod_input, only: dt, bc_xscale, bc_yscale, bc_zscale, bc_tscale, &
            sponge_top_coe, sponge_lateralx_coe, sponge_lateraly_coe, lsommerfeld, &
            x_boundary, y_boundary, z_boundary,  nelx, nely, nelz

        use mod_constants, only : pi, earth_radius

        !>    use mod_basis, only: ngl, nglx, ngly, nglz

        implicit none

        real(kind=r8), intent(inout) :: alpha(:), beta(:), errmask(:)
        integer,       intent(out) :: npoin_errmask
        logical,       intent(out) :: lsponge

        real(kind=r8) :: pio2, zs, xsl, xsr, ysl, ysr
        real(kind=r8) :: x, y, z, r, phi, lambda
        real(kind=r8) :: csidex, csidexl, csidexr, csidey, csideyl, csideyr
        real(kind=r8) :: ccorner1, ccorner2, ccorner3, ccorner4
        real(kind=r8) :: ccorner5, ccorner6, ccorner7, ccorner8
        real(kind=r8) :: ct, ctop, cs, csx, csy, csz, cside, dside, dblx, dbly, dblz, dbt
        real(kind=r8) :: alpha_coe, beta_coe, zdcorrect, zd, ztop, dsx, dsy, dsz, dxs, dys, dzs, zid, xid, yid
        real(kind=r8) :: ampx,ampy,ampz
        real(kind=r8) :: xl, xr, dbl
        real(kind=r8) :: abstaud, abstaudx, abstaudy, xbound
        integer :: i
        logical :: lcube

        pio2=0.5*pi

        if (geometry_type(1:4) == 'cube') then
            lcube=.true.
        else
            lcube=.false.
        end if

        !Initialize Sponges
        alpha(:)=1.0; beta(:)=0.0; errmask(:)=0.0; npoin_errmask=npoin

        !Check if sponges need to be used
        if(lcube) then
            lsponge = any(iboundary==6) .or. any(iboundary==8)
        else
            lsponge = any(iboundary(1:2)==6) .or. any(iboundary(1:2)==8)
        end if

        ! Return if there are no sponges
        if (.NOT. lsponge) return

        zs  = zmax - bc_zscale
        if(lcube) then
            xsl = xmin + bc_xscale
            xsr = xmax - bc_xscale
            ysl = ymin + bc_yscale
            ysr = ymax - bc_yscale
        end if

        if(sponge_type == 'default' .or. sponge_type == 'sine') then
       
            !
            !Sponge coefficients from input or default:
            !
            !X-lateral coefficient
            if(sponge_lateralx_coe < 0) then
                csx = 1.0
            else
                csx = sponge_lateralx_coe
            end if

            !Y-lateral coefficient
            if(sponge_lateraly_coe < 0) then
                csy = 1.0
            else
                csy = sponge_lateraly_coe
            end if

            !Z-top coefficient
            if(sponge_top_coe < 0) then
                ct = 1.0
            else
                ct = sponge_top_coe
            end if

       
            !loop thru the elements
            do i=1,npoin
                x=coord(1,i) ; y=coord(2,i) ; z=coord(3,i)

                if(lcube) then
                    !X-Lateral Sponge
                    if (x <= xsl .and. iboundary(5).eq.6) then
                        csidexl = csx*(sin(pio2*(x-xsl)/(xmin-xsl)) )**4
                       !beta(i) = beta(i) + csx*(sin(pio2*(x-xsl)/(xmin-xsl)) )**4
                    end if

                    if (x >= xsr .and. iboundary(6).eq.6) then
                        csidexr = csx*(sin(pio2*(x-xsr)/(xmax-xsr)) )**4
                       !beta(i) = beta(i) + csx*(sin(pio2*(x-xsr)/(xmax-xsr)) )**4
                    end if

                    !Y-Lateral Sponge
                    if (y <= ysl .and. iboundary(3).eq.6) then
                        csideyl = csy*(sin(pio2*(y-ysl)/(ymin-ysl)) )**4
                       !beta(i) = beta(i) + csy*(sin(pio2*(y-ysl)/(ymin-ysl)) )**4
                    end if

                    if (y >= ysr .and. iboundary(4).eq.6) then
                        csideyr = csy*(sin(pio2*(y-ysr)/(ymax-ysr)) )**4
                       !beta(i) = beta(i) + csy*(sin(pio2*(y-ysr)/(ymax-ysr)) )**4
                    end if
                else
                    call cart2sph(x,y,z,r,phi,lambda)
                    z = r - earth_radius
                end if
          
                !Top Sponge
                if (z >= zs .and. iboundary(2).eq.6) then
                    ctop = ct*(sin(pio2*(z-zs)/(zmax-zs)))**4 !ct*(dt/bc_tscale)*(sin(pio2*(z-zs)/(zmax-zs)))**4
                   !beta(i) = beta(i) + ct*(sin(pio2*(z-zs)/(zmax-zs)))**4 !ct*(dt/bc_tscale)*(sin(pio2*(z-zs)/(zmax-zs)))**4
                end if
                beta(i) = 1.0 - (1.0 - ctop)*(1.0-csidexl)*(1.0-csidexr)*(1.0-csideyl)*(1.0-csideyr)

                !Store Values
                beta(i) = min(beta(i), 1.0)
                alpha(i) = 1.0-beta(i)

                !Check if the point is in the sponge layer
                if(beta(i) > 0.0) then
                    errmask(i) = 1.0
                    npoin_errmask=npoin_errmask - 1
                endif

            end do !i

        else if(sponge_type == 'simple' ) then
            !
            !Sponge coefficients from input or default:
            !
            !X-lateral coefficient
            if(sponge_lateralx_coe < 0) then
                csx = 1.0
            else
                csx = sponge_lateralx_coe
            end if

            !Y-lateral coefficient
            if(sponge_lateraly_coe < 0) then
                csy = 1.0
            else
                csy = sponge_lateraly_coe
            end if

            !Z-top coefficient
            if(sponge_top_coe < 0) then
                ct = 1.0
            else
                ct = sponge_top_coe
            end if
       
            !
            ! loop thru the elements
            !
            do i=1,npoin
                z=coord(3,i)
          
                !Top Sponge
                if (z >= zs .and. iboundary(2).eq.6) then
                    ctop = ct*( (z-zs)/(zmax-zs) )**4
                end if

                !Store Values
                beta(i) = ctop
                beta(i) = min(beta(i), 1.0)
                alpha(i) = 1.0 - beta(i)

                !Check if the point is in the sponge layer
                if ( z >= zs) THEN
                    errmask(i) = 0.0
                else
                    errmask(i)    = 1.0
                    npoin_errmask = npoin_errmask - 1
                endif

            end do !i

       
        else if(sponge_type == 'original') then

            print*,' Sponge type: Original, from numa2dCGDG_AMR'
            beta  = 0.0
            alpha = 0.0

            !
            ! Sponge of Lilly and Klemp 1978
            !
            zd       =  17700.0!11500.0 ! m
            alpha_coe=0.0057

            !Sponge apmplitudes: ct (top), cs (side):
            ct=1.0
            cs=1.0

            !Constants
            ztop   = zmax
!            xbound = dxs
            xbound=(xmax-xmin)/8
            
            xl     = xmin + xbound
            xr     = xmax - xbound
       
            dsx       = (xmax-xmin)/(real(nelx)*(nglx-1)) ! equivalent grid spacing
            dsy       = (ymax-ymin)/(real(nely)*(ngly-1)) ! equivalent grid spacign
            dsz       = (ztop)/(real(nelz)*(nglz-1)) ! equivalent grid spacing
            zdcorrect = zd - 8*dsz

            !loop thru the elements
            npoin_errmask=0
            do i=1,npoin

                cside=0.0; ctop=0.0
          
                x=coord(1,i)
                y=coord(2,i)
                z=coord(3,i)

                !
                ! boundary damping
                !
                dbl   = min(x - xmin, xmax - x) ! distance from the boundary. xs in Restelli's thesis
                beta_coe =  1.0 - tanh(dbl/(3.0 * dsx))
                cside = cs*beta_coe

                !
                ! top damping
                ! first layer: damp lee waves
                !
                if(z <= zdcorrect) then
                    ctop = 0.0
                else
                    xid = (z-zd)/(ztop-zd) ! normalized coordinate
                    if(xid < 0.5) then
                        abstaud = 0.5*alpha_coe*(1.0 - cos(xid*pi))
                    else
                        abstaud = 0.5*alpha_coe*(1.0 + (xid - 0.5)*pi)
                    endif
                    ctop = ct * (dt*abstaud)/(1.0 + 0.5*dt*abstaud) ! Sasa
                endif
          
                !
                ! second layer: damp short waves
                !
                dbt = ztop - z! distance from the boundary
                beta_coe = 1.0 - tanh(dbt/(2.0 * dsz)) !sasa
                ctop = ctop + ct*beta_coe

                ! normalization
                ctop = min(ctop, 1.0)
                    
                ! boundary damping
                !Store Values
                bb(i) = min(ctop+cside,1.0)
                aa(i) = 1.0 - bb(i)

                !Check if the point is in the sponge layer
                !if ( (z >= zs) .or. (x <= xsl) .or. (x >= xsr) ) then
                if((dbl < 16.0*dsx) .or. (z > 0.9*zdcorrect)) then
                   errmask(i) = 0.0
                else
                   errmask(i) = 1.0
                   npoin_errmask=npoin_errmask + 1
                endif

            end do !i

            !Rewrite alpha and beta
            do i=1,npoin
               alpha(i)=aa(i)
               beta(i)=bb(i)
            end do

         else if(sponge_type == 'klemp'      .or. &
            sponge_type == 'LK78'       .or. &
            sponge_type == 'KL78'       .or. &
            sponge_type == 'lillyKlemp' .or. &
            sponge_type == 'klempLilly') then

            print*,' Sponge type: Lilly and Klemp 1978, added by Simone Marras'

            !
            ! Sponge of Lilly and Klemp 1978
            !
            zd       =  bc_zscale  !11500.0 ! m
            alpha(1) = 0.0057      ! s^-1
       
            beta  = 0.0
            alpha = 0.0
       
            !Sponge apmplitudes: ct (top), cs (side):
            ct=1.0
            cs=1.0

            !Constants
            ztop   = zmax
            xbound = dxs

            xl     = xmin + xbound
            xr     = xmax - xbound
       
            dsx       = (xmax-xmin)/(real(nelx)*(nglx-1)) ! equivalent grid spacing
            dsy       = (ymax-ymin)/(real(nely)*(ngly-1)) ! equivalent grid spacign
            dsz       = (ztop)/(real(nelz)*(nglz-1))      ! equivalent grid spacing
            zdcorrect = zd !- 8*dsz

            !loop thru the elements
            npoin_errmask=0
            do i=1,npoin

                cside=0.0; ctop=0.0
          
                x=coord(1,i)
                y=coord(2,i)
                z=coord(3,i)

                !
                ! boundary damping
                !
                dbl   = min(x - xmin, xmax - x) ! distance from the boundary. xs in Restelli's thesis
                beta(1)  =  1.0 - tanh(dbl/(3.0 * dsx))
                !cside = cs*beta(1)
                cside = 0.0

                !
                ! top damping
                ! first layer: damp lee waves
                !
                if(z <= zd) then
                    ctop = 0.0
                else
                    xid = (z-zd)/(ztop-zd) ! normalized coordinate
                    if(xid < 0.5) then
                        abstaud = 0.5*alpha(1)*(1.0 - cos(xid*pi))

                    else
                        abstaud = 0.5*alpha(1)*(1.0 + (xid - 0.5)*pi)
                    endif

                    ctop = ct * (dt*abstaud)/(1.0 + 0.5*dt*abstaud) ! Sasa
                   !ctop = ct*abstaud ! Simone
                endif
          
                !
                ! second layer: damp short waves
                !
                dbt = zmax - z + 100.0 ! distance from the boundary
                !beta(1) = (1.0 - tanh(dbt/dsz))/tanh((dbt)/dsz) !sm
                beta(1) = 1.0 - tanh(dbt/(2.0 * dsz)) !sasa

                ! normalization
                ctop = min(ctop, 1.0)
                    
                ! boundary damping
                !Store Values
                bb(i) = ctop !+ cside
                bb(i) = min(bb(i), 1.0)
                aa(i) = 1.0 - bb(i)

                !Check if the point is in the sponge layer
                !if ( (z >= zs) .or. (x <= xsl) .or. (x >= xsr) ) then
                if ( (z > 0.9*zdcorrect) .or. (dbl < 16.0*dsx)) then
                    errmask(i) = 0.0
                else
                    errmask(i) = 1.0
                    npoin_errmask=npoin_errmask + 1
                endif

            end do !i

        else if(sponge_type(1:7) == 'testing') then

            !Constants
            ztop = zmax

            !
            ! Sponge of Durran and Klemp 1983
            !
            zd = ztop - bc_zscale   ! m

            !
            !Sponge coefficients from input or default:
            !
            !X-lateral coefficient
            if(sponge_lateralx_coe < 0) then
                csx = 1.0
            else
                csx = sponge_lateralx_coe
            end if

            !Y-lateral coefficient
            if(sponge_lateraly_coe < 0) then
                csy = 1.0
            else
                csy = sponge_lateraly_coe
            end if

            !Z-top coefficient
            if(sponge_top_coe < 0) then
                ct = 1.0
            else
                ct = sponge_top_coe
            end if
       
            dsx    = bc_xscale
            dsy    = bc_yscale
            dsz    =  660.0
            zdcorrect = zd
            alpha_coe = 0.5

            !loop thru the elements
            npoin_errmask=0
            do i=1,npoin

                csidex=0.0; csidey=0.0; ctop=0.0

                x=coord(1,i)
                y=coord(2,i)
                z=coord(3,i)

                !
                ! top damping
                ! first layer: damp lee waves
                !
                ctop = 0.0
                if(z >= zd  .and. iboundary(2) == 6 .and. x >= xsl) then
                    zid = (z-zd)/(ztop-zd) ! normalized coordinate
                    if(zid >= 0.0 .and. zid <= 0.5) then
                        abstaud = alpha_coe*(1.0 - cos(zid*pi))

                    else
                        abstaud = alpha_coe*(1.0 + (zid - 0.5)*pi)

                    endif
                    ctop = ct*abstaud
                endif
          
                !X-Lateral Sponge
                if (x <= xsl .and. z < zd .and. iboundary(5) == 6) then
                    csidexl = csx*( (xsl - x)/(xsl - xmin) )**4
                else if (x >= xsr .and.  z < zd .and. iboundary(6) == 6) then
                    csidexr = csx*( (x - xsr)/(xmax - xsr) )**4
                end if
                if(lsommerfeld) then
                    csidexl = 0.0
                    csidexr = 0.0
                end if

                !Y-Lateral Sponge
                if (y <= ysl .and. z < zd .and. iboundary(3) == 6) then
                    csideyl = csy*( (ysl - y)/(ysl - ymin) )**4
                else if (y >= ysr .and. z < zd .and. iboundary(4) == 6) then
                    csideyr = csy*( (y - ysr)/(ymax - ysr) )**4
                end if
                if(lsommerfeld) then
                    csideyl = 0.0
                    csideyr = 0.0
                end if
                !
                ! second layer: damp short waves
                !
                !dbt = ztop - z ! distance from the boundary
                !if(dbt == 0.0) then
                !   beta_coe = 1.0
                !else
                !   beta_coe = (1.0 - tanh(dbt/dsz))/tanh(dbt/dsz)
                !end if
                !ctop = ctop + ct*beta_coe
          
                ! boundary damping
                !Store Values
                if(x <= xsl .and. z >= zd .and. iboundary(5) == 6) then
                    ccorner4 = max(ctop, csidexl)
                end if

                beta(i) = ctop + csidexl + csidexr + csideyl + csideyr + ccorner4
                beta(i) = min(beta(i), 1.0)
                alpha(i) = 1.0 - beta(i)

                !Check if the point is in the sponge layer
                if ( (z >= zs) .or. (x <= xsl) .or. (x >= xsr) .or. (y <= ysl) .or. (y >= ysr) ) then
                    errmask(i) = 0.0
                else
                    errmask(i) = 1.0
                    npoin_errmask=npoin_errmask + 1
                endif

            end do !ipoin

        else if(sponge_type(1:5) == 'mixed') then

            !Constants
            ztop = zmax

            !
            ! Sponge of Durran and Klemp 1983
            !
            zd = ztop - bc_zscale   ! m

            !
            !Sponge coefficients from input or default:
            !
            !Top coefficient
            if(sponge_top_coe < 0) then
                ct = 0.5
            else
                ct = sponge_top_coe
            end if

            !X-lateral coefficient
            if(sponge_lateralx_coe < 0) then
                csx = 0.0
            else
                csx = sponge_lateralx_coe
            end if

            !Y-lateral coefficient
            if(sponge_lateraly_coe < 0) then
                csy = 0.0
            else
                csy = sponge_lateraly_coe
            end if

            dsx    = bc_xscale
            dsy    = bc_yscale
            dsz    =  660.0
            zdcorrect = zd
            alpha_coe = 0.5

            !loop thru the elements
            npoin_errmask=0
            do i=1,npoin

                csidex = 0.0; csidexl=0.0; csidexr=0.0
                csidey=0.0;   csideyl=0.0; csideyr=0.0
                ctop=0.0

                x=coord(1,i)
                y=coord(2,i)
                z=coord(3,i)

                !
                ! top damping
                ! first layer: damp lee waves
                !
                if(z .le. zd) then
                    ctop = 0
                else
                    zid = (z-zd)/(ztop-zd) ! normalized coordinate
                    if(zid >= 0.0 .and. zid <= 0.5) then
                        abstaud = 0.5*alpha_coe*(1.0 - cos(zid*pi))

                    else
                        abstaud = 0.5*alpha_coe*(1.0 + (zid - 0.5)*pi)

                    endif
                    ctop = ct*abstaud
                endif

                !X-Lateral Sponge
                if (x <= xsl) then
                    csidexl = csx*( (xsl - x)/(xsl - xmin) )**4
                else if (x >= xsr) then
                    csidexr = csx*( (x - xsr)/(xmax - xsr) )**4
                end if
                if(lsommerfeld) csidexl = 0.0

                !Y-Lateral Sponge
                if (y <= ysl) then
                    csidey = csy*( (ysl - y)/(ysl - ymin) )**4
                else if (y >= ysr) then
                    csidey = csy*( (y - ysr)/(ymax - ysr) )**4
                end if

                !
                ! second layer: damp short waves
                !
                !dbt = ztop - z ! distance from the boundary
                !if(dbt == 0.0) then
                !   beta_coe = 1.0
                !else
                !   beta_coe = (1.0 - tanh(dbt/dsz))/tanh(dbt/dsz)
                !end if
                !ctop = ctop + ct*beta_coe

                ! boundary damping
                !Store Values
                beta(i) = ctop + csidexl + csidexr +csidey
                beta(i) = min(beta(i), 1.0)
                alpha(i) = 1.0 - beta(i)

                !Check if the point is in the sponge layer
                if ( (z >= zs) .or. (x <= xsl) .or. (x >= xsr) .or. (y <= ysl) .or. (y >= ysr) ) then
                    errmask(i) = 0.0
                else
                    errmask(i) = 1.0
                    npoin_errmask=npoin_errmask + 1
                endif

            end do !ipoin

        end if

    end subroutine create_sponge

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


    subroutine create_bc_list()!bc_list,bc_count)

      use mod_grid, only: npoin, mod_grid_get_face_ngl, nface, face, intma, coord
      
      use mod_face, only: imapl
      
      use mod_input, only: lfree_slip_exact

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
                  else if (bcflag(k) == -4 .and. .not.lfree_slip_exact) then !no-flux (velocity-neumann, pressure-dirichlet)
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
                  else if (bcflag(k) == -4 .and. .not.lfree_slip_exact) then !no-flux (velocity-neumann, pressure-dirichlet)
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

    !Subroutine to initialize sponge arrays read from gmsh
    subroutine mod_bc_init_restoring()

      use mod_input, only: lrestoring_sponge

      implicit none

      if(allocated(vc_el_type)) deallocate(vc_el_type)

      if(lrestoring_sponge) allocate(vc_el_type(nelem))

    end subroutine mod_bc_init_restoring

!------------------------------------------------------------------------------------
!>@brief Integrate sponge term to contribute to rhs
!>@author Michal A. Kopera
!------------------------------------------------------------------------------------

    subroutine mod_bc_apply_restoring(rhs,q)

      use mod_grid, only: intma

      use mod_initial, only: nvar, q_ref, q_exact, q_init

      use mod_input, only: restoring_time

      use mod_metrics, only: jac

      implicit none

      real, intent(inout) :: rhs(nvar,npoin)
      real, intent(in) :: q(nvar,npoin)

      real :: wq
      integer :: e, vct, i, j, k, m, ip
      integer, dimension(nvar) :: vcflag

      do e = 1,nelem
         vct = vc_el_type(e)
         vcflag = 0

         if(vct>0) then

            vcflag(1) = 0 !no sponge for pressure
            vcflag(2) = mod(vct,10) !velocity
            vcflag(3) = mod(vct,10) !the same
            vcflag(4) = mod(vct,10) !for all components
            vcflag(5) = mod(vct/10,10) !temperature
            vcflag(6) = mod(vct/100,10) !salinity

            do k=1,nglz
               do j=1,ngly
                  do i=1,nglx
                     ip = intma(i,j,k,e)
                     wq = jac(i,j,k,e)
                     do m=1,nvar
                        rhs(m,ip) = rhs(m,ip) + vcflag(m) * wq * 1.0/restoring_time * (q_init(m,ip)-q(m,ip))
                     end do

                  end do
               end do
            end do

         end if

      end do


    end subroutine mod_bc_apply_restoring

end module mod_bc
