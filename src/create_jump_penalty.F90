!-------------------------------------------------------------
!@>brief Calculate jump penalty flux
!-------------------------------------------------------------
subroutine compute_jump_penalty_weight()

    use mod_face, only: penalty, imapl, imapr, jac_face, normal_vector
    
    use mod_grid, only: nface, face, face_type, mod_grid_get_face_ngl
    
    use mod_basis, only: ngl, is_2d, dpsix, dpsiy, dpsiz, psix, psiy, psiz

    use mod_metrics, only: &
        ksi_x, ksi_y, ksi_z, &
        eta_x, eta_y, eta_z, &
        zeta_x, zeta_y, zeta_z, jac

    implicit none
    
    integer :: iface, itype, ier, er, iel, el, ilocl
    integer :: i, j, k, l, il, jl, kl, ii, jj, kk
    integer :: ngl_i, ngl_j, plane_ij
    real :: mul, mur, pconst, pen, nx, ny, nz
    real :: h_e, h_n, h_c, psi_e, psi_n, psi_c, psi_x, psi_y, psi_z

    !Shahbazi penalty constant
    if(is_2d) then
        pconst = (ngl + 1) * (ngl + 2) / 2.0
    else
        pconst = (ngl + 1) * (ngl + 3) / 3.0
    endif
        
    !Initialize
    penalty=0
        
    do iface=1,nface
        itype = face_type(iface)
        
        !left & right elements
        ilocl=face(5,iface)
        er=face(8,iface)
        el=face(7,iface)

        !get face ngls
        call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

        !Loop through Quadrature Points
        do j=1,ngl_j
            do i=1,ngl_i

                !normals
                nx=normal_vector(1,i,j,iface)
                ny=normal_vector(2,i,j,iface)
                nz=normal_vector(3,i,j,iface)
                
                !Store Left Pointers to Quadrature Points
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                
                !shahbazi penalty
                pen = pconst * jac_face(i,j,iface)/jac(il,jl,kl,el)
     
                !Compute Computational Derivatives
                psi_e=0; psi_n=0; psi_c=0
                do l=1,ngl_j
                    do k=1,ngl_i
                        ii=imapl(1,k,l,iface)
                        jj=imapl(2,k,l,iface)
                        kk=imapl(3,k,l,iface)
                        h_e=dpsix(ii,il)*psiy(jj,jl)*psiz(kk,kl)
                        h_n=psix(ii,il)*dpsiy(jj,jl)*psiz(kk,kl)
                        h_c=psix(ii,il)*psiy(jj,jl)*dpsiz(kk,kl)
                        psi_e=psi_e + h_e
                        psi_n=psi_n + h_n
                        psi_c=psi_c + h_c
                    end do !i
                end do !j

                !Compute Physical Derivatives
                psi_x=psi_e*ksi_x(il,jl,kl,el) + psi_n*eta_x(il,jl,kl,el) + psi_c*zeta_x(il,jl,kl,el)
                psi_y=psi_e*ksi_y(il,jl,kl,el) + psi_n*eta_y(il,jl,kl,el) + psi_c*zeta_y(il,jl,kl,el)
                psi_z=psi_e*ksi_z(il,jl,kl,el) + psi_n*eta_z(il,jl,kl,el) + psi_c*zeta_z(il,jl,kl,el)
                pen = pen + (nx*psi_x + ny*psi_y + nz*psi_z)
                
                !store penalty weight
                penalty(i,j,iface) = pen
            
            end do !i
        end do !j
        
    end do
    
end subroutine

!-------------------------------------------------------------
!@>brief Calculate jump penalty flux for gradient
!-------------------------------------------------------------
subroutine compute_gradient_penalty_flux(rhs,q,lpen)

    use mod_basis, only: ngl

    use mod_face, only: normal_vector, jac_face, imapl, imapr, penalty

    use mod_grid, only: npoin, intma, face_type, nface, face, mod_grid_get_face_ngl, is_cg_coupled

    use mod_input, only: icase, space_method, cgdg_method, visc, is_non_conforming_flg

    implicit none

    !global arrays
    real, intent(inout) :: rhs(3,npoin)
    real, intent(in)    :: q(npoin)
    logical, intent(in) :: lpen

    !local arrays
    real, dimension(ngl,ngl) :: ql, qr

    !local variables
    real unl, nx, ny, nz, wq, qt(3)
    real :: pen = 0
    integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er, itype
    integer iel, ier, ilocl, ilocr,  ifa, ip
    integer ngl_i, ngl_j, plane_ij, bcflag, start
    logical compute_boundary_flux, compute_mixed_flux

    if(space_method == 'cgc' .or. (space_method == 'cgd' .and. cgdg_method == 'separate')) return
    compute_boundary_flux = (space_method == 'cgd' .and. cgdg_method == 'unified')
    compute_mixed_flux    = (space_method == 'cgd' .and. cgdg_method == 'mixed')
   
    !Construct FVM-type Operators
    do iface=1,nface
        itype = face_type(iface)   !Type of face

        !Skip boundary faces
        if (itype == 2 .or. itype == 21 .or. itype == 12) then
            cycle
        end if

        !Skip periodic
        ier=face(8,iface)
        er=ier
        if (ier == -6) then
            cycle
        end if

        !Compute boundary flux
        if(compute_boundary_flux .and. ier > 0) then
            cycle
        end if

        !Skip non-conforming faces
        if(is_non_conforming_flg>0) then
            if(face(9,iface)>0) cycle
        end if

        !-------------------------------------
        !Store Left and Right Side Variables
        !-------------------------------------
        ilocl=face(5,iface)
        iel=face(7,iface)
        el=iel
        call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

        !----------------------------Left Element
        do j=1,ngl_j
            do i=1,ngl_i
                !Get Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                ip=intma(il,jl,kl,iel)

                ql(i,j)=q(ip)

            end do !i
        end do !j

        !----------------------------Right Element
        if (ier > 0 ) then
            ilocr=face(6,iface)
            do j=1,ngl_j
                do i=1,ngl_i

                    !Get Pointers
                    ir=imapr(1,i,j,iface)
                    jr=imapr(2,i,j,iface)
                    kr=imapr(3,i,j,iface)
                    ip=intma(ir,jr,kr,ier)

                    qr(i,j)=q(ip)
                end do !i
            end do !j

        !----------------------------Boundary element
        else if (ier <= 0) then
            bcflag = mod(ier,10)
            
            !Use Left Element IMAP and Metric Terms
            er=el
            ilocr=ilocl

            !default values
            qr=ql
            
        end if

        !---------------------------------
        !  Do Gauss-Lobatto Integration
        !---------------------------------
        !--------------Left Side------------------!
        do j=1,ngl_j
            do i=1,ngl_i
                !Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                wq=jac_face(i,j,iface)
                ip=intma(il,jl,kl,iel)
                if(.not.(compute_mixed_flux .and. is_cg_coupled(il,jl,kl,iel) .ne. 0 )) then
                    nx=normal_vector(1,i,j,iface)
                    ny=normal_vector(2,i,j,iface)
                    nz=normal_vector(3,i,j,iface)
                    
                    if(lpen) pen = penalty(i,j,iface)
                    rhs(1,ip)=rhs(1,ip) + wq * nx * (pen + 0.5) * (qr(i,j) - ql(i,j))
                    rhs(2,ip)=rhs(2,ip) + wq * ny * (pen + 0.5) * (qr(i,j) - ql(i,j))
                    rhs(3,ip)=rhs(3,ip) + wq * nz * (pen + 0.5) * (qr(i,j) - ql(i,j))
                endif

            end do !i
        end do !j

        !--------------Right Side------------------!
        if (ier > 0) then
            do j=1,ngl_j
                do i=1,ngl_i
                    !Pointers
                    ir=imapr(1,i,j,iface)
                    jr=imapr(2,i,j,iface)
                    kr=imapr(3,i,j,iface)
                    wq=jac_face(i,j,iface)
                    ip=intma(ir,jr,kr,ier)
                    if(.not.(compute_mixed_flux .and. is_cg_coupled(ir,jr,kr,ier) .ne. 0 )) then
                        nx=-normal_vector(1,i,j,iface)
                        ny=-normal_vector(2,i,j,iface)
                        nz=-normal_vector(3,i,j,iface)
                    
                        if(lpen) pen = penalty(i,j,iface)
                        rhs(1,ip)=rhs(1,ip) + wq * nx * (pen + 0.5) * (ql(i,j) - qr(i,j))
                        rhs(2,ip)=rhs(2,ip) + wq * ny * (pen + 0.5) * (ql(i,j) - qr(i,j))
                        rhs(3,ip)=rhs(3,ip) + wq * nz * (pen + 0.5) * (ql(i,j) - qr(i,j))
                    endif
                end do !i
            end do !j
        end if !ier

    end do !iface
    
end subroutine

!-------------------------------------------------------------
!@>brief Calculate jump penalty flux for gradient on nc faces
!-------------------------------------------------------------
subroutine compute_gradient_penalty_flux_nc(rhs,q,lpen)

    use mod_basis, only: ngl, is_2d, FACE_CHILDREN

    use mod_face, only: normal_vector, jac_face, imapl, imapr, face_nc, nface_nc, penalty

    use mod_grid, only: npoin, intma, face_type, nface, face, mod_grid_get_face_ngl, is_cg_coupled

    use mod_initial, only: nvar

    use mod_input, only: icase, space_method, cgdg_method

    use mod_p4est, only: scatter_element_2d, gather_element_2d, plist

    implicit none

    !global arrays
    real, intent(inout) :: rhs(3,npoin)
    real, intent(in)    :: q(npoin)
    logical, intent(in) :: lpen

    !local arrays
    real, dimension(ngl,ngl) :: ql
    real, dimension(ngl,ngl) :: flux_ql, flux_lq
    real, dimension(ngl,ngl) :: nxl, nyl, nzl, nxr, nyr, nzr
    real, dimension(ngl,ngl,FACE_CHILDREN) :: qr, qlp
    real, dimension(ngl,ngl,FACE_CHILDREN) :: flux_q, flux_qr
    real, dimension(ngl,ngl,FACE_CHILDREN)   :: jac_face_r

    !local variables
    real unl, nx, ny, nz
    real:: wq, pen = 0
    integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er, itype, jface
    integer iel, ier, ilocl, ilocr,  ifa, ip, ic, icc, ivar
    integer ngl_i, ngl_j, plane_ij
    logical compute_boundary_flux, compute_mixed_flux

    if(space_method == 'cgc' .or. (space_method == 'cgd' .and. cgdg_method == 'separate')) return
    compute_boundary_flux = (space_method == 'cgd' .and. cgdg_method == 'unified')
    compute_mixed_flux    = (space_method == 'cgd' .and. cgdg_method == 'mixed')

    !initialize local arrays
    flux_q = 0
    flux_lq = 0
    flux_ql = 0
    flux_qr = 0
    ql = 0
    qr = 0

    !Construct FVM-type Operators
    do jface=1,nface_nc
        iface = face_nc(jface)
        itype = face_type(iface)   !Type of face

        !Skip boundary faces
        if (itype == 2 .or. itype == 21) then
            cycle
        end if

        !Skip periodic
        ier=face(8,iface)
        if (ier == -6) then
            cycle
        end if

        !Skip Interior faces for CGd
        if (compute_boundary_flux .and. ier > 0) then
            cycle
        end if

        !-------------------------------------
        !Store Left and Right Side Variables
        !-------------------------------------
        ilocl=face(5,iface)
        iel=face(7,iface)
        call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

        !----------------------------Left Element
        do j=1,ngl_j
            do i=1,ngl_i
                !Get Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                ip=intma(il,jl,kl,iel)

                ql(i,j)=q(ip)

                !Normal Vector
                nxl(i,j)=normal_vector(1,i,j,iface)
                nyl(i,j)=normal_vector(2,i,j,iface)
                nzl(i,j)=normal_vector(3,i,j,iface)
            end do !j
        end do !i

        !----------------------------Right Element
        do ic=1,FACE_CHILDREN !four children
            ilocr=face(6,iface)
            icc = ic !reorder_nc_face(ic,ilocr)
            ier = face(7+icc, iface)

            if (ier > 0 ) then
                do j=1,ngl_j
                    do i=1,ngl_i

                        !Get Pointers
                        ir=imapr(1,i,j,iface)
                        jr=imapr(2,i,j,iface)
                        kr=imapr(3,i,j,iface)
                        ip=intma(ir,jr,kr,ier)

                        qr(i,j,ic)=q(ip)
                    end do !i
                end do !j
            end if
        end do


        !Right Normals
        !     nxr=-nxl
        !     nyr=-nyl
        !     nzr=-nzl

        !---------------------------------------------
        ! Project left element data to right elements
        !---------------------------------------------
        call scatter_element_2d(ql(:,:),qlp(:,:,:), ngl_i, ngl_j, plane_ij)

        !----------------------------------
        !  Compute face and element fluxes
        !----------------------------------
        flux_q = 0; flux_lq = 0; flux_ql = 0; flux_qr = 0

        do ic=1,FACE_CHILDREN

            icc = ic !reorder_nc_face(ic,ilocr)
            ier = face(7+icc,iface)

            if (ier>0) then
                !Get normals (redundant) and jacobian for the right faces
                call compute_normals_element(nxr,nyr,nzr,jac_face_r(:,:,ic),ier,ilocr)

                !Make normals outward-facing wrt left element
                nxr = -nxr
                nyr = -nyr
                nzr = -nzr

                flux_ql(:,:) = qlp(:,:,ic)
                flux_qr(:,:,ic) = qr(:,:,ic)
                flux_q(:,:,ic) = 0.5 * (flux_ql(:,:) + flux_qr(:,:,ic))

            end if
        end do

        !------------------------------------------
        ! Project flux data back to left element
        !------------------------------------------
        call gather_element_2d(flux_q(:,:,:),flux_lq(:,:), ngl_i, ngl_j, plane_ij)

        !----------------------------------
        !  Compute element flux for left element
        !----------------------------------
        flux_ql = ql

        !---------------------------------
        !  Do Gauss-Lobatto Integration
        !---------------------------------

        !--------------Left Side------------------!
        do j=1,ngl_j
            do i=1,ngl_i
                !Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                wq=jac_face(i,j,iface)
                ip=intma(il,jl,kl,iel)
                if (.not.(compute_mixed_flux .and. is_cg_coupled(il,jl,kl,iel) .ne. 0)) then
                    if(lpen) pen = penalty(i,j,iface)
                    unl = (flux_ql(i,j) - flux_lq(i,j)) + pen*(flux_lq(i,j) - flux_ql(i,j)) 
                    rhs(1,ip)=rhs(1,ip) + wq * unl * nxl(i,j)
                    rhs(2,ip)=rhs(2,ip) + wq * unl * nyl(i,j)
                    rhs(3,ip)=rhs(3,ip) + wq * unl * nzl(i,j)
                endif
            end do !i
        end do !j

        !--------------Right Side------------------!
        do ic=1,FACE_CHILDREN
            ier = face(7+ic,iface)
            if (ier > 0) then
                do j=1,ngl_j
                    do i=1,ngl_i
                        !Pointers
                        ir=imapr(1,i,j,iface)
                        jr=imapr(2,i,j,iface)
                        kr=imapr(3,i,j,iface)
                        wq=jac_face_r(i,j,ic)
                        ip=intma(ir,jr,kr,ier)
                        if(.not.(compute_mixed_flux .and. is_cg_coupled(ir,jr,kr,ier) .ne. 0)) then
                            if(lpen) pen = penalty(i,j,iface)
                            unl = (flux_qr(i,j,ic) - flux_q(i,j,ic)) + pen*(flux_qr(i,j,ic) - flux_ql(i,j)) 
                            rhs(1,ip)=rhs(1,ip) + wq * unl * nxr(i,j)
                            rhs(2,ip)=rhs(2,ip) + wq * unl * nyr(i,j)
                            rhs(3,ip)=rhs(3,ip) + wq * unl * nzr(i,j)
                        endif
                    end do !i
                end do !j
            end if !ier
        end do !ic

    end do !jface
end subroutine 

!-------------------------------------------------------------
!@>brief Calculate jump penalty flux for gradient on ip faces
!-------------------------------------------------------------
subroutine compute_gradient_penalty_flux_ip(rhs,q_send,q_recv,lpen)

    use mod_basis, only: ngl, FACE_CHILDREN

    use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send, penalty

    use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_ngl, face_type

    use mod_metrics, only: jac

    use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list

    use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv

    use mod_ref, only: nmessage

    use mod_input, only: space_method, cgdg_method

    implicit none

    !global arrays
    real, intent(inout) :: rhs(3,npoin)
    real,intent(in):: q_send(nmessage,ngl,ngl,nboun)
    real,intent(in):: q_recv(nmessage,ngl,ngl,nboun)
    logical, intent(in) :: lpen

    !local variables
    real, dimension(ngl,ngl) :: flux_q, flux_ql, flux_qr, flux_lq
    real, dimension(ngl,ngl) :: ql, qr, qlp
    real, dimension(ngl,ngl) :: jac_face_l
    real, dimension(ngl,ngl) :: nxl, nyl, nzl, nxr, nyr, nzr
    real, dimension(ngl,ngl,4) :: nxrp, nyrp, nzrp

    real:: wq, pen = 0, qt(3), unl
    integer iface, i, j, k, l, ip, il, jl, kl, ir, jr, kr
    integer iel, ier, ilocl, ilocr
    integer ifaceb, ngl_i, ngl_j, plane_ij

    integer isub, ic, jj, im, imm, inbh, ib, kk, ivar
    integer ftype, pface, subface, imulti
    logical compute_boundary_flux, compute_mixed_flux

    if(space_method == 'cgc' .or. (space_method == 'cgd' .and. cgdg_method == 'separate')) return
    compute_boundary_flux = (space_method == 'cgd' .and. cgdg_method == 'unified')
    compute_mixed_flux    = (space_method == 'cgd' .and. cgdg_method == 'mixed')

    !Constants
    pface = 0
    isub = 1

    jj=1
    kk=1
    imm = 0

    do inbh = 1, num_nbh
        do ib=1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)

            ftype = face_type(iface)

            do im = 1,imulti

                !---------------------
                ! determine face type
                !---------------------
                if (ftype==2) then
                    ilocl=face(5,iface)
                    iel=face(7,iface)
                else if (ftype==21) then
                    ilocl=face(6,iface)
                    i=1
                    ic=0
                    do while(ic<im)
                        if (face(7+i,iface)>0) then
                            ic=ic+1
                            iel=face(7+i,iface)
                            subface = i
                        end if
                        i=i+1
                    end do
                else if (ftype==12) then
                    ilocl=face(5,iface)
                    i=1
                    ic=0
                    do while(ic<im)
                        if (nbh_send_recv_half((jj-1)*FACE_CHILDREN+i)==1) then
                            ic=ic+1
                            subface=i
                        end if
                        i=i+1
                    end do
                    iel=face(7,iface)
                end if

                !-------------------------------------
                !Store Left/Right Side Variables
                !-------------------------------------
                call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

                ql(:,:)=q_send(1,:,:,kk)
                qr(:,:)=q_recv(1,:,:,kk)

                !project parent face
                if (ftype==12) then
                    qlp=ql
                    call scatter_element_2d_subface(qlp(:,:),ql(:,:),subface, ngl_i, ngl_j, plane_ij)
                else if (ftype==21) then
                    qlp=qr
                    call scatter_element_2d_subface(qlp(:,:),qr(:,:),subface, ngl_i, ngl_j, plane_ij)
                end if

                !------------------
                ! Store normals
                !------------------
                if (ftype==2) then
                    do j=1,ngl_j
                        do i=1,ngl_i
                            !Left Normals
                            nxl(i,j)=normal_vector(1,i,j,iface)
                            nyl(i,j)=normal_vector(2,i,j,iface)
                            nzl(i,j)=normal_vector(3,i,j,iface)
                            jac_face_l(i,j) = jac_face(i,j,iface)
                        end do !i
                    end do !j
                else if (ftype==21) then
                    call compute_normals_element(nxl,nyl,nzl,jac_face_l,iel,ilocl)
                else if (ftype==12) then
                    do i=1,ngl_i
                        do j=1,ngl_j
                            !Left Normals
                            nxl(i,j)=normal_vector(1,i,j,iface)
                            nyl(i,j)=normal_vector(2,i,j,iface)
                            nzl(i,j)=normal_vector(3,i,j,iface)
                            jac_face_l(i,j) = jac_face(i,j,iface)
                        end do !j
                    end do !i
                    !------------------------------------------
                    ! Project normals
                    !------------------------------------------
                    call scatter_element_2d(nxl,nxrp, ngl_i, ngl_j, plane_ij)
                    call scatter_element_2d(nyl,nyrp, ngl_i, ngl_j, plane_ij)
                    call scatter_element_2d(nzl,nzrp, ngl_i, ngl_j, plane_ij)

                    nxl = nxrp(:,:,subface)
                    nyl = nyrp(:,:,subface)
                    nzl = nzrp(:,:,subface)

                end if
                
                nxr=-nxl
                nyr=-nyl
                nzr=-nzl
                
                !----------------------------------
                !  Compute face and element fluxes
                !----------------------------------
                flux_ql = ql
                flux_qr = qr
                flux_q = 0.5 * (flux_ql + flux_qr)

                !---------------------------------
                !  Do Gauss-Lobatto Integration
                !---------------------------------
                if (ftype==12) then
                    !need to backward-project the flux
                    call gather_element_2d_subface(flux_q(:,:),flux_lq(:,:),subface, ngl_i, ngl_j, plane_ij)
                    flux_q=flux_lq

                    !elemental flux needs to be computed differently for face 12 - will do it in the regular nc subroutine
                    flux_ql = 0
                end if
        
                do j=1,ngl_j
                    do i=1,ngl_i
                        !Pointers
                        if (ftype==21) then
                            il=imapr(1,i,j,iface)
                            jl=imapr(2,i,j,iface)
                            kl=imapr(3,i,j,iface)
                        else
                            il=imapl(1,i,j,iface)
                            jl=imapl(2,i,j,iface)
                            kl=imapl(3,i,j,iface)
                        end if
                        wq=jac_face_l(i,j)
                        ip=intma(il,jl,kl,iel)

                        !Store RHS
                        if(lpen) pen = penalty(i,j,iface)
                        unl = (flux_q(i,j) - flux_ql(i,j)) + pen*(flux_qr(i,j) - flux_ql(i,j)) 
                        rhs(1,ip)=rhs(1,ip) + wq * unl * nxl(i,j)
                        rhs(2,ip)=rhs(2,ip) + wq * unl * nyl(i,j)
                        rhs(3,ip)=rhs(3,ip) + wq * unl * nzl(i,j)

                    end do !i
                end do !j

                kk=kk+1
            end do
            jj=jj+1
        end do

    end do !iface

end subroutine

!-------------------------------------------------------------
!@>brief Calculate jump penalty flux for divergence
!-------------------------------------------------------------
subroutine compute_divergence_penalty_flux(rhs,q,lpen)

    use mod_basis, only: ngl

    use mod_face, only: normal_vector, jac_face, imapl, imapr, penalty

    use mod_grid, only: npoin, intma, face_type, nface, face, mod_grid_get_face_ngl, is_cg_coupled

    use mod_input, only: is_non_conforming_flg, icase, space_method, cgdg_method

    implicit none

    !global arrays
    real, intent(inout) :: rhs(npoin)
    real, intent(in) :: q(3,npoin)
    logical, intent(in) :: lpen

    !local arrays
    real, dimension(3,ngl,ngl) :: ql, qr
    real, dimension(ngl,ngl) :: nxl, nyl, nzl, nxr, nyr, nzr

    !local variables
    real unl, wq, qt(3)
    real :: pen = 0
    integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er, itype
    integer iel, ier, ilocl, ilocr,  ifa, ip
    integer ngl_i, ngl_j, plane_ij, bcflag
    logical compute_boundary_flux, compute_mixed_flux

    if(space_method == 'cgc' .or. (space_method == 'cgd' .and. cgdg_method == 'separate')) return
    compute_boundary_flux = (space_method == 'cgd' .and. cgdg_method == 'unified')
    compute_mixed_flux    = (space_method == 'cgd' .and. cgdg_method == 'mixed')

    !Construct FVM-type Operators
    do iface=1,nface
        itype = face_type(iface)   !Type of face

        !Skip boundary faces
        if (itype == 2 .or. itype == 21 .or. itype == 12) then
            cycle
        end if

        !Skip periodic faces
        ier=face(8,iface)
        if (ier == -6) then
            cycle
        end if

        !Skip Interior Faces for CGD
        if (compute_boundary_flux .and. ier > 0) then
            cycle
        end if

        !Skip non-conforming faces
        if (is_non_conforming_flg>0) then
            if (face(9,iface)>0) cycle
        end if

        !-------------------------------------
        !Store Left and Right Side Variables
        !-------------------------------------
        ilocl=face(5,iface)
        iel=face(7,iface)
        call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

        !----------------------------Left Element
        do j=1,ngl_j
            do i=1,ngl_i
                !Get Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                ip=intma(il,jl,kl,iel)

                ql(1:3,i,j)=q(1:3,ip)

                !Normal Vector
                nxl(i,j)=normal_vector(1,i,j,iface)
                nyl(i,j)=normal_vector(2,i,j,iface)
                nzl(i,j)=normal_vector(3,i,j,iface)
            end do !j
        end do !i

        !----------------------------Right Element
        if (ier > 0 ) then
            ilocr=face(6,iface)
            do j=1,ngl_j
                do i=1,ngl_i
                    !Get Pointers
                    ir=imapr(1,i,j,iface)
                    jr=imapr(2,i,j,iface)
                    kr=imapr(3,i,j,iface)
                    ip=intma(ir,jr,kr,ier)

                    qr(1:3,i,j)=q(1:3,ip)
                end do !i
            end do !j

        !------------------------------
        !  Apply boundary conditions
        !------------------------------
        else if (ier <= 0) then
            bcflag = mod(ier, 10)
            
            !default to qr=ql
            qr=ql
            
            if (bcflag == -4) then    !no-flux
                do j=1,ngl_j
                    do i=1,ngl_i
                        unl=nxl(i,j)*ql(1,i,j) + nyl(i,j)*ql(2,i,j) + nzl(i,j)*ql(3,i,j)
                        qr(1,i,j)=ql(1,i,j) - 2*unl*nxl(i,j)
                        qr(2,i,j)=ql(2,i,j) - 2*unl*nyl(i,j)
                        qr(3,i,j)=ql(3,i,j) - 2*unl*nzl(i,j)
                    end do !i
                end do !j
            else if (bcflag == -6) then     !outflow
                do j=1,ngl_j
                    do i=1,ngl_i
                        unl=nxl(i,j)*ql(1,i,j) + nyl(i,j)*ql(2,i,j) + nzl(i,j)*ql(3,i,j)
                        qr(1,i,j)=2*ql(1,i,j) - 2*unl*nxl(i,j)
                        qr(2,i,j)=2*ql(2,i,j) - 2*unl*nyl(i,j)
                        qr(3,i,j)=2*ql(3,i,j) - 2*unl*nzl(i,j)
                    end do !i
                end do !j
            end if
            
        end if

        !Right Normals
        nxr=-nxl
        nyr=-nyl
        nzr=-nzl

        !---------------------------------
        !  Do Gauss-Lobatto Integration
        !---------------------------------

        !--------------Left Side------------------!
        do j=1,ngl_j
            do i=1,ngl_i
                !Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                wq=jac_face(i,j,iface)
                ip=intma(il,jl,kl,iel)
                if (.not.(compute_mixed_flux .and. is_cg_coupled(il,jl,kl,iel) .ne. 0)) then
                    if(lpen) pen = penalty(i,j,iface)
                    qt(1:3) = (pen + 0.5) * (qr(1:3,i,j) - ql(1:3,i,j))
                    unl = nxl(i,j)*qt(1) + nyl(i,j)*qt(2) + nzl(i,j)*qt(3)
                    rhs(ip)=rhs(ip) + wq*unl
                end if
            end do !i
        end do !j

        !--------------Right Side------------------!
        if (ier > 0) then
            do j=1,ngl_j
                do i=1,ngl_i
                    !Pointers
                    ir=imapr(1,i,j,iface)
                    jr=imapr(2,i,j,iface)
                    kr=imapr(3,i,j,iface)
                    wq=jac_face(i,j,iface)
                    ip=intma(ir,jr,kr,ier)
                    if (.not.(compute_mixed_flux .and. is_cg_coupled(ir,jr,kr,ier) .ne. 0)) then
                        if(lpen) pen = penalty(i,j,iface)
                        qt(1:3) = (pen + 0.5) * (ql(1:3,i,j) - qr(1:3,i,j))
                        unl = nxr(i,j)*qt(1) + nyr(i,j)*qt(2) + nzr(i,j)*qt(3)
                        rhs(ip)=rhs(ip) + wq*unl
                    end if
                end do !i
            end do !j
        end if !ier

    end do !iface

end subroutine

!-------------------------------------------------------------
!@>brief Calculate jump penalty flux for divergence on nc faces
!-------------------------------------------------------------
subroutine compute_divergence_penalty_flux_nc(rhs,q,lpen)

    use mod_basis, only: ngl, is_2d, FACE_CHILDREN

    use mod_face, only: normal_vector, jac_face, imapl, imapr, face_nc, nface_nc, penalty

    use mod_grid, only: npoin, intma, face_type, nface, face, mod_grid_get_face_ngl, is_cg_coupled

    use mod_initial, only: nvar

    use mod_input, only: icase, space_method, cgdg_method

    use mod_p4est, only: scatter_element_2d, gather_element_2d, plist

    implicit none

    !global arrays
    real, intent(inout) :: rhs(npoin)
    real, intent(in)    :: q(3,npoin)
    logical, intent(in) :: lpen

    !local arrays
    real, dimension(3,ngl,ngl) :: ql
    real, dimension(3,ngl,ngl) :: flux_ql, flux_lq
    real, dimension(ngl,ngl) :: nxl, nyl, nzl, nxr, nyr, nzr
    real, dimension(3,ngl,ngl,FACE_CHILDREN) :: qr, qlp
    real, dimension(3,ngl,ngl,FACE_CHILDREN) :: flux_q, flux_qr
    real, dimension(ngl,ngl,FACE_CHILDREN)   :: jac_face_r

    !local variables
    real unl, nx, ny, nz
    real:: wq, pen = 0, qt(3)
    integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er, itype, jface
    integer iel, ier, ilocl, ilocr,  ifa, ip, ic, icc, ivar
    integer ngl_i, ngl_j, plane_ij
    logical compute_boundary_flux, compute_mixed_flux

    if(space_method == 'cgc' .or. (space_method == 'cgd' .and. cgdg_method == 'separate')) return
    compute_boundary_flux = (space_method == 'cgd' .and. cgdg_method == 'unified')
    compute_mixed_flux    = (space_method == 'cgd' .and. cgdg_method == 'mixed')

    !initialize local arrays
    flux_q = 0
    flux_lq = 0
    flux_ql = 0
    flux_qr = 0
    ql = 0
    qr = 0

    !Construct FVM-type Operators
    do jface=1,nface_nc
        iface = face_nc(jface)
        itype = face_type(iface)   !Type of face

        !Skip boundary faces
        if (itype == 2 .or. itype == 21) then
            cycle
        end if

        !Skip periodic
        ier=face(8,iface)
        if (ier == -6) then
            cycle
        end if

        !Skip Interior faces for CGd
        if (compute_boundary_flux .and. ier > 0) then
            cycle
        end if

        !-------------------------------------
        !Store Left and Right Side Variables
        !-------------------------------------
        ilocl=face(5,iface)
        iel=face(7,iface)
        call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

        !----------------------------Left Element
        do j=1,ngl_j
            do i=1,ngl_i
                !Get Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                ip=intma(il,jl,kl,iel)

                ql(:,i,j)=q(:,ip)

                !Normal Vector
                nxl(i,j)=normal_vector(1,i,j,iface)
                nyl(i,j)=normal_vector(2,i,j,iface)
                nzl(i,j)=normal_vector(3,i,j,iface)
            end do !j
        end do !i

        !----------------------------Right Element
        do ic=1,FACE_CHILDREN !four children
            ilocr=face(6,iface)
            icc = ic !reorder_nc_face(ic,ilocr)
            ier = face(7+icc, iface)

            if (ier > 0 ) then
                do j=1,ngl_j
                    do i=1,ngl_i

                        !Get Pointers
                        ir=imapr(1,i,j,iface)
                        jr=imapr(2,i,j,iface)
                        kr=imapr(3,i,j,iface)
                        ip=intma(ir,jr,kr,ier)

                        qr(:,i,j,ic)=q(:,ip)
                    end do !i
                end do !j
            end if
        end do


        !Right Normals
        !     nxr=-nxl
        !     nyr=-nyl
        !     nzr=-nzl

        !---------------------------------------------
        ! Project left element data to right elements
        !---------------------------------------------
        do ivar = 1,3
            call scatter_element_2d(ql(ivar,:,:),qlp(ivar,:,:,:), ngl_i, ngl_j, plane_ij)
        enddo

        !----------------------------------
        !  Compute face and element fluxes
        !----------------------------------
        flux_q = 0; flux_lq = 0; flux_ql = 0; flux_qr = 0

        do ic=1,FACE_CHILDREN

            icc = ic !reorder_nc_face(ic,ilocr)
            ier = face(7+icc,iface)

            if (ier>0) then
                !Get normals (redundant) and jacobian for the right faces
                call compute_normals_element(nxr,nyr,nzr,jac_face_r(:,:,ic),ier,ilocr)

                !Make normals outward-facing wrt left element
                nxr = -nxr
                nyr = -nyr
                nzr = -nzr

                flux_ql(:,:,:) = qlp(:,:,:,ic)
                flux_qr(:,:,:,ic) = qr(:,:,:,ic)
                flux_q(:,:,:,ic) = 0.5 * (flux_ql(:,:,:) + flux_qr(:,:,:,ic))

            end if
        end do

        !------------------------------------------
        ! Project flux data back to left element
        !------------------------------------------
        call gather_element_2d(flux_q(:,:,:,:),flux_lq(:,:,:), ngl_i, ngl_j, plane_ij)

        !----------------------------------
        !  Compute element flux for left element
        !----------------------------------
        flux_ql = ql

        !---------------------------------
        !  Do Gauss-Lobatto Integration
        !---------------------------------

        !--------------Left Side------------------!
        do j=1,ngl_j
            do i=1,ngl_i
                !Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                wq=jac_face(i,j,iface)
                ip=intma(il,jl,kl,iel)
                if (.not.(compute_mixed_flux .and. is_cg_coupled(il,jl,kl,iel) .ne. 0)) then
                    if(lpen) pen = penalty(i,j,iface)
                    qt(1:3) = (flux_ql(:,i,j) - flux_lq(:,i,j)) + pen*(flux_lq(:,i,j) - flux_ql(:,i,j))
                    unl = nxl(i,j)*qt(1) + nyl(i,j)*qt(2) + nzl(i,j)*qt(3)
                    rhs(ip)=rhs(ip) + wq*unl
                endif
            end do !i
        end do !j

        !--------------Right Side------------------!
        do ic=1,FACE_CHILDREN
            ier = face(7+ic,iface)
            if (ier > 0) then
                do j=1,ngl_j
                    do i=1,ngl_i
                        !Pointers
                        ir=imapr(1,i,j,iface)
                        jr=imapr(2,i,j,iface)
                        kr=imapr(3,i,j,iface)
                        wq=jac_face_r(i,j,ic)
                        ip=intma(ir,jr,kr,ier)
                        if(.not.(compute_mixed_flux .and. is_cg_coupled(ir,jr,kr,ier) .ne. 0)) then
                            if(lpen) pen = penalty(i,j,iface)
                            qt(1:3) = (flux_qr(:,i,j,ic) - flux_q(:,i,j,ic)) + pen*(flux_qr(:,i,j,ic) - flux_ql(:,i,j))
                            unl = nxr(i,j)*qt(1) + nyr(i,j)*qt(2) + nzr(i,j)*qt(3)
                            rhs(ip)=rhs(ip) + wq*unl
                        endif
                    end do !i
                end do !j
            end if !ier
        end do !ic

    end do !jface
end subroutine 

!-------------------------------------------------------------
!@>brief Calculate jump penalty flux for divergence on ip faces
!-------------------------------------------------------------
subroutine compute_divergence_penalty_flux_ip(rhs,q_send,q_recv,lpen)

    use mod_basis, only: ngl, FACE_CHILDREN

    use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send, penalty

    use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_ngl, face_type

    use mod_metrics, only: jac

    use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list

    use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv

    use mod_ref, only: nmessage

    use mod_input, only: space_method, cgdg_method

    implicit none

    !global arrays
    real, intent(inout) :: rhs(npoin)
    real,intent(in):: q_send(nmessage,ngl,ngl,nboun)
    real,intent(in):: q_recv(nmessage,ngl,ngl,nboun)
    logical, intent(in) :: lpen

    !local variables
    real, dimension(3,ngl,ngl) :: flux_q, flux_ql, flux_qr, flux_lq
    real, dimension(3,ngl,ngl) :: ql, qr, qlp
    real, dimension(ngl,ngl) :: jac_face_l
    real, dimension(ngl,ngl) :: nxl, nyl, nzl, nxr, nyr, nzr
    real, dimension(ngl,ngl,4) :: nxrp, nyrp, nzrp

    real:: wq, pen = 0, qt(3), unl
    integer iface, i, j, k, l, ip, il, jl, kl, ir, jr, kr
    integer iel, ier, ilocl, ilocr
    integer ifaceb, ngl_i, ngl_j, plane_ij

    integer isub, ic, jj, im, imm, inbh, ib, kk, ivar
    integer ftype, pface, subface, imulti
    logical compute_boundary_flux, compute_mixed_flux

    if(space_method == 'cgc' .or. (space_method == 'cgd' .and. cgdg_method == 'separate')) return
    compute_boundary_flux = (space_method == 'cgd' .and. cgdg_method == 'unified')
    compute_mixed_flux    = (space_method == 'cgd' .and. cgdg_method == 'mixed')

    !Constants
    pface = 0
    isub = 1

    jj=1
    kk=1
    imm = 0

    do inbh = 1, num_nbh
        do ib=1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)

            ftype = face_type(iface)

            do im = 1,imulti

                !---------------------
                ! determine face type
                !---------------------
                if (ftype==2) then
                    ilocl=face(5,iface)
                    iel=face(7,iface)
                else if (ftype==21) then
                    ilocl=face(6,iface)
                    i=1
                    ic=0
                    do while(ic<im)
                        if (face(7+i,iface)>0) then
                            ic=ic+1
                            iel=face(7+i,iface)
                            subface = i
                        end if
                        i=i+1
                    end do
                else if (ftype==12) then
                    ilocl=face(5,iface)
                    i=1
                    ic=0
                    do while(ic<im)
                        if (nbh_send_recv_half((jj-1)*FACE_CHILDREN+i)==1) then
                            ic=ic+1
                            subface=i
                        end if
                        i=i+1
                    end do
                    iel=face(7,iface)
                end if

                !-------------------------------------
                !Store Left/Right Side Variables
                !-------------------------------------
                call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

                ql(:,:,:)=q_send(1:3,:,:,kk)
                qr(:,:,:)=q_recv(1:3,:,:,kk)

                !project parent face
                if (ftype==12) then
                    qlp=ql
                    do ivar = 1,3
                        call scatter_element_2d_subface(qlp(ivar,:,:),ql(ivar,:,:),subface, ngl_i, ngl_j, plane_ij)
                    enddo
                else if (ftype==21) then
                    qlp=qr
                    do ivar = 1,3
                        call scatter_element_2d_subface(qlp(ivar,:,:),qr(ivar,:,:),subface, ngl_i, ngl_j, plane_ij)
                    enddo
                end if

                !------------------
                ! Store normals
                !------------------
                if (ftype==2) then
                    do j=1,ngl_j
                        do i=1,ngl_i
                            !Left Normals
                            nxl(i,j)=normal_vector(1,i,j,iface)
                            nyl(i,j)=normal_vector(2,i,j,iface)
                            nzl(i,j)=normal_vector(3,i,j,iface)
                            jac_face_l(i,j) = jac_face(i,j,iface)
                        end do !i
                    end do !j
                else if (ftype==21) then
                    call compute_normals_element(nxl,nyl,nzl,jac_face_l,iel,ilocl)
                else if (ftype==12) then
                    do i=1,ngl_i
                        do j=1,ngl_j
                            !Left Normals
                            nxl(i,j)=normal_vector(1,i,j,iface)
                            nyl(i,j)=normal_vector(2,i,j,iface)
                            nzl(i,j)=normal_vector(3,i,j,iface)
                            jac_face_l(i,j) = jac_face(i,j,iface)
                        end do !j
                    end do !i
                    !------------------------------------------
                    ! Project normals
                    !------------------------------------------
                    call scatter_element_2d(nxl,nxrp, ngl_i, ngl_j, plane_ij)
                    call scatter_element_2d(nyl,nyrp, ngl_i, ngl_j, plane_ij)
                    call scatter_element_2d(nzl,nzrp, ngl_i, ngl_j, plane_ij)

                    nxl = nxrp(:,:,subface)
                    nyl = nyrp(:,:,subface)
                    nzl = nzrp(:,:,subface)

                end if
                
                nxr=-nxl
                nyr=-nyl
                nzr=-nzl
                
                !----------------------------------
                !  Compute face and element fluxes
                !----------------------------------
                flux_ql = ql
                flux_qr = qr
                flux_q = 0.5 * (flux_ql + flux_qr)

                !---------------------------------
                !  Do Gauss-Lobatto Integration
                !---------------------------------
                if (ftype==12) then
                    !need to backward-project the flux
                    do ivar = 1,3
                        call gather_element_2d_subface(flux_q(ivar,:,:),flux_lq(ivar,:,:),subface, ngl_i, ngl_j, plane_ij)
                    enddo
                    flux_q=flux_lq

                    !elemental flux needs to be computed differently for face 12 - will do it in the regular nc subroutine
                    flux_ql = 0
                end if
        
                do j=1,ngl_j
                    do i=1,ngl_i
                        !Pointers
                        if (ftype==21) then
                            il=imapr(1,i,j,iface)
                            jl=imapr(2,i,j,iface)
                            kl=imapr(3,i,j,iface)
                        else
                            il=imapl(1,i,j,iface)
                            jl=imapl(2,i,j,iface)
                            kl=imapl(3,i,j,iface)
                        end if
                        wq=jac_face_l(i,j)
                        ip=intma(il,jl,kl,iel)

                        !Store RHS
                        if(lpen) pen = penalty(i,j,iface)
                        qt(1:3) = (flux_q(:,i,j) - flux_ql(:,i,j)) + pen*(flux_qr(:,i,j) - flux_ql(:,i,j)) 
                        unl = nxl(i,j)*qt(1) + nyl(i,j)*qt(2) + nzl(i,j)*qt(3)
                        rhs(ip)=rhs(ip) + wq*unl

                    end do !i
                end do !j

                kk=kk+1
            end do
            jj=jj+1
        end do

    end do !iface

end subroutine

!-------------------------------------------------------------
!@>brief Calculate jump penalty flux for laplacian
!-------------------------------------------------------------
subroutine compute_laplacian_penalty_flux(rhs,q,lpen)

    use mod_basis, only: ngl

    use mod_face, only: normal_vector, jac_face, imapl, imapr, penalty

    use mod_grid, only: npoin, intma, face_type, nface, face, mod_grid_get_face_ngl, is_cg_coupled

    use mod_input, only: icase, space_method, cgdg_method, visc, is_non_conforming_flg

    implicit none

    !global arrays
    real, intent(inout) :: rhs(npoin)
    real, intent(in)    :: q(npoin)
    logical, intent(in) :: lpen
    
    !local arrays
    real, dimension(ngl,ngl) :: ql, qr

    !local variables
    real unl, wq
    real :: pen = 0
    integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er, itype
    integer iel, ier, ilocl, ilocr,  ifa, ip
    integer ngl_i, ngl_j, plane_ij, bcflag
    logical compute_boundary_flux, compute_mixed_flux

    if(space_method == 'cgc' .or. (space_method == 'cgd' .and. cgdg_method == 'separate')) return
    compute_boundary_flux = (space_method == 'cgd' .and. cgdg_method == 'unified')
    compute_mixed_flux    = (space_method == 'cgd' .and. cgdg_method == 'mixed')
  
    !Construct FVM-type Operators
    do iface=1,nface
        itype = face_type(iface)   !Type of face

        !Skip boundary faces
        if (itype == 2 .or. itype == 21 .or. itype == 12) then
            cycle
        end if

        !Skip periodic
        ier=face(8,iface)
        er=ier
        if (ier == -6) then
            cycle
        end if

        !Compute boundary flux
        if(compute_boundary_flux .and. ier > 0) then
            cycle
        end if

        !Skip non-conforming faces
        if(is_non_conforming_flg>0) then
            if(face(9,iface)>0) cycle
        end if

        !-------------------------------------
        !Store Left and Right Side Variables
        !-------------------------------------
        ilocl=face(5,iface)
        iel=face(7,iface)
        el=iel
        call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

        !----------------------------Left Element
        do j=1,ngl_j
            do i=1,ngl_i
                !Get Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                ip=intma(il,jl,kl,iel)

                ql(i,j)=q(ip)
            end do !i
        end do !j

        !----------------------------Right Element
        if (ier > 0 ) then
            ilocr=face(6,iface)
            do j=1,ngl_j
                do i=1,ngl_i

                    !Get Pointers
                    ir=imapr(1,i,j,iface)
                    jr=imapr(2,i,j,iface)
                    kr=imapr(3,i,j,iface)
                    ip=intma(ir,jr,kr,ier)

                    qr(i,j)=q(ip)
                end do !i
            end do !j

        !----------------------------Boundary element
        else if (ier <= 0) then
            bcflag = mod(ier,10)
            
            !Use Left Element IMAP and Metric Terms
            er=el
            ilocr=ilocl

            !default values
            qr=ql
            
        end if

        !---------------------------------
        !  Do Gauss-Lobatto Integration
        !---------------------------------
        !--------------Left Side------------------!
        do j=1,ngl_j
            do i=1,ngl_i
                !Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                wq=jac_face(i,j,iface)
                ip=intma(il,jl,kl,iel)
                if(.not.(compute_mixed_flux .and. is_cg_coupled(il,jl,kl,iel) .ne. 0 )) then
                    if(lpen) pen = penalty(i,j,iface)
                    rhs(ip)=rhs(ip) + wq * (pen + 0.5) * (qr(i,j) - ql(i,j))
                endif

            end do !i
        end do !j

        !--------------Right Side------------------!
        if (ier > 0) then
            do j=1,ngl_j
                do i=1,ngl_i
                    !Pointers
                    ir=imapr(1,i,j,iface)
                    jr=imapr(2,i,j,iface)
                    kr=imapr(3,i,j,iface)
                    wq=jac_face(i,j,iface)
                    ip=intma(ir,jr,kr,ier)
                    if(.not.(compute_mixed_flux .and. is_cg_coupled(ir,jr,kr,ier) .ne. 0 )) then
                        if(lpen) pen = penalty(i,j,iface)
                        rhs(ip)=rhs(ip) + wq * (pen + 0.5) * (ql(i,j) - qr(i,j))
                    endif
                end do !i
            end do !j
        end if !ier

    end do !iface
    
end subroutine

!-------------------------------------------------------------
!@>brief Calculate jump penalty flux for laplacian on nc faces
!-------------------------------------------------------------
subroutine compute_laplacian_penalty_flux_nc(rhs,q,lpen)

    use mod_basis, only: ngl, is_2d, FACE_CHILDREN

    use mod_face, only: normal_vector, jac_face, imapl, imapr, face_nc, nface_nc, penalty

    use mod_grid, only: npoin, intma, face_type, nface, face, mod_grid_get_face_ngl, is_cg_coupled

    use mod_initial, only: nvar

    use mod_input, only: icase, space_method, cgdg_method

    use mod_p4est, only: scatter_element_2d, gather_element_2d, plist

    implicit none

    !global arrays
    real, intent(inout) :: rhs(npoin)
    real, intent(in)    :: q(npoin)
    logical, intent(in) :: lpen

    !local arrays
    real, dimension(ngl,ngl) :: ql
    real, dimension(ngl,ngl) :: flux_ql, flux_lq
    real, dimension(ngl,ngl) :: nxl, nyl, nzl, nxr, nyr, nzr
    real, dimension(ngl,ngl,FACE_CHILDREN) :: qr, qlp
    real, dimension(ngl,ngl,FACE_CHILDREN) :: flux_q, flux_qr
    real, dimension(ngl,ngl,FACE_CHILDREN)   :: jac_face_r

    !local variables
    real unl, nx, ny, nz
    real:: wq, pen = 0
    integer iface, i, j, k, il, jl, kl, ir, jr, kr, el, er, itype, jface
    integer iel, ier, ilocl, ilocr,  ifa, ip, ic, icc, ivar
    integer ngl_i, ngl_j, plane_ij
    logical compute_boundary_flux, compute_mixed_flux

    if(space_method == 'cgc' .or. (space_method == 'cgd' .and. cgdg_method == 'separate')) return
    compute_boundary_flux = (space_method == 'cgd' .and. cgdg_method == 'unified')
    compute_mixed_flux    = (space_method == 'cgd' .and. cgdg_method == 'mixed')

    !initialize local arrays
    flux_q = 0
    flux_lq = 0
    flux_ql = 0
    flux_qr = 0
    ql = 0
    qr = 0

    !Construct FVM-type Operators
    do jface=1,nface_nc
        iface = face_nc(jface)
        itype = face_type(iface)   !Type of face

        !Skip boundary faces
        if (itype == 2 .or. itype == 21) then
            cycle
        end if

        !Skip periodic
        ier=face(8,iface)
        if (ier == -6) then
            cycle
        end if

        !Skip Interior faces for CGd
        if (compute_boundary_flux .and. ier > 0) then
            cycle
        end if

        !-------------------------------------
        !Store Left and Right Side Variables
        !-------------------------------------
        ilocl=face(5,iface)
        iel=face(7,iface)
        call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

        !----------------------------Left Element
        do j=1,ngl_j
            do i=1,ngl_i
                !Get Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                ip=intma(il,jl,kl,iel)

                ql(i,j)=q(ip)

                !Normal Vector
                nxl(i,j)=normal_vector(1,i,j,iface)
                nyl(i,j)=normal_vector(2,i,j,iface)
                nzl(i,j)=normal_vector(3,i,j,iface)
            end do !j
        end do !i

        !----------------------------Right Element
        do ic=1,FACE_CHILDREN !four children
            ilocr=face(6,iface)
            icc = ic !reorder_nc_face(ic,ilocr)
            ier = face(7+icc, iface)

            if (ier > 0 ) then
                do j=1,ngl_j
                    do i=1,ngl_i

                        !Get Pointers
                        ir=imapr(1,i,j,iface)
                        jr=imapr(2,i,j,iface)
                        kr=imapr(3,i,j,iface)
                        ip=intma(ir,jr,kr,ier)

                        qr(i,j,ic)=q(ip)
                    end do !i
                end do !j
            end if
        end do


        !Right Normals
        !     nxr=-nxl
        !     nyr=-nyl
        !     nzr=-nzl

        !---------------------------------------------
        ! Project left element data to right elements
        !---------------------------------------------
        call scatter_element_2d(ql(:,:),qlp(:,:,:), ngl_i, ngl_j, plane_ij)

        !----------------------------------
        !  Compute face and element fluxes
        !----------------------------------
        flux_q = 0; flux_lq = 0; flux_ql = 0; flux_qr = 0

        do ic=1,FACE_CHILDREN

            icc = ic !reorder_nc_face(ic,ilocr)
            ier = face(7+icc,iface)

            if (ier>0) then
                !Get normals (redundant) and jacobian for the right faces
                call compute_normals_element(nxr,nyr,nzr,jac_face_r(:,:,ic),ier,ilocr)

                !Make normals outward-facing wrt left element
                nxr = -nxr
                nyr = -nyr
                nzr = -nzr

                flux_ql(:,:) = qlp(:,:,ic)
                flux_qr(:,:,ic) = qr(:,:,ic)
                flux_q(:,:,ic) = 0.5 * (flux_ql(:,:) + flux_qr(:,:,ic))

            end if
        end do

        !------------------------------------------
        ! Project flux data back to left element
        !------------------------------------------
        call gather_element_2d(flux_q(:,:,:),flux_lq(:,:), ngl_i, ngl_j, plane_ij)

        !----------------------------------
        !  Compute element flux for left element
        !----------------------------------
        flux_ql = ql

        !---------------------------------
        !  Do Gauss-Lobatto Integration
        !---------------------------------

        !--------------Left Side------------------!
        do j=1,ngl_j
            do i=1,ngl_i
                !Pointers
                il=imapl(1,i,j,iface)
                jl=imapl(2,i,j,iface)
                kl=imapl(3,i,j,iface)
                wq=jac_face(i,j,iface)
                ip=intma(il,jl,kl,iel)
                if (.not.(compute_mixed_flux .and. is_cg_coupled(il,jl,kl,iel) .ne. 0)) then
                    if(lpen) pen = penalty(i,j,iface)
                    rhs(ip)=rhs(ip) + wq*( (flux_ql(i,j) - flux_lq(i,j)) + pen*(flux_lq(i,j) - flux_ql(i,j)) )
                endif
            end do !i
        end do !j

        !--------------Right Side------------------!
        do ic=1,FACE_CHILDREN
            ier = face(7+ic,iface)
            if (ier > 0) then
                do j=1,ngl_j
                    do i=1,ngl_i
                        !Pointers
                        ir=imapr(1,i,j,iface)
                        jr=imapr(2,i,j,iface)
                        kr=imapr(3,i,j,iface)
                        wq=jac_face_r(i,j,ic)
                        ip=intma(ir,jr,kr,ier)
                        if(.not.(compute_mixed_flux .and. is_cg_coupled(ir,jr,kr,ier) .ne. 0)) then
                            if(lpen) pen = penalty(i,j,iface)
                            rhs(ip)=rhs(ip) + wq*( (flux_qr(i,j,ic) - flux_q(i,j,ic)) + pen*(flux_qr(i,j,ic) - flux_ql(i,j)) )
                        endif
                    end do !i
                end do !j
            end if !ier
        end do !ic

    end do !jface
end subroutine 

!-------------------------------------------------------------
!@>brief Calculate jump penalty flux for laplacian on ip faces
!-------------------------------------------------------------
subroutine compute_laplacian_penalty_flux_ip(rhs,q_send,q_recv,lpen)

    use mod_basis, only: ngl, FACE_CHILDREN

    use mod_face, only: normal_vector, jac_face, imapl, imapr, face_send, penalty

    use mod_grid, only: nelem, npoin, intma, face, nboun, mod_grid_get_face_ngl, face_type

    use mod_metrics, only: jac

    use mod_p4est, only: scatter_element_2d, scatter_element_2d_subface, gather_element_2d_subface, plist, lev_list

    use mod_parallel, only: nbh_send_recv, nbh_send_recv_multi, nbh_send_recv_half, num_nbh, num_send_recv

    use mod_ref, only: nmessage

    use mod_input, only: space_method, cgdg_method

    implicit none

    !global arrays
    real, intent(inout) :: rhs(npoin)
    real,intent(in):: q_send(nmessage,ngl,ngl,nboun)
    real,intent(in):: q_recv(nmessage,ngl,ngl,nboun)
    logical, intent(in) :: lpen

    !local variables
    real, dimension(ngl,ngl) :: flux_q, flux_ql, flux_qr, flux_lq
    real, dimension(ngl,ngl) :: ql, qr, qlp
    real, dimension(ngl,ngl) :: jac_face_l

    real:: wq, pen = 0
    integer iface, i, j, k, l, ip, il, jl, kl, ir, jr, kr
    integer iel, ier, ilocl, ilocr
    integer ifaceb, ngl_i, ngl_j, plane_ij

    integer isub, ic, jj, im, imm, inbh, ib, kk, ivar
    integer ftype, pface, subface, imulti
    logical compute_boundary_flux, compute_mixed_flux

    if(space_method == 'cgc' .or. (space_method == 'cgd' .and. cgdg_method == 'separate')) return
    compute_boundary_flux = (space_method == 'cgd' .and. cgdg_method == 'unified')
    compute_mixed_flux    = (space_method == 'cgd' .and. cgdg_method == 'mixed')

    !Constants
    pface = 0
    isub = 1

    jj=1
    kk=1
    imm = 0

    do inbh = 1, num_nbh
        do ib=1,num_send_recv(inbh)
            iface = nbh_send_recv(jj)
            imulti = nbh_send_recv_multi(jj)

            ftype = face_type(iface)

            do im = 1,imulti

                !---------------------
                ! determine face type
                !---------------------
                if (ftype==2) then
                    ilocl=face(5,iface)
                    iel=face(7,iface)
                else if (ftype==21) then
                    ilocl=face(6,iface)
                    i=1
                    ic=0
                    do while(ic<im)
                        if (face(7+i,iface)>0) then
                            ic=ic+1
                            iel=face(7+i,iface)
                            subface = i
                        end if
                        i=i+1
                    end do
                else if (ftype==12) then
                    ilocl=face(5,iface)
                    i=1
                    ic=0
                    do while(ic<im)
                        if (nbh_send_recv_half((jj-1)*FACE_CHILDREN+i)==1) then
                            ic=ic+1
                            subface=i
                        end if
                        i=i+1
                    end do
                    iel=face(7,iface)
                end if

                !-------------------------------------
                !Store Left/Right Side Variables
                !-------------------------------------
                call mod_grid_get_face_ngl(ilocl, ngl_i, ngl_j, plane_ij)

                ql(:,:)=q_send(1,:,:,kk)
                qr(:,:)=q_recv(1,:,:,kk)

                !project parent face
                if (ftype==12) then
                    qlp=ql
                    call scatter_element_2d_subface(qlp(:,:),ql(:,:),subface, ngl_i, ngl_j, plane_ij)
                else if (ftype==21) then
                    qlp=qr
                    call scatter_element_2d_subface(qlp(:,:),qr(:,:),subface, ngl_i, ngl_j, plane_ij)
                end if

                !------------------
                ! Store normals
                !------------------
                if (ftype==2) then
                    do j=1,ngl_j
                        do i=1,ngl_i
                            jac_face_l(i,j) = jac_face(i,j,iface)
                        end do !i
                    end do !j
                else if (ftype==21) then
                else if (ftype==12) then
                    do i=1,ngl_i
                        do j=1,ngl_j
                            jac_face_l(i,j) = jac_face(i,j,iface)
                        end do !j
                    end do !i

                end if

                !----------------------------------
                !  Compute face and element fluxes
                !----------------------------------
                flux_ql = ql
                flux_qr = qr
                flux_q = 0.5 * (flux_ql + flux_qr)

                !---------------------------------
                !  Do Gauss-Lobatto Integration
                !---------------------------------
                if (ftype==12) then
                    !need to backward-project the flux
                    call gather_element_2d_subface(flux_q(:,:),flux_lq(:,:),subface, ngl_i, ngl_j, plane_ij)
                    flux_q=flux_lq

                    !elemental flux needs to be computed differently for face 12 - will do it in the regular nc subroutine
                    flux_ql = 0
                end if

                do j=1,ngl_j
                    do i=1,ngl_i
                        !Pointers
                        if (ftype==21) then
                            il=imapr(1,i,j,iface)
                            jl=imapr(2,i,j,iface)
                            kl=imapr(3,i,j,iface)
                        else
                            il=imapl(1,i,j,iface)
                            jl=imapl(2,i,j,iface)
                            kl=imapl(3,i,j,iface)
                        end if
                        wq=jac_face_l(i,j)
                        ip=intma(il,jl,kl,iel)

                        !Store RHS
                        if(lpen) pen = penalty(i,j,iface)
                        rhs(ip)=rhs(ip) + wq*( (flux_q(i,j) - flux_ql(i,j)) + pen*(flux_qr(i,j) - flux_ql(i,j)) )

                    end do !i
                end do !j

                kk=kk+1
            end do
            jj=jj+1
        end do

    end do !iface

end subroutine