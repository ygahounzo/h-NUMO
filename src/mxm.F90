!-----------------------------------------------------------------------
!Original code from Paul Fischer from NEK5000
!Modified by F.X. Giraldo to hardwire loop unrolling
!-----------------------------------------------------------------------
subroutine mxm(a,n1,b,n2,c,n3)
    !-----------------------------------------------------------------------
    !
    !     NOTE -- 4x4 stencil call
    !
    !     Matrix-vector product routine.
    !     NOTE: Use assembly coded routine if available.
    !
    !     Computes C(N1,N3) = A(N1,N2)*B(N2,N3)
    !----------------------------------------------------------------------
    REAL A(N1,N2),B(N2,N3),C(N1,N3)
    integer n1, n2, n3

    !Use F90 Intrinsic Function
    !      c=matmul(a,b)
    !      return

    !Use BLAS Matrix-Matrix Multiply => similar to MATMUL
    !      alpha=1
    !      beta=0
    !      call dgemm('N','N',n1,n3,n1,alpha,a,n1,b,n2,beta,c,n1)
    !      return

    !Use Paul Fischer's original codes from NEK5000
    !      if (n2==2) then
    !         call mxm44_2(a,n1,b,n2,c,n3)
    !      else
    !         call mxm44_0(a,n1,b,n2,c,n3)
    !      endif
    !      return

    !Use ad hoc MXM by F.X. Giraldo which does hardwire unrolling.
    if (n2 <= 17) then
        select case (n2)
            case (2)
                call mxm02(a,n1,b,n2,c,n3)
            case (3)
                call mxm03(a,n1,b,n2,c,n3)
            case (4)
                call mxm04(a,n1,b,n2,c,n3)
            case (5)
                call mxm05(a,n1,b,n2,c,n3)
            case (6)
                call mxm06(a,n1,b,n2,c,n3)
            case (7)
                call mxm07(a,n1,b,n2,c,n3)
            case (8)
                call mxm08(a,n1,b,n2,c,n3)
            case (9)
                call mxm09(a,n1,b,n2,c,n3)
            case (10)
                call mxm10(a,n1,b,n2,c,n3)
            case (11)
                call mxm11(a,n1,b,n2,c,n3)
            case (12)
                call mxm12(a,n1,b,n2,c,n3)
            case (13)
                call mxm13(a,n1,b,n2,c,n3)
            case (14)
                call mxm14(a,n1,b,n2,c,n3)
            case (15)
                call mxm15(a,n1,b,n2,c,n3)
            case (16)
                call mxm16(a,n1,b,n2,c,n3)
            case (17)
                call mxm17(a,n1,b,n2,c,n3)
        end select
    else if (n2 > 17) then
        call mxm44_0(a,n1,b,n2,c,n3)
    end if

end subroutine mxm
!-----------------------------------------------------------------------
subroutine mxm44_2(a, m, b, k, c, n)
    real a(m,2), b(2,n), c(m,n)
    integer m, k, n, nresid, n1, j, i

    nresid = iand(n,3)
    n1 = n - nresid + 1

    do j=1,n-nresid,4
        do i=1,m
            c(i,j) = a(i,1)*b(1,j) &
                + a(i,2)*b(2,j)
            c(i,j+1) = a(i,1)*b(1,j+1) &
                + a(i,2)*b(2,j+1)
            c(i,j+2) = a(i,1)*b(1,j+2) &
                + a(i,2)*b(2,j+2)
            c(i,j+3) = a(i,1)*b(1,j+3) &
                + a(i,2)*b(2,j+3)
        enddo
    enddo
    if (nresid == 0) then
        return
    elseif (nresid == 1) then
        do i=1,m
            c(i,n) = a(i,1)*b(1,n) &
                + a(i,2)*b(2,n)
        enddo
    elseif (nresid == 2) then
        do i=1,m
            c(i,n-1) = a(i,1)*b(1,n-1) &
                + a(i,2)*b(2,n-1)
            c(i,n) = a(i,1)*b(1,n) &
                + a(i,2)*b(2,n)
        enddo
    else
        do i=1,m
            c(i,n-2) = a(i,1)*b(1,n-2) &
                + a(i,2)*b(2,n-2)
            c(i,n-1) = a(i,1)*b(1,n-1) &
                + a(i,2)*b(2,n-1)
            c(i,n) = a(i,1)*b(1,n) &
                + a(i,2)*b(2,n)
        enddo
    endif

end subroutine mxm44_2
!-----------------------------------------------------------------------
subroutine mxm02(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j)
        end do !i
    end do !j

end subroutine mxm02
!-----------------------------------------------------------------------
subroutine mxm03(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j)
        end do !i
    end do !j

end subroutine mxm03
!-----------------------------------------------------------------------
subroutine mxm04(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j)
        end do !i
    end do !j

end subroutine mxm04
!-----------------------------------------------------------------------
subroutine mxm05(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j)
        end do !i
    end do !j

end subroutine mxm05
!-----------------------------------------------------------------------
subroutine mxm06(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j)
        end do !i
    end do !j

end subroutine mxm06
!-----------------------------------------------------------------------
subroutine mxm07(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j) &
                + a(i,7)*b(7,j)
        end do !i
    end do !j

end subroutine mxm07
!-----------------------------------------------------------------------
subroutine mxm08(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j) &
                + a(i,7)*b(7,j) + a(i,8)*b(8,j)
        end do !i
    end do !j

end subroutine mxm08
!-----------------------------------------------------------------------
subroutine mxm09(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j) &
                + a(i,7)*b(7,j) + a(i,8)*b(8,j) &
                + a(i,9)*b(9,j)
        end do !i
    end do !j

end subroutine mxm09
!-----------------------------------------------------------------------
subroutine mxm10(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j) &
                + a(i,7)*b(7,j) + a(i,8)*b(8,j) &
                + a(i,9)*b(9,j) + a(i,10)*b(10,j)
        end do !i
    end do !j

end subroutine mxm10
!-----------------------------------------------------------------------
subroutine mxm11(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j) &
                + a(i,7)*b(7,j) + a(i,8)*b(8,j) &
                + a(i,9)*b(9,j) + a(i,10)*b(10,j) &
                + a(i,11)*b(11,j)
        end do !i
    end do !j

end subroutine mxm11
!-----------------------------------------------------------------------
subroutine mxm12(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j) &
                + a(i,7)*b(7,j) + a(i,8)*b(8,j) &
                + a(i,9)*b(9,j) + a(i,10)*b(10,j) &
                + a(i,11)*b(11,j) + a(i,12)*b(12,j)
        end do !i
    end do !j

end subroutine mxm12
!-----------------------------------------------------------------------
subroutine mxm13(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j) &
                + a(i,7)*b(7,j) + a(i,8)*b(8,j) &
                + a(i,9)*b(9,j) + a(i,10)*b(10,j) &
                + a(i,11)*b(11,j) + a(i,12)*b(12,j) &
                + a(i,13)*b(13,j)
        end do !i
    end do !j

end subroutine mxm13
!-----------------------------------------------------------------------
subroutine mxm14(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j) &
                + a(i,7)*b(7,j) + a(i,8)*b(8,j) &
                + a(i,9)*b(9,j) + a(i,10)*b(10,j) &
                + a(i,11)*b(11,j) + a(i,12)*b(12,j) &
                + a(i,13)*b(13,j) + a(i,14)*b(14,j)
        end do !i
    end do !j

end subroutine mxm14
!-----------------------------------------------------------------------
subroutine mxm15(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j) &
                + a(i,7)*b(7,j) + a(i,8)*b(8,j) &
                + a(i,9)*b(9,j) + a(i,10)*b(10,j) &
                + a(i,11)*b(11,j) + a(i,12)*b(12,j) &
                + a(i,13)*b(13,j) + a(i,14)*b(14,j) &
                + a(i,15)*b(15,j)
        end do !i
    end do !j

end subroutine mxm15
!-----------------------------------------------------------------------
subroutine mxm16(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j) &
                + a(i,7)*b(7,j) + a(i,8)*b(8,j) &
                + a(i,9)*b(9,j) + a(i,10)*b(10,j) &
                + a(i,11)*b(11,j) + a(i,12)*b(12,j) &
                + a(i,13)*b(13,j) + a(i,14)*b(14,j) &
                + a(i,15)*b(15,j) + a(i,16)*b(16,j)
        end do !i
    end do !j

end subroutine mxm16
!-----------------------------------------------------------------------
subroutine mxm17(a, n1, b, n2, c, n3)
    real a(n1,n2), b(n2,n3), c(n1,n3)
    integer n1, n2, n3, i, j

    do j=1,n3
        do i=1,n1
            c(i,j)=a(i,1)*b(1,j) + a(i,2)*b(2,j) &
                + a(i,3)*b(3,j) + a(i,4)*b(4,j) &
                + a(i,5)*b(5,j) + a(i,6)*b(6,j) &
                + a(i,7)*b(7,j) + a(i,8)*b(8,j) &
                + a(i,9)*b(9,j) + a(i,10)*b(10,j) &
                + a(i,11)*b(11,j) + a(i,12)*b(12,j) &
                + a(i,13)*b(13,j) + a(i,14)*b(14,j) &
                + a(i,15)*b(15,j) + a(i,16)*b(16,j) &
                + a(i,17)*b(17,j)
        end do !i
    end do !j

end subroutine mxm17
!-----------------------------------------------------------------------
subroutine mxm44_0(a, m, b, k, c, n)
    real a(m,k), b(k,n), c(m,n)
    real s11, s12, s13, s14, s21, s22, s23, s24
    real s31, s32, s33, s34, s41, s42, s43, s44
    integer k, m, n, i, j, mresid, nresid, m1, n1, l

    !matrix multiply with a 4x4 stencil
    mresid = iand(m,3)
    nresid = iand(n,3)
    m1 = m - mresid + 1
    n1 = n - nresid + 1

    do i=1,m-mresid,4
        do j=1,n-nresid,4
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0
            s43 = 0.0d0
            s14 = 0.0d0
            s24 = 0.0d0
            s34 = 0.0d0
            s44 = 0.0d0
            do l=1,k
                s11 = s11 + a(i,l)*b(l,j)
                s12 = s12 + a(i,l)*b(l,j+1)
                s13 = s13 + a(i,l)*b(l,j+2)
                s14 = s14 + a(i,l)*b(l,j+3)

                s21 = s21 + a(i+1,l)*b(l,j)
                s22 = s22 + a(i+1,l)*b(l,j+1)
                s23 = s23 + a(i+1,l)*b(l,j+2)
                s24 = s24 + a(i+1,l)*b(l,j+3)

                s31 = s31 + a(i+2,l)*b(l,j)
                s32 = s32 + a(i+2,l)*b(l,j+1)
                s33 = s33 + a(i+2,l)*b(l,j+2)
                s34 = s34 + a(i+2,l)*b(l,j+3)

                s41 = s41 + a(i+3,l)*b(l,j)
                s42 = s42 + a(i+3,l)*b(l,j+1)
                s43 = s43 + a(i+3,l)*b(l,j+2)
                s44 = s44 + a(i+3,l)*b(l,j+3)
            enddo
            c(i,j)     = s11 
            c(i,j+1)   = s12
            c(i,j+2)   = s13
            c(i,j+3)   = s14

            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 

            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42

            c(i+1,j+2) = s23
            c(i+2,j+2) = s33
            c(i+3,j+2) = s43

            c(i+1,j+3) = s24
            c(i+2,j+3) = s34
            c(i+3,j+3) = s44
        enddo
        ! Residual when n is not multiple of 4
        if (nresid /= 0) then
            if (nresid == 1) then
                s11 = 0.0d0
                s21 = 0.0d0
                s31 = 0.0d0
                s41 = 0.0d0
                do l=1,k
                    s11 = s11 + a(i,l)*b(l,n)
                    s21 = s21 + a(i+1,l)*b(l,n)
                    s31 = s31 + a(i+2,l)*b(l,n)
                    s41 = s41 + a(i+3,l)*b(l,n)
                enddo
                c(i,n)     = s11
                c(i+1,n)   = s21
                c(i+2,n)   = s31
                c(i+3,n)   = s41
            elseif (nresid == 2) then
                s11 = 0.0d0
                s21 = 0.0d0
                s31 = 0.0d0
                s41 = 0.0d0
                s12 = 0.0d0
                s22 = 0.0d0
                s32 = 0.0d0
                s42 = 0.0d0
                do l=1,k
                    s11 = s11 + a(i,l)*b(l,j)
                    s12 = s12 + a(i,l)*b(l,j+1)

                    s21 = s21 + a(i+1,l)*b(l,j)
                    s22 = s22 + a(i+1,l)*b(l,j+1)

                    s31 = s31 + a(i+2,l)*b(l,j)
                    s32 = s32 + a(i+2,l)*b(l,j+1)

                    s41 = s41 + a(i+3,l)*b(l,j)
                    s42 = s42 + a(i+3,l)*b(l,j+1)
                enddo
                c(i,j)     = s11
                c(i,j+1)   = s12

                c(i+1,j)   = s21
                c(i+2,j)   = s31
                c(i+3,j)   = s41

                c(i+1,j+1) = s22
                c(i+2,j+1) = s32
                c(i+3,j+1) = s42
            else
                s11 = 0.0d0
                s21 = 0.0d0
                s31 = 0.0d0
                s41 = 0.0d0
                s12 = 0.0d0
                s22 = 0.0d0
                s32 = 0.0d0
                s42 = 0.0d0
                s13 = 0.0d0
                s23 = 0.0d0
                s33 = 0.0d0
                s43 = 0.0d0
                do l=1,k
                    s11 = s11 + a(i,l)*b(l,j)
                    s12 = s12 + a(i,l)*b(l,j+1)
                    s13 = s13 + a(i,l)*b(l,j+2)

                    s21 = s21 + a(i+1,l)*b(l,j)
                    s22 = s22 + a(i+1,l)*b(l,j+1)
                    s23 = s23 + a(i+1,l)*b(l,j+2)

                    s31 = s31 + a(i+2,l)*b(l,j)
                    s32 = s32 + a(i+2,l)*b(l,j+1)
                    s33 = s33 + a(i+2,l)*b(l,j+2)

                    s41 = s41 + a(i+3,l)*b(l,j)
                    s42 = s42 + a(i+3,l)*b(l,j+1)
                    s43 = s43 + a(i+3,l)*b(l,j+2)
                enddo
                c(i,j)     = s11
                c(i+1,j)   = s21
                c(i+2,j)   = s31
                c(i+3,j)   = s41
                c(i,j+1)   = s12
                c(i+1,j+1) = s22
                c(i+2,j+1) = s32
                c(i+3,j+1) = s42
                c(i,j+2)   = s13
                c(i+1,j+2) = s23
                c(i+2,j+2) = s33
                c(i+3,j+2) = s43
            endif
        endif
    enddo

    ! Residual when m is not multiple of 4
    if (mresid == 0) then
        return
    elseif (mresid == 1) then
        do j=1,n-nresid,4
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            do l=1,k
                s11 = s11 + a(m,l)*b(l,j)
                s12 = s12 + a(m,l)*b(l,j+1)
                s13 = s13 + a(m,l)*b(l,j+2)
                s14 = s14 + a(m,l)*b(l,j+3)
            enddo
            c(m,j)     = s11
            c(m,j+1)   = s12
            c(m,j+2)   = s13
            c(m,j+3)   = s14
        enddo
        ! mresid is 1, check nresid
        if (nresid == 0) then
            return
        elseif (nresid == 1) then
            s11 = 0.0d0
            do l=1,k
                s11 = s11 + a(m,l)*b(l,n)
            enddo
            c(m,n) = s11
            return
        elseif (nresid == 2) then
            s11 = 0.0d0
            s12 = 0.0d0
            do l=1,k
                s11 = s11 + a(m,l)*b(l,n-1)
                s12 = s12 + a(m,l)*b(l,n)
            enddo
            c(m,n-1) = s11
            c(m,n) = s12
            return
        else
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            do l=1,k
                s11 = s11 + a(m,l)*b(l,n-2)
                s12 = s12 + a(m,l)*b(l,n-1)
                s13 = s13 + a(m,l)*b(l,n)
            enddo
            c(m,n-2) = s11
            c(m,n-1) = s12
            c(m,n) = s13
            return
        endif          
    elseif (mresid == 2) then
        do j=1,n-nresid,4
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s21 = 0.0d0
            s22 = 0.0d0
            s23 = 0.0d0
            s24 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-1,l)*b(l,j)
                s12 = s12 + a(m-1,l)*b(l,j+1)
                s13 = s13 + a(m-1,l)*b(l,j+2)
                s14 = s14 + a(m-1,l)*b(l,j+3)

                s21 = s21 + a(m,l)*b(l,j)
                s22 = s22 + a(m,l)*b(l,j+1)
                s23 = s23 + a(m,l)*b(l,j+2)
                s24 = s24 + a(m,l)*b(l,j+3)
            enddo
            c(m-1,j)   = s11
            c(m-1,j+1) = s12
            c(m-1,j+2) = s13
            c(m-1,j+3) = s14
            c(m,j)     = s21
            c(m,j+1)   = s22
            c(m,j+2)   = s23
            c(m,j+3)   = s24
        enddo
        ! mresid is 2, check nresid
        if (nresid == 0) then
            return
        elseif (nresid == 1) then
            s11 = 0.0d0
            s21 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-1,l)*b(l,n)
                s21 = s21 + a(m,l)*b(l,n)
            enddo
            c(m-1,n) = s11
            c(m,n) = s21
            return
        elseif (nresid == 2) then
            s11 = 0.0d0
            s21 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-1,l)*b(l,n-1)
                s12 = s12 + a(m-1,l)*b(l,n)
                s21 = s21 + a(m,l)*b(l,n-1)
                s22 = s22 + a(m,l)*b(l,n)
            enddo
            c(m-1,n-1) = s11
            c(m-1,n)   = s12
            c(m,n-1)   = s21
            c(m,n)     = s22
            return
        else
            s11 = 0.0d0
            s21 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-1,l)*b(l,n-2)
                s12 = s12 + a(m-1,l)*b(l,n-1)
                s13 = s13 + a(m-1,l)*b(l,n)
                s21 = s21 + a(m,l)*b(l,n-2)
                s22 = s22 + a(m,l)*b(l,n-1)
                s23 = s23 + a(m,l)*b(l,n)
            enddo
            c(m-1,n-2) = s11
            c(m-1,n-1) = s12
            c(m-1,n)   = s13
            c(m,n-2)   = s21
            c(m,n-1)   = s22
            c(m,n)     = s23
            return
        endif
    else
        ! mresid is 3
        do j=1,n-nresid,4
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0

            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0

            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0

            s14 = 0.0d0
            s24 = 0.0d0
            s34 = 0.0d0

            do l=1,k
                s11 = s11 + a(m-2,l)*b(l,j)
                s12 = s12 + a(m-2,l)*b(l,j+1)
                s13 = s13 + a(m-2,l)*b(l,j+2)
                s14 = s14 + a(m-2,l)*b(l,j+3)

                s21 = s21 + a(m-1,l)*b(l,j)
                s22 = s22 + a(m-1,l)*b(l,j+1)
                s23 = s23 + a(m-1,l)*b(l,j+2)
                s24 = s24 + a(m-1,l)*b(l,j+3)

                s31 = s31 + a(m,l)*b(l,j)
                s32 = s32 + a(m,l)*b(l,j+1)
                s33 = s33 + a(m,l)*b(l,j+2)
                s34 = s34 + a(m,l)*b(l,j+3)
            enddo
            c(m-2,j)   = s11
            c(m-2,j+1) = s12
            c(m-2,j+2) = s13
            c(m-2,j+3) = s14

            c(m-1,j)   = s21
            c(m-1,j+1) = s22
            c(m-1,j+2) = s23
            c(m-1,j+3) = s24

            c(m,j)     = s31
            c(m,j+1)   = s32
            c(m,j+2)   = s33
            c(m,j+3)   = s34
        enddo
        ! mresid is 3, check nresid
        if (nresid == 0) then
            return
        elseif (nresid == 1) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-2,l)*b(l,n)
                s21 = s21 + a(m-1,l)*b(l,n)
                s31 = s31 + a(m,l)*b(l,n)
            enddo
            c(m-2,n) = s11
            c(m-1,n) = s21
            c(m,n) = s31
            return
        elseif (nresid == 2) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-2,l)*b(l,n-1)
                s12 = s12 + a(m-2,l)*b(l,n)
                s21 = s21 + a(m-1,l)*b(l,n-1)
                s22 = s22 + a(m-1,l)*b(l,n)
                s31 = s31 + a(m,l)*b(l,n-1)
                s32 = s32 + a(m,l)*b(l,n)
            enddo
            c(m-2,n-1) = s11
            c(m-2,n)   = s12
            c(m-1,n-1) = s21
            c(m-1,n)   = s22
            c(m,n-1)   = s31
            c(m,n)     = s32
            return
        else
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0
            do l=1,k
                s11 = s11 + a(m-2,l)*b(l,n-2)
                s12 = s12 + a(m-2,l)*b(l,n-1)
                s13 = s13 + a(m-2,l)*b(l,n)
                s21 = s21 + a(m-1,l)*b(l,n-2)
                s22 = s22 + a(m-1,l)*b(l,n-1)
                s23 = s23 + a(m-1,l)*b(l,n)
                s31 = s31 + a(m,l)*b(l,n-2)
                s32 = s32 + a(m,l)*b(l,n-1)
                s33 = s33 + a(m,l)*b(l,n)
            enddo
            c(m-2,n-2) = s11
            c(m-2,n-1) = s12
            c(m-2,n)   = s13
            c(m-1,n-2) = s21
            c(m-1,n-1) = s22
            c(m-1,n)   = s23
            c(m,n-2)   = s31
            c(m,n-1)   = s32
            c(m,n)     = s33
            return
        endif
    endif

end subroutine mxm44_0
!-----------------------------------------------------------------------
