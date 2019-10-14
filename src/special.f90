module special

! Special functions:
! This module offers some special functions such as
!- spherical Bessel functions
!- zeros of real Bessel functions
!- Bessel and Hankel functions of complex argument (wrapping the AMOS library
!   from http://netlib.org/amos)
!- Airy functions Ai and Bi

! Note on the AMOS wrappers: the AMOS routines generally offer the possibility
! of exponential scaling to remove exponential asymptotic behaviour.
! We do not expose the corresponding switches to the user, however, it is
! straight forward to add those.

! All references labelled with "NIST" refer to the NIST Handbook of
! Mathematical Functions (2010)

use types, only: dp
use constants, only: pi, i_
use utils, only: stop_error
use optimize, only: bisect
use amos
implicit none

private
public bessel_jn_zeros, spherical_bessel_jn, spherical_bessel_yn, &
       spherical_bessel_jn_zeros, &
       besselj, bessely, hankel1, hankel2, besseli, besselk, &  ! note the missing underscores in Bessel functions
       airyai, dairyai, airybi, dairybi  ! Airy functions Ai and Bi and their derivatives

! Bessel J function (first kind)
interface besselj
    module procedure besseljn_real
    module procedure besseljv_real
    module procedure besseljv_complex
end interface besselj

! Bessel Y function (second kind)
interface bessely
    module procedure besselyn_real
    module procedure besselyv_real
    module procedure besselyv_complex
end interface bessely

! Hankel H function of first kind
interface hankel1
    module procedure hankel1n_real
    module procedure hankel1v_real
    module procedure hankel1v_complex
end interface hankel1

! Hankel H function of second kind
interface hankel2
    module procedure hankel2n_real
    module procedure hankel2v_real
    module procedure hankel2v_complex
end interface hankel2

! modified Bessel function of first kind
interface besseli
    module procedure besseliv_real
    module procedure besseliv_complex
end interface besseli

! modified Bessel function of second kind
interface besselk
    module procedure besselkv_real
    module procedure besselkv_complex
end interface besselk

contains

function besseljn_real(order, x) result(res)
    integer, intent(in) :: order
    real(dp), intent(in) :: x
    real(dp) :: res

    ! convenience function that calls Bessel functions of integer order
    ! (intrinsic to F2008)

    if(abs(order) == 0) then
        res = bessel_j0(x)
    elseif(abs(order) == 1) then
        res = bessel_j1(x)
    else
        res = bessel_jn(abs(order), x)
    endif

    ! apply negative order NIST eq. (10.4.1) if needed:
    res = reflect(order)*res
end function besseljn_real


function besseljv_real(order, x) result(res)
    real(dp), intent(in) :: order, x
    real(dp) :: res
    complex(dp) :: resCmplx

    ! here, we do rely in an informed user: if the user wanted Bessel functions
    ! of integer order n, he/she should call bessel_jn() (as of F2008) directly.
    ! Hence, here we just convert the argument to a complex number and call the
    ! AMOS routine
    ! TODO: on the long run, wrap the cephes routines!

    resCmplx = besseljv_complex(order, cmplx(x, 0.0_dp, kind=dp))  ! negative order is handled in the routine called
    res = real(resCmplx)  ! Bessel functions of real arument are real
end function besseljv_real


function besseljv_complex(order, z) result(res)
    real(dp), intent(in) :: order
    complex(dp), intent(in) :: z
    complex(dp) :: res

    integer, parameter :: scaling=1, length=1
    real(dp) :: resReal(length), resImag(length), callorder
    integer :: underflows, ierr

    ! make sure AMOS zbesj() is called with postive order
    callorder = abs(order)

    call zbesj(real(z), aimag(z), callorder, scaling, length, resReal, resImag, underflows, ierr)
    if(ierr > 0) then
        if(ierr == 1) then
            print *, "IERR=1, INPUT ERROR - NO COMPUTATION"
        elseif(ierr == 2) then
            print *, "IERR=2, OVERFLOW - NO COMPUTATION, AIMAG(Z) TOO LARGE ON KODE=1"
        elseif(ierr == 3) then
            print *, "IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE"
            print *, "BUT LOSSES OF SIGNIFCANCE BY ARGUMENT"
            print *, "REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY"
        elseif(ierr == 4) then
            print *, "IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION BECAUSE"
            print *, "OF COMPLETE LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION"
        elseif(ierr == 5) then
            print *, "IERR=5, ERROR - NO COMPUTATION ALGORITHM TERMINATION CONDITION NOT MET"
        endif
        call stop_error("besselj: zbesj error")
    endif

    ! take care of negative order if needed:
    if(order < 0.0_dp) then  ! use NIST eq. (10.4.7):
        res = exp(-pi*i_*callorder)*(resReal(1) + i_*resImag(1)) &
            & + i_*sin(pi*callorder)*hankel1(callorder, z)
    else
        res = resReal(1) +i_*resImag(1)
    endif
end function besseljv_complex


function besselyn_real(order, x) result(res)
    integer, intent(in) :: order
    real(dp), intent(in) :: x
    real(dp) :: res

    ! convenience function that calls Bessel functions of integer order
    ! (intrinsic to F2008)

    if(abs(order) == 0) then
        res = bessel_y0(x)
    elseif(abs(order) == 1) then
        res = bessel_y1(x)
    else
        res = bessel_yn(abs(order), x)
    endif
    res = reflect(order)*res  ! apply negative order eq. (10.4.1) from NIST
end function besselyn_real


function besselyv_real(order, x) result(res)
    real(dp), intent(in) :: order, x
    real(dp) :: res
    complex(dp) :: resCmplx

    resCmplx = besselyv_complex(order, x+0*i_)  ! neg. order handled in the routine called
    res = real(resCmplx)
end function besselyv_real


function besselyv_complex(order, z) result(res)
    real(dp), intent(in) :: order
    complex(dp), intent(in) :: z
    complex(dp) :: res

    integer, parameter :: scaling=1, length=1
    real(dp) :: workr(length), worki(length), resReal(length), resImag(length), callorder
    integer :: underflows, ierr

    ! make sure AMOS zbesy() is called with positive order:
    callorder = abs(order)

    call zbesy(real(z), aimag(z), callorder, scaling, length, resReal, resImag, underflows, workr, worki, ierr)
    if(ierr > 0) then
        if(ierr == 1) then
            print *, "IERR=1, INPUT ERROR - NO COMPUTATION"
        elseif(ierr == 2) then
            print *, "IERR=2, OVERFLOW - NO COMPUTATION, FNU IS TOO LARGE OR CABS(Z)"
            print *, "IS TOO SMALL OR BOTH"
        elseif(ierr == 3) then
            print *, "IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE BUT LOSSES OF SIGNIFCANCE"
            print *, "BY ARGUMENT REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY"
        elseif(ierr == 4) then
            print *, "IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION BECAUSE OF COMPLETE"
            print *, "LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION"
        elseif(ierr == 5) then
            print *, "IERR=5, ERROR - NO COMPUTATION ALGORITHM TERMINATION CONDITION NOT MET"
        endif
        call stop_error("bessely: zbesy error")
    endif

    if(order < 0.0_dp) then  ! apply NIST eq. (10.4.7)
        res = exp(-callorder*pi*i_)*(resReal(1) + i_*resImag(1)) + &
            & sin(pi*callorder)*hankel1(callorder, z)
    else
        res = resReal(1) + i_*resImag(1)
    endif
end function besselyv_complex


function hankel1n_real(order, x) result(res)
    integer, intent(in) :: order
    real(dp), intent(in) :: x
    complex(dp) :: res

    ! There are two linearly independent solutions to Bessel's equation.
    ! Thus, the two Hankel functions can be written as superpositions of
    ! J and Y, which are taken to be the two linearly independent ones.
    ! NIST Handbook (10.4.3)

    ! negative order is taken care of according to NIST eq. (10.4.2)
    if(order == 0) then
        res = bessel_j0(x) + i_*bessel_y0(x)
    elseif(abs(order) == 1) then
        res = reflect(order)*(bessel_j1(x) + i_*bessel_y1(x))
    else
        res = reflect(order)*(bessel_jn(order, x) + i_*bessel_yn(order, x))
    endif
end function hankel1n_real


function hankel1v_real(order, x) result(res)
    real(dp), intent(in) :: order, x
    complex(dp) :: res

    res = hankel1v_complex(order, cmplx(x, 0.0_dp, kind=dp))
end function hankel1v_real


function hankel1v_complex(order, z) result(res)
    real(dp), intent(in) :: order
    complex(dp), intent(in) :: z
    complex(dp) :: res

    integer, parameter :: scaling=1, length=1, hankelkind=1
    real(dp) :: resReal(length), resImag(length), callorder
    integer :: underflows, ierr

    callorder = abs(order)

    call zbesh(real(z), aimag(z), callorder, scaling, hankelkind, length, resReal, resImag, underflows, ierr)
    if(ierr > 0) then
        if(ierr == 1) then
            print *, "IERR=1, INPUT ERROR - NO COMPUTATION"
        elseif(ierr == 2) then
            print *, "IERR=2, OVERFLOW - NO COMPUTATION, FNU TOO LARGE OR CABS(Z) TOO SMALL OR BOTH"
        elseif(ierr == 3) then
            print *, "IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE BUT LOSSES OF SIGNIFCANCE"
            print *, "BY ARGUMENT REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY"
        elseif(ierr == 4) then
            print *, "IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION BECAUSE OF COMPLETE"
            print *, "LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION"
        elseif(ierr == 5) then
            print *, "IERR=5, ERROR - NO COMPUTATION ALGORITHM TERMINATION CONDITION NOT MET"
        endif
        call stop_error("hankel1: zbesh error")
    endif

    ! take care of negative orders if needed:
    res = Hrotate(1, order)*(resReal(1) + i_*resImag(1))
end function hankel1v_complex


function hankel2n_real(order, x) result(res)
    integer, intent(in) :: order
    real(dp), intent(in) :: x
    complex(dp) :: res

    ! convenience function that calls Bessel functions of integer order
    ! (intrinsic to F2008)

    if(abs(order) == 0) then
        res = bessel_j0(x) - i_*bessel_y0(x)
    elseif(abs(order) == 1) then
        res = reflect(order)*(bessel_j1(x) - i_*bessel_y1(x))
    else
        res = reflect(order)*(bessel_jn(order, x) - i_*bessel_yn(order, x))
    endif
end function hankel2n_real


function hankel2v_real(order, x) result(res)
    real(dp), intent(in) :: order, x
    complex(dp) :: res

    res = hankel2v_complex(order, cmplx(x, 0.0_dp, kind=dp))
end function hankel2v_real


function hankel2v_complex(order, z) result(res)
    real(dp), intent(in) :: order
    complex(dp), intent(in) :: z
    complex(dp) :: res

    integer, parameter :: scaling=1, length=1, hankelkind=2
    real(dp) :: resReal(length), resImag(length), callorder
    integer :: underflows, ierr

    callorder = abs(order)

    call zbesh(real(z), aimag(z), callorder, scaling, hankelkind, length, resReal, resImag, underflows, ierr)
    if(ierr > 0) then
        if(ierr == 1) then
            print *, "IERR=1, INPUT ERROR - NO COMPUTATION"
        elseif(ierr == 2) then
            print *, "IERR=2, OVERFLOW - NO COMPUTATION, FNU TOO LARGE OR CABS(Z) TOO SMALL OR BOTH"
        elseif(ierr == 3) then
            print *, "IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE BUT LOSSES OF SIGNIFCANCE"
            print *, "BY ARGUMENT REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY"
        elseif(ierr == 4) then
            print *, "IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION BECAUSE OF COMPLETE"
            print *, "LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION"
        elseif(ierr == 5) then
            print *, "IERR=5, ERROR - NO COMPUTATION ALGORITHM TERMINATION CONDITION NOT MET"
        endif
        call stop_error("hankel2: zbesh error")
    endif

    ! take care of neg. orders by Hrotate() if needed:
    res = Hrotate(2, order)*(resReal(1) + i_*resImag(1))
end function hankel2v_complex


function besseliv_real(order, x) result(res)
    real(dp), intent(in) :: order, x
    real(dp) :: res
    complex(dp) :: resCmplx

    resCmplx = besseliv_complex(order, x+0*i_)
    res = real(resCmplx)
end function besseliv_real


function besseliv_complex(order, z) result(res)
    real(dp), intent(in) :: order
    complex(dp), intent(in) :: z
    complex(dp) :: res

    integer, parameter :: scaling=1, length=1
    real(dp) :: resReal(length), resImag(length), callorder
    integer :: underflows, ierr

    callorder = abs(order)

    call zbesi(real(z), aimag(z), callorder, scaling, length, resReal, resImag, underflows, ierr)
    if(ierr > 0) then
        if(ierr == 1) then
            print *, "IERR=1, INPUT ERROR - NO COMPUTATION"
        elseif(ierr == 2) then
            print *, "IERR=2, OVERFLOW - NO COMPUTATION, REAL(Z) TOO LARGE ON KODE=1"
        elseif(ierr == 3) then
            print *, "IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE BUT LOSSES OF SIGNIFCANCE"
            print *, "BY ARGUMENT REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY"
        elseif(ierr == 4) then
            print *, "IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION BECAUSE OF COMPLETE"
            print *, "LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION"
        elseif(ierr == 5) then
            print *, "IERR=5, ERROR - NO COMPUTATION ALGORITHM TERMINATION CONDITION NOT MET"
        endif
        call stop_error("besseli: zbesi error")
    endif

    ! negative order from NIST Handbook, (10.27.2)
    if(order < 0.0_dp) then
        res = (resReal(1) + i_*resImag(1)) &
            & + 2.0_dp/pi*sin(callorder*pi)*besselk(callorder, z)
    else
        res = resReal(1) + i_*resImag(1)
    endif
end function besseliv_complex

function besselkv_real(order, x) result(res)
    real(dp), intent(in) :: order, x
    real(dp) :: res
    complex(dp) :: resCmplx

    resCmplx = besselkv_complex(order, x+0*i_)
    res = real(resCmplx)
end function besselkv_real


function besselkv_complex(order, z) result(res)
    real(dp), intent(in) :: order
    complex(dp), intent(in) :: z
    complex(dp) :: res

    integer, parameter :: scaling=1, length=1
    real(dp) :: resReal(length), resImag(length), callorder
    integer :: underflows, ierr

    ! K_{-v} = K_v; NIST Handbook (10.27.3)
    callorder = abs(order)

    call zbesk(real(z), aimag(z), callorder, scaling, length, resReal, resImag, underflows, ierr)
    if(ierr > 0) then
        if(ierr == 1) then
            print *, "IERR=1, INPUT ERROR - NO COMPUTATION"
        elseif(ierr == 2) then
            print *, "IERR=2, OVERFLOW - NO COMPUTATION, FNU IS TOO LARGE OR CABS(Z) IS TOO SMALL"
            print *, "OR BOTH"
        elseif(ierr == 3) then
            print *, "IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE BUT LOSSES OF SIGNIFCANCE"
            print *, "BY ARGUMENT REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY"
        elseif(ierr == 4) then
            print *, "IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION BECAUSE OF COMPLETE"
            print *, "LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION"
        elseif(ierr == 5) then
            print *, "IERR=5, ERROR - NO COMPUTATION ALGORITHM TERMINATION CONDITION NOT MET"
        endif
        call stop_error("besselk: zbesk error")
    endif
    res = resReal(1) + i_*resImag(1)
end function besselkv_complex


function bessel_j0_zeros(nzeros, eps) result(zeros)
! Calculates zeros of bessel_j0().
! This function is then used in bessel_jn_zeros() to calculate zeros of j_n(x)
! for n > 0 by simply using the zeros of j_{n-1}(x), as they must lie between.
integer, intent(in) :: nzeros
real(dp), intent(in) :: eps
! zeros(i) is the i-th zero of bessel_j0(x)
real(dp) :: zeros(nzeros)
real(dp) :: points(0:nzeros)
integer :: n, j
! The zeros of j0(x) must lie between the 'points', as can be shown by:
!
!   In [97]: n = 100000; arange(1, n+1) * pi - scipy.special.jn_zeros(0, n)
!   Out[97]:
!   array([ 0.7367671 ,  0.7631072 ,  0.77105005, ...,  0.78539777,
!           0.78539777,  0.78539777])
!
! For large "x", the asymptotic is j0(x) ~ cos(x - ...), so there it holds and
! for small 'x', it just happens to be the case numerically (see above).

points = [ (n*pi, n = 0, nzeros) ]
do j = 1, nzeros
    zeros(j) = bisect(f, points(j-1), points(j), eps)
end do

contains

real(dp) function f(x)
real(dp), intent(in) :: x
f = bessel_j0(x)
end function

end function


function bessel_jn_zeros(nmax, nzeros, eps) result(zeros)
! Calculates 'nzeros' zeros of bessel_jn() for all n=0, 1, ..., nmax
! It uses the fact that zeros of j_{n-1}(x) lie between zeros of j_n(x) and
! uses bisection to calculate them.
integer, intent(in) :: nmax, nzeros
real(dp), intent(in) :: eps
! zeros(i, n) is the i-th zero of bessel_jn(n, x)
real(dp) :: zeros(nzeros, 0:nmax)
! points holds all zeros of j_{n-1} needed for j_n
real(dp) :: points(nmax+nzeros)
integer :: n, j
points = bessel_j0_zeros(nmax + nzeros, eps)
zeros(:, 0) = points(:nzeros)
do n = 1, nmax
    do j = 1, nmax + nzeros - n
        points(j) = bisect(f, points(j), points(j+1), eps)
    end do
    zeros(:, n) = points(:nzeros)
end do

contains

real(dp) function f(x)
real(dp), intent(in) :: x
f = bessel_jn(n, x)
end function

end function


real(dp) function spherical_bessel_jn(n, x) result(r)
integer, intent(in) :: n
real(dp), intent(in) :: x
integer :: nm
real(dp) :: sj(0:n), dj(0:n)
call sphj(n, x, nm, sj, dj)
if (nm /= n) call stop_error("spherical_bessel_jn: sphj didn't converge")
r = sj(n)
end function

real(dp) function spherical_bessel_yn(n, x) result(r)
integer, intent(in) :: n
real(dp), intent(in) :: x
integer :: nm
real(dp) :: sy(0:n), dy(0:n)
call sphy(n, x, nm, sy, dy)
if (nm /= n) call stop_error("spherical_bessel_yn: sphy didn't converge")
r = sy(n)
end function


function spherical_bessel_jn_zeros(nmax, nzeros, eps) result(zeros)
! Calculates 'nzeros' zeros of spherical_bessel_jn() for all n=0, 1, ..., nmax
! It uses the fact that zeros of j_{n-1}(x) lie between zeros of j_n(x) and
! uses bisection to calculate them.
integer, intent(in) :: nmax, nzeros
real(dp), intent(in) :: eps
! zeros(i, n) is the i-th zero of spherical_bessel_jn(n, x)
real(dp) :: zeros(nzeros, 0:nmax)
! points holds all zeros of j_{n-1} needed for j_n
real(dp) :: points(nmax+nzeros)
integer :: n, j
points = [ (n*pi, n = 1, nmax + nzeros) ]
zeros(:, 0) = points(:nzeros)
do n = 1, nmax
    do j = 1, nmax + nzeros - n
        points(j) = bisect(f, points(j), points(j+1), eps)
    end do
    zeros(:, n) = points(:nzeros)
end do

contains

    real(dp) function f(x)
        real(dp), intent(in) :: x
        f = spherical_bessel_jn(n, x)
    end function f

end function


! The SPHJ, SPHY, MSTA1, MSTA2 routines below are taken from SciPy's specfun.f.
! Authors: Shanjie Zhang and Jianming Jin
! Copyrighted but permission granted to use code in programs.
        SUBROUTINE SPHJ(N,X,NM,SJ,DJ)
!       =======================================================
!       Purpose: Compute spherical Bessel functions jn(x) and
!                their derivatives
!       Input :  x --- Argument of jn(x)
!                n --- Order of jn(x)  ( n = 0,1,… )
!       Output:  SJ(n) --- jn(x)
!                DJ(n) --- jn'(x)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       =======================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT integer (I-N)
        DIMENSION SJ(0:N),DJ(0:N)
        NM=N
        IF (DABS(X).LT.1.0D-100) THEN
           DO 10 K=0,N
              SJ(K)=0.0D0
              DJ(K)=0.0D0
10         END DO
           SJ(0)=1.0D0
           IF (N.GT.0) THEN
              DJ(1)=.3333333333333333D0
           ENDIF
           RETURN
        ENDIF
        SJ(0)=DSIN(X)/X
        DJ(0)=(DCOS(X)-DSIN(X)/X)/X
        IF (N.LT.1) THEN
           RETURN
        ENDIF
        SJ(1)=(SJ(0)-DCOS(X))/X
        IF (N.GE.2) THEN
           SA=SJ(0)
           SB=SJ(1)
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F=0.0D0
           F0=0.0D0
           F1=1.0D0-100
           DO 15 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F1/X-F0
              IF (K.LE.NM) SJ(K)=F
              F0=F1
              F1=F
15         END DO
           CS=0.0D0
           IF (DABS(SA).GT.DABS(SB)) CS=SA/F
           IF (DABS(SA).LE.DABS(SB)) CS=SB/F0
           DO 20 K=0,NM
              SJ(K)=CS*SJ(K)
20         END DO
        ENDIF
        DO 25 K=1,NM
           DJ(K)=SJ(K-1)-(K+1.0D0)*SJ(K)/X
25      END DO
        RETURN
        END

        SUBROUTINE SPHY(N,X,NM,SY,DY)
!       ======================================================
!       Purpose: Compute spherical Bessel functions yn(x) and
!                their derivatives
!       Input :  x --- Argument of yn(x) ( x ≥ 0 )
!                n --- Order of yn(x) ( n = 0,1,… )
!       Output:  SY(n) --- yn(x)
!                DY(n) --- yn'(x)
!                NM --- Highest order computed
!       ======================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT integer (I-N)
        DIMENSION SY(0:N),DY(0:N)
        NM=N
        IF (X.LT.1.0D-60) THEN
           DO 10 K=0,N
              SY(K)=-1.0D+300
              DY(K)=1.0D+300
10         END DO
           RETURN
        ENDIF
        SY(0)=-DCOS(X)/X
        F0=SY(0)
        DY(0)=(DSIN(X)+DCOS(X)/X)/X
        IF (N.LT.1) THEN
           RETURN
        ENDIF
        SY(1)=(SY(0)-DSIN(X))/X
        F1=SY(1)
        DO 15 K=2,N
           F=(2.0D0*K-1.0D0)*F1/X-F0
           SY(K)=F
           IF (DABS(F).GE.1.0D+300) GO TO 20
           F0=F1
           F1=F
15      END DO
20      NM=K-1
        DO 25 K=1,NM
           DY(K)=SY(K-1)-(K+1.0D0)*SY(K)/X
25      END DO
        RETURN
        END

        INTEGER FUNCTION MSTA1(X,MP)
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that the magnitude of
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point
!       ===================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT integer (I-N)
        A0=DABS(X)
        N0=INT(1.1D0*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20
           NN=int(N1-(N1-N0)/(1.0D0-F0/F1))
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
           F1=F
 10     END DO
 20     MSTA1=NN
        RETURN
        END

        INTEGER FUNCTION MSTA2(X,N,MP)
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT integer (I-N)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1D0*A0)+1
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=int(N1-(N1-N0)/(1.0D0-F0/F1))
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
           F1=F
10      END DO
20      MSTA2=NN+10
        RETURN
        END

        real(dp) function envj(n, x) result(r)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        r = log10(6.28_dp*n)/2 - n*log10(1.36_dp*x/n)
        end function


! Ai(z)
function airyai(z) result(res)
    implicit none
    complex(dp), intent(in) :: z
    complex(dp) :: res

    integer :: underflow, ierr
    real(dp) :: resReal, resImag
    integer, parameter :: deriv=0, scaling=1  ! compute unscaled Ai(z)

    call zairy(real(z), aimag(z), deriv, scaling, resReal, resImag, underflow, ierr)
    if(ierr > 0) then
        if(ierr == 1) then
            print *, "IERR=1, INPUT ERROR - NO COMPUTATION"
        elseif(ierr == 2) then
            print *, "IERR=2, OVERFLOW - NO COMPUTATION, REAL(ZTA) TOO LARGE ON KODE=1"
        elseif(ierr == 3) then
            print *, "IERR=3, CABS(Z) LARGE - COMPUTATION DONE BUT LOSSES OF SIGNIFCANCE"
            print *, "BY ARGUMENT REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY"
        elseif(ierr == 4) then
            print *, "IERR=4, CABS(Z) TOO LARGE - NO COMPUTATION BECAUSE OF COMPLETE"
            print *, "LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION"
        elseif(ierr == 5) then
            print *, "IERR=5, ERROR - NO COMPUTATION ALGORITHM TERMINATION CONDITION NOT MET"
        endif
        call stop_error("airyai: zairy error")
    endif
    res = resReal + i_*resImag
end function airyai


! d(Ai(z))/dz
function dairyai(z) result(res)
    implicit none
    complex(dp), intent(in) :: z
    complex(dp) :: res

    integer :: underflow, ierr
    real(dp) :: resReal, resImag
    integer, parameter :: deriv=1, scaling=1  ! compute unscaled d(Ai(z))/dz

    call zairy(real(z), aimag(z), deriv, scaling, resReal, resImag, underflow, ierr)
    if(ierr > 0) then
        if(ierr == 1) then
            print *, "IERR=1, INPUT ERROR - NO COMPUTATION"
        elseif(ierr == 2) then
            print *, "IERR=2, OVERFLOW - NO COMPUTATION, REAL(ZTA) TOO LARGE ON KODE=1"
        elseif(ierr == 3) then
            print *, "IERR=3, CABS(Z) LARGE - COMPUTATION DONE BUT LOSSES OF SIGNIFCANCE"
            print *, "BY ARGUMENT REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY"
        elseif(ierr == 4) then
            print *, "IERR=4, CABS(Z) TOO LARGE - NO COMPUTATION BECAUSE OF COMPLETE"
            print *, "LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION"
        elseif(ierr == 5) then
            print *, "IERR=5, ERROR - NO COMPUTATION ALGORITHM TERMINATION CONDITION NOT MET"
        endif
        call stop_error("dairyai: zairy error")
    endif
    res = resReal + i_*resImag
end function dairyai


! Bi(z)
function airybi(z) result(res)
    implicit none
    complex(dp), intent(in) :: z
    complex(dp) :: res

    integer :: ierr
    real(dp) :: resReal, resImag
    integer, parameter :: deriv=0, scaling=1  ! compute unscaled Bi(z)

    call zbiry(real(z), aimag(z), deriv, scaling, resReal, resImag, ierr)
    if(ierr > 0) then
        if(ierr == 1) then
            print *, "IERR=1, INPUT ERROR - NO COMPUTATION"
        elseif(ierr == 2) then
            print *, "IERR=2, OVERFLOW - NO COMPUTATION, REAL(Z) TOO LARGE ON KODE=1"
        elseif(ierr == 3) then
            print *, "IERR=3, CABS(Z) LARGE - COMPUTATION DONE BUT LOSSES OF SIGNIFCANCE"
            print *, "BY ARGUMENT REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY"
        elseif(ierr == 4) then
            print *, "IERR=4, CABS(Z) TOO LARGE - NO COMPUTATION BECAUSE OF COMPLETE"
            print *, "LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION"
        elseif(ierr == 5) then
            print *, "IERR=5, ERROR - NO COMPUTATION ALGORITHM TERMINATION CONDITION NOT MET"
        endif
        call stop_error("dairyai: zairy error")
    endif
    res = resReal + i_*resImag
end function airybi


! d(Bi(z))/dz
function dairybi(z) result(res)
    implicit none
    complex(dp), intent(in) :: z
    complex(dp) :: res

    integer :: ierr
    real(dp) :: resReal, resImag
    integer, parameter :: deriv=1, scaling=1  ! compute unscaled d(Bi(z))/dz

    call zbiry(real(z), aimag(z), deriv, scaling, resReal, resImag, ierr)
    if(ierr > 0) then
        if(ierr == 1) then
            print *, "IERR=1, INPUT ERROR - NO COMPUTATION"
        elseif(ierr == 2) then
            print *, "IERR=2, OVERFLOW - NO COMPUTATION, REAL(Z) TOO LARGE ON KODE=1"
        elseif(ierr == 3) then
            print *, "IERR=3, CABS(Z) LARGE - COMPUTATION DONE BUT LOSSES OF SIGNIFCANCE"
            print *, "BY ARGUMENT REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY"
        elseif(ierr == 4) then
            print *, "IERR=4, CABS(Z) TOO LARGE - NO COMPUTATION BECAUSE OF COMPLETE"
            print *, "LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION"
        elseif(ierr == 5) then
            print *, "IERR=5, ERROR - NO COMPUTATION ALGORITHM TERMINATION CONDITION NOT MET"
        endif
        call stop_error("dairyai: zairy error")
    endif
    res = resReal + i_*resImag
end function dairybi

! rotate Hankel functions of negative order:
function Hrotate(hkind, order) result(factor)
    integer, intent(in) :: hkind
    real(dp), intent(in) :: order
    complex(dp) :: factor

    ! NIST Handbook, 10.4.6; note the signs in the exponents
    ! (NIST gives formulas for order -v with positive v, our 'order' here is negative)

    if(order < 0.0_dp) then
        if(hkind == 1) then
            factor = exp(-pi*order*i_)
        elseif(hkind == 2) then
            factor = exp(pi*order*i_)
        endif
    else
        factor = (1.0_dp, 0.0_dp)
    endif

end function Hrotate

! 'reflect' negative integer order Bessel functions if needed
function reflect(order) result(res)
    integer, intent(in) :: order
    integer :: res

    if(order < 0) then
        res = (-1)**order
    else
        res = 1
    endif
end function

end module
