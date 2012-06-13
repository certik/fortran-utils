module special

! Special functions

use types, only: dp
use constants, only: pi
use utils, only: stop_error
use optimize, only: bisect
implicit none
private
public bessel_jn_zeros

contains

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

end module
