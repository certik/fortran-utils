program test_bessel_zeros
use types, only: dp
use utils, only: assert
use constants, only: pi
use special, only: bessel_jn_zeros, spherical_bessel_jn, &
    spherical_bessel_jn_zeros
implicit none

real(dp), allocatable :: zeros(:, :)
integer :: n, nzeros

n = 5
nzeros = 10
allocate(zeros(nzeros, 0:n))
zeros = bessel_jn_zeros(n, nzeros, 1e-13_dp)

! The following numbers were checked against SciPy's:
!   scipy.special.jn_zeros(n, 10)
! and they all agree to all printed digits (with correct rounding)

call assert(abs(zeros( 1, 0) -  2.40482556_dp) < 1e-8_dp)
call assert(abs(zeros( 2, 0) -  5.52007811_dp) < 1e-8_dp)
call assert(abs(zeros(10, 0) - 30.63460647_dp) < 1e-8_dp)

call assert(abs(zeros( 1, 1) -  3.83170597_dp) < 1e-8_dp)
call assert(abs(zeros( 2, 1) -  7.01558667_dp) < 1e-8_dp)
call assert(abs(zeros(10, 1) - 32.18967991_dp) < 1e-8_dp)

call assert(abs(zeros( 1, 4) -  7.58834243_dp) < 1e-8_dp)
call assert(abs(zeros( 2, 4) - 11.06470949_dp) < 1e-8_dp)
call assert(abs(zeros(10, 4) - 36.69900113_dp) < 1e-8_dp)

call assert(abs(zeros( 1, 5) -  8.77148382_dp) < 1e-8_dp)
call assert(abs(zeros( 2, 5) - 12.33860420_dp) < 1e-8_dp)
call assert(abs(zeros(10, 5) - 38.15986856_dp) < 1e-8_dp)


! Test that things work for the special case n=0:

n = 0
nzeros = 10
deallocate(zeros)
allocate(zeros(nzeros, 0:n))
zeros = bessel_jn_zeros(n, nzeros, 1e-13_dp)

call assert(abs(zeros( 1, 0) -  2.40482556_dp) < 1e-8_dp)
call assert(abs(zeros( 2, 0) -  5.52007811_dp) < 1e-8_dp)
call assert(abs(zeros(10, 0) - 30.63460647_dp) < 1e-8_dp)



! Test spherical Bessel functions:
n = 5
nzeros = 10
deallocate(zeros)
allocate(zeros(nzeros, 0:n))
zeros = spherical_bessel_jn_zeros(n, nzeros, 1e-13_dp)

call assert(abs(zeros( 1, 0) -    pi) < 1e-8_dp)
call assert(abs(zeros( 2, 0) -  2*pi) < 1e-8_dp)
call assert(abs(zeros(10, 0) - 10*pi) < 1e-8_dp)

call assert(abs(zeros( 1, 1) -  4.49340946_dp) < 1e-8_dp)
call assert(abs(zeros( 2, 1) -  7.72525184_dp) < 1e-8_dp)
call assert(abs(zeros(10, 1) - 32.95638904_dp) < 1e-8_dp)

call assert(abs(zeros( 1, 4) -  8.18256145_dp) < 1e-8_dp)
call assert(abs(zeros( 2, 4) - 11.70490715_dp) < 1e-8_dp)
call assert(abs(zeros(10, 4) - 37.43173677_dp) < 1e-8_dp)

call assert(abs(zeros( 1, 5) -  9.35581211_dp) < 1e-8_dp)
call assert(abs(zeros( 2, 5) - 12.96653017_dp) < 1e-8_dp)
call assert(abs(zeros(10, 5) - 38.88363096_dp) < 1e-8_dp)

end program
