program test_svdvals
use types, only: dp
use utils, only: assert
use linalg, only: svdvals, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: B(3, 2), s(2), sc(2)
complex(dp) :: AC(3, 2)

! test a general complex matrix
AC = reshape([1+i_, i_, 0*i_, -i_, 2*i_, 2+0*i_], shape=[3,2], order=[2,1])
sc = svdvals(AC)
! svdvals from SciPy: 3.0269254467476365, 1.6845540477620835
call assert(maxval(abs(sc - [3.0269254467476365_dp, 1.6845540477620835_dp])) < eps)

B = reshape([1, 2, 3, 4, 5, 6], shape=[3,2], order=[2,1])
s = svdvals(B)
! svdvals from SciPy: 9.5255180915651092, 0.51430058065864448
call assert(maxval(abs(s - [9.5255180915651092_dp, 0.51430058065864448_dp])) < eps)

! note that SciPy uses the same LAPACK routines. thus, comparing with 
! SciPy's results does not test correctness of the results, instead, the
! wrapper is tested.
end program
