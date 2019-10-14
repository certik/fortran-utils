program test_eigvals
use types, only: dp
use utils, only: assert
use linalg, only: eigvals, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: B(5, 5)
complex(dp) :: AC(2, 2), lam(2), lamb(5)

! test a general matrix
AC = reshape([1.0_dp+0*i_, 2*i_, 3*i_, -4.0_dp+0*i_], shape=[2,2])
lam = eigvals(AC)
call assert(all(abs(lam - [-1.0_dp, -2.0_dp]) < eps))

! test a multiple of the unit matrix:
B = 3*eye(5)
lamb = eigvals(B)
call assert(maxval(abs(lamb - 3.0_dp)) < eps)  ! all eigenvalues are 3

! TODO: add more tests

end program
