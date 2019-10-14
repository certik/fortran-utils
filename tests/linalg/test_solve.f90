program test_solve
use types, only: dp
use utils, only: assert
use linalg, only: solve, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(10,10), b(10), x(10)
complex(dp) :: C(3,3), d(3), y(3)

! test the eye*x = 0 with solution x = 0
A = eye(10)
b = 0.0_dp
x = solve(A,b)
call assert(maxval(x) < eps)

! test eye*x = 1, with solution x(:) = 1
b = 1.0_dp
x = solve(A,b)
call assert(maxval(abs(x-1.0_dp)) < eps)

! test i*eye*y = -1 with solution y(:) = i
d = cmplx(-1.0_dp + 0*i_)
C = i_*eye(3)
y = solve(C, d)
call assert(maxval(abs(y - i_)) < eps)

! TODO: implement complex tests
end program
