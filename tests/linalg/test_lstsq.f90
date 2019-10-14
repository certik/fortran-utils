program test_lstsq
use types, only: dp
use utils, only: assert
use linalg, only: lstsq, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(10,10), b(10), x(10)
real(dp) :: AA(5,3), bb(5), xx(3)
complex(dp) :: AC(5,3), bc(5), xc(3)
complex(dp) :: C(3,3), d(3), y(3)

! solve square systems of equations, recycle the solve() tests:

! test the eye*x = 0 with solution x = 0
A = eye(10)
b = 0.0_dp
x = lstsq(A,b)
call assert(maxval(x) < eps)

! test eye*x = 1, with solution x(:) = 1
b = 1.0_dp
x = lstsq(A,b)
call assert(maxval(abs(x-1.0_dp)) < eps)

! test i*eye*y = -1 with solution y(:) = i
d = cmplx(-1.0_dp + 0*i_)
C = i_*eye(3)
y = lstsq(C, d)
call assert(maxval(abs(y - i_)) < eps)

! perform least squares computations:

! test overdetermined systems:

! use model y = ax^2 + bx + c; choose 5 points (x=0,\pm 1,\pm 2) exactly on the parabola with a=b=c=1
! expected solution: [a, b, c] = [1, 1, 1]
AA = reshape([0, 0, 1, 1, -1, 1, 1, 1, 1, 4, -2, 1, 4, 2, 1], order=[2,1], shape=[5,3])
bb = [1, 1, 3, 3, 7]  ! RHS of eq. for x=0,-1,+1,-2,+2
xx = lstsq(AA, bb)
call assert(maxval(abs(xx - [1, 1, 1])) < eps)

! use same model above, but now use complex data: a=b=1, c=i
AC = reshape([0, 0, 1, 1, -1, 1, 1, 1, 1, 4, -2, 1, 4, 2, 1], order=[2,1], shape=[5,3])
bc = [i_, i_, 2+i_, 2+i_, 6+i_]  ! RHS of eq. for x=0,-1,+1,-2,+2
xc = lstsq(AC, bc)
call assert(maxval(abs(xc - [1+0*i_, 1+0*i_, i_])) < eps)

end program
