program test_cholesky
use types, only: dp
use utils, only: assert
use linalg, only: cholesky
implicit none

real(dp), parameter :: eps = 1e-12_dp
real(dp) :: A(3,3), A2(3,3), L(3,3)

A = reshape([4, 12, -16, 12, 37, -43, -16, -43, 98], shape=[3,3], order=[2,1])
L = cholesky(A)
A2 = matmul(L, transpose(L))
print *, A
print *, L
print *, A2
call assert(maxval(abs(A - A2)) < eps)

end program
