program test_det
use types, only: dp
use utils, only: assert
use linalg, only: det, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(3,3), Adet
complex(dp) :: B(3,3), Bdet

A = real(reshape([ 1, 2, 3, 4, 5, 6, 7, 8, -9 ], shape=[3,3]))
Adet = det(A)
call assert(abs(Adet - 54.0_dp) < eps)

B = 0*i_
B = A
B(1,1) = i_
! B has det(B) = 147 - 93i (SymPy)
Bdet = det(B)
call assert(abs(Bdet - (147.0_dp, -93.0_dp)) < eps)

end program
