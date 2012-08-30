program test_svd
use types, only: dp
use utils, only: assert
use linalg, only: svd, eye, diag
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: B(3, 2), s(2), sc(2), U(3, 3), Vtransp(2, 2), sigma(3, 2)
complex(dp) :: AC(3, 2), Uc(3, 3), Vtranspc(2, 2)

! test if the matrix' reconstruction from the SVD is faithful

AC = reshape([1+i_, i_, 0*i_, -i_, 2*i_, 2+0*i_], shape=[3,2], order=[2,1])
call svd(AC, sc, Uc, Vtranspc)
sigma = 0.0_dp
sigma(:2,:2) = diag(sc)
call assert(maxval(abs(AC - matmul(Uc, matmul(sigma, Vtranspc)))) < eps)

B = reshape([1, 2, 3, 4, 5, 6], shape=[3,2], order=[2,1])
call svd(B, s, U, Vtransp)
sigma(:2,:2) = diag(s)
call assert(maxval(abs(B - matmul(U, matmul(sigma, Vtransp)))) < eps)

end program
