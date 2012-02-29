module linalg
use types, only: dp
use lapack, only: dsygvd
use utils, only: stop_error
implicit none
private
public eigh

contains

subroutine eigh(Am, Bm, lam, c)
! solves generalized eigen value problem for all eigenvalues and eigenvectors
! Am must by symmetric, Bm symmetric positive definite.
! Only the lower triangular part of Am and Bm is used.
real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
real(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam Bm c
real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
integer n
! lapack variables
integer lwork, liwork, info
integer, allocatable:: iwork(:)
real(dp), allocatable:: Amt(:,:), Bmt(:,:), work(:)

! solve
n = size(Am,1)
lwork = 1 + 6*n + 2*n**2
liwork = 3 + 5*n
allocate(Amt(n,n), Bmt(n,n), work(lwork), iwork(liwork))
Amt = Am; Bmt = Bm  ! Amt,Bmt temporaries overwritten by dsygvd
call dsygvd(1,'V','L',n,Amt,n,Bmt,n,lam,work,lwork,iwork,liwork,info)
if (info /= 0) then
    print *, "dsygvd returned info =", info
    if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
    else if (info <= n) then
        print *, "the algorithm failed to compute an eigenvalue while working"
        print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
        print *, "through", mod(info, n+1)
    else
        print *, "The leading minor of order ", info-n, &
            "of B is not positive definite. The factorization of B could ", &
            "not be completed and no eigenvalues or eigenvectors were computed."
    end if
    call stop_error('eigh: dsygvd error')
end if
c = Amt
end subroutine

end module
