module linalg
use types, only: dp
use lapack, only: dsygvd, ilaenv, zgetri, zgetrf
use utils, only: stop_error
implicit none
private
public eigh, inv, eye

contains

subroutine eigh(Am, Bm, lam, c)
! solves generalized eigen value problem for all eigenvalues and eigenvectors
! Am must by symmetric, Bm symmetric positive definite.
! Only the lower triangular part of Am and Bm is used.
real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
real(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam Bm c
real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
integer :: n
! lapack variables
integer :: lwork, liwork, info
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

function inv(Am) result(Bm)
! Inverts the general complex matrix Am
complex(dp), intent(in) :: Am(:,:)   ! Matrix to be inverted
complex(dp) :: Bm(size(Am, 1), size(Am, 2))   ! Bm = inv(Am)
integer :: n, nb
! lapack variables
integer :: lwork, info
complex(dp), allocatable:: Amt(:,:), work(:)
integer, allocatable:: ipiv(:)

n = size(Am, 1)
nb = ilaenv(1, 'ZGETRI', "UN", n, -1, -1, -1)
if (nb < 1) nb = max(1, n)
lwork = n*nb
allocate(Amt(n,n), ipiv(n), work(lwork))
Amt = Am
call zgetrf(n, n, Amt, n, ipiv, info)
if (info /= 0) then
    print *, "zgetrf returned info =", info
    if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
    else
        print *, "U(", info, ",", info, ") is exactly zero; The factorization"
        print *, "has been completed, but the factor U is exactly"
        print *, "singular, and division by zero will occur if it is used"
        print *, "to solve a system of equations."
    end if
    call stop_error('inv: zgetrf error')
end if
call zgetri(n, Amt, n, ipiv, work, lwork, info)
if (info /= 0) then
    print *, "zgetri returned info =", info
    if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
    else
        print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
        print *, "singular and its inverse could not be computed."
    end if
    call stop_error('inv: zgetri error')
end if
Bm = Amt
end function

function eye(n) result(A)
! Returns the identity matrix of size n x n and type real.
integer, intent(in) :: n
real(dp) :: A(n, n)
integer :: i
A = 0
do i = 1, n
    A(i, i) = 1
end do
end function

end module
