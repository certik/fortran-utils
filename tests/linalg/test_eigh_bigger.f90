program test_eig_bigger
! You can adapt this test to test speed for bigger matrices.
use types, only: dp
use utils, only: assert
use linalg, only: eigh
implicit none

real(dp), parameter :: eps = 1e-9_dp
integer, parameter :: n = 10
real(dp) :: A(n, n), lam(n), c(n, n), r(n), norm
real(dp) :: t1, t2
integer :: co1, co2, comax, rate
integer :: i
A = 0
forall(i = 1:n) A(i, i) = 2
forall(i = 1:n-1) A(i, i+1) = -1
forall(i = 1:n-1) A(i+1, i) = -1
!To print the array (if n=10):
!print "(10(f6.2))", A
call cpu_time(t1)
call system_clock(co1, rate, comax)
call eigh(A, lam, c)
call cpu_time(t2)
call system_clock(co2, rate, comax)
print *, "Time (cpu):", t2-t1
print *, "Time (all):", (co2-co1) / real(rate, dp)
print *, "First 10 eigenvalues:"
print *, lam(:10)
do i = 1, n
    ! Test that c(:, i) are eigenvectors:
    r = matmul(A, c(:, i)) - lam(i) * c(:, i)
    call assert(sqrt(dot_product(r, r)) < eps)
    ! Test that eigenvectors are properly normalized:
    norm = dot_product(c(:, i), c(:, i))
    call assert(abs(norm - 1) < eps)
end do

end program
