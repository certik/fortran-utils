program test_splines
use types, only: dp
use constants, only: pi
use utils, only: assert
use mesh, only: meshexp, linspace
use splines, only: spline3, iix, iixmin, spline3ders
implicit none

real(dp), parameter :: eps = 1e-12
real(dp), allocatable :: x(:), y(:), x2(:), y2(:)
integer :: N, N2, i, ip
real(dp) :: t1, t2

! Basics:
allocate(x(5), y(5))
x = [0._dp, 1._dp, 2._dp, 3._dp, 4._dp]
y = [0._dp, 1._dp, 2._dp, 3._dp, 4._dp]
call assert(all(abs(spline3(x, y, [0._dp]) - [0._dp]) < eps))
call assert(all(abs(spline3(x, y, [0._dp, 0.5_dp, 1._dp, 10._dp]) - &
        [0._dp, 0.5_dp, 1._dp, 10._dp]) < eps))
call assert(all(abs(spline3(x, y, [10._dp, 1._dp, 0._dp, 0.5_dp, 1._dp, &
        10._dp]) - [10._dp, 1._dp, 0._dp, 0.5_dp, 1._dp, 10._dp]) < eps))

x = linspace(0.0_dp, pi, 5)
y = sin(x)
allocate(x2(100))
x2 = linspace(0.0_dp, pi, 100)
call assert(all(abs(spline3(x, y, x2) - sin(x2)) < 1.2e-2_dp))
deallocate(x, y, x2)

! Test iix correctness:
N = 100000
N2 = 1000
allocate(x(N+1), x2(N2+1))
x = meshexp(1e-7_dp, 50._dp, 1e5_dp, N)
x2 = meshexp(0._dp, 51._dp, 200._dp, N2)
call cpu_time(t1)
do i = 1, N2+1
    ip = iix(x2(i), x)
    if (x2(i) < x(1)) then
        call assert(ip == 1)
    elseif (x2(i) > x(size(x))) then
        call assert(ip == size(x)-1)
    else
        call assert(x(ip) <= x2(i) .and. x2(i) <= x(ip+1))
    end if
end do
call cpu_time(t2)
print *, "Pure iix:", t2-t1

! Test faster indexing:
call cpu_time(t1)
ip = 0
do i = 1, N2+1
    ip = iixmin(x2(i), x, ip)
    if (x2(i) < x(1)) then
        call assert(ip == 1)
    elseif (x2(i) > x(size(x))) then
        call assert(ip == size(x)-1)
    else
        call assert(x(ip) <= x2(i) .and. x2(i) <= x(ip+1))
    end if
end do
call cpu_time(t2)
print *, "Reducing size iix:", t2-t1
deallocate(x, x2)

! Correctness. This tests that BC at each end is correct:
N = 10
N2 = 1000
allocate(x(N+1), y(N+1), x2(N2+1), y2(N2+1))
x = meshexp(0._dp, pi/2, 1._dp, N)
y = sin(x)
x2 = meshexp(-pi/2, 3*pi/4, 1._dp, N2)
y2 = spline3(x, y, x2)
call assert(abs(y2(1) - sin(x2(1))) < 5e-2_dp)
call assert(abs(y2(size(y2)) - sin(x2(size(x2)))) < 5e-2_dp)
! 1st derivative
call spline3ders(x, y, x2, dynew=y2)
call assert(abs(y2(1) - cos(x2(1))) < 0.2_dp)
call assert(abs(y2(size(y2)) - cos(x2(size(x2)))) < 0.2_dp)
! 2nd derivative
call spline3ders(x, y, x2, d2ynew=y2)
call assert(abs(y2(1) + sin(x2(1))) < 0.52_dp)
call assert(abs(y2(size(y2)) + sin(x2(size(x2)))) < 0.52_dp)


x = meshexp(0._dp, 1._dp, 1._dp, N)
y = x**2
x2 = meshexp(-10._dp, 10._dp, 1._dp, N2)
y2 = spline3(x, y, x2)
call assert(abs(y2(1) - x2(1)**2) < 1e-10_dp)
call assert(abs(y2(size(y2)) - x2(size(x2))**2) < 1e-10_dp)
! 1st derivative
call spline3ders(x, y, x2, dynew=y2)
call assert(abs(y2(1) - 2*x2(1)) < 1e-10_dp)
call assert(abs(y2(size(y2)) - 2*x2(size(x2))) < 1e-10_dp)
! 2nd derivative
call spline3ders(x, y, x2, d2ynew=y2)
call assert(maxval(abs(y2 - 2)) < 1e-10_dp)
deallocate(x, y, x2, y2)

! Test speed
N  = 1000
N2 = 10000
allocate(x(N+1), y(N+1), x2(N2+1), y2(N2+1))
x = meshexp(1e-7_dp, 50._dp, 1e5_dp, N)
x2 = meshexp(0._dp, 51._dp, 200._dp, N2)
y = sin(x)
call cpu_time(t1)
y2 = spline3(x, y, x2)
call cpu_time(t2)
print *, "Spline time:", t2-t1

end program
