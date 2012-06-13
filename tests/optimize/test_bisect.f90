program test_bisect
use types, only: dp
use utils, only: assert
use optimize, only: bisect
implicit none

real(dp) :: x

x = bisect(f1, -4._dp, 4._dp, 1e-12_dp)
print *, x
call assert(abs(x - 0.73908513_dp) < 1e-8_dp)

x = bisect(f2, -4._dp, 4._dp, 1e-12_dp)
print *, x
call assert(abs(x - 0.86547403_dp) < 1e-8_dp)

x = bisect(f3, 0._dp, 3._dp, 1e-12_dp)
print *, x
call assert(abs(x - 2.40482556_dp) < 1e-8_dp)

x = bisect(f3, 3._dp, 6._dp, 1e-12_dp)
print *, x
call assert(abs(x - 5.52007811_dp) < 1e-8_dp)

x = bisect(f4, -4._dp, 4._dp, 1e-12_dp)
print *, x
call assert(abs(x) < 1e-8_dp)

contains

real(dp) function f1(x)
real(dp), intent(in) :: x
f1 = cos(x) - x
end function

real(dp) function f2(x)
real(dp), intent(in) :: x
f2 = cos(x) - x**3
end function

real(dp) function f3(x)
real(dp), intent(in) :: x
f3 = bessel_jn(0, x)
end function

real(dp) function f4(x)
real(dp), intent(in) :: x
f4 = x
end function

end program
