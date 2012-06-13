module optimize

! Optimization algorithms

use types, only: dp
use utils, only: stop_error
implicit none
private
public bisect

interface
    real(dp) function func(x)
    import :: dp
    implicit none
    real(dp), intent(in) :: x
    end function
end interface

contains

real(dp) function bisect(f, a, b, tol) result(c)
! Solves f(x) = 0 on the interval [a, b] using the bisection method
procedure(func) :: f
real(dp), intent(in) :: a, b, tol
real(dp) :: a_, b_, fa, fb, fc
a_ = a; b_ = b
fa = f(a_)
fb = f(b_)
if (fa * fb >= 0) then
    call stop_error("bisect: f(a) and f(b) must have opposite signs")
end if
do while (b_ - a_ > tol)
    c = (a_ + b_) / 2
    fc = f(c)
    if (fc == 0) return   ! We need to make sure f(c) is not zero below
    if (fa * fc < 0) then
        b_ = c
        fb = fc
    else
        a_ = c
        fa = fc
    end if
end do
c = (a_ + b_)/2
end function

end module
