module optimize

! Optimization algorithms

use types, only: dp
use utils, only: stop_error
implicit none
private
public bisect,secant

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
    if (abs(fc) < tiny(1.0_dp)) return   ! We need to make sure f(c) is not zero below
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

real(dp) function secant(f, a,b,tol, maxiter) result (c)
!Solves f(x) = 0 on using the a and b as starting values
procedure(func) :: f
real(dp),intent(in) :: a,b,tol
integer,optional :: maxiter

integer :: maxiter_
real(dp) :: xstart, xnext
real(dp) :: fstart, fnext
integer ::  i

if(present(maxiter)) then
    maxiter_ = maxiter
else
    maxiter_ = 100
endif


xstart = a
xnext = b
fstart = f(xstart)
if(abs(fstart ) < tiny(1.0_dp)) then
    c = xstart
    return
endif
fnext = f(xnext)
if(abs(fnext ) < tiny(1.0_dp)) then
    c = xnext
    return
endif

do i = 1, maxiter_
    
    if( abs(fnext - fstart) < tiny(1.0_dp)) then
       call stop_error("secant: division by zero")
    endif
    c = (xstart*fnext - xnext*fstart)/(fnext - fstart)
    if (abs(c - xnext) < tol) return ! root found 
    !update variables
    xstart = xnext
    fstart = fnext
    xnext = c
    fnext = f(c)
enddo
!max iterations number reached
call stop_error("secant: method did not converge")
end function secant
    

end module
