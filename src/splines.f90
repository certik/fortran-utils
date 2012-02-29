module splines

! Splines are fully specified by the interpolation points, except that
! at the ends, we have the freedom to prescribe the second derivatives.
! If we know a derivative at an end (exactly), then best is to impose that.
! Otherwise, it is better to use the "consistent" end conditions: the second
! derivative is determined such that it is smooth.
!
! High level API: spline3, spline3ders.
! Low level API: the rest of public soubroutines.
!
! Use the high level API to obtain cubic spline fit with consistent boundary
! conditions and optionally the derivatives. Use the low level API if more fine
! grained control is needed.
!
! This module is based on a code written by John E. Pask, LLNL.

use types, only: dp
use lapack, only: dgesv, dgbsv
use utils, only: stop_error
implicit none
private
public spline3pars, spline3valder, iix, iixmin, iixun, iixexp, poly3, dpoly3, &
    d2poly3, spline3, spline3ders

contains

function spline3(x, y, xnew) result(ynew)
! Takes the function values 'y' on the grid 'x' and returns new values 'ynew'
! at the given grid 'xnew' using cubic splines interpolation with such
! boundary conditions so that the 2nd derivative is consistent with the
! interpolating cubic.
real(dp), intent(in) :: x(:), y(:), xnew(:)
real(dp) :: ynew(size(xnew))
real(dp) :: c(0:4, size(x)-1)
integer :: i, ip
! get spline parameters: 2nd derivs at ends determined by cubic interpolation
call spline3pars(x, y, [2, 2], [0._dp, 0._dp], c)

ip = 0
do i = 1, size(xnew)
    ip = iixmin(xnew(i), x, ip)
    ynew(i) = poly3(xnew(i), c(:, ip))
end do
end function

subroutine spline3ders(x, y, xnew, ynew, dynew, d2ynew)
! Just like 'spline', but also calculate 1st and 2nd derivatives
real(dp), intent(in) :: x(:), y(:), xnew(:)
real(dp), intent(out), optional :: ynew(:), dynew(:), d2ynew(:)
real(dp) :: c(0:4, size(x)-1)
integer :: i, ip
call spline3pars(x, y, [2, 2], [0._dp, 0._dp], c)

ip = 0
do i = 1, size(xnew)
    ip = iixmin(xnew(i), x, ip)
    if (present(  ynew))   ynew(i) =   poly3(xnew(i), c(:, ip))
    if (present( dynew))  dynew(i) =  dpoly3(xnew(i), c(:, ip))
    if (present(d2ynew)) d2ynew(i) = d2poly3(xnew(i), c(:, ip))
end do
end subroutine

subroutine spline3pars(xi,yi,bctype,bcval,c)
! Returns parameters c defining cubic spline interpolating x-y data xi, yi, with
! boundary conditions specified by bcytpe, bcvals
real(dp), intent(in):: xi(:)        ! x values of data
real(dp), intent(in):: yi(:)        ! y values of data
integer, intent(in):: bctype(2)     ! type of boundary condition at each end:
   ! bctype(1) = type at left end, bctype(2) = type at right end.
   ! 1 = specified 2nd derivative, 2 = 2nd derivative consistent with interpolating cubic.
real(dp), intent(in):: bcval(2)     ! boundary condition values at each end:
   ! bcval(1) = value at left end, bcval(2) = value at right end
real(dp), intent(out):: c(0:,:)     ! parameters defining spline: c(i,j) = ith parameter of jth
   ! spline polynomial, p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data pts.
   ! dimensions: c(0:4,1:n-1)
real(dp) As(5,2*size(c,2))             ! spline eq. matrix -- LAPACK band form
real(dp) bs(2*size(c,2))               ! spline eq. rhs vector
real(dp) cs(2*size(c,2))               ! spline eq. solution vector
real(dp) hi(size(c,2))                 ! spline intervals
real(dp) Ae(4,4)                       ! end-cubic eq. matrix
real(dp) be(4)                         ! end-cubic eq. rhs vector
real(dp) ce(4)                         ! end-cubic eq. solution vector
real(dp) xe(4),ye(4)                   ! x,y values at ends
real(dp) d2p1,d2pn                     ! 2nd derivatives at ends
real(dp) x0                            ! expansion center
real(dp) c1,c2,c3,c4                   ! expansion coefficients
integer n                              ! number of data points
integer i,j,i2
! lapack variables
integer ipiv(4),ipiv2(2*size(c,2))
real(dp) bemat(4,1),bmat(2*size(c,2),1)
integer info

! check input parameters
if (bctype(1) < 1 .or. bctype(1) > 2) call stop_error("spline3pars error: bctype /= 1 or 2.")
if (bctype(2) < 1 .or. bctype(2) > 2) call stop_error("spline3pars error: bctype /= 1 or 2.")
if (size(c,1) /= 5) call stop_error("spline3pars error: size(c,1) /= 5.")
if (size(c,2) /= size(xi)-1) call stop_error("spline3pars error: size(c,2) /= size(xi)-1.")
if (size(xi) /= size(yi)) call stop_error("spline3pars error: size(xi) /= size(yi)")

! To get rid of compiler warnings:
d2p1 = 0
d2pn = 0

! initializations
n=size(xi)
do i=1,n-1
   hi(i)=xi(i+1)-xi(i)
end do

! compute interpolating-cubic 2nd derivs at ends, if required
   ! left end
if(bctype(1)==2) then
   if (n < 4) call stop_error("spline3pars error: n < 4")
   xe=xi(1:4)
   ye=yi(1:4)
   x0=xe(1) ! center at end
   do i=1,4
      do j=1,4
         Ae(i,j) = (xe(i)-x0)**(j-1)
      end do
   end do
   Ae(:,1) = 1    ! set 0^0 = 1
   be=ye; bemat(:,1)=be
   call dgesv(4, 1, Ae, 4, ipiv, bemat, 4, info)
   if (info /= 0) call stop_error("spline3pars error: dgesv error.")
   ce=bemat(:,1)
   d2p1=2*ce(3)
end if
   ! right end
if(bctype(2)==2) then
   if (n < 4) call stop_error("spline3pars error: n < 4")
   xe=xi(n-3:n)
   ye=yi(n-3:n)
   x0=xe(4) ! center at end
   do i=1,4
      do j=1,4
         Ae(i,j) = (xe(i)-x0)**(j-1)
      end do
   end do
   Ae(:,1) = 1    ! set 0^0 = 1
   be=ye; bemat(:,1)=be
   call dgesv(4, 1, Ae, 4, ipiv, bemat, 4, info)
   if (info /= 0) call stop_error("spline3pars error: dgesv error.")
   ce=bemat(:,1)
   d2pn=2*ce(3)
end if

! set 2nd derivs at ends
if(bctype(1)==1) d2p1=bcval(1)
if(bctype(2)==1) d2pn=bcval(2)
!write(*,*) d2p1,d2pn

! construct spline equations -- LAPACK band form
! basis: phi1 = -(x-x_i)/h_i, phi2 = (x-x_{i+1})/h_i, phi3 = phi1^3-phi1, phi4 = phi2^3-phi2
! on interval [x_i,x_{i+1}] of length h_i = x_{i+1}-x_i
!A=0  ! full matrix
As=0
   ! left end condition
!A(1,1)=6/hi(1)**2   ! full matrix
As(4,1)=6/hi(1)**2
bs(1)=d2p1
   ! internal knot conditions
do i=2,n-1
   i2=2*(i-1)
!   A(i2,i2-1) = 1/hi(i-1)    ! full matrix ...
!   A(i2,i2)   = 2/hi(i-1)
!   A(i2,i2+1) = 2/hi(i)
!   A(i2,i2+2) = 1/hi(i)
!   A(i2+1,i2) = 1/hi(i-1)**2
!   A(i2+1,i2+1) = -1/hi(i)**2
   As(5,i2-1) = 1/hi(i-1)
   As(4,i2)   = 2/hi(i-1)
   As(3,i2+1) = 2/hi(i)
   As(2,i2+2) = 1/hi(i)
   As(5,i2)   = 1/hi(i-1)**2
   As(4,i2+1) = -1/hi(i)**2
   bs(i2) = (yi(i+1) - yi(i))/hi(i) - (yi(i) - yi(i-1))/hi(i-1)
   bs(i2+1) = 0
end do
   ! right end condition   
!A(2*(n-1),2*(n-1))=6/hi(n-1)**2 ! full matrix
As(4,2*(n-1))=6/hi(n-1)**2
bs(2*(n-1))=d2pn

! solve spline equations -- full matrix
!bmat(:,1)=bs
!call dgesv(2*(n-1), 1, A, 2*(n-1), ipiv2, bmat, 2*(n-1), info)
!if (info /= 0) call stop_error("spline3pars error: dgesv error.")
!cs=bmat(:,1)

! solve spline equations -- LAPACK band form
bmat(:,1)=bs
call dgbsv(2*(n-1), 1, 2, 1, As, 5, ipiv2, bmat, 2*(n-1), info)
if (info /= 0) call stop_error("spline3pars error: dgbsv error.")
cs=bmat(:,1)
!write(*,*) cs(1:6)
!write(*,*) cs(2*(n-1)-5:2*(n-1))

! transform to (x-x0)^(i-1) basis and return
do i=1,n-1
   ! coefficients in spline basis:
   c1=yi(i)
   c2=yi(i+1)
   c3=cs(2*i-1)
   c4=cs(2*i)
   ! coefficients in (x-x0)^(i-1) basis
   c(0,i)=xi(i)
   c(1,i)=c1
   c(2,i)=-(c1-c2+2*c3+c4)/hi(i)
   c(3,i)=3*c3/hi(i)**2
   c(4,i)=(-c3+c4)/hi(i)**3
end do
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine spline3valder(x,xi,c,val,der)
! Returns value and 1st derivative of spline defined by knots xi and parameters c
! returned by spline3pars
real(dp), intent(in):: x            ! point at which to evaluate spline
real(dp), intent(in):: xi(:)        ! spline knots (x values of data)
real(dp), intent(in):: c(0:,:)      ! spline parameters: c(i,j) = ith parameter of jth
   ! spline polynomial, p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data pts.
   ! dimensions: c(0:4,1:n-1)
real(dp), intent(out):: val         ! value of spline at x
real(dp), intent(out):: der         ! 1st derivative of spline at x
integer n                           ! number of knots
integer i1

! initialize, check input parameters
n=size(xi)
if (size(c,1) /= 5) call stop_error("spline3 error: size(c,1) /= 5.")
if (size(c,2) /= size(xi)-1) call stop_error("spline3 error: size(c,2) /= size(xi)-1.")
! find interval containing x
i1=iix(x,xi)
! return value and derivative
val=poly3(x,c(:,i1))
der=dpoly3(x,c(:,i1))
end subroutine

!--------------------------------------------------------------------------------------------------!

integer function iix(x, xi) result(i1)
! Returns index i of interval [xi(i),xi(i+1)] containing x in mesh xi,
! with intervals indexed by left-most points.
! N.B.: x outside [x1,xn] are indexed to nearest end.
! Uses bisection, except if "x" lies in the first or second elements (which is
! often the case)
real(dp), intent(in) :: x            ! target value
real(dp), intent(in) :: xi(:)        ! mesh, xi(i) < xi(i+1)
integer n                            ! number of mesh points
integer i2, ic

n = size(xi)
i1 = 1
if (n < 2) then
    call stop_error("error in iix: n < 2")
elseif (n == 2) then
    i1 = 1
elseif (n == 3) then
    if (x <= xi(2)) then ! first element
        i1 = 1
    else
        i1 = 2
    end if
elseif (x <= xi(1)) then ! left end
    i1 = 1
elseif (x <= xi(2)) then ! first element
    i1 = 1
elseif (x <= xi(3)) then ! second element
    i1 = 2
elseif (x >= xi(n)) then  ! right end
    i1 = n-1
else
    ! bisection: xi(i1) <= x < xi(i2)
    i1 = 3; i2 = n
    do
        if (i2 - i1 == 1) exit
        ic = i1 + (i2 - i1)/2
        if (x >= xi(ic)) then
            i1 = ic
        else
            i2 = ic
        endif
    end do
end if
end function

integer function iixmin(x, xi, i_min) result(ip)
! Just like iix, but assumes that x >= xi(i_min)
real(dp), intent(in) :: x, xi(:)
integer, intent(in) :: i_min
if (i_min >= 1 .and. i_min <= size(xi)-1) then
    ip = iix(x, xi(i_min:)) + i_min - 1
else
    ip = iix(x, xi)
end if
end function

!--------------------------------------------------------------------------------------------------!

function iixun(x,n,x1,xn)
! Returns index i of interval [x(i),x(i+1)] containing x in uniform mesh defined by
!   x(i) = x1 + (i-1)/(n-1)*(xn-x1), i = 1 .. n,
! with intervals indexed by left-most points.
! N.B.: x outside [x1,xn] are indexed to nearest end.
integer iixun                       ! index i of interval [x(i),x(i+1)] containing x
real(dp), intent(in):: x            ! target value
integer, intent(in):: n             ! number of mesh points
real(dp), intent(in):: x1           ! initial point of mesh
real(dp), intent(in):: xn           ! final point of mesh
integer i

! compute index
i=int((x-x1)/(xn-x1)*(n-1))+1
! reset if ouside 1..n
if (i<1) i=1
if (i>n-1) i=n-1
iixun=i
end function

!--------------------------------------------------------------------------------------------------!

function iixexp(x,n,x1,alpha,beta)
! Returns index i of interval [x(i),x(i+1)] containing x in exponential mesh defined by
!   x(i) = x1 + alpha [ exp(beta(i-1)) - 1 ], i = 1 .. n,
! where alpha = (x(n) - x(1))/[ exp(beta(n-1)) - 1 ],
! beta = log(r)/(n-2), r = (x(n)-x(n-1))/(x(2)-x(1)) = ratio of last to first interval,
! and intervals indexed by left-most points.
! N.B.: x outside [x1,xn] are indexed to nearest end.
integer iixexp                      ! index i of interval [x(i),x(i+1)] containing x
real(dp), intent(in):: x            ! target value
integer, intent(in):: n             ! number of mesh points
real(dp), intent(in):: x1           ! initial point of mesh
real(dp), intent(in):: alpha        ! mesh parameter:
!   x(i) = x1 + alpha [ exp(beta(i-1)) - 1 ], i = 1 .. n,
! where alpha = (x(n) - x(1))/[ exp(beta(n-1)) - 1 ],
! beta = log(r)/(n-2), r = (x(n)-x(n-1))/(x(2)-x(1)) = ratio of last to first interval,
real(dp), intent(in):: beta         ! mesh parameter
integer i

! compute index
i=int(log((x-x1)/alpha + 1)/beta) + 1
! reset if outside 1..n
if (i<1) i=1
if (i>n-1) i=n-1
iixexp=i
end function

!--------------------------------------------------------------------------------------------------!

function poly3(x,c)
! returns value of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) poly3
real(dp), intent(in):: x      ! point at which to evaluate polynomial
real(dp), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dx
dx=x-c(0)
poly3=c(1)+c(2)*dx+c(3)*dx**2+c(4)*dx**3
end function

!--------------------------------------------------------------------------------------------------!

function dpoly3(x,c)
! returns 1st derivative of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dpoly3
real(dp), intent(in):: x      ! point at which to evaluate polynomial
real(dp), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dx
dx=x-c(0)
dpoly3=c(2)+2*c(3)*dx+3*c(4)*dx**2
end function

!--------------------------------------------------------------------------------------------------!

function d2poly3(x,c)
! returns 2nd derivative of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) d2poly3
real(dp), intent(in):: x      ! point at which to evaluate polynomial
real(dp), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dx
dx=x-c(0)
d2poly3=2*c(3)+6*c(4)*dx
end function

end module
