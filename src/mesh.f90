module mesh

use types, only: dp
use utils, only: stop_error

implicit none

private
public meshexp, meshexp_der, get_meshexp_pars, meshexp_der2, &
    linspace, meshgrid

contains

function meshexp(rmin, rmax, a, N) result(mesh)
! Generates exponential mesh of N elements on [rmin, rmax]
!
! Arguments
! ---------
!
! The domain [rmin, rmax], the mesh will contain both endpoints:
real(dp), intent(in) :: rmin, rmax
!
! The ratio of the rightmost to leftmost element lengths in the mesh (for a > 1
! this means the "largest/smallest"); The only requirement is a > 0. For a == 1
! a uniform mesh will be returned:
real(dp), intent(in) :: a
!
! The number of elements in the mesh:
integer, intent(in) :: N
!
! Returns
! -------
!
! The generated mesh:
real(dp) :: mesh(N+1)
!
! Note: Every exponential mesh is fully determined by the set of parameters
! (rmin, rmax, a, N). Use the get_meshexp_pars() subroutine to obtain them
! from the given mesh.
!
! Example
! -------
!
! real(dp) :: r(11)
! r = meshexp(0._dp, 50._dp, 1e9_dp, 10)

integer :: i
real(dp) :: alpha, beta
if (a < 0) then
    call stop_error("meshexp: a > 0 required")
else if (a == 1) then
    alpha = (rmax - rmin) / N
    do i = 1, N+1
        mesh(i) = alpha * (i-1.0_dp) + rmin
    end do
else
    if (N > 1) then
        beta = log(a) / (N-1)
        alpha = (rmax - rmin) / (exp(beta*N) - 1)
        do i = 1, N+1
            mesh(i) = alpha * (exp(beta*(i-1)) - 1) + rmin
        end do
    else if (N == 1) then
        mesh(1) = rmin
        mesh(2) = rmax
    else
        call stop_error("meshexp: N >= 1 required")
    end if
end if
end function

function meshexp_der(rmin, rmax, a, N) result(Rp)
! Generates dR/dt where R(t) is the mesh returned by meshexp()
!
! Input parameters the same as for meshexp(). The variable "t" is defined by:
! t = 1, 2, ..., N+1
! So it describes a uniform mesh, with a step size 1, and the corresponding
! physical points are given by the R(t) array.
!
! Output parameters:
!     Rp(N+1) ....... dR/dt
real(dp), intent(in) :: rmin
real(dp), intent(in) :: rmax
real(dp), intent(in) :: a
integer, intent(in) :: N
real(dp) :: Rp(N+1)

integer :: i
real(dp) :: alpha, beta
if (a < 0) then
    call stop_error("meshexp_der: a > 0 required")
else if (a == 1) then
    call stop_error("meshexp_der: a == 1 not implemented")
else
    if (N > 1) then
        beta = log(a)/(N-1)
        alpha = (rmax - rmin) / (exp(beta*N) - 1)
        do i = 1, N+1
            Rp(i) = alpha * beta * exp(beta*(i-1))
        end do
    else
        call stop_error("meshexp_der: N > 1 required")
    end if
end if
end function

function meshexp_der2(rmin, rmax, a, N) result(Rpp)
! Generates d^R/dt^2 where R(t) is the mesh returned by meshexp()
!
! Input parameters the same as for meshexp(). The variable "t" is defined by:
! t = 1, 2, ..., N+1
! So it describes a uniform mesh, with a step size 1, and the corresponding
! physical points are given by the R(t) array.
!
! Output parameters:
!     Rp(N+1) ....... d^2R/dt^2
real(dp), intent(in) :: rmin
real(dp), intent(in) :: rmax
real(dp), intent(in) :: a
integer, intent(in) :: N
real(dp) :: Rpp(N+1)

integer :: i
real(dp) :: alpha, beta
if (a < 0) then
    call stop_error("meshexp_der2: a > 0 required")
else if (a == 1) then
    call stop_error("meshexp_der2: a == 1 not implemented")
else
    if (N > 1) then
        beta = log(a)/(N-1)
        alpha = (rmax - rmin) / (exp(beta*N) - 1)
        do i = 1, N+1
            Rpp(i) = alpha * beta**2 * exp(beta*(i-1))
        end do
    else
        call stop_error("meshexp_der2: N > 1 required")
    end if
end if
end function

subroutine get_meshexp_pars(R, rmin, rmax, a, N)
! Given any exponential mesh R, it determines the get_mesh()'s parameters
!
! This only looks at the number of elements, the leftmost and the rightmost
! elements (so the middle elements are not checked/taken into account).
real(dp), intent(in) :: R(:)
real(dp), intent(out) :: rmin, rmax, a
integer, intent(out) :: N
rmin = R(1)
rmax = R(size(R))
a = (R(size(R)) - R(size(R)-1)) / (R(2) - R(1))
N = size(R) - 1
end subroutine

function linspace(a, b, n) result(s)
real(dp), intent(in) :: a, b
integer, intent(in) :: n
real(dp) :: s(n)
s = meshexp(a, b, 1.0_dp, n-1)
end function

subroutine meshgrid(x, y, x2, y2)
real(dp), intent(in) :: x(:), y(:)
real(dp), intent(out) :: x2(:, :), y2(:, :)
x2 = spread(x, 1, size(y))
y2 = spread(y, 2, size(x))
end subroutine

end module
