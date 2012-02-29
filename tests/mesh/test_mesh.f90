program test_mesh
use types, only: dp
use utils, only: assert
use mesh, only: meshexp, meshexp_der, meshexp_der2, get_meshexp_pars
implicit none

integer, parameter :: N = 100
real(dp), parameter :: rmin = 0._dp, rmax = 50._dp, a = 1e4_dp
real(dp), parameter :: eps = 1e-10_dp
real(dp) :: R(N+1), Rp(N+1), Rpp(N+1)
real(dp) :: rmin2, rmax2, a2
integer :: N2
R   = meshexp     (rmin, rmax, a, N)
Rp  = meshexp_der (rmin, rmax, a, N)
Rpp = meshexp_der2(rmin, rmax, a, N)
call get_meshexp_pars(R, rmin2, rmax2, a2, N2)
call assert(abs(rmin-rmin2) < eps)
call assert(abs(rmax-rmax2) < eps)
call assert(abs(a-a2) < eps)
call assert(N == N2)

R = meshexp(rmin, rmax, 1/a, N)
call get_meshexp_pars(R, rmin2, rmax2, a2, N2)
call assert(abs(rmin-rmin2) < eps)
call assert(abs(rmax-rmax2) < eps)
call assert(abs(1/a-a2) < eps)
call assert(N == N2)

R = meshexp(rmin, rmax, 1.0_dp, N)
call get_meshexp_pars(R, rmin2, rmax2, a2, N2)
call assert(abs(rmin-rmin2) < eps)
call assert(abs(rmax-rmax2) < eps)
call assert(abs(1-a2) < eps)
call assert(N == N2)

end program
