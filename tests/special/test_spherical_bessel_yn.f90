program test_spherical_bessel_yn
  use types, only: dp
  use utils, only: assert
  use special, only: spherical_bessel_yn
  implicit none

  ! check a few value against the ones computed with SciPy:
  ! note we use a rather 'large' eps - SciPy computes some the values we compare
  ! to with CEPHES whereas we use AMOS

  integer :: i
  real(dp), parameter :: eps = 1e-13_dp
  integer, parameter :: ordersReal(7) = [0,  1,  2,  3,  4,  5,  6]
  real(dp),  parameter :: xReal(7) = [1.1_dp,  3.0_dp,  10.0_dp,  2.1_dp,  5.5_dp,  1.1_dp,  2.0_dp]
  real(dp), parameter :: correctResultsRealArg(7) = [-0.4123601103868884_dp, 0.06295916360231597_dp, &
       & -0.0650693049937348_dp, -1.2845701875051379_dp, -0.10423608124826822_dp, -570.9016520477184_dp, -97.79165768518729_dp]

  ! test real argument with integer order:
  do i = 1,7
     call assert(abs(spherical_bessel_yn(ordersReal(i), xReal(i)) - correctResultsRealArg(i)) < eps)
  end do

end program test_spherical_bessel_yn
