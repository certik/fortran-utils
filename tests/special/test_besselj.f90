program test_besselj
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: besselj
  implicit none

  ! check a few value against the ones computed with SciPy:
  ! note we use a rather 'large' eps - SciPy computes some the values we compare
  ! to with CEPHES whereas we use AMOS

  integer :: i
  real(dp), parameter :: eps = 1d-13
  real(dp), parameter :: ordersReal(7) = [0.0_dp,  1.0_dp,  2.0_dp,  3.3_dp,  7.0_dp,  -1.5_dp,  -4.0_dp]
  real(dp),  parameter :: xReal(7) = [1.1_dp,  3.0_dp,  10.0_dp,  2.1_dp,  5.5_dp,  1.1_dp,  2.0_dp]
  real(dp),  parameter :: xImag(7) = [0.1_dp,  -0.3_dp,  2.5_dp,  -1.8_dp,  2.0_dp,  0.5_dp,  1.5_dp]
  real(dp), parameter :: correctResultsRealArg(7) = [0.71962201852751084_dp,  0.33905895852593654_dp, &
       0.25463031368512057_dp,   0.10198066675126187_dp,  0.086601225791612571_dp, -0.991692967169576_dp, 0.033995719807568_dp]
  complex(dp), parameter :: correctResultsComplexArg(7) = &
       & [(0.721080505182738_dp, -0.047148055544687_dp), &
       & (0.347065164074845_dp, 0.113385532812354_dp), &
       & (1.464459020326502_dp, -0.136902071794694_dp), &
       & (-0.099772866153918_dp, -0.296926484013360_dp), &
       & (-0.017883657462811_dp, 0.158730532048265_dp), &
       & (-0.829111259236468_dp, 0.394003201107883_dp), &
       (-0.060418553266548_dp, 0.071850227467280_dp)]


  ! test real argument with real order:
  do i = 1,7
     call assert(abs(besselj(ordersReal(i), xReal(i)) - correctResultsRealArg(i)) < eps)
  end do

  ! test real argument with integer order against the Fortran intrinsic functions:
  ! test those to machine precision
  call assert(abs(besselj(0, xReal(1)) - bessel_j0(xReal(1))) <= epsilon(1.0_dp))
  call assert(abs(besselj(1, xReal(2)) - bessel_j1(xReal(2))) <= epsilon(1.0_dp))
  do i = 2, 7
     call assert(abs(besselj(i, xReal(i)) - bessel_jn(i, xReal(i))) <= epsilon(1.0_dp))
  end do

  ! test complex data with real order:
  do i = 1,7
     call assert(abs(besselj(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)) < eps)
  end do

end program test_besselj
