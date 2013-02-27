program test_besselk
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: besselk
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-12   ! ask for 12 digits precison
  real(dp), parameter :: ordersReal(7) = [0.000000000000000_dp, &
       & 1.000000000000000_dp, 2.000000000000000_dp, 3.300000000000000_dp, &
       & 7.000000000000000_dp, -1.500000000000000_dp, -4.000000000000000_dp]
  real(dp), parameter :: xReal(7) = [1.100000000000000_dp, &
       & 3.000000000000000_dp, 10.000000000000000_dp, 2.100000000000000_dp, &
       & 5.500000000000000_dp, 1.100000000000000_dp, 2.000000000000000_dp]
  real(dp), parameter :: xImag(7) = [0.100000000000000_dp, -0.300000000000000_dp, &
       & 2.500000000000000_dp, -1.800000000000000_dp, 2.000000000000000_dp, &
       & 0.500000000000000_dp, 1.500000000000000_dp]
  real(dp), parameter :: correctResultsRealArg(7) = [0.365602391543186_dp, &
       & 0.040156431128194_dp, 0.000021509817007_dp, 0.745790230249143_dp, &
       & 0.096540163952931_dp, 0.759392450729773_dp, 2.195915927411959_dp]
  complex(dp), parameter :: correctResultsComplexArg(7) = &
       & [(0.361475928611182_dp, -0.050696525861221_dp), &
       & (0.037467476635317_dp, 0.014072530525557_dp), &
       & (-0.000018623842733_dp, -0.000009669896370_dp), &
       & (-0.364002441456099_dp, 0.034749972300503_dp), &
       & (-0.068914908728911_dp, 0.005086289483904_dp), &
       & (0.418185067467488_dp, -0.533740167576888_dp), &
       & (-1.020776727329708_dp, -0.110685051872821_dp)]

  ! test real argument with real order:
  do i = 1,5
     call assert(abs(besselk(ordersReal(i), xReal(i)) - correctResultsRealArg(i)) < eps)
  end do

  ! test complex data with real order:
  do i = 1,5
     call assert(abs(besselk(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)) < eps)
  end do

end program test_besselk
