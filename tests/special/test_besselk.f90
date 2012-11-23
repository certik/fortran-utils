program test_besselk
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: besselk
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-7  ! FIXME: the results agree to only ~single precision. why?
  real(dp), parameter :: ordersReal(5) = [0.0, 1.0, 2.0, 3.3, 7.0]
  real(dp), parameter :: xReal(5) = [1.1, 3.0, 10.0, 2.1, 5.5]
  real(dp), parameter :: xImag(5) = [0.1, -0.3, 2.5, -1.8, 2.0]
  real(dp), parameter :: correctResultsRealArg(5) = [0.365602391543186_dp, 0.040156431128194_dp,&
       0.000021509817007_dp, 0.745790230249143_dp, 0.096540163952931_dp]

  complex(dp), parameter :: correctResultsComplexArg(5) = [(0.361475928611182_dp, -0.050696525861221_dp),&
       (0.037467476635317_dp, 0.014072530525557_dp), (-0.000018623842733_dp, -0.000009669896370_dp),&
       (-0.364002441456099_dp, 0.034749972300503_dp), (-0.068914908728911_dp, 0.005086289483904_dp)]


  ! test real argument with real order:
  do i = 1,5
     call assert(abs(besselk(ordersReal(i), xReal(i)) - correctResultsRealArg(i)) < eps)
  end do

  ! test complex data with real order:
  do i = 1,5
     call assert(abs(besselk(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)) < eps)
  end do

end program test_besselk
