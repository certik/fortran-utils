program test_dairyai
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: dairybi
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-6  ! FIXME: the results agree to only ~single precision. why?
  real(dp), parameter :: xReal(5) = [1.1, 3.0, 10.0, 2.1, 5.5]
  real(dp), parameter :: xImag(5) = [0.1, -0.3, 2.5, -1.8, 2.0]

  complex(dp), parameter :: correctResults(5) = [(1.052054289401527_dp, 0.143160933990618_dp),&
       (19.320522812668564_dp, -11.875745094026080_dp), (-118511775.639862522482872_dp, 872764689.469506621360779_dp),&
       (-3.142357771618667_dp, -0.870151474552866_dp), (296.691403478894586_dp, -3076.714959917741453_dp)]

  do i = 1,5
     call assert(abs(dairybi(xReal(i)+i_*xImag(i)) - correctResults(i)) < eps)
  end do

end program test_dairyai
