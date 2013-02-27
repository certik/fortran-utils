program test_dairybi
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: dairybi
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-14
  real(dp), parameter :: xReal(5) = [1.1_dp, 3.0_dp, 10.0_dp, 2.1_dp, 5.5_dp]
  real(dp), parameter :: xImag(5) = [0.1_dp, -0.3_dp, 2.5_dp, -1.8_dp, 2.0_dp]

  complex(dp), parameter :: correctResults(5) = [(1.052054289401527_dp, 0.143160933990618_dp),&
       (19.320522812668564_dp, -11.875745094026080_dp), (-118511775.639862522482872_dp, 872764689.469506621360779_dp),&
       (-3.142357771618667_dp, -0.870151474552866_dp), (296.691403478894586_dp, -3076.714959917741453_dp)]

  do i = 1,5
     call assert(abs(dairybi(xReal(i)+i_*xImag(i)) - correctResults(i)) < eps)
  end do

end program test_dairybi
