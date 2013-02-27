program test_dairyai
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: dairyai
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-14
  real(dp), parameter :: xReal(5) = [1.1_dp, 3.0_dp, 10.0_dp, 2.1_dp, 5.5_dp]
  real(dp), parameter :: xImag(5) = [0.1_dp, -0.3_dp, 2.5_dp, -1.8_dp, 2.0_dp]

  complex(dp), parameter :: correctResults(5) = [(-0.145563472763267_dp, 0.013229788969267_dp),&
       (-0.010610911598755_dp, -0.005772105021642_dp), (0.000000000007158_dp, 0.000000000580016_dp),&
       (0.067991533120730_dp, -0.048089068490876_dp), (0.000009851131022_dp, -0.000124105490600_dp)]

  do i = 1,5
     call assert(abs(dairyai(xReal(i)+i_*xImag(i)) - correctResults(i)) < eps)
  end do

end program test_dairyai
