program test_airyai
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: airyai
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-7  ! FIXME: the results agree to only ~single precision. why?
  real(dp), parameter :: xReal(5) = [1.1, 3.0, 10.0, 2.1, 5.5]
  real(dp), parameter :: xImag(5) = [0.1, -0.3, 2.5, -1.8, 2.0]

  complex(dp), parameter :: correctResults(5) = [(0.119388546481840_dp, -0.014569895804212_dp),&
       (0.005713324781179_dp, 0.003443301475109_dp), (-0.000000000023655_dp, -0.000000000177859_dp),&
       (-0.046447326316385_dp, 0.014458111694061_dp), (0.000004409556792_dp, 0.000050512526029_dp)]

  do i = 1,5
     call assert(abs(airyai(xReal(i)+i_*xImag(i)) - correctResults(i)) < eps)
  end do

end program test_airyai
