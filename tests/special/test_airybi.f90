program test_airybi
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: airybi
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-14
  real(dp), parameter :: xReal(5) = [1.1_dp, 3.0_dp, 10.0_dp, 2.1_dp, 5.5_dp]
  real(dp), parameter :: xImag(5) = [0.1_dp, -0.3_dp, 2.5_dp, -1.8_dp, 2.0_dp]

  complex(dp), parameter :: correctResults(5) = [(1.299900779689054930443603552703D+00, 1.060290624984279589382296649092D-01),&
       (1.219954982972653745321167662041D+01, -6.511549752437751692468737019226D+00),&
       (-2.663057456595905125141143798828D+06, 2.762998408228511214256286621094D+08),&
       (-1.531578027758900617527615395375D+00, -1.157033930943324406825922778808D+00),&
       (-1.140718548954242521631385898218D+02, -1.292974440780106533566140569746D+03)]

  do i = 1,5
     call assert(abs(airybi(xReal(i)+i_*xImag(i)) - correctResults(i)) < eps)
  end do

end program test_airybi
