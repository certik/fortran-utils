program test_besseli
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: besseli
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-7  ! FIXME: the results agree to only ~single precision. why?
  real(dp), parameter :: ordersReal(5) = [0.0, 1.0, 2.0, 3.3, 7.0]
  real(dp), parameter :: xReal(5) = [1.1, 3.0, 10.0, 2.1, 5.5]
  real(dp), parameter :: xImag(5) = [0.1, -0.3, 2.5, -1.8, 2.0]
  real(dp), parameter :: correctResultsRealArg(5) = [1.326160183712653_dp, 3.953370217402610_dp,&
       2281.518967726004121_dp, 0.170417364461405_dp, 0.581060275636005_dp]

  complex(dp), parameter :: correctResultsComplexArg(5) = [(1.322429481157360_dp, 0.063667990066242_dp),&
       (3.810081291811203_dp, -1.055482824433039_dp), (-1717.580096406979010_dp, 1493.230141931700473_dp),&
       (-0.334007574377871_dp, -0.129806527897784_dp), (-0.814112770843113_dp, 0.055555790888214_dp)]


  ! test real argument with real order:
  do i = 1,5
     call assert(abs(besseli(ordersReal(i), xReal(i)) - correctResultsRealArg(i)) < eps)
  end do

  ! test complex data with real order:
  do i = 1,5
     call assert(abs(besseli(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)) < eps)
  end do

end program test_besseli
