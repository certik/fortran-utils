program test_besseli
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: besseli
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-14
  real(dp), parameter :: ordersReal(7) = [0.000000000000000_dp, &
       & 1.000000000000000_dp, 2.000000000000000_dp, 3.300000000000000_dp, &
       & 7.000000000000000_dp, -1.500000000000000_dp, -4.000000000000000_dp]
  real(dp), parameter :: xReal(7) = [1.100000000000000_dp, &
       & 3.000000000000000_dp, 10.000000000000000_dp, 2.100000000000000_dp, &
       & 5.500000000000000_dp, 1.100000000000000_dp, 2.000000000000000_dp]
  real(dp), parameter :: xImag(7) = [0.100000000000000_dp, -0.300000000000000_dp, &
       & 2.500000000000000_dp, -1.800000000000000_dp, 2.000000000000000_dp, &
       & 0.500000000000000_dp, 1.500000000000000_dp]
  real(dp), parameter :: correctResultsRealArg(7) = [1.326160183712653_dp,&
       & 3.953370217402610_dp, 2281.518967726004121_dp, &
       & 0.170417364461405_dp, 0.581060275636005_dp, &
       & -0.137839008500075_dp, 0.050728569979180_dp]
  complex(dp), parameter :: correctResultsComplexArg(7) = &
    & [(1.322429481157360_dp, 0.063667990066242_dp), &
    & (3.810081291811203_dp, -1.055482824433039_dp), &
    & (-1717.580096406979010_dp, 1493.230141931700473_dp), &
    & (-0.334007574377871_dp, -0.129806527897784_dp), &
    & (-0.814112770843113_dp, 0.055555790888214_dp), &
    & (0.019170564094263_dp, 0.604067083160980_dp), &
    & (-0.107633719213177_dp, 0.030049276381438_dp)]

  ! test real argument with real order:
  do i = 1,7
     call assert(abs(besseli(ordersReal(i), xReal(i)) - correctResultsRealArg(i)) < eps)
  end do

  ! test complex data with real order:
  do i = 1,7
     call assert(abs(besseli(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)) < eps)
  end do

end program test_besseli
