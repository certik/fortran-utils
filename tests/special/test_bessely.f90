program test_bessely
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: bessely
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
  real(dp), parameter :: correctResultsRealArg(7) = [0.162163202926887_dp, &
       & 0.324674424791800_dp, -0.005868082442209_dp, -1.267804061151291_dp, &
       & -0.874921069456376_dp, -0.271278757000343_dp, -2.765943226330601_dp]
  complex(dp), parameter :: correctResultsComplexArg(7) = &
       [(0.166138524472016_dp, 0.069712187633448_dp), &
       & (0.341776321522289_dp, -0.081072612962834_dp), &
       & (0.133466627775052_dp, 1.443253639912400_dp), &
       & (-0.296969290084257_dp, -0.184899578288515_dp), &
       & (-0.252094851286872_dp, 0.315973789539520_dp), &
       & (-0.277544249205699_dp, -0.161447123409619_dp), &
       & (0.333070350408092_dp, 0.690575991440879_dp)]

  ! test real argument with real order:
  do i = 1,7
     print *, abs(bessely(ordersReal(i), xReal(i)) - correctResultsRealArg(i))
     call assert(abs(bessely(ordersReal(i), xReal(i)) - correctResultsRealArg(i)) < eps)
  end do

  ! test real argument with integer order:
  ! test those to machine precision
  call assert(abs(bessely(0, xReal(1)) - bessel_y0(xReal(1))) <= epsilon(1.0_dp))
  call assert(abs(bessely(1, xReal(2)) - bessel_y1(xReal(2))) <= epsilon(1.0_dp))
  call assert(abs(bessely(2, xReal(3)) - bessel_yn(2, xReal(3))) <= epsilon(1.0_dp))
  call assert(abs(bessely(3, xReal(4)) - bessel_yn(3, xReal(4))) <= epsilon(1.0_dp))
  call assert(abs(bessely(4, xReal(5)) - bessel_yn(4, xReal(5))) <= epsilon(1.0_dp))

  ! test complex data with real order:
  do i = 1,7
     print *, bessely(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)
     call assert(abs(bessely(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)) < eps)
  end do

end program test_bessely
