program test_bessely
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: bessely
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-6   ! ask for only 6 digits precison
  real(dp), parameter :: ordersReal(5) = [0.0, 1.0, 2.0, 3.3, 7.0]
  real(dp), parameter :: xReal(5) = [1.1, 3.0, 10.0, 2.1, 5.5]
  real(dp), parameter :: xImag(5) = [0.1, -0.3, 2.5, -1.8, 2.0]
  real(dp), parameter :: correctResultsRealArg(5) = [0.162163202926887_dp, 0.324674424791800_dp,&
       -0.005868082442209_dp, -1.267804061151291_dp, -0.874921069456376_dp]

  complex(dp), parameter :: correctResultsComplexArg(5) = [(0.166138524472016_dp, 0.069712187633448_dp),&
       (0.341776321522289_dp, -0.081072612962834_dp), (0.133466627775052_dp, 1.443253639912400_dp),&
       (-0.296969290084257_dp, -0.184899578288515_dp), (-0.252094851286872_dp, 0.315973789539520_dp)]

  ! test real argument with real order:
  do i = 1,5
     call assert(abs(bessely(ordersReal(i), xReal(i)) - correctResultsRealArg(i)) < eps)
  end do

  ! test real argument with integer order:
  ! test those to machine precision
  call assert(abs(bessely(0, xReal(1)) - bessel_y0(xReal(1))) <= epsilon(1.0_dp))
  call assert(abs(bessely(1, xReal(2)) - bessel_y1(xReal(2))) <= epsilon(1.0_dp))
  call assert(abs(bessely(2, xReal(3)) - bessel_yn(2, xReal(3))) <= epsilon(1.0_dp))
  call assert(abs(bessely(3, xReal(4)) - bessel_yn(3, xReal(4))) <= epsilon(1.0_dp))
  call assert(abs(bessely(4, xReal(5)) - bessel_yn(4, xReal(5))) <= epsilon(1.0_dp))

  !  ! test complex data with real order:
  do i = 1,5
     print *, "Bessely", bessely(ordersReal(i), xReal(i)+i_*xImag(i))
     print *, "Scipy", correctResultsComplexArg(i)
     print *, bessely(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)
     call assert(abs(bessely(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)) < eps)
  end do

end program test_bessely
