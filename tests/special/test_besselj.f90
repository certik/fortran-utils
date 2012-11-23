program test_besselj
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: besselj
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-7  ! FIXME: the results agree to only ~single precision. why?
  real(dp), parameter :: ordersReal(5) = [0.0, 1.0, 2.0, 3.3, 7.0]
  real(dp), parameter :: xReal(5) = [1.1, 3.0, 10.0, 2.1, 5.5]
  real(dp), parameter :: xImag(5) = [0.1, -0.3, 2.5, -1.8, 2.0]
  real(dp), parameter :: correctResultsRealArg(5) = [0.71962201852751084, 0.33905895852593654,&
       0.25463031368512057,  0.10198066675126187, 0.086601225791612571]
  complex(dp), parameter :: correctResultsComplexArg(5) = [(0.72108050518273825,-0.047148055544687213),&
       (0.34706516407484478,+0.1133855328123541),&
       (1.4644590203265018,-0.13690207179469438),&
       (-0.099772866153917547,-0.29692648401336014),&
       (-0.01788365746281147,+0.15873053204826548)]

  ! test real argument with real order:
  do i = 1,5
     call assert(abs(besselj(ordersReal(i), xReal(i)) - correctResultsRealArg(i)) < eps)
  end do

  ! test real argument with integer order:
  ! test those to machine precision
  call assert(abs(besselj(0, xReal(1)) - bessel_j0(xReal(1))) <= epsilon(1.0_dp))
  call assert(abs(besselj(1, xReal(2)) - bessel_j1(xReal(2))) <= epsilon(1.0_dp))
  call assert(abs(besselj(2, xReal(3)) - bessel_jn(2, xReal(3))) <= epsilon(1.0_dp))
  call assert(abs(besselj(3, xReal(4)) - bessel_jn(3, xReal(4))) <= epsilon(1.0_dp))
  call assert(abs(besselj(4, xReal(5)) - bessel_jn(4, xReal(5))) <= epsilon(1.0_dp))

  ! test complex data with real order:
  do i = 1,5
     call assert(abs(besselj(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)) < eps)
  end do

end program test_besselj
