program test_hankel1
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: hankel1
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-7  ! FIXME: the results agree to only ~single precision. why?
  real(dp), parameter :: ordersReal(5) = [0.0, 1.0, 2.0, 3.3, 7.0]
  real(dp), parameter :: xReal(5) = [1.1, 3.0, 10.0, 2.1, 5.5]
  real(dp), parameter :: xImag(5) = [0.1, -0.3, 2.5, -1.8, 2.0]
  complex(dp), parameter :: correctResultsRealArg(5) = [(0.719622018527511_dp, 0.162163202926887_dp), &
       (0.339058958525936_dp, 0.324674424791800_dp), (0.254630313685121_dp, -0.005868082442209_dp),&
       (0.101980666751263_dp, -1.267804061151291_dp), (0.086601225791613_dp, -0.874921069456376_dp)]

  complex(dp), parameter :: correctResultsComplexArg(5) = [(0.651368317549290_dp, 0.118990468927329_dp),&
       (0.428137777037678_dp, 0.455161854334643_dp), (0.021205380414102_dp, -0.003435444019643_dp),&
       (0.085126712134597_dp, -0.593895774097617_dp), (-0.333857447002332_dp, -0.093364319238606_dp)]


  ! test real argument with real order:
  do i = 1,5
     call assert(abs(hankel1(ordersReal(i), xReal(i)) - correctResultsRealArg(i)) < eps)
  end do

  ! test real argument with integer order:
  ! test those to machine precision
  call assert(abs(hankel1(0, xReal(1)) - hankel1FromDef(0, xReal(1))) <= epsilon(1.0_dp))
  call assert(abs(hankel1(1, xReal(2)) - hankel1FromDef(1, xReal(2))) <= epsilon(1.0_dp))
  call assert(abs(hankel1(2, xReal(3)) - hankel1FromDef(2, xReal(3))) <= epsilon(1.0_dp))
  call assert(abs(hankel1(3, xReal(4)) - hankel1FromDef(3, xReal(4))) <= epsilon(1.0_dp))
  call assert(abs(hankel1(4, xReal(5)) - hankel1FromDef(4, xReal(5))) <= epsilon(1.0_dp))

  ! test complex data with real order:
  do i = 1,5
     call assert(abs(hankel1(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)) < eps)
  end do

contains

  function hankel1FromDef(order, x) result(z)
    integer, intent(in) :: order
    real(dp), intent(in) :: x
    complex(dp) :: z

    z = bessel_jn(order, x) + i_*bessel_yn(order, x)
  end function hankel1FromDef

end program test_hankel1
