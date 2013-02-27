program test_hankel1
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: hankel1
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-14
  real(dp), parameter :: ordersReal(7) = [0.000000000000000_dp, 1.000000000000000_dp, 2.000000000000000_dp,&
       3.300000000000000_dp, 7.000000000000000_dp,&
       -1.500000000000000_dp, -4.000000000000000_dp]
  real(dp), parameter :: xReal(7) = [1.100000000000000_dp, 3.000000000000000_dp, 10.000000000000000_dp,&
       2.100000000000000_dp, 5.500000000000000_dp, &
       1.100000000000000_dp, 2.000000000000000_dp]
  real(dp), parameter :: xImag(7) = [0.100000000000000_dp, -0.300000000000000_dp, 2.500000000000000_dp, &
       -1.800000000000000_dp, 2.000000000000000_dp, 0.500000000000000_dp, 1.500000000000000_dp]

  complex(dp), parameter :: correctResultsRealArg(7) = [(0.719622018527511_dp, 0.162163202926887_dp), &
       (0.339058958525936_dp, 0.324674424791800_dp), (0.254630313685121_dp, -0.005868082442209_dp),&
       (0.101980666751263_dp, -1.267804061151291_dp), (0.086601225791613_dp, -0.874921069456376_dp),&
       (-0.991692967169576_dp, -0.271278757000343_dp), (0.033995719807567_dp, -2.765943226330601_dp)]

  complex(dp), parameter :: correctResultsComplexArg(7) = [(0.651368317549290_dp, 0.118990468927329_dp),&
       (0.428137777037678_dp, 0.455161854334643_dp), (0.021205380414102_dp, -0.003435444019643_dp),&
       (0.085126712134597_dp, -0.593895774097617_dp), (-0.333857447002332_dp, -0.093364319238606_dp),&
       (-0.667664135826850_dp, 0.116458951902184_dp), (-0.750994544707426_dp, 0.404920577875372_dp)]

  ! test real argument with real order:
  do i = 1,7
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
  do i = 1,7
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
