program test_hankel2
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: hankel2
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-12
  real(dp), parameter :: ordersReal(7) = [0.0_dp,  1.0_dp,  2.0_dp,  3.3_dp,  7.0_dp,  -1.5_dp,  -4.0_dp]
  real(dp), parameter :: xReal(7) = [1.1_dp,  3.0_dp,  10.0_dp,  2.1_dp,  5.5_dp,  1.1_dp,  2.0_dp]
  real(dp), parameter :: xImag(7) = [0.1_dp,  -0.3_dp,  2.5_dp,  -1.8_dp,  2.0_dp,  0.5_dp,  1.5_dp]
  complex(dp), parameter :: correctResultsRealArg(7) =[(0.719622018527511_dp, -0.162163202926887_dp),&
       (0.339058958525936_dp, -0.324674424791800_dp), (0.254630313685121_dp, 0.005868082442209_dp),&
       (0.101980666751263_dp, 1.267804061151291_dp), (0.086601225791613_dp, 0.874921069456376_dp),&
       (-0.991692967169576_dp, 0.271278757000343_dp), (0.033995719807567_dp, 2.765943226330601_dp)]

  complex(dp), parameter :: correctResultsComplexArg(7) = [(0.790792692816186_dp, -0.213286580016703_dp),&
       (0.265992551112011_dp, -0.228390788709935_dp), (2.907712660238902_dp, -0.270368699569746_dp),&
       (-0.284672444442432_dp, 0.000042806070897_dp), (0.298090132076709_dp, 0.410825383335137_dp),&
       (-0.990558382646087_dp, 0.671547450313583_dp), (0.630157438174331_dp, -0.261220122940811_dp)]

  ! test real argument with real order:
  do i = 1,7
     call assert(abs(hankel2(ordersReal(i), xReal(i)) - correctResultsRealArg(i)) < eps)
  end do

  ! test real argument with integer order:
  ! test those to machine precision
  call assert(abs(hankel2(0, xReal(1)) - hankel2FromDef(0, xReal(1))) <= epsilon(1.0_dp))
  call assert(abs(hankel2(1, xReal(2)) - hankel2FromDef(1, xReal(2))) <= epsilon(1.0_dp))
  call assert(abs(hankel2(2, xReal(3)) - hankel2FromDef(2, xReal(3))) <= epsilon(1.0_dp))
  call assert(abs(hankel2(3, xReal(4)) - hankel2FromDef(3, xReal(4))) <= epsilon(1.0_dp))
  call assert(abs(hankel2(4, xReal(5)) - hankel2FromDef(4, xReal(5))) <= epsilon(1.0_dp))

  ! test complex data with real order:
  do i = 1,7
     call assert(abs(hankel2(ordersReal(i), xReal(i)+i_*xImag(i)) - correctResultsComplexArg(i)) < eps)
  end do

contains

  function hankel2FromDef(order, x) result(z)
    integer, intent(in) :: order
    real(dp), intent(in) :: x
    complex(dp) :: z

    z = bessel_jn(order, x) - i_*bessel_yn(order, x)
  end function hankel2FromDef

end program test_hankel2
