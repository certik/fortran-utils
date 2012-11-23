program test_hankel2
  use types, only: dp
  use constants, only: i_
  use utils, only: assert
  use special, only: hankel2
  implicit none

  ! check a few value against the ones computed with SciPy:
  integer :: i
  real(dp), parameter :: eps = 1d-7  ! FIXME: the results agree to only ~single precision. why?
  real(dp), parameter :: ordersReal(5) = [0.0, 1.0, 2.0, 3.3, 7.0]
  real(dp), parameter :: xReal(5) = [1.1, 3.0, 10.0, 2.1, 5.5]
  real(dp), parameter :: xImag(5) = [0.1, -0.3, 2.5, -1.8, 2.0]
  complex(dp), parameter :: correctResultsRealArg(5) =[(0.719622018527511_dp, -0.162163202926887_dp),&
       (0.339058958525936_dp, -0.324674424791800_dp), (0.254630313685121_dp, 0.005868082442209_dp),&
       (0.101980666751263_dp, 1.267804061151291_dp), (0.086601225791613_dp, 0.874921069456376_dp)]


  complex(dp), parameter :: correctResultsComplexArg(5) = [(0.790792692816186_dp, -0.213286580016703_dp),&
       (0.265992551112011_dp, -0.228390788709935_dp), (2.907712660238902_dp, -0.270368699569746_dp),&
       (-0.284672444442432_dp, 0.000042806070897_dp), (0.298090132076709_dp, 0.410825383335137_dp)]

  ! test real argument with real order:
  do i = 1,5
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
  do i = 1,5
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
