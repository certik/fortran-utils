program test_str
use types, only: dp
use utils, only: assert, str
implicit none

call assert(str(5) == "5")
call assert(str(5) /= " 5")
call assert(str(12) == "12")
call assert(str(12345) == "12345")
call assert(str(12.0_dp) == "12.000000")
call assert(str(12.48977_dp, 1) == "12.5")
call assert(str(12.48977_dp, 2) == "12.49")
call assert(str(12.48977_dp, 3) == "12.490")
call assert(str(12.48977_dp, 4) == "12.4898")
call assert(str(12.48977_dp, 5) == "12.48977")
call assert(str(12.48977_dp, 6) == "12.489770")
call assert(str(12.48977_dp, 7) == "12.4897700")
call assert(str(12.48977_dp, 8) == "12.48977000")
call assert(str(12.48977_dp, 12) == "12.489770000000")

end program
