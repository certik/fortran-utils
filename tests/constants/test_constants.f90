program test_constants
use types, only: dp
use constants, only: pi, e_, i_
use utils, only: assert
implicit none

! Euler's identity:
call assert(abs(e_**(pi*i_) + 1) < 1e-15_dp)

end program
