program test_assert_fail1
use utils, only: assert
implicit none

! Must fail
call assert(.false.)

end program
