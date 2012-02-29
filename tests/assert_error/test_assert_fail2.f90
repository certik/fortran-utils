program test_assert_fail2
use utils, only: assert
implicit none

integer :: a(5)

! Must fail
a = [1, 2, 3, 4, 5]
call assert(all(a == [2, 2, 3, 4, 5]))

end program
