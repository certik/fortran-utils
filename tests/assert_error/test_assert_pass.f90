program test_assert_pass
use types, only: dp
use utils, only: assert
implicit none

integer :: i, a(5)
real(dp) :: b(5)

! All these must pass
call assert(.true.)
call assert(5 == 5)

i = 5
call assert(i == 5)
call assert(i /= 6)

a = [1, 2, 3, 4, 5]
call assert(all(a == [1, 2, 3, 4, 5]))
call assert(.not. all(a == [2, 2, 3, 4, 5]))
b = [1, 2, 3, 4, 5]
call assert(all(abs(b - [1, 2, 3, 4, 5]) < epsilon(1._dp)))
call assert(.not. all(abs(b - [2, 2, 3, 4, 5]) < epsilon(1._dp)))

end program
