program test_sort
use types, only: dp
use utils, only: assert
use sorting, only: sort, sortpairs
implicit none

integer :: a(5)
real(dp) :: b(5), c(5), vec(2, 5)
a = [4, 3, 2, 1, 5]
call sort(a)
call assert(all(a == [1, 2, 3, 4, 5]))
a = [5, 4, 3, 2, 1]
call sort(a)
call assert(all(a == [1, 2, 3, 4, 5]))

b = [5, 4, 3, 2, 1]
call sort(b)
call assert(all(abs(b - [1, 2, 3, 4, 5]) < epsilon(1._dp)))

! ---------------------

a = [5, 3, 4, 2, 1]
vec = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [2, 5])
call sortpairs(a, vec)
call assert(all(a == [1, 2, 3, 4, 5]))
call assert(all(abs(vec - reshape([9, 10, 7, 8, 3, 4, 5, 6, 1, 2], [2, 5])) &
    < epsilon(1._dp)))

b = [5, 3, 4, 2, 1]
vec = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [2, 5])
call sortpairs(b, vec)
call assert(all(abs(b - [1, 2, 3, 4, 5]) < epsilon(1._dp)))
call assert(all(abs(vec - reshape([9, 10, 7, 8, 3, 4, 5, 6, 1, 2], [2, 5])) &
    < epsilon(1._dp)))

b = [5, 3, 4, 2, 1]
c = [1, 2, 3, 4, 5]
call sortpairs(b, c)
call assert(all(abs(b - [1, 2, 3, 4, 5]) < epsilon(1._dp)))
call assert(all(abs(c - [5, 4, 2, 3, 1]) < epsilon(1._dp)))

end program
