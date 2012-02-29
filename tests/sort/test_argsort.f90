program test_argsort
use types, only: dp
use utils, only: assert
use sorting, only: argsort
implicit none

call assert(all(argsort([4, 3, 2, 1, 5]) == [4, 3, 2, 1, 5]))
call assert(all(argsort([10, 9, 8, 7, 6]) == [5, 4, 3, 2, 1]))
call assert(all(argsort([1, -1]) == [2, 1]))
call assert(all(argsort([1, 2, 2, 2, 3]) == [1, 2, 3, 4, 5]))
call assert(all(argsort([2, 2, 2, 3, 1]) == [5, 2, 3, 1, 4]))

call assert(all(argsort([4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp, 5.0_dp]) == &
        [4, 3, 2, 1, 5]))
call assert(all(argsort([10.0_dp, 9.0_dp, 8.0_dp, 7.0_dp, 6.0_dp]) == &
        [5, 4, 3, 2, 1]))
call assert(all(argsort([1.0_dp, -1.0_dp]) == [2, 1]))
call assert(all(argsort([1.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 3.0_dp]) == &
        [1, 2, 3, 4, 5]))
call assert(all(argsort([2.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 1.0_dp]) == &
        [5, 2, 3, 1, 4]))

call assert(all(argsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) == &
        [4, 3, 2, 1, 5]))

end program
