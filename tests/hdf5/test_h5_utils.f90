program test_h5_utils
! This test shows how to use hdf5 in a compact way

use types, only: dp
use h5_utils, only: h5_open, h5_create_group, h5_write_array, h5_close, h5_file
implicit none

type(h5_file) :: f
integer :: i
integer :: B(5)
real(dp) :: C(10)

do i = 1, 5
    B(i) = i
end do
do i = 1, 10
    C(i) = 1.0_dp / i
end do

f = h5_open("dsetf2.h5")
call h5_create_group(f, "/g")
call h5_create_group(f, "/g/f")
call h5_write_array(f, "/g/f/test5", B)
call h5_write_array(f, "/test6", C)
call h5_close(f)

end program
