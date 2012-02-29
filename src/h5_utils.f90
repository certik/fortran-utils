module h5_utils

use hdf5, only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, &
    H5F_ACC_TRUNC_F, h5gcreate_f, h5gclose_f, h5screate_simple_f, &
    h5sclose_f, h5dcreate_f, h5dclose_f, h5dwrite_f, h5fcreate_f, &
    h5fclose_f, h5close_f, h5open_f
use types, only: dp
use utils, only: stop_error
implicit none
private
public h5_open, h5_close, h5_write_array, h5_file, h5_create_group

type h5_file
    integer(HID_T) :: file_id
end type

interface h5_write_array
    module procedure h5_write_array_int, h5_write_array_real
end interface

integer :: number_of_open_files = 0

contains

subroutine check(error)
integer, intent(in) :: error
if (error /= 0) call stop_error("Error when calling HDF5.")
end subroutine

type(h5_file) function h5_open(filename)
character(len=*), intent(in) :: filename
integer :: error
integer(HID_T) :: file_id
if (number_of_open_files == 0) then
    ! If this is the first call to the HDF5, we need to initialize the
    ! interface:
    call h5open_f(error)
    call check(error)
end if
number_of_open_files = number_of_open_files + 1
call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
call check(error)
h5_open%file_id = file_id
end function

subroutine h5_close(self)
type(h5_file), intent(in) :: self
integer :: error
call h5fclose_f(self%file_id, error)
call check(error)
number_of_open_files = number_of_open_files - 1
if (number_of_open_files == 0) then
    ! All files are closed, so we can close the HDF5 interface:
    call h5close_f(error)
    call check(error)
end if
end subroutine

subroutine h5_write_array_int(self, a_name, A)
type(h5_file), intent(in) :: self
character(LEN=*), intent(in) :: a_name
integer, intent(in) :: A(:)

integer, parameter :: rank = 1
integer(HID_T) :: dset_id, dspace_id
integer(HSIZE_T) :: dims(1)
integer :: error

dims(1) = size(A)

call h5screate_simple_f(rank, dims, dspace_id, error)
call check(error)
call h5dcreate_f(self%file_id, a_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
call check(error)

call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, A, dims, error)
call check(error)

call h5dclose_f(dset_id, error)
call check(error)
call h5sclose_f(dspace_id, error)
call check(error)
end subroutine

subroutine h5_write_array_real(self, a_name, A)
type(h5_file), intent(in) :: self
character(LEN=*), intent(in) :: a_name
real(dp), intent(in) :: A(:)

integer, parameter :: rank = 1
integer(HID_T) :: dset_id, dspace_id
integer(HSIZE_T) :: dims(1)
integer :: error

dims(1) = size(A)

call h5screate_simple_f(rank, dims, dspace_id, error)
call check(error)
call h5dcreate_f(self%file_id, a_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
call check(error)

call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, A, dims, error)
call check(error)

call h5dclose_f(dset_id, error)
call check(error)
call h5sclose_f(dspace_id, error)
call check(error)
end subroutine

subroutine h5_create_group(self, g_name)
type(h5_file), intent(in) :: self
character(LEN=*), intent(in) :: g_name

integer(HID_T) :: group_id
integer :: error

call h5gcreate_f(self%file_id, g_name, group_id, error)
call check(error)
call h5gclose_f(group_id, error)
call check(error)
end subroutine

end module
