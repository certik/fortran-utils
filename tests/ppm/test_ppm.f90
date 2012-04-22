program test_ppm
use ppm, only: loadppm, saveppm
use utils, only: assert
use convolution, only: apply
implicit none

integer, allocatable :: img(:, :, :), img2(:, :, :)
integer :: w, h

! Test loadppm() by testing that the corners and the center of the image have
! the correct RGB colors
call loadppm("colors.ppm", img)
w = size(img, 2)
h = size(img, 3)
call assert(all(img(:, 1, 1) == [0, 0, 0]))            ! black
call assert(all(img(:, 1, h) == [255, 0, 0]))          ! red
call assert(all(img(:, w, h) == [0, 255, 0]))          ! green
call assert(all(img(:, w, 1) == [0, 0, 255]))          ! blue
call assert(all(img(:, w/2, h/2) == [255, 255, 255]))  ! white

! Test saveppm() by saving and loading the image and comparing the array:
call saveppm("tmp.ppm", img)
call loadppm("tmp.ppm", img2)
call assert(all(img == img2))

! Apply some transformation to the image:
call apply(img)
call saveppm("tmp2.ppm", img)

end program
