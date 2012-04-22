module convolution

! Sample module to show how to calculate image convolution

implicit none
private
public apply

contains

subroutine apply(C)
integer, intent(inout) :: C(:, :, :)
integer :: C2(3, size(C, 2), size(C, 3))
integer :: k(3, 3)
integer :: i, j, o, w, h
! Filter:
!k = reshape([1, 1, 1, 1, 2, 1, 1, 1, 1], [3, 3])
k = reshape([-4, -2, 0, -2, 1, 2, 0, 2, 4], [3, 3])

w = size(C, 2)
h = size(C, 3)
C2 = C
do j = 2, h-1
    do i = 2, w-1
        do o = 1, 3
            C2(o, i, j) = sum(C(o, i-1:i+1, j-1:j+1)*k) / sum(k)
        end do
    end do
end do
C = C2
where (C < 0) C = 0
where (C > 255) C = 255
end subroutine

end module
