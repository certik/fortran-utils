module constants
use types, only: dp
implicit none
private
public pi, e_, i_

! Constants contain more digits than double precision, so that
! they are rounded correctly. Single letter constants contain underscore so
! that they do not clash with user variables ("e" and "i" are frequently used as
! loop variables)
real(dp), parameter :: pi    = 3.1415926535897932384626433832795_dp
real(dp), parameter :: e_    = 2.7182818284590452353602874713527_dp
complex(dp), parameter :: i_ = (0, 1)

end module
