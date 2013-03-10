module amos
  implicit none

  ! This module contains interfaces for the amos routines (in unmodified form
  ! found on http://netlib.org/amos).
  ! The slightly modified sources taken from SciPy's "legacy/amos/special"
  ! module are linked against.

  ! The original code was written by Donald E. Amos at Sandia National Labratories
  ! See also "Portable package for Bessel functions of a complex argument and
  ! non-negative order" by D. E. Amos in Trans. Math. Software (1986)

  ! Note that we made all routines pure routines. We do this to be able to
  ! use these routines in OpenMP loops and forall statements. The F77 routines
  ! are, of course, not declared as pure. However, they behave as such - the do
  ! not have side effects. Also, the intent of every argument is clear.
  ! Technically, there is one exception: the ?1mach.f90 routines (used in
  ! finding machine constants) each contain a write statement in case the query
  ! fails. We ignore this possibility (as we trust the queries made by amos).
  ! The author of the wrapper uses this without any problems in production code.

  interface

     ! Bessel J function:
     pure subroutine zbesj(re, im, order, scaling, length, zOut_r, zOut_i, &
          underflows, ierr)
       implicit none
       double precision, intent(in) :: re, im, order
       integer, intent(in) :: scaling, length
       double precision, intent(out) :: zOut_r(length), zOut_i(length)
       integer, intent(out) :: underflows, ierr
     end subroutine zbesj

     ! Bessel Y function of second kind:
     pure subroutine zbesy(re, im, order, scaling, length, zOut_r, zOut_i, &
          underflows, workr, worki, ierr)
       implicit none
       double precision, intent(in) :: re, im, order
       integer, intent(in) :: scaling, length
       double precision, intent(out) :: workr(length), worki(length), &
            zOut_r(length), zOut_i(length)
       integer, intent(out) :: underflows, ierr
     end subroutine zbesy

     ! modified Bessel I function of first kind:
     pure subroutine zbesi(re, im, order, scaling, length, zOut_r, zOut_i, &
          underflows, ierr)
       implicit none
       double precision, intent(in) :: re, im, order
       integer, intent(in) :: scaling, length
       double precision, intent(out) :: zOut_r(length), zOut_i(length)
       integer, intent(out) :: underflows, ierr
     end subroutine zbesi

     ! modified Bessel K function of second kind:
     pure subroutine zbesk(re, im, order, scaling, length, zOut_r, zOut_i, &
          underflows, ierr)
       implicit none
       double precision, intent(in) :: re, im, order
       integer, intent(in) :: scaling, length
       double precision, intent(out) :: zOut_r(length), zOut_i(length)
       integer, intent(out) :: underflows, ierr
     end subroutine zbesk

     ! Hankel H functions:
     pure subroutine zbesh(re, im, order, scaling, hankelkind, length, zOut_r, &
          zOut_i, underflows, ierr)
       implicit none
       double precision, intent(in) :: re, im, order
       integer, intent(in) :: scaling, hankelkind, length
       double precision, intent(out) :: zOut_r(length), zOut_i(length)
       integer, intent(out) :: underflows, ierr
     end subroutine zbesh

     ! Airy function Ai (or its derivative):
     pure subroutine zairy(re, im, deriv, scaling, zOut_r, zOut_i, underflow, ierr)
       implicit none
       double precision, intent(in) :: re, im
       integer, intent(in) :: deriv, scaling
       double precision, intent(out) :: zOut_r, zOut_i
       integer, intent(out) :: underflow, ierr
     end subroutine zairy

     ! Airy function Bi (or its derivative)
     pure subroutine zbiry(re, im, deriv, scaling, zOut_r, zOut_i, ierr)
       implicit none
       double precision, intent(in) :: re, im
       integer, intent(in) :: deriv, scaling
       double precision, intent(out) :: zOut_r, zOut_i
       integer, intent(out) :: ierr
     end subroutine zbiry

  end interface

end module amos
