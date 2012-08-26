module linalg
  use types, only: dp
  use lapack, only: dsyevd, dsygvd, ilaenv, zgetri, zgetrf, zheevd, &
       dgeev, zgeev, zhegvd, dgesv, zgesv, dgetrf, dgetri
  use utils, only: stop_error
  use constants, only: i_
  implicit none
  private
  public eig, eigvals, eigh, inv, solve, eye, det

  ! eigenvalue/-vector problem for general matrices:
  interface eig
     module procedure deig
     module procedure zeig
  end interface eig

  ! eigenvalue/-vector problem for real symmetric/complex hermitian matrices:
  interface eigh
     module procedure deigh_generalized
     module procedure deigh_simple
     module procedure zeigh_generalized
     module procedure zeigh_simple
  end interface eigh

  ! eigenvalues for general matrices:
  interface eigvals
     module procedure deigvals
     module procedure zeigvals
  end interface eigvals

  ! matrix inversion for real/complex matrices:
  interface inv
     module procedure dinv
     module procedure zinv
  end interface inv

  ! solution to linear systems of equation with real/complex coefficients:
  interface solve
     module procedure dsolve
     module procedure zsolve
  end interface solve

  interface det
     module procedure ddet
     module procedure zdet
  end interface det

contains

  ! TODO: add optional switch for left or right eigenvectors in deig() and zeig()?
  subroutine deig(A, lam, c)
    real(dp), intent(in) :: A(:, :)  ! matrix for eigenvalue compuation
    complex(dp), intent(out) :: lam(:)  ! eigenvalues: A c = lam c
    complex(dp), intent(out) :: c(:, :)  ! eigenvectors: A c = lam c; c(i,j) = ith component of jth vec.
    ! LAPACK variables for DGEEV:
    real(dp), allocatable ::  At(:,:), vl(:,: ), vr(:,:), wi(:), work(:), wr(:)
    integer :: info, lda, ldvl, ldvr, lwork, n, i

    lda = size(A(:,1))
    n = size(A(1,:))
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    allocate(At(lda,n), wr(n), wi(n), vl(ldvl,n), vr(ldvr,n), work(lwork))
    At = A

    call dgeev('N', 'V', n, At, lda, wr, wi, vl, ldvl, vr, ldvr, &
         work, lwork, info)
    if(info /= 0) then
       print *, "dgeev returned info = ", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "the QR algorithm failed to compute all the"
          print *, "eigenvalues, and no eigenvectors have been computed;"
          print *, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print *, "have converged."
       end if
       call stop_error('eig: dgeev error')
    end if

    lam = wr + i_*wi
    ! as DGEEV has a rather complicated way of returning the eigenvectors,
    ! it is necessary to build the complex array of eigenvectors from
    ! two real arrays:
    do i = 1,n
       if(wi(i) > 0.0) then  ! first of two conjugate eigenvalues
          c(:, i) = vr(:, i) + i_*vr(:, i+1)
       elseif(wi(i) < 0.0_dp) then  ! second of two conjugate eigenvalues
          c(:, i) = vr(:, i-1) - i_*vr(:, i)
       else
          c(:, i) = vr(:, i)
       end if
    end do
  end subroutine deig

  subroutine zeig(A, lam, c)
    complex(dp), intent(in) :: A(:, :)  ! matrix to solve eigenproblem for
    complex(dp), intent(out) :: lam(:)  ! eigenvalues: A c = lam c
    complex(dp), intent(out) :: c(:,:)  ! eigenvectors: A c = lam c; c(i,j) = ith component of jth vec.
    ! LAPACK variables:
    integer :: info, lda, ldvl, ldvr, lwork, n, lrwork
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: vl(:,:), vr(:,:), work(:)

    lda = size(A(:,1))
    n = size(A(1,:))
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    lrwork = 2*n
    allocate(vl(ldvl,n), vr(ldvr,n), work(lwork), rwork(lrwork))
    c = A
    call zgeev('N', 'V', n, c, lda, lam, vl, ldvl, vr, ldvr, work, &
         lwork, rwork, info)
    if(info /= 0) then
       print *, "zgeev returned info = ", info
       if(info < 0) then
          print *, "the ",-info, "-th argument had an illegal value."
       else
          print *, "the QR algorithm failed to compute all the"
          print *, "eigenvalues, and no eigenvectors have been computed;"
          print *, "elements and ", info+1, ":", n, " of W contain eigenvalues which have"
          print *, "converged."
       end if
       call stop_error('eig: zgeev error')
    end if
    c = vr
  end subroutine zeig

  function deigvals(A) result(lam)
    real(dp), intent(in) :: A(:, :)  ! matrix for eigenvalue compuation
    complex(dp), allocatable :: lam(:)  ! eigenvalues: A c = lam c
    ! LAPACK variables for DGEEV:
    real(dp), allocatable ::  At(:,:), vl(:,: ), vr(:,:), wi(:), work(:), wr(:)
    integer :: info, lda, ldvl, ldvr, lwork, n

    lda = size(A(:,1))
    n = size(A(1,:))
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    allocate(At(lda,n), wr(n), wi(n), vl(ldvl,n), vr(ldvr,n), work(lwork), lam(n))
    At = A

    call dgeev('N', 'N', n, At, lda, wr, wi, vl, ldvl, vr, ldvr, &
         work, lwork, info)
    if(info /= 0) then
       print *, "dgeev returned info = ", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "the QR algorithm failed to compute all the"
          print *, "eigenvalues, and no eigenvectors have been computed;"
          print *, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print *, "have converged."
       end if
       call stop_error('eigvals: dgeev error')
    end if

    lam = wr + i_*wi
  end function deigvals

  function zeigvals(A) result(lam)
    complex(dp), intent(in) :: A(:, :)  ! matrix to solve eigenproblem for
    complex(dp), allocatable :: lam(:)  ! eigenvalues: A c = lam c
    ! LAPACK variables:
    integer :: info, lda, ldvl, ldvr, lwork, n, lrwork
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: At(:,:), vl(:,:), vr(:,:), work(:)

    lda = size(A(:,1))
    n = size(A(1,:))
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    lrwork = 2*n
    allocate(At(lda,n), vl(ldvl,n), vr(ldvr,n), work(lwork), rwork(lrwork), lam(n))
    At = A
    call zgeev('N', 'N', n, At, lda, lam, vl, ldvl, vr, ldvr, work, &
         lwork, rwork, info)
    if(info /= 0) then
       print *, "zgeev returned info = ", info
       if(info < 0) then
          print *, "the ",-info, "-th argument had an illegal value."
       else
          print *, "the QR algorithm failed to compute all the"
          print *, "eigenvalues, and no eigenvectors have been computed;"
          print *, "elements and ", info+1, ":", n, " of W contain eigenvalues which have"
          print *, "converged."
       end if
       call stop_error('eig: zgeev error')
    end if
  end function zeigvals

  subroutine deigh_generalized(Am, Bm, lam, c)
    ! solves generalized eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric, Bm symmetric positive definite.
    ! Only the lower triangular part of Am and Bm is used.
    real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
    real(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam Bm c
    real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
    integer :: n
    ! lapack variables
    integer :: lwork, liwork, info
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: Bmt(:,:), work(:)

    ! solve
    n = size(Am,1)
    lwork = 1 + 6*n + 2*n**2
    liwork = 3 + 5*n
    allocate(Bmt(n,n), work(lwork), iwork(liwork))
    c = Am; Bmt = Bm  ! Bmt temporaries overwritten by dsygvd
    call dsygvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "dsygvd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else if (info <= n) then
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          print *, "through", mod(info, n+1)
       else
          print *, "The leading minor of order ", info-n, &
               "of B is not positive definite. The factorization of B could ", &
               "not be completed and no eigenvalues or eigenvectors were computed."
       end if
       call stop_error('eigh: dsygvd error')
    end if
  end subroutine deigh_generalized

  subroutine deigh_simple(Am, lam, c)
    ! solves eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric
    ! Only the lower triangular part of Am is used.
    real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam c
    real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam c; c(i,j) = ith component of jth vec.
    integer :: n
    ! lapack variables
    integer :: lwork, liwork, info
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: work(:)

    ! solve
    n = size(Am,1)
    lwork = 1 + 6*n + 2*n**2
    liwork = 3 + 5*n
    allocate(work(lwork), iwork(liwork))
    c = Am
    call dsyevd('V','L',n,c,n,lam,work,lwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "dsyevd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          print *, "through", mod(info, n+1)
       end if
       call stop_error('eigh: dsyevd error')
    end if
  end subroutine deigh_simple

  subroutine zeigh_generalized(Am, Bm, lam, c)
    ! solves generalized eigen value problem for all eigenvalues and eigenvectors
    ! Am must by hermitian, Bm hermitian positive definite.
    ! Only the lower triangular part of Am and Bm is used.
    complex(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
    complex(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
    real(dp), intent(out) :: lam(:)      ! eigenvalues: Am c = lam Bm c
    complex(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
    ! lapack variables
    integer :: info, liwork, lrwork, lwork, n
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: Bmt(:,:), work(:)

    ! solve
    n = size(Am,1)
    lwork = 2*n + n**2
    lrwork = 1 + 5*N + 2*n**2
    liwork = 3 + 5*n
    allocate(Bmt(n,n), work(lwork), rwork(lrwork), iwork(liwork))
    c = Am; Bmt = Bm  ! Bmt temporary overwritten by zhegvd
    call zhegvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,rwork,lrwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "zhegvd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else if (info <= n) then
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          print *, "through", mod(info, n+1)
       else
          print *, "The leading minor of order ", info-n, &
               "of B is not positive definite. The factorization of B could ", &
               "not be completed and no eigenvalues or eigenvectors were computed."
       end if
       call stop_error('eigh: zhegvd error')
    end if
  end subroutine zeigh_generalized

  subroutine zeigh_simple(Am, lam, c)
    ! solves eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric
    ! Only the lower triangular part of Am is used.
    complex(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam c
    complex(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam c; c(i,j) = ith component of jth vec.
    ! LAPACK variables:
    integer :: info, lda, liwork, lrwork, lwork, n
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)

    ! use LAPACK's zheevd routine
    n = size(Am, 1)
    lda = max(1, n)
    lwork = 2*n + n**2
    lrwork = 1 + 5*n + 2*n**2
    liwork = 3 + 5*n
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    c = Am
    call zheevd("V", "L", n, c, lda, lam, work, lwork, rwork, lrwork, &
         iwork, liwork, info)
    if (info /= 0) then
       print *, "zheevd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "through the submatrix lying in rows and columns through"
          print *, info/(n+1), " through ", mod(info, n+1)
       end if
       call stop_error('eigh: zheevd error')
    end if
  end subroutine zeigh_simple

  function dinv(Am) result(Bm)
    real(dp), intent(in) :: Am(:,:)  ! matrix to be inverted
    real(dp) :: Bm(size(Am, 1), size(Am, 2))   ! Bm = inv(Am)
    real(dp), allocatable :: Amt(:,:), work(:)  ! temporary work arrays

    ! LAPACK variables:
    integer ::  info, lda, n, lwork, nb
    integer, allocatable :: ipiv(:)

    ! use LAPACK's dgetrf and dgetri
    n = size(Am(1, :))
    lda = n
    nb = ilaenv(1, 'DGETRI', "UN", n, -1, -1, -1)
    lwork = n*nb
    if (nb < 1) nb = max(1, n)
    allocate(Amt(n,n), work(lwork), ipiv(n))
    Amt = Am
    call dgetrf(n, n, Amt, lda, ipiv, info)
    if(info /= 0) then
       print *, "dgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       call stop_error('inv: dgetrf error')
    end if

    call dgetri(n, Amt, n, ipiv, work, lwork, info)
    if (info /= 0) then
       print *, "dgetri returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
          print *, "singular and its inverse could not be computed."
       end if
       call stop_error('inv: dgetri error')
    end if
    Bm = Amt

  end function dinv

  function zinv(Am) result(Bm)
    ! Inverts the general complex matrix Am
    complex(dp), intent(in) :: Am(:,:)   ! Matrix to be inverted
    complex(dp) :: Bm(size(Am, 1), size(Am, 2))   ! Bm = inv(Am)
    integer :: n, nb
    ! lapack variables
    integer :: lwork, info
    complex(dp), allocatable:: Amt(:,:), work(:)
    integer, allocatable:: ipiv(:)

    n = size(Am, 1)
    nb = ilaenv(1, 'ZGETRI', "UN", n, -1, -1, -1)
    if (nb < 1) nb = max(1, n)
    lwork = n*nb
    allocate(Amt(n,n), ipiv(n), work(lwork))
    Amt = Am
    call zgetrf(n, n, Amt, n, ipiv, info)
    if (info /= 0) then
       print *, "zgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       call stop_error('inv: zgetrf error')
    end if
    call zgetri(n, Amt, n, ipiv, work, lwork, info)
    if (info /= 0) then
       print *, "zgetri returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
          print *, "singular and its inverse could not be computed."
       end if
       call stop_error('inv: zgetri error')
    end if
    Bm = Amt
  end function zinv

  function dsolve(A, b) result(x)
    ! solves a system of equations A x = b with one right hand side
    real(dp), intent(in) :: A(:,:)  ! coefficient matrix A
    real(dp), intent(in) :: b(:)  ! right-hand-side A x = b
    real(dp), allocatable :: x(:)
    ! LAPACK variables:
    real(dp), allocatable :: At(:,:), bt(:,:)
    integer :: n, info, lda
    integer, allocatable :: ipiv(:)

    n = size(A(1,:))
    lda = size(A(:, 1))
    allocate(At(lda,n), bt(n,1), ipiv(n), x(n))
    At = A
    bt(:,1) = b(:)
    call dgesv(n, 1, At, lda, ipiv, bt, n, info)
    if(info /= 0) then
       print *, "dgesv returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, so the solution could not be computed."
       end if
       call stop_error('inv: dgesv error')
    endif
    x = bt(:,1)
  end function dsolve

  function zsolve(A, b) result(x)
    ! solves a system of equations A x = b with one right hand side
    complex(dp), intent(in) :: A(:,:)  ! coefficient matrix A
    complex(dp), intent(in) :: b(:)  ! right-hand-side A x = b
    complex(dp), allocatable :: x(:)
    ! LAPACK variables:
    complex(dp), allocatable :: At(:,:), bt(:,:)
    integer :: n, info, lda
    integer, allocatable :: ipiv(:)

    n = size(A(1,:))
    lda = size(A(:, 1))
    allocate(At(lda,n), bt(n,1), ipiv(n), x(n))
    At = A
    bt(:,1) = b(:)
    call zgesv(n, 1, At, lda, ipiv, bt, n, info)
    if(info /= 0) then
       print *, "zgesv returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, so the solution could not be computed."
       end if
       call stop_error('inv: zgesv error')
    endif
    x = bt(:,1)
  end function zsolve

  function eye(n) result(A)
    ! Returns the identity matrix of size n x n and type real.
    integer, intent(in) :: n
    real(dp) :: A(n, n)
    integer :: i

    A = 0
    do i = 1, n
       A(i, i) = 1
    end do
  end function eye

  function ddet(A) result(x)
    ! compute the determinant of a real matrix using an LU factorization
    real(dp), intent(in) :: A(:, :)
    real(dp) :: x
    integer :: i
    ! LAPACK variables:
    integer :: info, n
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: At(:,:)

    n = size(A(1,:))
    allocate(At(n,n), ipiv(n))
    At = A
    call dgetrf(n, n, At, n, ipiv, info)
    if(info /= 0) then
       print *, "dgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       call stop_error('det: dgetrf error')
    end if

    ! At now contains the LU of the factorization A = PLU
    ! as L has unit diagonal entries, the determinant can be computed
    ! from the product of U's diagonal entries. Additional sign changes
    ! stemming from the permutations P have to be taken into account as well.
    x = 1.0_dp
    do i = 1,n
       if(ipiv(i) /= i) then  ! additional sign change
          x = -x*At(i,i)
       else
          x = x*At(i,i)
       endif
    end do
  end function ddet

  function zdet(A) result(x)
    ! compute the determinant of a real matrix using an LU factorization
    complex(dp), intent(in) :: A(:, :)
    complex(dp) :: x
    integer :: i
    ! LAPACK variables:
    integer :: info, n
    integer, allocatable :: ipiv(:)
    complex(dp), allocatable :: At(:,:)

    n = size(A(1,:))
    allocate(At(n,n), ipiv(n))
    At = A
    call zgetrf(n, n, At, n, ipiv, info)
    if(info /= 0) then
       print *, "zgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       call stop_error('inv: zgetrf error')
    end if

    ! for details on the computation, compare the comment in ddet().
    x = 1.0_dp + 0*i_
    do i = 1,n
       if(ipiv(i) /= i) then  ! additional sign change
          x = -x*At(i,i)
       else
          x = x*At(i,i)
       endif
    end do
  end function zdet

end module linalg
