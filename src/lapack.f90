module lapack_precision
implicit none
private
public dp

! This is the precision that LAPACK "d" routines were compiled with (typically
! double precision, unless a special compiler option was used while compiling
! LAPACK). This "dp" is only used in lapack.f90
! The "d" routines data type is defined as "double precision", so
! we make "dp" the same kind as 0.d0 ("double precision"), so
! as long as LAPACK and this file were compiled with the same compiler options,
! it will be consistent. (If for example all double precision is promoted to
! quadruple precision, it will be promoted both in LAPACK and here.)
integer, parameter:: dp=kind(0.d0)

end module




module lapack
implicit none

interface

    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    use lapack_precision, only: dp
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
    use lapack_precision, only: dp
    INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           AB( LDAB, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
    use lapack_precision, only: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, &
        LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, &
        IWORK, INFO )
    use lapack_precision, only: dp
    CHARACTER          FACT, UPLO
    INTEGER            INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS
    REAL(dp)           RCOND
    INTEGER            IPIV( * ), IWORK( * )
    REAL(dp)           A( LDA, * ), AF( LDAF, * ), B( LDB, * ), &
                       BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
    END SUBROUTINE

    SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, &
                       VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
                       LWORK, IWORK, IFAIL, INFO )
    use lapack_precision, only: dp
    CHARACTER          JOBZ, RANGE, UPLO
    INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
    REAL(dp)           ABSTOL, VL, VU
    INTEGER            IFAIL( * ), IWORK( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * ), &
                       Z( LDZ, * )
    END SUBROUTINE

    SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, &
        BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
    use lapack_precision, only: dp
    CHARACTER          JOBVL, JOBVR
    INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
    REAL(dp)           A( LDA, * ), ALPHAI( * ), ALPHAR( * ), &
                       B( LDB, * ), BETA( * ), VL( LDVL, * ), &
                       VR( LDVR, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
        LWORK, IWORK, LIWORK, INFO )
    use lapack_precision, only: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
    END SUBROUTINE

    REAL(dp) FUNCTION DLAMCH( CMACH )
    use lapack_precision, only: dp
    CHARACTER          CMACH
    END FUNCTION

    INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
    CHARACTER*( * )    NAME, OPTS
    INTEGER            ISPEC, N1, N2, N3, N4
    END FUNCTION

    SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
    use lapack_precision, only: dp
    INTEGER            INFO, LDA, M, N
    INTEGER            IPIV( * )
    COMPLEX(dp)        A( LDA, * )
    END SUBROUTINE

    SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
    use lapack_precision, only: dp
    INTEGER            INFO, LDA, LWORK, N
    INTEGER            IPIV( * )
    COMPLEX(dp)        A( LDA, * ), WORK( * )
    END SUBROUTINE

end interface

contains

end module
