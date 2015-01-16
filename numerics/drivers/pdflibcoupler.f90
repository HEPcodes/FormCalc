! pdflibcoupler.f90
!	This file supplies coupler routines for some (double precision)
!	PDFLIB functions. These are needed e.g. when doing computations
!	in quadruple precision.
!	last modified 27 May 99 th


real*16 function ALPHAS2quad(S)
real*16 :: S

real*8 :: ALPHAS2, S8
external ALPHAS2

S8 = S
ALPHAS2quad = ALPHAS2(S8)
end function ALPHAS2quad

!------------------------------------------------------------------------

subroutine ZHEEVquad(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO)
use qcomplex
implicit none
character :: JOBZ, UPLO
integer :: INFO, LDA, LWORK, N
real*16 :: RWORK(*), W(N)
type(complex32) :: A(LDA, N), WORK(*)

complex*16 :: A8(LDA, N)
real*8 :: W8(N)
integer :: i, j

do i = 1, N
  W8(i) = W(i)
  do j = 1, LDA
    A8(j, i) = A(j, i)
  enddo
enddo
call ZHEEV(JOBZ, UPLO, N, A8, LDA, W8, WORK, LWORK, RWORK, INFO)
do i = 1, N
  W(i) = W8(i)
  do j = 1, LDA
    A(j, i) = A8(j, i)
  enddo
enddo
end subroutine ZHEEVquad

!------------------------------------------------------------------------

subroutine ZGESVDquad(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, &
                      WORK, LWORK, RWORK, INFO)
use qcomplex
implicit none
character :: JOBU, JOBVT
integer :: INFO, LDA, LDU, LDVT, LWORK, M, N
real*16 :: RWORK(*), S(*)
type(complex32) :: A(LDA, *), U(LDU, *), VT(LDVT, *), WORK(*)

complex*16 :: A8(LDA, N), U8(LDU, N), VT8(LDVT, N)
real*8 :: S8(N)
integer :: i, j

do i = 1, N
  S8(i) = S(i)
  do j = 1, LDA
    A8(j, i) = A(j, i)
  enddo
  do j = 1, LDU
    U8(j, i) = U(j, i)
  enddo
  do j = 1, LDVT
    VT8(j, i) = VT(j, i)
  enddo
enddo
call ZGESVD(JOBU, JOBVT, M, N, A8, LDA, S8, U8, LDU, VT8, LDVT, &
  WORK, LWORK, RWORK, INFO)
do i = 1, N
  S(i) = S8(i)
  do j = 1, LDA
    A(j, i) = A8(j, i)
  enddo
  do j = 1, LDU
    U(j, i) = U8(j, i)
  enddo
  do j = 1, LDVT
    VT(j, i) = VT8(j, i)
  enddo
enddo
end subroutine ZGESVDquad

