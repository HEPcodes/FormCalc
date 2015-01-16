* CTCcut.F
* the FormCalc wrapper for the CutTools three-point function
* this file is part of FormCalc
* last modified 6 Oct 12 th

#include "types.h"


	subroutine CTCmas(mi, p1, p2, p1p2, m1, m2, m3)
	implicit none
	ComplexType mi(*)
	RealType p1, p2, p1p2, m1, m2, m3
	end


************************************************************************

	ComplexType function CTCcut(mi, hel, rank, num, numeps,
     &    k1, k2, m1, m2, m3)
	implicit none
	ComplexType mi(*)
	integer hel, rank
	external num, numeps
	integer*8 k1, k2
	RealType m1, m2, m3

#include "CT.h"

	RealType kj(0:3,3)
	ComplexType mj(3)
	logical stable
	ComplexType res(0:2)
	external CTmpnum

	CTCcut = 0
	if( hel .eq. 0 ) return

	kj(0,1) = 0
	kj(1,1) = 0
	kj(2,1) = 0
	kj(3,1) = 0
	call CTMomDecode(k1, kj(0,2))
	call CTMomDecode(k2, kj(0,3))

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3

c	PRINT *, "CTCcut:", rank
c	PRINT *, "k1=", kj(0,2)
c	PRINT *, "   ", kj(1,2)
c	PRINT *, "   ", kj(2,2)
c	PRINT *, "   ", kj(3,2)
c	PRINT *, "k2=", kj(0,3)
c	PRINT *, "   ", kj(1,3)
c	PRINT *, "   ", kj(2,3)
c	PRINT *, "   ", kj(3,3)
c	PRINT *, "m1=", m1
c	PRINT *, "m2=", m2
c	PRINT *, "m3=", m3

	call ctsxcut(imode, rootsvalue, muscale, 3, num, CTmpnum,
     &    rank, kj, mj, res(0), res(1), res(2), stable)
	CTCcut = res(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setscaloop(4)
	  call ctsxcut(imode, rootsvalue, muscale, 3, numeps, CTmpnum,
     &      rank, kj, mj, res(0), res(1), res(2), stable)
	  CTCcut = CTCcut - 2*res(1)
	  call setscaloop(1)
	endif
c	PRINT *, "Ccut=", CTCcut
	end

