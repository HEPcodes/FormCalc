* CTBcut.F
* the FormCalc wrapper for the CutTools two-point function
* this file is part of FormCalc
* last modified 6 Oct 12 th

#include "types.h"


	subroutine CTBmas(mi, p, m1, m2)
	implicit none
	ComplexType mi(*)
	RealType p, m1, m2
	end


************************************************************************

	ComplexType function CTBcut(mi, hel, rank, num, numeps,
     &    k, m1, m2)
	implicit none
	ComplexType mi(*)
	integer hel, rank
	external num, numeps
	integer*8 k
	RealType m1, m2

#include "CT.h"

	RealType kj(0:3,2)
	ComplexType mj(2)
	logical stable
	ComplexType res(0:2)
	external CTmpnum

	CTBcut = 0
	if( hel .eq. 0 ) return

	kj(0,1) = 0
	kj(1,1) = 0
	kj(2,1) = 0
	kj(3,1) = 0
	call CTMomDecode(k, kj(0,2))

	mj(1) = m1
	mj(2) = m2

c	PRINT *, "CTBcut:", rank
c	PRINT *, "k1=", kj(0,2)
c	PRINT *, "   ", kj(1,2)
c	PRINT *, "   ", kj(2,2)
c	PRINT *, "   ", kj(3,2)
c	PRINT *, "m1=", m1
c	PRINT *, "m2=", m2

	call ctsxcut(imode, rootsvalue, muscale, 2, num, CTmpnum,
     &    rank, kj, mj, res(0), res(1), res(2), stable)
	CTBcut = res(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setscaloop(4)
	  call ctsxcut(imode, rootsvalue, muscale, 2, numeps, CTmpnum,
     &      rank, kj, mj, res(0), res(1), res(2), stable)
	  CTBcut = CTBcut - 2*res(1)
	  call setscaloop(1)
	endif
c	PRINT *, "CTBcut=", CTBcut
	end

