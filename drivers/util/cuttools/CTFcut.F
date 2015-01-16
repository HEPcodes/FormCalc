* CTFcut.F
* the FormCalc wrapper for the CutTools six-point function
* this file is part of FormCalc
* last modified 14 Feb 12 th

#include "types.h"


	ComplexType function CTFcut(rank, num, numeps, heldelta,
     &    p1, p2, p3, p4, p5, m1, m2, m3, m4, m5, m6)
	implicit none
	integer rank, heldelta
	external num, numeps
	integer*8 p1, p2, p3, p4, p5
	RealType m1, m2, m3, m4, m5, m6

#include "CT.h"

	RealType pj(0:3,6)
	ComplexType mj(6)
	logical stable
	ComplexType res(0:2)
	external CTmpnum

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0
	call CTMomDecode(p1, pj(0,2))
	call CTMomDecode(p2, pj(0,3))
	call CTMomDecode(p3, pj(0,4))
	call CTMomDecode(p4, pj(0,5))
	call CTMomDecode(p5, pj(0,6))

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4
	mj(5) = m5
	mj(6) = m6

	call ctsxcut(imode, rootsvalue, muscale, 6, num, CTmpnum,
     &    rank, pj, mj, res(0), res(1), res(2), stable)
	CTFcut = res(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setscaloop(4)
	  call ctsxcut(imode, rootsvalue, muscale, 6, numeps, CTmpnum,
     &      rank, pj, mj, res(0), res(1), res(2), stable)
	  CTFcut = CTFcut - 2*res(0)
	  call setscaloop(1)
	endif
c	PRINT *, "CTFcut=", CTFcut
	end

