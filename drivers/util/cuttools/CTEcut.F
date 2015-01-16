* CTEcut.F
* the FormCalc wrapper for the CutTools five-point function
* this file is part of FormCalc
* last modified 14 Feb 12 th

#include "types.h"


	ComplexType function CTEcut(rank, num, numeps, heldelta,
     &    p1, p2, p3, p4, m1, m2, m3, m4, m5)
	implicit none
	integer rank, heldelta
	external num, numeps
	integer*8 p1, p2, p3, p4
	RealType m1, m2, m3, m4, m5

#include "CT.h"

	RealType pj(0:3,5)
	ComplexType mj(5)
	logical stable
	ComplexType res(0:2)
	external CTmpnum

	CTEcut = 0
	if( heldelta .eq. 0 ) return

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0
	call CTMomDecode(p1, pj(0,2))
	call CTMomDecode(p2, pj(0,3))
	call CTMomDecode(p3, pj(0,4))
	call CTMomDecode(p4, pj(0,5))

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4
	mj(5) = m5

	call ctsxcut(imode, rootsvalue, muscale, 5, num, CTmpnum,
     &    rank, pj, mj, res(0), res(1), res(2), stable)
	CTEcut = res(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setscaloop(4)
	  call ctsxcut(imode, rootsvalue, muscale, 5, numeps, CTmpnum,
     &      rank, pj, mj, res(0), res(1), res(2), stable)
	  CTEcut = CTEcut - 2*res(1)
	  call setscaloop(1)
	endif
	PRINT *, "CTEcut=", CTEcut
	end

