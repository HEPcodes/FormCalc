* CTDcut.F
* the FormCalc wrapper for the CutTools four-point function
* this file is part of FormCalc
* last modified 14 Feb 12 th

#include "types.h"


	ComplexType function CTDcut(rank, num, numeps, heldelta,
     &    p1, p2, p3, m1, m2, m3, m4)
	implicit none
	integer rank, heldelta
	external num, numeps
	integer*8 p1, p2, p3
	RealType m1, m2, m3, m4

#include "CT.h"

	RealType pj(0:3,4)
	ComplexType mj(4)
	logical stable
	ComplexType res(0:2)
	external CTmpnum

	CTDcut = 0
	if( heldelta .eq. 0 ) return

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0
	call CTMomDecode(p1, pj(0,2))
	call CTMomDecode(p2, pj(0,3))
	call CTMomDecode(p3, pj(0,4))

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4

c	PRINT *, "CTDcut:", rank
c	PRINT *, "p1=", pj(0,2)
c	PRINT *, "   ", pj(1,2)
c	PRINT *, "   ", pj(2,2)
c	PRINT *, "   ", pj(3,2)
c	PRINT *, "p2=", pj(0,3)
c	PRINT *, "   ", pj(1,3)
c	PRINT *, "   ", pj(2,3)
c	PRINT *, "   ", pj(3,3)
c	PRINT *, "p3=", pj(0,4)
c	PRINT *, "   ", pj(1,4)
c	PRINT *, "   ", pj(2,4)
c	PRINT *, "   ", pj(3,4)
c	PRINT *, "m1=", m1
c	PRINT *, "m2=", m2
c	PRINT *, "m3=", m3
c	PRINT *, "m4=", m4

	call ctsxcut(imode, rootsvalue, muscale, 4, num, CTmpnum,
     &    rank, pj, mj, res(0), res(1), res(2), stable)
	CTDcut = res(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setscaloop(4)
	  call ctsxcut(imode, rootsvalue, muscale, 4, numeps, CTmpnum,
     &      rank, pj, mj, res(0), res(1), res(2), stable)
	  CTDcut = CTDcut - 2*res(1)
	  call setscaloop(1)
	endif
c	PRINT *, "CTDcut=", CTDcut
	end

