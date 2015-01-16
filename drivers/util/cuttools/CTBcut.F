* CTBcut.F
* the FormCalc wrapper for the CutTools two-point function
* this file is part of FormCalc
* last modified 14 Feb 12 th

#include "types.h"


	ComplexType function CTBcut(rank, num, numeps, heldelta,
     &    p, m1, m2)
	implicit none
	integer rank, heldelta
	external num, numeps
	integer*8 p
	RealType m1, m2

#include "CT.h"

	RealType pj(0:3,2)
	ComplexType mj(2)
	logical stable
	ComplexType res(0:2)
	external CTmpnum

	CTBcut = 0
	if( heldelta .eq. 0 ) return

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0
	call CTMomDecode(p, pj(0,2))

	mj(1) = m1
	mj(2) = m2

c	PRINT *, "CTBcut:", rank
c	PRINT *, "p1=", pj(0,2)
c	PRINT *, "   ", pj(1,2)
c	PRINT *, "   ", pj(2,2)
c	PRINT *, "   ", pj(3,2)
c	PRINT *, "m1=", m1
c	PRINT *, "m2=", m2

	call ctsxcut(imode, rootsvalue, muscale, 2, num, CTmpnum,
     &    rank, pj, mj, res(0), res(1), res(2), stable)
	CTBcut = res(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setscaloop(4)
	  call ctsxcut(imode, rootsvalue, muscale, 2, numeps, CTmpnum,
     &      rank, pj, mj, res(0), res(1), res(2), stable)
	  CTBcut = CTBcut - 2*res(1)
	  call setscaloop(1)
	endif
c	PRINT *, "CTBcut=", CTBcut
	end

