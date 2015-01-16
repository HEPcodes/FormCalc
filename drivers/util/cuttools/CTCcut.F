* CTCcut.F
* the FormCalc wrapper for the CutTools three-point function
* this file is part of FormCalc
* last modified 14 Feb 12 th

#include "types.h"


	ComplexType function CTCcut(rank, num, numeps, heldelta,
     &    p1, p2, m1, m2, m3)
	implicit none
	integer rank, heldelta
	external num, numeps
	integer*8 p1, p2
	RealType m1, m2, m3

#include "CT.h"

	RealType pj(0:3,3)
	ComplexType mj(3)
	logical stable
	ComplexType res(0:2)
	external CTmpnum

	CTCcut = 0
	if( heldelta .eq. 0 ) return

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0
	call CTMomDecode(p1, pj(0,2))
	call CTMomDecode(p2, pj(0,3))

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3

c	PRINT *, "CTCcut:", rank
c	PRINT *, "p1=", pj(0,2)
c	PRINT *, "   ", pj(1,2)
c	PRINT *, "   ", pj(2,2)
c	PRINT *, "   ", pj(3,2)
c	PRINT *, "p2=", pj(0,3)
c	PRINT *, "   ", pj(1,3)
c	PRINT *, "   ", pj(2,3)
c	PRINT *, "   ", pj(3,3)
c	PRINT *, "m1=", m1
c	PRINT *, "m2=", m2
c	PRINT *, "m3=", m3

	call ctsxcut(imode, rootsvalue, muscale, 3, num, CTmpnum,
     &    rank, pj, mj, res(0), res(1), res(2), stable)
	CTCcut = res(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setscaloop(4)
	  call ctsxcut(imode, rootsvalue, muscale, 3, numeps, CTmpnum,
     &      rank, pj, mj, res(0), res(1), res(2), stable)
	  CTCcut = CTCcut - 2*res(1)
	  call setscaloop(1)
	endif
c	PRINT *, "Ccut=", CTCcut
	end

