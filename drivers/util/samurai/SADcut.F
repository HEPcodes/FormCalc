* SADcut.F
* the FormCalc wrapper for the Samurai four-point function
* this file is part of FormCalc
* last modified 13 Feb 12 th

#include "types.h"


	ComplexType function SADcut(rank, num, numeps, heldelta,
     &    p1, p2, p3, m1, m2, m3, m4)
	implicit none
	integer rank, heldelta
	ComplexType num, numeps
	external num, numeps
	integer*8 p1, p2, p3
	RealType m1, m2, m3, m4

#include "SA.h"

	RealType pj(4,4), mj(4)
	ComplexType tot(-2:0), totr
	logical ok

	SADcut = 0
	if( heldelta .eq. 0 ) return

	pj(1,1) = 0
	pj(1,2) = 0
	pj(1,3) = 0
	pj(1,4) = 0
	call SAMomDecode(p1, pj(2,1), 4)
	call SAMomDecode(p2, pj(3,1), 4)
	call SAMomDecode(p3, pj(4,1), 4)

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4

c	PRINT *, "SADcut:", rank
c	PRINT *, "p1=", pj(2,4)
c	PRINT *, "   ", pj(2,1)
c	PRINT *, "   ", pj(2,2)
c	PRINT *, "   ", pj(2,3)
c	PRINT *, "p2=", pj(3,4)
c	PRINT *, "   ", pj(3,1)
c	PRINT *, "   ", pj(3,2)
c	PRINT *, "   ", pj(3,3)
c	PRINT *, "p3=", pj(4,4)
c	PRINT *, "   ", pj(4,1)
c	PRINT *, "   ", pj(4,2)
c	PRINT *, "   ", pj(4,3)
c	PRINT *, "m1=", m1
c	PRINT *, "m2=", m2
c	PRINT *, "m3=", m3
c	PRINT *, "m4=", m4

	call samurai(num, tot, totr, pj, mj, 4, rank, 1, mudim, ok)
	SADcut = tot(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setisca(5)
	  call samurai(numeps, tot, totr, pj, mj, 4, rank, 1, mudim, ok)
	  SADcut = SADcut - 2*tot(-1)
	  call setisca(4)
	endif
c	PRINT *, "SADcut=", SADcut
	end

