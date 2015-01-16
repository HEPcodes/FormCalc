* SABcut.F
* the FormCalc wrapper for the Samurai two-point function
* this file is part of FormCalc
* last modified 13 Feb 12 th

#include "types.h"


	ComplexType function SABcut(rank, num, numeps, heldelta,
     &    p, m1, m2)
	implicit none
	integer rank, heldelta
	ComplexType num, numeps
	external num, numeps
	integer*8 p
	RealType m1, m2

#include "SA.h"

	RealType pj(2,4), mj(2)
	ComplexType tot(-2:0), totr
	logical ok

	SABcut = 0
	if( heldelta .eq. 0 ) return

	pj(1,1) = 0
	pj(1,2) = 0
	pj(1,3) = 0
	pj(1,4) = 0
	call SAMomDecode(p, pj(2,1), 2)

	mj(1) = m1
	mj(2) = m2

c	PRINT *, "SABcut:", rank
c	PRINT *, "p1=", pj(2,4)
c	PRINT *, "   ", pj(2,1)
c	PRINT *, "   ", pj(2,2)
c	PRINT *, "   ", pj(2,3)
c	PRINT *, "m1=", m1
c	PRINT *, "m2=", m2

	call samurai(num, tot, totr, pj, mj, 2, rank, 1, mudim, ok)
	SABcut = tot(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setisca(5)
	  call samurai(numeps, tot, totr, pj, mj, 2, rank, 1, mudim, ok)
	  SABcut = SABcut - 2*tot(-1)
	  call setisca(4)
	endif
c	PRINT *, "SABcut=", SABcut
	end

