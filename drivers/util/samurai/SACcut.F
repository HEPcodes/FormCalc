* SACcut.F
* the FormCalc wrapper for the Samurai three-point function
* this file is part of FormCalc
* last modified 13 Feb 12 th

#include "types.h"


	ComplexType function SACcut(rank, num, numeps, heldelta,
     &    p1, p2, m1, m2, m3)
	implicit none
	integer rank, heldelta
	ComplexType num, numeps
	external num, numeps
	integer*8 p1, p2
	RealType m1, m2, m3

#include "SA.h"

	RealType pj(3,4), mj(3)
	ComplexType tot(-2:0), totr
	logical ok

	SACcut = 0
	if( heldelta .eq. 0 ) return

	pj(1,1) = 0
	pj(1,2) = 0
	pj(1,3) = 0
	pj(1,4) = 0
	call SAMomDecode(p1, pj(2,1), 3)
	call SAMomDecode(p2, pj(3,1), 3)

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3

c	PRINT *, "SACcut:", rank
c	PRINT *, "p1=", pj(2,4)
c	PRINT *, "   ", pj(2,1)
c	PRINT *, "   ", pj(2,2)
c	PRINT *, "   ", pj(2,3)
c	PRINT *, "p2=", pj(3,4)
c	PRINT *, "   ", pj(3,1)
c	PRINT *, "   ", pj(3,2)
c	PRINT *, "   ", pj(3,3)
c	PRINT *, "m1=", m1
c	PRINT *, "m2=", m2
c	PRINT *, "m3=", m3

	call samurai(num, tot, totr, pj, mj, 3, rank, 1, mudim, ok)
	SACcut = tot(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setisca(5)
	  call samurai(numeps, tot, totr, pj, mj, 3, rank, 1, mudim, ok)
	  SACcut = SACcut - 2*tot(-1)
	  call setisca(4)
	endif

c	PRINT *, "SACcut=", SACcut
	end

