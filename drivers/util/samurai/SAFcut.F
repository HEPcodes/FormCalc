* SAFcut.F
* the FormCalc wrapper for the Samurai six-point function
* this file is part of FormCalc
* last modified 31 Mar 11 th

#include "types.h"


	ComplexType function SAFcut(rank, num, numeps, heldelta,
     &    p1, p2, p3, p4, p5, m1, m2, m3, m4, m5, m6)
	implicit none
	integer rank, heldelta
	ComplexType num, numeps
	external num, numeps
	integer*8 p1, p2, p3, p4, p5
	RealType m1, m2, m3, m4, m5, m6

#include "SA.h"

	RealType pj(6,4), mj(6)
	ComplexType tot(-2:0), totr
	logical ok

	SAFcut = 0
	if( heldelta .eq. 0 ) return

	pj(1,1) = 0
	pj(1,2) = 0
	pj(1,3) = 0
	pj(1,4) = 0
	call SAMomDecode(p1, pj(2,1), 6)
	call SAMomDecode(p2, pj(3,1), 6)
	call SAMomDecode(p3, pj(4,1), 6)
	call SAMomDecode(p4, pj(5,1), 6)
	call SAMomDecode(p5, pj(6,1), 6)

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4
	mj(5) = m5
	mj(6) = m6

	call samurai(num, tot, totr, pj, mj, 6, rank, 1, mudim, ok)
	SAFcut = tot(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setisca(5)
	  call samurai(numeps, tot, totr, pj, mj, 6, rank, 1, mudim, ok)
	  SAFcut = SAFcut - 2*tot(-1)
	  call setisca(4)
	endif
	end

