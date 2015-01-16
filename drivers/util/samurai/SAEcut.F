* SAEcut.F
* the FormCalc wrapper for the Samurai five-point function
* this file is part of FormCalc
* last modified 6 Oct 12 th

#include "types.h"


	subroutine SAEmas(mi, p1, p2, p3, p4, p5,
     &    p1p2, p2p3, p3p4, p4p5, p5p1,
     &    m1, m2, m3, m4, m5)
	implicit none
	ComplexType mi(*)
	RealType p1, p2, p3, p4, p5
	RealType p1p2, p2p3, p3p4, p4p5, p5p1
	RealType m1, m2, m3, m4, m5
	end


************************************************************************

	ComplexType function SAEcut(mi, hel, rank, num, numeps,
     &    k1, k2, k3, k4, m1, m2, m3, m4, m5)
	implicit none
	ComplexType mi(*), num, numeps
	external num, numeps
	integer hel, rank
	integer*8 k1, k2, k3, k4
	RealType m1, m2, m3, m4, m5

#include "SA.h"

	RealType kj(5,4), mj(5)
	ComplexType tot(-2:0), totr
	logical ok

	SAEcut = 0
	if( hel .eq. 0 ) return

	kj(1,1) = 0
	kj(1,2) = 0
	kj(1,3) = 0
	kj(1,4) = 0
	call SAMomDecode(k1, kj(2,1), 5)
	call SAMomDecode(k2, kj(3,1), 5)
	call SAMomDecode(k3, kj(4,1), 5)
	call SAMomDecode(k4, kj(5,1), 5)

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4
	mj(5) = m5

	call samurai(num, tot, totr, kj, mj, 5, rank, 1, mudim, ok)
	SAEcut = tot(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setisca(5)
	  call samurai(numeps, tot, totr, kj, mj, 5, rank, 1, mudim, ok)
	  SAEcut = SAEcut - 2*tot(-1)
	  call setisca(4)
	endif
	end

