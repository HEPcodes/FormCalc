* SABcut.F
* the FormCalc wrapper for the Samurai two-point function
* this file is part of FormCalc
* last modified 6 Oct 12 th

#include "types.h"


	subroutine SABmas(mi, p, m1, m2)
	implicit none
	ComplexType mi(*)
	RealType p, m1, m2
	end


************************************************************************

	ComplexType function SABcut(mi, hel, rank, num, numeps,
     &    k, m1, m2)
	implicit none
	ComplexType mi(*), num, numeps
	external num, numeps
	integer hel, rank
	integer*8 k
	RealType m1, m2

#include "SA.h"

	RealType kj(2,4), mj(2)
	ComplexType tot(-2:0), totr
	logical ok

	SABcut = 0
	if( hel .eq. 0 ) return

	kj(1,1) = 0
	kj(1,2) = 0
	kj(1,3) = 0
	kj(1,4) = 0
	call SAMomDecode(k, kj(2,1), 2)

	mj(1) = m1
	mj(2) = m2

c	PRINT *, "SABcut:", rank
c	PRINT *, "k1=", kj(2,4)
c	PRINT *, "   ", kj(2,1)
c	PRINT *, "   ", kj(2,2)
c	PRINT *, "   ", kj(2,3)
c	PRINT *, "m1=", m1
c	PRINT *, "m2=", m2

	call samurai(num, tot, totr, kj, mj, 2, rank, 1, mudim, ok)
	SABcut = tot(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setisca(5)
	  call samurai(numeps, tot, totr, kj, mj, 2, rank, 1, mudim, ok)
	  SABcut = SABcut - 2*tot(-1)
	  call setisca(4)
	endif
c	PRINT *, "SABcut=", SABcut
	end

