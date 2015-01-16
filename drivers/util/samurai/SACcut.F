* SACcut.F
* the FormCalc wrapper for the Samurai three-point function
* this file is part of FormCalc
* last modified 6 Oct 12 th

#include "types.h"


	subroutine SACmas(mi, p1, p2, p1p2, m1, m2, m3)
	implicit none
	ComplexType mi(*)
	RealType p1, p2, p1p2, m1, m2, m3
	end


************************************************************************

	ComplexType function SACcut(mi, hel, rank, num, numeps,
     &    k1, k2, m1, m2, m3)
	implicit none
	ComplexType mi(*), num, numeps
	external num, numeps
	integer hel, rank
	integer*8 k1, k2
	RealType m1, m2, m3

#include "SA.h"

	RealType kj(3,4), mj(3)
	ComplexType tot(-2:0), totr
	logical ok

	SACcut = 0
	if( hel .eq. 0 ) return

	kj(1,1) = 0
	kj(1,2) = 0
	kj(1,3) = 0
	kj(1,4) = 0
	call SAMomDecode(k1, kj(2,1), 3)
	call SAMomDecode(k2, kj(3,1), 3)

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3

c	PRINT *, "SACcut:", rank
c	PRINT *, "k1=", kj(2,4)
c	PRINT *, "   ", kj(2,1)
c	PRINT *, "   ", kj(2,2)
c	PRINT *, "   ", kj(2,3)
c	PRINT *, "k2=", kj(3,4)
c	PRINT *, "   ", kj(3,1)
c	PRINT *, "   ", kj(3,2)
c	PRINT *, "   ", kj(3,3)
c	PRINT *, "m1=", m1
c	PRINT *, "m2=", m2
c	PRINT *, "m3=", m3

	call samurai(num, tot, totr, kj, mj, 3, rank, 1, mudim, ok)
	SACcut = tot(epscoeff)
	if( epscoeff .eq. 0 ) then
	  call setisca(5)
	  call samurai(numeps, tot, totr, kj, mj, 3, rank, 1, mudim, ok)
	  SACcut = SACcut - 2*tot(-1)
	  call setisca(4)
	endif

c	PRINT *, "SACcut=", SACcut
	end

