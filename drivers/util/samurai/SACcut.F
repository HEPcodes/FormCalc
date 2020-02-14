* SACcut.F
* the FormCalc wrapper for the Samurai three-point function
* this file is part of FormCalc
* last modified 18 Mar 19 th

#include "types.h"


	subroutine SACmas(mi, p1, p2, p1p2, m1, m2, m3)
	implicit none
	ComplexType mi(*)
	RealType p1, p2, p1p2, m1, m2, m3
	end


************************************************************************

	ComplexType function SACcut(mi, hdelta, rank, num, numeps,
     &    k1, k2, m1, m2, m3)
	implicit none
	ComplexType mi(*), num, numeps
	external num, numeps
	integer hdelta, rank
	integer*8 k1, k2
	RealType m1, m2, m3

#include "SA.h"

	RealType kj(3,4), mj(3)

	SACcut = 0
	if( hdelta == 0 ) return

	kj(1,:) = 0
	call SAMomDecode(k1, kj(2,1),3)
	call SAMomDecode(k2, kj(3,1),3)

	mj = [m1, m2, m3]

#if SADEBUG >= 3
	print *, "  rank =", rank
	print '("   k1: ",Z8.8)', k1
	print *, "  k1 =", kj(2,4)
	print *, "      ", kj(2,1)
	print *, "      ", kj(2,2)
	print *, "      ", kj(2,3)
	print '("   k2: ",Z8.8)', k2
	print *, "  k2 =", kj(3,4)
	print *, "      ", kj(3,1)
	print *, "      ", kj(3,2)
	print *, "      ", kj(3,3)
	print *, "  m1 =", mj(1)
	print *, "  m2 =", mj(2)
	print *, "  m3 =", mj(3)
	call flush(6)
#endif

	call fcsamurai(num, numeps, 3, rank, kj, mj, SACcut)

#if SADEBUG
	print *, "SACcut =", SACcut
	call flush(6)
#endif
	end

