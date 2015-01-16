* SADcut.F
* the FormCalc wrapper for the Samurai four-point function
* this file is part of FormCalc
* last modified 29 Aug 14 th

#include "types.h"


	subroutine SADmas(mi, p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	ComplexType mi(*)
	RealType p1, p2, p3, p4, p1p2, p2p3
	RealType m1, m2, m3, m4
	end


************************************************************************

	ComplexType function SADcut(mi, hel, rank, num, numeps,
     &    k1, k2, k3, m1, m2, m3, m4)
	implicit none
	ComplexType mi(*), num, numeps
	external num, numeps
	integer hel, rank
	integer*8 k1, k2, k3
	RealType m1, m2, m3, m4

#include "SA.h"

	RealType kj(4,4), mj(4)
	ComplexType res(-2:0), resr
	logical ok

	SADcut = 0
	if( hel .eq. 0 ) return

	kj(1,1) = 0
	kj(1,2) = 0
	kj(1,3) = 0
	kj(1,4) = 0
	call SAMomDecode(k1, kj(2,1), 4)
	call SAMomDecode(k2, kj(3,1), 4)
	call SAMomDecode(k3, kj(4,1), 4)

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4

#if SADEBUG >= 4
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
	print '("   k3: ",Z8.8)', k3
	print *, "  k3 =", kj(4,4)
	print *, "      ", kj(4,1)
	print *, "      ", kj(4,2)
	print *, "      ", kj(4,3)
	print *, "  m1 =", mj(1)
	print *, "  m2 =", mj(2)
	print *, "  m3 =", mj(3)
	print *, "  m4 =", mj(4)
	call flush(6)
#endif

	call fcsamurai(num, numeps, 4, rank, kj, mj, SADcut)

#if SADEBUG
	print *, "SADcut =", SADcut
	call flush(6)
#endif
	end

