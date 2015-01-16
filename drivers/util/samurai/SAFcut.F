* SAFcut.F
* the FormCalc wrapper for the Samurai six-point function
* this file is part of FormCalc
* last modified 29 Aug 14 th

#include "types.h"


	subroutine SAFmas(mi, p1, p2, p3, p4, p5, p6,
     &    p1p2, p2p3, p3p4, p4p5, p5p6, p6p1,
     &    p1p2p3, p2p3p4, p3p4p5,
     &    m1, m2, m3, m4, m5, m6)
	implicit none
	ComplexType mi(*)
	RealType p1, p2, p3, p4, p5, p6
	RealType p1p2, p2p3, p3p4, p4p5, p5p6, p6p1
	RealType p1p2p3, p2p3p4, p3p4p5
	RealType m1, m2, m3, m4, m5, m6
	end


************************************************************************

	ComplexType function SAFcut(mi, hel, rank, num, numeps,
     &    k1, k2, k3, k4, k5, m1, m2, m3, m4, m5, m6)
	implicit none
	ComplexType mi(*), num, numeps
	external num, numeps
	integer hel, rank
	integer*8 k1, k2, k3, k4, k5
	RealType m1, m2, m3, m4, m5, m6

#include "SA.h"

	RealType kj(6,4), mj(6)
	ComplexType res(-2:0), resr
	logical ok

	SAFcut = 0
	if( hel .eq. 0 ) return

	kj(1,1) = 0
	kj(1,2) = 0
	kj(1,3) = 0
	kj(1,4) = 0
	call SAMomDecode(k1, kj(2,1),6)
	call SAMomDecode(k2, kj(3,1),6)
	call SAMomDecode(k3, kj(4,1),6)
	call SAMomDecode(k4, kj(5,1),6)
	call SAMomDecode(k5, kj(6,1),6)

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4
	mj(5) = m5
	mj(6) = m6

#if SADEBUG >= 6
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
	print '("   k4: ",Z8.8)', k4
	print *, "  k4 =", kj(5,4)
	print *, "      ", kj(5,1)
	print *, "      ", kj(5,2)
	print *, "      ", kj(5,3)
	print '("   k5: ",Z8.8)', k5
	print *, "  k5 =", kj(6,4)
	print *, "      ", kj(6,1)
	print *, "      ", kj(6,2)
	print *, "      ", kj(6,3)
	print *, "  m1 =", mj(1)
	print *, "  m2 =", mj(2)
	print *, "  m3 =", mj(3)
	print *, "  m4 =", mj(4)
	print *, "  m5 =", mj(5)
	call flush(6)
#endif

	call fcsamurai(num, numeps, 6, rank, kj, mj, SAFcut)

#if SADEBUG
	print *, "SAFcut =", SAFcut
	call flush(6)
#endif
	end

