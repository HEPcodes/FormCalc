* SAEcut.F
* the FormCalc wrapper for the Samurai five-point function
* this file is part of FormCalc
* last modified 18 Mar 19 th

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

	ComplexType function SAEcut(mi, hdelta, rank, num, numeps,
     &    k1, k2, k3, k4, m1, m2, m3, m4, m5)
	implicit none
	ComplexType mi(*), num, numeps
	external num, numeps
	integer hdelta, rank
	integer*8 k1, k2, k3, k4
	RealType m1, m2, m3, m4, m5

#include "SA.h"

	RealType kj(5,4), mj(5)

	SAEcut = 0
	if( hdelta == 0 ) return

	kj(1,:) = 0
	call SAMomDecode(k1, kj(2,1),5)
	call SAMomDecode(k2, kj(3,1),5)
	call SAMomDecode(k3, kj(4,1),5)
	call SAMomDecode(k4, kj(5,1),5)

	mj = [m1, m2, m3, m4, m5]

#if SADEBUG >= 5
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
	print *, "  m1 =", mj(1)
	print *, "  m2 =", mj(2)
	print *, "  m3 =", mj(3)
	print *, "  m4 =", mj(4)
	print *, "  m5 =", mj(5)
	call flush(6)
#endif

	call fcsamurai(num, numeps, 5, rank, kj, mj, SAEcut)

#if SADEBUG
	print *, "SAEcut =", SAEcut
	call flush(6)
#endif
	end

