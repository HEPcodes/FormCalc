* SABcut.F
* the FormCalc wrapper for the Samurai two-point function
* this file is part of FormCalc
* last modified 16 Jun 16 th

#include "types.h"


	subroutine SABmas(mi, p, m1, m2)
	implicit none
	ComplexType mi(*)
	RealType p, m1, m2
	end


************************************************************************

	ComplexType function SABcut(mi, hdelta, rank, num, numeps,
     &    k, m1, m2)
	implicit none
	ComplexType mi(*), num, numeps
	external num, numeps
	integer hdelta, rank
	integer*8 k
	RealType m1, m2

#include "SA.h"

	RealType kj(2,4), mj(2)

	SABcut = 0
	if( hdelta .eq. 0 ) return

	kj(1,1) = 0
	kj(1,2) = 0
	kj(1,3) = 0
	kj(1,4) = 0
	call SAMomDecode(k, kj(2,1), 2)

	mj(1) = m1
	mj(2) = m2

#if SADEBUG >= 2
	print *, "  rank =", rank
	print '("   k1: ",Z8.8)', k1
	print *, "  k1 =", kj(2,4)
	print *, "      ", kj(2,1)
	print *, "      ", kj(2,2)
	print *, "      ", kj(2,3)
	print *, "  m1 =", mj(1)
	print *, "  m2 =", mj(2)
	call flush(6)
#endif

	call fcsamurai(num, numeps, 2, rank, kj, mj, SABcut)

#if SADEBUG
	print *, "SABcut =", SABcut
	call flush(6)
#endif
	end

