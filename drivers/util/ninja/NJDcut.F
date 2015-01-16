* NJDcut.F
* the FormCalc wrapper for the Ninja four-point function
* this file is part of FormCalc
* last modified 29 Aug 14 th

#include "types.h"


	subroutine NJDmas(mi, p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	ComplexType mi(*)
	RealType p1, p2, p3, p4, p1p2, p2p3
	RealType m1, m2, m3, m4
	end


************************************************************************

	ComplexType function NJDcut(mi, hel, rank,
     &    num, muexp, t3exp, t2exp,
     &	  k1, k2, k3, m1, m2, m3, m4)
	implicit none
	ComplexType mi(*)
	external num, muexp, t3exp, t2exp
	integer hel, rank
	integer*8 k1, k2, k3
	RealType m1, m2, m3, m4

#include "NJ.h"

	RealType kj(0:3,4), mj(0:3)
	ComplexType res(-2:0)

	NJDcut = 0
	if( hel .eq. 0 ) return

	kj(0,1) = 0
	kj(1,1) = 0
	kj(2,1) = 0
	kj(3,1) = 0
	call NJMomDecode(k1, kj(0,2))
	call NJMomDecode(k2, kj(0,3))
	call NJMomDecode(k3, kj(0,4))

	mj(0) = m1
	mj(1) = m2
	mj(2) = m3
	mj(3) = m4

#if NJDEBUG >= 4
	print *, "  rank =", rank
	print '("   k1: ",Z8.8)', k1
	print *, "  k1 =", kj(0,2)
	print *, "      ", kj(1,2)
	print *, "      ", kj(2,2)
	print *, "      ", kj(3,2)
	print '("   k2: ",Z8.8)', k2
	print *, "  k2 =", kj(0,3)
	print *, "      ", kj(1,3)
	print *, "      ", kj(2,3)
	print *, "      ", kj(3,3)
	print '("   k3: ",Z8.8)', k3
	print *, "  k3 =", kj(0,4)
	print *, "      ", kj(1,4)
	print *, "      ", kj(2,4)
	print *, "      ", kj(3,4)
	print *, "  m1 =", mj(1)
	print *, "  m2 =", mj(2)
	print *, "  m3 =", mj(3)
	print *, "  m4 =", mj(4)
	call flush(6)
#endif

	call ninjaformcalc(num, muexp, t3exp, t2exp, 4, rank,
     &    kj, mj, mudim, res)
	NJDcut = res(-epsi)

#if NJDEBUG
	print *, "NJDcut =", NJDcut
	call flush(6)
#endif
	end

