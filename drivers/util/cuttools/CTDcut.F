* CTDcut.F
* the FormCalc wrapper for the CutTools four-point function
* this file is part of FormCalc
* last modified 4 Jul 14 th

#include "types.h"


	subroutine CTDmas(mi, p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	ComplexType mi(*)
	RealType p1, p2, p3, p4, p1p2, p2p3
	RealType m1, m2, m3, m4
	end


************************************************************************

	ComplexType function CTDcut(mi, hel, rank, num, numeps,
     &    k1, k2, k3, m1, m2, m3, m4)
	implicit none
	ComplexType mi(*)
	integer hel, rank
	external num, numeps
	integer*8 k1, k2, k3
	RealType m1, m2, m3, m4

#include "CT.h"

	RealType kj(0:3,4)
	ComplexType mj(4)
	logical stable
	ComplexType res(0:2)
	external CTmpnum

	CTDcut = 0
	if( hel .eq. 0 ) return

	kj(0,1) = 0
	kj(1,1) = 0
	kj(2,1) = 0
	kj(3,1) = 0
	call CTMomDecode(k1, kj(0,2))
	call CTMomDecode(k2, kj(0,3))
	call CTMomDecode(k3, kj(0,4))

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4

#if CTDEBUG >= 4
	print *, "  rank =", rank
	print *, "  k1 =", kj(0,2)
	print *, "      ", kj(1,2)
	print *, "      ", kj(2,2)
	print *, "      ", kj(3,2)
	print *, "  k2 =", kj(0,3)
	print *, "      ", kj(1,3)
	print *, "      ", kj(2,3)
	print *, "      ", kj(3,3)
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

	call ctsxcut(imode, rootsvalue, muscale, 4, num, CTmpnum,
     &    rank, kj, mj, res(0), res(1), res(2), stable)
	CTDcut = res(epsi)
	if( epsi .eq. 0 ) then
	  call setscaloop(4)
	  call ctsxcut(imode, rootsvalue, muscale, 4, numeps, CTmpnum,
     &      rank, kj, mj, res(0), res(1), res(2), stable)
	  CTDcut = CTDcut - 2*res(1)
	  call setscaloop(1)
	endif

#if CTDEBUG
	print *, "CTDcut =", CTDcut
	call flush(6)
#endif
	end

