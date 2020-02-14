* CTFcut.F
* the FormCalc wrapper for the CutTools six-point function
* this file is part of FormCalc
* last modified 18 Mar 19 th

#include "types.h"


	subroutine CTFmas(mi, p1, p2, p3, p4, p5, p6,
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

	ComplexType function CTFcut(mi, hdelta, rank, num, numeps,
     &    k1, k2, k3, k4, k5, m1, m2, m3, m4, m5, m6)
	implicit none
	ComplexType mi(*)
	integer hdelta, rank
	external num, numeps
	integer*8 k1, k2, k3, k4, k5
	RealType m1, m2, m3, m4, m5, m6

#include "CT.h"

	RealType kj(0:3,6)
	ComplexType mj(6)
	logical stable
	ComplexType res(0:2)
	external CTmpnum

	CTFcut = 0
	if( hdelta == 0 ) return

	kj(:,1) = 0
	call CTMomDecode(k1, kj(0,2))
	call CTMomDecode(k2, kj(0,3))
	call CTMomDecode(k3, kj(0,4))
	call CTMomDecode(k4, kj(0,5))
	call CTMomDecode(k5, kj(0,6))

	mj = [m1, m2, m3, m4, m5, m6]

#if CTDEBUG >= 2
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
	print *, "  k4 =", kj(0,5)
	print *, "      ", kj(1,5)
	print *, "      ", kj(2,5)
	print *, "      ", kj(3,5)
	print *, "  k5 =", kj(0,6)
	print *, "      ", kj(1,6)
	print *, "      ", kj(2,6)
	print *, "      ", kj(3,6)
	print *, "  m1 =", mj(1)
	print *, "  m2 =", mj(2)
	print *, "  m3 =", mj(3)
	print *, "  m4 =", mj(4)
	print *, "  m5 =", mj(5)
	print *, "  m6 =", mj(6)
	call flush(6)
#endif

	call ctsxcut(imode, rootsvalue, muscale, 6, num, CTmpnum,
     &    rank, kj, mj, res(0), res(1), res(2), stable)
	CTFcut = res(epsi)
	if( epsi == 0 ) then
	  call setscaloop(4)
	  call ctsxcut(imode, rootsvalue, muscale, 6, numeps, CTmpnum,
     &      rank, kj, mj, res(0), res(1), res(2), stable)
	  CTFcut = CTFcut - 2*res(1)
	  call setscaloop(1)
	endif

#if CTDEBUG
	print *, "CTFcut =", CTFcut
	call flush(6)
#endif
	end

