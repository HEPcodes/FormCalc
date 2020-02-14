* CTBcut.F
* the FormCalc wrapper for the CutTools two-point function
* this file is part of FormCalc
* last modified 18 Mar 19 th

#include "types.h"


	subroutine CTBmas(mi, p, m1, m2)
	implicit none
	ComplexType mi(*)
	RealType p, m1, m2
	end


************************************************************************

	ComplexType function CTBcut(mi, hdelta, rank, num, numeps,
     &    k, m1, m2)
	implicit none
	ComplexType mi(*)
	integer hdelta, rank
	external num, numeps
	integer*8 k
	RealType m1, m2

#include "CT.h"

	RealType kj(0:3,2)
	ComplexType mj(2)
	logical stable
	ComplexType res(0:2)
	external CTmpnum

	CTBcut = 0
	if( hdelta == 0 ) return

	kj(:,1) = 0
	call CTMomDecode(k, kj(0,2))

	mj = [m1, m2]

#if CTDEBUG >= 2
	print *, "  rank =", rank
	print *, "  k1 =", kj(0,2)
	print *, "      ", kj(1,2)
	print *, "      ", kj(2,2)
	print *, "      ", kj(3,2)
	print *, "  m1 =", mj(1)
	print *, "  m2 =", mj(2)
	call flush(6)
#endif

	call ctsxcut(imode, rootsvalue, muscale, 2, num, CTmpnum,
     &    rank, kj, mj, res(0), res(1), res(2), stable)
	CTBcut = res(epsi)
	if( epsi == 0 ) then
	  call setscaloop(4)
	  call ctsxcut(imode, rootsvalue, muscale, 2, numeps, CTmpnum,
     &      rank, kj, mj, res(0), res(1), res(2), stable)
	  CTBcut = CTBcut - 2*res(1)
	  call setscaloop(1)
	endif

#if CTDEBUG
	print *, "CTBcut =", CTBcut
	call flush(6)
#endif
	end

