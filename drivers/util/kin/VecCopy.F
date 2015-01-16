* VecCopy.F
* set up the correct polarization vectors and spinors
* for a particular helicity combination
* this file is part of FormCalc
* last modified 8 May 13 th

#include "distrib.h"
#include "types.h"
#include "util.h"


* v: the SIMD-vector index
* n: the number of external particles
* Hel: the helicities

	subroutine VecCopy(v, n, Hel)
	implicit none
	integer v, n, Hel(*)

#if SIMD > 0

#include "util.h"

	integer i, a, b

	do i = 1, n
	  a = k0(i)
	  b = k(i)
	  ves(v,1,1,b) = vec(1,1,a)
	  ves(v,2,1,b) = vec(2,1,a)
	  ves(v,1,2,b) = vec(1,2,a)
	  ves(v,2,2,b) = vec(2,2,a)

	  a = s0(i)
	  b = s(i)
	  ves(v,1,1,b) = vec(1,1,a)
	  ves(v,2,1,b) = vec(2,1,a)
	  ves(v,1,2,b) = vec(1,2,a)
	  ves(v,2,2,b) = vec(2,2,a)

	  a = e0(i)
	  b = e(i)
	  ves(v,1,1,b) = vec(1,1,a)
	  ves(v,2,1,b) = vec(2,1,a)
	  ves(v,1,2,b) = vec(1,2,a)
	  ves(v,2,2,b) = vec(2,2,a)

	  a = ec0(i)
	  b = ec(i)
	  ves(v,1,1,b) = vec(1,1,a)
	  ves(v,2,1,b) = vec(2,1,a)
	  ves(v,1,2,b) = vec(1,2,a)
	  ves(v,2,2,b) = vec(2,2,a)

	  a = Spinor0(i,-1,1)
	  b = Spinor(i,-1,1)
	  ves(v,1,1,b) = vec(1,1,a)
	  ves(v,2,1,b) = vec(2,1,a)
	  ves(v,1,2,b) = vec(1,2,a)
	  ves(v,2,2,b) = vec(2,2,a)

	  a = Spinor0(i,-1,2)
	  b = Spinor(i,-1,2)
	  ves(v,1,1,b) = vec(1,1,a)
	  ves(v,2,1,b) = vec(2,1,a)
	  ves(v,1,2,b) = vec(1,2,a)
	  ves(v,2,2,b) = vec(2,2,a)

	  a = Spinor0(i,+1,1)
	  b = Spinor(i,+1,1)
	  ves(v,1,1,b) = vec(1,1,a)
	  ves(v,2,1,b) = vec(2,1,a)
	  ves(v,1,2,b) = vec(1,2,a)
	  ves(v,2,2,b) = vec(2,2,a)

	  a = Spinor0(i,+1,2)
	  b = Spinor(i,+1,2)
	  ves(v,1,1,b) = vec(1,1,a)
	  ves(v,2,1,b) = vec(2,1,a)
	  ves(v,1,2,b) = vec(1,2,a)
	  ves(v,2,2,b) = vec(2,2,a)
	enddo

#endif
	end

