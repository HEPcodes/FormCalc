* VecCopy.F
* set up the correct polarization vectors and spinors
* for a particular helicity combination
* this file is part of FormCalc
* last modified 11 Jan 16 th

#include "distrib.h"
#include "types.h"
#include "util.h"


* v: the SIMD-vector index
* n: the number of external particles
* Hel: the helicities

	subroutine VecCopy(v, n)
	implicit none
	integer v, n

#if SIMD > 1

#include "util.h"

	integer i, a, b

	do i = 1, n
	  Hel(v,i) = Hel0(i)

	  a = k0(i)
	  b = k(i)
	  vec(v,1,1,b) = vec0(1,1,a)
	  vec(v,2,1,b) = vec0(2,1,a)
	  vec(v,1,2,b) = vec0(1,2,a)
	  vec(v,2,2,b) = vec0(2,2,a)

	  a = s0(i)
	  b = s(i)
	  vec(v,1,1,b) = vec0(1,1,a)
	  vec(v,2,1,b) = vec0(2,1,a)
	  vec(v,1,2,b) = vec0(1,2,a)
	  vec(v,2,2,b) = vec0(2,2,a)

	  a = e0(i)
	  b = e(i)
	  vec(v,1,1,b) = vec0(1,1,a)
	  vec(v,2,1,b) = vec0(2,1,a)
	  vec(v,1,2,b) = vec0(1,2,a)
	  vec(v,2,2,b) = vec0(2,2,a)

	  a = ec0(i)
	  b = ec(i)
	  vec(v,1,1,b) = vec0(1,1,a)
	  vec(v,2,1,b) = vec0(2,1,a)
	  vec(v,1,2,b) = vec0(1,2,a)
	  vec(v,2,2,b) = vec0(2,2,a)

	  a = Spinor0(i,-1,1)
	  b = Spinor(i,-1,1)
	  vec(v,1,1,b) = vec0(1,1,a)
	  vec(v,2,1,b) = vec0(2,1,a)
	  vec(v,1,2,b) = vec0(1,2,a)
	  vec(v,2,2,b) = vec0(2,2,a)

	  a = Spinor0(i,-1,2)
	  b = Spinor(i,-1,2)
	  vec(v,1,1,b) = vec0(1,1,a)
	  vec(v,2,1,b) = vec0(2,1,a)
	  vec(v,1,2,b) = vec0(1,2,a)
	  vec(v,2,2,b) = vec0(2,2,a)

	  a = Spinor0(i,+1,1)
	  b = Spinor(i,+1,1)
	  vec(v,1,1,b) = vec0(1,1,a)
	  vec(v,2,1,b) = vec0(2,1,a)
	  vec(v,1,2,b) = vec0(1,2,a)
	  vec(v,2,2,b) = vec0(2,2,a)

	  a = Spinor0(i,+1,2)
	  b = Spinor(i,+1,2)
	  vec(v,1,1,b) = vec0(1,1,a)
	  vec(v,2,1,b) = vec0(2,1,a)
	  vec(v,1,2,b) = vec0(1,2,a)
	  vec(v,2,2,b) = vec0(2,2,a)
	enddo

#endif
	end

