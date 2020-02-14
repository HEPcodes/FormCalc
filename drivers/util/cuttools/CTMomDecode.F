* CTMomDecode.F
* form a linear combination of momenta in CutTools conventions
* and return it as a Minkowski four-vector
* the momenta's (integer*1) prefactors are encoded as the
* bytes of an integer*8
* this file is part of FormCalc
* last modified 22 Aug 19 th

#include "types.h"
#include "util.h"


	subroutine CTMomDecode(code, p)
	implicit none
	integer*8 code
	RealType p(0:3)

#include "util.h"

#define K(x,y, i) vec0(x,y, k0(i))

	integer i
	integer*8 c
	ComplexType a(2,2)

	a = 0
	c = code
	i = 0

	do while( c /= 0 )
	  i = i + 1
	  a = a + iand(c, 255_8)*K(:,:,i)
	  c = ishft(c, -8)
	enddo

	p(0) = Re(a(1,1) + a(2,2))/2
	p(3) = Re(a(1,1) - a(2,2))/2
	p(1) = Re(a(2,1) + a(1,2))/2
	p(2) = Im(a(2,1) - a(1,2))/2
	end

