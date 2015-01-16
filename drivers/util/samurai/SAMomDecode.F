* SAMomDecode.F
* form a linear combination of momenta in Samurai conventions
* and return it as a Minkowski four-vector
* the momenta's (integer*1) prefactors are encoded as the
* bytes of an integer*8
* this file is part of FormCalc
* last modified 8 May 13 th

#include "types.h"
#include "util.h"


	subroutine SAMomDecode(code, p,ldp)
	implicit none
	integer*8 code
	integer ldp
	RealType p(ldp,4)

#include "util.h"

#define K(x,y, i) vec(x,y, k0(i))

	integer i
	ComplexType a(2,2)

	integer*8 c
	integer*1 f

	a(1,1) = 0
	a(2,1) = 0
	a(1,2) = 0
	a(2,2) = 0

	c = code
	i = 0

	do while( c .ne. 0 )
	  i = i + 1
	  f = iand(c, 255)
	  a(1,1) = a(1,1) + f*K(1,1,i)
	  a(2,1) = a(2,1) + f*K(2,1,i)
	  a(1,2) = a(1,2) + f*K(1,2,i)
	  a(2,2) = a(2,2) + f*K(2,2,i)
	  c = ishft(c, -8)
	enddo

	p(1,4) = .5D0*Re(a(1,1) + a(2,2))
	p(1,3) = .5D0*Re(a(1,1) - a(2,2))
	p(1,1) = .5D0*Re(a(2,1) + a(1,2))
	p(1,2) = .5D0*Im(a(2,1) - a(1,2))
	end

