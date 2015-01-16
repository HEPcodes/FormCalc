* Chain.F
* compute a spinor chain specified through an opcode string
* this file is part of FormCalc
* last modified 27 Feb 12 th

#include "types.h"


	ComplexType function Chain(t, opcodes)
	implicit none
	integer t
	integer*8 opcodes

	integer b, e
	ComplexType s1, s2, x1

	integer*8 ops
	integer*1 n
	equivalence (n, ops)

	ComplexType vec(0:1,0:1, 1)
	common /vectors/ vec

	ComplexType spi(0:1, 0:1, 1)
	equivalence (spi, vec)

#define sgn(i) (1-2*i)

	ops = opcodes
	b = iand(t, 7) - 6
* must correspond to EpsL in util.h:
	e = ibits(t, 3, 1)

*     a = { a1, a2}
* a.eps = {-a2, a1}
	s1 = sgn(e)*spi(e,ieor(e,1-b),n)
	s2 = spi(1-e,ieor(e,1-b),n)

	ops = ishft(ops, -8)
	do while( ishft(ops, -8) .ne. 0 )
	  x1 = vec(b,b,n)*s1 + sgn(b)*vec(1,0,n)*s2
	  s2 = sgn(b)*vec(0,1,n)*s1 + vec(1-b,1-b,n)*s2
	  s1 = x1

*  SxV(1) =  A11 b1 + A21 b2
*  SxV(2) =  A12 b1 + A22 b2
*
*  SxB(1) =  A22 b1 - A21 b2
*  SxB(2) = -A12 b1 + A11 b2

	  b = 1 - b
	  ops = ishft(ops, -8)
	enddo

* must correspond to EpsR in util.h:
	e = ibits(t, 4, 1)
	b = ieor(b, e)
* SxS = a1*b1 + a2*b2
* SeS = a1*b2 - a2*b1
	Chain = spi(e,b,n)*s1 + sgn(e)*spi(1-e,b,n)*s2
	end

