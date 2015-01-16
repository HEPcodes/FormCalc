* Eps.F
* calculates -I times the Levi-Civita tensor contracted with
* four four-vectors
* this file is part of FormCalc
* last modified 14 Feb 12 th

#include "types.h"


* Note: the arguments a, b, c, d specify the vector's position in
* the array vec, NOT the vector itself.  The macros k, e, ec, s
* (defined in util.h) convert from 1...LEGS to a, b, c, d.

	ComplexType function Eps(a, b, c, d)
	implicit none
	integer a, b, c, d

	ComplexType vec(2,2, 1)
	common /vectors/ vec

#define A(x,y) vec(x,y, a)
#define B(x,y) vec(x,y, b)

	ComplexType diag, off, e
	diag(a,b) = A(1,1)*B(2,2) - A(2,2)*B(1,1)
	off(a,b) = A(2,1)*B(1,2) - A(1,2)*B(2,1)
	e(a,b,c,d) = diag(a,b)*off(c,d) + diag(c,d)*off(a,b)

	Eps = .25D0*(e(a,b,c,d) - e(a,c,b,d) + e(a,d,b,c))
	end

