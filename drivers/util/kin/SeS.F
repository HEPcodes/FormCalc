* SeS.F
* computes Spinor.eps.Spinor
* this file is part of FormCalc
* last modified 9 Aug 11 th

#include "types.h"


	Complex function SeS(a, b)
	implicit none
	integer a, b

	Complex spi(2, 1)
	common /vectors/ spi

#define A(x) spi(x, a)
#define B(x) spi(x, b)

	SeS = A(1)*B(2) - A(2)*B(1)
	end

