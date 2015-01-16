* SxS.F
* computes Spinor.Spinor
* this file is part of FormCalc
* last modified 9 Aug 11 th

#include "types.h"


	Complex function SxS(a, b)
	implicit none
	integer a, b

	Complex spi(2, 1)
	common /vectors/ spi

#define A(x) spi(x, a)
#define B(x) spi(x, b)

	SxS = A(1)*B(1) + A(2)*B(2)
	end

