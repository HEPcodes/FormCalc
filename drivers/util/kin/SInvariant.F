* SInvariant.F
* calculates an s-type invariant, i.e. s_{ij} = (k_i + k_j)^2
* this file is part of FormCalc
* last modified 9 Aug 11 th

#include "types.h"


	Real function SInvariant(i, j)
	implicit none
	integer i, j

	Complex vec(2,2, 8, 0:1)
	common /vectors/ vec

#define Ki(x,y) vec(x,y, 1, i)
#define Kj(x,y) vec(x,y, 1, j)

	Real S11, S22
	Complex S12

	S11 = Re(Ki(1,1)) + Re(Kj(1,1))
	S22 = Re(Ki(2,2)) + Re(Kj(2,2))
	S12 = Ki(1,2) + Kj(1,2)

	SInvariant = S11*S22 - Sq(S12)
	end

