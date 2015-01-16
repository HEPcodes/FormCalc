* CTmpnum.F
* a dummy multiprecision numerator function
* this file is part of FormCalc
* last modified 13 Feb 12 th

#include "types.h"


	subroutine CTmpnum(q1in, res)
	implicit none
	RealType q1in(0:3)
	ComplexType res

	integer ini
	data ini /0/
	save ini

	if( ini .eq. 0 ) then
	  print *, "MP version not yet available"
	  ini = 1
	endif

	res = 0
	end

