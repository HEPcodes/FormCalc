* lumi_parton.F
* "luminosity" for a partonic process
* this file is part of FormCalc
* last modified 29 Mar 06 th


#define MINVAR SQRTS


************************************************************************

	subroutine LumiDefaults
	implicit none

	end


************************************************************************

	subroutine LumiIni(fail)
	implicit none
	integer fail

#include "xsection.h"

	integer i

	do i = 1, 10
	  call SetEnergy(fail, Var(SQRTS))
	  if( fail .ne. 1 ) return
	  Var(SQRTS) = threshold + .01D0
	  INFO "Adjusting sqrtS to above-threshold value ", Var(SQRTS)
	enddo

	WARN "Calculation of threshold did not converge"
	end


************************************************************************

	subroutine Luminosity(fac)
	implicit none
	double precision fac

	integer fail

#include "xsection.h"

	call SetEnergy(fail, Var(SQRTS))

	if( fail .ne. 0 ) fac = 0
	end

