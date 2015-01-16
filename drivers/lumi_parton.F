* lumi_parton.F
* "luminosity" for a partonic process
* this file is part of FormCalc
* last modified 22 May 08 th


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
#ifdef FORCE_ONSHELL
	  if( abs(mass_in - Var(SQRTS)) .lt. 1D-9 ) then
	    if( fail .eq. 1 ) INFO "Below threshold"
	    return
	  endif
	  Var(SQRTS) = mass_in
#else
	  if( fail .ne. 1 ) return
	  Var(SQRTS) = threshold + .01D0
	  INFO "Adjusting sqrtS to above-threshold value ", Var(SQRTS)
#endif
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

