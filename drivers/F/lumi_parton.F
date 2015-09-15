* lumi_parton.F
* "luminosity" for a partonic process
* this file is part of FormCalc
* last modified 11 Jun 15 th


#define MINVAR SQRTS


************************************************************************

	subroutine LumiDefaults(argc, argv)
	implicit none
	integer argc
	character*128 argv(*)

	end


************************************************************************

	subroutine LumiIni(fail, sqme)
	implicit none
	integer fail
	external sqme

#include "xsection.h"

	integer i

	do i = 1, 10
	  call SetEnergy(fail, Lower(SQRTS), sqme)
#ifdef FORCE_ONSHELL
	  if( abs(mass_in - Lower(SQRTS)) .lt. 1D-9 ) then
	    if( mass_in .lt. mass_out ) INFO "Below threshold"
	    return
	  endif
	  Lower(SQRTS) = mass_in
	  Upper(SQRTS) = mass_in
#else
	  if( fail .ne. 1 ) return
	  Lower(SQRTS) = minthreshold + .01D0
c	  INFO "Adjusting sqrtS to above-threshold value ", Lower(SQRTS)
#endif
	enddo

	Warning("Calculation of threshold did not converge")
	end


************************************************************************

	subroutine Luminosity(fac, sqme)
	implicit none
	RealType fac
	external sqme

	integer fail

#include "xsection.h"

	call SetEnergy(fail, Var(SQRTS), sqme)

	if( fail .ne. 0 ) fac = 0
	end

