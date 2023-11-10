* 1to2.F
* the kinematics for a 1 -> 2 process
* this file is part of FormCalc
* last modified 30 Aug 22 th


#define LEGS 3
#define LEGS_IN 1
#define LEGS_OUT 2

#define MASS_IN mass(1)
#define MASS_OUT mass(2) + mass(3)

* this defines the order of integrations:
#define XCOSTH1 1
#define MAXVAR XCOSTH1

#include LUMI


************************************************************************

	subroutine KinIni(fail)
	implicit none
	integer fail

#include "xsection.h"

	Lower(XCOSTH1) = 0
	Upper(XCOSTH1) = 1
	if( Lower(SQRTS) > Upper(SQRTS) ) then
	  Step(XCOSTH1) = XCOSTH1STEP
	else
	  Step(XCOSTH1) = 0
	endif

	fail = 0
	end


************************************************************************

	subroutine Kinematics(fac)
	implicit none
	RealType fac

#include "xsection.h"

	RealType p, ex,ey,ez, msum, minv
	integer i, key

#ifdef FORCE_ONSHELL
	p = 0
#else
	integer warn
	save warn

	p = dim(sqrtS, mass(1))
	if( warn == 0 .and. p > eps_sqrtS ) then
	  warn = 1
	  INFO "Computing off-shell decay"
	  INFO "Warning: Results correct only if decaying particle"
	  INFO "was chosen off-shell in CalcFeynAmp"
	endif
	p = sqrt(p*(sqrtS + mass(1)))
#endif

	i = FromPerm(crossing, 1)
	call VecSet(i, mass(i), p, 0D0,0D0,1D0)

	p = 0
	ex = 0
	ey = 0
	ez = 1
	minv = sqrtS
	msum = MASS_OUT

	key = 0
#ifdef COSTH1CUT
	CutMin(XCOSTH1) = -(1 - COSTH1CUT)
	CutMax(XCOSTH1) = +(1 - COSTH1CUT)
	key = key + Cut(CUT_COSTH, CUT_MIN + CUT_MAX)
#endif

	i = FromPerm(crossing, 3)
	call Split(i, mass(i), p, ex,ey,ez, minv, msum,
     &    fac, key, Var(FIXED), Var(XCOSTH1), Var(TRIVIAL))

	i = FromPerm(crossing, 2)
	call VecSet(i, mass(i), p, ex,ey,ez)
	end

