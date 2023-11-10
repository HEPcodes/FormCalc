* 2to2.F
* the kinematics for a 2 -> 2 process
* this file is part of FormCalc
* last modified 30 Aug 22 th


#define LEGS 4
#define LEGS_IN 2
#define LEGS_OUT 2

#define MASS_IN mass(1) + mass(2)
#define MASS_OUT mass(3) + mass(4)

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

#include "inline.h"

	p = ThreeMom(sqrtS, mass(1), mass(2))

	i = FromPerm(crossing, 1)
	call VecSet(i, mass(i), p, 0D0,0D0,1D0)
	i = FromPerm(crossing, 2)
	call VecSet(i, mass(i), p, 0D0,0D0,-1D0)

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

	i = FromPerm(crossing, 4)
	call VecSet(i, mass(i), p, ex,ey,ez)
	end

