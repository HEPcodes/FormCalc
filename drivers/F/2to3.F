* 2to3.F
* the kinematics for a 2 -> 3 process
* this file is part of FormCalc
* last modified 18 Mar 19 th


#define LEGS 5
#define LEGS_IN 2
#define LEGS_OUT 3

#define MASS_IN mass(1) + mass(2)
#define MASS_OUT mass(3) + mass(4) + mass(5)

* this defines the order of integrations:
#define XMREM1 1
#define XCOSTH1 2
#define XCOSTH2 3
#define XPHI2 4
#define MAXVAR XPHI2

#include LUMI


************************************************************************

	subroutine KinIni(fail)
	implicit none
	integer fail

#include "xsection.h"

	Lower(XMREM1) = 0
	Upper(XMREM1) = 1
	if( Lower(SQRTS) > Upper(SQRTS) ) then
	  Lower(XMREM1) = XMREM1STEP
	  Step(XMREM1) = XMREM1STEP
	else
	  Step(XMREM1) = 0
	endif

	Lower(XCOSTH1) = 0
	Upper(XCOSTH1) = 1
	Step(XCOSTH1) = 0

	Lower(XCOSTH2) = 0
	Upper(XCOSTH2) = 1
	Step(XCOSTH2) = 0

	Lower(XPHI2) = 0
	Upper(XPHI2) = 1
	Step(XPHI2) = 0

	fail = 0
	end


************************************************************************

	subroutine Kinematics(fac)
	implicit none
	RealType fac

#include "xsection.h"

	RealType p, ex,ey,ez, msum, minv
	integer key

#include "inline.h"

	p = ThreeMom(sqrtS, mass(1), mass(2))

	call VecSet(1, mass(1), p, 0D0,0D0,1D0)
	call VecSet(2, mass(2), p, 0D0,0D0,-1D0)

	p = 0
	ex = 0
	ey = 0
	ez = 1
	minv = sqrtS
	msum = MASS_OUT

	key = 0
#ifdef E1MIN
	CutMin(XMREM1) = E1MIN
	key = key + Cut(CUT_MREM_E, CUT_MIN)
#endif
#ifdef COSTH1CUT
	CutMin(XCOSTH1) = -(1 - COSTH1CUT)
	CutMax(XCOSTH1) = +(1 - COSTH1CUT)
	key = key + Cut(CUT_COSTH, CUT_MIN + CUT_MAX)
#endif

	call Split(5, mass(5), p, ex,ey,ez, minv, msum,
     &    fac, key, Var(XMREM1), Var(XCOSTH1), Var(TRIVIAL))
	Show(XMREM1) = momspec(SPEC_E, 5)
	if( fac == 0 ) return

	call Split(4, mass(4), p, ex,ey,ez, minv, msum,
     &    fac, 0, Var(FIXED), Var(XCOSTH2), Var(XPHI2))

	call VecSet(3, mass(3), p, ex,ey,ez)
	end

