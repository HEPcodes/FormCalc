* 1to3.F
* the kinematics for a 1 -> 3 process
* this file is part of FormCalc
* last modified 15 May 15 th


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
	if( Lower(SQRTS) .gt. Upper(SQRTS) ) then
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

	subroutine InitialState(m1)
	implicit none
	RealType m1

#include "xsection.h"

	RealType p

#ifdef FORCE_ONSHELL
	p = 0
#else
	p = dim(sqrtS, m1)
	if( p .gt. 1D-6 ) then
	  INFO "Computing off-shell decay"
	  INFO "Warning: Results correct only if decaying particle"
	  INFO "was chosen off-shell in CalcFeynAmp"
	endif
	p = sqrt(p*(sqrtS + m1))
#endif

	call VecSet(1, m1, p, 0D0, 0D0, 1D0)
	end


************************************************************************

	subroutine FinalState(fac, m2, m3, m4)
	implicit none
	RealType fac, m2, m3, m4

#include "xsection.h"

	RealType p, ex,ey,ez, msum, minv
	integer key

	p = 0
	ex = 0
	ey = 0
	ez = 1
	minv = sqrtS
	msum = m2 + m3 + m4

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

	call Split(4, m4, p, ex,ey,ez, minv, msum,
     &    fac, key, Var(XMREM1), Var(XCOSTH1), Var(TRIVIAL))
	Show(XMREM1) = momspec(SPEC_E, 4)
	if( fac .eq. 0 ) return

	call Split(3, m3, p, ex,ey,ez, minv, msum,
     &    fac, 0, Var(FIXED), Var(XCOSTH2), Var(XPHI2))

	call VecSet(2, m2, p, ex,ey,ez)
	end

