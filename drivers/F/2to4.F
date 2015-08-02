* 2to4.F
* the kinematics for a 2 -> 4 process
* this file is part of FormCalc
* last modified 15 May 15 th


* this defines the order of integrations:
#define XMREM1 1
#define XCOSTH1 2
#define XMREM2 3
#define XCOSTH2 4
#define XPHI2 5
#define XCOSTH3 6
#define XPHI3 7
#define MAXVAR XPHI3

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

	Lower(XMREM2) = 0
	Upper(XMREM2) = 1
	Step(XMREM2) = 0

	Lower(XCOSTH1) = 0
	Upper(XCOSTH1) = 1
	Step(XCOSTH1) = 0

	Lower(XCOSTH2) = 0
	Upper(XCOSTH2) = 1
	Step(XCOSTH2) = 0

	Lower(XPHI2) = 0
	Upper(XPHI2) = 1
	Step(XPHI2) = 0

	Lower(XCOSTH3) = 0
	Upper(XCOSTH3) = 1
	Step(XCOSTH3) = 0

	Lower(XPHI3) = 0
	Upper(XPHI3) = 1
	Step(XPHI3) = 0

	fail = 0
	end


************************************************************************

	subroutine InitialState(m1, m2)
	implicit none
	RealType m1, m2

#include "xsection.h"

	RealType p

#include "inline.h"

	p = ThreeMom(sqrtS, m1, m2)

	call VecSet(1, m1, p, 0D0, 0D0, 1D0)
	call VecSet(2, m2, p, 0D0, 0D0, -1D0)
	end


************************************************************************

	subroutine FinalState(fac, m3, m4, m5, m6)
	implicit none
	RealType fac, m3, m4, m5, m6

#include "xsection.h"

	RealType p, ex,ey,ez, msum, minv
	integer key

	p = 0
	ex = 0
	ey = 0
	ez = 1
	minv = sqrtS
	msum = m3 + m4 + m5 + m6

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

	call Split(6, m6, p, ex,ey,ez, minv, msum,
     &    fac, key, Var(XMREM1), Var(XCOSTH1), Var(TRIVIAL))
	Show(XMREM1) = momspec(SPEC_E, 6)
	if( fac .eq. 0 ) return

	call Split(5, m5, p, ex,ey,ez, minv, msum,
     &    fac, 0, Var(XMREM2), Var(XCOSTH2), Var(XPHI2))
	Show(XMREM2) = momspec(SPEC_E, 5)
	if( fac .eq. 0 ) return

	call Split(4, m4, p, ex,ey,ez, minv, msum,
     &    fac, 0, Var(FIXED), Var(XCOSTH3), Var(XPHI3))

	call VecSet(3, m3, p, ex,ey,ez)
	end

