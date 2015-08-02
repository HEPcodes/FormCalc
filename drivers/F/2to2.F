* 2to2.F
* the kinematics for a 2 -> 2 process
* this file is part of FormCalc
* last modified 6 May 13 th


* some definitions for xsection.F
#define MASS_IN MASS1 + MASS2
#define MASS_OUT MASS3 + MASS4
#define TYPES TYPE1, TYPE2, TYPE3, TYPE4

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
	if( Lower(SQRTS) .gt. Upper(SQRTS) ) then
	  Step(XCOSTH1) = XCOSTH1STEP
	else
	  Step(XCOSTH1) = 0
	endif

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

	subroutine FinalState(fac, m3, m4)
	implicit none
	RealType fac, m3, m4

#include "xsection.h"

	RealType p, ex,ey,ez, msum, minv
	integer key

	p = 0
	ex = 0
	ey = 0
	ez = 1
	minv = sqrtS
	msum = m3 + m4

	key = 0
#ifdef COSTH1CUT
	CutMin(XCOSTH1) = -(1 - COSTH1CUT)
	CutMax(XCOSTH1) = +(1 - COSTH1CUT)
	key = key + Cut(CUT_COSTH, CUT_MIN + CUT_MAX)
#endif

	call Split(3, m3, p, ex,ey,ez, minv, msum,
     &    fac, key, Var(FIXED), Var(XCOSTH1), Var(TRIVIAL))

	call VecSet(4, m4, p, ex,ey,ez)
	end

