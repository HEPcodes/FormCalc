* 2to2.F
* the kinematics for a 2 -> 2 process
* this file is part of FormCalc
* last modified 3 Jun 15 th


#define LEGS 4
#define LEGS_IN 2
#define LEGS_OUT 2

#define MASS_IN(pid) mass(1,pid) + mass(2,pid)
#define MASS_OUT(pid) mass(3,pid) + mass(4,pid)

* this defines the order of integrations:
#define XCOSTH1 1
#define MAXVAR XCOSTH1

#include LUMI


************************************************************************

	subroutine SetArray(ar, a1, a2, a3, a4)
	implicit none
	RealType ar(*), a1, a2, a3, a4

	ar(1) = a1
	ar(2) = a2
	ar(3) = a3
	ar(4) = a4
	end


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

	subroutine Kinematics(fac)
	implicit none
	RealType fac

#include "xsection.h"

	RealType p, ex,ey,ez, msum, minv
	integer key

#include "inline.h"

	p = ThreeMom(sqrtS, mass(1,pid), mass(2,pid))

	call VecSet(1, mass(1,pid), p, 0D0, 0D0, 1D0)

	call VecSet(2, mass(2,pid), p, 0D0, 0D0, -1D0)

	p = 0
	ex = 0
	ey = 0
	ez = 1
	minv = sqrtS
	msum = MASS_OUT(pid)

	key = 0
#ifdef COSTH1CUT
	CutMin(XCOSTH1) = -(1 - COSTH1CUT)
	CutMax(XCOSTH1) = +(1 - COSTH1CUT)
	key = key + Cut(CUT_COSTH, CUT_MIN + CUT_MAX)
#endif

	call Split(3, mass(3,pid), p, ex,ey,ez, minv, msum,
     &    fac, key, Var(FIXED), Var(XCOSTH1), Var(TRIVIAL))

	call VecSet(4, mass(4,pid), p, ex,ey,ez)
	end

