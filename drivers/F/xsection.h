* xsection.h
* common blocks for xsection.F
* this file is part of FormCalc
* last modified 11 Jun 15 th


#include "decl.h"

#ifndef SQRTS
#define SQRTS 0
#define FIXED MAXVAR+1
#define TRIVIAL MAXVAR+2
#define Var(v) var(1,v)
#define Show(v) var(2,v)
#define Lower(v) var(3,v)
#define Upper(v) var(4,v)
#define Step(v) var(5,v)
#define CutMin(v) var(6,v)
#define CutMax(v) var(7,v)
#endif

	integer nvars
	parameter (nvars = MAXVAR - (MINVAR) + 1)

	RealType var(8,MINVAR:TRIVIAL)
	RealType mass(LEGS,NPID), charge(LEGS,NPID)
	RealType avgfac(NPID), threshold(NPID), minthreshold
	RealType sqrtS, mass_in, mass_out, fscale
	integer*8 helmask, hel
	integer type, pid, parton1, parton2
	integer sqrtSinvalid, flags

	common /xsection/ var,
     &    mass, charge,
     &    avgfac, threshold, minthreshold,
     &    sqrtS, mass_in, mass_out, fscale,
     &    helmask, hel,
     &    type, pid, parton1, parton2,
     &    sqrtSinvalid, flags

