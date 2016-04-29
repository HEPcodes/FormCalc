* xsection.h
* common blocks for xsection.F
* this file is part of FormCalc
* last modified 18 Mar 16 th


#include "decl.h"

#ifndef SQRTS
#define SQRTS 0
#define FIXED MAXVAR+1
#define TRIVIAL MAXVAR+2
#endif

	integer nvars
	parameter (nvars = MAXVAR - (MINVAR) + 1)

	RealType var(8,MINVAR:TRIVIAL)
	RealType mass(LEGS,NPID)
	RealType charge(LEGS,NPID)
	RealType colorcharge(LEGS,NPID)
	RealType avgfac(NPID), threshold(NPID), minthreshold
	RealType sqrtS, mass_in, mass_out, fscale
	integer*8 helmask, helicities
	integer type, pid, parton1, parton2, serial
	integer sqrtSinvalid, flags

	common /xsection/ var,
     &    mass, charge, colorcharge,
     &    avgfac, threshold, minthreshold,
     &    sqrtS, mass_in, mass_out, fscale,
     &    helmask, helicities,
     &    type, pid, parton1, parton2, serial,
     &    sqrtSinvalid, flags

