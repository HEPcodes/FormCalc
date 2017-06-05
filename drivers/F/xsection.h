* xsection.h
* common blocks for xsection.F
* this file is part of FormCalc
* last modified 20 Feb 17 th


#include "decl.h"

#ifndef SQRTS
#define SQRTS 0
#define FIXED MAXVAR+1
#define TRIVIAL MAXVAR+2
#endif

	RealType eps_sqrtS
	parameter (eps_sqrtS = 1D-9)

	integer nvars
	parameter (nvars = MAXVAR - (MINVAR) + 1)

	RealType var(8,MINVAR:TRIVIAL)
	RealType mass(LEGS)
	RealType charge(LEGS,NPID)
	RealType colorcharge(LEGS,NPID)
	RealType sqrtS, mass_in, mass_out, fscale
	RealType threshold, avgfac
	integer*8 helmask
	integer type, pid, parton1, parton2, serial
	integer sqrtSinvalid, flags

	common /xsection/ var,
     &    mass, charge, colorcharge,
     &    sqrtS, mass_in, mass_out, fscale,
     &    threshold, avgfac,
     &    helmask,
     &    type, pid, parton1, parton2, serial,
     &    sqrtSinvalid, flags

