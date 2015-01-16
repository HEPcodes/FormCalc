* 2to2.h
* common blocks for 2to2.F and num.F
* this file is part of FormCalc
* last modified 24 Jul 01 th

#include "looptools.h"
#include "model.h"

	double complex vec(0:3, -4*4:4*4)
	double precision mass(4)
	integer bpol(4), epol(4)
	double precision avgfac, sqrtS, flux, Pout
	logical reset

	common /global/
     +    vec, mass, bpol, epol,
     +    avgfac, sqrtS, flux, Pout,
     +    reset

