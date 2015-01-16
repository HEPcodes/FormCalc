* 2to3.h
* common blocks for 2to3.F and num.F
* this file is part of FormCalc
* last modified 18 Sep 01 th

#include "looptools.h"
#include "model.h"

	double complex vec(0:3, -4*5:4*5)
	double precision mass(5)
	integer bpol(5), epol(5)
	double precision upper(NDIM), lower(NDIM), var(NDIM)
	double precision avgfac, sqrtS, flux
	logical reset

	common /global/
     +    vec, mass, bpol, epol,
     +    upper, lower, var,
     +    avgfac, sqrtS, flux,
     +    reset

