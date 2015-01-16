* 2to3.h
* common blocks for 2to3.F
* this file is part of FormCalc
* last modified 25 Apr 04 th


#include "model.h"

	double precision upper(MAXDIM), lower(MAXDIM), var(MAXDIM)
	double precision preflux, flux, sqrtS
	integer helicities
	logical reset
	common /var2to3/ upper, lower, var, preflux, flux, sqrtS,
     &    helicities, reset

