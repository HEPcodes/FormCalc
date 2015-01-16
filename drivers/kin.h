* kin.h
* common blocks with kinematical variables for num.F
* this file is part of FormCalc
* last modified 14 Jun 01 th

#include "looptools.h"
#include "model.h"

	double complex vec(0:3, -16:16)
	double precision mass(4)
	integer bpol(4), epol(4)
	double precision avgfac, beta
	logical reset

	common /numFglobal/ vec, mass, bpol, epol,
     +    avgfac, beta, reset

