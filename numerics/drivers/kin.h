* kin.h
* common blocks with kinematical variables for num.F
* this file is part of FormCalc
* last modified 2 May 00 th

#include "looptools.h"
#include "model.h"

	double complex vec(4, -16:16)
	double precision mass(4)
	integer bpol(4), epol(4)
	double precision avgfac, beta, sinth
	logical reset

	common /numFglobal/ vec, mass, bpol, epol,
     +    avgfac, beta, sinth, reset

	double precision degree
	parameter (degree = pi/180D0)

