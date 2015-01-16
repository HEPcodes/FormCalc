* 2to2.h
* common blocks for 2to2.F and num.F
* this file is part of FormCalc
* last modified 10 Jan 03 th

#include "looptools.h"
#include "model.h"

	double complex vec(0:3, -4*LEGS:4*LEGS)
	double precision mass2(LEGS)
	integer bpol(LEGS), epol(LEGS)
	double precision avgfac, sqrtS, flux, Pout
	integer cpus
	logical reset
	character*200 outfile

	common /global/
     &    vec, mass2, bpol, epol,
     &    avgfac, sqrtS, flux, Pout,
     &    cpus,
     &    reset,
     &    outfile

