* num.h
* declarations for the computation of the numerators
* this file is part of FormCalc
* last modified 4 Mar 08 th


#include "util.h"

#define q1 2
#define Q1(x,y) vec(x,y, q1)

	double complex vec(2,2, 1)
	common /vectors/ vec

	Q1(1,1) = q1in(0) + q1in(3)
	Q1(2,2) = q1in(0) - q1in(3)
	Q1(2,1) = q1in(1) + cI*q1in(2)
	Q1(1,2) = q1in(1) - cI*q1in(2)

