* num.h
* headers for the computation of the numerators
* this file is part of FormCalc
* last modified 20 Dec 10 th


#ifndef NUM_H
#define NUM_H

#define q1 2
#define Q1(x,y) vec(x,y, q1, 0)

#ifdef SAMURAI
#define NumeratorFunction(f) function f(ncut, q1in, MuTildeSq)
#define Result(f) f
#else
#define NumeratorFunction(f) subroutine f(res, q1in)
#define Result(f) res
#endif

#else

#ifdef SAMURAI
	integer ncut
	double complex q1in(4)
	double precision MuTildeSq
	Q1(1,1) = q1in(4) + q1in(3)
	Q1(2,2) = q1in(4) - q1in(3)
	Q1(2,1) = q1in(1) + cI*q1in(2)
	Q1(1,2) = q1in(1) - cI*q1in(2)
#else
	double complex q1in(0:3)
	Q1(1,1) = q1in(0) + q1in(3)
	Q1(2,2) = q1in(0) - q1in(3)
	Q1(2,1) = q1in(1) + cI*q1in(2)
	Q1(1,2) = q1in(1) - cI*q1in(2)
#endif

#endif

