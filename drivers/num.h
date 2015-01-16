* num.h
* the prototypes of the some functions from num.F needed by the
* FormCalc-generated code
* this file is part of FormCalc
* last modified 27 Nov 02 th

	double precision ThreeMom
	integer Delta
	double precision SInvariant, TInvariant
	double complex Pair, Eps
	double precision Li2

	external ThreeMom
	external Delta
	external SInvariant, TInvariant
	external Pair, Eps
	external Li2

#ifndef k
#define k(i) (4*i - 3)
#define s(i) (4*i - 2)
#define e(i) (4*i - 2 + Hel(i))
#define ec(i) -e(i)
#endif

