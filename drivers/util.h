* util.h
* prototypes for the functions in util.a
* this file is part of FormCalc
* last modified 13 Apr 06 th


	double precision ThreeMom
	double precision SInvariant, TInvariant
	double complex Pair, Eps
	double complex SxS, SeS
	integer VxS, VeS, BxS, BeS

	external ThreeMom
	external SInvariant, TInvariant
	external Pair, Eps
	external SxS, SeS
	external VxS, VeS, BxS, BeS

#ifdef LEGS
	double complex vec(2,2, 8, 0:LEGS)
	common /vectors/ vec

	double precision momspec(8, LEGS)
	common /momenta/ momspec
#endif

#ifndef Delta
#define Delta(i,j) ibits(ibits(i-(j),0,9)-1,9,1)

#define k(i) (8*i+1)
#define s(i) (8*i+3)
#define e(i) (8*i+3+Hel(i))
#define ec(i) (8*i+3-Hel(i))
#define Spinor(i,s,om) (s*2*Hel(i)+16*i+om+5)
#define DottedSpinor(i,s,om) (s*2*Hel(i)+16*i+om+7)

#define SPEC_M 1
#define SPEC_E 2
#define SPEC_K 3
#define SPEC_ET 4
#define SPEC_KT 5
#define SPEC_RAP 6
#define SPEC_PRAP 7
#define SPEC_DELTAK 8

#define Cut(c, m) (m)*(c)

#define CUT_MIN 1
#define CUT_MAX 2

#define CUT_COSTH 4
#define CUT_COSTHCMS 16
#define CUT_COSTH_E 64
#define CUT_COSTH_K 65
#define CUT_MREM 256
#define CUT_MREM_E 1024
#define CUT_MREM_K 1025
#define CUT_MREM_ET 4096
#define CUT_MREM_KT 4097
#define CUT_MREM_RAP 16384
#define CUT_MREM_PRAP 16385
#endif

