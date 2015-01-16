* opp.h
* declarations for the OPP routines
* this file is part of FormCalc
* last modified 7 Jul 14 th


#if SIMD > 0
#error SIMD must be set to 0 in distrib.h for OPP
#endif


#ifndef Bcut
* for now:
#define Mbb 1
#define Mcc 1
#define Mdd 1
#define Mee 1
#define Mff 1

#ifdef NINJA
#define Bcut NJBcut
#define Ccut NJCcut
#define Dcut NJDcut
#define Ecut NJEcut
#define Fcut NJFcut
#define Bmas NJBmas
#define Cmas NJCmas
#define Dmas NJDmas
#define Emas NJEmas
#define Fmas NJFmas
#elif defined SAMURAI
#define Bcut SABcut
#define Ccut SACcut
#define Dcut SADcut
#define Ecut SAEcut
#define Fcut SAFcut
#define Bmas SABmas
#define Cmas SACmas
#define Dmas SADmas
#define Emas SAEmas
#define Fmas SAFmas
#elif defined CUTTOOLS
#define Bcut CTBcut
#define Ccut CTCcut
#define Dcut CTDcut
#define Ecut CTEcut
#define Fcut CTFcut
#define Bmas CTBmas
#define Cmas CTCmas
#define Dmas CTDmas
#define Emas CTEmas
#define Fmas CTFmas
#else
#error No OPP method (NINJA, SAMURAI, CUTTOOLS) defined in user.h
#endif

#endif

	ComplexType Bcut, Ccut, Dcut, Ecut, Fcut
	external Bcut, Ccut, Dcut, Ecut, Fcut

#ifdef NINJA
	external None
#endif

