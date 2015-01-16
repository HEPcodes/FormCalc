* opp.h
* declarations for the OPP routines
* this file is part of FormCalc
* last modified 8 Oct 12 th


#ifndef Bcut
* for now:
#define Mbb 1
#define Mcc 1
#define Mdd 1
#define Mee 1
#define Mff 1

#ifdef SAMURAI
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
#error Neither SAMURAI nor CUTTOOLS defined
#endif

#endif

	ComplexType Bcut, Ccut, Dcut, Ecut, Fcut
	external Bcut, Ccut, Dcut, Ecut, Fcut

