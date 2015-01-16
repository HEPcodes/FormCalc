* opp.h
* declarations for the OPP routines
* this file is part of FormCalc
* last modified 6 Dec 10 th


#ifndef Acut
#ifdef SAMURAI
#define Acut SAAcut
#define Bcut SABcut
#define Ccut SACcut
#define Dcut SADcut
#define Ecut SAEcut
#define Fcut SAFcut
#else
#define Acut CTAcut
#define Bcut CTBcut
#define Ccut CTCcut
#define Dcut CTDcut
#define Ecut CTEcut
#define Fcut CTFcut
#endif
#endif

	double complex Acut, Bcut, Ccut, Dcut, Ecut, Fcut
	external Acut, Bcut, Ccut, Dcut, Ecut, Fcut

