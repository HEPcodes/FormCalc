#if 0
	const.h
	model-independent constants
	this file is part of FormCalc
	last modified 23 Apr 13 th
#endif


#define Pi 3.1415926535897932384626433832795029
#define sqrt2 1.4142135623730950488016887242096981
#define cI I

#if NOUNDERSCORE
#define renorm_ renorm
#endif

struct renorm_ {
  RealType Divergence, mudim, lambda, muscale;
  integer epscoeff, Finite;
} renorm_;

#define Divergence renorm_.Divergence
#define mudim renorm_.mudim
#define lambda renorm_.lambda
#define muscale renorm_.muscale
#define epscoeff renorm_.epscoeff
#define Finite renorm_.Finite

