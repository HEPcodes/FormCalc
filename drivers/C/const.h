#if 0
	const.h
	model-independent constants
	this file is part of FormCalc
	last modified 28 Jan 14 th
#endif


#define Pi 3.1415926535897932384626433832795029
#define sqrt2 1.4142135623730950488016887242096981
#define cI I

#if NOUNDERSCORE
#define renorm_ renorm
#endif

struct renorm_ {
  RealType Divergence, mudim, lambda, muscale;
  integer epsi;
} renorm_;

#define Divergence renorm_.Divergence
#define mudim renorm_.mudim
#define lambda renorm_.lambda
#define muscale renorm_.muscale
#define epsi renorm_.epsi

#define Finite (1>>epsi)
#define PVC(i) i+epsi

