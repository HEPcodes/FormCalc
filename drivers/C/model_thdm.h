#if 0
	model_thdm.h
	declarations for model_thdm.F
	this file is part of FormCalc
	last modified 9 Mar 13 th
#endif


#include "model_sm.h"

#if NOUNDERSCORE
#define thdmpara_ thdmpara
#endif

struct thdmpara_ {
  RealType Lambda5;
  RealType Yuk1, Yuk2, Yuk3;
} thdmpara_;

#define Lambda5 thdmpara_.Lambda5
#define Yuk1 thdmpara_.Yuk1
#define Yuk2 thdmpara_.Yuk2
#define Yuk3 thdmpara_.Yuk3

