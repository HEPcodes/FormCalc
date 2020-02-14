#if 0
	THDM.ch
	C declarations for THDM.F
	this file is part of FormCalc
	last modified 6 Oct 19 th
#endif


#include "SM.ch"

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

