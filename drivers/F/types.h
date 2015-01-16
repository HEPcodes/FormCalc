* types.h
* real-based type declarations
* this file is part of FormCalc
* last modified 14 May 13 th


#ifndef RealType
#define RealType double precision
#define ComplexType double complex
#define Re DBLE
#define Im DIMAG
#define Conjugate DCONJG
#define ToComplex DCMPLX
#define Sq(c) Re((c)*Conjugate(c))
#endif

#ifndef marker
#define marker integer
#endif

