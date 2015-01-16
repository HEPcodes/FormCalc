* types.h
* real-based type declarations
* this file is part of FormCalc
* last modified 18 Mar 13 th


#ifndef RealType
#define RealType double precision
#define ComplexType double complex
#define HelType ComplexType
#define Re DBLE
#define Im DIMAG
#define Conjugate DCONJG
#define ToComplex DCMPLX
#define Sq(c) Re((c)*Conjugate(c))
#endif

#ifndef marker
#define marker double precision
#endif

