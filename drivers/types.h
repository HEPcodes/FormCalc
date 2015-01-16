* types.h
* real-based type declarations
* this file is part of FormCalc
* last modified 9 Aug 11 th


#ifndef TYPES_H
#define TYPES_H

#ifdef QUAD
#define Real real*16
#define Complex complex*32
#define Re QEXT
#define Im QIMAG
#define Conjugate QCONJG
#define Cmplx QCMPLX
#else
#define Real double precision
#define Complex double complex
#define Re DBLE
#define Im DIMAG
#define Conjugate DCONJG
#define Cmplx DCMPLX
#endif

#define Sq(c) Re((c)*Conjugate(c))

#endif

