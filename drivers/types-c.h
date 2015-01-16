/*
  types-c.h
  real-based type declarations
  this file is part of FormCalc
  last modified 22 May 11 th
*/


#ifndef TYPES_H
#define TYPES_H

#ifdef QUAD
#define RealType long double
#define Re creall
#define Im cimagl
#define Conjugate conjl
#else
#define RealType double
#define Re creal
#define Im cimag
#define Conjugate conj
#endif

typedef int integer;
typedef RealType complex ComplexType;

#define ToComplex(r,i) ((r) + I*(i))
#define Sq(c) Re((c)*Conjugate(c))

#endif

