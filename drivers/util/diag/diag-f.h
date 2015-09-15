* diag-f.h
* global declarations for the Diag routines
* this file is part of Diag
* last modified 3 Sep 15 th


#include "types.h"


#if UCOLS
#define UL(i,j) U(j,i)
#define VL(i,j) V(j,i)
#define WL(i,j) W(j,i)
#else
#define UL(i,j) U(i,j)
#define VL(i,j) V(i,j)
#define WL(i,j) W(i,j)
#endif


* The maximum dimension of a matrix, needed for allocating internal
* memory, i.e. the routines handle at most MAXDIM-by-MAXDIM matrices.

#define MAXDIM 16


* A matrix is considered diagonal if the sum of the squares
* of the off-diagonal elements is less than EPS.  SYM_EPS is
* half of EPS since only the upper triangle is counted for
* symmetric matrices.
* 52 bits is the mantissa length for IEEE double precision.

#define EPS 2D0**(-102)

#define SYM_EPS 2D0**(-103)

#define DBL_EPS 2D0**(-52)

