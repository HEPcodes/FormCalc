* config.h
* global declarations for the Diag routines
* this file is part of Diag
* last modified 13 Dec 06 th


* The maximum dimension of a matrix, needed for allocating internal
* memory, i.e. the routines handle at most MAXDIM-by-MAXDIM matrices.

#define MAXDIM 16


* A matrix is considered diagonal if the sum of the absolute values
* of the off-diagonal elements is less than EPS.

#define EPS 5D-16


* The transposed versions are needed for C, which has row-major
* matrix access.

#ifdef TRANSPOSE

#define Element(A,i,j) A(j,i)
#define HEigensystem HEigensystemT
#define SEigensystem SEigensystemT
#define CEigensystem CEigensystemT
#define TakagiFactor TakagiFactorT
#define SVD SVDT

#else

#define Element(A,i,j) A(i,j)

#endif

