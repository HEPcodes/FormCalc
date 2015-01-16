/*
	vegas-f.c
		Fortran interface for Vegas
		this file is part of Vegas
		last modified 9 Feb 05 th
*/


#include "decl.h"

Extern char vegasstate_[MAXSTATESIZE];

#ifdef HAVE_UNDERSCORE
#define vegas vegas_
#endif


Extern void Vegas(ccount ndim, ccount ncomp, Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, ccount mineval, ccount maxeval,
  ccount nstart, ccount nincrease,
  count *pneval, int *pfail,
  real *integral, real *error, real *prob);


Extern void vegas(ccount *pndim, ccount *pncomp, Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, ccount *pmineval, ccount *pmaxeval,
  ccount *pnstart, ccount *pnincrease, 
  count *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  /* make sure the filename is null-terminated */
  if( *vegasstate_ ) {
    char *p;
    vegasstate_[sizeof(vegasstate_) - 1] = 0;
    if( (p = strchr(vegasstate_, ' ')) ) *p = 0;
  }

  Vegas(*pndim, *pncomp, integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pnstart, *pnincrease,
    pneval, pfail,
    integral, error, prob);
}

