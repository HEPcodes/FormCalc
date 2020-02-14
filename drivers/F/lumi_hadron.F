* lumi_hadron.F
* parton luminosity for a hadron collider
* this file is part of FormCalc
* last modified 18 Mar 19 th

* This code links to the LHAPDF library.
* The following variables are used:
*
* - parton1, parton2 are the PDG parton codes,
* - fscale determines the factorization scale,
* - PDFSET specifies the file name of the PDF-set to be used.


#define PARTON_T -1
#define PARTON_X -2
#define MINVAR PARTON_X

#define HAVE_ALPHASPDF

#define HAVE_BOOST

************************************************************************

	subroutine LumiDefaults(argc, argv)
	implicit none
	integer argc
	character*128 argv(*)

#include "xsection.h"

	integer partonsign
	common /hadron/ partonsign

	if( Upper(SQRTS) > 0 .and.
     &      Lower(SQRTS) /= Upper(SQRTS) )
     &    stop "Cannot handle a range of beam energies"

	if( Lower(SQRTS) == 2000 ) then
	  partonsign = -1	! FNAL p-pbar
	else if( Lower(SQRTS) >= 7000 .and.
     &           Lower(SQRTS) <= 14000 ) then
	  partonsign = 1	! LHC p-p
	else
	  stop "Unknown collider"
	endif

	call InitPDFsetByName(PDFSET)
	call InitPDF(PDFMEM)
	end


************************************************************************

	subroutine LumiIni(fail, sqme)
	implicit none
	integer fail
	external sqme

#include "xsection.h"

	if( parton1 == -999 .or. parton2 == -999 )
     &    stop "parton1,2 not set"

#ifdef M_INV
        call SetEnergy(fail, Re(M_INV), sqme)
        if( fail /= 0 ) return

        Lower(PARTON_T) = M_INV/Lower(SQRTS)
        Upper(PARTON_T) = Lower(PARTON_T)
        Step(PARTON_T) = 1
#else
	call SetEnergy(fail, Lower(SQRTS), sqme)
	if( fail /= 0 ) return

	Lower(PARTON_T) = threshold/sqrtS
	Upper(PARTON_T) = 1
	Step(PARTON_T) = 0
#endif

	Lower(PARTON_X) = 0
	Upper(PARTON_X) = 1
	Step(PARTON_X) = 0
	end


************************************************************************

	subroutine Luminosity(fac, sqme)
	implicit none
	RealType fac
	external sqme

#include "xsection.h"

	RealType xPDF(-6:6), tPDF(-6:6)
	RealType tau, x, xrange
	integer fail

	integer partonsign
	common /hadron/ partonsign

	RealType gam
        common /boost/ gam

	logical tooSmall

#include "inline.h"

c	tooSmall(x) = x < eps_sqrtS .or. x > 1 - eps_sqrtS
	tooSmall(x) = abs(1/2D0 - x) > 1/2D0 - eps_sqrtS

#ifndef M_INV
	call SetEnergy(fail, Var(PARTON_T)*Var(SQRTS), sqme)
	if( fail /= 0 ) then
	  fac = 0
	  return
	endif
#endif

	tau = Var(PARTON_T)**2
	xrange = (1 - Var(PARTON_T))*(1 + Var(PARTON_T))
	x = tau + xrange*Var(PARTON_X)

	if( tooSmall(x) .or. tooSmall(tau/x) ) then
	  fac = 0
	  return
	endif

	gam = x/Var(PARTON_T)

	call evolvePDF(x, fscale, xPDF)
	call evolvePDF(tau/x, fscale, tPDF)

#ifdef M_INV
#define DEN M_INV
#else
#define DEN Var(PARTON_T)
#endif

	fac = fac*xrange/(DEN*x)*
     &    (2 - IndexDelta(parton1, parton2))*
     &    ( tPDF(partonsign*parton1)*xPDF(parton2) +
     &      xPDF(partonsign*parton2)*tPDF(parton1) )
	end

