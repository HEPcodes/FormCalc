* lumi_hadron.F
* parton luminosity for a hadron collider
* this file is part of FormCalc
* last modified 29 Mar 16 th

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

	if( Upper(SQRTS) .gt. 0 .and.
     &      Lower(SQRTS) .ne. Upper(SQRTS) )
     &    stop "Cannot handle a range of beam energies"

	if( Lower(SQRTS) .eq. 2000 ) then
	  partonsign = -1
	else if( Lower(SQRTS) .ge. 7000 .and.
     &           Lower(SQRTS) .le. 14000 ) then
	  partonsign = 1
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

	if( parton1 .eq. -999 .or. parton2 .eq. -999 )
     &    stop "parton1,2 not set"

#ifdef M_INV
        call SetEnergy(fail, Re(M_INV), sqme)
        if( fail .ne. 0 ) return

        Lower(PARTON_T) = M_INV/Lower(SQRTS)
        Upper(PARTON_T) = Lower(PARTON_T)
        Step(PARTON_T) = 1
#else
	call SetEnergy(fail, Lower(SQRTS), sqme)
	if( fail .ne. 0 ) return

	Lower(PARTON_T) = minthreshold/sqrtS
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

c	tooSmall(x) = x .lt. SQRTS_EPS .or. x .gt. 1 - SQRTS_EPS
	tooSmall(x) = abs(.5D0 - x) .gt. .5D0 - SQRTS_EPS

#ifndef M_INV
	call SetEnergy(fail, Var(PARTON_T)*Var(SQRTS), sqme)
	if( fail .ne. 0 ) then
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

