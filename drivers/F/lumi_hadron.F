* lumi_hadron.F
* parton luminosity for a hadron collider
* this file is part of FormCalc
* last modified 11 Jun 15 th

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

	if( ior(parton1, parton2) .lt. 0 ) stop "parton1,2 not set"

	call SetEnergy(fail, Lower(SQRTS), sqme)
	if( fail .ne. 0 ) return

	Lower(PARTON_T) = minthreshold/sqrtS
	Upper(PARTON_T) = 1
	Step(PARTON_T) = 0

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

	call SetEnergy(fail, Var(PARTON_T)*Var(SQRTS), sqme)
	if( fail .ne. 0 ) then
	  fac = 0
	  return
	endif

	tau = Var(PARTON_T)**2
	xrange = (1 - Var(PARTON_T))*(1 + Var(PARTON_T))
	x = tau + xrange*Var(PARTON_X)

	if( x .lt. 1D-10 .or. x .gt. (1 - 1D-10) ) then
	  fac = 0
	  return
	endif

	call evolvePDF(x, fscale, xPDF)
	call evolvePDF(tau/x, fscale, tPDF)

	fac = fac*xrange/(Var(PARTON_T)*x)*
     &    (2 - IndexDelta(parton1, parton2))*
     &    ( tPDF(partonsign*parton1)*xPDF(parton2) +
     &      xPDF(partonsign*parton2)*tPDF(parton1) )
	end

