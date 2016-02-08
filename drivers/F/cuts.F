* cuts.F
* routines to compute cuts
* this file is part of FormCalc
* Christian Schappacher, 2 Dec 13
* last modified 18 Jan 16 cs


************************************************************************
** Cuts calculates standard cuts via some simple veto functions

	subroutine Cuts(fac)
        implicit none
	RealType fac

#include "xsection.h"

	RealType acosmin, cosdeltaalpha
	RealType deltay, deltaeta, deltaphi
	RealType sepR, seprho, invmass
	RealType x
	integer i, j

* Rapidity cuts:

#ifdef HAVE_BOOST
	RealType gam, y, omx2, eta
	common /boost/ gam
	y(i) = momspec(SPEC_RAP,i) - log(gam)
	omx2(x) = (1 - x)*(1 + x)
	eta(i) = atanh(tanh(y(i))/
     &    sqrt(omx2( 2*gam*momspec(SPEC_M,i)/
     &      (momspec(SPEC_E,i)*(1 + gam**2 +
     &        omx2(gam)*tanh(momspec(SPEC_RAP,i)))) )))
#else
#define y(i) momspec(SPEC_RAP,i)
#define eta(i) momspec(SPEC_PRAP,i)
#endif

#define CUT_y(i) abs(y(i))
#define CUT_eta(i) abs(eta(i))

* Energy and momentum cuts:

#define CUT_E(i) momspec(SPEC_E,i)
#define CUT_k(i) momspec(SPEC_K,i)
#define CUT_ET(i) momspec(SPEC_ET,i)
#define CUT_kT(i) momspec(SPEC_KT,i)

* Collinearity cuts:

	acosmin(x) = pi/2 - abs(acos(x) - pi/2)

	cosdeltaalpha(i,j) =
     &    momspec(SPEC_EX,i)*momspec(SPEC_EX,j) +
     &    momspec(SPEC_EY,i)*momspec(SPEC_EY,j) +
     &    momspec(SPEC_EZ,i)*momspec(SPEC_EZ,j)

#define CUT_deltatheta(i) acosmin(momspec(SPEC_EZ,i))
#define CUT_cosdeltatheta(i) abs(momspec(SPEC_EZ,i))
#define CUT_deltaalpha(i,j) acos(cosdeltaalpha(i,j))
#define CUT_cosdeltaalpha(i,j) abs(cosdeltaalpha(i,j))

* Rapidity-gap cuts:

	deltay(i,j) = y(j) - y(i)
	deltaeta(i,j) = eta(j) - eta(i)

#define CUT_deltay(i,j) abs(deltay(i,j))
#define CUT_deltaeta(i,j) abs(deltaeta(i,j))

* Separation cuts:

	deltaphi(i,j) = pi -
     &    abs(abs(momspec(SPEC_PHI,j) - momspec(SPEC_PHI,i)) - pi)

	sepR(i,j) = sqrt(deltay(i,j)**2 + deltaphi(i,j)**2)
	seprho(i,j) = sqrt(deltaeta(i,j)**2 + deltaphi(i,j)**2)

#define CUT_R(i,j) sepR(i,j)
#define CUT_rho(i,j) seprho(i,j)

* Opposite-hemisphere cuts:

#define CUT_yprod(i,j) y(i)*y(j)
#define CUT_etaprod(i,j) eta(i)*eta(j)

* Dijet invariant-mass cuts:

	invmass(i,j) = sqrt(
     &    (momspec(SPEC_E,i) + momspec(SPEC_E,j))**2 -
     &    (momspec(SPEC_K,i)*momspec(SPEC_EX,i) +
     &     momspec(SPEC_K,j)*momspec(SPEC_EX,j))**2 -
     &    (momspec(SPEC_K,i)*momspec(SPEC_EY,i) +
     &     momspec(SPEC_K,j)*momspec(SPEC_EY,j))**2 -
     &    (momspec(SPEC_K,i)*momspec(SPEC_EZ,i) +
     &     momspec(SPEC_K,j)*momspec(SPEC_EZ,j))**2 )

#define CUT_invmass(i,j) invmass(i,j)

#ifndef CUT1
#define CUT1
#endif
#ifndef CUT2
#define CUT2
#endif
#ifndef CUT3
#define CUT3
#endif
#ifndef CUT4
#define CUT4
#endif
#ifndef CUT5
#define CUT5
#endif
#ifndef CUT6
#define CUT6
#endif
#ifndef CUT7
#define CUT7
#endif
#ifndef CUT8
#define CUT8
#endif
#ifndef CUT9
#define CUT9
#endif
#ifndef CUT10
#define CUT10
#endif
#ifndef CUT11
#define CUT11
#endif
#ifndef CUT12
#define CUT12
#endif
#ifndef CUT13
#define CUT13
#endif
#ifndef CUT14
#define CUT14
#endif
#ifndef CUT15
#define CUT15
#endif
#ifndef CUT16
#define CUT16
#endif
#ifndef CUT17
#define CUT17
#endif
#ifndef CUT18
#define CUT18
#endif
#ifndef CUT19
#define CUT19
#endif
#ifndef CUT20
#define CUT20
#endif

	if(
     &    CUT1
     &    CUT2
     &    CUT3
     &    CUT4
     &    CUT5
     &    CUT6
     &    CUT7
     &    CUT8
     &    CUT9
     &    CUT10
     &    CUT11
     &    CUT12
     &    CUT13
     &    CUT14
     &    CUT15
     &    CUT16
     &    CUT17
     &    CUT18
     &    CUT19
     &    CUT20
     &  ) return

	fac = 0
	end

