* softphoton.F
* routines for the treatment of soft-photon radiation
* this file is part of FormCalc
* last modified 19 Nov 13 th


************************************************************************
** SoftPhotonFactor adds up the contributions from soft-photon emission
** off the different legs.  The result then has to be multiplied with
** the Born cross-section.

	RealType function SoftPhotonFactor()
	implicit none

#include "xsection.h"

	RealType PhaseSpaceIntegral
	external PhaseSpaceIntegral

	integer i, j

	RealType loglambda
	common /infrared/ loglambda

	RealType charge(LEGS), r

#ifndef Charge1
#define Charge1 0
#endif

#ifndef Charge2
#define Charge2 0
#endif

#ifndef Charge3
#define Charge3 0
#endif

#ifndef Charge4
#define Charge4 0
#endif

#ifndef Charge5
#define Charge5 0
#endif

#ifndef Charge6
#define Charge6 0
#endif

* not a very nice solution, the following, but necessary
* because there may be no "/" in a constant in the data
* statement (such as in 2/3D0)
#if LEGS == 3
#define CHARGES ch1, ch2, ch3
	RealType CHARGES
	parameter (ch1 = Charge1)
	parameter (ch2 = Charge2)
	parameter (ch3 = Charge3)
#elif LEGS == 4
#define CHARGES ch1, ch2, ch3, ch4
	RealType CHARGES
	parameter (ch1 = Charge1)
	parameter (ch2 = Charge2)
	parameter (ch3 = Charge3)
	parameter (ch4 = Charge4)
#elif LEGS == 5
#define CHARGES ch1, ch2, ch3, ch4, ch5
	RealType CHARGES
	parameter (ch1 = Charge1)
	parameter (ch2 = Charge2)
	parameter (ch3 = Charge3)
	parameter (ch4 = Charge4)
	parameter (ch5 = Charge5)
#elif LEGS == 6
#define CHARGES ch1, ch2, ch3, ch4, ch5, ch6
	RealType CHARGES
	parameter (ch1 = Charge1)
	parameter (ch2 = Charge2)
	parameter (ch3 = Charge3)
	parameter (ch4 = Charge4)
	parameter (ch5 = Charge5)
	parameter (ch6 = Charge6)
#else
#error No definition for LEGS external legs
#endif
	data charge /CHARGES/

	loglambda = log(4*(ESOFTMAX)**2/getlambda())

	SoftPhotonFactor = 0
	do i = 1, LEGS
	  do j = i, LEGS
	    r = charge(i)*charge(j)
	    if( r .ne. 0 ) SoftPhotonFactor = SoftPhotonFactor +
     &        sign(1, LEGS_IN - i)*sign(1, LEGS_IN - j)*
     &        r*PhaseSpaceIntegral(i, j)
	  enddo
	enddo
	SoftPhotonFactor = -4*pi*Alfa/(2*pi)**3*SoftPhotonFactor
	end


************************************************************************
** PhaseSpaceIntegral calculates the integral over the soft photon
** phase space.  See A. Denner, Fortschr. d. Physik, 41 (1993) 4.

	RealType function PhaseSpaceIntegral(i, j)
	implicit none
	integer i, j

#include "xsection.h"

	RealType mi, ki0, ki, kidelta, kilog
	RealType mj, kj0, kj, kjdelta, kjlog
	RealType kikj, alpha, d, vi

	RealType loglambda
	common /infrared/ loglambda

	ComplexType Pair
	Pair(i, j) = .5D0*(
     &    vec(1,1,i)*vec(2,2,j) + vec(2,2,i)*vec(1,1,j) -
     &    vec(1,2,i)*vec(2,1,j) - vec(2,1,i)*vec(1,2,j) )

	mi = momspec(SPEC_M, i)
	ki0 = momspec(SPEC_E, i)
	ki = momspec(SPEC_K, i)
	kidelta = momspec(SPEC_DELTAK, i)
	kilog = log(kidelta/(ki0 + ki))

	if( i .eq. j ) then
	  if( ki .lt. 1D-14 ) then
	    PhaseSpaceIntegral = pi*(loglambda - 2)
	  else
	    PhaseSpaceIntegral = pi*(loglambda + ki0/ki*kilog)
	  endif
	  return
	endif

	mj = momspec(SPEC_M, j)
	kj0 = momspec(SPEC_E, j)
	kj = momspec(SPEC_K, j)
	kjdelta = momspec(SPEC_DELTAK, j)
	kjlog = log(kjdelta/(kj0 + kj))

	kikj = Re(Pair(k0(i), k0(j)))

	d = sqrt((kikj - mi*mj)*(kikj + mi*mj))
	alpha = (kikj + d)/mi**2
	if( alpha*ki0 .lt. kj0 ) alpha = (kikj - d)/mi**2

	d = 1/((alpha*mi - mj)*(alpha*mi + mj))
	vi = 2*(alpha*ki0 - kj0)*d

	PhaseSpaceIntegral = 4*pi*alpha*kikj*d * Re(
     &    log(mi/mj*abs(alpha))*loglambda +
     &    .25D0*(kilog - kjlog)*(kilog + kjlog) +
     &    Li2omx((ki0 + ki)*alpha*vi) + Li2omx(kidelta*alpha*vi) -
     &    Li2omx((kj0 + kj)*vi) - Li2omx(kjdelta*vi) )
	end

