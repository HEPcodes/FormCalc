* softphoton.F
* routines for the treatment of soft-photon radiation
* this file is part of FormCalc
* last modified 23 Jun 04 th


************************************************************************
** SoftPhotonFactor adds up the contributions from soft-photon emission
** off the different legs.  The result then has to be multiplied with
** the Born cross-section.

	double precision function SoftPhotonFactor()
	implicit none

#include HEADER
#include "looptools.h"

	double precision PhaseSpaceIntegral
	external PhaseSpaceIntegral

	integer i, j

	double precision loglambda
	common /infrared/ loglambda

	double precision charge(LEGS), r

* not a very nice solution, the following, but necessary
* because there may be no "/" in a constant in the data
* statement (such as in 2/3D0)
#if LEGS == 3
	double precision ch1, ch2, ch3
	parameter (ch1 = CHARGE1)
	parameter (ch2 = CHARGE2)
	parameter (ch3 = CHARGE3)
	data charge / ch1, ch2, ch3 /
#elif LEGS == 4
	double precision ch1, ch2, ch3, ch4
	parameter (ch1 = CHARGE1)
	parameter (ch2 = CHARGE2)
	parameter (ch3 = CHARGE3)
	parameter (ch4 = CHARGE4)
	data charge / ch1, ch2, ch3, ch4 /
#elif LEGS == 5
	double precision ch1, ch2, ch3, ch4, ch5
	parameter (ch1 = CHARGE1)
	parameter (ch2 = CHARGE2)
	parameter (ch3 = CHARGE3)
	parameter (ch4 = CHARGE4)
	parameter (ch5 = CHARGE5)
	data charge / ch1, ch2, ch3, ch4, ch5 /
#else
	error: No definition for LEGS external legs
#endif

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
** PhaseSpaceIntegral calculates the integral over the soft photon's
** phase space.  See A. Denner, Fortschr. d. Physik, 41 (1993) 4.

	double precision function PhaseSpaceIntegral(i, j)
	implicit none
	integer i, j

#include "model.h"
#include "util.h"

	double precision ki2, ki0, ki, kidelta, kilog
	double precision kj2, kj0, kj, kjdelta, kjlog
	double precision kikj, alpha, d, vi

	double precision loglambda
	common /infrared/ loglambda

	ki2 = momspec(1, i)
	ki0 = momspec(2, i)
	ki = momspec(3, i)
	kidelta = momspec(4, i)
	kilog = log(kidelta/(ki0 + ki))

	if( i .eq. j ) then
	  if( ki .lt. 1D-14 ) then
	    PhaseSpaceIntegral = pi*(loglambda - 2)
	  else
	    PhaseSpaceIntegral = pi*(loglambda + ki0/ki*kilog)
	  endif
	  return
	endif

	kj2 = momspec(1, j)
	kj0 = momspec(2, j)
	kj = momspec(3, j)
	kjdelta = momspec(4, j)
	kjlog = log(kjdelta/(kj0 + kj))

	kikj = dble(Pair(k(i), k(j)))

	d = sqrt(kikj**2 - ki2*kj2)
	alpha = (kikj + d)/ki2
	if( sign(kj0, alpha*ki0 - kj0) .ne. kj0 ) alpha = (kikj - d)/ki2

	d = 1/(alpha**2*ki2 - kj2)
	vi = 2*(alpha*ki0 - kj0)*d

	PhaseSpaceIntegral = 4*pi*alpha*kikj*d * (
     &    .5D0*log(ki2/kj2*alpha**2)*loglambda +
     &    .25D0*(kilog**2 - kjlog**2) +
     &    Li2(1 - (ki0 + ki)*alpha*vi) + Li2(1 - kidelta*alpha*vi) -
     &    Li2(1 - (kj0 + kj)*vi) - Li2(1 - kjdelta*vi) )

	end

