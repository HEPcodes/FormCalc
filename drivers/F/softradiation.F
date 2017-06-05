* softradiation.F
* routines for the treatment of soft-photon or soft-gluon radiation
* this file is part of FormCalc
* last modified 18 Jun 16 th


************************************************************************
** SoftRadiationFactor adds up the contributions from the emission of
** a soft particle off the different legs.  The result then has to be
** multiplied with the Born cross-section.

	subroutine SoftRadiationFactor(fac, chg, weight)
	implicit none
	RealType fac, chg(LEGS), weight(LEGS,LEGS)

#include "xsection.h"

	integer i, j
	RealType loglambda, chsi, chsj, psint
	RealType mi, ki0, ki, kidelta, kilog
	RealType mj, kj0, kj, kjdelta, kjlog
	RealType kikj, alpha, d, vi

#include "inline.h"

	loglambda = log(4*(ESOFTMAX)**2/getlambda())

* psint = integral over the soft-photon phase space,
* from A. Denner, Fortschr. d. Physik, 41 (1993) 4.

	fac = 0
	do i = 1, LEGS
	  if( chg(i) .eq. 0 ) goto 1
	  chsi = sign(1, LEGS_IN - i)*chg(i)

	  mi = momspec(SPEC_M,i)
	  ki0 = momspec(SPEC_E,i)
	  ki = momspec(SPEC_K,i)
	  kidelta = momspec(SPEC_DELTAK,i)
	  kilog = log(kidelta/(ki0 + ki))

	  if( ki .lt. 1D-14 ) then
	    psint = pi*(loglambda - 2)
	  else
	    psint = pi*(loglambda + ki0/ki*kilog)
	  endif
	  fac = fac + chsi**2*psint

	  do j = i + 1, LEGS
	    if( chg(j) .eq. 0 ) goto 2
	    chsj = sign(1, LEGS_IN - j)*chg(j)

	    mj = momspec(SPEC_M,j)
	    kj0 = momspec(SPEC_E,j)
	    kj = momspec(SPEC_K,j)
	    kjdelta = momspec(SPEC_DELTAK,j)
	    kjlog = log(kjdelta/(kj0 + kj))

	    kikj = Re(Pair0(k0(i), k0(j)))

	    d = sqrt((kikj - mi*mj)*(kikj + mi*mj))
	    alpha = (kikj + d)/mi**2
	    if( alpha*ki0 .lt. kj0 ) alpha = (kikj - d)/mi**2

	    d = 1/((alpha*mi - mj)*(alpha*mi + mj))
	    vi = 2*(alpha*ki0 - kj0)*d

	    psint = 4*pi*alpha*kikj*d * Re(
     &        log(mi/mj*abs(alpha))*loglambda +
     &        .25D0*(kilog - kjlog)*(kilog + kjlog) +
     &        Li2omx((ki0 + ki)*alpha*vi) - Li2omx((kj0 + kj)*vi) +
     &        Li2omx(kidelta*alpha*vi) - Li2omx(kjdelta*vi) )
	    fac = fac + chsi*chsj*psint*weight(i,j)

2	    continue
	  enddo

1	  continue
	enddo
	end

