* softphoton.F
* routines for the treatment of soft-photon radiation
* this file is part of FormCalc
* last modified 8 Jun 15 th


************************************************************************
** SoftPhotonFactor adds up the contributions from soft-photon emission
** off the different legs.  The result then has to be multiplied with
** the Born cross-section.

	subroutine SoftPhotonFactor(fac)
	implicit none
	RealType fac

#include "xsection.h"

	integer i, j
	RealType r, loglambda, chsi, chsj, psint
	RealType mi, ki0, ki, kidelta, kilog
	RealType mj, kj0, kj, kjdelta, kjlog
	RealType kikj, alpha, d, vi

	ComplexType Pair
	Pair(i, j) = .5D0*(
     &    vec(1,1,i)*vec(2,2,j) + vec(2,2,i)*vec(1,1,j) -
     &    vec(1,2,i)*vec(2,1,j) - vec(2,1,i)*vec(1,2,j) )

	loglambda = log(4*(ESOFTMAX)**2/getlambda())

* psint = integral over the soft-photon phase space,
* from A. Denner, Fortschr. d. Physik, 41 (1993) 4.

	fac = 0
	do i = 1, LEGS
	  if( charge(i,pid) .eq. 0 ) goto 1
	  chsi = sign(1, LEGS_IN - i)*charge(i,pid)

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
	    if( charge(j,pid) .eq. 0 ) goto 2
	    chsj = sign(1, LEGS_IN - j)*charge(j,pid)

	    mj = momspec(SPEC_M,j)
	    kj0 = momspec(SPEC_E,j)
	    kj = momspec(SPEC_K,j)
	    kjdelta = momspec(SPEC_DELTAK,j)
	    kjlog = log(kjdelta/(kj0 + kj))

	    kikj = Re(Pair(k0(i), k0(j)))

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
	    fac = fac + chsi*chsj*psint

2	    continue
	  enddo

1	  continue
	enddo

	fac = -4*pi*Alfa/(2*pi)**3*fac
	end

