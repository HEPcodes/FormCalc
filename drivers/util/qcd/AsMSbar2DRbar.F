* AsMSbar2DRbar.F
* turn an MSbar into a DRbar alpha_s via Eq. (18) in arXiv:0706.2953
* this file is part of FormCalc
* last modified 24 Dec 15 th

#include "types.h"


	subroutine AsMSbartoDRbar(as, nf, Q, Qsusy, nloop)
	implicit none
	RealType as, Q, Qsusy
	integer nf, nloop

#include "qcd.h"

	RealType Lq, Lt, z1
	RealType Ltq, tq, z2
	integer iloop

	Lq = 2*log(Q/Qsusy)

	if( nf .eq. 6 ) then
	  z1 = -Lq
	  z2 = Lq*(Lq - 5/2D0) - 65/32D0
	else if( nf .eq. 5 ) then
	  Lt = 2*log(Q/Mquark(6))
	  z1 = -Lt/6D0 - Lq
	  tq = Mquark(6)/Qsusy
	  Ltq = 2*log(tq)
	  z2 = -215/96D0 - 19/24D0*Lt - 5/2D0*Lq + z1**2 +
     &      tq**2*(5/48D0 + 3/8D0*Ltq +
     &        tq*(-7*pi/36D0 +
     &          tq*(881/7200D0 - 1/80D0*Ltq + 7*pi/288D0*tq)))
	else
	  print *, "DRbar conversion not available for nf < 5"
	  z1 = 0
	  z2 = 0
	endif

	iloop = encLoop(nloop)

	as = as*(1 +
     &    ifLoop(iloop,1)*as/pi*(.25D0 - z1 +
     &      ifLoop(iloop,2)*as/pi*(11/8D0 - nf/12D0 -
     &        2*z1*(.25D0 - z1) - z2)))
        end

