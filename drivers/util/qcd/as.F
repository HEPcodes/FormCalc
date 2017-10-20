* as.F
* command-line interface to AlphaS
* compile with e.g.
* gfortran -Idrivers/F -o as drivers/util/qcd/as.F MacOSX-x86-64/util.a
* last modified 22 Sep 17 th

#include "types.h"

	program as
	implicit none

#include "const.h"
#include "model_sm.h"

	RealType Q
	integer nloop, nf, fail

#if U77EXT
	integer iargc
	external iargc
#endif

	integer argc
	character argv*100

	MU = exp(11/5D0*log(MZ) - 4/5D0*log(MC) -
     &    1/5D0*(log(MS) + log(MB)) - 55/30D0 -
     &    9/5D0*pi/2D0*DeltaAlfa5Had/Alfa)
	MD = MU

	call SetQCDPara(MZ, AlfasMZ, MU, MD, MS, MC, MBatMB, MT)

	nloop = 3

	argc = iargc() 
	if( argc .ge. 1 .and. argc .le. 2 ) goto 1

999	print *, "Usage:  as scale [nloop]"
	print *, "computes alpha_s(scale) at nloop order"
	call exit(1)

1	call getarg(1, argv)
	read(argv, *, iostat=fail, err=999) Q
	if( argc .gt. 1 ) then
	  call getarg(2, argv)
	  read(argv, *, iostat=fail, err=999) nloop
	endif

	call AlphaS(Alfas, nf, Q, nloop)

	print *, "scale      =", Q
	print *, "loop order =", nloop
	print *, "alpha_s    =", Alfas
	print *, "nf         =", nf
	end

