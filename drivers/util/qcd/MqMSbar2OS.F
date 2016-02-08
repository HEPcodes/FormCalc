* MqMSbar2OS.F
* compute the conversion factor from mq^MSbar(Q) -> mq^OS,
* code is transcribed from RunDec's mMS2mOS function.
* this file is part of FormCalc
* last modified 16 Jan 16 th

#include "types.h"


	subroutine MqMSbartoOS(m, mMS, Q, nloop)
	implicit none
	RealType m, mMS, Q
	integer nloop

#include "qcd.h"

	RealType as, t
	integer nf, iloop

* the following constants are from RunDec's
*   osfromms + api^3 osfromms3
* note: c is indexed for nf (3:6) but computed with nl = nf - 1

	RealType c(7,3:6)
	data c /
* nf = 3, nl = 2
     &    1.151075764833548879586528900462204D0,
     &      .6205922498093188500937617121595818D0,
     &        .1646469234187988785963041277158074D0,
     &    4.511843079333330161746044341013300D0,
     &      2.111166219183803548840895807401627D0,
     &        .5640658990973118993817148590010253D0,
     &          .09608269633224014486192403196801266D0,
* nf = 4, nl = 3
     &    1.045563236787674006166282776240741D0,
     &      .5840040446051413215168052393338468D0,
     &        .1562034914486040643093141724483301D0,
     &    3.757425572585101249486168152595144D0,
     &      1.727979654979010705731150972919183D0,
     &        .4914252902929065684177272264409039D0,
     &          .08563081017333290299660207044857072D0,
* nf = 5, nl = 4
     &    .9400507087417991327460366520192781D0,
     &      .5474158394009637929398487665081118D0,
     &        .1477600594784092500223242171808528D0,
     &    3.045108622173340689987629587919559D0,
     &      1.365441563609897623719760376122610D0,
     &        .4226668106332382130037163224451466D0,
     &          .07577617465207750352358422101595404D0,
* nf = 6, nl = 5
     &    .8345381806959242593257905277978152D0,
     &      .5108276341967862643628922936823768D0,
     &        .1393166275082144357353342619133755D0,
     &    2.374892228098048483250428646986546D0,
     &      1.023551945076464302806724017011909D0,
     &        .3577904601183068331396821470137534D0,
     &          .06651878976847394644287048367016261D0 /

	RealType dc1, dc2, dc3
	parameter (dc1 = 4/(3*pi**2)*pi**2/8D0)
	parameter (dc2 = -4/(3*pi**2)*.597D0)
	parameter (dc3 = 4/(3*pi**2)*.230D0)

	RealType Delta, DLite, x
	integer i
	Delta(x) = (sign(.5D0, 1 - x) + .5D0)*x*(dc1 + x*(dc2 + x*dc3))
	DLite(i) = Delta(Mquark(i)/mMS)

	call AlphaS(as, nf, Q, nloop)

	iloop = encLoop(nloop)
	t = 2*log(Q/mMS)

	m = m*(1 + ifLoop(iloop,1)*as*(4/3D0/pi + 1/pi*t +
     &    ifLoop(iloop,2)*as*(
     &        DLite(1) + DLite(2) + DLite(3) +
     &        DLite(4) + DLite(5) + DLite(6) +
     &        c(1,nf) + t*(c(2,nf) + t*c(3,nf)) +
     &      ifLoop(iloop,3)*as*(
     &        c(4,nf) + t*(c(5,nf) + t*(c(6,nf) + t*c(7,nf))) ) )))
	end

