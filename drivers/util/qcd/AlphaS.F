* AlphaS.F
* compute the strong coupling constant at a given scale,
* code is transcribed from RunDec's AsRunDec,
* see Eqs. (5) & (25) in hep-ph/0004189
* this file is part of FormCalc
* last modified 22 Dec 15 th

#include "types.h"


	subroutine AlphaS(as, nf, Q, nloop)
	implicit none
	RealType as, Q
	integer nf, nloop

#include "qcd.h"

	integer iloop
	RealType t0, t1
	external asfun

	RealType eps, h
	parameter (eps = 1D-7)
	parameter (h = .1D0)

	RealType d2, d30, d31
	parameter (d2 = 7/24D0/pi**2)
	parameter (d30 = (58933/124416D0 + 2/9D0*zeta2*(3 + log2) +
     &    80507/27648D0*zeta3)/pi**3)
	parameter (d31 = (2479/31104D0 + zeta2/9D0)/pi**3)

	RealType asthr(0:1,3:5)
	common /asthresh/ asthr
	integer asx
	save asx

	RealType decoup
	integer s, nl
	decoup(as, s, nl) = 1 +
     &    ifLoop(iloop,2)*as**2*s*(d2 +
     &      ifLoop(iloop,3)*as*(d30 - nl*d31))

	if( Q .lt. Mquark(3) )
     &    stop "AlphaS not implemented below s-quark mass"

	iloop = encLoop(nloop)

	as = AlfasMZ
	t0 = log(MZ**2)
	t1 = t0

	if( Q .ge. Mquark(6) ) then
	  t0 = log(Mquark(6)**2)
	  if( .not. btest(asx, 5) ) then
	    asx = ibset(asx, 5)
	    call ODESolve(1, asfun, as, t1, t0, eps, h, iloop + 5)
	    asthr(0,5) = as
	    asthr(1,5) = as*decoup(as, +1, 5)
	  endif
	  as = asthr(1,5)
	  nf = 6
	  goto 1
	endif

	nf = 5
	if( Q .ge. Mquark(5) ) goto 1

	t0 = log(Mquark(5)**2)
	if( .not. btest(asx, 4) ) then
	  asx = ibset(asx, 4)
	  call ODESolve(1, asfun, as, t1, t0, eps, h, iloop + 5)
	  asthr(1,4) = as
	  asthr(0,4) = as*decoup(as, -1, 4)
	endif
	as = asthr(0,4)
	nf = 4
	if( Q .ge. Mquark(4) ) goto 1

	t1 = t0
	t0 = log(Mquark(4)**2)
	if( .not. btest(asx, 3) ) then
	  asx = ibset(asx, 3)
	  call ODESolve(1, asfun, as, t1, t0, eps, h, iloop + 4)
	  asthr(1,3) = as
	  asthr(0,3) = as*decoup(as, -1, 3)
	endif
	as = asthr(0,3)
	nf = 3

1	call ODESolve(1, asfun, as, t0, log(Q**2), eps, h, iloop + nf)
	end


************************************************************************

	subroutine asfun(das_dt, as, t, iloop)
	implicit none
	RealType das_dt, as, t
	integer iloop

#include "qcd.h"

#if 0
the following constants were produced with

b0[nf_] := (11 - 2/3 nf)/4;
b10[nf_] := (102 - 38/3 nf)/(16 b0[nf]);
b20[nf_] := (2857/2 - 5033/18 nf + 325/54 nf^2)/(64 b0[nf]);
b30[nf_] := (149753/6 + 3564 zeta3 -
  (1078361/162 + 6508/27 zeta3) nf +
  (50065/162 + 6472/81 zeta3) nf^2 +
  1093/729 nf^3)/(256 b0[nf]);

Table[Through[{b0,b10,b20,b30}[n]], {n, 3, 6}] //Expand

#endif

	integer nf

	RealType b3_0, b3_10, b3_20, b3_30
	RealType b4_0, b4_10, b4_20, b4_30
	RealType b5_0, b5_10, b5_20, b5_30
	RealType b6_0, b6_10, b6_20, b6_30

	parameter (b3_0 = 9/4D0/pi)
	parameter (b3_10 = 16/9D0/pi)
	parameter (b3_20 = 3863/864D0/pi**2)
	parameter (b3_30 = (140599/10368D0 +
     &                      445/72D0*zeta3)/pi**3)

	parameter (b4_0 = 25/12D0/pi)
	parameter (b4_10 = 77/50D0/pi)
	parameter (b4_20 = 21943/7200D0/pi**2)
	parameter (b4_30 = (4918247/777600D0 +
     &                      15707/2160D0*zeta3)/pi**3)

	parameter (b5_0 = 23/12D0/pi)
	parameter (b5_10 = 29/23D0/pi)
	parameter (b5_20 = 9769/6624D0/pi**2)
	parameter (b5_30 = (-26017/31104D0 +
     &                      11027/1242D0*zeta3)/pi**3)

	parameter (b6_0 = 7/4D0/pi)
	parameter (b6_10 = 13/14D0/pi)
	parameter (b6_20 = -65/224D0/pi**2)
	parameter (b6_30 = (-63559/8064D0 +
     &                      11237/1008D0*zeta3)/pi**3)

	RealType b_n(0:3,3:6)
	data b_n /
     &    b3_0, b3_10, b3_20, b3_30,
     &    b4_0, b4_10, b4_20, b4_30,
     &    b5_0, b5_10, b5_20, b5_30,
     &    b6_0, b6_10, b6_20, b6_30 /

	nf = nFlavors(iloop)

	das_dt = -as**2*b_n(0,nf)*(1 +
     &    ifLoop(iloop,1)*as*(b_n(1,nf) +
     &      ifLoop(iloop,2)*as*(b_n(2,nf) +
     &        ifLoop(iloop,3)*as*b_n(3,nf))))
	end

