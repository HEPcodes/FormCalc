* MqRun.F
* compute the factor that runs a quark mass from Q_from to Q_to,
* code is transcribed from RunDec's mL2mH and mH2mL,
* see Eqs. (10) & (27) in hep-ph/0004189.
* last modified 24 Dec 15 th

#include "types.h"


	subroutine MqRun(m, Q_from, Q_to, nloop)
	implicit none
	RealType m, Q_from, Q_to
	integer nloop

#include "qcd.h"

	RealType as_from, as_to, r
	integer nf_from, nf_to, nf, iloop

	RealType eps
	parameter (eps = 1D-12)

#if 0
the following constants were produced with

B4 = 16 PolyLog[4, 1/2] + 2/3 Log[2]^4 - 2/3 Pi^2 Log[2]^2 - 13/180 Pi^4

c0[nf_] := 1/b0[nf];
c1[nf_] := (202/3 - nf 20/9)/(16 b0[nf]);
c2[nf_] := (1249 - nf (2216/27 + 160/3 zeta3) - nf^2 140/81)/(64 b0[nf]);
c3[nf_] := (4603055/162 + 135680/27 zeta3 - 8800 zeta5 +
  nf (-91723/3 - 34192 zeta3 + 7920 zeta4 + 18400 zeta5)/9 +
  nf^2 (5242/81 + 800/3 zeta3 - 160 zeta4)/3 +
  nf^3 (-332/9 + 64 zeta3)/27)/(256 b0[nf]);

A1[nf_] := (c1[nf] - b10[nf] c0[nf]);
A2[nf_] := (c2[nf] - b10[nf] c1[nf] + b10[nf]^2 c0[nf] - b20[nf] c0[nf]);
A3[nf_] := (c3[nf] - b10[nf] c2[nf] + b10[nf]^2 c1[nf] - b20[nf] c1[nf] -
  b10[nf]^3 c0[nf] + 2 b10[nf] b20[nf] c0[nf] - b30[nf] c0[nf]);

C1[nf_] := A1[nf] /Pi;
C2[nf_] := (A1[nf]^2 + A2[nf])/2 /Pi^2;
C3[nf_] := (A1[nf]^3/6 + A1[nf] A2[nf]/2 + A3[nf]/3) /Pi^3;

N[Table[Through[{C1,C2,C3}[nf]], {nf, 3, 6}] //Expand, 34]

#endif

	RealType B4, d2, d30, d31
	parameter (B4 = -1.762800087073770864061897634679819D0)
	parameter (d2 = 89/432D0/pi**2)
	parameter (d30 = (1871/2916D0 - B4/36D0 - 407/864D0*zeta3 +
     &    5/4D0*zeta4)/pi**3)
	parameter (d31 = (1327/11664D0 - 2/27D0*zeta3)/pi**3)

	RealType C_n(3,3:6)
	data C_n /
     &    .2849069968928990578578783418396862D0,
     &      .1389552059662087337299131752037630D0,
     &        .06294478093723101666667034039241161D0,
     &    .3228086659085215796955013077896918D0,
     &      .1407560845491501961210593657036523D0,
     &        .03517169812300769371007893557460479D0,
     &    .3741695606023071819494047392204481D0,
     &      .1520533792527639666058513363910428D0,
     &        .005562681911861937206888259088747334D0,
     &    .4449842286446869591905525628986626D0,
     &      .1817174007626149992000585582983386D0,
     &        -.02204174785022961834832249362399668D0 /

	RealType asthr(0:1,3:5)
	common /asthresh/ asthr

	RealType cfun, decoup, as
	integer s, nl
	cfun(as, nl) = (as/pi)**(4/(11 - 2/3D0*nl))*(1 +
     &    ifLoop(iloop,1)*as*(C_n(1,nl) +
     &      ifLoop(iloop,2)*as*(C_n(2,nl) +
     &        ifLoop(iloop,3)*as*C_n(3,nl))))
	decoup(as, s, nl) = 1 -
     &    ifLoop(iloop,2)*as**2*s*(d2 +
     &      ifLoop(iloop,3)*as*(d30 + nl*d31))

	call AlphaS(as_from, nf_from, Q_from, nloop)
	call AlphaS(as_to, nf_to, Q_to, nloop)

	iloop = encLoop(nloop)

	r = cfun(as_to, nf_to)/cfun(as_from, nf_from)

* nf_to > nf_from
	do nf = nf_from, nf_to - 1
	  r = r*cfun(asthr(0,nf), nf)/cfun(asthr(1,nf), nf + 1)*
     &      decoup(asthr(0,nf), +1, nf)
	enddo

* nf_to < nf_from
	do nf = nf_from - 1, nf_to, -1
	  r = r*cfun(asthr(1,nf), nf + 1)/cfun(asthr(0,nf), nf)*
     &      decoup(asthr(1,nf), -1, nf)
	enddo

	m = m*r
	end

