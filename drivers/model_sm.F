* model_sm.F
* all parameters for the Standard Model are initialized here
* this file is part of FormCalc
* last modified 5 May 11 th


	subroutine ModelDefaults(argc, argv)
	implicit none
	integer argc
	character*128 argv(*)

#include "decl.h"

#ifdef NO_CKM

	CKM(1,1) = 1
	CKM(2,2) = 1
	CKM(3,3) = 1
	CKM(2,1) = 0
	CKM(3,1) = 0
	CKM(1,2) = 0
	CKM(3,2) = 0
	CKM(1,3) = 0
	CKM(2,3) = 0

#else

	double precision s12, s23, s13
	double precision c12, c23, c13
	double complex s13c

	s12 = CKMlambda
	c12 = sqrt((1 - s12)*(1 + s12))

	s23 = s12**2*CKMA
	c23 = sqrt((1 - s23)*(1 + s23))

	s13c = s12*s23*DCMPLX(CKMrhobar, CKMetabar)
	s13c = c23/c12*s13c/(1 - s23/s12*s13c)
	s13 = abs(s13c)
	c13 = sqrt((1 - s13)*(1 + s13))

	CKM(1,1) = c12*c13
	CKM(2,1) = -s12*c23 - c12*s23*s13c
	CKM(3,1) = s12*s23 - c12*c23*s13c
	CKM(1,2) = s12*c13
	CKM(2,2) = c12*c23 - s12*s23*s13c
	CKM(3,2) = -c12*s23 - s12*c23*s13c
	CKM(1,3) = DCONJG(s13c)
	CKM(2,3) = s23*c13
	CKM(3,3) = c23*c13

#endif

	Mf(1,1) = 0
	Mf(2,1) = ME
	Mf(3,1) = MU
	Mf(4,1) = MD
	Mf(1,2) = 0
	Mf(2,2) = MM
	Mf(3,2) = MC
	Mf(4,2) = MS
	Mf(1,3) = 0
	Mf(2,3) = ML
	Mf(3,3) = MT
	Mf(4,3) = MB

	Mf2(1,1) = 0
	Mf2(2,1) = ME2
	Mf2(3,1) = MU2
	Mf2(4,1) = MD2
	Mf2(1,2) = 0
	Mf2(2,2) = MM2
	Mf2(3,2) = MC2
	Mf2(4,2) = MS2
	Mf2(1,3) = 0
	Mf2(2,3) = ML2
	Mf2(3,3) = MT2
	Mf2(4,3) = MB2

	EL = sqrt(4*pi*Alfa)
	SW = sqrt(SW2)
	end


************************************************************************

	subroutine ModelConstIni(fail)
	implicit none
	integer fail

#include "decl.h"

	fail = 0

	MH2 = MH**2
	end


************************************************************************

	subroutine ModelVarIni(fail, scale)
	implicit none
	integer fail
	double precision scale

#include "decl.h"

	double precision alphasPDF
	external alphasPDF

	fail = 0

	Alfas = alphasPDF(scale)
	Alfas2 = Alfas**2
	GS = sqrt(4*pi*Alfas)
	end


************************************************************************

	subroutine ModelDigest
	implicit none

#include "decl.h"

	INFO "============================"
	INFO "SM parameters"
	INFO "----------------------------"
	INFO "1/Alfa =", 1/Alfa
	INFO "Alfas =", Alfas
	INFO "GF =", GF
	INFO "SW2 =", SW2
	INFO "----------------------------"
	INFO "Gauge-boson masses:"
	INFO "MW =", MW
	INFO "MZ =", MZ
	INFO "----------------------------"
	INFO "Fermion masses:"
	INFO "ME =", ME
	INFO "MM =", MM
	INFO "ML =", ML
	INFO "MU =", MU
	INFO "MC =", MC
	INFO "MT =", MT
	INFO "MD =", MD
	INFO "MS =", MS
	INFO "MB =", MB
	INFO "----------------------------"
	INFO "Higgs mass:"
	INFO "MH =", MH
	INFO "============================"
	end


#ifndef HAVE_ALPHASPDF

************************************************************************
** alphasPDF (same name as in LHAPDF) computes the strong coupling
** constant at a given scale.  Code is transcribed from Java applet by
** I. Hinchliffe at http://www-theory.lbl.gov/~ianh/alpha/alpha.html.

	double precision function alphasPDF(scale)
	implicit none
	double precision scale

#include "decl.h"

* the betas and deltas are calculated as follows:
*
* beta0[nf_] := 11 - 2/3 nf
* beta1[nf_] := 51 - 19/3 nf
* beta2[nf_] := 2857 - 5033/9 nf + 325/27 nf^2
*
* delta0[nf_] := (beta0[nf] - beta0[nf + 1])/(2 beta0[nf])
* delta10[nf_] := beta1[nf]/beta0[nf]^2
* delta11[nf_] := (beta1[nf]/beta0[nf] -
*   beta1[nf + 1]/beta0[nf + 1])/beta0[nf]
* delta20[nf_] := (28/9 -
*   8 (beta1[nf + 1]/beta0[nf + 1])^2 + beta2[nf + 1]/beta0[nf + 1] +
*   8 (beta1[nf]/beta0[nf])^2 - beta2[nf]/beta0[nf])/
*   (4 beta0[nf + 1] beta0[nf])
* delta21[nf_] := 2 beta1[nf + 1]/beta0[nf + 1]^2 delta11[nf]

	double precision beta3_0, beta3_1, beta3_2
	parameter (beta3_0 = 9)
	parameter (beta3_1 = 32)
	parameter (beta3_2 = 3863/3D0)

	double precision beta4_0, beta4_1, beta4_2
	parameter (beta4_0 = 25/3D0)
	parameter (beta4_1 = 77/3D0)
	parameter (beta4_2 = 21943/27D0)

	double precision beta5_0, beta5_1, beta5_2
	parameter (beta5_0 = 23/3D0)
	parameter (beta5_1 = 58/3D0)
	parameter (beta5_2 = 9769/27D0)

	double precision beta6_0, beta6_1, beta6_2
	parameter (beta6_0 = 7)
	parameter (beta6_1 = 13)
	parameter (beta6_2 = -65)

	double precision delta54_0, delta54_10, delta54_11
	double precision delta54_20, delta54_21
	parameter (delta54_0  =  .040000000000000000000D0)
	parameter (delta54_10 = -.030817842663873271183D0)
	parameter (delta54_11 =  .066991304347826086957D0)
	parameter (delta54_20 = -.086878340429029341662D0)
	parameter (delta54_21 =  .044069893975507520342D0)

	double precision delta43_0, delta43_10, delta43_11
	double precision delta43_20, delta43_21
	parameter (delta43_0  =  .037037037037037037037D0)
	parameter (delta43_10 = -.030404361930322301228D0)
	parameter (delta43_11 =  .052839506172839506173D0)
	parameter (delta43_20 = -.057313053497942386831D0)
	parameter (delta43_21 =  .039058962962962962963D0)

	double precision beta0, beta1, beta2, lam, t, b, c, l

	double precision lam5, asMZ
	save lam5, asMZ
	data lam5 /204D-3/, asMZ /.1176D0/
* the 5-flavour scale lam5 reproduces alpha_s(MZ) = .1176

	do while( abs(asMZ - AlfasMZ) .gt. 1D-4 )
	  t = 2*log(MZ/lam5)
	  c = 2*beta5_1/beta5_0**2/t
	  l = c*log(t)
	  b = c**2*(beta5_2*beta5_0/(8*beta5_1**2) - 1)
	  asMZ = 4*pi/beta5_0/t*(1 + l*(l - c - 1) + b)
* Newton step is (as - AlfasMZ)/as'
	  lam5 = lam5*(1 - beta5_0/(8*pi)*t**2*(asMZ - AlfasMZ)/
     &      (1 + c*(c + 1) + l*(3*l - 5*c - 2) + 3*b))
	enddo

	if( scale .lt. .5D0 ) then
	  alphasPDF = 1D10
	  return
	endif

	if( scale .gt. MT ) then
c	  nf = 6
	  beta0 = beta6_0
	  beta1 = beta6_1
	  beta2 = beta6_2
	  t = MT/lam5
	  lam = lam5*t**(-2/21D0)*(2*log(t))**(-107/1127D0)
	else if( scale .gt. MB ) then
c	  nf = 5
	  beta0 = beta5_0
	  beta1 = beta5_1
	  beta2 = beta5_2
	  lam = lam5
	else
	  t = 2*log(MB/lam5)
	  l = log(t)
	  lam = lam5*exp(delta54_0*t +
     &      delta54_11*l + delta54_10 +
     &      (delta54_21*l + delta54_20)/t)
	  if( scale .gt. MC ) then
c	    nf = 4
	    beta0 = beta4_0
	    beta1 = beta4_1
	    beta2 = beta4_2
	  else
c	    nf = 3
	    beta0 = beta3_0
	    beta1 = beta3_1
	    beta2 = beta3_2
	    t = 2*log(MC/lam)
	    l = log(t)
	    lam = lam*exp(delta43_0*t +
     &        delta43_11*l + delta43_10 +
     &        (delta43_21*l + delta43_20)/t)
	  endif
	endif

	t = 2*log(scale/lam)
	c = 2*beta1/(beta0**2*t)
	l = c*log(t)
	b = c**2*(beta2*beta0/(8*beta1**2) - 1)
	alphasPDF = 4*pi/(beta0*t)*(1 + l*(l - c - 1) + b)
	end

#endif

