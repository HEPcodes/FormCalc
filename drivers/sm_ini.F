* sm_ini.F
* all parameters for the Standard Model are initialized here
* this file is part of FormCalc
* last modified 12 Feb 03 th


	subroutine ModelConstIni(*)
	implicit none

#include "model.h"

* Warning: With the PDG values used here, the CKM matrix is
* not strictly unitary.  In cases where GIM cancellations are
* important, this parametrization may give wrong results!

	data CKM /
     &    .97495D0, .2225D0,  .0035D0,
     &    .222D0,   .97415D0, .040D0,
     &    .009D0,   .039D0,   .99915D0 /

	data MLE  / ME,  MM,  ML  /
	data MLE2 / ME2, MM2, ML2 /

	data MQU  / MU,  MC,  MT  /
	data MQU2 / MU2, MC2, MT2 /

	data MQD  / MD,  MS,  MB  /
	data MQD2 / MD2, MS2, MB2 /

	EL = sqrt(4*pi*Alfa)
	MH2 = MH**2
	SW = sqrt(SW2)
	end


************************************************************************

	subroutine ModelVarIni(sqrtS, *)
	implicit none
	double precision sqrtS

#include "model.h"

	double precision ALPHAS2
	external ALPHAS2

	Alfas = ALPHAS2(sqrtS)
	Alfas2 = Alfas**2
	GS = sqrt(4*pi*Alfas)
	end


************************************************************************

	subroutine ModelDigest
	implicit none

#include "model.h"

	print *, "============================"
	print *, "SM parameters"
	print *, "----------------------------"
	print *, "1/Alfa =", 1/Alfa
	print *, "Alfas =", Alfas
	print *, "GF =", GF
	print *, "SW2 =", SW2
	print *, "----------------------------"
	print *, "Gauge-boson masses:"
	print *, "MW =", MW
	print *, "MZ =", MZ
	print *, "----------------------------"
	print *, "Fermion masses:"
	print *, "ME =", ME
	print *, "MM =", MM
	print *, "ML =", ML
	print *, "MU =", MU
	print *, "MC =", MC
	print *, "MT =", MT
	print *, "MD =", MD
	print *, "MS =", MS
	print *, "MB =", MB
	print *, "----------------------------"
	print *, "Higgs mass:"
	print *, "MH =", MH
	print *, "============================"
	end


#ifndef CERNLIB

************************************************************************
** ALPHAS2 (same name as CERNLIB routine) computes the strong coupling
** constant at a given scale.  Code is transcribed from I. Hinchliffe's
** Java applet at http://www-theory.lbl.gov/~ianh/alpha/alpha.html.

	double precision function ALPHAS2(scale)
	implicit none
	double precision scale

#include "model.h"

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
	parameter (delta54_0  =  0.040000000000000000000D0)
	parameter (delta54_10 = -0.030817842663873271183D0)
	parameter (delta54_11 =  0.066991304347826086957D0)
	parameter (delta54_20 = -0.086878340429029341662D0)
	parameter (delta54_21 =  0.044069893975507520342D0)

	double precision delta43_0, delta43_10, delta43_11
	double precision delta43_20, delta43_21
	parameter (delta43_0  =  0.037037037037037037037D0)
	parameter (delta43_10 = -0.030404361930322301228D0)
	parameter (delta43_11 =  0.052839506172839506173D0)
	parameter (delta43_20 = -0.057313053497942386831D0)
	parameter (delta43_21 =  0.039058962962962962963D0)

	double precision lambda5
	parameter (lambda5 = 199D-3)
* The 5-flavour scale lambda5 is chosen to reproduce the
* "world average" alpha_s(MZ) = .1172.  The 1-sigma interval
* on alpha_s(MZ) corresponds to 177 MeV < lambda5 < 223 MeV.

	double precision b0, b1, b2, lambda, t, logt, c

	if( scale .gt. MT ) then
	  b0 = beta6_0
	  b1 = beta6_1
	  b2 = beta6_2
	  t = MT/lambda5
	  lambda = lambda5*t**(-2/21D0)*(2*log(t))**(-107/1127D0)
	else if( scale .gt. MB ) then
	  b0 = beta5_0
	  b1 = beta5_1
	  b2 = beta5_2
	  lambda = lambda5
	else
	  t = 2*log(MB/lambda5)
	  logt = log(t)
	  lambda = lambda5*exp(delta54_0*t +
     &      delta54_11*logt + delta54_10 +
     &      (delta54_21*logt + delta54_20)/t)
	  if( scale .gt. MC ) then
	    b0 = beta4_0
	    b1 = beta4_1
	    b2 = beta4_2
	  else
	    b0 = beta3_0
	    b1 = beta3_1
	    b2 = beta3_2
	    t = 2*log(MC/lambda)
	    logt = log(t)
	    lambda = lambda*exp(delta43_0*t +
     &        delta43_11*logt + delta43_10 +
     &        (delta43_21*logt + delta43_20)/t)
	  endif
	endif

	t = 2*log(scale/lambda)
	logt = log(t)
	c = 2*b1/b0**2/t

	ALPHAS2 = 4*pi/(t*b0)*(1 - logt*c +
     &    ((logt - .5D0)**2 + b2*b0/(8*b1**2) - 5/4D0)*c**2)
	end

#endif

