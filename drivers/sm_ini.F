* sm_ini.F
* all parameters for the Standard Model are initialized here
* this file is part of FormCalc
* last modified 14 Jun 01 th


	subroutine model_ini(Ecms, reset, *)
	implicit none
	double precision Ecms
	logical reset

	double precision ALPHAS2
	external ALPHAS2

#include "model.h"


* Warning: With the PDG values used here, the CKM matrix is
* not strictly unitary. In cases where GIM cancellations are
* important, this parametrization may give wrong results!

	data CKM /
     +    .97495D0, .2225D0,  .0035D0,
     +    .222D0,   .97415D0, .040D0,
     +    .009D0,   .039D0,   .99915D0 /

	data MLE  / ME,  MM,  ML  /
	data MLE2 / ME2, MM2, ML2 /

	data MQU  / MU,  MC,  MT  /
	data MQU2 / MU2, MC2, MT2 /

	data MQD  / MD,  MS,  MB  /
	data MQD2 / MD2, MS2, MB2 /

	Alfas = ALPHAS2(Ecms)
	Alfas2 = Alfas**2
	GS = sqrt(4*pi*Alfas)

	if(.not. reset) return

	EL = sqrt(4*pi*Alfa)

* boson masses

	MH2 = MH**2
	SW = sqrt(SW2)

#ifdef DEBUG

	print *, "============================"
	print *, "SM parameters"
	print *, "----------------------------"
	print *, "1/alpha =", 1/Alfa
	print *, "alpha_s(s) =", Alfas
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

#endif

	end

