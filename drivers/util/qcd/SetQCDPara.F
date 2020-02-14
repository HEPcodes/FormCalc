* SetQCDPara.F
* set the relevant QCD parameters: (MZ, alpha_s(MZ)), quark masses
* this file is part of FormCalc
* last modified 25 Mar 19 th

#include "types.h"


	subroutine SetQCDPara(mx, asmx, m1, m2, m3, m4, m5, m6)
	implicit none
	RealType mx, asmx, m1, m2, m3, m4, m5, m6

#include "qcd.h"

	MZ = mx
	AlfasMZ = asmx

	Mquark = [m1, m2, m3, m4, m5, m6]
	end

