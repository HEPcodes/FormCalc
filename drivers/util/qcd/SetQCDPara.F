* SetQCDPara.F
* set the relevant QCD parameters: (MZ, alpha_s(MZ)), quark masses
* this file is part of FormCalc
* last modified 22 Dec 15 th

#include "types.h"


	subroutine SetQCDPara(mx, asmx, m1, m2, m3, m4, m5, m6)
	implicit none
	RealType mx, asmx, m1, m2, m3, m4, m5, m6

#include "qcd.h"

	MZ = mx
	AlfasMZ = asmx

	Mquark(1) = m1
	Mquark(2) = m2
	Mquark(3) = m3
	Mquark(4) = m4
	Mquark(5) = m5
	Mquark(6) = m6
	end

