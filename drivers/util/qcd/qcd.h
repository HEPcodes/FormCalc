* qcd.h
* constants concerning QCD running
* this file is part of FormCalc
* last modified 22 Dec 15 th


#define encLoop(n) (ibset(0,n) - 1)*8
#define nFlavors(i) iand(i,7)
#define ifLoop(i,n) ibits(i,n+2,1)

#include "const.h"

	RealType Mquark(6), MZ, AlfasMZ
	common /qcdvar/ Mquark, MZ, AlfasMZ

