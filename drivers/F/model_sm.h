* model_sm.h
* declarations for model_sm.F
* this file is part of FormCalc
* last modified 06 Jun 17 cs


#ifndef CKMC
#define CKMC(i,j) Conjugate(CKM(i,j))
#endif

	RealType CKMlambda, CKMA, CKMrhobar, CKMetabar
	parameter (CKMlambda = .22535D0)
	parameter (CKMA = .811D0)
	parameter (CKMrhobar = .131D0)
	parameter (CKMetabar = .345D0)

	RealType MZ, MZ2, MW, MW2, CW, CW2, SW2
	RealType GammaW, GammaZ
	parameter (MZ = 91.1876D0, MZ2 = MZ**2)
	parameter (MW = 80.385D0, MW2 = MW**2)
	parameter (CW = MW/MZ, CW2 = CW**2)
	parameter (SW2 = (1 - CW)*(1 + CW))
	parameter (GammaW = 2.1180D0)   
	parameter (GammaZ = 2.4952D0)

	RealType GF, Alfa0, AlfaGF, AlfaMZ, Alfa, Alfa2
	RealType DeltaAlfa5Had, AlfasMZ
	parameter (GF = 1.1663787D-5)
	parameter (Alfa0 = 1/137.035999074D0)
	parameter (AlfaGF = sqrt2/pi*GF*MW2*SW2)
	parameter (AlfaMZ = 1/127.944D0)
	parameter (Alfa = Alfa0, Alfa2 = Alfa**2)
	parameter (DeltaAlfa5Had = .027547D0)
	parameter (AlfasMZ = .1184D0)

	RealType MH_
	parameter (MH_ = 125)

	RealType ME, ME2, MM, MM2, ML, ML2
	parameter (ME = .510998928D-3, ME2 = ME**2)
	parameter (MM = 105.6583715D-3, MM2 = MM**2)
	parameter (ML = 1776.82D-3, ML2 = ML**2)

	RealType MC, MC2, MT, MT2
	parameter (MC = 1.275D0, MC2 = MC**2)
	parameter (MT = 173.21D0, MT2 = MT**2)

	RealType MS, MS2, MB1S, MBatMB
	parameter (MS = 95D-3, MS2 = MS**2)
	parameter (MB1S = 4.66D0)
	parameter (MBatMB = 4.18D0)

	ComplexType CKM(3,3)
	RealType Mf(4,3), Mf2(4,3), MH, MH2
	RealType EL, ELMZ, GS, Alfas, Alfas2, SW

	common /smpara/ CKM
	common /smpara/ Mf, Mf2, MH, MH2
	common /smpara/ EL, ELMZ, GS, Alfas, Alfas2, SW

	RealType MU, MU2, MD, MD2, MB, MB2
	equivalence (Mf(3,1), MU)
	equivalence (Mf2(3,1), MU2)
	equivalence (Mf(4,1), MD)
	equivalence (Mf2(4,1), MD2)
	equivalence (Mf(4,3), MB)
	equivalence (Mf2(4,3), MB2)

