* model.h
* common blocks for the model parameters
* this file is part of FormCalc
* last modified 16 Nov 06 th


	double precision pi, degree, sqrt2, hbar_c2

	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (degree = pi/180D0)
	parameter (sqrt2 = 1.41421356237309504880168872421D0)

	parameter (hbar_c2 = 3.8937966D8)
*         = hbar c^2 in picobarn

	double complex bogus, cI

	parameter (bogus = (-1D123, -2D123))
*	  some weird number likely to noticeably distort the final result;
*	  used for initializing arrays to check that all components
*	  have been calculated

	parameter (cI = (0D0, 1D0))

#ifndef WARN
#define WARN print *,
#define INFO print *,
#define Digit(i) char(i+48)
#define Polar(r,theta) r*exp(cI*degree*theta)
#endif

	double precision Divergence
	common /renorm/ Divergence


* SM parameters

	double precision MZ, MZ2, MW, MW2, CW, CW2, SW2
	parameter (MZ = 91.1875D0, MZ2 = MZ**2)
	parameter (MW = 80.450D0, MW2 = MW**2)
	parameter (CW = MW/MZ, CW2 = CW**2)
	parameter (SW2 = 1 - CW2)

	double precision GF, Alfa, Alfa2, AlfaMZ
	parameter (GF = 1.16639D-5)
	parameter (Alfa = 1/137.0359895D0, Alfa2 = Alfa**2)
c	parameter (Alfa = sqrt2/pi*GF*MW2*SW2, Alfa2 = Alfa**2)
	parameter (AlfaMZ = 1/127.934D0)

	double precision ME, ME2, MM, MM2, ML, ML2
	parameter (ME = .51099907D-3, ME2 = ME**2)
	parameter (MM = 105.658389D-3, MM2 = MM**2)
	parameter (ML = 1777D-3, ML2 = ML**2)

	double precision MU, MU2, MC, MC2, MT, MT2
	parameter (MU = 53.8D-3, MU2 = MU**2)
	parameter (MC = 1.50D0, MC2 = MC**2)
	parameter (MT = 171.4D0, MT2 = MT**2)

	double precision MD, MD2, MS, MS2, MB, MB2, MBatMB
	parameter (MD = 53.8D-3, MD2 = MD**2)
	parameter (MS = 150D-3, MS2 = MS**2)
	parameter (MB = 4.7D0, MB2 = MB**2)
	parameter (MBatMB = 4.25D0)

	double complex CKM(3,3)
	double precision Mf(4,3), Mf2(4,3)
	double precision MH, MH2, MG0, MG02, MGp, MGp2
	double precision EL, GS, Alfas, Alfas2, AlfasMT, SW
	logical sm_digest

	common /sm_para/ CKM
	common /sm_para/ Mf, Mf2
	common /sm_para/ MH, MH2, MG0, MG02, MGp, MGp2
	common /sm_para/ EL, GS, Alfas, Alfas2, AlfasMT, SW
	common /sm_para/ sm_digest


* MSSM parameters

	double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
	double complex USf(2,2,4,3), Af(2:4,3,3), Xf(2:4,3,3)
	double complex Atau, At, Ab, MUE, M_1, M_2, M_3
	double precision MCha(2), MCha2(2), MNeu(4), MNeu2(4)
	double precision MSS(2,2:4,3), MSS2(2,2:4,3), DSf(2,4)
	double precision MSf(2,4,3), MSf2(2,4,3), MSusy, MGl, MGl2
	double precision Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2
	double precision Mh02tree, MHH2tree, MHp2tree
	double precision CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	double precision CA, SA, CA2, SA2, C2A, S2A
	double precision CAB, SAB, CBA, SBA
	logical mssm_digest

	common /mssm_para/ UCha, VCha, ZNeu
	common /mssm_para/ USf, Af, Xf
	common /mssm_para/ Atau, At, Ab, MUE, M_1, M_2, M_3
	common /mssm_para/ MCha, MCha2, MNeu, MNeu2
	common /mssm_para/ MSS, MSS2, DSf
	common /mssm_para/ MSf, MSf2, MSusy, MGl, MGl2
	common /mssm_para/ Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2
	common /mssm_para/ Mh02tree, MHH2tree, MHp2tree
	common /mssm_para/ CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	common /mssm_para/ CA, SA, CA2, SA2, C2A, S2A
	common /mssm_para/ CAB, SAB, CBA, SBA
	common /mssm_para/ mssm_digest

	double precision Af_flat(3*3*3), Xf_flat(3*3*3)
	equivalence (Af, Af_flat)
	equivalence (Xf, Xf_flat)

	double precision ReImAtau(2), ReAtau, ImAtau
	equivalence (Atau, ReImAtau)
	equivalence (ReImAtau(1), ReAtau), (ReImAtau(2), ImAtau)

	double precision ReImAt(2), ReAt, ImAt
	equivalence (At, ReImAt)
	equivalence (ReImAt(1), ReAt), (ReImAt(2), ImAt)

	double precision ReImAb(2), ReAb, ImAb
	equivalence (Ab, ReImAb)
	equivalence (ReImAb(1), ReAb), (ReImAb(2), ImAb)

	double precision ReImMUE(2), ReMUE, ImMUE
	equivalence (MUE, ReImMUE)
	equivalence (ReImMUE(1), ReMUE), (ReImMUE(2), ImMUE)

	double precision ReImM_1(2), ReM_1, ImM_1
	equivalence (M_1, ReImM_1)
	equivalence (ReImM_1(1), ReM_1), (ReImM_1(2), ImM_1)

	double precision ReImM_2(2), ReM_2, ImM_2
	equivalence (M_2, ReImM_2)
	equivalence (ReImM_2(1), ReM_2), (ReImM_2(2), ImM_2)

	double precision ReImM_3(2), ReM_3, ImM_3
	equivalence (M_3, ReImM_3)
	equivalence (ReImM_3(1), ReM_3), (ReImM_3(2), ImM_3)

#ifndef CKMC
#define CKMC(a,b) DCONJG(CKM(a,b))
#define USfC(a,b,t,g) DCONJG(USf(a,b,t,g))
#define VChaC(a,b) DCONJG(VCha(a,b))
#define UChaC(a,b) DCONJG(UCha(a,b))
#define ZNeuC(a,b) DCONJG(ZNeu(a,b))
#endif


* flavour-violating parameters

	double complex deltaSf(3:4,6,6)
	double complex UASf(6,6,3:4)
	double precision MASf(6,3:4), MASf2(6,3:4)

	common /fv_para/ UASf, MASf, MASf2, deltaSf

	double complex deltaSf_flat(2*6*6)
	equivalence (deltaSf, deltaSf_flat)

	double precision ReImdeltaSf(2,3:4,6,6)
	equivalence (deltaSf, ReImdeltaSf)

	double complex deltaLLuc, deltaLRuc
	double complex deltaRLucC, deltaRRuc
	equivalence (deltaSf(3,1  ,2  ), deltaLLuc)
	equivalence (deltaSf(3,1  ,2+3), deltaLRuc)
	equivalence (deltaSf(3,2  ,1+3), deltaRLucC)
	equivalence (deltaSf(3,1+3,2+3), deltaRRuc)

	double precision RedeltaLLuc, ImdeltaLLuc
	double precision RedeltaLRuc, ImdeltaLRuc
	double precision RedeltaRLucC, ImdeltaRLucC
	double precision RedeltaRRuc, ImdeltaRRuc
	equivalence (ReImdeltaSf(1,3,1  ,2  ), RedeltaLLuc)
	equivalence (ReImdeltaSf(2,3,1  ,2  ), ImdeltaLLuc)
	equivalence (ReImdeltaSf(1,3,1  ,2+3), RedeltaLRuc)
	equivalence (ReImdeltaSf(2,3,1  ,2+3), ImdeltaLRuc)
	equivalence (ReImdeltaSf(1,3,2  ,1+3), RedeltaRLucC)
	equivalence (ReImdeltaSf(2,3,2  ,1+3), ImdeltaRLucC)
	equivalence (ReImdeltaSf(1,3,1+3,2+3), RedeltaRRuc)
	equivalence (ReImdeltaSf(2,3,1+3,2+3), ImdeltaRRuc)

	double complex deltaLLct, deltaLRct
	double complex deltaRLctC, deltaRRct
	equivalence (deltaSf(3,2  ,3  ), deltaLLct)
	equivalence (deltaSf(3,2  ,3+3), deltaLRct)
	equivalence (deltaSf(3,3  ,2+3), deltaRLctC)
	equivalence (deltaSf(3,2+3,3+3), deltaRRct)

	double precision RedeltaLLct, ImdeltaLLct
	double precision RedeltaLRct, ImdeltaLRct
	double precision RedeltaRLctC, ImdeltaRLctC
	double precision RedeltaRRct, ImdeltaRRct
	equivalence (ReImdeltaSf(1,3,2  ,3  ), RedeltaLLct)
	equivalence (ReImdeltaSf(2,3,2  ,3  ), ImdeltaLLct)
	equivalence (ReImdeltaSf(1,3,2  ,3+3), RedeltaLRct)
	equivalence (ReImdeltaSf(2,3,2  ,3+3), ImdeltaLRct)
	equivalence (ReImdeltaSf(1,3,3  ,2+3), RedeltaRLctC)
	equivalence (ReImdeltaSf(2,3,3  ,2+3), ImdeltaRLctC)
	equivalence (ReImdeltaSf(1,3,2+3,3+3), RedeltaRRct)
	equivalence (ReImdeltaSf(2,3,2+3,3+3), ImdeltaRRct)

	double complex deltaLLut, deltaLRut
	double complex deltaRLutC, deltaRRut
	equivalence (deltaSf(3,1  ,3  ), deltaLLut)
	equivalence (deltaSf(3,1  ,3+3), deltaLRut)
	equivalence (deltaSf(3,3  ,1+3), deltaRLutC)
	equivalence (deltaSf(3,1+3,3+3), deltaRRut)

	double precision RedeltaLLut, ImdeltaLLut
	double precision RedeltaLRut, ImdeltaLRut
	double precision RedeltaRLutC, ImdeltaRLutC
	double precision RedeltaRRut, ImdeltaRRut
	equivalence (ReImdeltaSf(1,3,1  ,3  ), RedeltaLLut)
	equivalence (ReImdeltaSf(2,3,1  ,3  ), ImdeltaLLut)
	equivalence (ReImdeltaSf(1,3,1  ,3+3), RedeltaLRut)
	equivalence (ReImdeltaSf(2,3,1  ,3+3), ImdeltaLRut)
	equivalence (ReImdeltaSf(1,3,3  ,1+3), RedeltaRLutC)
	equivalence (ReImdeltaSf(2,3,3  ,1+3), ImdeltaRLutC)
	equivalence (ReImdeltaSf(1,3,1+3,3+3), RedeltaRRut)
	equivalence (ReImdeltaSf(2,3,1+3,3+3), ImdeltaRRut)

	double complex deltaLLds, deltaLRds
	double complex deltaRLdsC, deltaRRds
	equivalence (deltaSf(4,1  ,2  ), deltaLLds)
	equivalence (deltaSf(4,1  ,2+3), deltaLRds)
	equivalence (deltaSf(4,2  ,1+3), deltaRLdsC)
	equivalence (deltaSf(4,1+3,2+3), deltaRRds)

	double precision RedeltaLLds, ImdeltaLLds
	double precision RedeltaLRds, ImdeltaLRds
	double precision RedeltaRLdsC, ImdeltaRLdsC
	double precision RedeltaRRds, ImdeltaRRds
	equivalence (ReImdeltaSf(1,4,1  ,2  ), RedeltaLLds)
	equivalence (ReImdeltaSf(2,4,1  ,2  ), ImdeltaLLds)
	equivalence (ReImdeltaSf(1,4,1  ,2+3), RedeltaLRds)
	equivalence (ReImdeltaSf(2,4,1  ,2+3), ImdeltaLRds)
	equivalence (ReImdeltaSf(1,4,2  ,1+3), RedeltaRLdsC)
	equivalence (ReImdeltaSf(2,4,2  ,1+3), ImdeltaRLdsC)
	equivalence (ReImdeltaSf(1,4,1+3,2+3), RedeltaRRds)
	equivalence (ReImdeltaSf(2,4,1+3,2+3), ImdeltaRRds)

	double complex deltaLLsb, deltaLRsb
	double complex deltaRLsbC, deltaRRsb
	equivalence (deltaSf(4,2  ,3  ), deltaLLsb)
	equivalence (deltaSf(4,2  ,3+3), deltaLRsb)
	equivalence (deltaSf(4,3  ,2+3), deltaRLsbC)
	equivalence (deltaSf(4,2+3,3+3), deltaRRsb)

	double precision RedeltaLLsb, ImdeltaLLsb
	double precision RedeltaLRsb, ImdeltaLRsb
	double precision RedeltaRLsbC, ImdeltaRLsbC
	double precision RedeltaRRsb, ImdeltaRRsb
	equivalence (ReImdeltaSf(1,4,2  ,3  ), RedeltaLLsb)
	equivalence (ReImdeltaSf(2,4,2  ,3  ), ImdeltaLLsb)
	equivalence (ReImdeltaSf(1,4,2  ,3+3), RedeltaLRsb)
	equivalence (ReImdeltaSf(2,4,2  ,3+3), ImdeltaLRsb)
	equivalence (ReImdeltaSf(1,4,3  ,2+3), RedeltaRLsbC)
	equivalence (ReImdeltaSf(2,4,3  ,2+3), ImdeltaRLsbC)
	equivalence (ReImdeltaSf(1,4,2+3,3+3), RedeltaRRsb)
	equivalence (ReImdeltaSf(2,4,2+3,3+3), ImdeltaRRsb)

	double complex deltaLLdb, deltaLRdb
	double complex deltaRLdbC, deltaRRdb
	equivalence (deltaSf(4,1  ,3  ), deltaLLdb)
	equivalence (deltaSf(4,1  ,3+3), deltaLRdb)
	equivalence (deltaSf(4,3  ,1+3), deltaRLdbC)
	equivalence (deltaSf(4,1+3,3+3), deltaRRdb)

	double precision RedeltaLLdb, ImdeltaLLdb
	double precision RedeltaLRdb, ImdeltaLRdb
	double precision RedeltaRLdbC, ImdeltaRLdbC
	double precision RedeltaRRdb, ImdeltaRRdb
	equivalence (ReImdeltaSf(1,4,1  ,3  ), RedeltaLLdb)
	equivalence (ReImdeltaSf(2,4,1  ,3  ), ImdeltaLLdb)
	equivalence (ReImdeltaSf(1,4,1  ,3+3), RedeltaLRdb)
	equivalence (ReImdeltaSf(2,4,1  ,3+3), ImdeltaLRdb)
	equivalence (ReImdeltaSf(1,4,3  ,1+3), RedeltaRLdbC)
	equivalence (ReImdeltaSf(2,4,3  ,1+3), ImdeltaRLdbC)
	equivalence (ReImdeltaSf(1,4,1+3,3+3), RedeltaRRdb)
	equivalence (ReImdeltaSf(2,4,1+3,3+3), ImdeltaRRdb)

#ifndef UASfC
#define UASfC(a,b,t) DCONJG(UASf(a,b,t))
#endif


* THDM parameters

	double precision Lambda5
	double precision Yuk1, Yuk2, Yuk3
	logical thdm_digest

	common /thdm_para/ Lambda5
	common /thdm_para/ Yuk1, Yuk2, Yuk3
	common /thdm_para/ thdm_digest

