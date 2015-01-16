* model.h
* common blocks for the model parameters
* this file is part of FormCalc
* last modified 16 Dec 04 th


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

#ifndef polar
#define polar(r, theta) r*exp(cI*degree*theta)
#endif

	double precision Divergence
	common /renorm/ Divergence


* SM parameters

	double precision MZ, MZ2, MW, MW2, CW, CW2, SW2
	parameter (MZ = 91.1875D0, MZ2 = MZ**2)
	parameter (MW = 80.450D0, MW2 = MW**2)
	parameter (CW = MW/MZ, CW2 = CW**2)
	parameter (SW2 = 1 - CW2)

	double precision GF, Alfa, Alfa2
	parameter (GF = 1.16639D-5)
	parameter (Alfa = 1/137.0359895D0, Alfa2 = Alfa**2)
c	parameter (Alfa = sqrt2/pi*GF*MW2*SW2, Alfa2 = Alfa**2)

	double precision ME, ME2, MM, MM2, ML, ML2
	parameter (ME = .51099907D-3, ME2 = ME**2)
	parameter (MM = 105.658389D-3, MM2 = MM**2)
	parameter (ML = 1777D-3, ML2 = ML**2)

	double precision MU, MU2, MC, MC2, MT, MT2
	parameter (MU = 53.8D-3, MU2 = MU**2)
	parameter (MC = 1.50D0, MC2 = MC**2)
	parameter (MT = 178D0, MT2 = MT**2)

	double precision MD, MD2, MS, MS2, MB, MB2
	parameter (MD = 53.8D-3, MD2 = MD**2)
	parameter (MS = 150D-3, MS2 = MS**2)
	parameter (MB = 4.7D0, MB2 = MB**2)

	double complex CKM(3,3)
	double precision Mf(4,3), Mf2(4,3)
	double precision MH, MH2, MG0, MG02, MGp, MGp2
	double precision EL, GS, Alfas, Alfas2, AlfasMT, SW

	common /sm_para/
     &    CKM, Mf, Mf2, MH, MH2, MG0, MG02, MGp, MGp2,
     &    EL, GS, Alfas, Alfas2, AlfasMT, SW


* MSSM parameters

	double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
	double complex USf(2,2,4,3), Af(2:4,3), Xf(2:4,3)
	double complex Atau, At, Ab, MUE
	double precision MCha(2), MCha2(2), MNeu(4), MNeu2(4)
	double precision MSS(2,2:4,3), MSS2(2,2:4,3), DSf(2,4)
	double precision MSf(2,4,3), MSf2(2,4,3), MSusy
	double precision Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2
	double precision Mh02tree, MHH2tree, MHp2tree
	double precision M_1, M_2, MGl, MGl2
	double precision CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	double precision CA, SA, CA2, SA2, C2A, S2A
	double precision CAB, SAB, CBA, SBA

	common /mssm_para/
     &    UCha, VCha, ZNeu,
     &    USf, Af, Xf, Atau, At, Ab, MUE,
     &    MCha, MCha2, MNeu, MNeu2,
     &    MSS, MSS2, DSf, MSf, MSf2, MSusy,
     &    Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2,
     &    Mh02tree, MHH2tree, MHp2tree,
     &    M_1, M_2, MGl, MGl2,
     &    CB, SB, TB, CB2, SB2, TB2, C2B, S2B,
     &    CA, SA, CA2, SA2, C2A, S2A,
     &    CAB, SAB, CBA, SBA

	double precision reimMUE(2), reMUE, imMUE
	equivalence (MUE, reimMUE)
	equivalence (reimMUE(1), reMUE), (reimMUE(2), imMUE)

	double precision reimAtau(2), reAtau, imAtau
	equivalence (Atau, reimAtau)
	equivalence (reimAtau(1), reAtau), (reimAtau(2), imAtau)

	double precision reimAt(2), reAt, imAt
	equivalence (At, reimAt)
	equivalence (reimAt(1), reAt), (reimAt(2), imAt)

	double precision reimAb(2), reAb, imAb
	equivalence (Ab, reimAb)
	equivalence (reimAb(1), reAb), (reimAb(2), imAb)

#ifndef CKMC
#define CKMC(a,b) dconjg(CKM(a,b))
#define USfC(a,b,t,g) dconjg(USf(a,b,t,g))
#define VChaC(a,b) dconjg(VCha(a,b))
#define UChaC(a,b) dconjg(UCha(a,b))
#define ZNeuC(a,b) dconjg(ZNeu(a,b))
#endif


* THDM parameters

	double precision Lambda5
	double precision Yuk1, Yuk2, Yuk3

	common /thdm_para/
     &    Lambda5,
     &    Yuk1, Yuk2, Yuk3

