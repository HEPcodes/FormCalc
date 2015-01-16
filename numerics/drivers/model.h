* model.h
* common blocks for the model parameters
* this file is part of FormCalc
* last modified 3 Apr 00 th


	double precision pi
	parameter (pi = 3.1415926535897932384626433832795029D0)

* SM parameters

	double precision EL, Alfa, Alfa2, Gmue, GS, Alfas, Alfas2
	double precision MW, MW2, MZ, MZ2, MH, MH2
	double precision SW, SW2, CW, CW2
	double precision ME, ME2, MM, MM2, ML, ML2, MLE(3), MLE2(3)
	double precision MU, MU2, MC, MC2, MT, MT2, MQU(3), MQU2(3)
	double precision MD, MD2, MS, MS2, MB, MB2, MQD(3), MQD2(3)

	parameter (Alfa = 1/137.0359895D0, Alfa2 = Alfa**2)
	parameter (Gmue = 1.16639D-5)

	parameter (MW = 80.39D0, MW2 = MW**2)
	parameter (MZ = 91.1867D0, MZ2 = MZ**2)
	parameter (SW2 = (MZ2 - MW2)/MZ2)
	parameter (CW = MW/MZ, CW2 = CW**2)

	parameter (ME = .51099907D-3, ME2 = ME**2)
	parameter (MM = 105.658389D-3, MM2 = MM**2)
	parameter (ML = 1777D-3, ML2 = ML**2)

	parameter (MU = 46.4D-3, MU2 = MU**2)
	parameter (MC = 1.50D0, MC2 = MC**2)

	parameter (MD = 46.5D-3, MD2 = MD**2)
	parameter (MS = 150D-3, MS2 = MS**2)
	parameter (MB = 4.7D0, MB2 = MB**2)

	common /sm_para/
     +    MLE, MLE2, MQU, MQU2, MQD, MQD2,
     +    EL, GS, Alfas, Alfas2, SW, MT, MT2, MH, MH2


* MSSM parameters

	double complex USf(2, 2, 4, 3)
	double complex UCha(2, 2), VCha(2, 2), ZNeu(4, 4)
	double complex Af(4, 3), MUE
	double precision M_1, M_2
	double precision MSf(2, 4, 3), MSf2(2, 4, 3)
	double precision MCha(2), MNeu(4), MCha2(2), MNeu2(4)
	double precision MSNE(3), MSLE1(3), MSQU1(3), MSQD1(3)
	double precision MSLE2(3), MSQU2(3), MSQD2(3)
	double precision Mh0, MHH, MA0, MG0, MHp, MGp, MGl, MSusy
	double precision Mh02, MHH2, MA02, MG02, MHp2, MGp2, MGl2
	double precision CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	double precision CA, SA, CA2, SA2, C2A, S2A
	double precision CAB, SAB, CBA, SBA

	common /mssm_para/
     +    USf, UCha, VCha, ZNeu,
     +    Af, MUE,
     +    M_1, M_2,
     +    MSf, MSf2,
     +    MCha, MNeu, MCha2, MNeu2,
     +    MSNE, MSLE1, MSQU1, MSQD1,
     +    MSLE2, MSQU2, MSQD2,
     +    Mh0, MHH, MA0, MG0, MHp, MGp, MGl, MSusy,
     +    Mh02, MHH2, MA02, MG02, MHp2, MGp2, MGl2,
     +    CB, SB, TB, CB2, SB2, TB2, C2B, S2B,
     +    CA, SA, CA2, SA2, C2A, S2A,
     +    CAB, SAB, CBA, SBA

#ifndef USfC
#define USfC(a, b, t, g) dconjg(USf(a, b, t, g))
#define VChaC(a, b) dconjg(VCha(a, b))
#define UChaC(a, b) dconjg(UCha(a, b))
#define ZNeuC(a, b) dconjg(ZNeu(a, b))
#endif

