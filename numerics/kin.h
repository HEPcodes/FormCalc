*	kin.h
*	common blocks with kinematical variables for num.F
*	(yeah, and some other variables, too, so we do not need so
*	many include files)
*	last modified 21 Jun 98 th

	double precision mass(4)
	integer type(4)
	common /process/ mass, type

	double complex resoc
#ifdef DYSON
	double complex reso, resoT, resoU
#else
	double precision reso, resoT, resoU
#endif
	double precision EE(4), Ecms, Pin, Pout, Pin2, Pout2
	double precision th, st, ct, S, T, U, kinf
	double precision bornamp, loopamp, bornsum, loopsum
	double precision bornAfb, loopAfb
	integer Cptr2, Dptr2
	common /kin/ 
     +    resoc, reso, resoT, resoU,
     +    EE, Ecms, Pin, Pout, Pin2, Pout2,
     +    th, st, ct, S, T, U, kinf,
     +    bornamp, loopamp, bornsum, loopsum,
     +    bornAfb, loopAfb,
     +    Cptr2, Dptr2

	double complex v4(4, 16)
	common /vec/ v4

	double precision avgfac
	integer bpol(4), epol(4)
	integer i1, i2, i3, i4
	character polstr*4
	common /pol/
     +    avgfac,
     +    bpol, epol,
     +    i1, i2, i3, i4,
     +    polstr

	double precision pi, ds2, degree, hbarc2
	parameter (pi = 3.1415926535897932384626433832795028841972D0)
	parameter (ds2 = .7071067811865475244008443621048490392848D0)
*	  = 1/sqrt(2)
	parameter (degree = pi/180D0)
	parameter (hbarc2 = 3.8937966D8)
*	  = \hbar c^2 in picobarn


* SM parameters
	double precision MLE(3), MQU(3), MQD(3)
	double precision MLE2(3), MQU2(3), MQD2(3)
	double precision MH, MH2, MT, MT2
	common /sm/ MLE, MQU, MQD,
     +    MLE2, MQU2, MQD2,
     +    MH, MH2, MT, MT2

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
	double precision CB, SB, TB, C2B, S2B
	double precision CA, SA, C2A, S2A
	double precision CAB, SAB, CBA, SBA
	common /mssm/ USf, UCha, VCha, ZNeu,
     +    Af, MUE,
     +    M_1, M_2,
     +    MSf, MSf2,
     +    MCha, MNeu, MCha2, MNeu2,
     +    MSNE, MSLE1, MSQU1, MSQD1,
     +    MSLE2, MSQU2, MSQD2,
     +    Mh0, MHH, MA0, MG0, MHp, MGp, MGl, MSusy,
     +    CB, SB, TB, C2B, S2B,
     +    CA, SA, C2A, S2A,
     +    CAB, SAB, CBA, SBA

#ifndef USfC
#define USfC(a, b, t, g) dconjg(USf(a, b, t, g))
#define VChaC(a, b) dconjg(VCha(a, b))
#define UChaC(a, b) dconjg(UCha(a, b))
#define ZNeuC(a, b) dconjg(ZNeu(a, b))
#endif

