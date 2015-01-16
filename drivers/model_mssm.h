* model_mssm.h
* declarations for model_mssm.F
* this file is part of FormCalc
* last modified 18 May 11 th


#include "model_sm.h"

	double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
	double complex USf(2,2,4,3), UCSf(3,4,2:4,3), UUSf(3,4,2:4,3)
	double complex XHiggs(3,3,2)
	double complex Af(2:4,3,3), Xf(2:4,3), MUETB(2:4)
	double complex Atau, At, Ab, MUE
	double complex Mino1, Mino2, Mino3, SqrtEGl
	double precision MCha(2), MCha2(2), MNeu(4), MNeu2(4)
	double precision MSS(2,2:4,3), MSS2(2,2:4,3), DSf(2,4)
	double precision MSf(2,4,3), MSf2(2,4,3), MSusy, MGl, MGl2
	double precision MHiggs(4), MHiggs2(4), MHiggstree2(4)
	double precision CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	double precision CA, SA, CA2, SA2, C2A, S2A
	double precision CAB, SAB, CBA, SBA, CBA2, SBA2
	double precision AlfasMT

	common /mssm_para/ UCha, VCha, ZNeu
	common /mssm_para/ USf, UCSf, UUSf
	common /mssm_para/ XHiggs
	common /mssm_para/ Af, Xf
	common /mssm_para/ Atau, At, Ab, MUE, MUETB
	common /mssm_para/ Mino1, Mino2, Mino3, SqrtEGl
	common /mssm_para/ MCha, MCha2, MNeu, MNeu2
	common /mssm_para/ MSS, MSS2, DSf
	common /mssm_para/ MSf, MSf2, MSusy, MGl, MGl2
	common /mssm_para/ MHiggs, MHiggs2, MHiggstree2
	common /mssm_para/ CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	common /mssm_para/ CA, SA, CA2, SA2, C2A, S2A
	common /mssm_para/ CAB, SAB, CBA, SBA, CBA2, SBA2
	common /mssm_para/ AlfasMT

#ifndef USfC
#define USfC(i,j,t,g) DCONJG(USf(i,j,t,g))
#define UCSfC(i,j,t,g) DCONJG(UCSf(i,j,t,g))
#define UUSfC(i,j,t,g) DCONJG(UUSf(i,j,t,g))
#define USf2(i,j,t,g) DBLE(UCSf(i,j,t,g))
#define VChaC(i,j) DCONJG(VCha(i,j))
#define UChaC(i,j) DCONJG(UCha(i,j))
#define ZNeuC(i,j) DCONJG(ZNeu(i,j))
#define UHiggsC(i,j) DCONJG(UHiggs(i,j))
#define ZHiggsC(i,j) DCONJG(ZHiggs(i,j))
#define AfC(t,g1,g2) DCONJG(Af(t,g1,g2))
#define Mino3C DCONJG(Mino3)
#define MUEC DCONJG(MUE)
#define SqrtEGlC DCONJG(SqrtEGl)
#endif

	double precision Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2
	equivalence (MHiggs(1), Mh0), (MHiggs2(1), Mh02)
	equivalence (MHiggs(2), MHH), (MHiggs2(2), MHH2)
	equivalence (MHiggs(3), MA0), (MHiggs2(3), MA02)
	equivalence (MHiggs(4), MHp), (MHiggs2(4), MHp2)

	double precision Mh0tree2, MHHtree2, MA0tree2, MHptree2
	equivalence (MHiggstree2(1), Mh0tree2)
	equivalence (MHiggstree2(2), MHHtree2)
	equivalence (MHiggstree2(3), MA0tree2)
	equivalence (MHiggstree2(4), MHptree2)

	double complex UHiggs(3,3), ZHiggs(3,3)
	equivalence (XHiggs(1,1,1), UHiggs)
	equivalence (XHiggs(1,1,2), ZHiggs)

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

	double precision ReImMino1(2), ReMino1, ImMino1
	double precision M_1, ReM_1, ImM_1
	equivalence (Mino1, M_1, ReImMino1)
	equivalence (ReImMino1(1), ReMino1, ReM_1)
	equivalence (ReImMino1(2), ImMino1, ImM_1)

	double precision ReImMino2(2), ReMino2, ImMino2
	double precision M_2, ReM_2, ImM_2
	equivalence (Mino2, M_2, ReImMino2)
	equivalence (ReImMino2(1), ReMino2, ReM_2)
	equivalence (ReImMino2(2), ImMino2, ImM_2)

	double precision ReImMino3(2), ReMino3, ImMino3
	double precision M_3, ReM_3, ImM_3
	equivalence (Mino3, M_3, ReImMino3)
	equivalence (ReImMino3(1), ReMino3, ReM_3)
	equivalence (ReImMino3(2), ImMino3, ImM_3)


* flavour-violating parameters

	double complex deltaSf(6,6,3:4)
	double complex UASf(6,6,3:4)
	double precision MASf(6,3:4), MASf2(6,3:4)

	common /fv_para/ UASf, MASf, MASf2, deltaSf

#ifndef UASfC
#define UASfC(i,j,t) DCONJG(UASf(i,j,t))

#define deltaSf_LL(i,j,t) deltaSf(i,j,t)
#define deltaSf_LR(i,j,t) deltaSf(i,j+3,t)
#define deltaSf_RL(i,j,t) deltaSf(j,i+3,t)
#define deltaSf_RR(i,j,t) deltaSf(i+3,j+3,t)

#define ImdeltaSf_LL(i,j,t) ReImdeltaSf(2,i,j,t)
#define ImdeltaSf_LR(i,j,t) ReImdeltaSf(2,i,j+3,t)
#define ImdeltaSf_RL(i,j,t) ReImdeltaSf(2,j,i+3,t)
#define ImdeltaSf_RR(i,j,t) ReImdeltaSf(2,i+3,j+3,t)
#endif

	double complex deltaSf_flat(6*6*2)
	equivalence (deltaSf, deltaSf_flat)

	double precision ReImdeltaSf(2,6,6,3:4)
	equivalence (deltaSf, ReImdeltaSf)

	double complex deltaLRuc
	double precision RedeltaLRuc, ImdeltaLRuc
	equivalence (deltaSf_LR(1,2,3), deltaLRuc, RedeltaLRuc)
	equivalence (ImdeltaSf_LR(1,2,3), ImdeltaLRuc)
	double complex deltaLRct
	double precision RedeltaLRct, ImdeltaLRct
	equivalence (deltaSf_LR(2,3,3), deltaLRct, RedeltaLRct)
	equivalence (ImdeltaSf_LR(2,3,3), ImdeltaLRct)
	double complex deltaLRut
	double precision RedeltaLRut, ImdeltaLRut
	equivalence (deltaSf_LR(1,3,3), deltaLRut, RedeltaLRut)
	equivalence (ImdeltaSf_LR(1,3,3), ImdeltaLRut)

	double complex deltaRLuc
	double precision RedeltaRLuc, ImdeltaRLuc
	equivalence (deltaSf_RL(1,2,3), deltaRLuc, RedeltaRLuc)
	equivalence (ImdeltaSf_RL(1,2,3), ImdeltaRLuc)
	double complex deltaRLct
	double precision RedeltaRLct, ImdeltaRLct
	equivalence (deltaSf_RL(2,3,3), deltaRLct, RedeltaRLct)
	equivalence (ImdeltaSf_RL(2,3,3), ImdeltaRLct)
	double complex deltaRLut
	double precision RedeltaRLut, ImdeltaRLut
	equivalence (deltaSf_RL(1,3,3), deltaRLut, RedeltaRLut)
	equivalence (ImdeltaSf_RL(1,3,3), ImdeltaRLut)

	double complex deltaRRuc
	double precision RedeltaRRuc, ImdeltaRRuc
	equivalence (deltaSf_RR(1,2,3), deltaRRuc, RedeltaRRuc)
	equivalence (ImdeltaSf_RR(1,2,3), ImdeltaRRuc)
	double complex deltaRRct
	double precision RedeltaRRct, ImdeltaRRct
	equivalence (deltaSf_RR(2,3,3), deltaRRct, RedeltaRRct)
	equivalence (ImdeltaSf_RR(2,3,3), ImdeltaRRct)
	double complex deltaRRut
	double precision RedeltaRRut, ImdeltaRRut
	equivalence (deltaSf_RR(1,3,3), deltaRRut, RedeltaRRut)
	equivalence (ImdeltaSf_RR(1,3,3), ImdeltaRRut)

	double complex deltaLL12
	double precision RedeltaLL12, ImdeltaLL12
	equivalence (deltaSf_LL(1,2,4), deltaLL12, RedeltaLL12)
	equivalence (ImdeltaSf_LL(1,2,4), ImdeltaLL12)
	double complex deltaLL23
	double precision RedeltaLL23, ImdeltaLL23
	equivalence (deltaSf_LL(2,3,4), deltaLL23, RedeltaLL23)
	equivalence (ImdeltaSf_LL(2,3,4), ImdeltaLL23)
	double complex deltaLL13
	double precision RedeltaLL13, ImdeltaLL13
	equivalence (deltaSf_LL(1,3,4), deltaLL13, RedeltaLL13)
	equivalence (ImdeltaSf_LL(1,3,4), ImdeltaLL13)

	double complex deltaLRds
	double precision RedeltaLRds, ImdeltaLRds
	equivalence (deltaSf_LR(1,2,4), deltaLRds, RedeltaLRds)
	equivalence (ImdeltaSf_LR(1,2,4), ImdeltaLRds)
	double complex deltaLRsb
	double precision RedeltaLRsb, ImdeltaLRsb
	equivalence (deltaSf_LR(2,3,4), deltaLRsb, RedeltaLRsb)
	equivalence (ImdeltaSf_LR(2,3,4), ImdeltaLRsb)
	double complex deltaLRdb
	double precision RedeltaLRdb, ImdeltaLRdb
	equivalence (deltaSf_LR(1,3,4), deltaLRdb, RedeltaLRdb)
	equivalence (ImdeltaSf_LR(1,3,4), ImdeltaLRdb)

	double complex deltaRLds
	double precision RedeltaRLds, ImdeltaRLds
	equivalence (deltaSf_RL(1,2,4), deltaRLds, RedeltaRLds)
	equivalence (ImdeltaSf_RL(1,2,4), ImdeltaRLds)
	double complex deltaRLsb
	double precision RedeltaRLsb, ImdeltaRLsb
	equivalence (deltaSf_RL(2,3,4), deltaRLsb, RedeltaRLsb)
	equivalence (ImdeltaSf_RL(2,3,4), ImdeltaRLsb)
	double complex deltaRLdb
	double precision RedeltaRLdb, ImdeltaRLdb
	equivalence (deltaSf_RL(1,3,4), deltaRLdb, RedeltaRLdb)
	equivalence (ImdeltaSf_RL(1,3,4), ImdeltaRLdb)

	double complex deltaRRds
	double precision RedeltaRRds, ImdeltaRRds
	equivalence (deltaSf_RR(1,2,4), deltaRRds, RedeltaRRds)
	equivalence (ImdeltaSf_RR(1,2,4), ImdeltaRRds)
	double complex deltaRRsb
	double precision RedeltaRRsb, ImdeltaRRsb
	equivalence (deltaSf_RR(2,3,4), deltaRRsb, RedeltaRRsb)
	equivalence (ImdeltaSf_RR(2,3,4), ImdeltaRRsb)
	double complex deltaRRdb
	double precision RedeltaRRdb, ImdeltaRRdb
	equivalence (deltaSf_RR(1,3,4), deltaRRdb, RedeltaRRdb)
	equivalence (ImdeltaSf_RR(1,3,4), ImdeltaRRdb)
