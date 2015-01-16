* model_mssm.h
* declarations for model_mssm.F
* this file is part of FormCalc
* last modified 8 Mar 13 th


#include "model_sm.h"

	ComplexType UCha(2,2), VCha(2,2), ZNeu(4,4)
	ComplexType USf(2,2,4,3), UCSf(3,4,2:4,3), UUSf(3,4,2:4,3)
	ComplexType XHiggs(3,3,2)
	ComplexType Af(2:4,3,3), Xf(2:4,3), MUETB(2:4)
	ComplexType Atau, At, Ab, MUE
	ComplexType Mino1, Mino2, Mino3, SqrtEGl
	RealType MCha(2), MCha2(2), MNeu(4), MNeu2(4)
	RealType MSS(2,2:4,3), MSS2(2,2:4,3), DSf(2,4)
	RealType MSf(2,4,3), MSf2(2,4,3), MSusy, MGl, MGl2
	RealType MHiggs(4), MHiggs2(4), MHiggstree2(4)
	RealType CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	RealType CA, SA, CA2, SA2, C2A, S2A
	RealType CAB, SAB, CBA, SBA, CBA2, SBA2
	RealType AlfasMT

	common /mssmpara/ UCha, VCha, ZNeu
	common /mssmpara/ USf, UCSf, UUSf
	common /mssmpara/ XHiggs
	common /mssmpara/ Af, Xf
	common /mssmpara/ Atau, At, Ab, MUE, MUETB
	common /mssmpara/ Mino1, Mino2, Mino3, SqrtEGl
	common /mssmpara/ MCha, MCha2, MNeu, MNeu2
	common /mssmpara/ MSS, MSS2, DSf
	common /mssmpara/ MSf, MSf2, MSusy, MGl, MGl2
	common /mssmpara/ MHiggs, MHiggs2, MHiggstree2
	common /mssmpara/ CB, SB, TB, CB2, SB2, TB2, C2B, S2B
	common /mssmpara/ CA, SA, CA2, SA2, C2A, S2A
	common /mssmpara/ CAB, SAB, CBA, SBA, CBA2, SBA2
	common /mssmpara/ AlfasMT

#ifndef USfC
#define USfC(i,j,t,g) Conjugate(USf(i,j,t,g))
#define UCSfC(i,j,t,g) Conjugate(UCSf(i,j,t,g))
#define UUSfC(i,j,t,g) Conjugate(UUSf(i,j,t,g))
#define USf2(i,j,t,g) Re(UCSf(i,j,t,g))
#define VChaC(i,j) Conjugate(VCha(i,j))
#define UChaC(i,j) Conjugate(UCha(i,j))
#define ZNeuC(i,j) Conjugate(ZNeu(i,j))
#define UHiggsC(i,j) Conjugate(UHiggs(i,j))
#define ZHiggsC(i,j) Conjugate(ZHiggs(i,j))
#define AfC(t,g1,g2) Conjugate(Af(t,g1,g2))
#define Mino3C Conjugate(Mino3)
#define MUEC Conjugate(MUE)
#define SqrtEGlC Conjugate(SqrtEGl)
#endif

	RealType Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2
	equivalence (MHiggs(1), Mh0), (MHiggs2(1), Mh02)
	equivalence (MHiggs(2), MHH), (MHiggs2(2), MHH2)
	equivalence (MHiggs(3), MA0), (MHiggs2(3), MA02)
	equivalence (MHiggs(4), MHp), (MHiggs2(4), MHp2)

	RealType Mh0tree2, MHHtree2, MA0tree2, MHptree2
	equivalence (MHiggstree2(1), Mh0tree2)
	equivalence (MHiggstree2(2), MHHtree2)
	equivalence (MHiggstree2(3), MA0tree2)
	equivalence (MHiggstree2(4), MHptree2)

	ComplexType UHiggs(3,3), ZHiggs(3,3)
	equivalence (XHiggs(1,1,1), UHiggs)
	equivalence (XHiggs(1,1,2), ZHiggs)

	RealType Af_flat(3*3*3), Xf_flat(3*3*3)
	equivalence (Af, Af_flat)
	equivalence (Xf, Xf_flat)

	RealType ReImAtau(2), ReAtau, ImAtau
	equivalence (Atau, ReImAtau)
	equivalence (ReImAtau(1), ReAtau), (ReImAtau(2), ImAtau)

	RealType ReImAt(2), ReAt, ImAt
	equivalence (At, ReImAt)
	equivalence (ReImAt(1), ReAt), (ReImAt(2), ImAt)

	RealType ReImAb(2), ReAb, ImAb
	equivalence (Ab, ReImAb)
	equivalence (ReImAb(1), ReAb), (ReImAb(2), ImAb)

	RealType ReImMUE(2), ReMUE, ImMUE
	equivalence (MUE, ReImMUE)
	equivalence (ReImMUE(1), ReMUE), (ReImMUE(2), ImMUE)

	RealType ReImMino1(2), ReMino1, ImMino1
	RealType M_1, ReM_1, ImM_1
	equivalence (Mino1, M_1, ReImMino1)
	equivalence (ReImMino1(1), ReMino1, ReM_1)
	equivalence (ReImMino1(2), ImMino1, ImM_1)

	RealType ReImMino2(2), ReMino2, ImMino2
	RealType M_2, ReM_2, ImM_2
	equivalence (Mino2, M_2, ReImMino2)
	equivalence (ReImMino2(1), ReMino2, ReM_2)
	equivalence (ReImMino2(2), ImMino2, ImM_2)

	RealType ReImMino3(2), ReMino3, ImMino3
	RealType M_3, ReM_3, ImM_3
	equivalence (Mino3, M_3, ReImMino3)
	equivalence (ReImMino3(1), ReMino3, ReM_3)
	equivalence (ReImMino3(2), ImMino3, ImM_3)


* flavour-violating parameters

	ComplexType deltaSf(6,6,2:4)
	ComplexType UASf(6,6,2:4)
	RealType MASf(6,2:4), MASf2(6,2:4)

	common /mssmpara/ UASf, MASf, MASf2, deltaSf

#ifndef UASfC
#define UASfC(i,j,t) Conjugate(UASf(i,j,t))

#define deltaSf_LL(i,j,t) deltaSf(i,j,t)
#define deltaSf_LR(i,j,t) deltaSf(i,j+3,t)
#define deltaSf_RL(i,j,t) deltaSf(j,i+3,t)
#define deltaSf_RR(i,j,t) deltaSf(i+3,j+3,t)

#define ImdeltaSf_LL(i,j,t) ReImdeltaSf(2,i,j,t)
#define ImdeltaSf_LR(i,j,t) ReImdeltaSf(2,i,j+3,t)
#define ImdeltaSf_RL(i,j,t) ReImdeltaSf(2,j,i+3,t)
#define ImdeltaSf_RR(i,j,t) ReImdeltaSf(2,i+3,j+3,t)
#endif

	ComplexType deltaSf_flat(6*6*2)
	equivalence (deltaSf, deltaSf_flat)

	RealType ReImdeltaSf(2,6,6,3:4)
	equivalence (deltaSf, ReImdeltaSf)

	ComplexType deltaULR12
	RealType RedeltaULR12, ImdeltaULR12
	equivalence (deltaSf_LR(1,2,3), deltaULR12, RedeltaULR12)
	equivalence (ImdeltaSf_LR(1,2,3), ImdeltaULR12)
	ComplexType deltaULR23
	RealType RedeltaULR23, ImdeltaULR23
	equivalence (deltaSf_LR(2,3,3), deltaULR23, RedeltaULR23)
	equivalence (ImdeltaSf_LR(2,3,3), ImdeltaULR23)
	ComplexType deltaULR13
	RealType RedeltaULR13, ImdeltaULR13
	equivalence (deltaSf_LR(1,3,3), deltaULR13, RedeltaULR13)
	equivalence (ImdeltaSf_LR(1,3,3), ImdeltaULR13)

	ComplexType deltaURL12
	RealType RedeltaURL12, ImdeltaURL12
	equivalence (deltaSf_RL(1,2,3), deltaURL12, RedeltaURL12)
	equivalence (ImdeltaSf_RL(1,2,3), ImdeltaURL12)
	ComplexType deltaURL23
	RealType RedeltaURL23, ImdeltaURL23
	equivalence (deltaSf_RL(2,3,3), deltaURL23, RedeltaURL23)
	equivalence (ImdeltaSf_RL(2,3,3), ImdeltaURL23)
	ComplexType deltaURL13
	RealType RedeltaURL13, ImdeltaURL13
	equivalence (deltaSf_RL(1,3,3), deltaURL13, RedeltaURL13)
	equivalence (ImdeltaSf_RL(1,3,3), ImdeltaURL13)

	ComplexType deltaURR12
	RealType RedeltaURR12, ImdeltaURR12
	equivalence (deltaSf_RR(1,2,3), deltaURR12, RedeltaURR12)
	equivalence (ImdeltaSf_RR(1,2,3), ImdeltaURR12)
	ComplexType deltaURR23
	RealType RedeltaURR23, ImdeltaURR23
	equivalence (deltaSf_RR(2,3,3), deltaURR23, RedeltaURR23)
	equivalence (ImdeltaSf_RR(2,3,3), ImdeltaURR23)
	ComplexType deltaURR13
	RealType RedeltaURR13, ImdeltaURR13
	equivalence (deltaSf_RR(1,3,3), deltaURR13, RedeltaURR13)
	equivalence (ImdeltaSf_RR(1,3,3), ImdeltaURR13)

	ComplexType deltaQLL12
	RealType RedeltaQLL12, ImdeltaQLL12
	equivalence (deltaSf_LL(1,2,4), deltaQLL12, RedeltaQLL12)
	equivalence (ImdeltaSf_LL(1,2,4), ImdeltaQLL12)
	ComplexType deltaQLL23
	RealType RedeltaQLL23, ImdeltaQLL23
	equivalence (deltaSf_LL(2,3,4), deltaQLL23, RedeltaQLL23)
	equivalence (ImdeltaSf_LL(2,3,4), ImdeltaQLL23)
	ComplexType deltaQLL13
	RealType RedeltaQLL13, ImdeltaQLL13
	equivalence (deltaSf_LL(1,3,4), deltaQLL13, RedeltaQLL13)
	equivalence (ImdeltaSf_LL(1,3,4), ImdeltaQLL13)

	ComplexType deltaDLR12
	RealType RedeltaDLR12, ImdeltaDLR12
	equivalence (deltaSf_LR(1,2,4), deltaDLR12, RedeltaDLR12)
	equivalence (ImdeltaSf_LR(1,2,4), ImdeltaDLR12)
	ComplexType deltaDLR23
	RealType RedeltaDLR23, ImdeltaDLR23
	equivalence (deltaSf_LR(2,3,4), deltaDLR23, RedeltaDLR23)
	equivalence (ImdeltaSf_LR(2,3,4), ImdeltaDLR23)
	ComplexType deltaDLR13
	RealType RedeltaDLR13, ImdeltaDLR13
	equivalence (deltaSf_LR(1,3,4), deltaDLR13, RedeltaDLR13)
	equivalence (ImdeltaSf_LR(1,3,4), ImdeltaDLR13)

	ComplexType deltaDRL12
	RealType RedeltaDRL12, ImdeltaDRL12
	equivalence (deltaSf_RL(1,2,4), deltaDRL12, RedeltaDRL12)
	equivalence (ImdeltaSf_RL(1,2,4), ImdeltaDRL12)
	ComplexType deltaDRL23
	RealType RedeltaDRL23, ImdeltaDRL23
	equivalence (deltaSf_RL(2,3,4), deltaDRL23, RedeltaDRL23)
	equivalence (ImdeltaSf_RL(2,3,4), ImdeltaDRL23)
	ComplexType deltaDRL13
	RealType RedeltaDRL13, ImdeltaDRL13
	equivalence (deltaSf_RL(1,3,4), deltaDRL13, RedeltaDRL13)
	equivalence (ImdeltaSf_RL(1,3,4), ImdeltaDRL13)

	ComplexType deltaDRR12
	RealType RedeltaDRR12, ImdeltaDRR12
	equivalence (deltaSf_RR(1,2,4), deltaDRR12, RedeltaDRR12)
	equivalence (ImdeltaSf_RR(1,2,4), ImdeltaDRR12)
	ComplexType deltaDRR23
	RealType RedeltaDRR23, ImdeltaDRR23
	equivalence (deltaSf_RR(2,3,4), deltaDRR23, RedeltaDRR23)
	equivalence (ImdeltaSf_RR(2,3,4), ImdeltaDRR23)
	ComplexType deltaDRR13
	RealType RedeltaDRR13, ImdeltaDRR13
	equivalence (deltaSf_RR(1,3,4), deltaDRR13, RedeltaDRR13)
	equivalence (ImdeltaSf_RR(1,3,4), ImdeltaDRR13)
