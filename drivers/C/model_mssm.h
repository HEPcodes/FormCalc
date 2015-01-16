#if 0
	model_mssm.h
	declarations for model_mssm.F
	this file is part of FormCalc
	last modified 9 Mar 13 th
#endif


#include "model_sm.h"

#if NOUNDERSCORE
#define mssmpara_ mssmpara
#define fvpara_ fvpara
#endif

struct mssmpara_ {
  ComplexType UCha[2][2], VCha[2][2], ZNeu[4][4];
  ComplexType USf[3][4][2][2], UCSf[3][4-1][4][3], UUSf[3][4-1][4][3];
  ComplexType XHiggs[2][3][3];
  ComplexType Af[3][3][4-1], Xf[3][4-1], MUETB[4-1];
  ComplexType Atau, At, Ab, MUE;
  ComplexType Mino1, Mino2, Mino3, SqrtEGl;
  RealType MCha[2], MCha2[2], MNeu[4], MNeu2[4];
  RealType MSS[3][4-1][2], MSS2[3][4-1][2], DSf[4][2];
  RealType MSf[3][4][2], MSf2[3][4][2], MSusy, MGl, MGl2;
  RealType MHiggs[4], MHiggs2[4], MHiggstree2[4];
  RealType CB, SB, TB, CB2, SB2, TB2, C2B, S2B;
  RealType CA, SA, CA2, SA2, C2A, S2A;
  RealType CAB, SAB, CBA, SBA, CBA2, SBA2;
  RealType AlfasMT;
} mssmpara_;

#define UCha(i,j) mssmpara_.UCha[j-1][i-1]
#define UChaC(i,j) Conjugate(UCha(i,j))
#define VCha(i,j) mssmpara_.VCha[j-1][i-1]
#define VChaC(i,j) Conjugate(VCha(i,j))
#define ZNeu(i,j) mssmpara_.ZNeu[j-1][i-1]
#define ZNeuC(i,j) Conjugate(ZNeu(i,j))
#define USf(i,j,t,g) mssmpara_.USf[g-1][t-1][j-1][i-1]
#define USfC(i,j,t,g) Conjugate(USf(i,j,t,g))
#define UCSf(i,j,t,g) mssmpara_.UCSfC[g-1][t-2][j-1][i-1]
#define UCSfC(i,j,t,g) Conjugate(UCSf(i,j,t,g))
#define UUSf(i,j,t,g) mssmpara_.UUSf[g-1][t-2][j-1][i-1]
#define UUSfC(i,j,t,g) Conjugate(UUSf(i,j,t,g))
#define USf2(i,j,t,g) Re(UCSf(i,j,t,g))
#define UHiggs(i,j) mssmpara_.XHiggs[0][j-1][i-1]
#define UHiggsC(i,j) Conjugate(UHiggs(i,j))
#define ZHiggs(i,j) mssmpara_.XHiggs[1][j-1][i-1]
#define ZHiggsC(i,j) Conjugate(ZHiggs(i,j))
#define Af(t,g1,g2) mssmpara_.Af[g2-1][g1-1][t-2]
#define AfC(t,g1,g2) Conjugate(Af(t,g1,g2))
#define Xf(t,g) mssmpara_.Xf[g-1][t-2]
#define MUETB(t) mssmpara_.MUETB[t-2]
#define Atau mssmpara_.Atau
#define At mssmpara_.At
#define Ab mssmpara_.Ab
#define MUE mssmpara_.MUE
#define MUEC Conjugate(MUE)
#define Mino1 mssmpara_.Mino1
#define Mino2 mssmpara_.Mino2
#define Mino3 mssmpara_.Mino3
#define Mino3C Conjugate(Mino3)
#define SqrtEGl mssmpara_.SqrtEGl
#define SqrtEGlC Conjugate(SqrtEGl)
#define MCha(i) mssmpara_.MCha[i-1]
#define MCha2(i) mssmpara_.MCha2[i-1]
#define MNeu(i) mssmpara_.MNeu[i-1]
#define MNeu2(i) mssmpara_.MNeu2[i-1]
#define MSS(n,t,g) mssmpara_.MSS[g-1][t-2][n-1]
#define MSS2(n,t,g) mssmpara_.MSS2[g-1][t-2][n-1]
#define DSf(i,t) mssmpara_.DSf[t-1][i-1]
#define MSf(i,t,g) mssmpara_.MSf[g-1][t-1][i-1]
#define MSf2(i,t,g) mssmpara_.MSf2[g-1][t-1][i-1]
#define MSusy mssmpara_.MSusy
#define MGl mssmpara_.MGl
#define MGl2 mssmpara_.MGl2
#define MHiggs(i) mssmpara_.MHiggs[i-1]
#define Mh0 MHiggs(1)
#define MHH MHiggs(2)
#define MA0 MHiggs(3)
#define MHp MHiggs(4)
#define MHiggs2(i) mssmpara_.MHiggs2[i-1]
#define Mh02 MHiggs2(1)
#define MHH2 MHiggs2(2)
#define MA02 MHiggs2(3)
#define MHp2 MHiggs2(4)
#define MHiggstree2(i) mssmpara_.MHiggstree2[i-1]
#define Mh0tree2 MHiggstree2(1)
#define MHHtree2 MHiggstree2(2)
#define MA0tree2 MHiggstree2(3)
#define MHptree2 MHiggstree2(4)
#define CB mssmpara_.CB
#define SB mssmpara_.SB
#define TB mssmpara_.TB
#define CB2 mssmpara_.CB2
#define SB2 mssmpara_.SB2
#define TB2 mssmpara_.TB2
#define C2B mssmpara_.C2B
#define S2B mssmpara_.S2B
#define CA mssmpara_.CA
#define SA mssmpara_.SA
#define CA2 mssmpara_.CA2
#define SA2 mssmpara_.SA2
#define C2A mssmpara_.C2A
#define S2A mssmpara_.S2A
#define CAB mssmpara_.CAB
#define SAB mssmpara_.SAB
#define CBA mssmpara_.CBA
#define SBA mssmpara_.SBA
#define CBA2 mssmpara_.CBA2
#define SBA2 mssmpara_.SBA2
#define AlfasMT mssmpara_.AlfasMT


// flavour-violating parameters

struct fvpara_ {
  ComplexType deltaSf[4-1][6][6];
  ComplexType UASf[4-1][6][6];
  RealType MASf[4-1][6], MASf2[4-1][6];
} fvpara_;

#define deltaSf(i,j,t) fvpara_.deltaSf[t-2][j-1][i-1]
#define deltaSf_LL(i,j,t) deltaSf(i,j,t)
#define deltaSf_LR(i,j,t) deltaSf(i,j+3,t)
#define deltaSf_RL(i,j,t) deltaSf(j,i+3,t)
#define deltaSf_RR(i,j,t) deltaSf(i+3,j+3,t)
#define UASf(i,j,t) fvpara_.UASf[t-2][j-1][i-1]
#define UASfC(i,j,t) Conjugate(UASf(i,j,t))
#define MASf(i,t) fvpara_.MASf[t-2][i-1]
#define MASf2(i,t) fvpara_.MASf[t-2][i-1]

#define deltaULR12 deltaSf_LR(1,2,3)
#define deltaULR23 deltaSf_LR(2,3,3)
#define deltaULR13 deltaSf_LR(1,3,3)

#define deltaURL12 deltaSf_RL(1,2,3)
#define deltaURL23 deltaSf_RL(2,3,3)
#define deltaURL13 deltaSf_RL(1,3,3)

#define deltaURR12 deltaSf_RR(1,2,3)
#define deltaURR23 deltaSf_RR(2,3,3)
#define deltaURR13 deltaSf_RR(1,3,3)

#define deltaQLL12 deltaSf_LL(1,2,4)
#define deltaQLL23 deltaSf_LL(2,3,4)
#define deltaQLL13 deltaSf_LL(1,3,4)

#define deltaDLR12 deltaSf_LR(1,2,4)
#define deltaDLR23 deltaSf_LR(2,3,4)
#define deltaDLR13 deltaSf_LR(1,3,4)

#define deltaDRL12 deltaSf_RL(1,2,4)
#define deltaDRL23 deltaSf_RL(2,3,4)
#define deltaDRL13 deltaSf_RL(1,3,4)

#define deltaDRR12 deltaSf_RR(1,2,4)
#define deltaDRR23 deltaSf_RR(2,3,4)
#define deltaDRR13 deltaSf_RR(1,3,4)

