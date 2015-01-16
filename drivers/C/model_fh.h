#if 0
	model_fh.h
	declarations for model_fh.F
	this file is part of FormCalc
	last modified 2 Oct 13 th
#endif


#if NOUNDERSCORE
#define smpara_ smpara
#define mssmpara_ mssmpara
#endif

struct smpara_ {
  ComplexType CKM[3][3];
  RealType Mf[3][4], Mf2[3][4];
  RealType MZ, MZ2, MW, MW2, MH, MH2, MBatMB;
  RealType CW, CW2, SW, SW2;
  RealType ELMZ, AlfaMZ, GF, AlfaGF, AlfasMZ;
  RealType EL, Alfa, Alfa2, GS, Alfas, Alfas2;
  RealType CKMlambda, CKMA, CKMrhobar, CKMetabar;
} smpara_;

#define CKM(i,j) smpara_.CKM[j-1][i-1]
#define CKMC(i,j) Conjugate(CKM(i,j))
#define Mf(t,g) smpara_.Mf[g-1][t-1]
#define Mf2(t,g) smpara_.Mf2[g-1][t-1]
#define MZ smpara_.MZ
#define MZ2 smpara_.MZ2
#define MW smpara_.MW
#define MW2 smpara_.MW2
#define MH smpara_.MH
#define MH2 smpara_.MH2
#define MBatMB smpara_.MBatMB
#define CW smpara_.CW
#define CW2 smpara_.CW2
#define SW smpara_.SW
#define SW2 smpara_.SW2
#define ELMZ smpara_.ELMZ
#define AlfaMZ smpara_.AlfaMZ
#define GF smpara_.GF
#define AlfaGF smpara_.AlfaGF
#define AlfasMZ smpara_.AlfasMZ
#define EL smpara_.EL
#define Alfa smpara_.Alfa
#define Alfa2 smpara_.Alfa2
#define GS smpara_.GS
#define Alfas smpara_.Alfas
#define Alfas2 smpara_.Alfas2
#define CKMlambda smpara_.CKMlambda
#define CKMA smpara_.CKMA
#define CKMrhobar smpara_.CKMrhobar
#define CKMetabar smpara_.CKMetabar

#define Alfa0 (1/137.0359895)

#define ME Mf(2,1)
#define MM Mf(2,2)
#define ML Mf(2,3)
#define MU Mf(3,1)
#define MC Mf(3,2)
#define MT Mf(3,3)
#define MD Mf(4,1)
#define MS Mf(4,2)
#define MB Mf(4,3)

#define ME2 Mf2(2,1)
#define MM2 Mf2(2,2)
#define ML2 Mf2(2,3)
#define MU2 Mf2(3,1)
#define MC2 Mf2(3,2)
#define MT2 Mf2(3,3)
#define MD2 Mf2(4,1)
#define MS2 Mf2(4,2)
#define MB2 Mf2(4,3)


struct mssmpara_ {
  ComplexType UCha[2][2], VCha[2][2], ZNeu[4][4];
  ComplexType XHiggs[2][3][3];
  ComplexType deltaSf[4][6][6], USf[3][5][2][2], UASf[5][6][6];
  ComplexType MSS2[5][3][3], Afd[3][4-2], Kf[4-2][3][3];
  ComplexType MUE, Mino1, Mino2, Mino3, SqrtEGl;
  RealType MCha[2], MCha2[2], MNeu[4], MNeu2[4];
  RealType MSS[3][5], MSf[3][5][2], MSf2[3][5][2];
  RealType MASf[5][6], MASf2[5][6];
  RealType MHiggs[4], MHiggs2[4], MHtree[4], MHtree2[4];
  RealType MGl, MGl2;
  RealType CB, SB, TB, CB2, SB2, TB2, C2B, S2B;
  RealType CA, SA, CA2, SA2, C2A, S2A;
  RealType CAB, SAB, CBA, SBA, CBA2, SBA2, SAeff;
  integer nmfv;
} mssmpara_;

#define UCha(i,j) mssmpara_.UCha[j-1][i-1]
#define UChaC(i,j) Conjugate(UCha(i,j))
#define VCha(i,j) mssmpara_.VCha[j-1][i-1]
#define VChaC(i,j) Conjugate(VCha(i,j))
#define ZNeu(i,j) mssmpara_.ZNeu[j-1][i-1]
#define ZNeuC(i,j) Conjugate(ZNeu(i,j))
#define deltaSf(i,j,t) mssmpara_.deltaSf[t-1][j-1][i-1]
#define deltaSf_LL(i,j,t) deltaSf(i,j,t)
#define deltaSf_LR(i,j,t) deltaSf(i,j+3,t)
#define deltaSf_RL(i,j,t) deltaSf(j,i+3,t)
#define deltaSf_RR(i,j,t) deltaSf(i+3,j+3,t)
#define UASf(i,j,t) mssmpara_.UASf[t-1][j-1][i-1]
#define UASfC(i,j,t) Conjugate(UASf(i,j,t))
#define USf(i,j,t,g) mssmpara_.USf[g-1][t-1][j-1][i-1]
#define USfC(i,j,t,g) Conjugate(USf(i,j,t,g))
#define MSS2(n,g1,g2) mssmpara_.MSS2[g2-1][g1-1][n-1]
#define Afd(t,g) mssmpara_.Afd[g-1][t-2]
#define Af(t,g1,g2) Mf(t,g1)*Kf(g1,g2,t)
#define AfC(t,g1,g2) Conjugate(Af(t,g1,g2))
#define Kf(g1,g2,t) mssmpara_.Kf[t-2][g2-1][g1-1]
#define KfC(g1,g2,t) Conjugate(Kf(g1,g2,t))
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
#define MSS(n,g) mssmpara_.MSS[g-1][n-1]
#define MSf(s,t,g) mssmpara_.MSf[g-1][t-1][s-1]
#define MSf2(s,t,g) mssmpara_.MSf2[g-1][t-1][s-1]
#define MASf(as,t) mssmpara_.MASf[t-1][as-1]
#define MASf2(as,t) mssmpara_.MASf[t-1][as-1]
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
#define MHtree(i) mssmpara_.MHtree[i-1]
#define Mh0tree MHtree(1)
#define MHHtree MHtree(2)
#define MA0tree MHtree(3)
#define MHptree MHtree(4)
#define MHtree2(i) mssmpara_.MHtree2[i-1]
#define Mh0tree2 MHtree2(1)
#define MHHtree2 MHtree2(2)
#define MA0tree2 MHtree2(3)
#define MHptree2 MHtree2(4)
#define MGl mssmpara_.MGl
#define MGl2 mssmpara_.MGl2
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
#define SAeff mssmpara_.SAeff
#define nmfv mssmpara_.nmfv

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

