#if 0
	FV.ch
	C declarations for extended flavor mixing
	this file is part of FormCalc
	last modified 6 Oct 19 th
#endif


#ifndef FV
#define FV 1
#endif

#if NOUNDERSCORE
#define fvpara_ fvpara
#endif

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

#define deltaLLL12 deltaSf_LL(1,2,2)
#define deltaLLL23 deltaSf_LL(2,3,2)
#define deltaLLL13 deltaSf_LL(1,3,2)

#define deltaELR12 deltaSf_LR(1,2,2)
#define deltaELR23 deltaSf_LR(2,3,2)
#define deltaELR13 deltaSf_LR(1,3,2)

#define deltaERL12 deltaSf_RL(1,2,2)
#define deltaERL23 deltaSf_RL(2,3,2)
#define deltaERL13 deltaSf_RL(1,3,2)

#define deltaERR12 deltaSf_RR(1,2,2)
#define deltaERR23 deltaSf_RR(2,3,2)
#define deltaERR13 deltaSf_RR(1,3,2)

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

