* model_fv.h
* flavor-violating extension for model_mssm.h and model_fh.h
* this file is part of FormCalc
* last modified 31 Jul 14 th


	ComplexType deltaSf(6,6,4)
	ComplexType UASf(6,6,5)
	RealType MASf(6,5), MASf2(6,5)

	common /fvpara/ deltaSf, UASf, MASf, MASf2

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

	ComplexType deltaSf_flat(6*6*4)
	equivalence (deltaSf, deltaSf_flat)

	RealType ReImdeltaSf(2,6,6,4)
	equivalence (deltaSf, ReImdeltaSf)

	RealType MASf_flat(6*5), MASf2_flat(6*5)
	equivalence (MASf, MASf_flat)
	equivalence (MASf2, MASf2_flat)

* The following named deltas are needed because Fortran
* does not accept an array element as loop variable:

	ComplexType deltaLLL12
	RealType RedeltaLLL12, ImdeltaLLL12
	equivalence (deltaSf_LL(1,2,2), deltaLLL12, RedeltaLLL12)
	equivalence (ImdeltaSf_LL(1,2,2), ImdeltaLLL12)
	ComplexType deltaLLL23
	RealType RedeltaLLL23, ImdeltaLLL23
	equivalence (deltaSf_LL(2,3,2), deltaLLL23, RedeltaLLL23)
	equivalence (ImdeltaSf_LL(2,3,2), ImdeltaLLL23)
	ComplexType deltaLLL13
	RealType RedeltaLLL13, ImdeltaLLL13
	equivalence (deltaSf_LL(1,3,2), deltaLLL13, RedeltaLLL13)
	equivalence (ImdeltaSf_LL(1,3,2), ImdeltaLLL13)

	ComplexType deltaELR12
	RealType RedeltaELR12, ImdeltaELR12
	equivalence (deltaSf_LR(1,2,2), deltaELR12, RedeltaELR12)
	equivalence (ImdeltaSf_LR(1,2,2), ImdeltaELR12)
	ComplexType deltaELR23
	RealType RedeltaELR23, ImdeltaELR23
	equivalence (deltaSf_LR(2,3,2), deltaELR23, RedeltaELR23)
	equivalence (ImdeltaSf_LR(2,3,2), ImdeltaELR23)
	ComplexType deltaELR13
	RealType RedeltaELR13, ImdeltaELR13
	equivalence (deltaSf_LR(1,3,2), deltaELR13, RedeltaELR13)
	equivalence (ImdeltaSf_LR(1,3,2), ImdeltaELR13)

	ComplexType deltaERL12
	RealType RedeltaERL12, ImdeltaERL12
	equivalence (deltaSf_RL(1,2,2), deltaERL12, RedeltaERL12)
	equivalence (ImdeltaSf_RL(1,2,2), ImdeltaERL12)
	ComplexType deltaERL23
	RealType RedeltaERL23, ImdeltaERL23
	equivalence (deltaSf_RL(2,3,2), deltaERL23, RedeltaERL23)
	equivalence (ImdeltaSf_RL(2,3,2), ImdeltaERL23)
	ComplexType deltaERL13
	RealType RedeltaERL13, ImdeltaERL13
	equivalence (deltaSf_RL(1,3,2), deltaERL13, RedeltaERL13)
	equivalence (ImdeltaSf_RL(1,3,2), ImdeltaERL13)

	ComplexType deltaERR12
	RealType RedeltaERR12, ImdeltaERR12
	equivalence (deltaSf_RR(1,2,2), deltaERR12, RedeltaERR12)
	equivalence (ImdeltaSf_RR(1,2,2), ImdeltaERR12)
	ComplexType deltaERR23
	RealType RedeltaERR23, ImdeltaERR23
	equivalence (deltaSf_RR(2,3,2), deltaERR23, RedeltaERR23)
	equivalence (ImdeltaSf_RR(2,3,2), ImdeltaERR23)
	ComplexType deltaERR13
	RealType RedeltaERR13, ImdeltaERR13
	equivalence (deltaSf_RR(1,3,2), deltaERR13, RedeltaERR13)
	equivalence (ImdeltaSf_RR(1,3,2), ImdeltaERR13)


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
