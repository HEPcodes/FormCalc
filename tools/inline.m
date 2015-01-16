SInvariant[a_, b_] =
  (Re[vec[1,1,a]] + Re[vec[1,1,b]])*
  (Re[vec[2,2,a]] + Re[vec[2,2,b]]) -
  Sq[vec[1,2,a] + vec[1,2,b]]

TInvariant[a_, b_] =
  (Re[vec[1,1,a]] - Re[vec[1,1,b]])*
  (Re[vec[2,2,a]] - Re[vec[2,2,b]]) -
  Sq[vec[1,2,a] - vec[1,2,b]]

Pair[a_, b_] = .5D0*(
  vec[1,1,a]*vec[2,2,b] + vec[2,2,a]*vec[1,1,b] -
  vec[1,2,a]*vec[2,1,b] - vec[2,1,a]*vec[1,2,b] )

eps[a_, b_, c_, d_] =
  (vec[1,1,a]*vec[2,2,b] - vec[2,2,a]*vec[1,1,b])*
  (vec[2,1,c]*vec[1,2,d] - vec[1,2,c]*vec[2,1,d])

Eps[a_, b_, c_, d_] = .25D0*(
  eps[a, b, c, d] + eps[c, d, a, b] -
  eps[a, c, b, d] - eps[b, d, a, c] +
  eps[a, d, b, c] + eps[b, c, a, d] )


bar[{{A11_, A12_}, {A21_, A22_}}] := {{A22, -A12}, {-A21, A11}}

SxS[l_, r_] = l.r

SxV[l_, A_] = l.A

SxB[l_, A_] = l.bar[A]

VxS[A_, r_] = A.r

BxS[A_, r_] = bar[A].r

SxVxB[l_, a_, b_] = SxB[SxV[l, a], b]

SxBxV[l_, a_, b_] = SxV[SxB[l, a], b]

BxVxS[b_, a_, r_] = BxS[b_, VxS[a, r]]

VxBxS[b_, a_, r_] = VxS[b_, BxS[a, r]]

VxBxS2[b_, a_, r1_,r2_] =
  VxS2[b_, BxS1[a_, r1_,r2_],BxS2[a_, r1_,r2_]]

SxVxBxV1[l1_,l2_, a_, b_, c_] =
  SxBxV1[SxV1[l1_,l2_, a_],SxV2[l1_,l2_, a_], b_, c_]
SxVxBxV2[l1_,l2_, a_, b_, c_] =
  SxBxV2[SxV1[l1_,l2_, a_],SxV2[l1_,l2_, a_], b_, c_]

SxBxVxB1[l1_,l2_, a_, b_, c_] =
  SxVxB1[SxB1[l1_,l2_, a_],SxB2[l1_,l2_, a_], b_, c_]
SxBxVxB2[l1_,l2_, a_, b_, c_] =
  SxVxB2[SxB1[l1_,l2_, a_],SxB2[l1_,l2_, a_], b_, c_]

VxBxVxS1[c_, b_, a_, r1_,r2_] =
  VxBxS1[c_, b_, VxS1[a_, r1_,r2_],VxS2[a_, r1_,r2_]]
VxBxVxS2[c_, b_, a_, r1_,r2_] =
  VxBxS2[c_, b_, VxS1[a_, r1_,r2_],VxS2[a_, r1_,r2_]]

BxVxBxS1[c_, b_, a_, r1_,r2_] =
  BxVxS1[c_, b_, BxS1[a_, r1_,r2_],BxS2[a_, r1_,r2_]]
BxVxBxS2[c_, b_, a_, r1_,r2_] =
  BxVxS2[c_, b_, BxS1[a_, r1_,r2_],BxS2[a_, r1_,r2_]]

#define SpiLV[iL,eL] [1-2*eL]*vec[1+eL,1+eL,iL], vec[2-eL,1+eL,iL]
#define SpiLB[iL,eL] [1-2*eL]*vec[1+eL,2-eL,iL], vec[2-eL,2-eL,iL]

#define SpiRV[eR,iR] vec[1+eR,1+eR,iR], [1-2*eR]*vec[2-eR,1+eR,iR]
#define SpiRB[eR,iR] vec[1+eR,2-eR,iR], [1-2*eR]*vec[2-eR,2-eR,iR]

ChainV0[Spinor[iL_,sL_,dL_], eL_, eR_, Spinor[iR_,sR_,dR_]] :=
  <<<

	ChainV0[iL_,eL_, eR_,iR_] = SxS[
          SpiLB[iL_,eL_],
          SpiRV[eR_,iR_] ]
	ChainB0[iL_,eL_, eR_,iR_] = SxS[
          SpiLV[iL_,eL_],
          SpiRB[eR_,iR_] ]

	ChainV1[iL_,eL_, a_, eR_,iR_] = SxS[
          SxV1[SpiLB[iL_,eL_], a_],
          SxV2[SpiLB[iL_,eL_], a_],
          SpiRB[eR_,iR_] ]
	ChainB1[iL_,eL_, a_, eR_,iR_] = SxS[
          SxB1[SpiLV[iL_,eL_], a_],
          SxB2[SpiLV[iL_,eL_], a_],
          SpiRV[eR_,iR_] ]

	ChainV2[iL_,eL_, a_, b_, eR_,iR_] = SxS[
          SxV1[SpiLB[iL_,eL_], a_],
          SxV2[SpiLB[iL_,eL_], a_],
          BxS1[b_, SpiRV[eR_,iR_]],
          BxS2[b_, SpiRV[eR_,iR_]] ]
	ChainB2[iL_,eL_, a_, b_, eR_,iR_] = SxS[
          SxB1[SpiLV[iL_,eL_], a_],
          SxB2[SpiLV[iL_,eL_], a_],
          VxS1[b_, SpiRB[eR_,iR_]],
          VxS2[b_, SpiRB[eR_,iR_]] ]

	ChainV3[iL_,eL_, a_, b_, c_, eR_,iR_] = SxS[
          SxVxB1[SpiLB[iL_,eL_], a_, b_],
          SxVxB2[SpiLB[iL_,eL_], a_, b_],
          VxS1[c_, SpiRB[eR_,iR_]],
          VxS2[c_, SpiRB[eR_,iR_]] ]
	ChainB3[iL_,eL_, a_, b_, c_, eR_,iR_] = SxS[
          SxBxV1[SpiLV[iL_,eL_], a_, b_],
          SxBxV2[SpiLV[iL_,eL_], a_, b_],
          BxS1[c_, SpiRV[eR_,iR_]],
          BxS2[c_, SpiRV[eR_,iR_]] ]

	ChainV4[iL_,eL_, a_, b_, c_, d_, eR_,iR_] = SxS[
          SxVxB1[SpiLB[iL_,eL_], a_, b_],
          SxVxB2[SpiLB[iL_,eL_], a_, b_],
          VxBxS1[c_, d_, SpiRV[eR_,iR_]],
          VxBxS2[c_, d_, SpiRV[eR_,iR_]] ]
	ChainB4[iL_,eL_, a_, b_, c_, d_, eR_,iR_] = SxS[
          SxBxV1[SpiLV[iL_,eL_], a_, b_],
          SxBxV2[SpiLV[iL_,eL_], a_, b_],
          BxVxS1[c_, d_, SpiRB[eR_,iR_]],
          BxVxS2[c_, d_, SpiRB[eR_,iR_]] ]

	ChainV5[iL_,eL_, a_, b_, c_, d_, e_, eR_,iR_] = SxS[
          SxVxBxV1[SpiLB[iL_,eL_], a_, b_, c_],
          SxVxBxV2[SpiLB[iL_,eL_], a_, b_, c_],
          BxVxS1[d_, e_, SpiRB[eR_,iR_]],
          BxVxS2[d_, e_, SpiRB[eR_,iR_]] ]
	ChainB5[iL_,eL_, a_, b_, c_, d_, e_, eR_,iR_] = SxS[
          SxBxVxB1[SpiLV[iL_,eL_], a_, b_, c_],
          SxBxVxB2[SpiLV[iL_,eL_], a_, b_, c_],
          VxBxS1[d_, e_, SpiRV[eR_,iR_]],
          VxBxS2[d_, e_, SpiRV[eR_,iR_]] ]

	ChainV6[iL_,eL_, a_, b_, c_, d_, e_, f_, eR_,iR_] = SxS[
          SxVxBxV1[SpiLB[iL_,eL_], a_, b_, c_],
          SxVxBxV2[SpiLB[iL_,eL_], a_, b_, c_],
          BxVxBxS1[d_, e_, f_, SpiRV[eR_,iR_]],
          BxVxBxS2[d_, e_, f_, SpiRV[eR_,iR_]] ]
	ChainB6[iL_,eL_, a_, b_, c_, d_, e_, f_, eR_,iR_] = SxS[
          SxBxVxB1[SpiLV[iL_,eL_], a_, b_, c_],
          SxBxVxB2[SpiLV[iL_,eL_], a_, b_, c_],
          VxBxVxS1[d_, e_, f_, SpiRB[eR_,iR_]],
          VxBxVxS2[d_, e_, f_, SpiRB[eR_,iR_]] ]

	IndexEps_[a_] = 2*signbit[a_] - signbit[ior[a_, -a_]]
	IndexEps[a_, b_, c_] =
          IndexEps_[b_ - a_]*IndexEps_[b_ - c_]*IndexEps_[c_ - a_]

	ThreeMom_[ma_, mb_] = [ma_ - mb_]*[ma_ + mb_]
	ThreeMom[sqrtS_, ma_, mb_] = sqrt[ThreeMom_[
          .5D0*[sqrtS_ - ThreeMom_[ma_, mb_]/sqrtS_], mb_ ]]

