* inline.h
* inline versions of the util functions
* this file is part of FormCalc
* last modified 26 Oct 15 th


#ifndef INLINE_H
#define INLINE_H

#define signbit(i) ibits(i,31,1)
#define IndexDelta(i,j) signbit(ieor(i,j)-1)

#else

	integer IndexSign, IndexEps
	RealType SqDiff, ThreeMom, SInvariant, TInvariant
	ComplexType Pair0, Eps0, Eps0_

	integer a_, b_, c_, d_
	RealType sqrtS_, ma_, mb_

	SInvariant(a_, b_) =
     &    (Re(vec0(1,1,k0(a_))) + Re(vec0(1,1,k0(b_))))*
     &    (Re(vec0(2,2,k0(a_))) + Re(vec0(2,2,k0(b_)))) -
     &    Sq(vec0(1,2,k0(a_)) + vec0(1,2,k0(b_)))

	TInvariant(a_, b_) =
     &    (Re(vec0(1,1,k0(a_))) - Re(vec0(1,1,k0(b_))))*
     &    (Re(vec0(2,2,k0(a_))) - Re(vec0(2,2,k0(b_)))) -
     &    Sq(vec0(1,2,k0(a_)) - vec0(1,2,k0(b_)))

	Pair0(a_, b_) = .5D0*(
     &    vec0(1,1,a_)*vec0(2,2,b_) + vec0(2,2,a_)*vec0(1,1,b_) -
     &    vec0(1,2,a_)*vec0(2,1,b_) - vec0(2,1,a_)*vec0(1,2,b_) )

	Eps0_(a_, b_, c_, d_) =
     &    (vec0(1,1,a_)*vec0(2,2,b_) - vec0(2,2,a_)*vec0(1,1,b_))*
     &    (vec0(2,1,c_)*vec0(1,2,d_) - vec0(1,2,c_)*vec0(2,1,d_))
	Eps0(a_, b_, c_, d_) = .25D0*(
     &    Eps0_(a_, b_, c_, d_) + Eps0_(c_, d_, a_, b_) -
     &    Eps0_(a_, c_, b_, d_) - Eps0_(b_, d_, a_, c_) +
     &    Eps0_(a_, d_, b_, c_) + Eps0_(b_, c_, a_, d_) )

	IndexSign(a_) = signbit(ior(a_, -a_)) - 2*signbit(a_)
	IndexEps(a_, b_, c_) =
     &    IndexSign(a_ - b_)*IndexSign(c_ - b_)*IndexSign(a_ - c_)

	SqDiff(ma_, mb_) = (ma_ - mb_)*(ma_ + mb_)
	ThreeMom(sqrtS_, ma_, mb_) = sqrt(SqDiff(
     &    .5D0*(sqrtS_ - SqDiff(ma_, mb_)/sqrtS_), mb_ ))

#endif

