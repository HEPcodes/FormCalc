* inline.h
* inline versions of the util functions
* this file is part of FormCalc
* last modified 12 Apr 13 th


#ifndef INLINE_H
#define INLINE_H

#define signbit(i) ibits(i,31,1)
#define IndexDelta(i,j) signbit(ieor(i,j)-1)

#else

	integer IndexSign, IndexEps
	RealType SqDiff, ThreeMom

	integer a_, b_, c_
	RealType sqrtS_, ma_, mb_

	IndexSign(a_) = signbit(ior(a_, -a_)) - 2*signbit(a_)
	IndexEps(a_, b_, c_) =
     &    IndexSign(a_ - b_)*IndexSign(c_ - b_)*IndexSign(a_ - c_)

	SqDiff(ma_, mb_) = (ma_ - mb_)*(ma_ + mb_)
	ThreeMom(sqrtS_, ma_, mb_) = sqrt(SqDiff(
     &    .5D0*(sqrtS_ - SqDiff(ma_, mb_)/sqrtS_), mb_ ))

#endif

