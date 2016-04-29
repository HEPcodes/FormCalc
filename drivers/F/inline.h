* inline.h
* inline util functions and computation of numerators
* this file is part of FormCalc
* last modified 26 Feb 16 th


#ifndef INLINE_H
#define INLINE_H

#ifdef NINJA

#define NumFunction(f) subroutine f(ncut, q1in, MuTildeSq, res)
#define Result(f) res
#define MuExpFunction(f) subroutine f(ncut, vTin, njcoeff)
#define T3ExpFunction(f) subroutine f(ncut, v0in, v3in, v4in, para, mindeg, njcoeff)
#define T3ExpCoeff(d,m) if( d .gt. mindeg ) return
* i.e. coefficient of tnj^(rank-d) MuTildeSq^m
#define T2ExpFunction(f) subroutine f(ncut, v1in, v2in, v3in, v4in, para, mindeg, njcoeff)
#define T2ExpCoeff(d,x,m) if( d .gt. mindeg ) return
* i.e. coefficient of tnj^(rank-d) xnj^x MuTildeSq^m
#define CNNum CNum
#define b0nj para(1)
#define b1nj para(2)
#define b2nj para(3)

#elif defined SAMURAI

#define NumFunction(f) function f(ncut, q1in, MuTildeSq)
#define Result(f) f
#define CSNum CNum

#elif defined CUTTOOLS

#define NumFunction(f) subroutine f(q1in, res)
#define Result(f) res
#define CCNum CNum

#endif

#if NUMDEBUG > 1
#define NumDebug(f) \
print '(A,/4G14.5,A/4G14.5/)', f, \
  Re(q1in), " q1", Im(q1in)
#define MuExpDebug(f) \
print '(A/4G14.5,A/4G14.5/)', f, \
  Re(vTin), " vT", Im(vTin)
#define T3ExpDebug(f) \
print '(A," mindeg ",I1,3(/4G14.5,A/4G14.5/))', f, mindeg, \
  Re(v0in), " v0", Im(v0in), \
  Re(v3in), " v3", Im(v3in), \
  Re(v4in), " v4", Im(v4in)
#define T2ExpDebug(f) \
print '(A," mindeg ",I1,4(/4G14.5,A/4G14.5/))', f, mindeg, \
  Re(v1in), " v1", Im(v1in), \
  Re(v2in), " v2", Im(v2in), \
  Re(v3in), " v3", Im(v3in), \
  Re(v4in), " v4", Im(v4in)
#elif NUMDEBUG == 1
#define NumDebug(f) print '(A)', f
#define MuExpDebug(f) print '(A)', f
#define T3ExpDebug(f) print '(A," mindeg ",I1)', f, mindeg
#define T2ExpDebug(f) print '(A," mindeg ",I1)', f, mindeg
#else
#define NumDebug(f)
#define MuExpDebug(f)
#define T3ExpDebug(f)
#define T2ExpDebug(f)
#endif

#else

CNNum	integer ncut
CNNum	ComplexType q1in(0:3)
CNNum	ComplexType MuTildeSq

CSNum	integer ncut
CSNum	ComplexType q1in(4)
CSNum	RealType MuTildeSq

CCNum	ComplexType q1in(0:3)

CMuExp	integer ncut
CMuExp	ComplexType vTin(0:3)

CT3Exp	integer ncut, mindeg
CT3Exp	ComplexType v0in(0:3)
CT3Exp	ComplexType v3in(0:3)
CT3Exp	ComplexType v4in(0:3)
CT3Exp	ComplexType para(*)

CT2Exp	integer ncut, mindeg
CT2Exp	ComplexType v1in(0:3)
CT2Exp	ComplexType v2in(0:3)
CT2Exp	ComplexType v3in(0:3)
CT2Exp	ComplexType v4in(0:3)
CT2Exp	ComplexType para(*)

	integer IndexDelta, IndexSign, IndexEps
	RealType Sq, SqDiff, ThreeMom, SInvariant, TInvariant
	ComplexType Pair0, Eps0, Eps0_

	integer a_, b_, c_, d_
	RealType sqrtS_, ma_, mb_
	ComplexType z_

	IndexDelta(a_, b_) = merge(1, 0, a_ .eq. b_)

	IndexSign(a_) = signbit(ior(a_, -a_)) - 2*signbit(a_)
	IndexEps(a_, b_, c_) =
     &    IndexSign(a_ - b_)*IndexSign(c_ - b_)*IndexSign(a_ - c_)

	Sq(z_) = Re(z_*Conjugate(z_))

	SqDiff(ma_, mb_) = (ma_ - mb_)*(ma_ + mb_)
	ThreeMom(sqrtS_, ma_, mb_) = sqrt(SqDiff(
     &    .5D0*(sqrtS_ - SqDiff(ma_, mb_)/sqrtS_), mb_ ))

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

CNNum	Vec(1,1,q1) = q1in(0) + q1in(3)
CNNum	Vec(2,2,q1) = q1in(0) - q1in(3)
CNNum	Vec(2,1,q1) = q1in(1) + cI*q1in(2)
CNNum	Vec(1,2,q1) = q1in(1) - cI*q1in(2)

CSNum	Vec(1,1,q1) = q1in(4) + q1in(3)
CSNum	Vec(2,2,q1) = q1in(4) - q1in(3)
CSNum	Vec(2,1,q1) = q1in(1) + cI*q1in(2)
CSNum	Vec(1,2,q1) = q1in(1) - cI*q1in(2)
CSNum	muscale = MuTildeSq

CCNum	Vec(1,1,q1) = q1in(0) + q1in(3)
CCNum	Vec(2,2,q1) = q1in(0) - q1in(3)
CCNum	Vec(2,1,q1) = q1in(1) + cI*q1in(2)
CCNum	Vec(1,2,q1) = q1in(1) - cI*q1in(2)

CMuExp	Vec(1,1,vTnj) = vTin(0) + vTin(3)
CMuExp	Vec(2,2,vTnj) = vTin(0) - vTin(3)
CMuExp	Vec(2,1,vTnj) = vTin(1) + cI*vTin(2)
CMuExp	Vec(1,2,vTnj) = vTin(1) - cI*vTin(2)

CT3Exp	Vec(1,1,v0nj) = v0in(0) + v0in(3)
CT3Exp	Vec(2,2,v0nj) = v0in(0) - v0in(3)
CT3Exp	Vec(2,1,v0nj) = v0in(1) + cI*v0in(2)
CT3Exp	Vec(1,2,v0nj) = v0in(1) - cI*v0in(2)
CT3Exp	Vec(1,1,v3nj) = v3in(0) + v3in(3)
CT3Exp	Vec(2,2,v3nj) = v3in(0) - v3in(3)
CT3Exp	Vec(2,1,v3nj) = v3in(1) + cI*v3in(2)
CT3Exp	Vec(1,2,v3nj) = v3in(1) - cI*v3in(2)
CT3Exp	Vec(1,1,v4nj) = v4in(0) + v4in(3)
CT3Exp	Vec(2,2,v4nj) = v4in(0) - v4in(3)
CT3Exp	Vec(2,1,v4nj) = v4in(1) + cI*v4in(2)
CT3Exp	Vec(1,2,v4nj) = v4in(1) - cI*v4in(2)

CT2Exp	Vec(1,1,v1nj) = v1in(0) + v1in(3)
CT2Exp	Vec(2,2,v1nj) = v1in(0) - v1in(3)
CT2Exp	Vec(2,1,v1nj) = v1in(1) + cI*v1in(2)
CT2Exp	Vec(1,2,v1nj) = v1in(1) - cI*v1in(2)
CT2Exp	Vec(1,1,v2nj) = v2in(0) + v2in(3)
CT2Exp	Vec(2,2,v2nj) = v2in(0) - v2in(3)
CT2Exp	Vec(2,1,v2nj) = v2in(1) + cI*v2in(2)
CT2Exp	Vec(1,2,v2nj) = v2in(1) - cI*v2in(2)
CT2Exp	Vec(1,1,v3nj) = v3in(0) + v3in(3)
CT2Exp	Vec(2,2,v3nj) = v3in(0) - v3in(3)
CT2Exp	Vec(2,1,v3nj) = v3in(1) + cI*v3in(2)
CT2Exp	Vec(1,2,v3nj) = v3in(1) - cI*v3in(2)
CT2Exp	Vec(1,1,v4nj) = v4in(0) + v4in(3)
CT2Exp	Vec(2,2,v4nj) = v4in(0) - v4in(3)
CT2Exp	Vec(2,1,v4nj) = v4in(1) + cI*v4in(2)
CT2Exp	Vec(1,2,v4nj) = v4in(1) - cI*v4in(2)

#endif

