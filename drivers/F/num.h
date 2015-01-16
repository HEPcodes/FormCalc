* num.h
* headers for the computation of the numerators
* this file is part of FormCalc
* last modified 10 Sep 14 th

#ifndef NUM_H
#define NUM_H

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
#define CNNum1 CNum1
#define CNNum2 CNum2
#define b0nj para(1)
#define b1nj para(2)
#define b2nj para(3)

#elif defined SAMURAI

#define NumFunction(f) function f(ncut, q1in, MuTildeSq)
#define Result(f) f
#define CSNum1 CNum1
#define CSNum2 CNum2

#else

#define NumFunction(f) subroutine f(q1in, res)
#define Result(f) res
#define CCNum1 CNum1
#define CCNum2 CNum2

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

CNNum1	integer ncut
CNNum1	ComplexType q1in(0:3)
CNNum2	Vec(1,1,q1) = q1in(0) + q1in(3)
CNNum2	Vec(2,2,q1) = q1in(0) - q1in(3)
CNNum2	Vec(2,1,q1) = q1in(1) + cI*q1in(2)
CNNum2	Vec(1,2,q1) = q1in(1) - cI*q1in(2)
CNNum1	ComplexType MuTildeSq

CMuExp1	integer ncut
CMuExp1	ComplexType vTin(0:3)
CMuExp2	Vec(1,1,vTnj) = vTin(0) + vTin(3)
CMuExp2	Vec(2,2,vTnj) = vTin(0) - vTin(3)
CMuExp2	Vec(2,1,vTnj) = vTin(1) + cI*vTin(2)
CMuExp2	Vec(1,2,vTnj) = vTin(1) - cI*vTin(2)

CT3Exp1	integer ncut, mindeg
CT3Exp1	ComplexType v0in(0:3)
CT3Exp2	Vec(1,1,v0nj) = v0in(0) + v0in(3)
CT3Exp2	Vec(2,2,v0nj) = v0in(0) - v0in(3)
CT3Exp2	Vec(2,1,v0nj) = v0in(1) + cI*v0in(2)
CT3Exp2	Vec(1,2,v0nj) = v0in(1) - cI*v0in(2)
CT3Exp1	ComplexType v3in(0:3)
CT3Exp2	Vec(1,1,v3nj) = v3in(0) + v3in(3)
CT3Exp2	Vec(2,2,v3nj) = v3in(0) - v3in(3)
CT3Exp2	Vec(2,1,v3nj) = v3in(1) + cI*v3in(2)
CT3Exp2	Vec(1,2,v3nj) = v3in(1) - cI*v3in(2)
CT3Exp1	ComplexType v4in(0:3)
CT3Exp2	Vec(1,1,v4nj) = v4in(0) + v4in(3)
CT3Exp2	Vec(2,2,v4nj) = v4in(0) - v4in(3)
CT3Exp2	Vec(2,1,v4nj) = v4in(1) + cI*v4in(2)
CT3Exp2	Vec(1,2,v4nj) = v4in(1) - cI*v4in(2)
CT3Exp1	ComplexType para(*)

CT2Exp1	integer ncut, mindeg
CT2Exp1	ComplexType v1in(0:3)
CT2Exp2	Vec(1,1,v1nj) = v1in(0) + v1in(3)
CT2Exp2	Vec(2,2,v1nj) = v1in(0) - v1in(3)
CT2Exp2	Vec(2,1,v1nj) = v1in(1) + cI*v1in(2)
CT2Exp2	Vec(1,2,v1nj) = v1in(1) - cI*v1in(2)
CT2Exp1	ComplexType v2in(0:3)
CT2Exp2	Vec(1,1,v2nj) = v2in(0) + v2in(3)
CT2Exp2	Vec(2,2,v2nj) = v2in(0) - v2in(3)
CT2Exp2	Vec(2,1,v2nj) = v2in(1) + cI*v2in(2)
CT2Exp2	Vec(1,2,v2nj) = v2in(1) - cI*v2in(2)
CT2Exp1	ComplexType v3in(0:3)
CT2Exp2	Vec(1,1,v3nj) = v3in(0) + v3in(3)
CT2Exp2	Vec(2,2,v3nj) = v3in(0) - v3in(3)
CT2Exp2	Vec(2,1,v3nj) = v3in(1) + cI*v3in(2)
CT2Exp2	Vec(1,2,v3nj) = v3in(1) - cI*v3in(2)
CT2Exp1	ComplexType v4in(0:3)
CT2Exp2	Vec(1,1,v4nj) = v4in(0) + v4in(3)
CT2Exp2	Vec(2,2,v4nj) = v4in(0) - v4in(3)
CT2Exp2	Vec(2,1,v4nj) = v4in(1) + cI*v4in(2)
CT2Exp2	Vec(1,2,v4nj) = v4in(1) - cI*v4in(2)
CT2Exp1	ComplexType para(*)

CSNum1	integer ncut
CSNum1	ComplexType q1in(4)
CSNum2	Vec(1,1,q1) = q1in(4) + q1in(3)
CSNum2	Vec(2,2,q1) = q1in(4) - q1in(3)
CSNum2	Vec(2,1,q1) = q1in(1) + cI*q1in(2)
CSNum2	Vec(1,2,q1) = q1in(1) - cI*q1in(2)
CSNum1	RealType MuTildeSq
CSNum2	muscale = MuTildeSq

CCNum1	ComplexType q1in(0:3)
CCNum2	Vec(1,1,q1) = q1in(0) + q1in(3)
CCNum2	Vec(2,2,q1) = q1in(0) - q1in(3)
CCNum2	Vec(2,1,q1) = q1in(1) + cI*q1in(2)
CCNum2	Vec(1,2,q1) = q1in(1) - cI*q1in(2)

#endif

