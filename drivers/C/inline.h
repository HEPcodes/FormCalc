#if 0
	inline.h
	inline versions of the util functions
	this file is part of FormCalc
	last modified 10 Mar 13 th
#endif


#ifndef INLINE_H
#define INLINE_H

#define ComplexFun static inline ComplexType
#define SpiType(s1,s2) cComplexType s1, cComplexType s2
#define SpiSpec(i,e) cint i, cint e
#define SpiLV(i,e) (1-2*e)*vec(1+e,1+e,i), vec(2-e,1+e,i)
#define SpiLB(i,e) (1-2*e)*vec(1+e,2-e,i), vec(2-e,2-e,i)
#define SpiRV(e,i) vec(1+e,1+e,i), (1-2*e)*vec(2-e,1+e,i)
#define SpiRB(e,i) vec(1+e,2-e,i), (1-2*e)*vec(2-e,2-e,i)

enum { MBbits = 8, MBmask = (1LL << MBbits) - 1 };
#define MomEncoding(f,i) (((f) & (JK-1)) << MBbits*(i-1))

#define Sqrt sqrt

#define Power(x, n) \
  (sizeof(x) == sizeof(ComplexType) ? CPower(x, n) : RPower(x, n))

static inline RealType RPower(RealType x, int n) {
  RealType result = 1;
  while( n ) {
    if( n & 1 ) result *= x;
    n >>= 1;
    x *= x;
  }
  return result;
}

static inline ComplexType CPower(ComplexType x, int n) {
  ComplexType result = 1;
  while( n ) {
    if( n & 1 ) result *= x;
    n >>= 1;
    x *= x;
  }
  return result;
}

static inline RealType Sq(cComplexType c) {
  return Re(c*Conjugate(c));
}

ComplexFun SInvariant(cint a, cint b) {
  return (Re(vec(1,1,a)) + Re(vec(1,1,b)))*
         (Re(vec(2,2,a)) + Re(vec(2,2,b))) -
         Sq(vec(1,2,a) + vec(1,2,b));
}

ComplexFun TInvariant(cint a, cint b) {
  return (Re(vec(1,1,a)) - Re(vec(1,1,b)))*
         (Re(vec(2,2,a)) - Re(vec(2,2,b))) -
         Sq(vec(1,2,a) - vec(1,2,b));
}

ComplexFun Pair(cint a, cint b) {
  return .5*(vec(1,1,a)*vec(2,2,b) + vec(2,2,a)*vec(1,1,b) -
             vec(1,2,a)*vec(2,1,b) - vec(2,1,a)*vec(1,2,b));
}

ComplexFun Eps_(cint a, cint b, cint c, cint d) {
  return (vec(1,1,a)*vec(2,2,b) - vec(2,2,a)*vec(1,1,b))*
         (vec(2,1,c)*vec(1,2,d) - vec(1,2,c)*vec(2,1,d));
}

ComplexFun Eps(cint a, cint b, cint c, cint d) {
  return .25*(
    Eps_(a, b, c, d) + Eps_(c, d, a, b) -
    Eps_(a, c, b, d) - Eps_(b, d, a, c) +
    Eps_(a, d, b, c) + Eps_(b, c, a, d) );
}

ComplexFun SxS(SpiType(l1,l2), SpiType(r1,r2)) {
  return l1*r1 + l2*r2;
}

ComplexFun SxV1(SpiType(l1,l2), cint a) {
  return l1*vec(1,1,a) + l2*vec(2,1,a);
}

ComplexFun SxV2(SpiType(l1,l2), cint a) {
  return l2*vec(2,2,a) + l1*vec(1,2,a);
}

ComplexFun SxB1(SpiType(l1,l2), cint a) {
  return l1*vec(2,2,a) - l2*vec(2,1,a);
}

ComplexFun SxB2(SpiType(l1,l2), cint a) {
  return l2*vec(1,1,a) - l1*vec(1,2,a);
}

ComplexFun VxS1(cint a, SpiType(r1,r2)) {
  return vec(1,1,a)*r1 + vec(1,2,a)*r2;
}

ComplexFun VxS2(cint a, SpiType(r1,r2)) {
  return vec(2,1,a)*r1 + vec(2,2,a)*r2;
}

ComplexFun BxS1(cint a, SpiType(r1,r2)) {
  return vec(2,2,a)*r1 - vec(1,2,a)*r2;
}

ComplexFun BxS2(cint a, SpiType(r1,r2)) {
  return vec(1,1,a)*r2 - vec(2,1,a)*r1;
}

ComplexFun SxVxB1(SpiType(l1,l2), cint a, cint b) {
  return SxB1(SxV1(l1,l2, a),SxV2(l1,l2, a), b);
}

ComplexFun SxVxB2(SpiType(l1,l2), cint a, cint b) {
  return SxB2(SxV1(l1,l2, a),SxV2(l1,l2, a), b);
}

ComplexFun SxBxV1(SpiType(l1,l2), cint a, cint b) {
  return SxV1(SxB1(l1,l2, a),SxB2(l1,l2, a), b);
}

ComplexFun SxBxV2(SpiType(l1,l2), cint a, cint b) {
  return SxV2(SxB1(l1,l2, a),SxB2(l1,l2, a), b);
}

ComplexFun BxVxS1(cint b, cint a, SpiType(r1,r2)) {
  return BxS1(b, VxS1(a, r1,r2),VxS2(a, r1,r2));
}

ComplexFun BxVxS2(cint b, cint a, SpiType(r1,r2)) {
  return BxS2(b, VxS1(a, r1,r2),VxS2(a, r1,r2));
}

ComplexFun VxBxS1(cint b, cint a, SpiType(r1,r2)) {
  return VxS1(b, BxS1(a, r1,r2),BxS2(a, r1,r2));
}

ComplexFun VxBxS2(cint b, cint a, SpiType(r1,r2)) {
  return VxS2(b, BxS1(a, r1,r2),BxS2(a, r1,r2));
}

ComplexFun SxVxBxV1(SpiType(l1,l2), cint a, cint b, cint c) {
  return SxBxV1(SxV1(l1,l2, a),SxV2(l1,l2, a), b, c);
}

ComplexFun SxVxBxV2(SpiType(l1,l2), cint a, cint b, cint c) {
  return SxBxV2(SxV1(l1,l2, a),SxV2(l1,l2, a), b, c);
}

ComplexFun SxBxVxB1(SpiType(l1,l2), cint a, cint b, cint c) {
  return SxVxB1(SxB1(l1,l2, a),SxB2(l1,l2, a), b, c);
}

ComplexFun SxBxVxB2(SpiType(l1,l2), cint a, cint b, cint c) {
  return SxVxB2(SxB1(l1,l2, a),SxB2(l1,l2, a), b, c);
}

ComplexFun VxBxVxS1(cint c, cint b, cint a, SpiType(r1,r2)) {
  return VxBxS1(c, b, VxS1(a, r1,r2),VxS2(a, r1,r2));
}

ComplexFun VxBxVxS2(cint c, cint b, cint a, SpiType(r1,r2)) {
  return VxBxS2(c, b, VxS1(a, r1,r2),VxS2(a, r1,r2));
}

ComplexFun BxVxBxS1(cint c, cint b, cint a, SpiType(r1,r2)) {
  return BxVxS1(c, b, BxS1(a, r1,r2),BxS2(a, r1,r2));
}

ComplexFun BxVxBxS2(cint c, cint b, cint a, SpiType(r1,r2)) {
  return BxVxS2(c, b, BxS1(a, r1,r2),BxS2(a, r1,r2));
}

ComplexFun ChainV0(SpiSpec(iL,eL), SpiSpec(eR,iR)) {
  return SxS(SpiLB(iL,eL), SpiRV(eR,iR));
}

ComplexFun ChainB0(SpiSpec(iL,eL), SpiSpec(eR,iR)) {
  return SxS(SpiLV(iL,eL),
             SpiRB(eR,iR));
}

ComplexFun ChainV1(SpiSpec(iL,eL), cint a, SpiSpec(eR,iR)) {
  return SxS(SxV1(SpiLB(iL,eL), a),
             SxV2(SpiLB(iL,eL), a),
             SpiRB(eR,iR));
}

ComplexFun ChainB1(SpiSpec(iL,eL), cint a, SpiSpec(eR,iR)) {
  return SxS(SxB1(SpiLV(iL,eL), a),
             SxB2(SpiLV(iL,eL), a),
             SpiRV(eR,iR));
}

ComplexFun ChainV2(SpiSpec(iL,eL), cint a, cint b, SpiSpec(eR,iR)) {
  return SxS(SxV1(SpiLB(iL,eL), a),
             SxV2(SpiLB(iL,eL), a),
             BxS1(b, SpiRV(eR,iR)),
             BxS2(b, SpiRV(eR,iR)));
}

ComplexFun ChainB2(SpiSpec(iL,eL), cint a, cint b, SpiSpec(eR,iR)) {
  return SxS(SxB1(SpiLV(iL,eL), a),
             SxB2(SpiLV(iL,eL), a),
             VxS1(b, SpiRB(eR,iR)),
             VxS2(b, SpiRB(eR,iR)));
}

ComplexFun ChainV3(SpiSpec(iL,eL), cint a, cint b,
    cint c, SpiSpec(eR,iR)) {
  return SxS(SxVxB1(SpiLB(iL,eL), a, b),
             SxVxB2(SpiLB(iL,eL), a, b),
             VxS1(c, SpiRB(eR,iR)),
             VxS2(c, SpiRB(eR,iR)));
}

ComplexFun ChainB3(SpiSpec(iL,eL), cint a, cint b,
    cint c, SpiSpec(eR,iR)) {
  return SxS(SxBxV1(SpiLV(iL,eL), a, b),
             SxBxV2(SpiLV(iL,eL), a, b),
             BxS1(c, SpiRV(eR,iR)),
             BxS2(c, SpiRV(eR,iR)));
}

ComplexFun ChainV4(SpiSpec(iL,eL), cint a, cint b,
    cint c, cint d, SpiSpec(eR,iR)) {
  return SxS(SxVxB1(SpiLB(iL,eL), a, b),
             SxVxB2(SpiLB(iL,eL), a, b),
             VxBxS1(c, d, SpiRV(eR,iR)),
             VxBxS2(c, d, SpiRV(eR,iR)));
}

ComplexFun ChainB4(SpiSpec(iL,eL), cint a, cint b,
    cint c, cint d, SpiSpec(eR,iR)) {
  return SxS(SxBxV1(SpiLV(iL,eL), a, b),
             SxBxV2(SpiLV(iL,eL), a, b),
             BxVxS1(c, d, SpiRB(eR,iR)),
             BxVxS2(c, d, SpiRB(eR,iR)));
}

ComplexFun ChainV5(SpiSpec(iL,eL), cint a, cint b, cint c,
    cint d, cint e, SpiSpec(eR,iR)) {
  return SxS(SxVxBxV1(SpiLB(iL,eL), a, b, c),
             SxVxBxV2(SpiLB(iL,eL), a, b, c),
             BxVxS1(d, e, SpiRB(eR,iR)),
             BxVxS2(d, e, SpiRB(eR,iR)));
}

ComplexFun ChainB5(SpiSpec(iL,eL), cint a, cint b, cint c,
    cint d, cint e, SpiSpec(eR,iR)) {
  return SxS(SxBxVxB1(SpiLV(iL,eL), a, b, c),
             SxBxVxB2(SpiLV(iL,eL), a, b, c),
             VxBxS1(d, e, SpiRV(eR,iR)),
             VxBxS2(d, e, SpiRV(eR,iR)));
}

ComplexFun ChainV6(SpiSpec(iL,eL), cint a, cint b, cint c,
    cint d, cint e, cint f, SpiSpec(eR,iR)) {
  return SxS(SxVxBxV1(SpiLB(iL,eL), a, b, c),
             SxVxBxV2(SpiLB(iL,eL), a, b, c),
             BxVxBxS1(d, e, f, SpiRV(eR,iR)),
             BxVxBxS2(d, e, f, SpiRV(eR,iR)));
}

ComplexFun ChainB6(SpiSpec(iL,eL), cint a, cint b, cint c,
    cint d, cint e, cint f, SpiSpec(eR,iR)) {
  return SxS(SxBxVxB1(SpiLV(iL,eL), a, b, c),
             SxBxVxB2(SpiLV(iL,eL), a, b, c),
             VxBxVxS1(d, e, f, SpiRB(eR,iR)),
             VxBxVxS2(d, e, f, SpiRB(eR,iR)));
}

static inline int IndexDelta(cint a, cint b) {
  enum { n = 8*sizeof a - 1 };
  return (unsigned)((a ^ b) - 1) >> n;
}

static inline int IndexSign(cint a) {
  enum { n = 8*sizeof a - 1 };
  return (a >> n) | (((unsigned)-a) >> n);
}

static inline int IndexEps(cint a, cint b, cint c) {
  return IndexSign(a - b)*IndexSign(c - b)*IndexSign(a - c);
}

static inline RealType SqDiff(cRealType ma, cRealType mb) {
  return (ma - mb)*(ma + mb);
}

static inline RealType ThreeMom(cRealType sqrtS,
    cRealType ma, cRealType mb) {
  return sqrt(SqDiff(.5*(sqrtS - SqDiff(ma, mb)/sqrtS), mb));
}

#endif

