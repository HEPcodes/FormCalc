#if 0
	inline.h
	inline versions of the util functions
	this file is part of FormCalc
	last modified 8 Jun 15 th
#endif


#ifndef INLINE_H
#define INLINE_H

#define HelFun static inline HelType
#define SpiType(s1,s2) cHelType s1, cHelType s2
#define SpiSpec(i,e) cinteger i, cinteger e
#define SpiLV(i,e) RxH(1-2*e,Vec(1+e,1+e,i)), Vec(2-e,1+e,i)
#define SpiLB(i,e) RxH(1-2*e,Vec(1+e,2-e,i)), Vec(2-e,2-e,i)
#define SpiRV(e,i) Vec(1+e,1+e,i), RxH(1-2*e,Vec(2-e,1+e,i))
#define SpiRB(e,i) Vec(1+e,2-e,i), RxH(1-2*e,Vec(2-e,2-e,i))

#define Sqrt sqrt

#define Power(x, n) ({ \
  typeof(x) _x = x, _r = 1; \
  typeof(n) _n = n; \
  while( _n ) { \
    if( _n & 1 ) _r *= _x; \
    _n >>= 1; \
    _x *= _x; \
  } \
  _r; \
})

static inline RealType Sq(cComplexType c) {
  return Re(c*Conjugate(c));
}

static inline ComplexType SInvariant(cinteger a, cinteger b) {
  return (Re(vec(1,1,k0(a))) + Re(vec(1,1,k0(b))))*
         (Re(vec(2,2,k0(a))) + Re(vec(2,2,k0(b)))) -
         Sq(vec(1,2,k0(a)) + vec(1,2,k0(b)));
}

static inline ComplexType TInvariant(cinteger a, cinteger b) {
  return (Re(vec(1,1,k0(a))) - Re(vec(1,1,k0(b))))*
         (Re(vec(2,2,k0(a))) - Re(vec(2,2,k0(b)))) -
         Sq(vec(1,2,k0(a)) - vec(1,2,k0(b)));
}

HelFun Pair(cinteger a, cinteger b) {
  return RxH(.5,
    HxH(Vec(1,1,a), Vec(2,2,b)) + HxH(Vec(2,2,a), Vec(1,1,b)) -
    HxH(Vec(1,2,a), Vec(2,1,b)) - HxH(Vec(2,1,a), Vec(1,2,b)));
}

HelFun Eps_(cinteger a, cinteger b, cinteger c, cinteger d) {
  return HxH(HxH(Vec(1,1,a), Vec(2,2,b)) -
             HxH(Vec(2,2,a), Vec(1,1,b)),
             HxH(Vec(2,1,c), Vec(1,2,d)) -
             HxH(Vec(1,2,c), Vec(2,1,d)));
}

HelFun Eps(cinteger a, cinteger b, cinteger c, cinteger d) {
  return RxH(.25,
    Eps_(a, b, c, d) + Eps_(c, d, a, b) -
    Eps_(a, c, b, d) - Eps_(b, d, a, c) +
    Eps_(a, d, b, c) + Eps_(b, c, a, d));
}

HelFun SxS(SpiType(l1,l2), SpiType(r1,r2)) {
  return HxH(l1, r1) + HxH(l2, r2);
}

HelFun SxV1(SpiType(l1,l2), cinteger a) {
  return HxH(l1, Vec(1,1,a)) + HxH(l2, Vec(2,1,a));
}

HelFun SxV2(SpiType(l1,l2), cinteger a) {
  return HxH(l2, Vec(2,2,a)) + HxH(l1, Vec(1,2,a));
}

HelFun SxB1(SpiType(l1,l2), cinteger a) {
  return HxH(l1, Vec(2,2,a)) - HxH(l2, Vec(2,1,a));
}

HelFun SxB2(SpiType(l1,l2), cinteger a) {
  return HxH(l2, Vec(1,1,a)) - HxH(l1, Vec(1,2,a));
}

HelFun VxS1(cinteger a, SpiType(r1,r2)) {
  return HxH(Vec(1,1,a), r1) + HxH(Vec(1,2,a), r2);
}

HelFun VxS2(cinteger a, SpiType(r1,r2)) {
  return HxH(Vec(2,1,a), r1) + HxH(Vec(2,2,a), r2);
}

HelFun BxS1(cinteger a, SpiType(r1,r2)) {
  return HxH(Vec(2,2,a), r1) - HxH(Vec(1,2,a), r2);
}

HelFun BxS2(cinteger a, SpiType(r1,r2)) {
  return HxH(Vec(1,1,a), r2) - HxH(Vec(2,1,a), r1);
}

HelFun SxVxB1(SpiType(l1,l2), cinteger a, cinteger b) {
  return SxB1(SxV1(l1,l2, a),SxV2(l1,l2, a), b);
}

HelFun SxVxB2(SpiType(l1,l2), cinteger a, cinteger b) {
  return SxB2(SxV1(l1,l2, a),SxV2(l1,l2, a), b);
}

HelFun SxBxV1(SpiType(l1,l2), cinteger a, cinteger b) {
  return SxV1(SxB1(l1,l2, a),SxB2(l1,l2, a), b);
}

HelFun SxBxV2(SpiType(l1,l2), cinteger a, cinteger b) {
  return SxV2(SxB1(l1,l2, a),SxB2(l1,l2, a), b);
}

HelFun BxVxS1(cinteger b, cinteger a, SpiType(r1,r2)) {
  return BxS1(b, VxS1(a, r1,r2),VxS2(a, r1,r2));
}

HelFun BxVxS2(cinteger b, cinteger a, SpiType(r1,r2)) {
  return BxS2(b, VxS1(a, r1,r2),VxS2(a, r1,r2));
}

HelFun VxBxS1(cinteger b, cinteger a, SpiType(r1,r2)) {
  return VxS1(b, BxS1(a, r1,r2),BxS2(a, r1,r2));
}

HelFun VxBxS2(cinteger b, cinteger a, SpiType(r1,r2)) {
  return VxS2(b, BxS1(a, r1,r2),BxS2(a, r1,r2));
}

HelFun SxVxBxV1(SpiType(l1,l2), cinteger a, cinteger b, cinteger c) {
  return SxBxV1(SxV1(l1,l2, a),SxV2(l1,l2, a), b, c);
}

HelFun SxVxBxV2(SpiType(l1,l2), cinteger a, cinteger b, cinteger c) {
  return SxBxV2(SxV1(l1,l2, a),SxV2(l1,l2, a), b, c);
}

HelFun SxBxVxB1(SpiType(l1,l2), cinteger a, cinteger b, cinteger c) {
  return SxVxB1(SxB1(l1,l2, a),SxB2(l1,l2, a), b, c);
}

HelFun SxBxVxB2(SpiType(l1,l2), cinteger a, cinteger b, cinteger c) {
  return SxVxB2(SxB1(l1,l2, a),SxB2(l1,l2, a), b, c);
}

HelFun VxBxVxS1(cinteger c, cinteger b, cinteger a, SpiType(r1,r2)) {
  return VxBxS1(c, b, VxS1(a, r1,r2),VxS2(a, r1,r2));
}

HelFun VxBxVxS2(cinteger c, cinteger b, cinteger a, SpiType(r1,r2)) {
  return VxBxS2(c, b, VxS1(a, r1,r2),VxS2(a, r1,r2));
}

HelFun BxVxBxS1(cinteger c, cinteger b, cinteger a, SpiType(r1,r2)) {
  return BxVxS1(c, b, BxS1(a, r1,r2),BxS2(a, r1,r2));
}

HelFun BxVxBxS2(cinteger c, cinteger b, cinteger a, SpiType(r1,r2)) {
  return BxVxS2(c, b, BxS1(a, r1,r2),BxS2(a, r1,r2));
}

HelFun SxVxBxVxB1(SpiType(l1,l2), cinteger a, cinteger b, cinteger c, cinteger d) {
  return SxBxVxB1(SxV1(l1,l2, a),SxV2(l1,l2, a), b, c, d);
}

HelFun SxVxBxVxB2(SpiType(l1,l2), cinteger a, cinteger b, cinteger c, cinteger d) {
  return SxBxVxB2(SxV1(l1,l2, a),SxV2(l1,l2, a), b, c, d);
}

HelFun SxBxVxBxV1(SpiType(l1,l2), cinteger a, cinteger b, cinteger c, cinteger d) {
  return SxVxBxV1(SxB1(l1,l2, a),SxB2(l1,l2, a), b, c, d);
}

HelFun SxBxVxBxV2(SpiType(l1,l2), cinteger a, cinteger b, cinteger c, cinteger d) {
  return SxVxBxV2(SxB1(l1,l2, a),SxB2(l1,l2, a), b, c, d);
}

HelFun BxVxBxVxS1(cinteger d, cinteger c, cinteger b, cinteger a, SpiType(r1,r2)) {
  return BxVxBxS1(d, c, b, VxS1(a, r1,r2),VxS2(a, r1,r2));
}

HelFun BxVxBxVxS2(cinteger d, cinteger c, cinteger b, cinteger a, SpiType(r1,r2)) {
  return BxVxBxS2(d, c, b, VxS1(a, r1,r2),VxS2(a, r1,r2));
}

HelFun VxBxVxBxS1(cinteger d, cinteger c, cinteger b, cinteger a, SpiType(r1,r2)) {
  return VxBxVxS1(d, c, b, BxS1(a, r1,r2),BxS2(a, r1,r2));
}

HelFun VxBxVxBxS2(cinteger d, cinteger c, cinteger b, cinteger a, SpiType(r1,r2)) {
  return VxBxVxS2(d, c, b, BxS1(a, r1,r2),BxS2(a, r1,r2));
}

HelFun ChainV0(SpiSpec(iL,eL), SpiSpec(eR,iR)) {
  return SxS(SpiLB(iL,eL), SpiRV(eR,iR));
}

HelFun ChainB0(SpiSpec(iL,eL), SpiSpec(eR,iR)) {
  return SxS(SpiLV(iL,eL),
             SpiRB(eR,iR));
}

HelFun ChainV1(SpiSpec(iL,eL), cinteger a, SpiSpec(eR,iR)) {
  return SxS(SxV1(SpiLB(iL,eL), a),
             SxV2(SpiLB(iL,eL), a),
             SpiRB(eR,iR));
}

HelFun ChainB1(SpiSpec(iL,eL), cinteger a, SpiSpec(eR,iR)) {
  return SxS(SxB1(SpiLV(iL,eL), a),
             SxB2(SpiLV(iL,eL), a),
             SpiRV(eR,iR));
}

HelFun ChainV2(SpiSpec(iL,eL), cinteger a, cinteger b, SpiSpec(eR,iR)) {
  return SxS(SxV1(SpiLB(iL,eL), a),
             SxV2(SpiLB(iL,eL), a),
             BxS1(b, SpiRV(eR,iR)),
             BxS2(b, SpiRV(eR,iR)));
}

HelFun ChainB2(SpiSpec(iL,eL), cinteger a, cinteger b, SpiSpec(eR,iR)) {
  return SxS(SxB1(SpiLV(iL,eL), a),
             SxB2(SpiLV(iL,eL), a),
             VxS1(b, SpiRB(eR,iR)),
             VxS2(b, SpiRB(eR,iR)));
}

HelFun ChainV3(SpiSpec(iL,eL), cinteger a, cinteger b,
    cinteger c, SpiSpec(eR,iR)) {
  return SxS(SxVxB1(SpiLB(iL,eL), a, b),
             SxVxB2(SpiLB(iL,eL), a, b),
             VxS1(c, SpiRB(eR,iR)),
             VxS2(c, SpiRB(eR,iR)));
}

HelFun ChainB3(SpiSpec(iL,eL), cinteger a, cinteger b,
    cinteger c, SpiSpec(eR,iR)) {
  return SxS(SxBxV1(SpiLV(iL,eL), a, b),
             SxBxV2(SpiLV(iL,eL), a, b),
             BxS1(c, SpiRV(eR,iR)),
             BxS2(c, SpiRV(eR,iR)));
}

HelFun ChainV4(SpiSpec(iL,eL), cinteger a, cinteger b,
    cinteger c, cinteger d, SpiSpec(eR,iR)) {
  return SxS(SxVxB1(SpiLB(iL,eL), a, b),
             SxVxB2(SpiLB(iL,eL), a, b),
             VxBxS1(c, d, SpiRV(eR,iR)),
             VxBxS2(c, d, SpiRV(eR,iR)));
}

HelFun ChainB4(SpiSpec(iL,eL), cinteger a, cinteger b,
    cinteger c, cinteger d, SpiSpec(eR,iR)) {
  return SxS(SxBxV1(SpiLV(iL,eL), a, b),
             SxBxV2(SpiLV(iL,eL), a, b),
             BxVxS1(c, d, SpiRB(eR,iR)),
             BxVxS2(c, d, SpiRB(eR,iR)));
}

HelFun ChainV5(SpiSpec(iL,eL), cinteger a, cinteger b, cinteger c,
    cinteger d, cinteger e, SpiSpec(eR,iR)) {
  return SxS(SxVxBxV1(SpiLB(iL,eL), a, b, c),
             SxVxBxV2(SpiLB(iL,eL), a, b, c),
             BxVxS1(d, e, SpiRB(eR,iR)),
             BxVxS2(d, e, SpiRB(eR,iR)));
}

HelFun ChainB5(SpiSpec(iL,eL), cinteger a, cinteger b, cinteger c,
    cinteger d, cinteger e, SpiSpec(eR,iR)) {
  return SxS(SxBxVxB1(SpiLV(iL,eL), a, b, c),
             SxBxVxB2(SpiLV(iL,eL), a, b, c),
             VxBxS1(d, e, SpiRV(eR,iR)),
             VxBxS2(d, e, SpiRV(eR,iR)));
}

HelFun ChainV6(SpiSpec(iL,eL), cinteger a, cinteger b, cinteger c,
    cinteger d, cinteger e, cinteger f, SpiSpec(eR,iR)) {
  return SxS(SxVxBxV1(SpiLB(iL,eL), a, b, c),
             SxVxBxV2(SpiLB(iL,eL), a, b, c),
             BxVxBxS1(d, e, f, SpiRV(eR,iR)),
             BxVxBxS2(d, e, f, SpiRV(eR,iR)));
}

HelFun ChainB6(SpiSpec(iL,eL), cinteger a, cinteger b, cinteger c,
    cinteger d, cinteger e, cinteger f, SpiSpec(eR,iR)) {
  return SxS(SxBxVxB1(SpiLV(iL,eL), a, b, c),
             SxBxVxB2(SpiLV(iL,eL), a, b, c),
             VxBxVxS1(d, e, f, SpiRB(eR,iR)),
             VxBxVxS2(d, e, f, SpiRB(eR,iR)));
}

HelFun ChainV7(SpiSpec(iL,eL), cinteger a, cinteger b, cinteger c,
    cinteger d, cinteger e, cinteger f, cinteger g, SpiSpec(eR,iR)) {
  return SxS(SxVxBxVxB1(SpiLB(iL,eL), a, b, c, d),
             SxVxBxVxB2(SpiLB(iL,eL), a, b, c, d),
             VxBxVxS1(e, f, g, SpiRB(eR,iR)),
             VxBxVxS2(e, f, g, SpiRB(eR,iR)));
}

HelFun ChainB7(SpiSpec(iL,eL), cinteger a, cinteger b, cinteger c,
    cinteger d, cinteger e, cinteger f, cinteger g, SpiSpec(eR,iR)) {
  return SxS(SxBxVxBxV1(SpiLV(iL,eL), a, b, c, d),
             SxBxVxBxV2(SpiLV(iL,eL), a, b, c, d),
             BxVxBxS1(e, f, g, SpiRV(eR,iR)),
             BxVxBxS2(e, f, g, SpiRV(eR,iR)));
}

HelFun ChainV8(SpiSpec(iL,eL), cinteger a, cinteger b, cinteger c,
    cinteger d, cinteger e, cinteger f, cinteger g,
    cinteger h, SpiSpec(eR,iR)) {
  return SxS(SxVxBxVxB1(SpiLB(iL,eL), a, b, c, d),
             SxVxBxVxB2(SpiLB(iL,eL), a, b, c, d),
             VxBxVxBxS1(e, f, g, h, SpiRV(eR,iR)),
             VxBxVxBxS2(e, f, g, h, SpiRV(eR,iR)));
}

HelFun ChainB8(SpiSpec(iL,eL), cinteger a, cinteger b, cinteger c,
    cinteger d, cinteger e, cinteger f, cinteger g,
    cinteger h, SpiSpec(eR,iR)) {
  return SxS(SxBxVxBxV1(SpiLV(iL,eL), a, b, c, d),
             SxBxVxBxV2(SpiLV(iL,eL), a, b, c, d),
             BxVxBxVxS1(e, f, g, h, SpiRB(eR,iR)),
             BxVxBxVxS2(e, f, g, h, SpiRB(eR,iR)));
}

static inline integer IndexDelta(cinteger a, cinteger b) {
  enum { n = 8*sizeof a - 1 };
  return (unsigned)((a ^ b) - 1) >> n;
}

static inline integer IndexSign(cinteger a) {
  enum { n = 8*sizeof a - 1 };
  return (a >> n) | (((unsigned)-a) >> n);
}

static inline integer IndexEps(cinteger a, cinteger b, cinteger c) {
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

