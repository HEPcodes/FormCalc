#if 0
	util.h
	prototypes for the util functions
	this file is part of FormCalc
	SIMD functions by J.-N. Lang
	last modified 8 Sep 14 th
#endif


#ifndef LEGS
#define LEGS 1
#endif

struct {
  RealType eps;
  int n;
} hsel_;

#define hseleps hsel_.eps
#define hseln hsel_.n

/* special vectors needed by num.h; overlaps intended */
enum {
  vTnj = 0,
  v0nj = -2,
  v1nj = -3,
  v2nj = -2,
  v3nj = -1,
  v4nj = 0,
  q1 = 0,
  minvec = -3 };

enum { nvec = 12 };

struct {
  ComplexType vec[LEGS*nvec-minvec+1][2][2];
} vec_;

#define vec(i,j,n) vec_.vec[n-minvec][j-1][i-1]

#define k0(i) (1+nvec*(i-1))
#define s0(i) (3+nvec*(i-1))
#define e0(i) (3+nvec*(i-1)+Hel(i))
#define ec0(i) (3+nvec*(i-1)-Hel(i))
#define Spinor0(i,af,d) (af*2+d+7+nvec*(i-1)+Hel(i))

#define SIMD_ONLY(x)
#define SIMD_MULT(x)
#define SIMD_CEIL(n) n

#if SIMD > 0 && !defined DEPS

#include <immintrin.h>

#undef SIMD_ONLY
#define SIMD_ONLY(x) x

#if SIMD == 2 && defined __AVX__

#undef SIMD_MULT
#define SIMD_MULT(x) x
#undef SIMD_CEIL
#define SIMD_CEIL(n) (n+SIMD-1)/SIMD

typedef __m128d ResType;
typedef const ResType cResType;
typedef __m256d HelType;
typedef const HelType cHelType;
#define HelZero _mm256_setzero_pd()
#define ResZero _mm_setzero_pd()

static inline HelType ToH(cComplexType a) {
  return _mm256_broadcast_pd((__m128d *)&a);
}

static inline HelType HxH(cHelType a, cHelType b) {
  HelType t1 = _mm256_movedup_pd(b);
  HelType t2 = _mm256_shuffle_pd(b, b, 0xF);
  t1 = _mm256_mul_pd(a, t1);
  t2 = _mm256_mul_pd(a, t2);
  t2 = _mm256_shuffle_pd(t2, t2, 0x5);
  return _mm256_addsub_pd(t1, t2);
}

#define SxH(a,b) (sizeof(a) == sizeof(ComplexType) ? CxH(a,b) : RxH(a,b))

static inline HelType CxH(cComplexType a, cHelType b) {
  return HxH(_mm256_broadcast_pd((__m128d *)&a), b);
}

static inline HelType RxH(cRealType a, cHelType b) {
  return _mm256_broadcast_sd(&a)*b;
}

static inline ResType ReHcH(cHelType a, cHelType b) {
  HelType ab = a*b;
  HelType t = _mm256_hadd_pd(ab, ab);
  return (ResType){t[0], t[2]};
}

static inline RealType HelSum(cResType a) {
  return _mm_hadd_pd(a, a)[0];
}

static inline ResType ResAbs(cResType a) {
  return (ResType){fabs(a[0]), fabs(a[1])};
}

static inline ResType ToRes(cRealType a) {
  return _mm_load1_pd(&a);
}

#elif SIMD == 1 && defined __SSE3__

typedef RealType ResType;
typedef const ResType cResType;
typedef __m128d HelType;
typedef const HelType cHelType;
#define HelZero _mm_setzero_pd()
#define ResZero 0
#define HelSum(a) (a)
#define ResAbs(a) fabs(a)
#define ToRes(a) (a)

static inline HelType ToH(cComplexType a) {
  return *(HelType *)&a;
}

static inline HelType HxH(cHelType a, cHelType b) {
  HelType t1 = _mm_movedup_pd(b);
  HelType t2 = _mm_shuffle_pd(b, b, 0x3);
  t1 = _mm_mul_pd(a, t1);
  t2 = _mm_mul_pd(a, t2);
  t2 = _mm_shuffle_pd(t2, t2, 0x1);
  return _mm_addsub_pd(t1, t2);
}

#define SxH(a,b) (sizeof(a) == sizeof(ComplexType) ? CxH(a,b) : RxH(a,b))

static inline HelType CxH(cComplexType a, cHelType b) {
  return HxH(*(HelType *)&a, b);
}

static inline HelType RxH(cRealType a, cHelType b) {
  return _mm_load1_pd(&a)*b;
}

static inline ResType ReHcH(cHelType a, cHelType b) {
  HelType ab = a*b;
  return _mm_hadd_pd(ab, ab)[0];
}

#else

#if SIMD != 1
#error This value of SIMD not implemented.
#endif

typedef RealType ResType;
typedef const ResType cResType;
typedef ComplexType HelType;
typedef const HelType cHelType;
#define HelZero 0
#define ResZero 0

#define ToH(a) (a)
#define HxH(a,b) (a)*(b)
#define SxH(a,b) (a)*(b)
#define CxH(a,b) (a)*(b)
#define RxH(a,b) (a)*(b)
#define ReHcH(a,b) Re(Conjugate(a)*(b))
#define HelSum(a) (a)
#define ResAbs(a) fabs(a)
#define ToRes(a) (a)

#endif

enum { nves = 8 };

struct {
  HelType ves[LEGS*nves+1][2][2];
} ves_;

#define ves(i,j,n) ves_.ves[n][j-1][i-1]


#define k(i) (1+nves*(i-1))
#define s(i) (2+nves*(i-1))
#define e(i) (3+nves*(i-1))
#define ec(i) (4+nves*(i-1))
#define Spinor(i,af,d) (af+d+5+nves*(i-1))

#define Vec(x,y,i) ves(x,y,i)
#define bVec ves_

#if NOUNDERSCORE
#define veccopy_ veccopy
#endif

extern void veccopy_(cinteger *v, cinteger *n, cinteger *hel);

#else

typedef RealType ResType;
typedef const ResType cResType;
typedef ComplexType HelType;
typedef const HelType cHelType;
#define HelZero 0
#define ResZero 0

#define ToH(a) (a)
#define HxH(a,b) (a)*(b)
#define SxH(a,b) (a)*(b)
#define CxH(a,b) (a)*(b)
#define RxH(a,b) (a)*(b)
#define ReHcH(a,b) Re(Conjugate(a)*(b))
#define HelSum(a) (a)
#define ResAbs(a) fabs(a)
#define ToRes(a) (a)

#define k k0
#define s s0
#define e e0
#define ec ec0
#define Spinor Spinor0
#define Vec(x,y,i) vec(x,y,i)
#define bVec vec_

#endif

#define DEB(a,x) printf(a " (%.13lg,%.13lg)\n", Re(x), Im(x))
#define LOOP(var,from,to,step) for( var = from; var <= to; var += step ) {
#define ENDLOOP(var) }
#define TEST(i,b) if( *(i) & (1 << (b)) ) {
#define ENDTEST(i,b) }

#define BIT_RESET 0
#define BIT_LOOP 1
#define MASK_HEL(i) (1ULL << (5*(LEGS-i)+Hel(i)+2))

#define INI_S(seq) clearcache()
#define INI_ANGLE(seq) markcache()
#define DEINI(seq) restorecache()

#if PARALLEL
#if NOUNDERSCORE
#define sqmeprep_ sqmeprep
#define sqmeexec_ sqmeexec
#define sqmesync_ sqmesync
#define sqmewait_ sqmewait
#endif

void sqmeprep_(void *v, void *ve, void *r, void *re,
  void *s, void *se, void *a, void *ae, void *h, void *he);
void sqmeexec_(void (*foo)(ResType *, cinteger *),
  ResType *res, cinteger *flags);
void sqmesync_();
void sqmewait_();

#define lim(x) &x, ((char *)&x + sizeof(x))
#define PAR_PREP(r, s, a, h) sqmeprep_(lim(bVec), lim(r), lim(s), lim(a), lim(h))
#define PAR_EXEC(f, res, flags) sqmeexec_(&f, res, flags)
#define PAR_SYNC() sqmesync_()
#else
#define PAR_PREP(r, s, a, h)
#define PAR_EXEC(f, res, flags) f(res, flags)
#define PAR_SYNC()
#endif

