#if 0
	util.h
	prototypes for the util functions
	this file is part of FormCalc
	SIMD functions by J.-N. Lang
	last modified 3 Nov 15 th
#endif


#ifndef LEGS
#define LEGS 1
#endif

struct {
  RealType eps;
  int n, seq[2];
} hsel_;

#define hseleps hsel_.eps
#define hseln hsel_.n
#define seq(i) hsel_.seq[i-1]

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

enum { nvec0 = 12 };

struct {
  ComplexType vec0[LEGS*nvec0-minvec+1][2][2];
} vec0_;

#define vec0(i,j,n) vec0_.vec0[n-minvec][j-1][i-1]

struct {
  integer Hel0[LEGS];
} hel0_;

#define Hel0(i) hel0_.Hel0[i-1]

#define k0(i) (1+nvec0*(i-1))
#define s0(i) (3+nvec0*(i-1))
#define e0(i) (3+nvec0*(i-1)+Hel0(i))
#define ec0(i) (3+nvec0*(i-1)-Hel0(i))
#define Spinor0(i,af,d) (af*2+d+7+nvec0*(i-1)+Hel0(i))

#define NaN (NAN + I*NAN)

#define QH 32LL
#define ldQK 8
#define QK (1<<ldQK)

#define MomEncoding(f,i) ((integer8)((f) & (QK-1)) << (i-1)*ldQK)

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
#ifdef __clang__
typedef integer HelInt __attribute__((ext_vector_type(SIMD)));
#else
typedef integer HelInt __attribute__((__vector_size__(SIMD*sizeof(integer))));
#endif
typedef const HelInt cHelInt;

#define ResZero _mm_setzero_pd()
#define HelZero _mm256_setzero_pd()
#define HelNaN {NAN, NAN, NAN, NAN}


#define StoH(a) (sizeof(a) == sizeof(ComplexType) ? CtoH(a) : RtoH(a))

static inline HelType CtoH(cComplexType a) {
  return (HelType){Re(a),Im(a), Re(a),Im(a)};
}

static inline HelType RtoH(cRealType a) {
  return (HelType){a,0, a,0};
}

static inline HelType ItoH(cHelInt a) {
  return (HelType){a[0],0, a[1],0};
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
  return HxH(CtoH(a), b);
}

static inline HelType RxH(cRealType a, cHelType b) {
  return (HelType){a,a, a,a}*b;
}

static inline HelType IxH(cHelInt a, cHelType b) {
  return (HelType){a[0],a[0], a[1],a[1]}*b;
}

#define SxI(a,b) (sizeof(a) == sizeof(ComplexType) ? CxI(a,b) : RxI(a,b))

static inline HelType CxI(cComplexType a, cHelInt b) {
  return CtoH(a)*(HelType){b[0],b[0], b[1],b[1]};
}

static inline HelType RxI(cRealType a, cHelInt b) {
  return RtoH(a)*(HelType){b[0],b[0], b[1],b[1]};
}


static inline HelType ConjugateH(cHelType a) {
  return (HelType){1,-1, 1,-1}*a;
}

static inline ResType ReH(cHelType a) {
  return (ResType){a[0], a[2]};
}

static inline ResType AbsRxRes(cRealType a, cResType b) {
  return (ResType){fabs(a*b[0]), fabs(a*b[1])};
}

static inline RealType HelSum(cResType a) {
  return _mm_hadd_pd(a, a)[0];
}

#elif SIMD == 1 && defined __SSE3__

typedef RealType ResType;
typedef const ResType cResType;
typedef __m128d HelType;
typedef const HelType cHelType;
typedef int HelInt;
typedef const HelInt cHelInt;

#define HelZero _mm_setzero_pd()
#define ResZero 0
#define HelNaN {NAN, NAN}


#define StoH(a) (sizeof(a) == sizeof(ComplexType) ? CtoH(a) : RtoH(a))

static inline HelType CtoH(cComplexType a) {
  return (HelType){Re(a),Im(a)};
}

static inline HelType RtoH(cRealType a) {
  return (HelType){a,0};
}

static inline HelType ItoH(cHelInt a) {
  return (HelType){a,0};
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
  return HxH(CtoH(a), b);
}

static inline HelType RxH(cRealType a, cHelType b) {
  return (HelType){a,a}*b;
}

static inline HelType IxH(cHelInt a, cHelType b) {
  return (HelType){a,a}*b;
}


#define SxI(a,b) (sizeof(a) == sizeof(ComplexType) ? CxI(a,b) : RxI(a,b))

static inline HelType CxI(cComplexType a, cHelInt b) {
  return CtoH(a*b);
}

static inline HelType RxI(cRealType a, cHelInt b) {
  return RtoH(a*b);
}

static inline HelType ConjugateH(cHelType a) {
  return (HelType){1,-1}*a;
}

static inline ResType ReH(cHelType a) {
  return a[0];
}

#define AbsRxRes(a,b) fabs((a)*(b))
#define HelSum(a) (a)

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
#define HelNaN NaN

#define StoH(a) (a)
#define ItoH(a) (a)
#define HxH(a,b) (a)*(b)
#define SxH(a,b) (a)*(b)
#define RxH(a,b) (a)*(b)
#define IxH(a,b) (a)*(b)
#define SxI(a,b) (a)*(b)
#define ConjugateH(a) Conjugate(a)
#define ReH(a) Re(a)
#define AbsRxRes(a,b) fabs((a)*(b))
#define HelSum(a) (a)

#endif

enum { nvec = 8 };

struct {
  HelType vec[LEGS*nvec-minvec+1][2][2];
} vec_;

#define vec(i,j,n) vec_.vec[n-minvec][j-1][i-1]

struct {
  HelInt Hel[LEGS];
} hel_;

#define Hel(i) hel_.Hel[i-1]

#define k(i) (1+nvec*(i-1))
#define s(i) (2+nvec*(i-1))
#define e(i) (3+nvec*(i-1))
#define ec(i) (4+nvec*(i-1))
#define Spinor(i,af,d) (af+d+5+nvec*(i-1))

#define Vec(x,y,i) vec(x,y,i)
#define bVec vec_

#if NOUNDERSCORE
#define veccopy_ veccopy
#endif

extern void veccopy_(cinteger *v, cinteger *n);

#else

typedef RealType ResType;
typedef const ResType cResType;
typedef ComplexType HelType;
typedef const HelType cHelType;
#define HelZero 0
#define ResZero 0
#define HelNaN NaN

#define StoH(a) (a)
#define ItoH(a) (a)
#define HxH(a,b) (a)*(b)
#define SxH(a,b) (a)*(b)
#define RxH(a,b) (a)*(b)
#define IxH(a,b) (a)*(b)
#define SxI(a,b) (a)*(b)
#define ConjugateH(a) Conjugate(a)
#define ReH(a) Re(a)
#define AbsRxRes(a,b) fabs((a)*(b))
#define HelSum(a) (a)

#define k k0
#define s s0
#define e e0
#define ec ec0
#define Spinor Spinor0
#define Hel Hel0

#define Vec(x,y,i) vec0(x,y,i)
#define bVec vec0_

#endif

#define DEBr " %.13lg"
#define DEBc " (%.13lg,%.13lg)"
#define DEBx(x,i) ((RealType *)&x)[i]
#define DEB(a,x) switch( sizeof(x) ) { \
case sizeof(int): \
  printf(a " %d\n", *((int *)&x)); \
  break; \
case sizeof(RealType): \
  printf(a DEBr "\n", DEBx(x,0)); \
  break; \
case sizeof(ComplexType): \
  printf(a DEBc "\n", DEBx(x,0),DEBx(x,1)); \
  break; \
case 2*sizeof(ComplexType): \
  printf(a DEBc DEBc "\n", DEBx(x,0),DEBx(x,1), DEBx(x,2),DEBx(x,3)); \
  break; \
default: \
  printf("Don't know how to print " #x "\n"); \
}

#define LOOP(var,from,to,step) for( var = from; var <= to; var += step ) {
#define ENDLOOP(var) }
#define TEST(i,b) if( *(i) & (1 << (b)) ) {
#define ENDTEST(i,b) }

#define BIT_SETMASS 0
#define BIT_RESET 1
#define BIT_LOOP 2

#define ARG_ID(i,x) x
#define ARG_RE(i,x) Re(x)
#define ARG_HEL(i,x) (1LL << (Hel0(i)+2))
#define JOIN_SEQ(a,b) a,b
#define JOIN_MUL(a,b) a*b
#define JOIN_OCT(a,b) b+8*(a)
#define JOIN_DEC(a,b) b+10*(a)
#define JOIN_HEL(a,b) b+QH*(a)

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

