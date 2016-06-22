* util.h
* prototypes for the util functions
* this file is part of FormCalc
* last modified 7 Jun 16 th


#ifndef UTIL_H
#define UTIL_H

#define signbit(i) ibits(i,31,1)

#define k0(i) (1+nvec0*(i-1))
#define s0(i) (3+nvec0*(i-1))
#define e0(i) (3+nvec0*(i-1)+Hel0(i))
#define ec0(i) (3+nvec0*(i-1)-Hel0(i))
#define Spinor0(i,af,d) (af*2+d+7+nvec0*(i-1)+Hel0(i))

#define Finite ishft(1,-epsi)
#define Epsi(i) i+epsi

#define NaN(n) n*bogus

#if SIMD > 1

#define SIMD_ONLY(x) x
#define SIMD_CEIL(n) (n+SIMD-1)/SIMD

#define ResType RealType, dimension(SIMD) ::
#define HelType ComplexType, dimension(SIMD) ::
#define HelDim(i) SIMD,i
#define HelAll(i) :,i
#define HelInd(v,i) v,i
#define HelLoop(x,y,v,vmax) (x,y, v = 1,vmax)
#define HelSum(x) sum(x)

#if SIMD == 2
#define HelNaN(n) n*bogus,n*bogus
#elif SIMD == 4
#define HelNaN(n) n*bogus,n*bogus,n*bogus,n*bogus
#endif

#define k(i) (1+nvec*(i-1))
#define s(i) (2+nvec*(i-1))
#define e(i) (3+nvec*(i-1))
#define ec(i) (4+nvec*(i-1))
#define Spinor(i,af,d) (af+d+5+nvec*(i-1))

#define Vec(x,y,i) vec(:,x,y,i)
#define bVec vec(1,1,1,1)
#define eVec vec_end

#else

#define SIMD_CEIL(n) n
#define SIMD_ONLY(x)

#define ResType RealType
#define HelType ComplexType
#define HelDim(i) i
#define HelAll(i) i
#define HelInd(v,i) i
#define HelLoop(x,y,v,vmax) x,y
#define HelSum(x) x
#define HelNaN(n) n*bogus

#define k k0
#define s s0
#define e e0
#define ec ec0
#define Spinor Spinor0
#define Hel Hel0

#define Vec(x,y,i) vec0(x,y,i)
#define bVec vec0(1,1,1)
#define eVec vec0_end

#endif

#define MomEncoding(f,i) iand(f,QK-1)*QK**(i-1)

#define Digit(i) char(i+48)
#define Polar(r,theta) r*exp(cI*degree*theta)
#define SafeLog(x) log(max(x,1D-300))

#define Error(err,msg) call m_(err, __LINE__, __FILE__, msg)
#define Warning(msg) call m_(0, 0, __FILE__, msg)
#define INFO print *,
#define LOOP(var,from,to,step) do var = from, to, step
#define ENDLOOP(var) enddo
#define TEST(i,b) if( btest(i,b) ) then
#define ENDTEST(i,b) endif

#define BIT_SETMASS 0
#define BIT_RESET 1
#define BIT_LOOP 2

#define ARG_ID(i,x,o) x
#define ARG_RE(i,x,o) Re(x)
#define ARG_HEL(i,x,o) ibset(0,Hel0(i)+2)
#define JOIN_SEQ(a,b) a,b
#define JOIN_MUL(a,b) a*b
#define JOIN_OCT(a,b) b+8*(a)
#define JOIN_DEC(a,b) b+10*(a)
#define JOIN_HEL(a,b) b+QH*(a)

#define DEB(tag,var) print *, tag, " =", var

#define CHK_INI(seq) chkyes = ior(seq(1), seq(2))
#define CHK_PRE(var) chkval = abs(var)
#define CHK_POST(tag,var) if( chkyes .ne. 0 .and. abs(abs(var) - chkval) .gt. 1D10 ) print *, tag, " differs"

#define INI_S() call clearcache
#define INI_A() call markcache
#define DEINI() call restorecache

#define Var(v) var(1,v)
#define Show(v) var(2,v)
#define Lower(v) var(3,v)
#define Upper(v) var(4,v)
#define Step(v) var(5,v)
#define CutVar(c,v) var(c+5,v)
#define CutMin(v) var(6,v)
#define CutMax(v) var(7,v)

#define Cut(c,m) (m)*(c)

#define CUT_MIN 1
#define CUT_MAX 2

#define CUT_COSTH 4
#define CUT_MREM 16
#define CUT_MREM_E 256
#define CUT_MREM_K 257
#define CUT_MREM_ET 1024
#define CUT_MREM_KT 1025
#define CUT_MREM_RAP 4096
#define CUT_MREM_PRAP 4097

#define SPEC_M 1
#define SPEC_K 2
#define SPEC_E 3
#define SPEC_KT 4
#define SPEC_ET 5
#define SPEC_PRAP 6
#define SPEC_RAP 7
#define SPEC_DELTAK 8
#define SPEC_PHI 9
#define SPEC_EX 10
#define SPEC_EY 11
#define SPEC_EZ 12

#define SOFT 8265681
#define HARD 8265682

#define GAUSS 7378841
#define PATTERSON 7378842
#define VEGAS 7378843
#define SUAVE 7378844
#define DIVONNE 7378845
#define CUHRE 7378846

#else

#ifndef LEGS
#define LEGS 1
#endif

* special vectors needed by num.h; overlaps intended
	integer vTnj, v0nj, v1nj, v2nj, v3nj, v4nj, q1, minvec
	parameter (vTnj = 0)
	parameter (v0nj = -2)
	parameter (v1nj = -3)
	parameter (v2nj = -2)
	parameter (v3nj = -1)
	parameter (v4nj = 0)
	parameter (q1 = 0)
	parameter (minvec = -3)

	integer nvec0
	parameter (nvec0 = 12)

	ComplexType vec0(2,2,minvec:nvec0*LEGS), vec0_end
	common /vec0/ vec0, vec0_end

	integer Hel0(LEGS)
	common /hel0/ Hel0

	RealType hseleps
	integer hseln
	common /hsel/ hseleps, hseln

	RealType chkval
	integer chkyes
	common /chktmp/ chkval, chkyes

#if SIMD > 1
	integer nvec
	parameter (nvec = 8)

	HelType vec(HelDim(2),2,minvec:nvec*LEGS), vec_end
	common /vec/ vec, vec_end

	integer Hel(HelDim(LEGS))
	common /hel/ Hel
#endif

	RealType momspec(12,LEGS)
	common /momspec/ momspec

* QH (QK) is encoding base for helicities (momenta)
	integer ldQH
	integer*8 QH, QK
	parameter (ldQH = 5, QH = 2**ldQH, QK = 256)
#endif

