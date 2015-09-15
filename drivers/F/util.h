* util.h
* prototypes for the util functions
* this file is part of FormCalc
* last modified 10 Jun 15 th


#ifndef UTIL_H
#define UTIL_H

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

#define k0(i) (1+nvec*(i-1))
#define s0(i) (3+nvec*(i-1))
#define e0(i) (3+nvec*(i-1)+Hel(i))
#define ec0(i) (3+nvec*(i-1)-Hel(i))
#define Spinor0(i,af,d) (af*2+d+7+nvec*(i-1)+Hel(i))

#define Finite ishft(1, -epsi)
#define PVC(i) i+epsi

#if SIMD > 0

#define ResType RealType, dimension(SIMD) ::
#define HelType ComplexType, dimension(SIMD) ::
#define HelDim(i) SIMD,i
#define HelInd(i) :,i
#define HelSum(x) sum(x)

#define k(i) (1+nves*(i-1))
#define s(i) (2+nves*(i-1))
#define e(i) (3+nves*(i-1))
#define ec(i) (4+nves*(i-1))
#define Spinor(i,af,d) (af+d+5+nves*(i-1))

#define Vec(x,y,i) ves(:,x,y,i)
#define bVec ves(1,1,1,1)
#define eVec ves_end

#define SIMD_CEIL(n) (n+SIMD-1)/SIMD
#define SIMD_ONLY(x) x
#if SIMD > 1
#define SIMD_MULT(x) x
#else
#define SIMD_MULT(x)
#endif

#else

#define ResType RealType
#define HelType ComplexType
#define HelDim(i) i
#define HelInd(i) i
#define HelSum(x) x

#define k k0
#define s s0
#define e e0
#define ec ec0
#define Spinor Spinor0

#define Vec(x,y,i) vec(x,y,i)
#define bVec vec(1,1,1)
#define eVec vec_end

#define SIMD_CEIL(n) n
#define SIMD_ONLY(x)
#define SIMD_MULT(x)
#endif

#define MomEncoding(f,i) iand(f,QK-1)*QK**(i-1)

#define Digit(i) char(i+48)
#define Polar(r,theta) r*exp(cI*degree*theta)

#define Error(err,msg) call m_(err, __LINE__, __FILE__, msg)
#define Warning(msg) call m_(0, 0, __FILE__, msg)
#define INFO print *,
#define DEB(a,x) print *, a, x
#define LOOP(var,from,to,step) do var = from, to, step
#define ENDLOOP(var) enddo
#define TEST(i,b) if( btest(i,b) ) then
#define ENDTEST(i,b) endif

#define BIT_SETMASS 0
#define BIT_RESET 1
#define BIT_LOOP 2

#define ARG_ID(i,x) x
#define ARG_RE(i,x) Re(x)
#define ARG_HEL(i,x) ibset(0,Hel(i)+2)
#define JOIN_SEQ(a,b) a,b
#define JOIN_MUL(a,b) a*b
#define JOIN_OCT(a,b) b+8*(a)
#define JOIN_DEC(a,b) b+10*(a)
#define JOIN_HEL(a,b) b+QH*(a)

#define INI_S(seq) call clearcache
#define INI_ANGLE(seq) call markcache
#define DEINI(seq) call restorecache

#if PARALLEL
#define PAR_PREP(r,re, s,se, a,ae, h,he) call sqmeprep(bVec,eVec, r,re, s,se, a,ae, h,he)
#define PAR_EXEC(f, res, flags) call sqmeexec(f, res, flags)
#define PAR_SYNC() call sqmesync()
#else
#define PAR_PREP(r,re, s,se, a,ae, h,he)
#define PAR_EXEC(f, res, flags) call f(res, flags)
#define PAR_SYNC()
#endif

#define Cut(c,m) (m)*(c)

#define CUT_MIN 1
#define CUT_MAX 2

#define CUT_COSTH 4
#define CUT_COSTHCMS 16
#define CUT_COSTH_E 64
#define CUT_COSTH_K 65
#define CUT_MREM 256
#define CUT_MREM_E 1024
#define CUT_MREM_K 1025
#define CUT_MREM_ET 4096
#define CUT_MREM_KT 4097
#define CUT_MREM_RAP 16384
#define CUT_MREM_PRAP 16385

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

	integer nvec
	parameter (nvec = 12)

	ComplexType vec(2,2,minvec:nvec*LEGS), vec_end
	common /vec/ vec, vec_end

	RealType hseleps
	integer hseln
	common /hsel/ hseleps, hseln

#if SIMD > 0
	integer nves
	parameter (nves = 8)

	HelType ves(HelDim(2),2,minvec:nves*LEGS), ves_end
	common /ves/ ves, ves_end
#endif

	RealType momspec(12,LEGS)
	common /momspec/ momspec

* QH (QK) is encoding base for helicities (momenta)
	integer ldQH
	integer*8 QH, QK
	parameter (ldQH = 5, QH = 2**ldQH, QK = 256)
#endif

