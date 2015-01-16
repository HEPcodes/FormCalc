/*
	util-c.h
	C prototypes for the functions in util.a
	this file is part of FormCalc
	last modified 20 Jul 12 th
*/


#ifndef LEGS
#define LEGS 1
#endif

enum { nvec = 10 };

struct {
  ComplexType vec[LEGS][nvec][2][2];
} vectors_;

/* encoding base for spinor chains (JC) and momenta (JK) */
enum {JC = 256LL, JK = 256LL};

	RealType momspec(16,LEGS)
	common /momenta/ momspec


#ifndef SPEC_M

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

#define k(i) (8*i+1)
#define s(i) (8*i+3)
#define e(i) (8*i+3+Hel(i))
#define ec(i) (8*i+3-Hel(i))
#define Spinor(i,s,d) (s*Hel(i)+8*i+d+5)

#define EpsL 8
#define EpsR 16

#define MomEncoding(f,i) iand(f,JK-1)*JK**(i-1)

#define signbit(i) ibits(i,31,1)
#define IndexDelta(i,j) signbit(ieor(i,j)-1)
#define Digit(i) char(i+48)
#define Polar(r,theta) r*exp(cI*degree*theta)

#define Error(err,msg) call m_(err, __LINE__, __FILE__, msg)
#define Warning(msg) call m_(0, 0, __FILE__, msg)
#define INFO print *,
#define DEB(a,x) print *, a, x
#define LOOP(var,from,to,step) for( var = from; var <= to; var += step ) {
#define ENDLOOP(var) }
#define TEST(i,b) if( *(i) & (1 << (b)) ) {
#define ENDTEST(i,b) }

#define Power(x, n) \
  ((n) == 2) ? (x)*(x) : \
  ((n) == 3) ? (x)*(x)*(x) : \
  ((n) == 4) ? (x)*(x)*(x)*(x) : \
  ((n) == 5) ? (x)*(x)*(x)*(x)*(x) : pow(x, n)

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

#endif

