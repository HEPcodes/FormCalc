#if 0
	util.h
	prototypes for the util functions
	this file is part of FormCalc
	last modified 9 Mar 13 th
#endif


#ifndef LEGS
#define LEGS 1
#endif

enum { nvec = 10 };

struct {
  ComplexType vec[LEGS*nvec+1][2][2];
} vectors_;

#define vec(i,j,n) vectors_.vec[n][j-1][i-1]

#define k(i) (nvec*(i-1)+1)
#define s(i) (nvec*(i-1)+3)
#define e(i) (nvec*(i-1)+3+Hel(i))
#define ec(i) (nvec*(i-1)+3-Hel(i))
#define Spinor(i,s,d) (s*Hel(i)+nvec*(i-1)+d+5)

#define DEB(a,x) printf(a " (%.13lg,%.13lg)\n", Re(x), Im(x))
#define LOOP(var,from,to,step) for( var = from; var <= to; var += step ) {
#define ENDLOOP(var) }
#define TEST(i,b) if( *(i) & (1 << (b)) ) {
#define ENDTEST(i,b) }

#define BIT_RESET 0
#define BIT_LOOP 1
#define BIT_HEL(i) (5*(LEGS-i)+Hel(i)+2)
#define LOOP_HEL(h) for( h = -2; h <= 2; ++h ) {
#define ENDLOOP_HEL(h) }

#define INI_S(seq) clearcache()
#define INI_ANGLE(seq) markcache()
#define DEINI(seq) restorecache()

#define PREP(h,he, v,ve, a,ae, s,se)
#define EXEC(f, res, flags) f(res, flags)
#define SYNC(res)

