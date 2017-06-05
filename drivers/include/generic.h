#if 0
	generic.h
	macros for the use with the functions of specs.h
	this file is part of FormCalc
	last modified 16 Jun 16 th


ARG-type functions:
arg 1 'f' = function
arg 2 'i' = particle no. (1..LEGS)
arg 3 'o' = incoming 1, outgoing 2

particle-type functions (used with Generic):
arg 1 't' = particle type (running number, for comparison etc)
arg 2 'h' = helicity bit pattern (e.g. 10 = B'01010' for photon)
arg 3 'f' = bit mask for fermions (B'11111'), 0 otherwise
arg 4 'b' = bit pattern used for testing helicities

#endif

#define ARG_Ptyp(f,i,o) f(ARG_ptyp,i,o)
#define ARG_ptyp(t,h,f,b) t
#define ARG_Phel(f,i,o) f(ARG_phel,i,o)
#define ARG_phel(t,h,f,b) h
#define ARG_Ferm(f,i,o) f(ARG_ferm,i,o)
#define ARG_ferm(t,h,f,b) f
#define ARG_Bhel(f,i,o) f(ARG_bhel,i,o)
#define ARG_bhel(t,h,f,b) b

#define SCALAR(f,i,o) f(1,4,0,0)
#if DIRACFERMIONS
#define FERMION(f,i,o) f(2,10,31,0)
#else
#define FERMION(f,i,o) f(2,10,0,HelSet(i))
#endif
#define PHOTON(f,i,o) f(3,10,0,HelSet(i))
#define VECTOR(f,i,o) f(4,14,0,HelSet(i))
#define GRAVITINO(f,i,o) f(5,27,0,HelSet(i))
#define GRAVITON(f,i,o) f(6,27,0,HelSet(i))
#define TENSOR(f,i,o) f(7,31,0,HelSet(i))

#define ARG_ID(x,i,o) x
#define ARG_RE(x,i,o) Re(x)

#define JOIN_SEQ(a,b) a,b
#define JOIN_ADD(a,b) a+b
#define JOIN_MUL(a,b) a*b
#define JOIN_OCT(a,b) b+8*(a)
#define JOIN_HEL(a,b) b+QH*(a)

