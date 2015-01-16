	double complex resoc
#ifdef DYSON
	double complex reso, resoT, resoU
#else
	double precision reso, resoT, resoU
#endif
	double precision S, T, U, MH2, ME2, MD2, MU2
	double precision psq, p, ps, EE, EEs, st, ct, kinf
	double precision bornamp, feynamp, bornsum, feynsum
	double precision Ecms, th
	integer gen, Cptr2, Dptr2
	common /kin/ 
     +    resoc, reso, resoT, resoU,
     +    S, T, U, MH2, ME2, MD2, MU2,
     +    psq, p, ps, EE, EEs, st, ct, kinf,
     +    bornamp, feynamp, bornsum, feynsum,
     +    Ecms, th,
     +    gen, Cptr2, Dptr2

	double complex v(4,16)
	common /vec/ v

	double precision avgfac
	integer bpol(4), epol(4)
	integer i1, i2, i3, i4
	character polstr*4
	common /pol/
     +    avgfac,
     +    bpol, epol,
     +    i1, i2, i3, i4,
     +    polstr

	double precision pi, ds2, degree
	parameter (pi = 3.1415926535897932384626433832795028841972D0)
	parameter (ds2 = .7071067811865475244008443621048490392848D0)
*	  = 1/sqrt(2)
	parameter (degree = pi/180D0)

#ifdef MBOS
	double precision m2
	parameter (m2 = MBOS**2)
#endif

