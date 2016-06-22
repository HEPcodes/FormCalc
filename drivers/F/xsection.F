* xsection.F
* routines to compute the cross-section
* this file is part of FormCalc
* last modified 8 Jun 16 th

#include "decl.h"
#include "inline.h"
#include "contains.h"

#include "process.h"

* particle-type functions
* arg 1 't' = particle type (running number, for comparison etc)
* arg 2 'h' = helicity bit pattern (e.g. 10 = B'01010' for photon)

#define ARG_Ptyp(i,x,o) x(ARG_Xtyp,o)
#define ARG_Xtyp(t,h,o) t
#define ARG_Phel(i,x,o) x(ARG_Xhel,o)
#define ARG_Xhel(t,h,o) h
#define ARG_Pavg(i,x,o) x(ARG_Xavg,o)
#define ARG_Xavg(i,x,o) o

* UFERMION is an unpolarized Dirac fermion, computed with Hel[.] = 0
#define SCALAR(f,o) f(1,4,1)
#define FERMION(f,o) f(2,10,1)
#define UFERMION(f,o) f(2,4,o)
#define PHOTON(f,o) f(3,10,1)
#define VECTOR(f,o) f(4,14,1)
#define GRAVITINO(f,o) f(5,27,1)
#define GRAVITON(f,o) f(6,27,1)
#define TENSOR(f,o) f(7,31,1)

#define SQRTS_EPS 1D-9

#include "partonic.h"
#include KIN
#include MODEL

#if PHOTONRADIATION == SOFT || GLUONRADIATION == SOFT
#include "softradiation.F"
#endif

#if PHOTONRADIATION == HARD || GLUONRADIATION == HARD
#if LEGS != 3
#error "Hard radiation not available for LEGS .ne. 3"
#endif
#include "hardradiation.F"
#endif

#include "cuts.F"


************************************************************************
** ProcessIni translates the polarization string into bit-encoded
** helicities and initializes the model defaults.

	subroutine ProcessIni(fail, pol,
     &    sqrtSfrom, sqrtSto, sqrtSstep,
     &    argc, argv)
#ifdef SAMURAI
	use msamurai
#endif
	implicit none
	integer fail, argc
	character*(*) pol
	character*128 argv(*)
	RealType sqrtSfrom, sqrtSto, sqrtSstep

#include "xsection.h"

#if U77EXT
#define len_trim lnblnk
	integer lnblnk
	external lnblnk
#endif

	integer i, c, bits

	if( len_trim(pol) .ne. LEGS ) then
	  fail = 1
	  return
	endif

	helmask = 0
	do i = 1, LEGS
	  c = ior(ichar(pol(i:i)), 32)
	  if( c .eq. ichar('u') ) then
	    bits = 16+8+4+2+1
	  else if( c .eq. ichar('t') ) then
	    bits = 16+8+2+1
	  else if( c .eq. ichar('p') ) then
	    bits = 16
	  else if( c .eq. ichar('+') ) then
	    bits = 8
	  else if( c .eq. ichar('l') ) then
	    bits = 4
	  else if( c .eq. ichar('-') ) then
	    bits = 2
	  else if( c .eq. ichar('m') ) then
	    bits = 1
	  else
	    Error(fail, "Invalid polarization for leg "/ /Digit(i))
	    return
	  endif
	  pol(i:i) = char(c)
	  helmask = helmask*32 + bits
	enddo

	Lower(SQRTS) = sqrtSfrom
	Upper(SQRTS) = sqrtSto
	Step(SQRTS) = sqrtSstep

	Var(FIXED) = 0
	Step(FIXED) = 1

	Var(TRIVIAL) = 0
	Step(TRIVIAL) = 0

	sqrtS = -1
	fscale = -1
	sqrtSinvalid = 1

	call ltini

#ifdef CUTTOOLS
* args are:
* 1. limitvalue = precision limit above which the mp routines are activated
* 2. scaloop = scalar integrals (1: LoopTools)
	rootsvalue = sqrtSfrom
	call ctsinit(1D-2, 1)
#endif

#ifdef SAMURAI
* args are:
* 1. imeth = "diag" (numerators) or "tree" (products of tree amps)
* 2. isca = scalar integrals (1: QCDloop, 2: OneLOop, 4: LoopTools)
* 3. verbosity = 0, 1, 2, 3
* 4. itest = 0 (none), 1 (powertest), 2 (nntest), 3 (lnntest)
	call initsamurai("diag", 4, 3, 1)
#endif

#ifdef NINJA
c	call NinjaFormCalcSetVerbosity(3)
c	call NinjaFormCalcDebug(1)
#endif

	call ModelDefaults(argc, argv)
	call LumiDefaults(argc, argv)

	fail = 0
	end


************************************************************************
** ProcessExi wraps up the calculation, e.g. prints a summary of
** messages, deallocates arrays etc.

	subroutine ProcessExi
#ifdef SAMURAI
	use msamurai
#endif
	implicit none

#ifdef SAMURAI
	call exitsamurai
#endif

	call ltexi
	end


************************************************************************
** ParameterScan performs the actual calculation.
** It scans over the parameter loops LOOP* declared in run.F.

	subroutine ParameterScan(dir,
     &    serialfrom, serialto, serialstep)
	implicit none
	character*(*) dir
	integer serialfrom, serialto, serialstep

#include "xsection.h"

	integer openlog
	external openlog

	integer next

	next = serialfrom
	serial = 0

#ifdef LOOP1
	LOOP1
#endif
#ifdef LOOP2
	LOOP2
#endif
#ifdef LOOP3
	LOOP3
#endif
#ifdef LOOP4
	LOOP4
#endif
#ifdef LOOP5
	LOOP5
#endif
#ifdef LOOP6
	LOOP6
#endif
#ifdef LOOP7
	LOOP7
#endif
#ifdef LOOP8
	LOOP8
#endif
#ifdef LOOP9
	LOOP9
#endif
#ifdef LOOP10
	LOOP10
#endif
#ifdef LOOP11
	LOOP11
#endif
#ifdef LOOP12
	LOOP12
#endif
#ifdef LOOP13
	LOOP13
#endif
#ifdef LOOP14
	LOOP14
#endif
#ifdef LOOP15
	LOOP15
#endif
#ifdef LOOP16
	LOOP16
#endif
#ifdef LOOP17
	LOOP17
#endif
#ifdef LOOP18
	LOOP18
#endif
#ifdef LOOP19
	LOOP19
#endif
#ifdef LOOP20
	LOOP20
#endif

	serial = serial + 1
	if( serial .lt. next ) goto 1

	call flush(6)
	if( openlog(dir, serial) .eq. 0 ) then
	  call TotalXS
	  call flush(6)
	  call closelog
	endif

	next = next + serialstep
	if( next .gt. serialto ) return

1	continue
	end


************************************************************************
** SetEnergy sets the energy for the partonic scattering process.
** All scale-dependent quantities are initialized at FSCALE.

	subroutine SetEnergy(fail, newsqrtS, sqme)
	implicit none
	integer fail
	RealType newsqrtS
	external sqme

#include "xsection.h"

	RealType oldfscale, oldmass_in, oldmass_out
	logical reset
	integer i

	fail = 0

	call sqme(mass(1,pid), 0, 2**BIT_SETMASS)

	reset = abs(newsqrtS - sqrtS) .gt. SQRTS_EPS
	sqrtS = newsqrtS
	oldfscale = max(Re(FSCALE), 1D0)

	if( sqrtSinvalid .eq. 1 .or.
     &      btest(flags, BIT_RESET) .or.
     &      abs(oldfscale - fscale) .gt. SQRTS_EPS ) then
	  oldmass_in = MASS_IN(pid)
	  oldmass_out = MASS_OUT(pid)

	  do i = 1, 10
	    call ModelVarIni(fail, oldfscale)
	    call sqme(mass(1,pid), 0, 2**BIT_SETMASS)

	    mass_in = MASS_IN(pid)
	    mass_out = MASS_OUT(pid)
	    fscale = max(Re(FSCALE), 1D0)

	    if( abs(fscale - oldfscale) +
     &          abs(mass_in - oldmass_in) +
     &          abs(mass_out - oldmass_out) .lt. SQRTS_EPS ) goto 1

	    oldfscale = fscale
	    oldmass_in = mass_in
	    oldmass_out = mass_out
	  enddo

1	  threshold(pid) = max(mass_in, mass_out)
	  minthreshold = min(minthreshold, threshold(pid))
	  sqrtSinvalid = fail
	  reset = .TRUE.

#ifdef MUDIM
	  call setmudim(Re(MUDIM))
#endif
#ifdef DELTA
	  call setdelta(Re(DELTA))
#endif
#ifdef LAMBDA
	  call setlambda(Re(LAMBDA))
#endif
	  Divergence = getdelta()
	  mudim = getmudim()
	  muscale = mudim
	  lambda = getlambda()
	  epsi = getepsi()
	endif

	if( reset ) flags = ibset(flags, BIT_RESET)
	end


************************************************************************
** TotalXS computes the total cross-section at one point in parameter
** space.  The cross-section is differential in all variables with 
** non-zero step, i.e. integration variables have step = 0.

	subroutine TotalXS
	implicit none

#include "xsection.h"

	integer fail, f, v, c, ndim
	integer nshow, nfixed, ishow(nvars), ifixed(nvars)
	RealType xs(NCOMP,2), show(nvars)

#define CPdecl
#include "partonic.h"
#undef CPdecl

	integer h, bc1, bc2, bc3, bc4, bitcount, spin_df
	integer*8 h8
	integer m1, m2, m3, m4
	parameter (m1 = Z'555', m2 = Z'333', m3 = Z'F0F', m4 = Z'3F')
	bc1(h) = h - iand(ishft(h, -1), m1)
	bc2(h) = iand(h, m2) + iand(ishft(h, -2), m2)
	bc3(h) = iand(h + ishft(h, -4), m3)
	bc4(h) = iand(h + ishft(h, -8), m4)
	bitcount(h) = bc4(bc3(bc2(bc1(h))))
	spin_df(h8) = bitcount(int(ibits(h8, (LEGS-1)*ldQH, ldQH)))
#if LEGS_IN == 2
     &               *bitcount(int(ibits(h8, (LEGS-2)*ldQH, ldQH)))
#endif

	INFO "code generated by ", FCVERSION

	call ModelConstIni(fail)
	if( fail .ne. 0 ) goto 999

	flags = 2**BIT_RESET
#if NCOMP > 1
     &    + 2**BIT_LOOP
#endif

	call KinIni(fail)
	if( fail .ne. 0 ) goto 999

	minthreshold = 1D300

#define PartonicIni(sqme) call LumiIni(fail, sqme)
#define CPini
#include "partonic.h"
#undef CPini

	call ModelDigest

#define SHOW print 100,
100	format("|# ", A, "=", F10.4, SP, F10.4, " I")

#ifdef MMA
	call MmaBeginPara(serial)
#endif

#ifdef PRINT1
	PRINT1
#endif
#ifdef PRINT2
	PRINT2
#endif
#ifdef PRINT3
	PRINT3
#endif
#ifdef PRINT4
	PRINT4
#endif
#ifdef PRINT5
	PRINT5
#endif
#ifdef PRINT6
	PRINT6
#endif
#ifdef PRINT7
	PRINT7
#endif
#ifdef PRINT8
	PRINT8
#endif
#ifdef PRINT9
	PRINT9
#endif
#ifdef PRINT10
	PRINT10
#endif
#ifdef PRINT11
	PRINT11
#endif
#ifdef PRINT12
	PRINT12
#endif
#ifdef PRINT13
	PRINT13
#endif
#ifdef PRINT14
	PRINT14
#endif
#ifdef PRINT15
	PRINT15
#endif
#ifdef PRINT16
	PRINT16
#endif
#ifdef PRINT17
	PRINT17
#endif
#ifdef PRINT18
	PRINT18
#endif
#ifdef PRINT19
	PRINT19
#endif
#ifdef PRINT20
	PRINT20
#endif

#ifdef MMA
	call MmaEndPara
#endif

	nshow = 0
	nfixed = 0
	do v = MINVAR, MAXVAR
	  if( Step(v) .ne. 0 ) then
	    nshow = nshow + 1
	    ishow(nshow) = v
	    if( Step(v) .ne. -999 ) then
	      nfixed = nfixed + 1
	      ifixed(nfixed) = v
	      Var(v) = Lower(v)
	    endif
	  endif
	enddo

	ndim = nvars - nshow

1	do c = 1, NCOMP
	  xs(c,1) = 0
	enddo
	do c = 1, NCOMP
	  xs(c,2) = 0
	enddo

#if JOIN_PARTONIC
#define PartonicXS(sqme) call IntegrateXS(xs, sqme, ndim)
#define CPxs
#include "partonic.h"
#undef CPxs
#else
* note: SQUAREDME_FUNC is a dummy here, used only because it is
* suitably declared so as not to give a compiler warning
	call IntegrateXS(xs, SQUAREDME_FUNC, ndim)
#endif

	do c = 1, NCOMP
	  xs(c,2) = sqrt(xs(c,2))
	enddo

	do f = 1, nshow
	  show(f) = Show(ishow(f))
	enddo

#ifdef MMA
	call MmaData(show, nshow, xs(1,1), xs(1,2), NCOMP)
#else
* Note: "real" data lines are tagged with "|" in the output.
101	format("| ", 10(4G19.10, :, /"|+"))
	print 101, (show(f), f = 1, nshow)

102	format("|+   ", NCOMP G24.15)
	print 102, xs

	call flush(6)
#endif

	do f = nfixed, 1, -1
	  v = ifixed(f)
	  Var(v) = Var(v) + Step(v)
	  if( (Var(v) - Upper(v))/Step(v) .lt. 1D-10 ) goto 1
	  Var(v) = Lower(v)
	enddo

#ifndef MMA
103	format("|"/"|"/)
	print 103
#endif

999	continue
	end


************************************************************************
** IntegrateXS is a chooser for the Cuba routines, with special cases
** for ndim = 0 (integrand evaluation) and ndim = 1 (non-Cuba
** integration).

	subroutine IntegrateXS(xs, sqme, ndim)
	implicit none
	RealType xs(NCOMP,2)
	external sqme
	integer ndim

#include "xsection.h"

	RealType res(NCOMP), err(NCOMP), prob(NCOMP)
	integer nregions, neval, fail, c

	integer novec
	parameter (novec = 1)

	integer DifferentialXS
	external DifferentialXS

#define FLAGS_LAST 4
#define FLAGS_PSEUDO 8

	if( ndim .eq. 0 ) then

	  c = DifferentialXS(0, res, NCOMP, res, sqme)
	  do c = 1, NCOMP
	    err(c) = 0
	  enddo
	  goto 1

	else if( ndim .eq. 1 ) then

#if METHOD == GAUSS
	  neval = 32
	  call Gauss(NCOMP, 0D0, 1D0, DifferentialXS, sqme,
     &      neval, res)
	  do c = 1, NCOMP
	    err(c) = 0
	  enddo
	  fail = 0
	  INFO "Gauss integration results:"
#else
	  call Patterson(NCOMP, 0D0, 1D0, DifferentialXS, sqme,
     &      Re(RELACCURACY), Re(ABSACCURACY),
     &      neval, fail, res, err)
	  nregions = 1
	  INFO "Patterson integration results:"
#endif
	  nregions = 1
	  do c = 1, NCOMP
	    prob(c) = -1
	  enddo

	else

#if METHOD == VEGAS
	  call vegas(ndim, NCOMP, DifferentialXS, sqme, novec,
     &      Re(RELACCURACY), Re(ABSACCURACY),
     &      VERBOSE, SEED, MINEVAL, MAXEVAL,
     &      NSTART, NINCREASE, NBATCH,
     &      GRIDNO, STATEFILE, SPIN,
     &      neval, fail, res, err, prob)
	  nregions = 1
	  INFO "Vegas integration results:"
#elif METHOD == SUAVE
	  call suave(ndim, NCOMP, DifferentialXS, sqme, novec,
     &      Re(RELACCURACY), Re(ABSACCURACY),
     &      VERBOSE + FLAGS_LAST, SEED, MINEVAL, MAXEVAL,
     &      NNEW, NMIN, Re(FLATNESS), STATEFILE, SPIN,
     &      nregions, neval, fail, res, err, prob)
	  INFO "Suave integration results:"
#elif METHOD == DIVONNE
	  call divonne(ndim, NCOMP, DifferentialXS, sqme, novec,
     &      Re(RELACCURACY), Re(ABSACCURACY),
     &      VERBOSE, SEED, MINEVAL, MAXEVAL,
     &      KEY1, KEY2, KEY3, MAXPASS,
     &      Re(BORDER), Re(MAXCHISQ), Re(MINDEVIATION),
     &      0, NDIM, 0, 0, 0, STATEFILE, SPIN,
     &      nregions, neval, fail, res, err, prob)
	  INFO "Divonne integration results:"
#else
	  call cuhre(ndim, NCOMP, DifferentialXS, sqme, novec,
     &      Re(RELACCURACY), Re(ABSACCURACY),
     &      VERBOSE, MINEVAL, MAXEVAL,
     &      KEY, STATEFILE, SPIN,
     &      nregions, neval, fail, res, err, prob)
	  INFO "Cuhre integration results:"
#endif

	endif

	INFO "nregions =", nregions
	INFO "neval    =", neval
	INFO "fail     =", fail

	if( fail .gt. 0 ) then
	  Warning("Failed to reach the desired accuracy")
	else if( fail .lt. 0 ) then
	  Error(fail, "Integration error")
	endif

100	format(I2, G24.15, " +- ", G24.15, "    p = ", F6.3)
	print 100, (c, res(c), err(c), prob(c), c = 1, NCOMP)

1	do c = 1, NCOMP
	  xs(c,1) = xs(c,1) + res(c)
	enddo
	do c = 1, NCOMP
	  xs(c,2) = xs(c,2) + err(c)**2
	enddo
	end


************************************************************************
** DifferentialXS computes the differential cross-section at x.
** For all integration variables (those with zero step) it factors in
** the Jacobian, too.

	integer function DifferentialXS(ndim, x, ncomp, xs, sqme)
	implicit none
	integer ndim, ncomp
	RealType x(ndim), xs(ncomp)
	external sqme

#include "xsection.h"

	RealType jac, span
	integer v, d, c

#define CPdecl
#include "partonic.h"
#undef CPdecl

	DifferentialXS = 0

	jac = 1
	d = 0
	do v = MINVAR, MAXVAR
	  if( Step(v) .eq. 0 ) then
	    span = Upper(v) - Lower(v)
	    d = d + 1
	    Var(v) = Lower(v) + x(d)*span
	    jac = jac*span
	  endif
	  Show(v) = Var(v)
	enddo

	do c = 1, NCOMP
	  xs(c) = 0
	enddo

#if JOIN_PARTONIC
	call EvaluateXS(xs, sqme)
#else
#define PartonicXS(sqme) call EvaluateXS(xs, sqme)
#define CPxs
#include "partonic.h"
#undef CPxs
#endif

	do c = 1, NCOMP
	  xs(c) = xs(c)*jac
	enddo
	end


************************************************************************

	subroutine EvaluateXS(xs, sqme)
	implicit none
	RealType xs(*)
	external sqme

#include "xsection.h"

	RealType res(NCOMP), fac, rad, flux
	integer c

#if PHOTONRADIATION == SOFT
	RealType weight1(LEGS,LEGS)
	integer sizeof_weight1
	parameter (sizeof_weight1 = LEGS*LEGS)
	data weight1 /sizeof_weight1*1/
#endif

	fac = avgfac(pid)

	call Luminosity(fac, sqme)
	if( fac .eq. 0 .or.
     &      sqrtS - threshold(pid) .lt. -SQRTS_EPS ) return

	call Kinematics(fac)
	if( fac .eq. 0 ) return

	call Cuts(fac)
	if( fac .eq. 0 ) return

	if( btest(flags, BIT_RESET) ) then
	  call clearcache
#ifndef NO_RENCONST
	  call RenConst
#endif
	endif

	call sqme(res, helicities, flags)

	flags = ibclr(flags, BIT_RESET)

#if PHOTONRADIATION == SOFT
	call SoftRadiationFactor(rad, charge(1,pid), weight1)
	res(2) = res(2) - 2/(2*pi)**2*Alfa*rad*res(1)
#endif

#if GLUONRADIATION == SOFT
	call SoftRadiationFactor(rad, colorcharge(1,pid), weight1)
	res(2) = res(2) - 2/(2*pi)**2*Alfas*rad*res(1)
#endif

#if LEGS_IN == 1
	flux = 2*sqrtS
#else
	flux = 4/hbar_c2*sqrtS*momspec(SPEC_K,1)
#endif
	fac = fac/((2*pi)**(3*LEGS_OUT - 4)*2*sqrtS*flux)

	do c = 1, NCOMP
	  if( .not. abs(res(c)) .lt. 1D16 ) then
	    Warning("Got strange values from SquaredME:")
	    INFO res
	    INFO "(Did you compute the colour matrix elements?)"
	    stop
	  endif
	  xs(c) = xs(c) + fac*res(c)
	enddo
	end


************************************************************************
** fcsamurai is the FormCalc wrapper for Samurai
** it is placed here so that util.a can be built without requiring
** the .mod files of Samurai

#ifdef SAMURAI
	subroutine fcsamurai(num, numeps, nleg, rank, kj, mj, res)
	use msamurai
	use options
	implicit none
	ComplexType num, numeps, res
	external num, numeps
	integer nleg, rank
	RealType kj(nleg,4), mj(nleg)

#include "decl.h"

	ComplexType tot(-2:0), totr
	integer isca_save
	logical ok

	call samurai(num, tot, totr, kj, mj, nleg, rank, 1, mudim, ok)
	res = tot(-epsi)
	if( epsi .eq. 0 ) then
	  isca_save = isca
	  isca = 5
	  call samurai(numeps, tot, totr, kj, mj, nleg, rank,
     &      1, mudim, ok)
	  res = res - 2*tot(-1)
	  isca = isca_save
	endif
	end
#endif

