* xsection.F
* routines to compute the cross-section
* this file is part of FormCalc
* last modified 15 Jan 07 th


#include "process.h"

#include MODEL

#ifdef BREMSSTRAHLUNG
#include "softphoton.F"
#endif

	block data xsection_h
	implicit none

#include "xsection.h"

	data sqrtS /-1/
	data threshold /-1/
	data scale /-1/
	data sqrtSinvalid /1/
	data Var(FIXED) /0/
	data Step(FIXED) /1/
	data Var(TRIVIAL) /0/
	data Step(TRIVIAL) /0/
	end


************************************************************************
** ProcessIni translates the polarization string into bit-encoded
** helicities, determines the averaging factor, and initializes the
** model defaults.

	subroutine ProcessIni(fail, pol,
     &    sqrtSfrom, sqrtSto, sqrtSstep)
	implicit none
	integer fail
	character*(*) pol
	double precision sqrtSfrom, sqrtSto, sqrtSstep

#include "xsection.h"

	integer lnblnk
	external lnblnk

	integer i, c, bits, df

#define SCALAR 0
#define FERMION 1
#define PHOTON 2
#define GLUON PHOTON
#define VECTOR 3

	integer type(LEGS)
	data type /TYPES/

	if( lnblnk(pol) .ne. LEGS ) then
	  fail = 1
	  return
	endif

	df = 2**(LEGS)

	do i = 1, LEGS
	  if( type(i) .eq. SCALAR ) then
	    c = ichar('S')
	    bits = 1
	  else
	    c = ichar(pol(i:i))
	    if( c .eq. ichar('+') ) then
	      bits = 4
	    else if( c .eq. ichar('-') ) then
	      bits = 1
	    else
	      c = iand(c, 223)
	      if( c .eq. ichar('U') ) then
	        bits = 7
	      else if( c .eq. ichar('T') ) then
	        bits = 5
	      else if( c .eq. ichar('L') ) then
	        bits = 2
	      else
	        bits = 0
	      endif

	      if( type(i) .ne. VECTOR ) then
	        bits = iand(bits, -3)
	        c = ior(c, 1)
	      endif

	      if( bits .eq. 0 ) then
	        WARN "Invalid polarization for leg ", Digit(i)
	        fail = 2
	        return
	      endif

	      if( i .le. LEGS_IN ) df = df*(ibits(bits, 0, 1) +
     &          ibits(bits, 1, 1) + ibits(bits, 2, 1))

#ifdef DIRACFERMIONS
	      if( type(i) .eq. FERMION ) then
	        df = df/2
	        bits = 2
	      endif
#endif
	    endif
	  endif

	  pol(i:i) = char(c)
	  helicities = helicities*8 + bits
	enddo

	avgfac = 2**(LEGS)*DBLE(COLOURFACTOR)*DBLE(IDENTICALFACTOR)/df

	Lower(SQRTS) = sqrtSfrom
	Upper(SQRTS) = sqrtSto
	Step(SQRTS) = sqrtSstep

	call ffini
	call ModelDefaults
	call LumiDefaults

	fail = 0
	end


************************************************************************
** SetEnergy sets the energy for the partonic scattering process.
** All scale-dependent quantities are initialized at SCALE.

	subroutine SetEnergy(fail, newsqrtS)
	implicit none
	integer fail
	double precision newsqrtS

#include "xsection.h"

	double precision newthreshold, newscale

	if( newsqrtS .eq. sqrtS ) then
	  fail = 2*sqrtSinvalid
	  return
	endif

	fail = 0
	sqrtS = newsqrtS

1	newthreshold = max(DBLE(MASS_IN), DBLE(MASS_OUT))
	newscale = max(DBLE(SCALE), 1D0)

	if( abs(threshold - newthreshold) +
     &      abs(scale - newscale) .gt. 1D-13 ) then
	  threshold = newthreshold
	  scale = newscale
	  call ModelVarIni(fail, scale)
	  goto 1
	endif

	if( fail .ne. 0 .or. sqrtS .lt. threshold ) then
	  fail = 1
	else
	  call InitialState
	  reset = .TRUE.
	endif

	sqrtSinvalid = fail
	end


************************************************************************
** ParameterScan performs the actual calculation.
** It scans over the parameter loops LOOP1..15 declared in run.F.

	subroutine ParameterScan(dir,
     &    serialfrom, serialto, serialstep)
	implicit none
	character*(*) dir
	integer serialfrom, serialto, serialstep

#include "util.h"
#include "xsection.h"

	integer openlog
	external openlog

	integer serial, next

	next = serialfrom
	serial = 0

	LOOP1
	LOOP2
	LOOP3
	LOOP4
	LOOP5
	LOOP6
	LOOP7
	LOOP8
	LOOP9
	LOOP10
	LOOP11
	LOOP12
	LOOP13
	LOOP14
	LOOP15

	serial = serial + 1
	if( serial .lt. next ) goto 1

	if( openlog(dir, serial) .eq. 0 ) then

#define SHOW print 100,
100	  format("|# ", A, "=", F10.4, SP, F10.4, " I")

#ifdef MMA
	  call MmaSetPara
#endif

	  PRINT1
	  PRINT2
	  PRINT3
	  PRINT4
	  PRINT5
	  PRINT6
	  PRINT7
	  PRINT8
	  PRINT9
	  PRINT10
	  PRINT11
	  PRINT12
	  PRINT13
	  PRINT14
	  PRINT15

#ifdef MMA
	  call MmaEndSet
#endif

	  call IntegratedCS
	  call closelog
	endif

	next = next + serialstep
	if( next .gt. serialto ) return

1	continue
	end


************************************************************************
** IntegratedCS computes the integrated cross-section at one point in
** parameter space.  The cross-section is differential in all variables
** with non-zero step, i.e. integration variables have step = 0.

	subroutine IntegratedCS
	implicit none

#include "looptools.h"
#include "xsection.h"

	integer nvars
	parameter (nvars = MAXVAR - (MINVAR) + 1)

	external DifferentialCS

	double precision result(NCOMP), error(NCOMP), show(nvars)
	integer fail, f, v, nfix, ndim, fix(nvars)

#ifdef MMA
	call MmaSetData
#endif

	call ModelConstIni(fail)
	if( fail .ne. 0 ) goto 999

	reset = .TRUE.

#ifdef MUDIM
	call setmudim(DBLE(MUDIM))
#endif
#ifdef DELTA
	call setdelta(DBLE(DELTA))
#endif
#ifdef LAMBDA
	call setlambda(DBLE(LAMBDA))
#endif
	Divergence = getdelta()

	call KinIni(fail)
	if( fail .ne. 0 ) goto 999

	nfix = 0
	do v = MINVAR, MAXVAR
	  if( Step(v) .ne. 0 ) then
	    nfix = nfix + 1
	    fix(nfix) = v
	    Var(v) = Lower(v)
	  endif
	enddo

	ndim = nvars - nfix

	call LumiIni(fail)
	if( fail .ne. 0 ) goto 999

1	call Cuba(ndim, DifferentialCS, result, error)

	do f = 1, nfix
	  show(f) = Show(fix(f))
	enddo

#ifdef MMA

	call MmaData(show, nfix, result, error, NCOMP)

#else

* Note: "real" data lines are tagged with "|" in the output.

100	format("| ", 10(4G19.10, :, /"|+"))
	print 100, (show(f), f = 1, nfix)

101	format("|+   ", NCOMP G24.15)
	print 101, result
	print 101, error

	call flush(6)

#endif

	do f = nfix, 1, -1
	  v = fix(f)
	  Var(v) = Var(v) + Step(v)
	  if( (Var(v) - Upper(v))/Step(v) .lt. 1D-10 ) goto 1
	  Var(v) = Lower(v)
	enddo

999	continue

#ifdef MMA

	call MmaEndSet

#else

102	format("|"/"|"/)
	print 102

#endif

	end


************************************************************************
** DifferentialCS computes the differential cross-section at x.
** For all integration variables (those with zero step) it factors in
** the Jacobian, too.

	subroutine DifferentialCS(ndim, x, ncomp, result)
	implicit none
	integer ndim, ncomp
	double precision x(ndim), result(ncomp)

#include "util.h"
#include "xsection.h"
#include "renconst.h"

#ifdef BREMSSTRAHLUNG
	double precision SoftPhotonFactor
	external SoftPhotonFactor
#endif

	double precision fac, range, flux
	integer v, d, c

	fac = avgfac
	d = 0
	do v = MINVAR, MAXVAR
	  if( Step(v) .eq. 0 ) then
	    range = Upper(v) - Lower(v)
	    d = d + 1
	    Var(v) = Lower(v) + x(d)*range
	    fac = fac*range
	  endif
	  Show(v) = Var(v)
	enddo

	do c = 1, ncomp
	  result(c) = 0
	enddo

	call Luminosity(fac)
	if( fac .eq. 0 ) return

	call FinalState(fac)
	if( fac .eq. 0 ) return

	if( reset ) then
	  call clearcache
	  call CalcRenConst
	endif

	call SquaredME(result, helicities, reset)

	reset = .FALSE.

#if LEGS_IN == 1
	flux = 2*sqrtS
#else
	flux = 4/hbar_c2*sqrtS*momspec(SPEC_K, 1)
#endif
	fac = fac/((2*pi)**(3*LEGS_OUT - 4)*2*sqrtS*flux)

	do c = 1, ncomp
	  if( .not. abs(result(c)) .lt. 1D9 ) then
	    WARN "Got strange values from SquaredME:"
	    WARN result
	    WARN "(Did you compute the colour matrix elements?)"
	    stop
	  endif
	  result(c) = result(c)*fac
	enddo

#ifdef BREMSSTRAHLUNG
	result(2) = result(2) + SoftPhotonFactor()*result(1)
#endif

#ifdef WF_RENORMALIZATION
	result(2) = result(2) + (WF_RENORMALIZATION)*result(1)
#endif
	end


************************************************************************
** Cuba is a chooser for the Cuba routines, with special cases
** for ndim = 0 (integrand evaluation) and ndim = 1 (Patterson
** integration).

	subroutine Cuba(ndim, integrand, result, error)
	implicit none
	integer ndim
	external integrand
	double precision result(NCOMP), error(NCOMP)

	integer nregions, neval, fail, c
	double precision prob(NCOMP)

#define GAUSS 1
#define PATTERSON 2
#define VEGAS 3
#define SUAVE 4
#define DIVONNE 5
#define CUHRE 6

#define FLAGS_LAST 4
#define FLAGS_PSEUDO 8

	if( ndim .eq. 0 ) then

	  call integrand(0, 0D0, NCOMP, result)
	  do c = 1, NCOMP
	    error(c) = 0
	  enddo
	  return

	else if( ndim .eq. 1 ) then

#if METHOD == GAUSS

	  neval = 32
	  call Gauss(NCOMP, 0D0, 1D0, integrand,
     &      neval, result)
	  do c = 1, NCOMP
	    error(c) = -1
	    prob(c) = -1
	  enddo
	  nregions = 1
	  fail = 0
	  INFO "Gauss integration results:"

#else

	  call Patterson(NCOMP, 0D0, 1D0, integrand,
     &      DBLE(RELACCURACY), DBLE(ABSACCURACY),
     &      neval, fail, result, error)
	  do c = 1, NCOMP
	    prob(c) = -1
	  enddo
	  nregions = 1
	  INFO "Patterson integration results:"

#endif

	else

#if METHOD == VEGAS

	  call vegas(ndim, NCOMP, integrand,
     &      DBLE(RELACCURACY), DBLE(ABSACCURACY),
     &      VERBOSE, MINEVAL, MAXEVAL,
     &      NSTART, NINCREASE,
     &      neval, fail, result, error, prob)
	  nregions = 1
	  INFO "Vegas integration results:"

#elif METHOD == SUAVE

	  call suave(ndim, NCOMP, integrand,
     &      DBLE(RELACCURACY), DBLE(ABSACCURACY),
     &      VERBOSE + FLAGS_LAST, MINEVAL, MAXEVAL,
     &      NNEW, DBLE(FLATNESS),
     &      nregions, neval, fail, result, error, prob)
	  INFO "Suave integration results:"

#elif METHOD == DIVONNE

	  call divonne(ndim, NCOMP, integrand,
     &      DBLE(RELACCURACY), DBLE(ABSACCURACY),
     &      VERBOSE, MINEVAL, MAXEVAL,
     &      KEY1, KEY2, KEY3, MAXPASS,
     &      DBLE(BORDER), DBLE(MAXCHISQ), DBLE(MINDEVIATION),
     &      0, NDIM, 0, 0, 0,
     &      nregions, neval, fail, result, error, prob)
	  INFO "Divonne integration results:"

#else

	  call cuhre(ndim, NCOMP, integrand,
     &      DBLE(RELACCURACY), DBLE(ABSACCURACY),
     &      VERBOSE + FLAGS_LAST, MINEVAL, MAXEVAL,
     &      KEY,
     &      nregions, neval, fail, result, error, prob)
	  INFO "Cuhre integration results:"

#endif

	endif

	INFO "nregions =", nregions
	INFO "neval    =", neval
	INFO "fail     =", fail

	if( fail .gt. 0 ) then
	  WARN "Failed to reach the desired accuracy."
	else if( fail .lt. 0 ) then
	  WARN "Integration error."
	endif

100	format(I2, G24.15, " +- ", G24.15, "    p = ", F6.3)
	print 100, (c, result(c), error(c), prob(c), c = 1, NCOMP)
	end

