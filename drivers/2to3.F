* 2to3.F
* the kinematics for a 2 -> 3 process
* this file is part of FormCalc
* last modified 18 Nov 04 th


* some definitions for main.F
#define LEGS 5
#define LEGS_IN 2
#define TYPES TYPE1, TYPE2, TYPE3, TYPE4, TYPE5

* this defines the order of integrations:
#define K50 4
#define COSTH 3
#define K30 2
#define ETA 1

#define MAXDIM 4

#define HEADER "2to3.h"


************************************************************************
** DoIntegrated computes the integrated cross section.

	subroutine DoIntegrated
	implicit none

#include HEADER

	double precision result(NCOMP), error(NCOMP)

	print *, "sqrt(s) =", sqrtS
	call flush(6)

	call IntegratedCS(result, error, MAXDIM)

100	format(TAG, G19.10, NCOMP G24.15, NCOMP G15.6)
	print 100, sqrtS, result, error
	call flush(6)
	end


************************************************************************
** DoDifferential computes the differential cross section.

	subroutine DoDifferential
	implicit none

#include HEADER

	double precision result(NCOMP), error(NCOMP), k50

* at k50 = MASS5 (i.e. k5 = 0) the phase-space volume is precisely
* zero, plus we have difficulties computing certain kinematical
* quantities, such as cosxi ~ 1/k5, hence we move a tiny bit away
* from k5 = 0:
	k50 = max(lower(K50), dble(MASS5) + .1D0)

	do k50 = k50, upper(K50), K50STEP
	  print *, "sqrt(s) =", sqrtS, "  k50 =", k50
	  call flush(6)

	  var(K50) = k50
	  call SetLimits
	  call IntegratedCS(result, error, MAXDIM - 1)

100	  format(TAG, G19.10, NCOMP G24.15, NCOMP G15.6)
	  print 100, k50, result, error
	  call flush(6)
	enddo
	end


************************************************************************
** SetProcess sets up process-dependent constants.

	subroutine SetProcess(hel, avgfac, sqrtSfrom, sqrtSto)
	implicit none
	integer hel
	double precision avgfac, sqrtSfrom, sqrtSto

#include HEADER

	helicities = hel
	preflux = hbar_c2/(8*(2*pi)**4*4)*avgfac

	upper(ETA) = 2*pi
	lower(ETA) = 0
	end


************************************************************************
** SetEnergy sets up energy-dependent constants.

	subroutine SetEnergy(newsqrtS, *)
	implicit none
	double precision newsqrtS

#include HEADER
#include "util.h"

	double precision threshold, m, Pin
	logical adjust

	adjust = .FALSE.
1	continue
	call ModelVarIni(newsqrtS, *999)
	m = (MASS3 + MASS4)**2 - (MASS5)**2
	sqrtS = newsqrtS
	lower(K50) = max(dble(MASS5), dble(K50CUT))
	threshold = max(
     &    dble(MASS1 + MASS2),
     &    dble(MASS3 + MASS4 + MASS5),
     &    lower(K50) + sqrt(lower(K50)**2 + m) )
	if( newsqrtS .le. threshold ) then
	  newsqrtS = threshold + .01D0
	  adjust = .TRUE.
	  goto 1
	endif

100	format(TAG,
     &    "# WARNING: adjusted sqrt(s) to above-threshold value ",
     &    F10.2)
	if( adjust ) print 100, newsqrtS

	upper(K50) = .5D0*(sqrtS - m/sqrtS)

	upper(COSTH) = cos(dble(THETACUT))
	lower(COSTH) = -upper(COSTH)

	Pin = ThreeMom(sqrtS, dble(MASS1), dble(MASS2))

	flux = preflux/(Pin*sqrtS)

	call VecSet(1, dble(MASS1), Pin, 0D0, 0D0, 1D0)
	call VecSet(2, dble(MASS2), Pin, 0D0, 0D0, -1D0)

	call PREFIX CalcRenConst
	reset = .TRUE.
	return

999	return 1
	end


************************************************************************
** SetLimits is called after a new value for var(K50) has been set and
** computes the limits of the k30 integration that depend on k50.

	subroutine SetLimits
	implicit none

#include HEADER

	double precision k50, k5, sig, tau, m1, m2

	k50 = var(K50)
	k5 = sqrt(k50**2 - (MASS5)**2)
	sig = sqrtS - k50
	tau = sqrtS*(sig - k50) + (MASS5)**2
	m1 = MASS3 + MASS4
	m2 = MASS3 - (MASS4)
	sig = sig*(1 + m1*m2/tau)
	tau = k5*sqrt((1 - m1**2/tau)*(1 - m2**2/tau))
	lower(K30) = .5D0*(sig - tau)
	upper(K30) = .5D0*(sig + tau)
	end


************************************************************************
** DifferentialCS actually invokes the generated function SquaredME.
** If necessary, soft-photon bremsstrahlung corrections are added.
** DifferentialCS = dsigma/(dk50 dk30 dcosth deta)

	subroutine DifferentialCS(result, k50, costh, k30, eta)
	implicit none
	double precision result(*), k50, costh, k30, eta

#include HEADER
#include "util.h"
#include "renconst.h"

#ifdef BREMSSTRAHLUNG
	double precision SoftPhotonFactor
	external SoftPhotonFactor
#endif

	double precision sinth, cosxi, sinxi, coseta, sineta
	double precision e3x, e3y, e3z
	double precision k3, k4, k5
	integer c

	k3 = sqrt(k30**2 - (MASS3)**2)
	k4 = sqrt((sqrtS - k30 - k50)**2 - (MASS4)**2)
	k5 = sqrt(k50**2 - (MASS5)**2)
	cosxi = (k4**2 - k3**2 - k5**2)/(2*k3*k5)
	sinxi = sqrt(1 - cosxi**2)
	sinth = sqrt(1 - costh**2)
	coseta = cos(eta)
	sineta = sin(eta)

	e3x = costh*coseta*sinxi + sinth*cosxi
	e3y = sineta*sinxi
	e3z = costh*cosxi - sinth*coseta*sinxi

	call VecSet(3, dble(MASS3), k3, e3x, e3y, e3z)
	call VecSet(4, dble(MASS4), k4,
     &    (-k3*e3x - k5*sinth)/k4, -k3*e3y/k4, (-k3*e3z - k5*costh)/k4)
	call VecSet(5, dble(MASS5), k5, sinth, 0D0, costh)

	call SquaredME(result, helicities, reset)

	do c = 1, NCOMP
	  result(c) = result(c)*flux
	enddo

#ifdef BREMSSTRAHLUNG
	result(2) = result(2) + SoftPhotonFactor()*result(1)
#endif

#ifdef WF_RENORMALIZATION
	result(2) = result(2) + (WF_RENORMALIZATION)*result(1)
#endif
	end


************************************************************************
** Sample samples the differential cross-section at x.  The arrays
** upper and lower contain the integration limits.

	subroutine Sample(ndim, x, ncomp, result)
	implicit none
	integer ndim, ncomp
	double precision x(*), result(*)

#include HEADER

	integer dim, c
	double precision range, jacobian

	jacobian = 1
	do dim = ndim, 1, -1
	  range = upper(dim) - lower(dim)
	  jacobian = jacobian*range
	  var(dim) = lower(dim) + x(dim)*range
	  if( dim .eq. K50 ) call SetLimits
	enddo

	call DifferentialCS(result,
     &    var(K50), var(COSTH), var(K30), var(ETA))

	do c = 1, ncomp
	  result(c) = result(c)*jacobian
	enddo
	end


************************************************************************
** IntegratedCS integrates the differential cross-section over the
** dimensions 1..ndim.

#define VEGAS 1
#define SUAVE 2
#define DIVONNE 3
#define CUHRE 4

#define FLAGS_LAST 4

	subroutine IntegratedCS(result, error, ndim)
	implicit none
	double precision result(*), error(*)
	integer ndim

#include HEADER

	integer nregions, neval, fail, c
	double precision prob(NCOMP)
	external Sample

#if METHOD == VEGAS

	call vegas(ndim, NCOMP, Sample,
     &    dble(RELACCURACY), dble(ABSACCURACY),
     &    VERBOSE, MINEVAL, MAXEVAL,
     &    NSTART, NINCREASE,
     &    neval, fail, result, error, prob)
	nregions = 1

#elif METHOD == SUAVE

	call suave(ndim, NCOMP, Sample,
     &    dble(RELACCURACY), dble(ABSACCURACY),
     &    VERBOSE + FLAGS_LAST, MINEVAL, MAXEVAL,
     &    NNEW, dble(FLATNESS),
     &    nregions, neval, fail, result, error, prob)

#elif METHOD == DIVONNE

	call divonne(ndim, NCOMP, Sample,
     &    dble(RELACCURACY), dble(ABSACCURACY),
     &    VERBOSE, MINEVAL, MAXEVAL,
     &    KEY1, KEY2, KEY3, MAXPASS,
     &    dble(BORDER), dble(MAXCHISQ), dble(MINDEVIATION),
     &    0, NDIM, 0, 0, 0,
     &    nregions, neval, fail, result, error, prob)

#else

	call cuhre(ndim, NCOMP, Sample,
     &    dble(RELACCURACY), dble(ABSACCURACY),
     &    VERBOSE + FLAGS_LAST, MINEVAL, MAXEVAL,
     &    KEY,
     &    nregions, neval, fail, result, error, prob)

#endif

	if( fail .gt. 0 ) then
	  print *, "Failed to reach the desired accuracy."
	  print *, "Results thus obtained:"
100	  format(I2, G24.15, " +- ", G24.15, "  p = ", F8.3)
	  print 100, (c, result(c), error(c), prob(c), c = 1, NCOMP)
	else if( fail .lt. 0 ) then
	  print *, "integration error"
	endif

	print *, neval, " evaluations, ", nregions, " regions"
	end

