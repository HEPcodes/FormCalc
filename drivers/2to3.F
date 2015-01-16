* 2to3.F
* this file supplies num.F with the main program and other
* functions that depend on the kinematics for a 2 -> 3 process
* this file is part of FormCalc
* last modified 18 Oct 01 th


* some definitions for num.F
#define LEGS 5
#define TYPES TYPE1, TYPE2, TYPE3, TYPE4, TYPE5
#define COMMON "2to3.h"

* this defines the order of integrations:
#define K50 4
#define COSTH 3
#define K30 2
#define ETA 1

#define NDIM 4
#define NCOMP 2


	program main
	implicit none

#include COMMON

	double precision Ecmsfrom, Ecmsto, Ecms, k50, tree, loop

	call ParseCommandline(Ecmsfrom, Ecmsto)

	call ffini
#ifdef MUDIM
	call setmudim(MUDIM)
#endif
#ifdef DELTA
	call setdelta(DELTA)
#endif
#ifdef LAMBDA
	call setlambda(LAMBDA)
#endif

	upper(COSTH) = cos(THETACUT)
	lower(COSTH) = -upper(COSTH)

	upper(ETA) = 2*pi
	lower(ETA) = 0

27	format(3G24.15)
28	format()

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

	Ecms = Ecmsfrom
	call SetEnergy(Ecms, .TRUE., *1)

	if(Ecmsto .ge. Ecmsfrom) then
18	  continue
	  print *, "sqrt(s) =", Ecms
	  call IntegratedCS(tree, loop, 4)
	  write(11, 27) Ecms, tree, loop
	  call flush(11)
	  if(Ecms .lt. Ecmsto) then
	    Ecms = Ecms + ECMSSTEP
	    call SetEnergy(Ecms, .FALSE., *1)
	    goto 18
	  endif
	else
	  do k50 = lower(K50), upper(K50), K50STEP
	    print *, "sqrt(s) =", Ecms, "  k50 =", k50
	    var(K50) = k50
	    call SetLimits
	    call IntegratedCS(tree, loop, 3)
	    write(11, 27) k50, tree, loop
	    call flush(11)
	  enddo
	endif

	write(11, 28)
1	continue

	close(11)
	call ffexi
	end


************************************************************************
** SetEnergy sets the CMS energy and initializes everything that depends
** on it.

	subroutine SetEnergy(Ecms, fullini, *)
	implicit none
	double precision Ecms
	logical fullini

#include COMMON

	double precision threshold, m, Pin
	logical adjust, ini

25	format("# WARNING: adjusted sqrt(s) to above-threshold value ",
     +    F10.2)
26	format("# ", A, F8.2)

	reset = fullini
	ini = fullini
	adjust = .FALSE.
1	continue
	call model_ini(Ecms, ini, *999)
	m = (MASS3 + MASS4)**2 - MASS5**2
	lower(K50) = max(dble(MASS5) + .01D0, dble(K50CUT))
	threshold = max(
     +    dble(MASS1 + MASS2),
     +    dble(MASS3 + MASS4 + MASS5),
     +    lower(K50) + sqrt(lower(K50)**2 + m) )
	if(Ecms .lt. threshold) then
	  Ecms = threshold + .01D0
	  ini = .FALSE.
	  adjust = .TRUE.
	  goto 1
	endif

	if(adjust) then
	  write(11, 25) Ecms
	  print 25, Ecms
	endif

	if(fullini) then
#define PRINT write(11, 26)
	  PRINT1
	  PRINT2
	  PRINT3
	  PRINT4
	  PRINT5
#undef PRINT

#define PRINT print 26,
	  PRINT1
	  PRINT2
	  PRINT3
	  PRINT4
	  PRINT5
#undef PRINT
	endif

	mass(1) = MASS1**2
	mass(2) = MASS2**2
	mass(3) = MASS3**2
	mass(4) = MASS4**2
	mass(5) = MASS5**2
	call calc_renconst

	upper(K50) = .5D0*(Ecms - m/Ecms)

	Pin = sqrt((.5D0*(Ecms - (mass(1) - mass(2))/Ecms))**2
     +    - mass(2))

	sqrtS = Ecms
	flux = hbar_c2/(8*(2*pi)**4*4)/(Pin*Ecms)*avgfac

	call VecSet(1, 1, dble(MASS1), Pin, 0D0, 0D0, 1D0)
	call VecSet(5, 1, dble(MASS2), Pin, 0D0, 0D0, -1D0)
	return

999	continue
	return 1
	end


************************************************************************
** SetLimits is called after a new value for var(K50) has been set and
** computes the limits of the k30 integration that depend on k50.

	subroutine SetLimits
	implicit none

#include COMMON

	double precision k50, k5, sig, tau, m1, m2

	k50 = var(K50)
	k5 = sqrt(k50**2 - mass(5))
	sig = sqrtS - k50
	tau = sqrtS*(sig - k50) + mass(5)
	m1 = MASS3 + MASS4
	m2 = MASS3 - MASS4
	sig = sig*(1 + m1*m2/tau)
	tau = k5*sqrt((1 - m1**2/tau)*(1 - m2**2/tau))
	lower(K30) = .5D0*(sig - tau)
	upper(K30) = .5D0*(sig + tau)
	end


************************************************************************
** DifferentialCS actually invokes the generated function squared_me.
** If necessary, soft-photon bremsstrahlung corrections are added.
** DifferentialCS = dsigma/(dk50 dk30 dcosth deta)

	subroutine DifferentialCS(tree, loop, k50, costh, k30, eta)
	implicit none
	double precision tree, loop, k50, costh, k30, eta

#include COMMON
#include "renconst.h"

#ifdef BREMSSTRAHLUNG
	double precision SoftPhotonFactor
	external SoftPhotonFactor
#endif

	double precision sinth, cosxi, sinxi, coseta, sineta
	double precision e3x, e3y, e3z
	double precision k3, k4, k5

	k3 = sqrt(k30**2 - mass(3))
	k4 = sqrt((sqrtS - k30 - k50)**2 - mass(4))
	k5 = sqrt(k50**2 - mass(5))
	cosxi = (k4**2 - k3**2 - k5**2)/(2*k3*k5)
	sinxi = sqrt(1 - cosxi**2)
	sinth = sqrt(1 - costh**2)
	coseta = cos(eta)
	sineta = sin(eta)

	e3x = costh*coseta*sinxi + sinth*cosxi
	e3y = sineta*sinxi
	e3z = costh*cosxi - sinth*coseta*sinxi

	call VecSet(9, -1, dble(MASS3), k3, e3x, e3y, e3z)
	call VecSet(13, -1, dble(MASS4), k4,
     +    (-k3*e3x - k5*sinth)/k4, -k3*e3y/k4, (-k3*e3z - k5*costh)/k4)
	call VecSet(17, -1, dble(MASS5), k5,
     +    sinth, 0D0, costh)

	call squared_me(tree, loop, sqrtS,
     +    bpol(1), epol(1),
     +    bpol(2), epol(2),
     +    bpol(3), epol(3),
     +    bpol(4), epol(4),
     +    bpol(5), epol(5), reset)

	tree = flux*tree
	loop = flux*loop + tree

#ifdef BREMSSTRAHLUNG
	loop = loop + tree*SoftPhotonFactor()
#endif

#ifdef WF_RENORMALIZATION
	loop = loop + tree*WF_RENORMALIZATION
#endif
	end


#if defined(VEGAS) || defined(DCUHRE)

#ifdef VEGAS
#include "vegas.F"
#endif

************************************************************************
** IntegratedCS integrates DifferentialCS over the scattering angle.

	subroutine IntegratedCS(tree, loop, ndim)
	implicit none
	double precision tree, loop
	integer ndim

#include COMMON

	double precision result(2)
	integer neval
	external sample

#ifdef DCUHRE
* this is for checking the integration with the DCUHRE algorithm, see
* http://www.acm.org/pubs/citations/journals/toms/1991-17-4/p452-berntsen

	double precision zero(4), one(4)
	data zero /0, 0, 0, 0/, one /1, 1, 1, 1/

	integer nwork, fail
	parameter (nwork = 25000)
	double precision work(nwork)
	double precision error(2)

	call dcuhre(ndim, NCOMP, zero, one, 0, NCALL*MAXITER,
     +    sample, 1D-15, ACCURACY, 0, nwork, 0,
     +    result, error, neval, fail, work)
	if(fail .gt. 1) print *, "integration error ", fail
	if(fail .eq. 1) then
	  print *,
     +      "Warning: DCUHRE failed to reach the desired accuracy."
	  print *, "Remaining relative errors:"
	  print *, error(1)/result(1), " (tree)"
	  print *, error(2)/result(2), " (loop)"
	endif
#else
	call vegas(sample, result, ACCURACY, ndim,
     +    NCALL, MAXITER, neval)
#endif
	print *, neval, " evaluations used"

	tree = result(1)
	loop = result(2)
	end


************************************************************************
** sample samples the differential cross-section at x. The arrays upper
** and lower contain the integration limits.

	subroutine sample(ndim, x, ncomp, result)
	implicit none
	integer ndim, ncomp
	double precision x(ndim), result(ncomp)

#include COMMON

	integer dim
	double precision range, jacobian

	jacobian = 1
	do dim = ndim, 1, -1
	  range = upper(dim) - lower(dim)
	  jacobian = jacobian*range
	  var(dim) = lower(dim) + x(dim)*range
	  if(dim .eq. K50) call SetLimits
	enddo

	call DifferentialCS(result(1), result(2),
     +    var(K50), var(COSTH), var(K30), var(ETA))
	result(1) = result(1)*jacobian
	result(2) = result(2)*jacobian
	end

#else

#include "multigauss.F"

************************************************************************
** IntegratedCS integrates DifferentialCS over the scattering angle.

	subroutine IntegratedCS(tree, loop, ndim)
	implicit none
	double precision tree, loop
	integer ndim

#include COMMON

	double precision result(2)
	integer d

	double precision para(-1:NCOMP, 0:MAXINTV, NDIM), accuracy
	integer intv(3, NDIM), dim, neval
	common /intpara/ para, accuracy, intv, dim, neval

	neval = 0
	accuracy = ACCURACY
	dim = ndim
	do d = 1, ndim
	  para(-1, 0, d) = lower(d)
	  para(-1, 1, d) = upper(d)
	  intv(1, d) = 1
	  intv(2, d) = 1
	  intv(3, d) = 0
	enddo

	call integrate(result)

	print *, "intervals:   ", (intv(2, d), d = 1, ndim)
	print *, "fluctuation: ", (intv(3, d), d = 1, ndim)
	print *, neval, " evaluations used"

	tree = result(1)
	loop = result(2)
	end


************************************************************************
** sample takes care of the multi-dimensional sampling of the
** differential cross-section.

	subroutine sample(result, x)
	implicit none
	double precision result(*), x

#include COMMON

	double precision para(-1:NCOMP, 0:MAXINTV, NDIM), accuracy
	integer intv(3, NDIM), dim, neval
	common /intpara/ para, accuracy, intv, dim, neval

	integer nodes(NDIM)
	data nodes(K50) /NODES_K50/
	data nodes(COSTH) /NODES_COSTH/
	data nodes(K30) /NODES_K30/
	data nodes(ETA) /NODES_ETA/

	var(dim) = x
	if(dim .eq. 1) then
	  call DifferentialCS(result(1), result(2),
     +      var(K50), var(COSTH), var(K30), var(ETA))
	  neval = neval + 1
	  if(mod(neval, 10000) .eq. 0) print *, neval, " evaluations"
	  return
	endif

	if(dim .eq. K50) then
	  call SetLimits
	  para(-1, 0, K30) = lower(K30)
	  para(-1, 1, K30) = upper(K30)
	  intv(1, K30) = 1
	endif
	dim = dim - 1
	accuracy = accuracy/GAIN

	entry integrate(result)

	call multigauss(result, accuracy,
     +    para(-1, 0, dim), intv(1, dim), MAXINTV, nodes(dim))
	accuracy = accuracy*GAIN
	dim = dim + 1
	end

#endif

