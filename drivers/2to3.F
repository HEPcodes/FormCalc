* 2to3.F
* this file supplies num.F with the main program and other
* functions that depend on the kinematics for a 2 -> 3 process
* this file is part of FormCalc
* last modified 14 Feb 03 th


* some definitions for num.F
#define LEGS 5
#define LEGS_IN 2
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

	double precision sqrtSfrom, sqrtSto, k50, tree, loop

#ifdef PARALLEL
	integer fork, wait, lnblnk
	external fork, wait, lnblnk

	integer seqnumber, processes, pid, i
	character*200 datafile, logfile
	data seqnumber /0/, processes /0/
#endif

	call ParseCommandline(sqrtSfrom, sqrtSto)

	call ffini
#ifdef MUDIM
	call setmudim(dble(MUDIM))
#endif
#ifdef DELTA
	call setdelta(dble(DELTA))
#endif
#ifdef LAMBDA
	call setlambda(dble(LAMBDA))
#endif

	Divergence = getdelta()

	upper(COSTH) = cos(dble(THETACUT))
	lower(COSTH) = -upper(COSTH)

	upper(ETA) = 2*pi
	lower(ETA) = 0

#ifndef PARALLEL
	open(11, file=outfile, status="UNKNOWN")
#endif

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

	call ModelConstIni(*1)

#ifdef PARALLEL
	if( processes .lt. cpus ) then
	  processes = processes + 1
	else
	  pid = wait(0)
	endif

	seqnumber = seqnumber + 1
	pid = fork()
	if( pid .gt. 0 ) goto 1
	if( pid .lt. 0 ) stop "Fork error"

	i = lnblnk(outfile)

24	format(A, ".log-", I5.5)
	write(logfile, 24) outfile(1:i), seqnumber
	open(6, file=logfile, status="UNKNOWN")

25	format(A, ".data-", I5.5)
	write(datafile, 25) outfile(1:i), seqnumber
	open(11, file=datafile, status="UNKNOWN")
#endif

	sqrtS = sqrtSfrom
	call SetEnergy(*1)

	call ModelDigest

26	format("# ", A, F10.4)

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

27	format(3G24.15)

	if( sqrtSto .ge. sqrtSfrom ) then
18	  continue
	    print *, "sqrt(s) =", sqrtS
	    call flush(6)
	    call IntegratedCS(tree, loop, 4)
	    write(11, 27) sqrtS, tree, tree + loop
	    call flush(11)
	  if( sqrtS .lt. sqrtSto ) then
	    sqrtS = sqrtS + SQRTSSTEP
	    call SetEnergy(*1)
	    goto 18
	  endif
	else
	  do k50 = lower(K50), upper(K50), K50STEP
	    print *, "sqrt(s) =", sqrtS, "  k50 =", k50
	    call flush(6)
	    var(K50) = k50
	    call SetLimits
	    call IntegratedCS(tree, loop, 3)
	    write(11, 27) k50, tree, tree + loop
	    call flush(11)
	  enddo
	endif

28	format(/)
	write(11, 28)

#ifdef PARALLEL
	close(6)
	close(11)
	call ffexi
	return
#endif

1	continue

#ifdef PARALLEL
	do i = 1, processes
	  pid = wait(0)
	enddo
#else
	close(11)
	call ffexi
#endif
	end


************************************************************************
** SetEnergy sets the CMS energy and initializes everything that depends
** on it.

	subroutine SetEnergy(*)
	implicit none

#include COMMON
#include "num.h"

	double precision threshold, m, Pin
	logical adjust

	adjust = .FALSE.
1	continue
	call ModelVarIni(sqrtS, *999)
	m = (MASS3 + MASS4)**2 - (MASS5)**2
	lower(K50) = max(dble(MASS5), dble(K50CUT))
	threshold = max(
     &    dble(MASS1 + MASS2),
     &    dble(MASS3 + MASS4 + MASS5),
     &    lower(K50) + sqrt(lower(K50)**2 + m) )
	if( sqrtS .lt. threshold ) then
	  sqrtS = threshold + .01D0
	  adjust = .TRUE.
	  goto 1
	endif

25	format("# WARNING: adjusted sqrt(s) to above-threshold value ",
     &    F10.2)

	if( adjust ) then
	  write(11, 25) sqrtS
	  print 25, sqrtS
	endif

	mass2(1) = (MASS1)**2
	mass2(2) = (MASS2)**2
	mass2(3) = (MASS3)**2
	mass2(4) = (MASS4)**2
	mass2(5) = (MASS5)**2
	call CalcRenConst

	upper(K50) = .5D0*(sqrtS - m/sqrtS)

	Pin = ThreeMom(sqrtS, mass2(1), mass2(2))

	flux = hbar_c2/(8*(2*pi)**4*4)/(Pin*sqrtS)*avgfac

	call VecSet(k(1), dble(MASS1), Pin, 0D0, 0D0, 1D0)
	call VecSet(k(2), dble(MASS2), Pin, 0D0, 0D0, -1D0)

	reset = .TRUE.
	return

999	return 1
	end


************************************************************************
** SetLimits is called after a new value for var(K50) has been set and
** computes the limits of the k30 integration that depend on k50.

	subroutine SetLimits
	implicit none

#include COMMON

	double precision k50, k5, sig, tau, m1, m2

	k50 = var(K50)
	k5 = sqrt(k50**2 - mass2(5))
	sig = sqrtS - k50
	tau = sqrtS*(sig - k50) + mass2(5)
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

	subroutine DifferentialCS(tree, loop, k50, costh, k30, eta)
	implicit none
	double precision tree, loop, k50, costh, k30, eta

#include COMMON
#include "num.h"
#include "renconst.h"

#ifdef BREMSSTRAHLUNG
	double precision SoftPhotonFactor
	external SoftPhotonFactor
#endif

	double precision sinth, cosxi, sinxi, coseta, sineta
	double precision e3x, e3y, e3z
	double precision k3, k4, k5

	k3 = sqrt(k30**2 - mass2(3))
	k4 = sqrt((sqrtS - k30 - k50)**2 - mass2(4))
	k5 = sqrt(k50**2 - mass2(5))
	cosxi = (k4**2 - k3**2 - k5**2)/(2*k3*k5)
	sinxi = sqrt(1 - cosxi**2)
	sinth = sqrt(1 - costh**2)
	coseta = cos(eta)
	sineta = sin(eta)

	e3x = costh*coseta*sinxi + sinth*cosxi
	e3y = sineta*sinxi
	e3z = costh*cosxi - sinth*coseta*sinxi

	call VecSet(k(3), dble(MASS3), k3, e3x, e3y, e3z)
	call VecSet(k(4), dble(MASS4), k4,
     &    (-k3*e3x - k5*sinth)/k4, -k3*e3y/k4, (-k3*e3z - k5*costh)/k4)
	call VecSet(k(5), dble(MASS5), k5, sinth, 0D0, costh)

	call SquaredME(tree, loop, bpol, epol, reset)

	tree = flux*tree
	loop = flux*loop

#ifdef BREMSSTRAHLUNG
	loop = loop + tree*SoftPhotonFactor()
#endif

#ifdef WF_RENORMALIZATION
	loop = loop + tree*(WF_RENORMALIZATION)
#endif
	end


#ifdef MULTIGAUSS
#include "multigauss.F"

************************************************************************
** IntegratedCS integrates DifferentialCS over the scattering angle.

	subroutine IntegratedCS(tree, loop, ndim)
	implicit none
	double precision tree, loop
	integer ndim

#include COMMON

	double precision result(NCOMP)
	integer d

	double precision para(2, 0:NCOMP, 0:MAXINTV, NDIM), relacc
	integer intv(3, NDIM), dim, neval
	common /intpara/ para, relacc, intv, dim, neval

	neval = 0
	relacc = RELACCURACY
	dim = ndim
	do d = 1, ndim
	  para(1, 0, 0, d) = lower(d)
	  para(1, 0, 1, d) = upper(d)
	  intv(1, d) = 1
	  intv(2, d) = 1
	  intv(3, d) = 0
	enddo

	call Integrate(result)

	print *, "intervals:   ", (intv(2, d), d = 1, ndim)
	print *, "fluctuation: ", (intv(3, d), d = 1, ndim)
	print *, neval, " evaluations used"

	tree = result(1)
	loop = result(2)
	end


************************************************************************
** Sample takes care of the multi-dimensional sampling of the
** differential cross-section.

	options /recursive

	subroutine Sample(result, x)
	implicit none
	double precision result(*), x

#include COMMON

	double precision para(2, 0:NCOMP, 0:MAXINTV, NDIM), relacc
	integer intv(3, NDIM), dim, neval
	common /intpara/ para, relacc, intv, dim, neval

	integer nodes(NDIM)
	data nodes(K50) /NODES_K50/
	data nodes(COSTH) /NODES_COSTH/
	data nodes(K30) /NODES_K30/
	data nodes(ETA) /NODES_ETA/

	var(dim) = x
	if( dim .eq. 1 ) then
	  call DifferentialCS(result(1), result(2),
     &      var(K50), var(COSTH), var(K30), var(ETA))
	  neval = neval + 1
	  if( mod(neval, 10000) .eq. 0 ) print *, neval, " evaluations"
	  return
	endif

	if( dim .eq. K50 ) then
	  call SetLimits
	  para(1, 0, 0, K30) = lower(K30)
	  para(1, 0, 1, K30) = upper(K30)
	  intv(1, K30) = 1
	endif
	dim = dim - 1
#ifdef GAIN
	relacc = relacc/(GAIN)
#endif

	entry Integrate(result)

	call MultiGauss(result, dble(ABSACCURACY), relacc,
     &    para(1, 0, 0, dim), intv(1, dim), MAXINTV, nodes(dim))

#ifdef GAIN
	relacc = relacc*(GAIN)
#endif
	dim = dim + 1
	end

#else

#ifdef VEGAS
#include "vegas.F"

************************************************************************
** IntegratedCS integrates DifferentialCS over the scattering angle.

	subroutine IntegratedCS(tree, loop, ndim)
	implicit none
	double precision tree, loop
	integer ndim

#include COMMON

	double precision result(NCOMP)
	integer neval
	external Sample

	call Vegas(Sample, result,
     &    dble(ABSACCURACY), dble(RELACCURACY), 
     &    ndim, NCALL, MAXITER, neval)

	print *, neval, " evaluations used"

	tree = result(1)
	loop = result(2)
	end

#else
	include "dcuhre.f"

************************************************************************
** IntegratedCS integrates DifferentialCS over the scattering angle.

	subroutine IntegratedCS(tree, loop, ndim)
	implicit none
	double precision tree, loop
	integer ndim

#include COMMON

	double precision result(NCOMP)
	integer neval
	external Sample

* this is for integration with the DCUHRE algorithm, see
* http://www.sci.wsu.edu/math/faculty/genz/papers/cuhre/cuhre.html
* DCUHRE is algorithm 698 from http://www.netlib.org/toms/

	double precision zero(4), one(4)
	data zero /0, 0, 0, 0/, one /1, 1, 1, 1/

* see dcuhre.f for the determination of NUM
* the first argument of max covers NDIM = 2 and 3
#define NUM(d) max(d*(7 - d)*31 - 245, 1 + 2*d*(7 + d*(3 + 2*d))/3 + 2**d)

#define NWORK(d, n) 2*((maxpts - n)/(2*n) + 1)*(d + NCOMP + 1) + 17*NCOMP + 1

	integer maxpts, dim1, dim2, num1, num2, nwork1, nwork2, nwork
	parameter (maxpts = NCALL*MAXITER)
	parameter (dim1 = NDIM - 1)
	parameter (dim2 = NDIM)
	parameter (num1 = NUM(dim1))
	parameter (num2 = NUM(dim2))
	parameter (nwork1 = NWORK(dim1, num1))
	parameter (nwork2 = NWORK(dim2, num2))
	parameter (nwork = max(nwork1, nwork2))

* this nwork is sufficient for dim1 <= ndim <= dim2

	double precision work(nwork)
	double precision error(NCOMP)
	integer fail

	call dcuhre(ndim, NCOMP, zero, one, 0, maxpts, Sample,
     &    dble(ABSACCURACY), dble(RELACCURACY), 0, nwork, 0,
     &    result, error, neval, fail, work)
	if( fail .gt. 1 ) print *, "integration error ", fail
	if( fail .eq. 1 ) then
	  print *, "DCUHRE failed to reach the desired accuracy."
	  print *, "Remaining relative errors:"
	  if( result(1) .gt. ABSACCURACY )
     &      print *, error(1)/result(1), " (tree)"
	  if( result(2) .gt. ABSACCURACY )
     &      print *, error(2)/result(2), " (loop)"
	endif

	print *, neval, " evaluations used"

	tree = result(1)
	loop = result(2)
	end

#endif VEGAS


************************************************************************
** Sample samples the differential cross-section at x.  The arrays
** upper and lower contain the integration limits.

	subroutine Sample(ndim, x, ncomp, result)
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
	  if( dim .eq. K50 ) call SetLimits
	enddo

	call DifferentialCS(result(1), result(2),
     &    var(K50), var(COSTH), var(K30), var(ETA))
	result(1) = result(1)*jacobian
	result(2) = result(2)*jacobian
	end

#endif MULTIGAUSS

