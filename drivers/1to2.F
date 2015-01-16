* 1to2.F
* this file supplies num.F with the main program and other
* functions that depend on the kinematics for a 1 -> 2 process
* this file is part of FormCalc
* last modified 15 Jan 03 th


* some definitions for num.F
#define LEGS 3
#define LEGS_IN 1
#define TYPES TYPE1, TYPE2, TYPE3
#define COMMON "1to2.h"


	program main
	implicit none

#include COMMON

	double precision sqrtSfrom, sqrtSto, theta, tree, loop

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
	  call IntegratedCS(tree, loop)
	  write(11, 27) dble(MASS1), tree, tree + loop
	  call flush(11)
	else
	  do theta = THETACUT, pi - (THETACUT), THETASTEP
	    print *, "theta =", theta
	    call flush(6)
	    call DifferentialCS(tree, loop, cos(theta))
	    write(11, 27) theta, tree, tree + loop
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

	double precision sqrtS

	sqrtS = MASS1
	if( sqrtS .le. 0 ) sqrtS = 1
1	continue
	call ModelVarIni(sqrtS, *999)
	if( abs(MASS1 - sqrtS) .gt. 1D-14 ) goto 1
	if( sqrtS .lt. MASS2 + MASS3 ) stop "Below threshold"

	mass2(1) = (MASS1)**2
	mass2(2) = (MASS2)**2
	mass2(3) = (MASS3)**2
	call CalcRenConst

	Pout = ThreeMom(sqrtS, mass2(2), mass2(3))

	flux = 1/(2*(2*pi)**2*4)/mass2(1)*Pout*avgfac

	call VecSet(k(1), dble(MASS1), 0D0, 0D0, 0D0, 1D0)

	reset = .TRUE.
	return

999	return 1
	end


************************************************************************
** DifferentialCS actually invokes the generated function SquaredME.
** If necessary, soft-photon bremsstrahlung corrections are added.

	subroutine DifferentialCS(tree, loop, costh)
	implicit none
	double precision tree, loop, costh

#include COMMON
#include "num.h"
#include "renconst.h"

#ifdef BREMSSTRAHLUNG
	double precision SoftPhotonFactor
	external SoftPhotonFactor
#endif

	double precision sinth, conv

	sinth = sqrt(1 - costh**2)

	call VecSet(k(2), dble(MASS2), Pout, sinth, 0D0, costh)
	call VecSet(k(3), dble(MASS3), Pout, -sinth, 0D0, -costh)

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


************************************************************************
** IntegratedCS integrates DifferentialCS over the scattering angle.

	subroutine IntegratedCS(tree, loop)
	implicit none
	double precision tree, loop

#include COMMON

	integer i
	double precision mid, x, w, t, l
	double precision treeF, treeB, loopF, loopB

#include "gausspoints.F"

* Integration is done in two intervals separately:
* [THETACUT, pi/2] and [pi/2, pi - THETACUT].
* In this way, we get the F-B asymmetry for free, even
* though the number of sampling points is the same,
* e.g. instead of 32 over the whole interval we can do
* with 16 over each half.

	mid = .5D0*cos(dble(THETACUT))

	treeF = 0
	loopF = 0
	treeB = 0
	loopB = 0

	do i = 1, GAUSSPOINTS/2
	  x = mid*gauss_x(i)
	  w = 2*pi*(IDENTICALFACTOR)*mid*gauss_w(i)
*	      ^^^^ from phi integration
	  call DifferentialCS(t, l, mid - x)
	  treeF = treeF + w*t
	  loopF = loopF + w*l
	  call DifferentialCS(t, l, mid + x)
	  treeF = treeF + w*t
	  loopF = loopF + w*l
	  call DifferentialCS(t, l, -mid - x)
	  treeB = treeB + w*t
	  loopB = loopB + w*l
	  call DifferentialCS(t, l, -mid + x)
	  treeB = treeB + w*t
	  loopB = loopB + w*l
	enddo

	tree = treeF + treeB
	loop = loopF + loopB

#ifdef FB_ASYMMETRY
	if( tree .eq. 0 ) then
	  tree = 999
	else
	  tree = (treeF - treeB)/tree
	endif

	if( loop .eq. 0 ) then
	  loop = 999
	else
	  loop = (loopF - loopB)/loop
	endif
#endif
	end

