* 2to2.F
* this file supplies num.F with the main program and other
* functions that depend on the kinematics for a 2 -> 2 process
* this file is part of FormCalc
* last modified 22 Aug 01 th


* some definitions for num.F
#define LEGS 4
#define TYPES TYPE1, TYPE2, TYPE3, TYPE4
#define COMMON "2to2.h"


	program main
	implicit none

#include COMMON

	double precision Ecmsfrom, Ecmsto, Ecms, theta, tree, loop

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
	  call IntegratedCS(tree, loop)
	  write(11, 27) Ecms, tree, loop
	  call flush(11)
	  if(Ecms .lt. Ecmsto) then
	    Ecms = Ecms + ECMSSTEP
	    call SetEnergy(Ecms, .FALSE., *1)
	    goto 18
	  endif
	else
	  do theta = THETACUT, pi - THETACUT, THETASTEP
	    print *, "sqrt(s) =", Ecms, "  theta =", theta
	    call DifferentialCS(tree, loop, cos(theta))
	    write(11, 27) theta, tree, loop
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

	double precision threshold, Pin
	logical adjust, ini

25	format("# WARNING: adjusted sqrt(s) to above-threshold value ",
     +    F10.2)
26	format("# ", A, F8.2)

	reset = fullini
	ini = fullini
	adjust = .FALSE.
1	continue
	call model_ini(Ecms, ini, *999)
	threshold = max(
     +    dble(MASS1 + MASS2),
     +    dble(MASS3 + MASS4) )
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
	call calc_renconst

	Pin = sqrt((.5D0*(Ecms - (mass(1) - mass(2))/Ecms))**2
     +    - mass(2))
	Pout = sqrt((.5D0*(Ecms - (mass(3) - mass(4))/Ecms))**2
     +    - mass(4))

	sqrtS = Ecms
	flux = hbar_c2/(4*(2*pi)**2*4)/(Pin*Ecms**2)*Pout*avgfac

	call VecSet(1, 1, dble(MASS1), Pin, 0D0, 0D0, 1D0)
	call VecSet(5, 1, dble(MASS2), Pin, 0D0, 0D0, -1D0)
	return

999	continue
	return 1
	end


************************************************************************
** DifferentialCS actually invokes the generated function squared_me.
** If necessary, soft-photon bremsstrahlung corrections are added.

	subroutine DifferentialCS(tree, loop, costh)
	implicit none
	double precision tree, loop, costh

#include COMMON
#include "renconst.h"

#ifdef BREMSSTRAHLUNG
	double precision SoftPhotonFactor
	external SoftPhotonFactor
#endif

	double precision sinth, conv

	sinth = sqrt(1 - costh**2)

	call VecSet(9, -1, dble(MASS3), Pout, sinth, 0D0, costh)
	call VecSet(13, -1, dble(MASS4), Pout, -sinth, 0D0, -costh)

	call squared_me(tree, loop, sqrtS,
     +    bpol(1), epol(1),
     +    bpol(2), epol(2),
     +    bpol(3), epol(3), 
     +    bpol(4), epol(4), reset)

	tree = flux*tree
	loop = flux*loop + tree

#ifdef BREMSSTRAHLUNG
	loop = loop + tree*SoftPhotonFactor()
#endif

#ifdef WF_RENORMALIZATION
	loop = loop + tree*WF_RENORMALIZATION
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

#include "gauss.F"

* Integration is done in two intervals separately:
* [THETACUT, pi/2] and [pi/2, pi - THETACUT].
* In this way, we get the F-B asymmetry for free, even
* though the number of sampling points is the same,
* e.g. instead of 32 over the whole interval we can do
* with 16 over each half.

	mid = .5D0*cos(THETACUT)

	treeF = 0
	loopF = 0
	treeB = 0
	loopB = 0

	do i = 1, GAUSSPOINTS/2
	  x = mid*gauss_x(i)
	  w = 2*pi*IDENTICALFACTOR*mid*gauss_w(i)
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
	if(tree .eq. 0) then
	  tree = 999
	else
	  tree = (treeF - treeB)/tree
	endif

	if(loop .eq. 0) then
	  loop = 999
	else
	  loop = (loopF - loopB)/loop
	endif
#endif
	end

