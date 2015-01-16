* num.F
* this is the driver program for a 2 -> 2 process
* for Fortran code produced by WriteSquaredME
* this file is part of FormCalc
* last modified 18 Jun 01 th


#include "prefix.h"
#include "process.h"

#define SCALAR 0
#define FERMION 1
#define PHOTON 2
#define VECTOR 3


************************************************************************
** THE MAIN PROGRAM STARTS HERE:

	program num
	implicit none

#include "kin.h"

	double precision Ecms, Ecmsfrom, Ecmsto, theta
	double precision treeres, loopres

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

	reset = .TRUE.

	if(Ecmsto .ge. Ecmsfrom) then
	  Ecms = Ecmsfrom
18	  continue
	  print *, "sqrt(s) =", Ecms
	  call IntegratedCS(treeres, loopres, Ecms, *1)
	  write(11, 27) Ecms, treeres, loopres
	  call flush(11)
	  if(Ecms .lt. Ecmsto) then
	    Ecms = Ecms + ECMSSTEP
	    goto 18
	  endif
	else
	  do theta = CUTOFF, pi - CUTOFF, THETASTEP
	    print *, "sqrt(s) =", Ecmsfrom, "  theta =", theta
	    call DifferentialCS(treeres, loopres, Ecmsfrom, theta, *1)
	    write(11, 27) theta, treeres, loopres
	    call flush(11)
	  enddo
	endif

	write(11, 28)
1	continue

	close(11)
	call ffexi
	end


************************************************************************
** ParseCommandline reads the polarizations and the energy range from
** the command line. From these values, it constructs a file name and
** opens the file for writing as file #11.

	subroutine ParseCommandline(Ecmsfrom, Ecmsto)
	implicit none
	double precision Ecmsfrom, Ecmsto

#include "kin.h"

	integer*4 iargc
	external iargc

	integer argc, i, j
	character*200 argv, outf
	character*4 polstr
	character c

	integer type(4)
	data type / TYPE1, TYPE2, TYPE3, TYPE4 /

	argc = iargc()
	if(argc .lt. 5 .or. argc .gt. 6) then
	  print *, "USAGE:"
	  print *, "to calculate a differential cross-section:"
	  print *, "  num p1 p2 p3 p4 sqrts"
	  print *, "to calculate an integrated cross-section:"
	  print *, "  num p1 p2 p3 p4 sqrts_from sqrts_to"
	  print *, "for polarizations p1 + p2 -> p3 + p4"
	  print *, "where p[1-4] may be one of"
	  print *, "  0  for longitudinal   polarization"
	  print *, "  +  for right-circular      ''"
	  print *, "  -  for left-circular       ''"
	  print *, "  t  for transverse          ''"
	  print *, "  u  for unpolarized particle"
	  stop
	endif

	avgfac = 1

	do i = 1, 4
	  if(type(i) .eq. SCALAR) then
	    c = 'S'
	    bpol(i) = 0
	    epol(i) = 0
	  else
	    call getarg(i, argv)
	    c = argv(1:1)
	    if(type(i) .eq. FERMION) then
	      if(c .eq. 'u') then
	        j = 0
	        if(i .gt. 2) avgfac = avgfac*2
	        c = 'U'
	      else if(c .eq. '+') then
	        j = 1
	      else if(c .eq. '-') then
	        j = -1
	      else
	        stop "Invalid polarization"
	      endif
	      bpol(i) = j
	      epol(i) = j
	    else
	      if(c .eq. 'u') then
	        if(type(i) .eq. PHOTON)
     +            stop "No longitudinal polarizations allowed"
	        bpol(i) = 0
	        epol(i) = 2
	        if(i .lt. 3) avgfac = avgfac/3D0
	        c = 'U'
	      else if(c .eq. 't') then
	        bpol(i) = 1
	        epol(i) = 2
	        if(i .lt. 3) avgfac = avgfac/2D0
	        c = 'T'
	      else
	        if(c .eq. 'l' .or. c .eq. '0') then
	          if(type(i) .eq. PHOTON)
     +              stop "No longitudinal polarizations allowed"
	          j = 0
	          c = 'L'
	        else if(c .eq. '+') then
	          j = 1
	        else if(c .eq. '-') then
	          j = 2
	        else
	          stop "Invalid polarization"
	        endif
	        bpol(i) = j
	        epol(i) = j
	      endif
	    endif
	  endif
	  polstr(i:i) = c
	enddo

22	format(E20.0)
23	format("diff.pol=", A4, ".E=", I5.5)
24	format("tot.pol=", A4, ".E=", I5.5, "-", I5.5)

	call getarg(5, argv)
	read(argv, 22) Ecmsfrom
	if(argc .eq. 5) then
	  Ecmsto = -1
	  write(outf, 23) polstr, int(Ecmsfrom)
	else
	  call getarg(6, argv)
	  read(argv, 22) Ecmsto
	  write(outf, 24) polstr, int(Ecmsfrom), int(Ecmsto)
	endif

	call getarg(0, argv)
	i = 0
12	continue
	i = i + 1
	if(argv(i:i) .ne. ' ') goto 12
	outf = argv(1:i - 1) // "-" // outf

	open(11, file=outf, status='UNKNOWN')
	print *, "output file: ", outf
	end


************************************************************************
** VecSet constructs the CMS representation of a momentum and the
** three corresponding polarization vectors. The arguments:
** n: the index of the momentum in the array vec
** dir: 1 for incoming, -1 for outgoing particle
** m, p: mass and three-momentum of the particle
** ex, ey, ez: the unit three-vector of the momentum

	subroutine VecSet(n, dir, m, p, ex, ey, ez)
	implicit none
	integer n, dir
	double precision m, p, ex, ey, ez

#include "kin.h"

	double precision E, sinth
	integer rho, n2

	E = sqrt(p**2 + m**2)

	vec(0, n) = E
	vec(1, n) = p*ex
	vec(2, n) = p*ey
	vec(3, n) = p*ez

	if(m .ne. 0) then
	  vec(0, n + 1) = p/m
	  vec(1, n + 1) = E/m*ex
	  vec(2, n + 1) = E/m*ey
	  vec(3, n + 1) = E/m*ez
	endif

	vec(0, n + 2) = 0
	sinth = 1 - ez**2
	if(sinth .ne. 0) then
	  sinth = sqrt(sinth)
	  vec(1, n + 2) = (ex*ez - dir*ey*cI)/sinth/sqrt2
	  vec(2, n + 2) = (ey*ez + dir*ex*cI)/sinth/sqrt2
	  vec(3, n + 2) = -sinth/sqrt2
	else
	  vec(1, n + 2) = 1/sqrt2
	  vec(2, n + 2) = dir*ez/sqrt2*cI
	  vec(3, n + 2) = 0
	endif

	do rho = 0, 3
	  vec(rho, n + 3) = dconjg(vec(rho, n + 2))
	enddo

	do n2 = n, n + 3
	  do rho = 0, 3
	    vec(rho, -n2) = dconjg(vec(rho, n2))
	  enddo
	enddo
	end


************************************************************************
** DifferentialCS actually invokes the generated function squared_me.
** If necessary, soft-photon bremsstrahlung corrections are added.

	subroutine DifferentialCS(treeres, loopres, Ecms, theta, *)
	implicit none
	double precision treeres, loopres, Ecms, theta

#include "kin.h"
#include "renconst.h"

#ifdef BREMSSTRAHLUNG
	double precision SoftPhotonFactor
	external SoftPhotonFactor
#endif

	double precision Ecms2, flux, Pin, Pout
	double precision threshold, sinth, costh
	integer component, n
	logical reini
	save Ecms2, flux, Pin, Pout

	double precision hbar_c2
	parameter (hbar_c2 = 3.8937966D8)
*	  = hbar c^2 in picobarn

	double precision prevEcms
	data prevEcms /-1/

27	format("# WARNING: adjusted sqrt(s) to above-threshold value ",
     +    F10.2)

	if(reset .or. Ecms .ne. prevEcms) then
	  reini = reset
1	  continue
	  call model_ini(Ecms, reini, *999)
	  reini = .FALSE.
	  threshold = max(
     +      dble(MASS1 + MASS2),
     +      dble(MASS3 + MASS4))
	  if(Ecms .lt. threshold) then
	    print *, "Warning: sqrt(s) =", Ecms, " is below threshold."
	    Ecms = threshold + .01D0
	    print *, "Using sqrt(s) = ", Ecms, " now."
	    write(11, 27) Ecms
	    goto 1
	  endif

	  if(reset) then
26	    format("# ", A, F8.2)
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

	  prevEcms = Ecms
	  Ecms2 = Ecms**2

	  Pin = sqrt((.5D0*(Ecms + (mass(2) - mass(1))/Ecms))**2
     +      - mass(2))
	  Pout = sqrt((.5D0*(Ecms + (mass(4) - mass(3))/Ecms))**2
     +      - mass(4))

	  flux = avgfac*COLOURFACTOR*
     +      hbar_c2/(64*pi**2*Ecms2)*Pout/Pin

* some apps need the relative velocity of the final-state particles
	  beta = Pout/(.5D0*Ecms)

	  call VecSet(1, 1, dble(MASS1), Pin, 0D0, 0D0, 1D0)
	  call VecSet(5, 1, dble(MASS2), Pin, 0D0, 0D0, -1D0)
	endif

	sinth = sin(theta)
	costh = cos(theta)
	call VecSet(9, -1, dble(MASS3), Pout, sinth, 0D0, costh)
	call VecSet(13, -1, dble(MASS4), Pout, -sinth, 0D0, -costh)

	call squared_me(treeres, loopres, Ecms2,
     +    bpol(1), epol(1),
     +    bpol(2), epol(2),
     +    bpol(3), epol(3), 
     +    bpol(4), epol(4), reset)

	treeres = treeres*flux
	loopres = loopres*flux + treeres

#ifdef BREMSSTRAHLUNG
	loopres = loopres + treeres*SoftPhotonFactor(Ecms)
#endif

#ifdef WF_RENORMALIZATION
	loopres = loopres + treeres*WF_RENORMALIZATION
#endif
	return

999	continue
	return 1
	end


************************************************************************
** IntegratedCS integrates DifferentialCS over the scattering angle.

	subroutine IntegratedCS(treeres, loopres, Ecms, *)
	implicit none
	double precision treeres, loopres, Ecms

#include "kin.h"

	double precision gauss_x(GAUSSPOINTS/2)
	double precision gauss_w(GAUSSPOINTS/2)

	integer i
	double precision x, w, w1, w2
	double precision treeamp, loopamp, treeF, treeB, loopF, loopB
	double precision lower, upper, mid1, mid2, range
	parameter ( lower = CUTOFF, upper = pi - CUTOFF,
     +    mid1 = .5D0*(lower + pi/2),
     +    mid2 = .5D0*(upper + pi/2),
     +    range = .25D0*(upper - lower) )

#include "gauss.F"

* Integration is done in two intervals separately:
* [CUTOFF, pi/2] and [pi/2, pi - CUTOFF].
* In this way, we get the F-B asymmetry for free, even
* though the number of sampling points is the same,
* e.g. instead of 32 over the whole interval we can do
* with 16 over each half.

	treeF = 0
	loopF = 0
	treeB = 0
	loopB = 0

	do i = 1, GAUSSPOINTS/2
	  x = range*gauss_x(i)
	  w = range*2*pi*IDENTICALFACTOR*gauss_w(i)
*			*2*pi from phi integration
*			possibly /2 because of identical particles
	  w1 = w*sin(mid1 - x)
	  w2 = w*sin(mid1 + x)
	  call DifferentialCS(treeamp, loopamp, Ecms, mid1 - x, *999)
	  treeF = treeF + w1*treeamp
	  loopF = loopF + w1*loopamp
	  call DifferentialCS(treeamp, loopamp, Ecms, mid1 + x, *999)
	  treeF = treeF + w2*treeamp
	  loopF = loopF + w2*loopamp
	  call DifferentialCS(treeamp, loopamp, Ecms, mid2 - x, *999)
	  treeB = treeB + w2*treeamp
	  loopB = loopB + w2*loopamp
	  call DifferentialCS(treeamp, loopamp, Ecms, mid2 + x, *999)
	  treeB = treeB + w1*treeamp
	  loopB = loopB + w1*loopamp
	enddo

	treeres = treeF + treeB
	loopres = loopF + loopB

#ifdef FB_ASYMMETRY
	if(treeres .eq. 0) then
	  treeres = 999D0
	else
	  treeres = (treeF - treeB)/treeres
	endif

	if(loopres .eq. 0) then
	  loopres = 999D0
	else
	  loopres = (loopF - loopB)/loopres
	endif
#endif
	return

999	continue
	return 1
	end


************************************************************************
** The Kronecker delta.

	integer function Delta(i, j)
	implicit none
	integer i, j

	if(i .eq. j) then
	  Delta = 1
	else
	  Delta = 0
	endif
	end


************************************************************************
** MomSquare is a function used to calculate Mandelstam variables
** which in general are the squared difference of two momenta. If the
** squared sum is needed, the second argument must be passed as a
** negative integer. Since MomSquare is called only on momentum
** vectors, all operations are done only on the real part.

	double precision function MomSquare(i, j)
	implicit none
	integer i, j

#include "kin.h"

	if(j .lt. 0) then
	  j = -j
	  MomSquare =
     +      (dble(vec(0, i)) + dble(vec(0, j)))**2 -
     +      (dble(vec(1, i)) + dble(vec(1, j)))**2 -
     +      (dble(vec(2, i)) + dble(vec(2, j)))**2 -
     +      (dble(vec(3, i)) + dble(vec(3, j)))**2
	else
	  MomSquare =
     +      (dble(vec(0, i)) - dble(vec(0, j)))**2 -
     +      (dble(vec(1, i)) - dble(vec(1, j)))**2 -
     +      (dble(vec(2, i)) - dble(vec(2, j)))**2 -
     +      (dble(vec(3, i)) - dble(vec(3, j)))**2
	endif
	end


************************************************************************
** Pair calculates the scalar product of two four-vectors.

	double complex function Pair(i, j)
	implicit none
	integer i, j

#include "kin.h"

	Pair =
     +    vec(0, i)*vec(0, j) -
     +    vec(1, i)*vec(1, j) -
     +    vec(2, i)*vec(2, j) -
     +    vec(3, i)*vec(3, j)
	end


************************************************************************
** Eps calculates -I times the Levi-Civita tensor contracted with
** four four-vectors.

	double complex function Eps(i, j, k, l)
	implicit none
	integer i, j, k, l

#include "kin.h"

	double complex x01, x02, x03, x12, x13, x23

	double complex Cross
	integer mu_, nu_

	Cross(mu_, nu_) =
     +    vec(mu_, k)*vec(nu_, l) - vec(nu_, k)*vec(mu_, l)

	x01 = Cross(0, 1)
	x02 = Cross(0, 2)
	x03 = Cross(0, 3)
	x12 = Cross(1, 2)
	x13 = Cross(1, 3)
	x23 = Cross(2, 3)

	Eps =
     +    vec(0, i)*(vec(1, j)*x23 - vec(2, j)*x13 + vec(3, j)*x12) -
     +    vec(1, i)*(vec(0, j)*x23 - vec(2, j)*x03 + vec(3, j)*x02) +
     +    vec(2, i)*(vec(0, j)*x13 - vec(1, j)*x03 + vec(3, j)*x01) -
     +    vec(3, i)*(vec(0, j)*x12 - vec(1, j)*x02 + vec(2, j)*x01)

* Note: although Eps is defined as -I*(Levi-Civita tensor) in FormCalc,
* we have to multiply with I (not -I) here because we've neglected
* the metric in the summation above.

	Eps = cI*Eps
	end


#ifdef BREMSSTRAHLUNG

************************************************************************
** SoftPhotonFactor adds up the contributions from soft-photon emission
** off the different legs. The result then has to be multiplied with
** the Born cross-section.

	double precision function SoftPhotonFactor(Ecms)
	implicit none
	double precision Ecms

#include "kin.h"

	double precision PhaseSpaceIntegral
	external PhaseSpaceIntegral

	integer i, j

	double precision loglambda
	common /infrared/ loglambda

	double precision ch1, ch2, ch3, ch4, relsign(4), r
	parameter (ch1 = CHARGE1, ch2 = CHARGE2,
     +    ch3 = -(CHARGE3), ch4 = -(CHARGE4))
	data relsign / ch1, ch2, ch3, ch4 /

	loglambda = log(4*(ESOFTMAX)**2/getlambda())

	SoftPhotonFactor = 0
	do i = 1, 4
	  do j = i, 4
	    r = relsign(i)*relsign(j)
	    if(r .ne. 0) SoftPhotonFactor = SoftPhotonFactor +
     +        r*PhaseSpaceIntegral(i, j)
	  enddo
	enddo
	SoftPhotonFactor = -4*pi*Alfa/(2*pi)**3*SoftPhotonFactor
	end


************************************************************************
** PhaseSpaceIntegral calculates the integral over the soft photon's
** phase space. See A. Denner, Fortschr. d. Physik, 41 (1993) 4.

	double precision function PhaseSpaceIntegral(i, j)
	implicit none
	integer i, j

#include "kin.h"

	double complex Pair
	double precision Li2
	external Pair, Li2

	integer ii, jj
	double precision ki0, ki, kj0, kj, kikj
	double precision alpha, vnum, v

	double precision loglambda
	common /infrared/ loglambda

	double precision AbsVec
	integer k_

	AbsVec(k_) = sqrt(dble(vec(1, k_))**2 +
     +    dble(vec(2, k_))**2 + dble(vec(3, k_))**2)

	ii = 4*i - 3
	ki0 = dble(vec(0, ii))
	ki = AbsVec(ii)

	if(i .eq. j) then
	  PhaseSpaceIntegral =
     +      pi*(loglambda + ki0/ki*log((ki0 - ki)/(ki0 + ki)))
	  return
	endif

	jj = 4*j - 3
	kj0 = dble(vec(0, jj))
	kj = AbsVec(jj)
	kikj = dble(Pair(ii, jj))

	v = sqrt(kikj**2 - mass(i)*mass(j))
	alpha = (kikj + v)/mass(i)
	if((alpha*ki0 - kj0)/kj0 .lt. 0) alpha = (kikj - v)/mass(i)

	vnum = alpha**2*mass(i) - mass(j)
	v = .5D0*vnum/(alpha*ki0 - kj0)

	PhaseSpaceIntegral = 4*pi*alpha*kikj/vnum * (
     +    .5D0*log(alpha**2*mass(i)/mass(j))*loglambda +
     +    .25D0*(log((ki0 - ki)/(ki0 + ki)) -
     +           log((kj0 - kj)/(kj0 + kj))) +
     +    Li2(1 - (ki0 + ki)*alpha/v) - Li2(1 - (kj0 + kj)/v) +
     +    Li2(1 - (ki0 - ki)*alpha/v) - Li2(1 - (kj0 - kj)/v) )

	end

#endif


************************************************************************
** Li2 calls the dilogarithm function of FF, ffzxdl

	double precision function Li2(x)
	implicit none
	double precision x

#include "kin.h"

	double precision pi12
	parameter (pi12 = pi**2/12D0)

	double complex cli2, dummy
	integer ier, ipi12

	ier = 0
	call ffzxdl(cli2, ipi12, dummy, x, 1, ier)
	Li2 = dble(cli2) + ipi12*pi12
	end
