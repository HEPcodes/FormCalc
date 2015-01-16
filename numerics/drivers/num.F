* num.F
* this is the driver program for a 2 -> 2 process
* for code produced by NumPrep.m
* this file is part of FormCalc
* last modified 14 May 00 th


#include "process.h"
#include "squared_me.F"

#define SCALAR 0
#define FERMION 1
#define PHOTON 2
#define VECTOR 3


* The momenta and polarization vectors are defined as
* (this is the order NumPrep expects them in):

* vec(1)     = k1     = (E1, 0, 0, -Pin)
* vec(2)     = e1(0)  = (-Pin, 0, 0, E1)/MASS1
* vec(3,4)   = e1(+-) = (0, -1, +-i, 0)/sqrt(2)

#define k1(i) vec(i, 1)
#define e1(pol, i) vec(i, 2 + pol)

* vec(5)     = k2     = (E2, 0, 0, Pin)
* vec(6)     = e2(0)  = (-Pin, 0, 0, -E2)/MASS2
* vec(7,8)   = e2(+-) = (0, 1, +-i, 0)/sqrt(2)

#define k2(i) vec(i, 5)
#define e2(pol, i) vec(i, 6 + pol)

* vec(9)     = k3     = (E3, -Pout sin(th), 0, -Pout cos(th))
* vec(10)    = e3(0)  = (-Pout, E3 sin(th), 0, E3 cos(th))/MASS3
* vec(11,12) = e3(+-) = (0, -cos(th), -+i, sin(th))/sqrt(2)

#define k3(i) vec(i, 9)
#define e3(pol, i) vec(i, 10 + pol)

* vec(13)    = k4     = (E4, Pout sin(th), 0, Pout cos(th))
* vec(14)    = e4(0)  = (-Pout, -E4 sin(th), 0, -E4 cos(th))/MASS4
* vec(15,16) = e4(+-) = (0, cos(th), -+i, -sin(th))/sqrt(2)

#define k4(i) vec(i, 13)
#define e4(pol, i) vec(i, 14 + pol)


************************************************************************
** THE MAIN PROGRAM STARTS HERE:

	program num
	implicit none

#include "kin.h"

	double precision Ecms, Ecmsfrom, Ecmsto, theta
	double precision treeres, loopres

	call ParseCommandline(Ecmsfrom, Ecmsto)

	call ffini

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

	double precision Ecms2, flux, EE(4), Pin, Pout
	double precision threshold, costh
	integer component, n
	logical reini
	save Ecms2, flux, EE, Pin, Pout

	double precision kallen, x1, x2

	double precision sqrt2, hbar_c2
	parameter (sqrt2 = 1.4142135623730950488016887242096981D0)
	parameter (hbar_c2 = 3.8937966D8)
*	  = \hbar c^2 in picobarn

	double precision prevEcms
	data prevEcms /-1/

	kallen(x1, x2) =
     +    sqrt(Ecms2*(Ecms2 - 2*(x1 + x2)) + (x1 - x2)**2)

26	format(VARFMT)
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
	    write(11, 26) VAROUT
	    print 26, VAROUT
	  endif

#ifdef MUDIM
	  call setmudim(MUDIM)
#endif
#ifdef DELTA
	  call setdelta(DELTA)
#endif
#ifdef LAMBDA
	  call setlambda(LAMBDA)
#endif

	  mass(1) = MASS1**2
	  mass(2) = MASS2**2
	  mass(3) = MASS3**2
	  mass(4) = MASS4**2
	  call calc_renconst

	  prevEcms = Ecms
	  Ecms2 = Ecms**2

	  flux = avgfac*COLOURFACTOR*
     +      hbar_c2/(64*pi**2*Ecms2)*
     +      kallen(mass(3), mass(4))/kallen(mass(1), mass(2))

	  Pin = (Ecms2 + mass(2) - mass(1))**2/4D0/Ecms2 - mass(2)
	  Pout = (Ecms2 + mass(4) - mass(3))**2/4D0/Ecms2 - mass(4)
	  EE(1) = sqrt(Pin + mass(1))
	  EE(2) = sqrt(Pin + mass(2))
	  EE(3) = sqrt(Pout + mass(3))
	  EE(4) = sqrt(Pout + mass(4))
	  Pin = sqrt(Pin)
	  Pout = sqrt(Pout)

* some apps need the relative velocity of the final-state particles
	  beta = Pout/(.5D0*Ecms)

	  k1(1) = EE(1)
	  k1(2) = 0
	  k1(3) = 0
	  k1(4) = -Pin

	  if(MASS1 .ne. 0) then
	    e1(0,1) = -Pin/MASS1
	    e1(0,2) = 0
	    e1(0,3) = 0
	    e1(0,4) = EE(1)/MASS1
	  endif

	  e1(1,1) = 0
	  e1(1,2) = -1/sqrt2
	  e1(1,3) = dcmplx(0D0, 1/sqrt2)
	  e1(1,4) = 0

	  e1(2,1) = 0
	  e1(2,2) = -1/sqrt2
	  e1(2,3) = -dcmplx(0D0, 1/sqrt2)
	  e1(2,4) = 0

	  k2(1) = EE(2)
	  k2(2) = 0
	  k2(3) = 0
	  k2(4) = Pin

	  if(MASS2 .ne. 0) then
	    e2(0,1) = -Pin/MASS2
	    e2(0,2) = 0
	    e2(0,3) = 0
	    e2(0,4) = -EE(2)/MASS2
	  endif

	  e2(1,1) = 0
	  e2(1,2) = 1/sqrt2
	  e2(1,3) = dcmplx(0D0, 1/sqrt2)
	  e2(1,4) = 0

	  e2(2,1) = 0
	  e2(2,2) = 1/sqrt2
	  e2(2,3) = -dcmplx(0D0, 1/sqrt2)
	  e2(2,4) = 0

	  do n = 1, 8
	    do component = 1, 4
	      vec(component, -n) = dconjg(vec(component, n))
	    enddo
	  enddo
	endif

	sinth = sin(theta)
	costh = cos(theta)

	k3(1) = EE(3)
	k3(2) = -Pout*sinth
	k3(3) = 0
	k3(4) = -Pout*costh

	if(MASS3 .ne. 0) then
	  e3(0,1) = -Pout/MASS3
	  e3(0,2) = EE(3)/MASS3*sinth
	  e3(0,3) = 0
	  e3(0,4) = EE(3)/MASS3*costh
	endif

	e3(1,1) = 0
	e3(1,2) = -costh/sqrt2
	e3(1,3) = -dcmplx(0D0, 1/sqrt2)
	e3(1,4) = sinth/sqrt2

	e3(2,1) = 0
	e3(2,2) = -costh/sqrt2
	e3(2,3) = dcmplx(0D0, 1/sqrt2)
	e3(2,4) = sinth/sqrt2

	k4(1) = EE(4)
	k4(2) = Pout*sinth
	k4(3) = 0
	k4(4) = Pout*costh

	if(MASS4 .ne. 0) then
	  e4(0,1) = -Pout/MASS4
	  e4(0,2) = -EE(4)/MASS4*sinth
	  e4(0,3) = 0
	  e4(0,4) = -EE(4)/MASS4*costh
	endif

	e4(1,1) = 0
	e4(1,2) = costh/sqrt2
	e4(1,3) = -dcmplx(0D0, 1/sqrt2)
	e4(1,4) = -sinth/sqrt2

	e4(2,1) = 0
	e4(2,2) = costh/sqrt2
	e4(2,3) = dcmplx(0D0, 1/sqrt2)
	e4(2,4) = -sinth/sqrt2

	do n = 9, 16
	  do component = 1, 4
	    vec(component, -n) = dconjg(vec(component, n))
	  enddo
	enddo

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
	double precision x, w
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
	  call DifferentialCS(treeamp, loopamp, Ecms, mid1 + x, *999)
	  treeF = treeF + sinth*w*treeamp
	  loopF = loopF + sinth*w*loopamp
	  call DifferentialCS(treeamp, loopamp, Ecms, mid1 - x, *999)
	  treeF = treeF + sinth*w*treeamp
	  loopF = loopF + sinth*w*loopamp
	  call DifferentialCS(treeamp, loopamp, Ecms, mid2 + x, *999)
	  treeB = treeB + sinth*w*treeamp
	  loopB = loopB + sinth*w*loopamp
	  call DifferentialCS(treeamp, loopamp, Ecms, mid2 - x, *999)
	  treeB = treeB + sinth*w*treeamp
	  loopB = loopB + sinth*w*loopamp
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

	integer function IndexDelta(i, j)
	implicit none
	integer i, j

	if(i .eq. j) then
	  IndexDelta = 1
	else
	  IndexDelta = 0
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
     +      (dble(vec(1,i)) + dble(vec(1,j)))**2 -
     +      (dble(vec(2,i)) + dble(vec(2,j)))**2 -
     +      (dble(vec(3,i)) + dble(vec(3,j)))**2 -
     +      (dble(vec(4,i)) + dble(vec(4,j)))**2
	else
	  MomSquare =
     +      (dble(vec(1,i)) - dble(vec(1,j)))**2 -
     +      (dble(vec(2,i)) - dble(vec(2,j)))**2 -
     +      (dble(vec(3,i)) - dble(vec(3,j)))**2 -
     +      (dble(vec(4,i)) - dble(vec(4,j)))**2
	endif
	end


************************************************************************
** Pair calculates the scalar product of two four-vectors.

	double complex function Pair(i, j)
	implicit none
	integer i, j

#include "kin.h"

	Pair =
     +    vec(1,i)*vec(1,j) -
     +    vec(2,i)*vec(2,j) -
     +    vec(3,i)*vec(3,j) -
     +    vec(4,i)*vec(4,j)
	end


************************************************************************
** Eps calculates -I times the Levi-Civita tensor contracted with
** four four-vectors.

	double complex function Eps(i, j, k, l)
	implicit none
	integer i, j, k, l

#include "kin.h"

	double complex kl12, kl13, kl14, kl23, kl24, kl34

	kl12 = vec(1,k)*vec(2,l) - vec(2,k)*vec(1,l)
	kl13 = vec(1,k)*vec(3,l) - vec(3,k)*vec(1,l)
	kl14 = vec(1,k)*vec(4,l) - vec(4,k)*vec(1,l)
	kl23 = vec(2,k)*vec(3,l) - vec(3,k)*vec(2,l)
	kl24 = vec(2,k)*vec(4,l) - vec(4,k)*vec(2,l)
	kl34 = vec(3,k)*vec(4,l) - vec(4,k)*vec(3,l)

	Eps =
     +    vec(1,i)*(vec(2,j)*kl34 - vec(3,j)*kl24 + vec(4,j)*kl23) -
     +    vec(2,i)*(vec(1,j)*kl34 - vec(3,j)*kl14 + vec(4,j)*kl13) +
     +    vec(3,i)*(vec(1,j)*kl24 - vec(2,j)*kl14 + vec(4,j)*kl12) -
     +    vec(4,i)*(vec(1,j)*kl23 - vec(2,j)*kl13 + vec(3,j)*kl12)

* Note: although Eps is defined as -I*(Levi-Civita tensor) in FormCalc,
* we have to multiply with I (not -I) here because we've neglected
* the metric in the summation above.

	Eps = dcmplx(0, 1)*Eps
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

	double precision ch1, ch2, ch3, ch4, relsign(4)
	parameter (ch1 = CHARGE1, ch2 = CHARGE2,
     +    ch3 = -(CHARGE3), ch4 = -(CHARGE4))
	data relsign / ch1, ch2, ch3, ch4 /

	loglambda = log(4*(ESOFTMAX)**2/getlambda())

	SoftPhotonFactor = 0
	do i = 1, 4
	  do j = i, 4
	    SoftPhotonFactor = SoftPhotonFactor +
     +        relsign(i)*relsign(j)*PhaseSpaceIntegral(i, j)
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
	integer k

	AbsVec(k) = sqrt(dble(vec(2,k))**2 +
     +    dble(vec(3,k))**2 + dble(vec(4,k))**2)

	ii = 4*i - 3
	ki0 = dble(vec(1,ii))
	ki = AbsVec(ii)

	if(i .eq. j) then
	  PhaseSpaceIntegral =
     +      pi*(loglambda + ki0/ki*log((ki0 - ki)/(ki0 + ki)))
	  return
	endif

	jj = 4*j - 3
	kj0 = dble(vec(1,jj))
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

#endif

