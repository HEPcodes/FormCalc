* num.F
* this is the driver program for Fortran code produced by
* WriteSquaredME
* this file is part of FormCalc
* last modified 22 Jan 03 th

* Note: this file actually supplies only utility functions.
* The main program and all other parts of the code that depend
* on the kinematics are included by process.h (e.g. 2to2.F).


#define SCALAR 0
#define FERMION 1
#define PHOTON 2
#define VECTOR 3

#include "prefix.h"
#include "process.h"

* In the function arguments here, vector indices
* i, j, etc. specify the vector itself (as in k_i, k_j, etc.),
* a, b, etc. specify the vector's position in the vec array.
* The macros k, e, ec, s (defined in num.h) convert from i, j, ...
* to a, b, ...


************************************************************************
** ParseCommandline reads the polarizations and the energy range from
** the command line.  From these values, it constructs a base file name.

	subroutine ParseCommandline(sqrtSfrom, sqrtSto)
	implicit none
	double precision sqrtSfrom, sqrtSto

#include COMMON

	integer iargc, lnblnk
	external iargc, lnblnk

	integer argc, i, j
	character*200 argv
	character*LEGS polstr
	character c

	integer type(LEGS)
	data type /TYPES/

	argc = iargc()
	if( argc .lt. LEGS + 1 .or. argc .gt. LEGS + 2 ) then
999	  print *, "USAGE:"
	  print *, "to calculate a differential cross-section:"
	  print *, "  num p1 p2 ... pn sqrts"
	  print *, "to calculate an integrated cross-section:"
	  print *, "  num p1 p2 ... pn sqrts_from sqrts_to"
	  print *, "where the polarizations p1...pn may be one of"
	  print *, "  0  for longitudinal   polarization"
	  print *, "  +  for right-circular      ''"
	  print *, "  -  for left-circular       ''"
	  print *, "  t  for transverse          ''"
	  print *, "  u  for unpolarized particle"
	  stop
	endif

	avgfac = COLOURFACTOR

	do i = 1, LEGS
	  if( type(i) .eq. SCALAR ) then
	    c = 'S'
	    bpol(i) = 0
	    epol(i) = 0
	  else
	    call getarg(i, argv)
	    c = argv(1:1)
	    if( type(i) .eq. FERMION ) then
	      if( c .eq. 'u' ) then
	        j = 0
	        if( i .gt. LEGS_IN ) avgfac = avgfac*2
	        c = 'U'
	      else if( c .eq. '+' ) then
	        j = 1
	      else if( c .eq. '-' ) then
	        j = -1
	      else
	        stop "Invalid polarization"
	      endif
	      bpol(i) = j
	      epol(i) = j
	    else
	      if( c .eq. 'u' ) then
	        if( type(i) .eq. PHOTON )
     &            stop "No longitudinal polarizations allowed"
	        bpol(i) = 0
	        epol(i) = 2
	        if( i .le. LEGS_IN ) avgfac = avgfac/3D0
	        c = 'U'
	      else if( c .eq. 't' ) then
	        bpol(i) = 1
	        epol(i) = 2
	        if( i .le. LEGS_IN ) avgfac = avgfac/2D0
	        c = 'T'
	      else
	        if( c .eq. 'l' .or. c .eq. '0' ) then
	          if( type(i) .eq. PHOTON )
     &              stop "No longitudinal polarizations allowed"
	          j = 0
	          c = 'L'
	        else if( c .eq. '+' ) then
	          j = 1
	        else if( c .eq. '-' ) then
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

	call getarg(LEGS + 1, argv)
	read(argv, 22, err=999) sqrtSfrom
	if( argc .gt. LEGS + 1 ) then
	  call getarg(LEGS + 2, argv)
	  read(argv, 22, err=999) sqrtSto
23	  format("tot.pol=", A, ".E=", I5.5, "-", I5.5)
	  write(outfile, 23) polstr, int(sqrtSfrom), int(sqrtSto)
	else
	  sqrtSto = -1
24	  format("diff.pol=", A, ".E=", I5.5)
	  write(outfile, 24) polstr, int(sqrtSfrom)
	endif

	call getarg(0, argv)
	outfile = argv(1:lnblnk(argv)) // "-" // outfile

#ifdef PARALLEL
25	format(I)
	call getenv("NCPUS", argv)
	read(argv, 25) cpus
	if( cpus .lt. 1 ) cpus = 1
#endif
	end


************************************************************************
** ThreeMom computes the length of \vec p_b in the frame in which
** \vec p_a + \vec p_b vanishes.  S = (p_a + p_b)^2.

	double precision function ThreeMom(sqrtS, ma2, mb2)
	implicit none
	double precision sqrtS, ma2, mb2

	ThreeMom = sqrt((.5D0*(sqrtS - (ma2 - mb2)/sqrtS))**2 - mb2)
	end


************************************************************************
** VecSet constructs the CMS representation of a momentum and the
** three corresponding polarization vectors.  The arguments:
** a: the index of the momentum
** m, p: mass and three-momentum of the particle
** ex, ey, ez: the unit three-vector of the momentum

	subroutine VecSet(a, m, p, ex, ey, ez)
	implicit none
	integer a
	double precision m, p, ex, ey, ez

#include COMMON

	double precision E, sinth
	integer b, rho

	E = sqrt(p**2 + m**2)

	b = a
	vec(0, b) = E
	vec(1, b) = p*ex
	vec(2, b) = p*ey
	vec(3, b) = p*ez

	b = b + 1
	if( m .ne. 0 ) then
	  vec(0, b) = p/m
	  vec(1, b) = E/m*ex
	  vec(2, b) = E/m*ey
	  vec(3, b) = E/m*ez
	endif

	b = b + 1
	vec(0, b) = 0
	sinth = 1 - ez**2
	if( sinth .ne. 0 ) then
	  sinth = sqrt(sinth)
	  vec(1, b) = (ex*ez - ey*cI)/sinth/sqrt2
	  vec(2, b) = (ey*ez + ex*cI)/sinth/sqrt2
	  vec(3, b) = -sinth/sqrt2
	else
	  vec(1, b) = 1/sqrt2
	  vec(2, b) = ez/sqrt2*cI
	  vec(3, b) = 0
	endif

	b = b + 1
	do rho = 0, 3
	  vec(rho, b) = dconjg(vec(rho, b - 1))
	enddo

	do b = a, a + 3
	  do rho = 0, 3
	    vec(rho, -b) = dconjg(vec(rho, b))
	  enddo
	enddo
	end


************************************************************************
** The Kronecker delta

	integer function Delta(i, j)
	implicit none
	integer i, j

	if( i .eq. j ) then
	  Delta = 1
	else
	  Delta = 0
	endif
	end


************************************************************************
** SInvariant computes an s-type invariant, i.e. s_{ij} = (k_i + k_j)^2.

	double precision function SInvariant(a, b)
	implicit none
	integer a, b

#include COMMON

	SInvariant =
     &    (dble(vec(0, a)) + dble(vec(0, b)))**2 -
     &    (dble(vec(1, a)) + dble(vec(1, b)))**2 -
     &    (dble(vec(2, a)) + dble(vec(2, b)))**2 -
     &    (dble(vec(3, a)) + dble(vec(3, b)))**2
	end


************************************************************************
** TInvariant computes an t-type invariant, i.e. t_{ij} = (k_i - k_j)^2.

	double precision function TInvariant(a, b)
	implicit none
	integer a, b

#include COMMON

	TInvariant =
     &    (dble(vec(0, a)) - dble(vec(0, b)))**2 -
     &    (dble(vec(1, a)) - dble(vec(1, b)))**2 -
     &    (dble(vec(2, a)) - dble(vec(2, b)))**2 -
     &    (dble(vec(3, a)) - dble(vec(3, b)))**2
	end


************************************************************************
** Pair calculates the scalar product of two four-vectors.

	double complex function Pair(a, b)
	implicit none
	integer a, b

#include COMMON

	Pair =
     &    vec(0, a)*vec(0, b) -
     &    vec(1, a)*vec(1, b) -
     &    vec(2, a)*vec(2, b) -
     &    vec(3, a)*vec(3, b)
	end


************************************************************************
** Eps calculates -I times the Levi-Civita tensor contracted with
** four four-vectors.

	double complex function Eps(a, b, c, d)
	implicit none
	integer a, b, c, d

#include COMMON

	double complex x01, x02, x03, x12, x13, x23

	x01 = vec(0, c)*vec(1, d) - vec(1, c)*vec(0, d)
	x02 = vec(0, c)*vec(2, d) - vec(2, c)*vec(0, d)
	x03 = vec(0, c)*vec(3, d) - vec(3, c)*vec(0, d)
	x12 = vec(1, c)*vec(2, d) - vec(2, c)*vec(1, d)
	x13 = vec(1, c)*vec(3, d) - vec(3, c)*vec(1, d)
	x23 = vec(2, c)*vec(3, d) - vec(3, c)*vec(2, d)

	Eps =
     &    vec(0, a)*(vec(1, b)*x23 - vec(2, b)*x13 + vec(3, b)*x12) -
     &    vec(1, a)*(vec(0, b)*x23 - vec(2, b)*x03 + vec(3, b)*x02) +
     &    vec(2, a)*(vec(0, b)*x13 - vec(1, b)*x03 + vec(3, b)*x01) -
     &    vec(3, a)*(vec(0, b)*x12 - vec(1, b)*x02 + vec(2, b)*x01)

* Note: although Eps is defined as -I*(Levi-Civita tensor) in FormCalc,
* we have to multiply with I (not -I) here because we've neglected
* the metric in the summation above.

	Eps = cI*Eps
	end


#ifdef BREMSSTRAHLUNG

************************************************************************
** SoftPhotonFactor adds up the contributions from soft-photon emission
** off the different legs.  The result then has to be multiplied with
** the Born cross-section.

	double precision function SoftPhotonFactor()
	implicit none

#include COMMON

	double precision PhaseSpaceIntegral
	external PhaseSpaceIntegral

	integer i, j

	double precision loglambda
	common /infrared/ loglambda

	double precision charge(LEGS), r

* not a very nice solution, the following, but necessary
* because there may be no "/" in a constant in the data
* statement (such as in 2/3D0)
#if LEGS == 3
	double precision ch1, ch2, ch3
	parameter (ch1 = CHARGE1, ch2 = CHARGE2, ch3 = CHARGE3)
	data charge / ch1, ch2, ch3 /
#elif LEGS == 4
	double precision ch1, ch2, ch3, ch4
	parameter (ch1 = CHARGE1, ch2 = CHARGE2,
     &             ch3 = CHARGE3, ch4 = CHARGE4)
	data charge / ch1, ch2, ch3, ch4 /
#elif LEGS == 5
	double precision ch1, ch2, ch3, ch4, ch5
	parameter (ch1 = CHARGE1, ch2 = CHARGE2,
     &             ch3 = CHARGE3, ch4 = CHARGE4, ch5 = CHARGE5)
	data charge / ch1, ch2, ch3, ch4, ch5 /
#else
	error: No definition for LEGS external legs
#endif

	loglambda = log(4*(ESOFTMAX)**2/getlambda())

	SoftPhotonFactor = 0
	do i = 1, LEGS
	  do j = i, LEGS
	    r = charge(i)*charge(j)
	    if( r .ne. 0 ) SoftPhotonFactor = SoftPhotonFactor +
     &        sign(1, LEGS_IN - i)*sign(1, LEGS_IN - j)*
     &        r*PhaseSpaceIntegral(i, j)
	  enddo
	enddo
	SoftPhotonFactor = -4*pi*Alfa/(2*pi)**3*SoftPhotonFactor
	end


************************************************************************
** PhaseSpaceIntegral calculates the integral over the soft photon's
** phase space.  See A. Denner, Fortschr. d. Physik, 41 (1993) 4.

	double precision function PhaseSpaceIntegral(i, j)
	implicit none
	integer i, j

#include COMMON
#include "num.h"

	double precision ki0, ki, ki2, kidelta, kilog
	double precision kj0, kj, kj2, kjdelta, kjlog
	double precision kikj, alpha, d, vi

	double precision loglambda
	common /infrared/ loglambda

	call GetMomSpecs(i, ki0, ki, kidelta, kilog)

	if( i .eq. j ) then
	  PhaseSpaceIntegral = pi*(loglambda + ki0/ki*kilog)
	  return
	endif

	call GetMomSpecs(j, kj0, kj, kjdelta, kjlog)
	kikj = dble(Pair(k(i), k(j)))
	ki2 = mass2(i)
	kj2 = mass2(j)

	d = sqrt(kikj**2 - ki2*kj2)
	alpha = (kikj + d)/ki2
	if( sign(kj0, alpha*ki0 - kj0) .ne. kj0 ) alpha = (kikj - d)/ki2

	d = 1/(alpha**2*ki2 - kj2)
	vi = 2*(alpha*ki0 - kj0)*d

	PhaseSpaceIntegral = 4*pi*alpha*kikj*d * (
     &    .5D0*log(ki2/kj2*alpha**2)*loglambda +
     &    .25D0*(kilog**2 - kjlog**2) +
     &    Li2(1 - (ki0 + ki)*alpha*vi) + Li2(1 - kidelta*alpha*vi) -
     &    Li2(1 - (kj0 + kj)*vi) - Li2(1 - kjdelta*vi) )

	end


************************************************************************
** GetMomSpecs computes the energy (ki0) and three-momentum (ki) for
** a given momentum vector i.  It also returns the difference (ki0 - ki),
** which needs to be computed using a series expansion if the mass is
** very small to achieve sufficient accuracy.

	subroutine GetMomSpecs(i, ki0, ki, kidelta, kilog)
	implicit none
	integer i
	double precision ki0, ki, kidelta, kilog

#include COMMON
#include "num.h"

	integer a

	a = k(i)
	ki0 = dble(vec(0, a))
	ki = sqrt( dble(vec(1, a))**2 + dble(vec(2, a))**2 +
     &    dble(vec(3, a))**2 )

	kidelta = ki0 - ki
	if( kidelta .lt. 1D-8*ki ) kidelta = .5D0*mass2(i)/ki

	kilog = log(kidelta/(ki0 + ki))
	end

#endif BREMSSTRAHLUNG


************************************************************************
** Li2 calls the dilogarithm function of FF, ffzxdl.

	double precision function Li2(x)
	implicit none
	double precision x

#include COMMON

	double precision pi12
	parameter (pi12 = pi**2/12D0)

	double complex cli2, dummy
	integer ier, ipi12

	ier = 0
	call ffzxdl(cli2, ipi12, dummy, x, 1, ier)
	Li2 = dble(cli2) + ipi12*pi12
	end

