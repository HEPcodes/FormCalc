* num.F
* this is the driver program for Fortran code
* produced by WriteSquaredME
* this file is part of FormCalc
* last modified 4 Mar 02 th

* Note: this file actually supplies only utility functions.
* The main program and all other parts of the code that depend
* on the kinematics are included by process.h (e.g. 2to2.F).


#define SCALAR 0
#define FERMION 1
#define PHOTON 2
#define VECTOR 3

#include "prefix.h"
#include "process.h"


************************************************************************
** ParseCommandline reads the polarizations and the energy range from
** the command line. From these values, it constructs a file name and
** opens the file for writing as file #11.

	subroutine ParseCommandline(Ecmsfrom, Ecmsto)
	implicit none
	double precision Ecmsfrom, Ecmsto

#include COMMON

	integer*4 iargc
	external iargc

	integer argc, i, j
	character*200 argv, outf
	character*LEGS polstr
	character c

	integer type(LEGS)
	data type /TYPES/

	argc = iargc()
	if(argc .lt. LEGS + 1 .or. argc .gt. LEGS + 2) then
	  print *, "USAGE:"
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
23	format("diff.pol=", A, ".E=", I5.5)
24	format("tot.pol=", A, ".E=", I5.5, "-", I5.5)

	call getarg(LEGS + 1, argv)
	read(argv, 22) Ecmsfrom
	if(argc .gt. LEGS + 1) then
	  call getarg(LEGS + 2, argv)
	  read(argv, 22) Ecmsto
	  write(outf, 24) polstr, int(Ecmsfrom), int(Ecmsto)
	else
	  Ecmsto = -1
	  write(outf, 23) polstr, int(Ecmsfrom)
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

#include COMMON

	double precision E, sinth
	integer rho, nn

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

	do nn = n, n + 3
	  do rho = 0, 3
	    vec(rho, -nn) = dconjg(vec(rho, nn))
	  enddo
	enddo
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

#include COMMON

	if(j .lt. 0) then
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

#include COMMON

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

#include COMMON

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

	double precision function SoftPhotonFactor()
	implicit none

#include COMMON

	double precision PhaseSpaceIntegral
	external PhaseSpaceIntegral

	integer i, j

	double precision loglambda
	common /infrared/ loglambda

	double precision Ecms
	equivalence (Ecms, sqrtS)

	double precision ch1, ch2, ch3, ch4, relsign(4), r
	parameter (ch1 = CHARGE1, ch2 = CHARGE2,
     +             ch3 = CHARGE3, ch4 = CHARGE4)
	data relsign / ch1, ch2, -ch3, -ch4 /

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

#include COMMON

	double complex Pair
	double precision Li2
	external Pair, Li2

	double precision ki0, ki, ki2, kidelta, kilog
	double precision kj0, kj, kj2, kjdelta, kjlog
	double precision kikj, alpha, d, vi

	double precision loglambda
	common /infrared/ loglambda

	call getkk0(i, ki0, ki, kidelta, kilog)

	if(i .eq. j) then
	  PhaseSpaceIntegral = pi*(loglambda + ki0/ki*kilog)
	  return
	endif

	call getkk0(j, kj0, kj, kjdelta, kjlog)
	kikj = dble(Pair(4*i - 3, 4*j - 3))
	ki2 = mass(i)
	kj2 = mass(j)

	d = sqrt(kikj**2 - ki2*kj2)
	alpha = (kikj + d)/ki2
	if(sign(kj0, alpha*ki0 - kj0) .ne. kj0) alpha = (kikj - d)/ki2

	d = 1/(alpha**2*ki2 - kj2)
	vi = 2*(alpha*ki0 - kj0)*d

	PhaseSpaceIntegral = 4*pi*alpha*kikj*d * (
     +    .5D0*log(ki2/kj2*alpha**2)*loglambda +
     +    .25D0*(kilog**2 - kjlog**2) +
     +    Li2(1 - (ki0 + ki)*alpha*vi) + Li2(1 - kidelta*alpha*vi) -
     +    Li2(1 - (kj0 + kj)*vi) - Li2(1 - kjdelta*vi) )

	end


************************************************************************
** getkk0 computes k and k0 for a given momentum. It also returns the
** difference (k0 - k), which needs to be computed using a series
** expansion if the mass is very small to achieve sufficient accuracy.

	subroutine getkk0(i, k0, k, kdelta, klog)
	implicit none
	integer i
	double precision k0, k, kdelta, klog

#include COMMON

	integer ii

	ii = 4*i - 3
	k0 = dble(vec(0, ii))
	k = sqrt( dble(vec(1, ii))**2 + dble(vec(2, ii))**2 +
     +    dble(vec(3, ii))**2 )

	kdelta = k0 - k
	if(kdelta .lt. 1D-8*k) kdelta = .5D0*mass(i)/k

	klog = log(kdelta/(k0 + k))
	end

#endif BREMSSTRAHLUNG


************************************************************************
** Li2 calls the dilogarithm function of FF, ffzxdl

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
