*	num.F
*	this is the driver program for a 2 -> 2 process
*	for code produced by NumPrep.m
*	last modified 23 Jun 99 th

#include "feyn.F"
#include "process.h"
#include "tools.F"

#define SCALAR 0
#define FERMION 1
#define PHOTON 2
#define VECTOR 3
#define FERMIONIC \
  (TYPE1 == FERMION || TYPE2 == FERMION || \
   TYPE3 == FERMION || TYPE4 == FERMION)

#ifdef DEBUG
#define THETASTEP (125D0*degree)
*#define THETASTEP (40D0*degree)
#else
#define THETASTEP (1D0*degree)
#endif

*#define ECMSSTEP max(5.00653422312312D0, dble(S - MH2)*2.000452034D-6)
*		        ^^^ some silly numbers to avoid S===MH2 ^^^
#define ECMSSTEP 10D0


*			   1-4: p1, p2, k1, k2
#define E1START 5
*			   5-7: e1(0, +, -)
#define E2START 8
*			  8-10: e2(0, +, -)
#define E3START 11
*			 11-13: e3(0, +, -)
#define E4START 14
*			 14-16: e4(0, +, -)

* the momenta and polarization vectors are taken as:

* p1     = (E, 0, 0, -p)
* e1(0)  = (-p, 0, 0, E)/MASS1
* e1(+-) = (0, -1, +-i, 0)/sqrt(2)

* p2     = (E, 0, 0, p)
* e2(0)  = (-p, 0, 0, -E)/MASS2
* e2(+-) = (0, 1, +-i, 0)/sqrt(2)

* k1     = (E, -p sin(th), 0, -p cos(th))
* e3(0)  = (-p, E sin(th), 0, E cos(th))/MASS3
* e3(+-) = (0, -cos(th), -+i, sin(th))/sqrt(2)

* k2     = (E, p sin(th), 0, p cos(th))
* e4(0)  = (-p, -E sin(th), 0, -E cos(th))/MASS4
* e4(+-) = (0, cos(th), -+i, -sin(th))/sqrt(2)

#define p1(i) v4(i, 1)
#define p2(i) v4(i, 2)
#define k1(i) v4(i, 3)
#define k2(i) v4(i, 4)
#define e1(pol, i) v4(i, E1START + pol)
#define e2(pol, i) v4(i, E2START + pol)
#define e3(pol, i) v4(i, E3START + pol)
#define e4(pol, i) v4(i, E4START + pol)


************************************************************************
** MAIN PROGRAM STARTS HERE:

	program num
	implicit none
#include "kin.h"
#include "tools.h"

	integer*4 iargc, argc, lnblnk
	external iargc, getarg, lnblnk
	character*40 argv, outf, basename
	double precision Ecmsfrom, Ecmsto

	data type / TYPE1, TYPE2, TYPE3, TYPE4 /

	argc = iargc()
	if(argc .lt. 5 .or. argc .gt. 6) then
	  print *, "USAGE:"
	  print *, "to calculate differential cross-section:"
	  print *, "	num i1 i2 i3 i4 ecms"
	  print *, "to calculate total cross-section:"
	  print *, "	num i1 i2 i3 i4 ecms_from ecms_to"
	  print *, "for polarization i1 + i2 -> i3 + i4"
	  print *, "where i[1-4] may be one of"
	  print *, "  0 or l  for longitudinal   polarization"
	  print *, "    +     for right-circular      ''"
	  print *, "    -     for left-circular       ''"
	  print *, "    t     for transverse          ''"
	  print *, "    u     for unpolarized particle"
	  stop
	endif

	avgfac = 1D0
	call getpol(1)
	call getpol(2)
	call getpol(3)
	call getpol(4)

22	format(E20.0)
23	format("diff.pol=", A4, ".E=", I5.5)
24	format("tot.pol=", A4, ".E=", I5.5, "-", I5.5)

	call getarg(5, argv)
	read(argv, 22) Ecmsfrom
	Ecms = Ecmsfrom
	if(argc .eq. 5) then
	  write(outf, 23) polstr, int(Ecms)
	else
	  call getarg(6, argv)
	  read(argv, 22) Ecmsto
	  write(outf, 24) polstr, int(Ecmsfrom), int(Ecmsto)
	endif
	call getarg(0, basename)
	outf = basename(1:lnblnk(basename)) // "-" // outf

	print *, "output file: ", outf
	print *, "angular cut: ", CUTOFF

	call setup
	open(11, file=outf, status='UNKNOWN')

26	format(3G24.15)
27	format(VARFMT)
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

	MODELINI
	mass(1) = MASS1**2
	mass(2) = MASS2**2
	mass(3) = MASS3**2
	mass(4) = MASS4**2
	call calc_renconst
	write(11, 27) VAROUT
	print 27, VAROUT

	Ecms = max(Ecmsfrom,
     +    MASS1 + MASS2 + .01D0,
     +    MASS3 + MASS4 + .01D0)
	if(Ecms .ne. Ecmsfrom) print *, "Warning: sqrt(s) = ",
     +    Ecmsfrom, " is below threshold."

	if(argc .eq. 5) then
	  call setenergy
	  do th = CUTOFF, pi - CUTOFF, THETASTEP
	    print *, "th =", th
	    call getamp
	    bornamp = bornamp*kinf
#ifdef oneloop
	    loopamp = 999D0
#else
	    loopamp = loopamp*kinf
#endif
	    write(11, 26) th, bornamp, loopamp
	    call flush(11)
	    Cptr = Cptr2
	    Dptr = Dptr2
	  enddo
	else
18	  continue
	  call intgauss
	  bornamp = bornsum*kinf
#ifdef oneloop
	  loopamp = 999D0
	  loopAfb = 999D0
#else
	  loopamp = loopsum*kinf
#endif
	  write(11, 26) Ecms, bornamp, loopamp
*	  write(11, 26) Ecms, bornAfb, loopAfb
	  call flush(11)
	  if(Ecms .lt. Ecmsto) then
	    Ecms = Ecms + ECMSSTEP
	    goto 18
	  endif
	endif

	write(11, 28)
1	continue

	close(11)
	call ffexi
	end

************************************************************************
	subroutine getpol(i)
	implicit none
	integer i
#include "kin.h"

	character c*1, argv*40
	integer j

	integer offset(4)
	data offset / E1START, E2START, E3START, E4START /

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
	      if(i .gt. 2) avgfac = avgfac*2D0
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
     +          stop "No longitudinal polarizations allowed"
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
     +            stop "No longitudinal polarizations allowed"
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
	    bpol(i) = bpol(i) + offset(i)
	    epol(i) = epol(i) + offset(i)
	  endif
	endif
	polstr(i:i) = c
	end

************************************************************************
	subroutine setup
	implicit none
#include "kin.h"

	double complex ds2I, ds2mI
	parameter (ds2I = dcmplx(0D0, ds2),
     +    ds2mI = dcmplx(0D0, -ds2))

	integer i, j

	do i = 1, 4
	  do j = 1, 16 
	    v4(i, j) = dcmplx(0D0)
	  enddo
	enddo

	e1(1, 2) = -ds2
	e1(1, 3) = ds2I
	e1(2, 2) = -ds2
	e1(2, 3) = ds2mI

	e2(1, 2) = ds2
	e2(1, 3) = ds2I
	e2(2, 2) = ds2
	e2(2, 3) = ds2mI

	e3(1, 3) = ds2mI
	e3(2, 3) = ds2I

	e4(1, 3) = ds2mI
	e4(2, 3) = ds2I

	call bcaini
	end

************************************************************************
	subroutine setenergy
	implicit none
#include "tools.h"
#include "kin.h"
#include "rcsm.h"

#ifdef DYSON
	double complex shh
	external shh
#endif

	double precision kallen, ma, mb
	kallen(ma, mb) =
     +    sqrt(S*(S - 2D0*(ma + mb)) + (ma - mb)**2)

	S = Ecms*Ecms
	kinf = hbarc2/(64D0*pi**2*S)*COLORFACTOR*avgfac*
     +    kallen(mass(3), mass(4))/kallen(mass(1), mass(2))
	Pin2 = (S + mass(2) - mass(1))**2/4D0/S - mass(2)
	Pout2 = (S + mass(4) - mass(3))**2/4D0/S - mass(4)
	EE(1) = sqrt(Pin2 + mass(1))
	EE(2) = sqrt(Pin2 + mass(2))
	EE(3) = sqrt(Pout2 + mass(3))
	EE(4) = sqrt(Pout2 + mass(4))
	Pin = sqrt(Pin2)
	Pout = sqrt(Pout2)

#ifdef DYSON
	reso = 1D0/(S - MH2 + shh(S))
#else
	reso = 1D0/(S - MH2)
	resoc = 1D0/(S - MH2 + dcmplx(0D0, 1D0)*GammaHMH)
#endif

	Cptr = 1
	Dptr = 1
	call calc_bca
	Cptr2 = Cptr
	Dptr2 = Dptr

	p1(1) = dcmplx(EE(1))
	p1(4) = dcmplx(-Pin)
	p2(1) = dcmplx(EE(2))
	p2(4) = dcmplx(Pin)
	k1(1) = dcmplx(EE(3))
	k2(1) = dcmplx(EE(4))

	if(MASS1 .ne. 0D0) then
	  e1(0, 1) = dcmplx(-Pin/MASS1)
	  e1(0, 4) = dcmplx(EE(1)/MASS1)
	endif
	if(MASS2 .ne. 0D0) then
	  e2(0, 1) = dcmplx(-Pin/MASS2)
	  e2(0, 4) = dcmplx(-EE(2)/MASS2)
	endif

	print *, "sqrt(s) =", Ecms
	end

************************************************************************
	subroutine settheta
	implicit none
#include "kin.h"

#ifdef DYSON
	double complex shh
	external shh
#endif

	st = sin(th)
	ct = cos(th)
	T = mass(1) + mass(3) - 2D0*(EE(1)*EE(3) - Pin*Pout*ct)
	U = mass(2) + mass(4) - 2D0*(EE(2)*EE(4) + Pin*Pout*ct)

#ifdef DYSON
	resoT = 1D0/(T - MH2 + shh(T))
	resoU = 1D0/(U - MH2 + shh(U))
#else
	resoT = 1D0/(T - MH2)
	resoU = 1D0/(U - MH2)
#endif

	call calc_bca2

	k1(2) = dcmplx(-Pout*st)
	k1(4) = dcmplx(-Pout*ct)
	k2(2) = dcmplx(Pout*st)
	k2(4) = dcmplx(Pout*ct)

	if(MASS3 .ne. 0D0) then
	  e3(0, 1) = dcmplx(-Pout/MASS3)
	  e3(0, 2) = dcmplx(EE(3)/MASS3*st)
	  e3(0, 4) = dcmplx(EE(3)/MASS3*ct)
	  e3(1, 2) = dcmplx(-ds2*ct)
	  e3(1, 4) = dcmplx(ds2*st)
	  e3(2, 2) = dcmplx(-ds2*ct)
	  e3(2, 4) = dcmplx(ds2*st)
	endif
	if(MASS4 .ne. 0D0) then
	  e4(0, 1) = dcmplx(-Pout/MASS4)
	  e4(0, 2) = dcmplx(-EE(4)/MASS4*st)
	  e4(0, 4) = dcmplx(-EE(4)/MASS4*ct)
	  e4(1, 2) = dcmplx(ds2*ct)
	  e4(1, 4) = dcmplx(-ds2*st)
	  e4(2, 2) = dcmplx(ds2*ct)
	  e4(2, 4) = dcmplx(-ds2*st)
	endif
	end

************************************************************************
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
	double complex function Pair(i, j)
	implicit none
	integer i, j
#include "kin.h"

	Pair = v4(1, i)*v4(1, j) -
     +    v4(2, i)*v4(2, j) -
     +    v4(3, i)*v4(3, j) -
     +    v4(4, i)*v4(4, j)
	end

************************************************************************
	double complex function Eps(i, j, k, l)
	implicit none
	integer i, j, k, l
#include "kin.h"

	integer r, ev(4, 12), od(4, 12)

	data ev /
     +    1,2,3,4,  1,3,4,2,
     +    1,4,2,3,  2,1,4,3,
     +    2,3,1,4,  2,4,3,1, 
     +    3,1,2,4,  3,2,4,1,
     +    3,4,1,2,  4,1,3,2,
     +    4,2,1,3,  4,3,2,1 /
	data od /
     +    1,2,4,3,  1,3,2,4,
     +    1,4,3,2,  2,1,3,4,
     +    2,3,4,1,  2,4,1,3,
     +    3,1,4,2,  3,2,1,4,
     +    3,4,2,1,  4,1,2,3,
     +    4,2,3,1,  4,3,1,2 /

	Eps = dcmplx(0D0, 0D0)
	do r = 1, 12
	  Eps = Eps +
     +      v4(ev(1, r), i)*v4(ev(2, r), j)*
     +        v4(ev(3, r), k)*v4(ev(4, r), l) -
     +      v4(od(1, r), i)*v4(od(2, r), j)*
     +        v4(od(3, r), k)*v4(od(4, r), l)
	enddo
	Eps = Eps*dcmplx(0D0, -1D0)
	end

************************************************************************
#ifdef BREMSSTRAHLUNG
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

************************************************************************
	double precision function norm(z)
	implicit none
	double complex z

	norm = dble(z)**2 + dimag(z)**2
	end

************************************************************************
	subroutine getamp
	implicit none
#include "kin.h"
#include "rcsm.h"

#if FERMIONIC
#  define RESTYPE double precision
#else
#  define RESTYPE double complex
	double precision norm
	external norm
#  ifdef BORN
	double complex BORN
	external BORN
#  endif
#endif

	RESTYPE Mborn, Mloop
	RESTYPE born
	external born

#ifndef oneloop
	RESTYPE oneloop
	external oneloop
#endif

#ifdef BREMSSTRAHLUNG
	double precision dSB
#  ifndef softphot
	double precision softphot
	external softphot
#  endif
#endif

	call settheta
#ifdef BREMSSTRAHLUNG
	dSB = softphot()
#  ifdef DEBUG
	print *, "dSB =", dSB
#  endif
#endif
	bornamp = 0D0
	loopamp = 0D0
	do i4 = bpol(4), epol(4)
	  do i3 = bpol(3), epol(3)
	    do i2 = bpol(2), epol(2)
	      do i1 = bpol(1), epol(1)
	        call calc_sme
	        call calc_abbr
#if FERMIONIC
		call calc_mat
		Mborn = born()
		bornamp = bornamp + Mborn
		Mloop = oneloop()
#  ifdef BREMSSTRAHLUNG
	        Mborn = Mborn*dSB
#  endif
		loopamp = loopamp + Mborn + Mloop
#else
	        Mborn = born()
	        bornamp = bornamp + norm(Mborn)
	        Mloop = oneloop()
#  ifdef BREMSSTRAHLUNG
	        Mloop = Mloop + Mborn*dSB
#  endif
	        loopamp = loopamp + 2D0*dble(dconjg(Mborn)*Mloop)
#  ifdef BORN
	        Mborn = BORN()
#  endif
	        loopamp = loopamp + norm(Mborn)
#endif

#ifdef DEBUG
	        print *, "  Mborn =", Mborn
	        print *, "  Mloop =", Mloop
#endif
	      enddo
	    enddo
	  enddo
	enddo
	end

************************************************************************
	subroutine intgauss
	implicit none
#include "tools.h"
#include "kin.h"

	double precision gauss_x(GAUSSPOINTS/2)
	double precision gauss_w(GAUSSPOINTS/2)

	integer i
	double precision x, w
	double precision bornsumF, bornsumB, loopsumF, loopsumB
	double precision lower, upper, mid1, mid2, range
	parameter ( lower = CUTOFF, upper = pi - CUTOFF,
     +    mid1 = .5D0*(lower + .5D0*pi),
     +    mid2 = .5D0*(upper + .5D0*pi),
     +    range = .25D0*(upper - lower) )

#include "gauss.F"

* Integration is done in two intervals separately:
* [CUTOFF, pi/2] and [pi/2, pi - CUTOFF].
* In this way, we get the F-B asymmetry for free, even
* though the number of sampling points is the same,
* e.g. instead of 32 over the whole interval we can do
* with 16 over each half.

	call setenergy
	bornsumF = 0D0
	loopsumF = 0D0
	bornsumB = 0D0
	loopsumB = 0D0
	do i = 1, GAUSSPOINTS/2
	  x = gauss_x(i)*range
	  w = gauss_w(i)*range*2D0*pi*IDENTICALFACTOR
*			*2*pi from phi integration
*			possibly /2 because of identical particles
	  th = mid1 + x
	  call getamp
	  bornsumF = bornsumF + st*w*bornamp
	  loopsumF = loopsumF + st*w*loopamp
	  th = mid1 - x
	  call getamp
	  bornsumF = bornsumF + st*w*bornamp
	  loopsumF = loopsumF + st*w*loopamp
	  th = mid2 + x
	  call getamp
	  bornsumB = bornsumB + st*w*bornamp
	  loopsumB = loopsumB + st*w*loopamp
	  th = mid2 - x
	  call getamp
	  bornsumB = bornsumB + st*w*bornamp
	  loopsumB = loopsumB + st*w*loopamp
	  Cptr = Cptr2
	  Dptr = Dptr2
	enddo
	bornsum = bornsumF + bornsumB
	loopsum = loopsumF + loopsumB
	if(bornsum .eq. 0D0) then
	  bornAfb = 999D0
	else
	  bornAfb = (bornsumF - bornsumB)/bornsum
	endif
	if(loopsum .eq. 0D0) then
	  loopAfb = 999D0
	else
	  loopAfb = (loopsumF - loopsumB)/loopsum
	endif
	end
