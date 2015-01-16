*	num.F
*	this is the master file for numerics
*	for a four-gauge-boson process
*	apr96 th; last modified jul98

*#define DEBUG

*#define BREMSSTRAHLUNG
#define MBOS MZ
*--> the external mass
#define ALLOW_LONGITUDINAL

#define ESOFT EE
*#define ESOFT (.1D0*EE)
*#define ESOFT MBOS
*--> soft photon cutoff

#ifdef DEBUG
#define STEP (40D0*degree)
#else
#define STEP (1D0*degree)
#endif

#define CUTOFF (10d0*degree)

#define MINSTEP 5.00653422312312D0
#define STEPDIV 2.00045203475398759D-6
*--> some silly number to avoid S===MH2

#define GAUSSPOINTS 32

*			   1-4: p1_,p2_,k1_,k2_
#define E1START 5
*			   5-7: e1_(0, +, -)
#define E2START 8
*			  8-10: e2_(0, +, -)
#define E3START 11
*			 11-13: e3_(0, +, -)
#define E4START 14
*			 14-16: e4_(0, +, -)

#define p1(i) v(i,1)
#define p2(i) v(i,2)
#define k1(i) v(i,3)
#define k2(i) v(i,4)
#define e1(pol,i) v(i,E1START + pol)
#define e2(pol,i) v(i,E2START + pol)
#define e3(pol,i) v(i,E3START + pol)
#define e4(pol,i) v(i,E4START + pol)

#include "feyn.F"
#include "tools.F"

**************************************************************
** MAIN PROGRAM STARTS HERE:

	program num
#include "kin.h"
#include "tools.h"

	integer*4 iargc, argc
	external iargc, getarg
	character*40 argv, outf
	double precision Ecmsto

	argc = iargc()
	if(argc .lt. 5 .or. argc .gt. 6) then
	  print *, "USAGE:"
	  print *, "to calculate differential cross-section:"
	  print *, "	num i1 i2 i3 i4 ecms"
	  print *, "to calculate total cross-section:"
	  print *, "	num i1 i2 i3 i4 ecms_from ecms_to"
	  print *, "for polarization i1+i2 -> i3+i4"
	  print *, "where i[1-4] may be one of"
	  print *, "  0 or l   for longitudinal   polarization"
	  print *, "    +      for right-circular      ''"
	  print *, "    -      for left-circular       ''"
	  print *, "    t      for transverse          ''"
	  print *, "    u      for unpolarized boson"
	  stop
	endif

	avgfac = 1D0
	call getpol(1)
	call getpol(2)
	call getpol(3)
	call getpol(4)

22	format(E20.0)
23	format("diff.MH=", I4.4, ".pol=", A4, ".E=", I5.5)
24	format("tot.MH=", I4.4, ".pol=", A4, ".E=", I5.5, "-", I5.5)

	call getarg(5, argv)
	read(argv, 22) Ecms
	if(argc .eq. 5) then
	  write(outf, 23) int(MH), polstr, int(Ecms)
	else
	  call getarg(6, argv)
	  read(argv, 22) Ecmsto
	  write(outf, 24) int(MH), polstr, int(Ecms), int(Ecmsto)
	endif
	print *, "output file: ", outf
	print *, "angular cut: ", CUTOFF

	call setup
	open(11, file=outf, status='UNKNOWN')

26	format(3G24.15)

	if(argc .eq. 5) then
	  call setenergy
	  do 17 th = CUTOFF, pi-CUTOFF, STEP
	    print *, "th =", th
	    call getamp
#ifdef oneloop
	    feynamp = 999D0
#else
	    feynamp = feynamp*kinf
#endif
	    write(11, 26) th, bornamp*kinf, feynamp
	    close(11)
	    open(11, file=outf, access='APPEND', status='OLD')
	    Cptr = Cptr2
	    Dptr = Dptr2
17	  continue
	else
18	  continue
	  call intgauss
#ifdef oneloop
	  feynsum = 999D0
#else
	  feynsum = feynsum*kinf
#endif
	  write(11, 26) Ecms, bornsum*kinf, feynsum
	  close(11)
	  open(11, file=outf, access='APPEND', status='OLD')
	  if(Ecms .lt. Ecmsto) then
	    Ecms = Ecms + max(MINSTEP, dble(S - MH2)*STEPDIV)
	    goto 18
	  endif
	endif

	close(11)
	call ffexi
	end

**************************************************************
	subroutine getpol(i)
	implicit none
	integer i
#include "kin.h"

	character c*1, argv*40
	integer j

	integer offset(4)
	data offset / E1START, E2START, E3START, E4START /

	call getarg(i, argv)
	c = argv(1:1)
	if(c .eq. 'u') then
#ifndef ALLOW_LONGITUDINAL
	  stop "No longitudinal polarizations allowed"
#else
	  bpol(i) = 0
	  epol(i) = 2
	  if(i .lt. 3) avgfac = avgfac*1D0/3D0
	  c = 'U'
#endif
	else if(c .eq. 't') then
	  bpol(i) = 1
	  epol(i) = 2
	  if(i .lt. 3) avgfac = avgfac*1D0/2D0
	  c = 'T'
	else
	  if(c .eq. 'l' .or. c .eq. '0') then 
#ifndef ALLOW_LONGITUDINAL
	    stop "No longitudinal polarizations allowed"
#else
	    j = 0
	    c = 'L'
#endif
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
	polstr(i:i) = c
	bpol(i) = bpol(i) + offset(i)
	epol(i) = epol(i) + offset(i)
	end

**************************************************************
	subroutine setup
	implicit none
#include "kin.h"

	double complex ds2I, ds2mI, zero
	parameter (ds2I = dcmplx(0D0, ds2),
     +    ds2mI = dcmplx(0D0, -ds2),
     +    zero = dcmplx(0d0, 0d0))

	integer i, j

	do 7 i = 1, 4
	  do 7 j = 1, 16 
	    v(i, j) = zero
7	continue

	e1(1, 2) = -ds2
	e1(1, 3) = ds2I
	e1(2, 2) = -ds2
	e1(2, 3) = ds2mI

	e2(1,2) = ds2
	e2(1,3) = ds2I
	e2(2,2) = ds2
	e2(2,3) = ds2mI

	e3(1,3) = ds2mI
	e3(2,3) = ds2I

	e4(1,3) = ds2mI
	e4(2,3) = ds2I

	MH2 = dble(MH)**2
	call bcaini
	call calc_renconst
	end

**************************************************************
	subroutine setenergy
	implicit none
#include "tools.h"
#include "kin.h"
#include "rcsm.h"

#ifdef DYSON
	double complex shh
	external shh
#endif DYSON

#ifdef BREMSSTRAHLUNG
	double precision Iij
	external Iij

	double precision delta
	common /ffcut/ delta

	double precision logEmin, brS
	common /brems/ logEmin, brS
#endif BREMSSTRAHLUNG

	EE = .5D0*Ecms
	S = Ecms*Ecms
	kinf = KINFAC*avgfac

	psq = EE*EE - m2
	p = sqrt(psq)

#ifdef ALLOW_LONGITUDINAL
	EEs = EE/MBOS
	ps = p/MBOS
#else
	EEs = 999D0
	ps = 999D0
#endif

#ifdef DYSON
	reso = 1D0/(S - MH2 + shh(S))
#else
	reso = 1D0/(S - MH2)
	resoc = 1D0/(S - MH2 + dcmplx(0D0, 1D0)*GammaHMH)
#endif DYSON

#ifdef BREMSSTRAHLUNG
	logEmin = log(4D0*ESOFT**2/delta)
	brS = logEmin + EE/p*log((EE - p)/(EE + p)) -
     +    (S - 2D0*m2)*Iij(S + sqrt(S*(S - 4D0*m2)))
#endif BREMSSTRAHLUNG

	Cptr = 1
	Dptr = 1
	call calc_bca
	Cptr2 = Cptr
	Dptr2 = Dptr

	p1(1) = dcmplx(EE)
	p1(4) = dcmplx(-p)
	p2(1) = dcmplx(EE)
	p2(4) = dcmplx(p)
	k1(1) = dcmplx(EE)
	k2(1) = dcmplx(EE)

	e1(0,1) = dcmplx(-ps)
	e1(0,4) = dcmplx(EEs)
	e2(0,1) = dcmplx(-ps)
	e2(0,4) = dcmplx(-EEs)
	e3(0,1) = dcmplx(ps)
	e4(0,1) = dcmplx(ps)

	print *, "Ecms =", Ecms
	end

**************************************************************
	subroutine settheta
	implicit none
#include "kin.h"

#ifdef DYSON
	double complex shh
	external shh
#endif DYSON

	st = sin(th)
	ct = cos(th)
	T = -2D0*psq*(1D0 - ct)
	U = -2D0*psq*(1D0 + ct)

#ifdef DYSON
	resoT = 1D0/(T - MH2 + shh(T))
	resoU = 1D0/(U - MH2 + shh(U))
#else
	resoT = 1D0/(T - MH2)
	resoU = 1D0/(U - MH2)
#endif DYSON

	call calc_bca2

	k1(2) = dcmplx(-p*st)
	k1(4) = dcmplx(-p*ct)
	k2(2) = dcmplx(p*st)
	k2(4) = dcmplx(p*ct)

	e3(0, 2) = dcmplx(-EEs*st)
	e3(0, 4) = dcmplx(-EEs*ct)
	e4(0, 2) = dcmplx(EEs*st)
	e4(0, 4) = dcmplx(EEs*ct)

	e3(1, 2) = dcmplx(-ds2*ct)
	e3(1, 4) = dcmplx(ds2*st)
	e4(1, 2) = dcmplx(ds2*ct)
	e4(1, 4) = dcmplx(-ds2*st)

	e3(2, 2) = dcmplx(-ds2*ct)
	e3(2, 4) = dcmplx(ds2*st)
	e4(2, 2) = dcmplx(ds2*ct)
	e4(2, 4) = dcmplx(-ds2*st)
	end

**************************************************************
	double complex function Pair(i, j)
	implicit none
	integer i, j
#include "kin.h"

	Pair = v(1,i)*v(1,j) -
     +    v(2,i)*v(2,j) -
     +    v(3,i)*v(3,j) -
     +    v(4,i)*v(4,j)
	end

**************************************************************
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
	do 10 r = 1, 12
	  Eps = Eps +
     +     v(ev(1, r), i)*v(ev(2, r), j)*
     +       v(ev(3, r), k)*v(ev(4, r), l) -
     +     v(od(1, r), i)*v(od(2, r), j)*
     +       v(od(3, r), k)*v(od(4, r), l)
10	continue
	Eps = Eps*dcmplx(0D0, -1D0)
	end

#ifdef BREMSSTRAHLUNG
**************************************************************
	double precision function Iij(x)
	implicit none
	double precision x
#include "kin.h"

	double precision logEmin, brS
	common /brems/ logEmin, brS

	double precision xx, h1, h2
	double complex tmp, li2a, li2b, li2c, li2d
	integer ier, ipa, ipb, ipc, ipd

	double precision pi12
	parameter (pi12 = pi**2/12D0)

	xx = x - 2D0*m2
	h1 = 2D0*Ecms*(EE - p)/x
	h2 = 2D0*Ecms*(EE + p)/x
	ier = 0
	call ffzxdl(li2a, ipa, tmp, dble(1D0 - h1*xx/2D0/m2), 1, ier)
	call ffzxdl(li2b, ipb, tmp, dble(1D0 - h2*xx/2D0/m2), 1, ier)
	call ffzxdl(li2c, ipc, tmp, dble(1D0 - h1), 1, ier)
	call ffzxdl(li2d, ipd, tmp, dble(1D0 - h2), 1, ier)

	Iij = 2D0*xx/(x*(x - 4D0*m2))*(logEmin*log(xx/2D0/m2) +
     +    dble(li2a) + dble(li2b) - dble(li2c) - dble(li2d) +
     +    (ipa + ipb - ipc - ipd)*pi12)
	end

**************************************************************
	double precision function softphot()
	implicit none
#include "kin.h"

	double precision Iij
	external Iij

	double precision logEmin, brS
	common /brems/ logEmin, brS

	softphot = -Alpha/pi*(brS -
     +    (2D0*m2 - T)*Iij(4D0*m2 - T + sqrt(T*(T - 4D0*m2))) +
     +    (2D0*m2 - U)*Iij(4D0*m2 - U + sqrt(U*(U - 4D0*m2))) )
	end
#endif BREMSSTRAHLUNG

**************************************************************
	subroutine getamp
	implicit none
#include "kin.h"
#include "rcsm.h"

	double complex bo, fe
	double precision dSB
	double complex born
	external born

#ifndef oneloop
	double complex oneloop
	external oneloop
#endif

#ifdef BORN
	double complex BORN
	external BORN
#endif

#ifdef BREMSSTRAHLUNG
	double precision softphot
	external softphot
#endif

	double precision norm
	double complex z
	norm(z) = dble(z)**2 + dimag(z)**2

	call settheta
#ifdef BREMSSTRAHLUNG
	dSB = softphot()
#  ifdef DEBUG
	print *, "dSB =", dSB
#  endif
#endif
	bornamp = 0D0
	feynamp = 0D0
	do 60 i4 = bpol(4), epol(4)
	  do 60 i3 = bpol(3), epol(3)
	    do 60 i2 = bpol(2), epol(2)
	      do 60 i1 = bpol(1), epol(1)
	        call calc_sme
	        call calc_abbr
	        bo = born()
	        fe = oneloop()
#ifdef DEBUG
	        print *, "  born             =", bo
#  ifdef BORN
	        print *, "  BORN             =", BORN()
#  endif BORN
	        print *, "  oneloop          =", fe
	        print *, "  fe+bo*dWFW       =", fe + bo*dWFW
#  ifdef BREMSSTRAHLUNG
	        print *, "  fe+bo*(dWFW+dSB) =", fe + bo*(dWFW + dSB)
#  endif BREMSSTRAHLUNG
#endif DEBUG
	        bornamp = bornamp + norm(bo)
#ifdef BREMSSTRAHLUNG
	        fe = fe + bo*(dWFW + dSB)
	        feynamp = feynamp + 2D0*dble(dconjg(bo)*fe)
#  ifdef BORN
	        bo = BORN()
#  endif BORN
	        feynamp = feynamp + norm(bo)
#else
* this is for ZZ -> ZZ
#  ifdef DYSON
	        feynamp = feynamp + norm((bo + fe)*(1D0 + dWFZ))
#  else
	        fe = fe + bo*dWFZ
#    ifdef BORN
	        bo = BORN()
#    endif BORN
	        feynamp = feynamp + norm(bo + fe)
#  endif DYSON
#endif BREMSSTRAHLUNG
60	continue
	end

**************************************************************
	subroutine intgauss
	implicit none
#include "tools.h"
#include "kin.h"

	double precision gauss_x(GAUSSPOINTS/2)
	double precision gauss_w(GAUSSPOINTS/2)

	integer i
	double precision w
	double precision lower, upper, mid, range
	parameter ( lower = CUTOFF, upper = pi - CUTOFF,
     +    mid = .5D0*(upper + lower),
     +    range = .5D0*(upper - lower) )

#include "gauss.F"

	call setenergy
	bornsum = 0D0
	feynsum = 0D0
	do 150 i = 1, GAUSSPOINTS/2
#ifdef BREMSSTRAHLUNG
	  w = gauss_w(i)*range*pi*2D0
#else
	  w = gauss_w(i)*range*pi
#endif
*			*2*pi from phi integration
*			/2 because of identical particles (ZZ -> ZZ)
	  th = mid + gauss_x(i)*range
	  call getamp
	  bornsum = bornsum + st*w*bornamp
	  feynsum = feynsum + st*w*feynamp
	  th = 2D0*mid - th
	  call getamp
	  bornsum = bornsum + st*w*bornamp
	  feynsum = feynsum + st*w*feynamp
	  Cptr = Cptr2
	  Dptr = Dptr2
150	continue
	end
