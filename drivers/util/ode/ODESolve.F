* ODESolve.F
* ODE solver a la Numerical Recipes
* this file is part of FormCalc
* last modified 10 Apr 19 th

#include "types.h"

c#define VERBOSE
c#define WARNINGS


	subroutine ODESolve(nvar, f, ystart, x1, x2, eps, h1, userarg)
	implicit none
	integer nvar
	external f, userarg
	RealType ystart(nvar), x1, x2, eps, h1

	RealType safety, pgrow, pshrink, errcon
	parameter (safety = .9D0)
	parameter (pgrow = -.2D0)
	parameter (pshrink = -.25D0)
	parameter (errcon = 1.89D-4)

	RealType A2, B21
	RealType A3, B31, B32
	RealType A4, B41, B42, B43
	RealType A5, B51, B52, B53, B54
	RealType A6, B61, B62, B63, B64, B65
	RealType C1, C3, C4, C6, DC1, DC3, DC4, DC5, DC6
	parameter (A2 = .2D0, B21 = .2D0)
	parameter (A3 = .3D0, B31 = 3/40D0, B32 = 9/40D0)
	parameter (A4 = .6D0, B41 = .3D0, B42 = -.9D0, B43 = 1.2D0)
	parameter (A5 = 1, B51 = -11/54D0, B52 = 2.5D0, B53 = -70/27D0,
     &    B54 = 35/27D0)
	parameter (A6 = .875D0, B61 = 1631/55296D0, B62 = 175/512D0,
     &    B63 = 575/13824D0, B64 = 44275/110592D0, B65 = 253/4096D0)
	parameter (C1 = 37/378D0, DC1 = C1 - 2825/27648D0)
	parameter (C3 = 250/621D0, DC3 = C3 - 18575/48384D0)
	parameter (C4 = 125/594D0, DC4 = C4 - 13525/55296D0)
	parameter (DC5 = -277/14336D0)
	parameter (C6 = 512/1771D0, DC6 = C6 - .25D0)

	integer maxstep
	parameter (maxstep = 10000)

	RealType x0, yy(nvar,-1:6)
	RealType h, h0, yerr, yscal, errmax
	integer i, step

#define y0	yy(1,-1)
#define Y0	yy(1:nvar,-1)
#define y	yy(1,0)
#define Y	yy(1:nvar,0)
#define dy(i)	yy(1,i)
#define DY(i)	yy(1:nvar,i)

	h = sign(h1, x2 - x1)
	x0 = x1
	Y0 = ystart

	do step = 1, maxstep
#ifdef VERBOSE
	  print *, "step,stepsize,x ", step, h, x0
	  print *, "y =", Y0
#endif
	  call f(dy(1), y0, x0, userarg)
	  h0 = h
	  if( (x0 + h - x2)*(x0 + h - x1) > 0 ) h = x2 - x0

* determine Cash-Karp step
1	  Y = Y0 + h*B21*DY(1)
	  call f(dy(2), y, x0 + A2*h, userarg)

	  Y = Y0 + h*(B31*DY(1) + B32*DY(2))
	  call f(dy(3), y, x0 + A3*h, userarg)

	  Y = Y0 + h*(B41*DY(1) + B42*DY(2) + B43*DY(3))
	  call f(dy(4), y, x0 + A4*h, userarg)

	  Y = Y0 + h*(B51*DY(1) + B52*DY(2) + B53*DY(3) + B54*DY(4))
	  call f(dy(5), y, x0 + A5*h, userarg)

	  Y = Y0 + h*(B61*DY(1) + B62*DY(2) + B63*DY(3) + B64*DY(4) + B65*DY(5))
	  call f(dy(6), y, x0 + A6*h, userarg)

	  errmax = maxval(
     &      abs(h*(DC1*DY(1) + DC3*DY(3) + DC4*DY(4) + DC5*DY(5) + DC6*DY(6)))/
     &      (abs(Y0) + abs(h0*DY(1)) + 1D-30)
     &    )/eps
	  if( errmax > 1 ) then
	    h = sign(max(abs(safety*h*errmax**pshrink), .1D0*abs(h)), h)
	    if( x0 + h /= x0 ) goto 1
#ifdef WARNINGS
	    print *, "stepsize underflow in ODESolve"
#endif
	  endif

	  Y0 = Y0 + h*(C1*DY(1) + C3*DY(3) + C4*DY(4) + C6*DY(6))
	  x0 = x0 + h
	  if( (x0 - x2)*(x2 - x1) >= 0 ) goto 9

	  if( errmax > errcon ) then
	    h = safety*h*errmax**pgrow
	  else
	    h = 5*h
	  endif
	enddo

#ifdef WARNINGS
	print *, "RGE solution did not converge"
#endif

9	ystart = Y0
	end

