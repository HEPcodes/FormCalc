* hardradiation.F
* routines for the treatment of hard- and soft-photon radiation
* this file is part of FormCalc
* Christian Schappacher, 24 Oct 13
* last modified 22 Jul 18 th


************************************************************************
** HardRadiationFactor computes the contributions from fully inclusive
** hard- and soft-photon emission off three different legs.
** Keep the following order for the (generic) particles: 
** SVV, SSV, SFF, FFS, FFV and VFF. 
** Massless particles (neutrinos) must always be the second particle.
** SSS, SFF and SSV are based loosely on code from SFOLD, arXiv:1104.2151.
** FFS and FFV are based on code by S. Liebler.
** The integrals are from A. Denner, Forts. Phys. 41 (1993) 307 [arXiv:0709.1075].

	subroutine HardRadiationFactors(chg)
	implicit none
	RealType chg(LEGS)

#include "xsection.h"

	integer i, j
	RealType m1, m2, m3, m12, m22, m32
	RealType d12, d13, d23, X, Y, Z
	RealType g1, g2, g3, g11, g12, g13, g22, g23, g33
	RealType k, logk, beta, fac, t
	RealType x1, x2, x3, logb1, logb2, logb3
	RealType l0, l2, li2_1, li2_2, li2_3
	RealType li2_12, li2_13, li2_23
	RealType I0, I1, I2, I00, I01, I02, I11, I12, I22
	RealType II, I0_1, I0_2, I1_0, I1_2, I2_0, I2_1
	RealType I00_11, I00_12, I00_22
	RealType I11_00, I11_02, I11_22
	RealType I22_00, I22_01, I22_11
	RealType G1m
	ComplexType cR, cL

#define VERT(a,b,c) Compose(ARG_Ptyp,JOIN_OCT,a,b,c)

	integer SSS, SSV, SVV, SFF, FFS, FFV, VFF
	parameter (SSS = VERT(SCALAR,SCALAR,SCALAR))
	parameter (SSV = VERT(SCALAR,SCALAR,VECTOR))
	parameter (SVV = VERT(SCALAR,VECTOR,VECTOR))
	parameter (SFF = VERT(SCALAR,FERMION,FERMION))
	parameter (FFS = VERT(FERMION,FERMION,SCALAR))
	parameter (FFV = VERT(FERMION,FERMION,VECTOR))
	parameter (VFF = VERT(VECTOR,FERMION,FERMION))

#define zero(m) (m .lt. 1D-6)
#define nonzero(m) (m .gt. 1D-6)

	g1 = chg(1)
	g11 = g1**2
	g2 = chg(2)
	g22 = g2**2
	g12 = g1*x2
	g3 = chg(3)
	g33 = g3**2
	g13 = g1*g3
	g23 = g1*g3

	m1 = mass(1)
	m12 = m1**2
	m2 = mass(2)
	m22 = m2**2
	m3 = mass(3)
	m32 = m3**2
	d12 = (m1 - m2)*(m1 + m2)
	d13 = (m1 - m3)*(m1 + m3)
	d23 = (m2 - m3)*(m2 + m3)
	X = d12 - m32
	Y = d12 + m32
	Z = d23 + m12

	k = sqrt(X**2 - 4*m22*m32)

	if( nonzero(m2) .and. nonzero(m3) ) then
	  logk = 2*log(k) - log(m12*m22*m32*getlambda())/2

	  x1 = (X + k)/(2*m2*m3)
	  logb1 = log(x1)
	  x2 = (Y - k)/(2*m1*m3)
	  logb2 = log(x2)
	  x3 = (Z - k)/(2*m1*m2)
	  logb3 = log(x3)

	  fac = .25D0/m12

	  t = k*(logk - 1)
	  I00 = fac*(t - d23*(logb2 - logb3) - logb1*m12)/m12
	  I11 = fac*(t - d13*(logb1 - logb3) - logb2*m22)/m22
	  I22 = fac*(t - d12*(logb1 - logb2) - logb3*m32)/m32

	  x1 = Re(Li2omx(x1**2))
	  x2 = Re(Li2omx(x2**2))
	  x3 = Re(Li2omx(x3**2))
	  li2_12 = x1 - x2 + (logb1 - logb2)*(logb1 + logb2)
	  li2_13 = x1 - x3 + (logb1 - logb3)*(logb1 + logb3)
	  li2_23 = x2 - x3 + (logb2 - logb3)*(logb2 + logb3)

	  I01 = fac*(2*logk*logb3 - li2_13 - li2_23)
	  I02 = fac*(2*logk*logb2 - li2_12 + li2_23)
	  I12 = fac*(2*logk*logb1 + li2_12 + li2_13)

	  II = fac*(k*(m12 + m22 + m32)/2 + 
     &      2*(logb1*m22*m32 + m12*(logb3*m22 + logb2*m32)))

	  I0 = -fac*(k + 2*(logb3*m22 + logb2*m32))
	  I1 = -fac*(k + 2*(logb3*m12 + logb1*m32))
	  I2 = -fac*(k + 2*(logb2*m12 + logb1*m22))

	  I0_1 = fac*(logb3*m22**2 - logb2*m32*(d12 + Y) +
     &      k*(2*(d23 - m32) - Y)/4)
	  I0_2 = fac*(logb2*m32**2 - logb3*m22*(d13 + Z) -
     &      k*(2*(d23 + m22) + Z)/4)
	  I1_0 = fac*(logb3*m12**2 + logb1*m32*(d12 + X) +
     &      k*(2*(d13 - m32) + X)/4)
	  I1_2 = fac*(logb1*m32**2 - logb3*m12*(d23 + Z) -
     &      k*(2*(d13 + m12) + Z)/4)
	  I2_0 = fac*(logb2*m12**2 + logb1*m22*(d13 + X) +
     &      k*(2*(d12 - m22) + X)/4)
	  I2_1 = fac*(logb1*m22**2 + logb2*m12*(d23 - Y) - 
     &      k*(2*(d12 + m12) + Y)/4)

	  t = k**3/(6*m12)
	  I00_12 = -fac*(t + logb3*m22**2 + logb2*m32**2 +
     &      k*(.5D0*m12 - .75D0*X))
	  I00_11 = fac*(t + 2*m32*(k + logb2*Y))
	  I00_22 = fac*(t + 2*m22*(k + logb3*Z))

	  t = k**3/(6*m22)
	  I11_02 = -fac*(t + logb3*m12**2 + logb1*m32**2 +
     &      k*(.5D0*m22 + .75D0*Y))
	  I11_00 = fac*(t + 2*m32*(k - logb1*X))
	  I11_22 = fac*(t + 2*m12*(k + logb3*Z))

	  t = k**3/(6*m32)
	  I22_01 = -fac*(t + logb2*m12**2 + logb1*m22**2 +
     &      k*(.5D0*m32 + .75D0*Z))
	  I22_00 = fac*(t + 2*m22*(k - logb1*X))
	  I22_11 = fac*(t + 2*m12*(k + logb2*Y))

* If particle two is massless (e.g. neutrino):
	else if( zero(m2) .and. nonzero(m3) ) then
	  logk = 2*log(X) - log(m12*m32*getlambda())/2

	  logb2 = log(m3/m1)

	  fac = .25D0/m12

	  I0 = -fac*(2*logb2*m32 + X)
	  I2 = -fac*(2*logb2*m12 + X)

	  II = fac*(m32*X + m12*(4*logb2*m32 + X))/2

	  I0_1 = -fac*(logb2*m32*(m32 + 2*m12) + (m12 + 5*m32)*X/4)
	  I2_1 = -fac*(logb2*m12*(m12 + 2*m32) + (m32 + 5*m12)*X/4)

	  I0_2 = fac*(logb2*m32**2 + (2*m32 - d13)*X/4)
	  I2_0 = fac*(logb2*m12**2 + (2*m12 + d13)*X/4)

	  t = X*(logk - 1)
	  l0 = log(X/m1)
	  l2 = log(X/m3)
	  I00 = fac*(t - m12*l2 + m32*log(m3*X/m12))/m12
	  I22 = fac*(t + m32*l0 - m12*log(m1*X/m32))/m32

	  I02 = -fac*((l0 - l2)*(l0 + l2) -
     &      2*(logb2*(logb2 + logk) + Re(Li2omx((m22 + m32)/m12))))

	  I00_12 = -fac*(X**3/(6*m12) +
     &      logb2*m32**2 - (m12 - 3*m32)*X/4)

	else
	  stop "Don't know how to handle this mass configuration"
	endif


	if( type .eq. SSS ) then
	  fac = -4*Alfa/Pi*m12/k
	  radii = fac*(g11*I00*m12 + g22*I11*m22 +
     &      g33*I22*m32 - g23*(I1 + I2 - I12*X) +
     &      g13*(I0 + I2 + I02*Y) + g12*(I0 + I1 + I01*Z))
	  radij = bogus

	else if( type .eq. SVV .and.
     &	    mass(2) .eq. MW .and. mass(3) .eq. MW ) then
	  beta = sqrt(1 - 4*m32/m12)
	  fac = Alfa/(4*Pi)*m32/
     &      (beta*(beta**2*m12**2 + 12*m32**2))
	  radii = -fac*
     &      (I11_00 + I22_00 + 46*(I11_02 + I22_01) +
     &        13*(I11_22 + I22_11) + 100*(I1_0 + I2_0) +
     &        80*(I1_2 + I2_1) + 134*II +
     &        (-42*(I1_0 + I1_2 + I2_0 + I2_1 + 2*II)*m12 +
     &           16*((m22 - d12)**2 + 8*m22**2)*
     &            (I1 + I2 - d12*I12 + (I11 + I12 + I22)*m22))/m22)
	  radij = bogus

	else if( type .eq. SSV ) then
	  fac = Alfa/Pi*m12*m32/k**3
	  radii = fac*((II*
     &         (d12**2*(g13 - g23 - g33) + 
     &           (-6*g13 + 5*(g11 + g22) - 2*(g12 + g23) - g33)*
     &            m32**2))/m32**2 - 
     &      (g11*(I00_11 + 2*I00_12 + I00_22)*m12 + 
     &         g22*(4*d12**2*I11 + I11_00 + 2*I11_02 + I11_22)*
     &          m22 + II*(3*(g11*m12 + g22*m22) - 
     &            4*g12*(m12 + m22) + g33*(m12 + 3*m22) - 
     &            g13*(5*m12 + 3*m22) + g23*(m12 + 7*m22)) + 
     &         4*d12*(d12**2*(g13*I02 + g23*I12) + 
     &            d12*g12*I01*(m12 + m22) + 
     &            I0*(2*g11*m12 - (g12 + g13)*(m12 + m22)) + 
     &            I1*(-2*g22*m22 + g12*(m12 + m22) - 
     &               g23*(m12 + m22))) - 
     &         2*g33*(I22_00 - 2*I22_01 + I22_11)*m32 + 
     &         2*(d12**2*((3*g13 - 3*g23 - g33)*I2 + 
     &               2*g11*I00*m12) + 
     &            I2_1*(d12*(g13 - g23 - g33) + 4*(g13 + g23)*m32)) -
     &           (g12 - g22 - g23)*I1_2*X - 
     &         (g11 - g12 - g13)*I0_2*Y)/m32 - 
     &      4*g33*I22*(2*m12*X - Y*Z) + 
     &      4*(I0*(2*g11*m12 + (g12 + g13)*(d23 + m22)) + 
     &         I1*(2*g22*m22 + (g12 - g23)*(d13 + m12)) + 
     &         g12*I01*(2*m12*m22 + 3*(d13*m12 + d23*m22) +  m32**2) +
     &         g13*I02*(d12*m12 + d13*m32 + 3*m22*Y) + 
     &         (g11*I00*m12 + g22*I11*m22)*(m12 + m22 + Z)) + 
     &      2*(I2*(g33*(2*Z - m32) + 
     &            (g13 - g23)*(2*(m12 + m22) + m32)) + 
     &         g23*I12*(m32*(X + Y) + X*(2*Y + 4*Z))))
	  radij = bogus

	else if( type .eq. SFF ) then
	  fac = Alfa/(16*Pi**2)/m1
	  radii = 2*fac*(
     &      g22*I1_2 + g33*I2_1 - (g11 - g12 - g13 + 2*g23)*II + 
     &      I0*((g12 + g13)*(m12 + m22 + m32) - g11*(2*m12 + X)) -
     &      2*(I1*(d13*(g12 - g23) - g22*m22) + 
     &         I2*(d12*(g13 - g23) - g33*m32) + 
     &         X*(g11*I00*m12 + g22*I11*m22 + g33*I22*m32 + 
     &            g23*I12*X + g13*I02*Y + g12*I01*Z)) )
	  radij = 2*fac*m2*m3*(
     &      g11*(I0 + 2*I00*m12) + g22*(I1 + 2*I11*m22) + 
     &      g33*(I2 + 2*I22*m32) - g23*(I1 + I2 - 2*I12*X) + 
     &      g13*(I0 + I2 + 2*I02*Y) + g12*(I0 + I1 + 2*I01*Z) )

	else if( type .eq. FFS ) then
	  fac = Alfa/(32*Pi**2)/m1
	  radii = 2*fac*(
     &      g13*I0_2 - g12*(I0_1 + I1_0) - g23*I1_2 + 
     &      (g11 - 4*g12 - 2*g13 + g22 + 2*(g23 + g33))*II - 
     &      I0*(g11*Y - (g12 + g13)*(Y - 2*Z)) - 
     &      2*(I2*((g13 - g23)*m32 - g33*(m32 - Z)) + 
     &         (g11*I00*m12 + g22*I11*m22 + g33*I22*m32 + 
     &          g23*I12*X + g13*I02*Y)*Z + g12*I01*Z**2) + 
     &      I1*(g22*X - (g12 - g23)*(X + 2*Z)) )
	  radij = -16*fac*m1*m2*(
     &      g11*I00*m12 + g22*I11*m22 + g33*(I2 + I22*m32) - 
     &      g23*(I1 - I12*X) + g13*(I0 + I02*Y) + 
     &      g12*(I0 + I1 + I01*Z) )

	else if( type .eq. FFV ) then
	  fac = Alfa/(16*Pi**2)/(m1*m32)
	  radii = fac*(
     &      4*((2*g12 - g22)*I22_01 + (4*g12 - 3*g22)*I2_1 - 
     &          g11*(I22_01 + I2_1))*m32 + 
     &      (g11*I0_2 - g12*(I0_1 + I0_2) - g22*I1_0)*
     &        (m12 + m22 + 2*m32) - 
     &      2*(g11*(I0 + I00*m12) + g22*(I1 + I11*m22) + 
     &          (g11 - 2*g12 + g22)*(I2 + I22*m32) + 
     &          I01*(2*g12*m22 + g22*X) + 
     &          I02*(g22*X + g11*Y - g12*(X + Y)))*(X*Y + m32*Z) +
     &      II*(-4*g22*m32 + g11*(m12 + m22 + 2*m32) - 
     &            g12*(m12 + m22 + 2*Z)) )
	  radij = 2*fac*m1*m2*(
     &      g12*I0_1 - (g11 - g12)*I0_2 + g22*I1_0 - 
     &      (g11 - 3*g12)*II - 
     &      6*m32*(g12*(2*(d12*I02 + I2) - 2*I01*m22 + 2*I22*m32) -
     &        g22*(I1 + I2 + I11*m22 + I22*m32 + (I01 + I02)*X) -
     &        g11*(I0 + I2 + I00*m12 + I22*m32 + I02*Y)) )

	else if( type .eq. VFF ) then
	  fac = 4*Alfa/Pi*m12/k
	  G1m = X + Y*Z/m12
	  radii = -fac*
     &      (g23*I12*X - ((g23 + g33)*(2*I0_1 - G1m*I02*Y) + 
     &           (g22 + g23)*(2*I0_2 - G1m*I01*Z))/G1m + 
     &        (g22 + 2*g23 + g33)*(I0 + I00*m12 - 
     &          (2*I00_12*m12 - II*(m22 + m32)/2)/(G1m*m12)) + 
     &        g22*(I1 + I11*m22 + 
     &          (II*m12 + I1_0*(m12 + (m22 + m32)/2))/(G1m*m12)) + 
     &        g33*(I2 + I22*m32 + 
     &          (II*m12 + I2_0*(m12 + (m22 + m32)/2))/(G1m*m12)))
	  radij = radii

	else
	  stop "Don't know how to handle this vertex"
	endif
	end

