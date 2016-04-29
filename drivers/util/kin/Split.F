* Split.F
* splits a parent momentum into two child momenta
* loosely based on code by T. Fritzsche
* this file is part of FormCalc
* last modified 5 Apr 16 th

#include "types.h"
#include "util.h"


* The first child momentum is stored in vector i, the second
* replaces the input variables such that repeated invocation
* builds up the entire phase-space.

* on entry:
*   minv = available invariant mass
*   msum = sum of particle masses 1..n
*   p (ex,ey,ez) = parent momentum
*   xmrem = the fraction of invariant mass that goes into
*     child particles 2..n
*   xcosth and xphi, the fractions of costh and phi which
*     specify the direction relative to (ex,ey,ez) in the CMS
*     frame of the child particles

* on exit:
*   vector i is assigned child particle 1
*   minv = invariant mass remaining for child particles 2..n
*   msum = sum of particle masses 2..n
*   p (ex,ey,ez) = sum of momenta of child particles 2..n


* to resolve Var, Show, Lower, Upper, CutMin, CutMax:
#define var(i,v) v(i)

#define Bits(c,m) iand(key, Cut(c,m))
#define Test(c,m) Bits(c,m) .ne. 0
#define Zero(c,m) key = iand(key, not(Cut(c,m)))
#define kQ(m) ibits(key, m-1, 1)

#define CUT_ANY CUT_MIN + CUT_MAX
#define CUT_MREM_INDIRECT CUT_MREM_E + CUT_MREM_ET + CUT_MREM_RAP

	subroutine Split(i, m, p, ex, ey, ez, minv, msum,
     &    cutfac, cutkey, xmrem, xcosth, xphi)
	implicit none
	integer i, cutkey
	RealType m, p, ex, ey, ez, minv, msum
	RealType cutfac, xmrem(*), xcosth(*), xphi(*)

#include "util.h"

	RealType costh, sinth, cosphi, sinphi
	RealType bmrem(2), bcosth(2), cmrem(2)
	RealType fxz, fyz
	RealType fx, fy, fz
	RealType nx, ny, nz
	RealType fac, cut, c0, c1, c2, d, k2
	RealType r, r0, r0costh, rz, r1, k2max
	RealType m2, mrem, k, k10, q, qe, qn
	integer key, c, s, j, n

	RealType rpart(0:1)
	equivalence (rpart(0), r0), (rpart(1), r)

	RealType pi, eps
	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (eps = 1D-7)

	RealType sqdiff, k2eq, a_, b_
	sqdiff(a_, b_) = (a_ - b_)*(a_ + b_)
	k2eq(k2) = 2*c1*sqrt(k2*(k2 + m2)) + c2*k2 - c0

	key = cutkey
	fac = cutfac
	cutfac = 0

* rotate the unit vector e = (ex,ey,ez) into the desired direction:
* rotate by theta in the e-n-plane and then by phi around e.

	if( Step(xphi) .eq. 0 ) fac = 2*pi*fac
	Show(xphi) = 2*pi*Var(xphi)
	cosphi = cos(Show(xphi))
	sinphi = sin(Show(xphi))

	fx = ex
	fy = ey
	fz = ez
	fxz = fx**2 + fz**2
	fyz = fy**2 + fz**2
	if( fyz .lt. fxz ) then
	  fxz = sqrt(fxz)
	  nx = (-cosphi*fx*fy - sinphi*fz)/fxz
	  ny = cosphi*fxz
	  nz = (-cosphi*fy*fz + sinphi*fx)/fxz
	else
	  fyz = sqrt(fyz)
	  nx = cosphi*fyz
	  ny = (-cosphi*fx*fy + sinphi*fz)/fyz
	  nz = (-cosphi*fx*fz - sinphi*fy)/fyz
	endif

	r = p/minv
	r0 = sqrt(r**2 + 1)
	m2 = m**2

	msum = msum - m
	bmrem(1) = msum
	bmrem(2) = minv - m

* direct cut on mrem
	if( Test(CUT_MREM, CUT_MIN) ) then
	  Zero(CUT_MREM, CUT_MIN)
	  bmrem(1) = max(bmrem(1), CutMin(xmrem))
	endif
	if( Test(CUT_MREM, CUT_MAX) ) then
	  Zero(CUT_MREM, CUT_MAX)
	  bmrem(2) = min(bmrem(2), CutMax(xmrem))
	endif

	if( bmrem(2) .le. bmrem(1) ) return

	bcosth(1) = -1
	bcosth(2) = 1

* direct cut on costh
	if( Test(CUT_COSTH, CUT_MIN) ) then
	  Zero(CUT_COSTH, CUT_MIN)
	  bcosth(1) = max(bcosth(1), CutMin(xcosth))
	endif
	if( Test(CUT_COSTH, CUT_MAX) ) then
	  Zero(CUT_COSTH, CUT_MAX)
	  bcosth(2) = min(bcosth(2), CutMax(xcosth))
	endif

	if( bcosth(2) .le. bcosth(1) ) return
	bcosth(2) = bcosth(2) - bcosth(1)
	if( Step(xcosth) .eq. 0 ) fac = fac*bcosth(2)
	costh = bcosth(1) + Var(xcosth)*bcosth(2)
	Show(xcosth) = costh

	r0costh = r0*costh
	sinth = sqrt(sqdiff(1D0, costh))
	rz = r*fz
	r1 = r**2*(1 + costh**2) + 1

* indirect cuts on mrem:
* - energy cut
* - momentum cut
* - transverse-energy cut
* - transverse-momentum cut
* - rapidity cut
* - pseudo-rapidity cut

	do c = CUT_MIN + signbit(Bits(CUT_MREM_INDIRECT, CUT_MIN) - 1),
     &         CUT_MAX - signbit(Bits(CUT_MREM_INDIRECT, CUT_MAX) - 1)
	  cut = CutVar(c,xmrem)
	  if( Test(CUT_MREM_E, c) ) then
	    c0 = sqdiff(cut, m*rpart(kQ(c)))
	    c1 = r*r0costh
	    c2 = r1
	    Zero(CUT_MREM_K, c)
	  else if( Test(CUT_MREM_ET, c) ) then
* caveat: the c0 here corresponds to the 'old' definition
* momspec(SPEC_ET,i) = sqrt(momspec(SPEC_KT,i)**2 + m2)
	    c0 = sqdiff(cut, m*rpart(kQ(c))) + m2*rz**2
	    c2 = r0costh*fz + sinth*nz
	    c1 = r*(r0costh - fz*c2)
	    c2 = sqdiff(sqrt(r1), c2) - rz**2
	    Zero(CUT_MREM_KT, c)
	  else if( Test(CUT_MREM_RAP, c) ) then
	    cut = tanh(cut)
	    c0 = m2*sqdiff(rpart(kQ(c))*cut, rz)
	    c2 = r0costh*fz + sinth*nz
	    c1 = r*(fz*c2 - r0costh*cut**2)
	    c2 = sqdiff(c2, cut*sqrt(r1)) + rz**2
	    Zero(CUT_MREM_PRAP, c)
	  else
	    goto 1
	  endif

	  k2max = .25D0*sqdiff(minv - m, bmrem(1))*
     &                  sqdiff(minv + m, bmrem(1))/minv**2

* -int(.) here because mrem and k2 move oppositely
	  j = ieor(c, ibits(-int(sign(1D0, k2eq(k2max))), 1, 2))

	  d = c1**2*(c0*(c0 + c2*m2) + (c1*m2)**2)
	  if( d .lt. 0 ) goto (999, 1) j
	  d = sqrt(d)

	  n = 0
	  do s = 1, 2
	    d = -d
	    k2 = (c0*c2 + 2*(c1**2*m2 + d))/sqdiff(c2, 2*c1)
	    if( k2 .gt. 0 .and. k2 .lt. k2max .and.
     &          abs(k2eq(k2)) .lt. eps ) then
	      n = n + 1
	      cmrem(n) = sqrt(m2 + minv*(minv - 2*sqrt(k2 + m2)))
	    endif
	  enddo

	  if( n .eq. 2 .and. abs(cmrem(2) - cmrem(1)) .gt. eps ) then
	    if( j .eq. 2 ) then
	      print *, "Cannot satisfy direct cut, use veto cut"
	      print *, "bnd mrem = ", bmrem
	      print *, "cut mrem = ", cmrem
	    else
	      bmrem(1) = max(bmrem(1), cmrem(1))
	      bmrem(2) = min(bmrem(2), cmrem(2))
	    endif
	  else if( n .ne. 0 ) then
	    if( (2*j - 3)*(bmrem(j) - cmrem(1)) .gt. 0 )
     &        bmrem(j) = cmrem(1)
	  else
	    goto (999, 1) j
	  endif

	  if( bmrem(2) .ge. bmrem(1) ) goto 1
999	  return

1	  continue
	enddo

	bmrem(2) = bmrem(2) - bmrem(1)
	if( Step(xmrem) .eq. 0 ) fac = fac*bmrem(2)
	mrem = bmrem(1) + Var(xmrem)*bmrem(2)
	Show(xmrem) = mrem

	k10 = .5D0*(minv - sqdiff(mrem, m)/minv)
	k = sqrt(sqdiff(k10, m))
	cutfac = .5D0*k*fac

* k is in the CMS of the decay minv -> m + mrem and must finally
* be boosted into the CMS of the entire process

	qe = k10*r + k*r0costh
	qn = k*sinth
	q = sqrt(qe**2 + qn**2)
	if( q .ne. 0 ) then
	  fx = (qe*fx + qn*nx)/q
	  fy = (qe*fy + qn*ny)/q
	  fz = (qe*fz + qn*nz)/q
	endif
	call VecSet(i, m, q, fx, fy, fz)

	fx = p*ex - q*fx
	fy = p*ey - q*fy
	fz = p*ez - q*fz
	p = sqrt(fx**2 + fy**2 + fz**2)
	if( p .ne. 0 ) then
	  ex = fx/p
	  ey = fy/p
	  ez = fz/p
	endif

	minv = mrem

	if( key .ne. 0 ) print *, "Warning: Superfluous bits in key"
	end

