* SVD.F
* singular value decomposition of an m-by-n matrix
* this file is part of the Diag library
* last modified 25 Mar 19 th

#include "diag-f.h"


************************************************************************
** SVD performs a singular value decomposition.
** Input:	m, n, A = m-by-n matrix.
** Output:	d = nm-vector of singular values, nm = min(m, n),
** for UCOLS=0:	V = nm-by-m left transformation matrix,
**		W = nm-by-n right transformation matrix,
** which fulfill d = V^* A W^+,  A = V^T d W,  V^* A = d W
** for UCOLS=1:	V = m-by-nm left transformation matrix,
**		W = n-by-nm right transformation matrix,
** which fulfill d = V^+ A W^*,  A = V d W^T,  A W^* = V d

	subroutine SVD(m, n, A,ldA, d, V,ldV, W,ldW, sort)
	implicit none
	integer m, n, ldA, ldV, ldW, sort
	ComplexType A(ldA,*), V(ldV,*), W(ldW,*)
	RealType d(*)

	integer nx, nm, p, q, px, qx, j, rev, pi(min(m,n))
	RealType red, off, thresh
	RealType t, dv, dw, xv, xw, invc
	ComplexType App, Apq, Aqp, Aqq
	ComplexType x, y, sv, sw, tv, tw, f
	ComplexType VW(max(m,n),max(m,n),0:2)

* note: for better cache efficiency, the V_, W_ arrays
* contain the *transpose* of the transformation matrices
#define V_(i,j) VW(i,j,0)
#define W_(i,j) VW(i,j,1)
#define A_(i,j) VW(i,j,2)

	integer sweep
	common /nsweeps/ sweep

	nx = max(m, n)

	VW = 0
	forall(p = 1:nx) VW(p,p,0:1) = 1

	rev = merge(1, 0, m < n)
	if( rev == 1 ) then
	  A_(1:n,1:m) = transpose(A(1:m,1:n))
	else
	  A_(1:m,1:n) = A(1:m,1:n)
	endif

	red = .01D0/nx**4

	do sweep = 1, 50
	  off = sum([(Sq(A_(1:q-1,q)) + Sq(A_(q,1:q-1)), q = 2,n)])
	  if( .not. off > EPS ) goto 1

	  thresh = 0
	  if( sweep < 4 ) thresh = off*red

	  do q = 2, nx
	    do p = 1, q - 1
	      px = p
	      qx = q
	      if( Sq(A_(p,p)) + Sq(A_(q,q)) < Sq(A_(p,q)) + Sq(A_(q,p)) ) then
	        px = q
	        qx = p
	      endif

	      App = A_(px,p)
	      Aqq = A_(qx,q)
	      Apq = A_(px,q)
	      Aqp = A_(qx,p)
	      off = Sq(Apq) + Sq(Aqp)
	      if( sweep > 4 .and.
     &            off < EPS*(Sq(App) + Sq(Aqq)) ) then
	        A_(px,q) = 0
	        A_(qx,p) = 0
	      else if( off > thresh ) then
	        xv = Re((App - Aqq)*Conjugate(App + Aqq))
	        xw = Re((Apq - Aqp)*Conjugate(Apq + Aqp))
	        dv = (xv + xw)/2
	        dw = (xv - xw)/2

	        tv = Conjugate(App)*Aqp + Aqq*Conjugate(Apq)
	        tw = Conjugate(App)*Apq + Aqq*Conjugate(Aqp)
c	        t = sqrt(dv**2 + Sq(tv))
	        t = sqrt(dw**2 + Sq(tw))

	        xv = min(abs(dv + t), abs(dw + t))
	        xw = min(abs(dv - t), abs(dw - t))
	        if( xv + xw > DBL_EPS ) then
	          t = sign(t, xv - xw)
	          tv = tv/(dv + t)
	          tw = tw/(dw + t)
	        else
	          tv = 0
	          tw = Apq/App
	        endif

	        invc = sqrt(1 + Sq(tv))
	        sv = tv/invc
	        tv = tv/(invc + 1)

	        invc = sqrt(1 + Sq(tw))
	        sw = tw/invc
	        tw = tw/(invc + 1)

	        do j = 1, nx
	          x = A_(j,p)
	          y = A_(j,q)
	          A_(j,p) = x + Conjugate(sw)*(y - tw*x)
	          A_(j,q) = y - sw*(x + Conjugate(tw)*y)
	          x = A_(px,j)
	          y = A_(qx,j)
	          A_(p,j) = x + Conjugate(sv)*(y - tv*x)
	          A_(q,j) = y - sv*(x + Conjugate(tv)*y)
	        enddo

	        A_(p,p) = invc*(App + Conjugate(sv)*(Aqp - tv*App))
	        A_(q,q) = invc*(Aqq - sv*(Apq + Conjugate(tv)*Aqq))
	        A_(q,p) = 0
	        A_(p,q) = 0

	        do j = 1, nx
	          x = V_(j,px)
	          y = V_(j,qx)
	          V_(j,p) = x + sv*(y - Conjugate(tv)*x)
	          V_(j,q) = y - Conjugate(sv)*(x + tv*y)
	        enddo

	        do j = 1, nx
	          x = W_(j,p)
	          y = W_(j,q)
	          W_(j,p) = x + sw*(y - Conjugate(tw)*x)
	          W_(j,q) = y - Conjugate(sw)*(x + tw*y)
	        enddo
	        goto 2
	      endif

	      if( p /= px ) then
	        do j = 1, nx
	          x = A_(p,j)
	          A_(p,j) = A_(q,j)
	          A_(q,j) = x
	        enddo

	        do j = 1, nx
	          x = V_(j,p)
	          V_(j,p) = V_(j,q)
	          V_(j,q) = x
	        enddo
	      endif

2	      continue
	    enddo
	  enddo
	enddo

	print *, "Bad convergence in SVD"

1	continue

	nm = min(m, n)

* make the diagonal elements nonnegative

	do p = 1, nm
	  d(p) = abs(A_(p,p))
	  if( d(p) > DBL_EPS .and. d(p) /= Re(A_(p,p)) ) then
	    f = A_(p,p)/d(p)
	    W_(1:nm,p) = W_(1:nm,p)*f
	  endif
	enddo

* sort the singular values

	forall(p = 1:nm) pi(p) = p

	do p = 1, nm
	  j = p
	  t = d(p)
	  if( sort /= 0 ) then
	    do q = p + 1, nm
	      if( sort*(t - d(q)) > 0 ) then
	        j = q
	        t = d(q)
	      endif
	    enddo
	  endif

	  d(j) = d(p)
	  d(p) = t

	  q = pi(j)
	  pi(j) = pi(p)

	  VL(p,1:m) = VW(1:m,q,rev)
	  WL(p,1:n) = VW(1:n,q,1-rev)
	enddo
	end

