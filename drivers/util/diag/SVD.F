* SVD.F
* singular value decomposition of an m-by-n matrix
* this file is part of the Diag library
* last modified 21 Aug 15 th

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

	integer nx, nm, p, q, px, qx, j, rev, pi(MAXDIM)
	RealType red, off, thresh
	RealType t, dv, dw, xv, xw, invc
	ComplexType App, Apq, Aqp, Aqq
	ComplexType x, y, sv, sw, tv, tw, f
	ComplexType VW(MAXDIM,MAXDIM,0:2)

* note: for better cache efficiency, the Vx, Wx arrays
* contain the *transpose* of the transformation matrices
	ComplexType V_(MAXDIM,MAXDIM)
	ComplexType W_(MAXDIM,MAXDIM)
	ComplexType A_(MAXDIM,MAXDIM)
	equivalence (VW(1,1,0), V_)
	equivalence (VW(1,1,1), W_)
	equivalence (VW(1,1,2), A_)

	integer sweep
	common /nsweeps/ sweep

	nx = max(m, n)

	if( nx .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, nx
	  do q = 1, nx
	    V_(q,p) = 0
	    W_(q,p) = 0
	    A_(q,p) = 0
	  enddo
	  V_(p,p) = 1
	  W_(p,p) = 1
	enddo

	rev = ibits(m - n, 15, 1)
	if( rev .eq. 1 ) then
	  do p = 1, n
	    do q = 1, m
	      A_(p,q) = A(q,p)
	    enddo
	  enddo
	else
	  do p = 1, n
	    do q = 1, m
	      A_(q,p) = A(q,p)
	    enddo
	  enddo
	endif

	red = .01D0/nx**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, nx
	    do p = 1, q - 1
	      off = off + Sq(A_(p,q)) + Sq(A_(q,p))
	    enddo
	  enddo
	  if( .not. off .gt. EPS ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, nx
	    do p = 1, q - 1
	      px = p
	      qx = q
	      if( Sq(A_(p,p)) + Sq(A_(q,q)) .lt.
     &            Sq(A_(p,q)) + Sq(A_(q,p)) ) then
	        px = q
	        qx = p
	      endif

	      App = A_(px,p)
	      Aqq = A_(qx,q)
	      Apq = A_(px,q)
	      Aqp = A_(qx,p)
	      off = Sq(Apq) + Sq(Aqp)
	      if( sweep .gt. 4 .and.
     &            off .lt. EPS*(Sq(App) + Sq(Aqq)) ) then
	        A_(px,q) = 0
	        A_(qx,p) = 0
	      else if( off .gt. thresh ) then
	        xv = Re((App - Aqq)*Conjugate(App + Aqq))
	        xw = Re((Apq - Aqp)*Conjugate(Apq + Aqp))
	        dv = .5D0*(xv + xw)
	        dw = .5D0*(xv - xw)

	        tv = Conjugate(App)*Aqp + Aqq*Conjugate(Apq)
	        tw = Conjugate(App)*Apq + Aqq*Conjugate(Aqp)
c	        t = sqrt(dv**2 + Sq(tv))
	        t = sqrt(dw**2 + Sq(tw))

	        xv = min(abs(dv + t), abs(dw + t))
	        xw = min(abs(dv - t), abs(dw - t))
	        if( xv + xw .gt. DBL_EPS ) then
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

	      if( p .ne. px ) then
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
	  if( d(p) .gt. DBL_EPS .and. d(p) .ne. Re(A_(p,p)) ) then
	    f = A_(p,p)/d(p)
	    do q = 1, nm
	      W_(q,p) = W_(q,p)*f
	    enddo
	  endif
	enddo

* sort the singular values

	do p = 1, nm
	  pi(p) = p
	enddo

	do p = 1, nm
	  j = p
	  t = d(p)
	  if( sort .ne. 0 ) then
	    do q = p + 1, nm
	      if( sort*(t - d(q)) .gt. 0 ) then
	        j = q
	        t = d(q)
	      endif
	    enddo
	  endif

	  d(j) = d(p)
	  d(p) = t

	  q = pi(j)
	  pi(j) = pi(p)

	  do j = 1, m
	    VL(p,j) = VW(j,q,rev)
	  enddo
	  do j = 1, n
	    WL(p,j) = VW(j,q,1-rev)
	  enddo
	enddo
	end

