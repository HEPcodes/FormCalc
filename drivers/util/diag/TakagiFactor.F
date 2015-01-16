* TakagiFactor.F
* computes the Takagi factorization of a complex symmetric matrix
* code adapted from the "Handbook" routines
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 22 Dec 06 th

#include "config.h"


************************************************************************
** TakagiFactor factorizes a complex symmetric n-by-n matrix
** Input: n, A = n-by-n matrix, complex symmetric
** (only the upper triangle of A needs to be filled).
** Output: d = vector of diagonal values, U = transformation matrix
** these fulfill diag(d) = U^* A U^+ with U unitary.

	subroutine TakagiFactor(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	double complex A(ldA,*), U(ldU,*)
	double precision d(*)

	integer sweep, p, q, j
	double precision test, threshold, absA, absB, absX
	double precision fsq, t, invc, s
	double complex f, x, y
	double complex ev(2,MAXDIM)

	if( n .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = Element(A,p,p)
	enddo

	do p = 1, n
	  do q = 1, n
	    Element(U,q,p) = 0
	  enddo
	  Element(U,p,p) = 1
	enddo

	do sweep = 1, 50
	  test = 0
	  do q = 2, n
	    do p = 1, q - 1
	      test = test + abs(Element(A,p,q))
	    enddo
	  enddo
	  if( test .lt. EPS ) goto 1

	  threshold = 0
	  if( sweep .lt. 4 ) threshold = test/(5*n**2)

	  do q = 2, n
	    do p = 1, q - 1
	      absA = abs(ev(2,p))
	      absB = abs(ev(2,q))
	      absX = abs(Element(A,p,q))
	      test = 100*absX
	      if( sweep .gt. 4 .and.
     &            absA + test .eq. absA .and.
     &            absB + test .eq. absB ) then
	        Element(A,p,q) = 0
	      else
	        if( absX .gt. threshold ) then
	          f = sign(1D0, absA - absB)*
     &              (ev(2,q)*DCONJG(Element(A,p,q)) +
     &               DCONJG(ev(2,p))*Element(A,p,q))
	          fsq = DBLE(f)**2 + DIMAG(f)**2
	          if( fsq .eq. 0 ) then
	            f = 1
	            fsq = 1
	          endif
	          t = .5D0*abs(absA - absB)*(absA + absB)
	          t = 1/(t + sqrt(t**2 + fsq))

	          ev(1,p) = ev(1,p) + t*Element(A,p,q)*DCONJG(f)
	          ev(2,p) = Element(A,p,p) + ev(1,p)
	          ev(1,q) = ev(1,q) - t*Element(A,p,q)*f
	          ev(2,q) = Element(A,q,q) + ev(1,q)

	          invc = sqrt(fsq*t**2 + 1)
	          s = t/invc
	          t = t*fsq/(invc + 1)

	          do j = 1, p - 1
	            x = Element(A,j,p)
	            y = Element(A,j,q)
	            Element(A,j,p) = x + s*(DCONJG(f)*y - t*x)
	            Element(A,j,q) = y - s*(f*x + t*y)
	          enddo

	          do j = p + 1, q - 1
	            x = Element(A,p,j)
	            y = Element(A,j,q)
	            Element(A,p,j) = x + s*(DCONJG(f)*y - t*x)
	            Element(A,j,q) = y - s*(f*x + t*y)
	          enddo

	          do j = q + 1, n
	            x = Element(A,p,j)
	            y = Element(A,q,j)
	            Element(A,p,j) = x + s*(DCONJG(f)*y - t*x)
	            Element(A,q,j) = y - s*(f*x + t*y)
	          enddo

	          do j = 1, n
	            x = Element(U,p,j)
	            y = Element(U,q,j)
	            Element(U,p,j) = x + s*(f*y - t*x)
	            Element(U,q,j) = y - s*(DCONJG(f)*x + t*y)
	          enddo

	          Element(A,p,q) = 0
	        endif
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    Element(A,p,p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in TakagiFactor"

1	continue

* make the diagonal elements nonnegative

	do p = 1, n
	  d(p) = abs(Element(A,p,p))
	  if( d(p) .gt. EPS ) then
	    f = sqrt(Element(A,p,p)/d(p))
	    do q = 1, n
	      Element(U,p,q) = Element(U,p,q)*f
	    enddo
	  endif
	enddo

	if( sort .eq. 0 ) return

* sort the eigenvalues

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(t - d(q)) .gt. 0 ) then
	      j = q
	      t = d(q)
	    endif
	  enddo

	  if( j .ne. p ) then
	    d(j) = d(p)
	    d(p) = t
	    do q = 1, n
	      x = Element(U,p,q)
	      Element(U,p,q) = Element(U,j,q)
	      Element(U,j,q) = x
	    enddo
	  endif
	enddo
	end

