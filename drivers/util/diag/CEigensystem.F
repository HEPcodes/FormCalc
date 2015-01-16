* CEigensystem.F
* diagonalization of a complex n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 22 Dec 06 th

#include "config.h"


************************************************************************
** CEigensystem diagonalizes a general complex n-by-n matrix.
** Input: n, A = n-by-n matrix
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^-1.

	subroutine CEigensystem(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	double complex A(ldA,*), U(ldU,*), d(*)

	integer sweep, p, q, j
	double precision test, threshold, absA, absB, absX
	double complex delta, t, invc, sx, sy, tx, ty
	double complex x, y
	double complex ev(2,MAXDIM)

	if( n .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = Element(A,p,p)
	  d(p) = ev(2,p)
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
	      test = test + abs(Element(A,p,q)) + abs(Element(A,q,p))
	    enddo
	  enddo
	  test = .5D0*test
	  if( test .lt. EPS ) goto 1

	  threshold = 0
	  if( sweep .lt. 4 ) threshold = test/(5*n**2)

	  do q = 2, n
	    do p = 1, q - 1
	      absA = abs(ev(2,p))
	      absB = abs(ev(2,q))
	      absX = .5D0*(abs(Element(A,p,q)) + abs(Element(A,q,p)))
	      test = 100*absX
	      if( sweep .gt. 4 .and.
     &            absA + test .eq. absA .and.
     &            absB + test .eq. absB ) then
	        Element(A,p,q) = 0
	        Element(A,q,p) = 0
	      else
	        if( absX .gt. threshold ) then
	          delta = Element(A,p,q)*Element(A,q,p)
	          x = .5D0*(ev(2,p) - ev(2,q))
	          y = sqrt(x**2 + delta)
	          t = x - y
	          x = x + y
	          if( abs(t) .lt. abs(x) ) t = x

	          t = 1/t
	          delta = delta*t
	          ev(1,p) = ev(1,p) + delta
	          ev(2,p) = d(p) + ev(1,p)
	          ev(1,q) = ev(1,q) - delta
	          ev(2,q) = d(q) + ev(1,q)

	          invc = sqrt(delta*t + 1)
	          x = t/invc
	          t = t/(invc + 1)
	          sx = x*Element(A,p,q)
	          ty = t*Element(A,p,q)
	          sy = x*Element(A,q,p)
	          tx = t*Element(A,q,p)

	          do j = 1, n
	            x = Element(A,j,p)
	            y = Element(A,j,q)
	            Element(A,j,p) = x + sy*(y - ty*x)
	            Element(A,j,q) = y - sx*(x + tx*y)
	            x = Element(A,p,j)
	            y = Element(A,q,j)
	            Element(A,p,j) = x + sx*(y - tx*x)
	            Element(A,q,j) = y - sy*(x + ty*y)
	          enddo

	          do j = 1, n
	            x = Element(U,p,j)
	            y = Element(U,q,j)
	            Element(U,p,j) = x + sx*(y - tx*x)
	            Element(U,q,j) = y - sy*(x + ty*y)
	          enddo

	          Element(A,p,q) = 0
	          Element(A,q,p) = 0
	        endif
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in CEigensystem"

1	if( sort .eq. 0 ) return

* sort the eigenvalues by their real part

	do p = 1, n - 1
	  j = p
	  x = d(p)
	  do q = p + 1, n
	    if( sort*(DBLE(x) - DBLE(d(q))) .gt. 0 ) then
	      j = q
	      x = d(q)
	    endif
	  enddo

	  if( j .ne. p ) then
	    d(j) = d(p)
	    d(p) = x
	    do q = 1, n
	      x = Element(U,p,q)
	      Element(U,p,q) = Element(U,j,q)
	      Element(U,j,q) = x
	    enddo
	  endif
	enddo
	end

