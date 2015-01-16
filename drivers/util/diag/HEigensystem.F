* HEigensystem.F
* diagonalization of a Hermitian n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 22 Dec 06 th

#include "config.h"


************************************************************************
** HEigensystem diagonalizes a Hermitian n-by-n matrix.
** Input: n, A = n-by-n matrix, Hermitian
** (only the upper triangle of A needs to be filled).
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^+ = U A U^-1 with U unitary.

	subroutine HEigensystem(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	double complex A(ldA,*), U(ldU,*)
	double precision d(*)

	integer sweep, p, q, j
	double precision test, threshold, absA, absB, absX
	double precision delta, t, invc, s
	double complex Apq, x, y
	double precision ev(2,MAXDIM)

	if( n .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = DBLE(Element(A,p,p))
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
	          t = .5D0*(ev(2,p) - ev(2,q))
	          t = 1/(t + sign(sqrt(t**2 + absX**2), t))

	          delta = t*absX**2
	          ev(1,p) = ev(1,p) + delta
	          ev(2,p) = d(p) + ev(1,p)
	          ev(1,q) = ev(1,q) - delta
	          ev(2,q) = d(q) + ev(1,q)

	          invc = sqrt(delta*t + 1)
	          s = t/invc
	          t = delta/(invc + 1)

	          Apq = Element(A,p,q)

	          do j = 1, p - 1
	            x = Element(A,j,p)
	            y = Element(A,j,q)
	            Element(A,j,p) = x + s*(DCONJG(Apq)*y - t*x)
	            Element(A,j,q) = y - s*(Apq*x + t*y)
	          enddo

	          do j = p + 1, q - 1
	            x = Element(A,p,j)
	            y = Element(A,j,q)
	            Element(A,p,j) = x + s*(Apq*DCONJG(y) - t*x)
	            Element(A,j,q) = y - s*(Apq*DCONJG(x) + t*y)
	          enddo

	          do j = q + 1, n
	            x = Element(A,p,j)
	            y = Element(A,q,j)
	            Element(A,p,j) = x + s*(Apq*y - t*x)
	            Element(A,q,j) = y - s*(DCONJG(Apq)*x + t*y)
	          enddo

	          do j = 1, n
	            x = Element(U,p,j)
	            y = Element(U,q,j)
	            Element(U,p,j) = x + s*(Apq*y - t*x)
	            Element(U,q,j) = y - s*(DCONJG(Apq)*x + t*y)
	          enddo

	          Element(A,p,q) = 0
	        endif
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in HEigensystem"

1	if( sort .eq. 0 ) return

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

