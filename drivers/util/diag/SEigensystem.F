* SEigensystem.F
* diagonalization of a complex symmetric n-by-n matrix using
* the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 25 Mar 19 th

#include "diag-f.h"


************************************************************************
** SEigensystem diagonalizes a complex symmetric n-by-n matrix.
** Input:	n, A = n-by-n matrix, complex symmetric
**		(only the upper triangle of A needs to be filled).
** Output:	d = vector of eigenvalues,
**		U = transformation matrix (U^-1 = U^T),
** these fulfill
**	d = U A U^T,  A = U^T d U,  U A = d U  (UCOLS=0),
**	d = U^T A U,  A = U d U^T,  A U = U d  (UCOLS=1).

	subroutine SEigensystem(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	ComplexType A(ldA,n), U(ldU,n), d(n)

	integer p, q, j
	RealType red, off, thresh
	ComplexType delta, t, invc, s
	ComplexType x, y
	ComplexType ev(2,n)

	integer sweep
	common /nsweeps/ sweep

	ev(1,1:n) = 0
	forall(p = 1:n) ev(2,p) = A(p,p)
	d(1:n) = ev(2,1:n)

	U(1:n,1:n) = 0
	forall(p = 1:n) U(p,p) = 1

	red = .04D0/n**4

	do sweep = 1, 50
	  off = sum([(Sq(A(1:q-1,q)), q = 2,n)])
	  if( .not. off > SYM_EPS ) goto 1

	  thresh = 0
	  if( sweep < 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      delta = A(p,q)
	      off = Sq(delta)
	      if( sweep > 4 .and.
     &            off < SYM_EPS*(Sq(ev(2,p)) + Sq(ev(2,q))) ) then
	        A(p,q) = 0
	      else if( off > thresh ) then
	        x = (ev(2,p) - ev(2,q))/2
	        y = sqrt(x**2 + delta**2)
	        t = x - y
	        s = x + y
	        if( Sq(t) < Sq(s) ) t = s

	        t = delta/t
	        delta = delta*t
	        ev(1,p) = ev(1,p) + delta
	        ev(2,p) = d(p) + ev(1,p)
	        ev(1,q) = ev(1,q) - delta
	        ev(2,q) = d(q) + ev(1,q)

	        invc = sqrt(t**2 + 1)
	        s = t/invc
	        t = t/(invc + 1)

	        do j = 1, p - 1
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + s*(y - t*x)
	          A(j,q) = y - s*(x + t*y)
	        enddo

	        do j = p + 1, q - 1
	          x = A(p,j)
	          y = A(j,q)
	          A(p,j) = x + s*(y - t*x)
	          A(j,q) = y - s*(x + t*y)
	        enddo

	        do j = q + 1, n
	          x = A(p,j)
	          y = A(q,j)
	          A(p,j) = x + s*(y - t*x)
	          A(q,j) = y - s*(x + t*y)
	        enddo

	        A(p,q) = 0

	        do j = 1, n
	          x = UL(p,j)
	          y = UL(q,j)
	          UL(p,j) = x + s*(y - t*x)
	          UL(q,j) = y - s*(x + t*y)
	        enddo
	      endif
	    enddo
	  enddo

	  ev(1,1:n) = 0
	  d(1:n) = ev(2,1:n)
	enddo

	print *, "Bad convergence in SEigensystem"

1	if( sort == 0 ) return

* sort the eigenvalues by their real part

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(Re(t) - Re(d(q))) > 0 ) then
	      j = q
	      t = d(q)
	    endif
	  enddo

	  if( j /= p ) then
	    d(j) = d(p)
	    d(p) = t
	    do q = 1, n
	      x = UL(p,q)
	      UL(p,q) = UL(j,q)
	      UL(j,q) = x
	    enddo
	  endif
	enddo
	end

