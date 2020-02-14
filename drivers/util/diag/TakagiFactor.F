* TakagiFactor.F
* computes the Takagi factorization of a complex symmetric matrix
* code adapted from the "Handbook" routines
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 29 Mar 19 th

#include "diag-f.h"


************************************************************************
** TakagiFactor factorizes a complex symmetric n-by-n matrix
** Input:	n, A = n-by-n matrix, complex symmetric
**		(only the upper triangle of A needs to be filled),
** Output:	d = vector of diagonal values,
**		U = transformation matrix, unitary (U^-1 = U^+),
** these fulfill
**	d = U^* A U^+,  A = U^T d U,  U^* A = d U  (UCOLS=0),
**	d = U^+ A U^*,  A = U d U^T,  A U^* = U d  (UCOLS=1).

	subroutine TakagiFactor(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	ComplexType A(ldA,n), U(ldU,n)
	RealType d(n)

	integer p, q, j
	RealType red, off, thresh
	RealType sqp, sqq, t, invc
	ComplexType f, x, y
	ComplexType ev(2,n)

	integer sweep
	common /nsweeps/ sweep

	ev(1,:) = 0
	forall(p = 1:n) ev(2,p) = A(p,p)

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
	      off = Sq(A(p,q))
	      sqp = Sq(ev(2,p))
	      sqq = Sq(ev(2,q))
	      if( sweep > 4 .and. off < SYM_EPS*(sqp + sqq) ) then
	        A(p,q) = 0
	      else if( off > thresh ) then
	        t = abs(sqp - sqq)/2
	        if( t > 0 ) then
	          f = sign(1D0, sqp - sqq)*
     &              (ev(2,q)*Conjugate(A(p,q)) +
     &               Conjugate(ev(2,p))*A(p,q))
	        else
	          f = 1
	          if( sqp /= 0 ) f = sqrt(ev(2,q)/ev(2,p))
	        endif
	        t = t + sqrt(t**2 + Sq(f))
	        f = f/t

	        ev(1,p) = ev(1,p) + A(p,q)*Conjugate(f)
	        ev(2,p) = A(p,p) + ev(1,p)
	        ev(1,q) = ev(1,q) - A(p,q)*f
	        ev(2,q) = A(q,q) + ev(1,q)

	        t = Sq(f)
	        invc = sqrt(t + 1)
	        f = f/invc
	        t = t/(invc*(invc + 1))

	        do j = 1, p - 1
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + (Conjugate(f)*y - t*x)
	          A(j,q) = y - (f*x + t*y)
	        enddo

	        do j = p + 1, q - 1
	          x = A(p,j)
	          y = A(j,q)
	          A(p,j) = x + (Conjugate(f)*y - t*x)
	          A(j,q) = y - (f*x + t*y)
	        enddo

	        do j = q + 1, n
	          x = A(p,j)
	          y = A(q,j)
	          A(p,j) = x + (Conjugate(f)*y - t*x)
	          A(q,j) = y - (f*x + t*y)
	        enddo

	        A(p,q) = 0

	        do j = 1, n
	          x = UL(p,j)
	          y = UL(q,j)
	          UL(p,j) = x + (f*y - t*x)
	          UL(q,j) = y - (Conjugate(f)*x + t*y)
	        enddo
	      endif
	    enddo
	  enddo

	  ev(1,:) = 0
	  forall(p = 1:n) A(p,p) = ev(2,p)
	enddo

	print *, "Bad convergence in TakagiFactor"

1	continue

* make the diagonal elements nonnegative

	do p = 1, n
	  d(p) = abs(A(p,p))
	  if( d(p) > DBL_EPS .and. d(p) /= Re(A(p,p)) ) then
	    f = sqrt(A(p,p)/d(p))
	    UL(p,1:n) = UL(p,1:n)*f
	  endif
	enddo

	if( sort == 0 ) return

* sort the eigenvalues

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(t - d(q)) > 0 ) then
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

