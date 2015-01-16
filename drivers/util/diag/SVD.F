* SVD.F
* singular value decomposition of an m-by-n matrix
* this file is part of the Diag library
* last modified 13 Dec 06 th

#include "config.h"


************************************************************************
** SVD performs a singular value decomposition.
** Input: m, n, A = m-by-n matrix (m > n).
** Output: d = n-vector of singular values,
** V = n-by-m left transformation matrix,
** W = n-by-n right transformation matrix,
** these fulfill diag(d) = V^* A W^+.

	subroutine SVD(m, n, A,ldA, d, V,ldV, W,ldW, sort)
	implicit none
	integer m, n, ldA, ldV, ldW, sort
	double complex A(ldA,*), V(ldV,*), W(ldW,*)
	double precision d(*)

	integer i, j, k, l, nz, iz(MAXDIM)
	logical valid(MAXDIM)
	double complex s
	double precision norm

	double complex AdagA(MAXDIM,MAXDIM)

	if( n .gt. m .or. m .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

* The singular values are the roots of the eigenvalues of A^+ A

	do k = 1, n
	  do i = 1, k
	    s = 0
	    do j = 1, m
	      s = s + DCONJG(Element(A,j,i))*Element(A,j,k)
	    enddo
	    Element(AdagA,i,k) = s
	  enddo
	enddo
	call HEigensystem(n, AdagA,MAXDIM, d, W,ldW, sort)

	nz = 0
	do k = 1, n
	  valid(k) = d(k) .gt. EPS
	  if( valid(k) ) then
	    d(k) = sqrt(d(k))
	    norm = 1/d(k)
	    do i = 1, m
	      s = 0
	      do j = 1, n
	        s = s + Element(A,i,j)*DCONJG(Element(W,k,j))
	      enddo
	      Element(V,k,i) = s*norm
	    enddo
	  else
	    nz = nz + 1
	    iz(nz) = k
	  endif
	enddo

* The vectors corresponding to the zero singular values can be chosen
* arbitrarily.  We choose u_i = delta_ij and use Gram-Schmidt to
* orthonormalize these.

#define Delta(i,j) ibits(ieor(i,j)-1,15,1)

	l = m
	do i = 1, nz
	  norm = 0

1	  do j = 1, m
	    s = 0
	    do k = 1, n
	      if( valid(k) )
     &          s = s + Element(V,k,l)*DCONJG(Element(V,k,j))
	    enddo
	    s = Delta(l, j) - s
	    Element(V,iz(i),j) = s
	    norm = norm + DBLE(s)**2 + DIMAG(s)**2
	  enddo
	  valid(iz(i)) = .TRUE.
	  l = l - 1
	  if( norm .lt. EPS ) goto 1

	  norm = 1/sqrt(norm)
	  k = iz(i)
	  do j = 1, m
	    Element(V,k,j) = norm*Element(V,k,j)
	  enddo
	enddo
	end

