* diag.F
* routines for the diagonalization and SVD of matrices
* based on code by M. Frank and G. Moortgat-Pick
* this file is part of FormCalc
* last modified 10 Jan 03 th

* if DESC is defined, Eigensystem and SingularValues will order
* the eigen- and singular values in descending order
c#define DESC


************************************************************************
* Eigen2x2: diagonalization of a Hermitian 2-by-2 matrix

	subroutine Eigen2x2(m, msq, U, M11, M22, M12, *)
	implicit none
	logical error
	double precision m(2), msq(2)
	double complex U(2, 2)
	double precision M11, M22
	double complex M12

	double precision s
	double complex r1, r2

	s = sqrt((M11 - M22)**2 + 4*(dble(M12)**2 + dimag(M12)**2))
	msq(1) = .5D0*(M11 + M22 - s)
	if( msq(1) .lt. 0 ) then
	  print *, "Negative mass-squares."
	  return 1
	endif
	msq(2) = .5D0*(M11 + M22 + s)

	m(1) = sqrt(msq(1))
	m(2) = sqrt(msq(2))

	if( s .eq. 0 ) then
	  U(1, 1) = 1
	  U(1, 2) = 0
	else if( M11 .eq. msq(2) ) then
	  U(1, 1) = 0
	  U(1, 2) = -1
	else
	  r1 = sqrt(dcmplx(M11 - msq(2)))
	  r2 = dcmplx(0D0, sqrt(s))
	  U(1, 1) = r1/r2
	  U(1, 2) = M12/(r1*r2)
	endif
	U(2, 1) = -dconjg(U(1, 2))
	U(2, 2) = U(1, 1)

	error = .FALSE.
	end


************************************************************************
* Eigensystem: diagonalization of a Hermitian n-by-n matrix
* (originally from EISPACK's htridi, tql2, and htribk)
* input: n, A = n-by-n matrix, Hermitian
* (only the lower triangle of A needs to be filled)
* output: d = vector of eigenvalues, U = transformation matrix
* these fulfill diag(d) = U A U^+

	subroutine Eigensystem(n, A, lda, d, U, ldu)
	implicit none
	integer n, lda, ldu
	double complex A(lda, *), U(ldu, *)
	double precision d(*)

	integer i, j, k, l
	double precision sqsum, r, f, s, c, p, diff, tdiff, ldiff
	double complex g, h

	integer nmax
	parameter (nmax = 4)
	double complex phase(0:nmax), tmp(nmax)
	double precision err(0:nmax), reU(nmax, nmax)

	double precision eps
	parameter (eps = 5D-16)

* reduce A to a real symmetric tridiagonal matrix (htridi)

	do i = 1, n
	  d(i) = dble(A(i, i))
	enddo

	phase(n) = 1

	do i = n, 1, -1
	  l = i - 1

	  sqsum = 0
	  do k = 1, l
	    s = dble(A(i, k))**2 + dimag(A(i, k))**2
	    sqsum = sqsum + s
	  enddo

	  if( sqsum .eq. 0 ) then
	    phase(l) = 1
	    err(l) = 0
	    goto 127
	  endif

	  r = sqrt(sqsum)
	  err(l) = r

	  if( s .eq. 0 ) then
	    phase(l) = -phase(i)
	    A(i, l) = r
	  else
	    s = sqrt(s)
	    phase(l) = -phase(i)*A(i, l)/s
	    A(i, l) = A(i, l)*(1 + r/s)
	    sqsum = sqsum + r*s
	    if( l .eq. 1 ) goto 127
	  endif

	  f = 0
	  do j = 1, l
	    g = 0
	    do k = 1, j
	      g = g + A(i, k)*dconjg(A(j, k))
	    enddo
	    do k = j + 1, l
	      g = g + A(i, k)*A(k, j)
	    enddo
	    tmp(j) = g/sqsum
	    f = f + dble(dconjg(tmp(j))*A(i, j))
	  enddo

	  f = f/(2*sqsum)
	  do j = 1, l
	    tmp(j) = tmp(j) - f*A(i, j)
	    g = dconjg(tmp(j))
	    h = dconjg(A(i, j))
	    do k = 1, j
	      A(j, k) = A(j, k) - h*tmp(k) - g*A(i, k)
	    enddo
	  enddo

127	  continue

	  f = dble(A(i, i))
	  A(i, i) = dcmplx(d(i), sqsum)
	  d(i) = f
	enddo

* diagonalize the reduced matrix using the QL method (tql2)

	do i = 1, n
	  do j = 1, n
	    reU(j, i) = 0
	  enddo
	  reU(i, i) = 1
	enddo

	do l = 1, n

210	  do j = l, n - 1
	    p = abs(d(j)) + abs(d(j + 1))
	    if( abs(err(j)) .lt. eps*p ) goto 211
	  enddo

211	  if( j .gt. l ) then
	    p = d(l)
	    diff = (d(l + 1) - p)/(2*err(l))
	    r = sqrt(diff**2 + 1)
	    diff = d(j) - p + err(l)/(diff + sign(r, diff))

	    s = 1
	    c = 1
	    p = 0

	    do i = j - 1, l, -1
	      tdiff = s*err(i)
	      ldiff = c*err(i)
	      if( abs(tdiff) .lt. abs(diff) ) then
	        s = tdiff/diff
	        r = sqrt(s**2 + 1)
	        err(i + 1) = diff*r
	        c = 1/r
	        s = s*c
	      else
	        c = diff/tdiff
	        r = sqrt(c**2 + 1)
	        err(i + 1) = tdiff*r
	        s = 1/r
	        c = c*s
	      endif
	      diff = d(i + 1) - p
	      r = (d(i) - diff)*s + 2*c*ldiff
	      p = s*r
	      d(i + 1) = diff + p
	      diff = c*r - ldiff
	      do k = 1, n
	        f = reU(i + 1, k)
	        reU(i + 1, k) = s*reU(i, k) + c*f
	        reU(i, k) = c*reU(i, k) - s*f
	      enddo
	    enddo

	    d(l) = d(l) - p
	    err(l) = diff
	    err(j) = 0
	    goto 210
	  endif
	enddo

* sort the eigenvalues

	do i = 1, n - 1
	  k = i
	  f = d(i)
	  do j = i + 1, n
#ifdef DESC
	    if( d(j) .gt. f ) then
#else
	    if( d(j) .lt. f ) then
#endif
	      k = j
	      f = d(j)
	    endif
	  enddo

	  if( k .ne. i ) then
	    d(k) = d(i)
	    d(i) = f
	    do j = 1, n
	      f = reU(i, j)
	      reU(i, j) = reU(k, j)
	      reU(k, j) = f
	    enddo
	  endif
	enddo

	do k = 1, n
	  g = phase(k)
	  do j = 1, n
	    U(j, k) = reU(j, k)*g
	  enddo
	enddo

* transform the eigenvectors of the real symmetric tridiagonal
* matrix to those of the Hermitian tridiagonal matrix (htribk)

	do i = 2, n
	  sqsum = dimag(A(i, i))
	  if( sqsum .ne. 0 ) then
	    l = i - 1
	    do j = 1, n
	      g = 0
	      do k = 1, l
	        g = g + U(j, k)*dconjg(A(i, k))
	      enddo
	      g = g/sqsum
	      do k = 1, l
	        U(j, k) = U(j, k) - g*A(i, k)
	      enddo
	    enddo
	  endif
	enddo
	end


************************************************************************
* SingularValues: singular-value decomposition of an m-by-n matrix
* input: m, n, A = m-by-n matrix (m > n)
* output: d = vector of singular values
* U = n-by-m left transformation matrix
* V = n-by-n right transformation matrix
* these fulfill diag(d) = U A V^+

	subroutine SingularValues(m, n, A, lda, d, U, ldu, V, ldv)
	implicit none
	integer m, n, lda, ldu, ldv
	double complex A(lda, *)
	double precision d(*)
	double complex U(ldu, *), V(ldv, *)

	integer i, j, k, l, r
	double complex s
	double precision x

	integer nmax
	parameter (nmax = 4)
	double complex AdagA(nmax, nmax)

	double precision eps
	parameter (eps = 5D-16)

* The singular values are the roots of the eigenvalues of A^+ A

	do k = 1, n
	  do i = 1, k
	    s = 0
	    do j = 1, m
	      s = s + dconjg(A(j, k))*A(j, i)
	    enddo
	    AdagA(k, i) = s
	  enddo
	enddo
	call Eigensystem(n, AdagA, nmax, d, V, ldv)

#ifdef DESC
	do k = 1, n
#else
	do k = n, 1, -1
#endif
	  if( d(k) .lt. eps ) goto 10
	  d(k) = sqrt(d(k))
	  x = 1/d(k)

	  do i = 1, m
	    s = 0
	    do j = 1, n
	      s = s + V(k, j)*dconjg(A(i, j))
	    enddo
	    U(k, i) = s*x
	  enddo
	enddo

10	continue

* The vectors corresponding to the zero singular values can be chosen
* arbitrarily.  We choose u_i = delta_ij and use Gram-Schmidt to
* orthonormalize these.

	j = 0
#ifdef DESC
	do i = k, n
#else
	do i = k, 1, -1
#endif
11	  j = j + 1
	  x = 0
	  do r = 1, m
	    if( r .eq. j ) then
	      s = 1
	    else
	      s = 0
	    endif
#ifdef DESC
	    do l = 1, i - 1
#else
	    do l = i + 1, k
#endif
	      s = s - dconjg(U(l, j))*U(l, r)
	    enddo
	    U(i, r) = s
	    x = x + dble(s)**2 + dimag(s)**2
	  enddo
	  if( x .le. eps ) goto 11

	  x = 1/sqrt(x)
	  do r = 1, m
	    U(i, r) = U(i, r)*x
	  enddo
	enddo
	end
