* vegas.F
* Vegas Monte-Carlo integration of a vector either with ordinary random
* numbers or with quasi-random numbers
* based on code by M. Martinez, J. Illana, J. Bossert, and A. Vicini
* this file is part of FormCalc
* last modified 24 Jan 03 th

* NDIM is the maximum number of dimensions the integrand may have.
* Note: if your Fortran compiler allows, replace NDIM by ndim
* (i.e. #define NDIM ndim) to dimension the arrays dynamically.
c#define NDIM 3

* NCOMP is the number of components of the integrand vector.
c#define NCOMP 2

* If DEBUG is defined, the error is printed out after each iteration
* so that one can tune the parameters for a particular integral.
#define DEBUG

* RNG determines which random-number generator is used:
* - 1 uses a subtractive Mitchell-Moore generator,
* - 2 uses a Faure quasi-random sequence,
* - 3 uses a Sobol quasi-random sequence (default).
#define RNG 3


************************************************************************
** Vegas integrates a vector in the ndim-dimensional unit hypercube
** using the Monte-Carlo method.  The NCOMP-dimensional integrand is
** invoked via the subroutine integrand(ndim, x, NCOMP, result).  After
** sampling ncall points, the grid is refined.  This is called one
** iteration.  The iterations loop terminates if the absolute value
** of the integral is below absacc or the relative error is below
** relacc, or after maxiter iterations.

	subroutine Vegas(integrand, result, absacc, relacc, ndim,
     &    ncall, maxiter, neval)
	implicit none
	double precision result(NCOMP), relacc, absacc
	integer ndim, ncall, maxiter, neval
	external integrand

	double precision fun(NCOMP), sfun(NCOMP), sfun2(NCOMP)
	double precision sint(NCOMP), sint2(NCOMP), sweight(NCOMP)
	double precision fun2, weight
	double precision r, dr, xo, xn, err
	integer iter, calls, dim, grid, g, c, cmax

	integer ngrid
	parameter (ngrid = 100)

	double precision xi(ngrid, NDIM), d(ngrid, NDIM, NCOMP)
	double precision x(NDIM), imp(ngrid), tmp(ngrid - 1)
	integer pos(NDIM)

	call IniRandom(maxiter*ncall, ndim)
	iter = 0
	neval = 0
	do c = 1, NCOMP
	  sint(c) = 0
	  sint2(c) = 0
	  sweight(c) = 0
	enddo

* define the initial distribution of intervals
	do grid = 1, ngrid
	  r = dble(grid)/ngrid
	  do dim = 1, ndim
	    xi(grid, dim) = r
	  enddo
	enddo

* iterations loop
1	continue
	iter = iter + 1

* initialize iteration variables
	do c = 1, NCOMP
	  sfun(c) = 0
	  sfun2(c) = 0
	  do dim = 1, ndim
	    do grid = 1, ngrid
	      d(grid, dim, c) = 0
	    enddo
	  enddo
	enddo

	do calls = 1, ncall
	  weight = 1D0/ncall

* compute the point position
	  call GetRandom(x)
	  do dim = 1, ndim
	    r = x(dim)*ngrid + 1
	    grid = int(r)
	    xo = 0
	    if( grid .gt. 1 ) xo = xi(grid - 1, dim)
	    xn = xi(grid, dim) - xo
	    x(dim) = xo + (r - grid)*xn
	    pos(dim) = grid
	    weight = weight*xn*ngrid
	  enddo

* compute the integrand
	  call integrand(ndim, x, NCOMP, fun)
	  do c = 1, NCOMP
	    fun2 = fun(c)*weight
	    sfun(c) = sfun(c) + fun2
	    fun2 = fun2**2
	    sfun2(c) = sfun2(c) + fun2
	    do dim = 1, ndim
	      d(pos(dim), dim, c) = d(pos(dim), dim, c) + fun2
	    enddo
	  enddo
	enddo
	neval = neval + ncall

* compute the integral and error values
	err = 0
	do c = 1, NCOMP
	  fun2 = sfun(c)**2
	  sint2(c) = sint2(c) + fun2
	  r = sfun2(c)*ncall - fun2
	  if( r .ne. 0 ) then
	    weight = fun2/abs(r)*(ncall - 1)
	    sweight(c) = sweight(c) + weight
	    sint(c) = sint(c) + sfun(c)*weight
	  endif
	  if( sweight(c) .eq. 0 ) then
	    result(c) = 0
	  else
	    r = sint(c)/sweight(c)
	    result(c) = r
* if the integrand is very close to zero, it is pointless (and costly)
* to insist on a certain relative accuracy
	    if( abs(r) .gt. absacc ) then
              r = sint2(c)/(sint(c)*r)
	      if( r .gt. err ) then
	        err = r
	        cmax = c
	      endif
	    endif
	  endif
	enddo
	err = sqrt(err/iter)

#ifdef DEBUG
	print *, "iteration ", iter, "  error ", err
#endif
	if( err .lt. relacc ) return
	if( iter .gt. maxiter ) then
	  print *, "VEGAS failed to reach the desired accuracy."
	  print *, "Remaining relative error: ", err
	  return
	endif

* redefine the grid (importance sampling)
* - smooth the f^2 value stored for each interval
	do dim = 1, ndim
	  xo = d(1, dim, cmax)
	  xn = d(2, dim, cmax)
	  d(1, dim, cmax) = .5D0*(xo + xn)
	  x(dim) = d(1, dim, cmax)
	  do grid = 2, ngrid - 1
	    r = xo + xn
	    xo = xn
	    xn = d(grid + 1, dim, cmax)
	    d(grid, dim, cmax) = (r + xn)/3D0
	    x(dim) = x(dim) + d(grid, dim, cmax)
	  enddo
	  d(ngrid, dim, cmax) = .5D0*(xo + xn)
	  x(dim) = x(dim) + d(ngrid, dim, cmax)
	enddo

* - compute the importance function of each interval
	do dim = 1, ndim
	  r = 0
	  do grid = 1, ngrid
	    imp(grid) = 0
	    if( d(grid, dim, cmax) .gt. 0 ) then
	      xo = x(dim)/d(grid, dim, cmax)
	      imp(grid) = ((xo - 1)/xo/log(xo))**1.5D0
	    endif
	    r = r + imp(grid)
	  enddo
	  r = r/ngrid

* - redefine the size of each interval
	  dr = 0
	  xn = 0
	  g = 0
	  do grid = 1, ngrid - 1
	    do while( dr .lt. r )
	      g = g + 1
	      dr = dr + imp(g)
	      xo = xn
	      xn = xi(g, dim)
	    enddo
	    dr = dr - r
	    tmp(grid) = xn - (xn - xo)*dr/imp(g)
	  enddo
	  do grid = 1, ngrid - 1
	    xi(grid, dim) = tmp(grid)
	  enddo
	  xi(ngrid, dim) = 1
	enddo

	goto 1
	end


#if RNG == 1

************************************************************************
** IniRandom sets up the random-number generator to produce at most
** max dims-dimensional quasi-random vectors

	subroutine IniRandom(max, dims)
	implicit none
	integer max, dims

	integer ndim
	common /rngdata/ ndim

	ndim = dims
	end


************************************************************************
** GetRandom is a subtractive Mitchell-Moore random-number generator.
** The algorithm is n(i) = (n(i - 24) - n(i - 55)) mod m, implemented
** as a circular array with n(i + 55) = n(i) and m = 2^30 in this
** version.  The array n has been initialized by setting n(i) = i and
** running the algorithm 100,000 times.  Code by Ronald Kleiss.

	subroutine GetRandom(array)
	implicit none
	double precision array(*)

	integer ndim
	common /rngdata/ ndim

	integer dim, a, b, j, m, n(55)
	parameter (m = 2**30)
	data a /55/, b /31/
	data n /
     &    980629335, 889272121, 422278310,1042669295, 531256381,
     &    335028099,  47160432, 788808135, 660624592, 793263632,
     &    998900570, 470796980, 327436767, 287473989, 119515078,
     &    575143087, 922274831,  21914605, 923291707, 753782759,
     &    254480986, 816423843, 931542684, 993691006, 343157264,
     &    272972469, 733687879, 468941742, 444207473, 896089285,
     &    629371118, 892845902, 163581912, 861580190,  85601059,
     &    899226806, 438711780, 921057966, 794646776, 417139730,
     &    343610085, 737162282,1024718389,  65196680, 954338580,
     &    642649958, 240238978, 722544540, 281483031,1024570269,
     &    602730138, 915220349, 651571385, 405259519, 145115737 /

	do dim = 1, ndim
	  a = mod(a, 55) + 1
	  b = mod(b, 55) + 1
	  j = n(b) - n(a)
	  if( j .lt. 0 ) j = j + m
	  n(a) = j
	  array(dim) = dble(j)/m
	enddo
	end

#elif RNG == 2

************************************************************************
** IniRandom sets up the random-number generator to produce a Faure
** sequence of at most max dims-dimensional quasi-random vectors.
** Adapted from ACM TOMS algorithm 659, see
** http://www.acm.org/pubs/citations/journals/toms/1988-14-1/p88-bratley

	subroutine IniRandom(max, dims)
	implicit none
	integer max, dims

	integer coeff(0:19, 0:19)
	integer ndim, prime, nextn, testn, digits
	common /rngdata/ coeff, ndim, prime, nextn, testn, digits

	integer i, j, h

	integer primes(40)
	save primes

	data primes /
     &    1, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11,
     &    13, 13, 17, 17, 17, 17, 19, 19,
     &    23, 23, 23, 23, 29, 29, 29, 29,
     &    29, 29, 31, 31, 37, 37, 37, 37,
     &    37, 37, 41, 41, 41 /

	ndim = dims
	prime = primes(dims)
	testn = prime**4
	nextn = testn - 1
	digits = 3

	h = nint(log(dble(max + testn))/log(dble(prime)))
	coeff(0, 0) = 1
	do j = 1, h
	  coeff(j, 0) = 1
	  coeff(j, j) = 1
	enddo
	do j = 1, h
	  do i = j + 1, h
	    coeff(i, j) =
     &        mod(coeff(i - 1, j) + coeff(i - 1, j - 1), prime)
	  enddo
	enddo
	end


************************************************************************
** GetRandom generates a vector of random numbers

	subroutine GetRandom(array)
	implicit none
	double precision array(*)

	integer coeff(0:19, 0:19)
	integer ndim, prime, nextn, testn, digits
	common /rngdata/ coeff, ndim, prime, nextn, testn, digits

	integer y(0:19), digit, dim, d, j, p
	double precision r

	p = testn
	j = nextn
	do digit = digits, 0, -1
	  p = p/prime
	  d = mod(j, p)
	  y(digit) = (j - d)/p
	  j = d
	enddo

	r = 0
	do digit = digits, 0, -1
	  r = (r + y(digit))/prime
	enddo
	array(1) = r

	do dim = 2, ndim
	  r = 0
	  p = 1
	  do digit = 0, digits
	    j = 0
	    do d = digit, digits
	      j = j + coeff(d, digit)*y(d)
	    enddo
	    y(digit) = mod(j, prime)
	    p = p*prime
	    r = r + dble(y(digit))/p
	  enddo
	  array(dim) = r
	enddo

	nextn = nextn + 1
	if( nextn .eq. testn ) then
	   testn = testn*prime
	   digits = digits + 1
	endif
	end

#else

************************************************************************
** IniRandom sets up the random-number generator to produce a Sobol
** sequence of at most max dims-dimensional quasi-random vectors.
** Adapted from ACM TOMS algorithm 659, see
** http://www.acm.org/pubs/citations/journals/toms/1988-14-1/p88-bratley

	subroutine IniRandom(max, dims)
	implicit none
	integer max, dims

	integer v(40, 30), lastq(40)
	integer ndim, count, norm
	common /rngdata/ v, lastq, ndim, count, norm

	integer bits, powers, degree, newv, dim, bit, deg, j

	integer poly(2:40), vinit(2:40, 1:8)
	save poly, vinit

	data poly / 3, 7, 11, 13, 19, 25, 37, 59, 47,
     &    61, 55, 41, 67, 97, 91, 109, 103, 115, 131,
     &    193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
     &    213, 191, 253, 203, 211, 239, 247, 285, 369, 299 /

	data (vinit(dim, 1), dim = 2, 40) / 39*1 /
	data (vinit(dim, 2), dim = 3, 40) /
     &          1, 3, 1, 3, 1, 3, 3, 1,
     &    3, 1, 3, 1, 3, 1, 1, 3, 1, 3,
     &    1, 3, 1, 3, 3, 1, 3, 1, 3, 1,
     &    3, 1, 1, 3, 1, 3, 1, 3, 1, 3 /
	data (vinit(dim, 3), dim = 4, 40) /
     &             7, 5, 1, 3, 3, 7, 5,
     &    5, 7, 7, 1, 3, 3, 7, 5, 1, 1,
     &    5, 3, 3, 1, 7, 5, 1, 3, 3, 7,
     &    5, 1, 1, 5, 7, 7, 5, 1, 3, 3 /
	data (vinit(dim, 4), dim = 6, 40) /
     &                  1, 7, 9, 13, 11,
     &    1, 3, 7, 9, 5, 13, 13, 11, 3, 15,
     &    5, 3, 15, 7, 9, 13, 9, 1, 11, 7,
     &    5, 15, 1, 15, 11, 5, 3, 1, 7, 9 /
	data (vinit(dim, 5), dim = 8, 40) /
     &                            9, 3, 27,
     &    15, 29, 21, 23, 19, 11, 25, 7, 13, 17,
     &    1, 25, 29, 3, 31, 11, 5, 23, 27, 19,
     &    21, 5, 1, 17, 13, 7, 15, 9, 31, 9 /
	data (vinit(dim, 6), dim = 14, 40) /
     &            37, 33, 7, 5, 11, 39, 63,
     &    27, 17, 15, 23, 29, 3, 21, 13, 31, 25,
     &    9, 49, 33, 19, 29, 11, 19, 27, 15, 25 /
	data (vinit(dim, 7), dim = 20, 40) /
     &                                  13,
     &    33, 115, 41, 79, 17, 29, 119, 75, 73, 105,
     &    7, 59, 65, 21, 3, 113, 61, 89, 45, 107 /
	data (vinit(dim, 8), dim = 38, 40) / 7, 23, 39 /

	j = max
	bits = 0
	do while( j .ne. 0 )
	  bits = bits + 1
	  j = ishft(j, -1)
	enddo

	do bit = 1, bits
	  v(1, bit) = 1
	enddo

	do dim = 2, dims
	  powers = poly(dim)

	  j = powers
	  degree = -1
	  do while( j .ne. 0 )
	    degree = degree + 1
	    j = ishft(j, -1)
	  enddo

	  do bit = 1, degree
	    v(dim, bit) = vinit(dim, bit)
	  enddo

	  do bit = degree + 1, bits
	    newv = v(dim, bit - degree)
	    j = powers
	    do deg = degree, 1, -1
	      if( btest(j, 0) )
     &          newv = ieor(newv, ishft(v(dim, bit - deg), deg))
	      j = ishft(j, -1)
	    enddo
	    v(dim, bit) = newv
	  enddo
	enddo

	do bit = 1, bits - 1
	  do dim = 1, dims
	    v(dim, bit) = ishft(v(dim, bit), bits - bit)
	  enddo
	enddo
	norm = ishft(1, bits)

	count = 0
	ndim = dims
	do dim = 1, dims
	  lastq(dim) = 0
	enddo
	end


************************************************************************
** GetRandom generates a vector of random numbers

	subroutine GetRandom(array)
	implicit none
	double precision array(*)

	integer v(40, 30), lastq(40)
	integer ndim, count, norm
	common /rngdata/ v, lastq, ndim, count, norm

	integer c, zerobit, dim

	c = count
	zerobit = 1
	do while( btest(c, 0) )
	  zerobit = zerobit + 1
	  c = ishft(c, -1)
	enddo

	do dim = 1, ndim
	  lastq(dim) = ieor(lastq(dim), v(dim, zerobit))
	  array(dim) = dble(lastq(dim))/norm
	enddo

	count = count + 1
	end

#endif

