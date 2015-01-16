* ndimgauss_vector.f
* this file is part of FormCalc
* last revision: 28 Feb 00 Thomas Hahn
*
* n-dimensional Gaussian integration for a vector of functions.
*	Note: this file works only with compilers that accept
*	      recursive functions!
*
* - The integration interval in each dimension is assumed to be
*   [0, 1] (similar to VEGAS).
*
* - The function to be integrated is calculated in the
*		subroutine func(result, intvar, from, to)
*		integer from, to
*		double precision result, intvar(10)
*		dimension res(to)
*   which receives the integration variables in the vector intvar, filled
*   with "maxdim" entries such that intvar(1) is the innermost integration
*   variable. func is expected to fill the elements "from" to "to" of the
*   vector "res".
*
* - gauss returns its results in the vector "result". Only the elements
*   "from" to "to" are filled.
*
* - gaussNDim.f contains the arrays of sampling points and weights
*   calculated by GaussPointsNDim.m. The number of sampling points and
*   can be set in the makefile.
*
* - "dim" is mainly for the internal recursion. The user calls the
*   subroutine with dim = maxdim.
*
* Integrating scalar functions can be treated as a special case of
* an n-dimensional integration with n = 1, for instance to calculate
* \int_0^1 dx dy dz sin(x) sin(y) sin(z) one might use the following
* sample invocation:
*
* 	subroutine integrand(result, X)
*	double precision result, X(10)
*	result = sin(X(3)) * sin(X(2)) * sin(X(1))
*	end
*
*	program main
*	double precision integral
*	call gauss(integral, integrand, 3, 3, 1, 1)
*	print *, integral
*	end
*
* Most Fortran compilers allow the arguments "from" and "to" to be omitted
* in the definition of the integrand function, i.e. the subroutine
* "integrand" defined above should really have four arguments:
* res, X, from, to.


	subroutine gauss(result, func, dim, maxdim, from, to)
	implicit none
	integer dim, maxdim, from, to
	double precision result
	dimension result(to)
	external func

	integer j, k
	double precision dx
	double precision left, right
	dimension left(to), right(to)
	double precision intvar(10)
	common /static/ intvar

	automatic j, dx, left, right

	include "GaussNDim.f"

	if(dim .eq. 0) then
	  call func(result, intvar, from, to)
	else 
	  do k = from, to
	    result(k) = 0
	  enddo
	  do j = 1, points(dim)
	    dx = .5D0*gauss_x(dim, j)

*		may want to delete this in final version:
*	    if(dim .ge. 4) print *, "lev ", dim, ": dx=", dx

	    intvar(maxdim + 1 - dim) = .5D0 + dx
	    call gauss(left, func, dim - 1, maxdim, from, to)
	    intvar(maxdim + 1 - dim) = .5D0 - dx
	    call gauss(right, func, dim - 1, maxdim, from, to)
	    do k = from, to
	      result(k) = result(k) +
     +          gauss_w(dim, j)*(left(k) + right(k))
	    enddo
	  enddo
	endif
	end
