* ndimgauss_scalar.f
* last revision: 26 Apr 99 Thomas Hahn
*
* n-dimensional Gaussian integration for a vector of functions.
*	Note: this file works only with compilers that accept
*	      recursive functions!
*
* - The integration interval in each dimension is assumed to be
*   [0, 1] (similar to VEGAS).
*
* - The function to be integrated is calculated in the
*		subroutine func(intvar, res)
*		double precision intvar(10), result
*   which receives the integration variables in the vector intvar, filled
*   with "maxdim" entries such that intvar(1) is the innermost integration
*   variable.
*
* - GaussNDim.f contains the arrays of sampling points and weights
*   calculated by GaussPointsNDim.m. The number of sampling points and
*   can be set in the makefile.
*
* - "dim" is mainly for the internal recursion. The user calls the
*   subroutine with dim = maxdim.
*
* For instance to calculate \int_0^1 dx dy dz sin(x) cos(y) exp(z) one
* might use the following sample invocation:
*
* 	double precision integrand(X)
*	double precision X(10)
*	integrand = sin(X(3)) * cos(X(2)) * exp(X(1))
*	end
*
*	program main
*	double precision integral
*	call gauss(integral, integrand, 3, 3)
*	print *, integral
*	end


	subroutine gauss(result, func, dim, maxdim)
	implicit none
	integer dim, maxdim
	double precision result, func
	external func

	integer j
	double precision dx
	double precision left, right
	double precision intvar(10)
	common /static/ intvar

	automatic j, dx, left, right

	include "GaussNDim.f"

	if(dim .eq. 0) then
	  result = func(intvar)
	else 
	  result = 0D0
	  do j = 1, points(dim)
	    dx = .5D0*gauss_x(dim, j)

*		may want to delete this in final version:
*	    if(dim .ge. 4) print *, "lev ", dim, " dx=", dx

	    intvar(maxdim + 1 - dim) = .5D0 + dx
	    call gauss(left, func, dim - 1, maxdim)
	    intvar(maxdim + 1 - dim) = .5D0 - dx
	    call gauss(right, func, dim - 1, maxdim)
	    result = result + gauss_w(dim, j)*(left + right)
	  enddo
	endif
	end

