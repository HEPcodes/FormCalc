	program test
	implicit none
	double precision res
	double precision integrand
	external integrand

	call gauss(res, integrand, 3, 3, 1, 1)
	print *, "integral = ", res
	end

************************************************************************
	subroutine integrand(res, X)
	implicit none
	double precision res, X(10)

	res = sin(X(1))*cos(X(2))*exp(X(3))
	end

************************************************************************
	include "ndimgauss_vector.f"

