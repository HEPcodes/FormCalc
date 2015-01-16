	program test
	implicit none
	double precision res
	double precision integrand
	external integrand

	call gauss(res, integrand, 3, 3)
	print *, "integral = ", res
	end

************************************************************************

	double precision function integrand(X)
	implicit none
	double precision X(10)

	integrand = sin(X(1))*cos(X(2))*exp(X(3))
	end

************************************************************************

	include "ndimgauss_scalar.f"

