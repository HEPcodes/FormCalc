* const.h
* model-independent constants
* this file is part of FormCalc
* last modified 19 Oct 15 th


	RealType pi, degree, sqrt2, sqrt3, hbar_c2

	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (degree = pi/180D0)
	parameter (sqrt2 = 1.4142135623730950488016887242096981D0)
	parameter (sqrt3 = 1.7320508075688772935274463415058724D0)

	parameter (hbar_c2 = 3.8937966D8)
*         = hbar c^2 in picobarn

	ComplexType bogus, cI

	parameter (bogus = (-1D123, -2D123))
*	  some weird number likely to noticeably distort the final result;
*	  used for initializing arrays to check that all components
*	  have been calculated

	parameter (cI = (0D0, 1D0))

	RealType C_F, C_A
	parameter (C_F = 4/3D0, C_A = 3)

	RealType Divergence, mudim, lambda, muscale
	integer epsi
	common /renorm/ Divergence, mudim, lambda, muscale
	common /renorm/ epsi

	RealType rootsvalue, imode
	common /cuttools_para/ rootsvalue, imode

