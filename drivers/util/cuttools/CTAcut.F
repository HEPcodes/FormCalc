* CTAcut.F
* the FormCalc wrapper for the CutTools one-point function
* this file is part of FormCalc
* last modified 8 Dec 10 th


	double complex function CTAcut(rank, num, m)
	implicit none
	integer rank
	external num
	double precision m

	double precision pj(0:3,1)
	external CTmpnum

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0

	PRINT *, "CTAcut:", rank
	PRINT *, "m1=", m

	call ctsxcut(1, num, CTmpnum, rank, pj, m, CTAcut)
c	call ctsxcut(1, num, CTmpnum, 1, pj, m, CTAcut)
	PRINT *, "CTAcut=", CTAcut
	end

