* CTBcut.F
* the FormCalc wrapper for the CutTools two-point function
* this file is part of FormCalc
* last modified 8 Dec 10 th


	double complex function CTBcut(rank, num, p, m1, m2)
	implicit none
	integer rank
	external num
	integer*8 p
	double precision m1, m2

	double precision pj(0:3,2), mj(2)
	external CTmpnum

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0
	call CTMomDecode(p, pj(0,2))

	mj(1) = m1
	mj(2) = m2

	PRINT *, "CTBcut:", rank
	PRINT *, "p1=", pj(0,2)
	PRINT *, "   ", pj(1,2)
	PRINT *, "   ", pj(2,2)
	PRINT *, "   ", pj(3,2)
	PRINT *, "m1=", m1
	PRINT *, "m2=", m2

	call ctsxcut(2, num, CTmpnum, rank, pj, mj, CTBcut)
c	call ctsxcut(2, num, CTmpnum, 2, pj, mj, CTBcut)
	PRINT *, "CTBcut=", CTBcut
	end

