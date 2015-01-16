* CTEcut.F
* the FormCalc wrapper for the CutTools five-point function
* this file is part of FormCalc
* last modified 8 Dec 10 th


	double complex function CTEcut(rank, num,
     &    p1, p2, p3, p4, m1, m2, m3, m4, m5)
	implicit none
	integer rank
	external num
	integer*8 p1, p2, p3, p4
	double precision m1, m2, m3, m4, m5

	double precision pj(0:3,5), mj(5)
	external CTmpnum

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0
	call CTMomDecode(p1, pj(0,2))
	call CTMomDecode(p2, pj(0,3))
	call CTMomDecode(p3, pj(0,4))
	call CTMomDecode(p4, pj(0,5))

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4
	mj(5) = m5

	call ctsxcut(5, num, CTmpnum, rank, pj, mj, CTEcut)
c	call ctsxcut(5, num, CTmpnum, 5, pj, mj, CTEcut)
	end

