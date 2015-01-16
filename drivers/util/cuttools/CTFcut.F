* CTFcut.F
* the FormCalc wrapper for the CutTools six-point function
* this file is part of FormCalc
* last modified 8 Dec 10 th


	double complex function CTFcut(rank, num,
     &    p1, p2, p3, p4, p5, m1, m2, m3, m4, m5, m6)
	implicit none
	integer rank
	external num
	integer*8 p1, p2, p3, p4, p5
	double precision m1, m2, m3, m4, m5, m6

	double precision pj(0:3,6), mj(6)
	external CTmpnum

	pj(0,1) = 0
	pj(1,1) = 0
	pj(2,1) = 0
	pj(3,1) = 0
	call CTMomDecode(p1, pj(0,2))
	call CTMomDecode(p2, pj(0,3))
	call CTMomDecode(p3, pj(0,4))
	call CTMomDecode(p4, pj(0,5))
	call CTMomDecode(p5, pj(0,6))

	mj(1) = m1
	mj(2) = m2
	mj(3) = m3
	mj(4) = m4
	mj(5) = m5
	mj(6) = m6

	call ctsxcut(6, num, CTmpnum, rank, pj, mj, CTFcut)
c	call ctsxcut(6, num, CTmpnum, 6, pj, mj, CTFcut)
	end

