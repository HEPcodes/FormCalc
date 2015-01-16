* rc_zero.F
* sets all renormalization constants to zero for debugging purposes
* this file is part of FormCalc
* last modified 23 Feb 00 th


	subroutine calc_renconst
	implicit none

#include "renconst.h"

	integer g1

	dMZsq1 = 0
	dMWsq1 = 0
	dMHsq1 = 0
	dSWsq1 = 0
	dCWsq1 = 0
	dSW1 = 0
	dZAA1 = 0
	dZAZ1 = 0
	dZZA1 = 0
	dZZZ1 = 0
	dZW1 = 0
	dZH1 = 0
	dZphi1 = 0
	dZchi1 = 0
	dWFZ1 = 0
	dWFW1 = 0
	dZe1 = 0
	dTad1 = 0

	do g1 = 1, 3
	  dMf1(2,g1) = 0
	  dZfL1(2,g1,g1) = 0
	  dZfR1(2,g1,g1) = 0
	  dMf1(3,g1) = 0
	  dZfL1(3,g1,g1) = 0
	  dZfR1(3,g1,g1) = 0
	  dMf1(4,g1) = 0
	  dZfL1(4,g1,g1) = 0
	  dZfR1(4,g1,g1) = 0
	enddo

	end

