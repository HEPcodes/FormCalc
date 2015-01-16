* This file sets all renormalization constants to zero
* for debugging purposes

	subroutine calc_renconst
	implicit none
#include "rcsm.h"

	integer g1

	do g1 = 1, 3
	  dMf1(2,g1) = 0D0
	  dZfL1(2,g1,g1) = 0D0
	  dZfR1(2,g1,g1) = 0D0
	  dMf1(3,g1) = 0D0
	  dZfL1(3,g1,g1) = 0D0
	  dZfR1(3,g1,g1) = 0D0
	  dMf1(4,g1) = 0D0
	  dZfL1(4,g1,g1) = 0D0
	  dZfR1(4,g1,g1) = 0D0
	enddo

	dMZsq1 = 0D0
	dMWsq1 = 0D0
	dMHsq1 = 0D0
	dSWsq1 = 0D0
	dCWsq1 = 0D0
	dSW1 = 0D0
	dZAA1 = 0D0
	dZAZ1 = 0D0
	dZZA1 = 0D0
	dZZZ1 = 0D0
	dZW1 = 0D0
	dZH1 = 0D0
	dZphi1 = 0D0
	dZchi1 = 0D0
	dWFZ1 = 0D0
	dWFW1 = 0D0
	dZe1 = 0D0
	dTad1 = 0D0
	GammaHMH = 1D0
	end

