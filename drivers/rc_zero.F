* rc_zero.F
* sets all renormalization constants to zero for debugging purposes
* this file is part of FormCalc
* last modified 21 Jun 01 th

#include "prefix.h"

	subroutine calc_renconst
	implicit none

#include "renconst.h"

	data rc /sizeof_rc*0/

	end

