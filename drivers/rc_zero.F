* rc_zero.F
* sets all renormalization constants to zero for debugging purposes
* this file is part of FormCalc
* last modified 14 Jan 03 th

#include "prefix.h"

	subroutine CalcRenConst
	implicit none

#include "renconst.h"

	data rc /sizeof_rc*0/

	end

