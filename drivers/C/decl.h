#if 0
	decl.h
	these declarations are included "everywhere"
	this file is part of FormCalc
	last modified 23 Apr 13 th
#endif


#ifndef DECL_H
#define DECL_H

// declarations for the whole file (e.g. preprocessor defs)

#include <math.h>
#include <complex.h>
#include <stdio.h>

#include "distrib.h"
#include "clooptools.h"
typedef int integer;
typedef const int cinteger;
typedef long long integer8;
typedef const integer8 cinteger8;

#if NOUNDERSCORE
#define CalcRenConst calcrenconst
#else
#define CalcRenConst calcrenconst_
#endif

#include "user.h"
#include "const.h"
#include "util.h"
#include "renconst.h"

#else

// declarations for every subroutine

#include "user.h"

#endif

