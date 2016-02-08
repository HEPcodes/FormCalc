#if 0
	decl.h
	these declarations are included "everywhere"
	this file is part of FormCalc
	last modified 6 Jan 16 th
#endif


#ifndef DECL_H
#define DECL_H

// declarations for the whole file (e.g. preprocessor defs)

#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <string.h>

#include "distrib.h"
#include "clooptools.h"
typedef int integer;
typedef const int cinteger;
typedef long long integer8;
typedef const integer8 cinteger8;

#if NOUNDERSCORE
#define RenConst renconst
#else
#define RenConst renconst_
#endif

#include "extra.h"
#include "const.h"

#else

// declarations for every subroutine

#endif

#include "user.h"
#include "util.h"

#ifdef SQUAREDME
#include "RenConst.h"
#endif

