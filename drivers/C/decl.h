#if 0
	decl.h
	these declarations are included "everywhere"
	this file is part of FormCalc
	last modified 18 Mar 13 th
#endif


#ifndef DECL_H
#define DECL_H

// declarations for the whole file (e.g. preprocessor defs)

#include <math.h>
#include <complex.h>
#include <stdio.h>

#include "clooptools.h"
typedef ComplexType HelType;
typedef int integer;
typedef const int cint;

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

