* decl.h
* these declarations are included "everywhere"
* this file is part of FormCalc
* last modified 10 Feb 17 th


#ifndef DECL_H
#define DECL_H

* declarations for the whole file (e.g. preprocessor defs)

#include "distrib.h"
#include "types.h"
#include "extra.h"
#include "generic.h"

#else

* declarations for every subroutine

#include "const.h"
#include "looptools.h"

#endif

#include "user.h"
#include "util.h"

#include "RenConst.h.F"
#include "MassShift.h.F"

