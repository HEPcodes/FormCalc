* decl.h
* these declarations are included "everywhere"
* this file is part of FormCalc
* last modified 23 Dec 15 th


#ifndef DECL_H
#define DECL_H

* declarations for the whole file (e.g. preprocessor defs)

#include "distrib.h"
#include "types.h"
#include "extra.h"

#else

* declarations for every subroutine

#include "const.h"
#include "looptools.h"

#endif

#include "user.h"
#include "util.h"

#ifdef SQUAREDME
#include "RenConst.h"
#endif

