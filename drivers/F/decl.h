* decl.h
* these declarations are included "everywhere"
* this file is part of FormCalc
* last modified 19 Nov 13 th


#ifndef DECL_H
#define DECL_H

* declarations for the whole file (e.g. preprocessor defs)

#include "distrib.h"
#include "types.h"
#include "user.h"
#include "extra.h"
#include "util.h"

#else

* declarations for every subroutine

#include "const.h"
#include "user.h"
#include "util.h"
#include "looptools.h"

#ifndef DRIVER
#include "renconst.h"
#endif

#endif

