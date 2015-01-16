* user.h
* the model declarations should be inserted here
* as well as any kind of user definition
* this file is part of FormCalc
* last modified 4 Jul 14 th


#ifndef USER_H
#define USER_H
* declarations for the whole file (e.g. preprocessor defs)

#define NINJA
c#define SAMURAI
c#define CUTTOOLS

* Possibly some wave-function renormalization
* (e.g. if calculating in the background-field method)

c#define WF_RENORMALIZATION (nW*dWFW1 + nZ*dWFZ1)

#else
* declarations for every subroutine

c#include "opp.h"
#include "model_mssm.h"

#endif

