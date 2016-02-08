* partonic.h
* defines how to build the total cross-section from partonic ones
* this file is part of FormCalc
* last modified 20 Jan 16

* To add a partonic process, copy and adapt the part between
* "BEGIN PARTONIC PROCESS" and "END PARTONIC PROCESS":
*
* PARTON1 and PARTON2 identify the incoming partons in
* lumi_hadron.F by their PDG code, where
* 0 = gluon
* 1 = down   3 = strange   5 = bottom
* 2 = up     4 = charm     6 = top


* NPID is the number of partonic processes.

#ifndef NPID
#define NPID 1
#endif


*** BEGIN PARTONIC PROCESS

#define PID 1
#define PARTON1 -999
#define PARTON2 -999
#include "squaredme/specs.h"
#include "parton.h"

*** END PARTONIC PROCESS

