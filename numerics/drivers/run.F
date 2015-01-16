* run.F
* defines parameter settings for one run, then includes num.F
* this file is part of FormCalc
* last modified 5 Apr 00 th


* Whether to run in debugging mode

#define DEBUG

* Uncomment the following to check UV- and IR-finiteness

*#define DELTA 1D10
*#define MUDIM 1D100
*#define LAMBDA 1D10

* The variables in VAROUT are written out as comments into the
* data file using the format VARFMT

#define VAROUT MH
#define VARFMT "# MH=", F10.2

* You can define LOOP[1..10] for scanning the parameter space. LOOP1
* will be the outermost loop. Do not remove #defines which are not
* needed. Note that all loops terminate at label 1, e.g.
*	#define LOOP1 do 1 MH = 100, 500, 200
* For a fixed parameter, use something like #define LOOP1 MH = 300.

#define LOOP1 MT = 174
#define LOOP2 MH = 100
#define LOOP3
#define LOOP4
#define LOOP5
#define LOOP6
#define LOOP7
#define LOOP8
#define LOOP9
#define LOOP10

* The angular cutoff

#define CUTOFF (10*degree)

* The increment in angle and energy

#ifdef DEBUG
#define THETASTEP (125*degree)
#else
#define THETASTEP (1*degree)
#endif

#define ECMSSTEP 10

* The number of sampling points for the Gaussian integration

#define GAUSSPOINTS 16


* Now go for it:

#include "num.F"

