*	run.F
*	defines parameter settings for one run, then includes num.F
*	last modified 24 Aug 99 th

* Whether to run in debugging mode

#define DEBUG

* The variables in VAROUT are written out as comments into the
* data file using the format VARFMT

#define VAROUT MH
#define VARFMT "# MH=", F10.2

* You can define LOOP[1..10] for scanning the parameter space. LOOP1
* will be the outermost loop. Do not remove #defines which are not
* needed. Note that all loops terminate at label 1, e.g.
*	#define LOOP1 do 1 MH = 100D0, 500D0, 200D0
* For a fixed parameter, use something like #define LOOP1 MH = 300D0.

#define LOOP1 MT = 175D0
#define LOOP2 MT2 = MT**2
#define LOOP3 MH = 100D0
#define LOOP4 MH2 = MH**2
#define LOOP5
#define LOOP6
#define LOOP7
#define LOOP8
#define LOOP9
#define LOOP10

* The angular cutoff

#define CUTOFF (10d0*degree)

* The number of sampling points for the Gaussian integration

#define GAUSSPOINTS 16


* Now go for it:

#include "num.F"

