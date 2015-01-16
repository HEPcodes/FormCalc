* run.F
* defines parameter settings for one run, then includes num.F
* this file is part of FormCalc
* last modified 19 Jun 01 th


* Whether to run in debugging mode

#define DEBUG

* Uncomment the following to check UV- and IR-finiteness

*#define DELTA 1D10
*#define MUDIM 1D100
*#define LAMBDA 1D10

* The model initialization code is included here. Options for the
* XXX_ini.F file (e.g. NO_SQUARK_MIXING for mssm_ini.F) must be
* defined before the #include line.

*#define NO_SQUARK_MIXING
*#define COMPLEX_PARAMETERS
*#define SM_ONLY
*#define NO_EXCLUSION_LIMITS

#include "sm_ini.F"

* You can define LOOP[1..10] for scanning the parameter space. LOOP1
* will be the outermost loop. Do not remove #defines which are not
* needed. Note that all loops terminate at label 1, e.g.
*	#define LOOP1 do 1 MH = 100, 500, 200
* For a fixed parameter, use something like #define LOOP1 MH = 300.

#define LOOP1 MH = 100
#define LOOP2
#define LOOP3
#define LOOP4
#define LOOP5
#define LOOP6
#define LOOP7
#define LOOP8
#define LOOP9
#define LOOP10
#define LOOP11
#define LOOP12
#define LOOP13
#define LOOP14
#define LOOP15

* PRINT[1..5] defines output statements with which changes in 
* parameters (typically the ones included in the scan) can be 
* written to the output file and screen. For example,
*	PRINT "MH=", MH, "MT=", MT
* writes out the values of MH and MT in each pass of the loop.
* Only real parameters are allowed (use dble and dimag for complex
* ones) and a string identifier must be given as in the example above.

#define PRINT1 PRINT "MH=", MH
#define PRINT2
#define PRINT3
#define PRINT4
#define PRINT5

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

