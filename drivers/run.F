* run.F
* defines parameter settings for one run, then includes num.F
* this file is part of FormCalc
* last modified 10 Oct 01 th


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
*#define TREE_LEVEL_MH

#include "sm_ini.F"

* You can define LOOP[1..15] for scanning the parameter space. LOOP1
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

* Cuts and increments. May be defined in terms of Ecms, the CMS energy.

#define THETACUT (10*degree)
#ifdef DEBUG
#define THETASTEP (125*degree)
#else
#define THETASTEP (1*degree)
#endif

#define K50CUT ESOFTMAX
#define K50STEP .05D0*Ecms

#define ECMSSTEP 10

* The integration parameters:
* For 2 -> 2: the number of nodes for the Gaussian quadrature.

#define GAUSSPOINTS 16

* For 2 -> 3: the desired accuracy,
* when using the adaptive Gauss-Kronrod quadrature:
* - the number of nodes used in each dimension (must be even
*   integers between 4 and 12, e.g. 4 = 4/9 Gauss-Kronrod rule),
* - the maximum number of subintervals the integration interval
*   may be split into, and
* - the gain in accuracy from one dimension to the next, i.e. the
*   integral in dimension (n - 1) is demanded to be GAIN times
*   more precise than the integral in dimension n,
* when using the Monte-Carlo integration:
* - the VEGAS switch,
* - the number of calls per iteration, and
* - the maximum number of iterations for VEGAS,
* when using the DCUHRE algorithm (link with -ldcuhre):
* - the DCUHRE switch,
* - the maximum number of integrand evaluations is NCALL*MAXITER,
*   where NCALL and MAXITER are the VEGAS parameters.

#define ACCURACY 1D-3

#define NODES_K50 8
#define NODES_COSTH 4
#define NODES_K30 6
#define NODES_ETA 4
#define MAXINTV 25
#define GAIN 1

*#define VEGAS
*#define DCUHRE
#define NCALL 1000
#define MAXITER 50

* Now go for it:

#include "num.F"

