* run.F
* defines parameter settings for one run, then includes num.F
* this file is part of FormCalc
* last modified 16 Dec 04 th


* Whether to run in debugging mode

#define DEBUG


* Uncomment the following to check UV- and IR-finiteness

c#define DELTA 1D7
c#define MUDIM 1D100
c#define LAMBDA 1D10

* The following sets all ren. constants to zero for debugging

c#define CalcRenConst ZeroRenConst


* The model initialization code is included here.  Options for the
* XXX_ini.F file (e.g. NO_SQUARK_MIXING for mssm_ini.F) must be
* defined before the #include line.

* options for sm_ini.F:
c#define NO_CKM

* options for mssm_ini.F:
c#define NO_SQUARK_MIXING
c#define HIGGS_MASSES TREE
c#define HIGGS_MASSES SIMPLE
c#define HIGGS_MASSES FEYNHIGGS
c#define SM_HIGGS
c#define NO_EXCLUSION_LIMITS

* options for thdm_ini.F:
c#define MODEL_TYPE_I
c#define MODEL_TYPE_II
c#define NO_EXCLUSION_LIMITS

#include "sm_ini.F"


* You can define LOOP[1..15] for scanning the parameter space.
* LOOP1 will be the outermost loop.  Do not remove #defines which
* are not needed.  Note that all loops terminate at label 1, e.g.
*	#define LOOP1 do 1 MH = 100, 500, 200
* For a fixed parameter, use something like #define LOOP1 MH = 300.

#define LOOP1 MH = 100
#define LOOP2
#define LOOP3 TB = 5
#define LOOP4 MA0 = 500
#define LOOP5 MSusy = 1000
#define LOOP6 MUE = 200
#define LOOP7 M_2 = 200
#define LOOP8 MGl = .8D0*MSusy
#define LOOP9 At = MUE/TB + 2*MSusy
#define LOOP10 Ab = At
#define LOOP11 Atau = At
#define LOOP12
#define LOOP13
#define LOOP14
#define LOOP15


* PRINT[1..5] defines output statements with which changes in 
* parameters (typically the ones included in the scan) can be 
* written to the output file and screen.  For example,
*	PRINT "MH=", MH, "MT=", MT
* writes out the values of MH and MT in each pass of the loop.
* Only real parameters are allowed (use dble and dimag for complex
* ones) and a string identifier must be given as in the example above.

#define PRINT1 PRINT "MH=", MH
#define PRINT2
#define PRINT3
#define PRINT4
#define PRINT5


* Cuts and increments.  May be defined in terms of sqrtS, the CMS energy.

#define THETACUT (10*degree)
#ifdef DEBUG
#define THETASTEP (125*degree)
#else
#define THETASTEP (1*degree)
#endif

#define K50CUT ESOFTMAX
#define K50STEP .05D0*sqrtS

#define SQRTSSTEP 10


* The integration parameters:
* - the desired relative and absolute accuracy

#define RELACCURACY 1D-3
#define ABSACCURACY 1D-10

* The multidimensional integration routines (2 -> 3) also need

#define METHOD DIVONNE
#define VERBOSE 1
#define MINEVAL 0
#define MAXEVAL 50000

* for Vegas:
#define NSTART 1000
#define NINCREASE 500

* for Suave:
#define NNEW 1000
#define FLATNESS 50

* for Divonne:
#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 1D-6
#define MAXCHISQ 10
#define MINDEVIATION .25D0

* for Cuhre:
#define KEY 0

* Now go for it:

#include "main.F"

