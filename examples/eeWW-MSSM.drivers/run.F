* run.F
* defines parameter settings for one run, then includes num.F
* this file is part of FormCalc
* last modified 8 Apr 04 th


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

* options for mssm_ini.F:
c#define NO_SQUARK_MIXING
c#define COMPLEX_PARAMETERS
c#define SM_ONLY
c#define NO_EXCLUSION_LIMITS
c#define TREE_LEVEL_MH

* options for thdm_ini.F:
c#define MODEL_TYPE_I
c#define MODEL_TYPE_II
c#define NO_EXCLUSION_LIMITS

#include "mssm_ini.F"


* You can define LOOP[1..15] for scanning the parameter space.
* LOOP1 will be the outermost loop.  Do not remove #defines which
* are not needed.  Note that all loops terminate at label 1, e.g.
*	#define LOOP1 do 1 MH = 100, 500, 200
* For a fixed parameter, use something like #define LOOP1 MH = 300.

#define LOOP1 TB = 1.5D0
#define LOOP2 MA0 = 700
#define LOOP3 MSusy = 1000
#define LOOP4 MUE = -1000
#define LOOP5 M_2 = 100
#define LOOP6 At = MSusy
#define LOOP7 Ab = At
#define LOOP8 Atau = At
#define LOOP9
#define LOOP10
#define LOOP11
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

#define PRINT1 PRINT "TB=", TB, "MA0=", MA0
#define PRINT2 PRINT "M_2=", M_2, "MUE=", reMUE
#define PRINT3 PRINT "MSusy=", MSusy, "At=", reAt
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
#define ABSACCURACY 1D-7

* The multidimensional integration routines (2 -> 3) also need

#define METHOD DIVONNE
#define FLAGS 1
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

