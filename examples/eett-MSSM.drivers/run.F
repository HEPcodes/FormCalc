* run.F
* defines parameter settings for one run, then includes main.F
* this file is part of FormCalc
* last modified 28 Mar 22 th


* Whether to run in testing mode

#ifndef TESTRUN
#define TESTRUN 1
#endif

* Uncomment the following to check UV- and IR-finiteness

#ifndef DELTA
c#define DELTA 1D7
#endif

#ifndef MUDIM
c#define MUDIM 1D100
#endif

#ifndef LAMBDA
c#define LAMBDA 1D10
#endif


* Options for the model.F file (e.g. NO_SQUARK_MIXING for MSSM.F).

* options for SM.F:
c#define NO_CKM
#define MB_SCHEME ONSHELL
#define ALPHAS_ORDER 3

* options for MSSM.F:
c#define NO_SQUARK_MIXING
c#define HIGGS_MASSES TREE
c#define HIGGS_MASSES SIMPLE
c#define HIGGS_MASSES FEYNHIGGS
c#define SM_HIGGS
c#define FV 3
c#define NO_EXCLUSION_LIMITS
c#define NO_GUT_RELATION

* options for THDM.F:
c#define MODEL_TYPE_I
c#define MODEL_TYPE_II
c#define NO_EXCLUSION_LIMITS


* The LOOPn set up the model parameters for the calculation, either
* - as fixed declarations, e.g.
*       TB = 12.3D0
* - or as do-loops terminating on label 1, e.g.
*       do 1 MA0 = 100, 500, 10

#define LOOP1 TB = 5
#define LOOP2 MA0 = 750
#define LOOP3 MSusy = 1000
#define LOOP4 MUE = 200
#define LOOP5 Mino2 = 200
#define LOOP6 MGl = .8D0*MSusy
#define LOOP7 At = MUE/TB + 2*MSusy
#define LOOP8 Ab = At
#define LOOP9 Atau = At
#define LOOP10
#define LOOP11
#define LOOP12
#define LOOP13
#define LOOP14
#define LOOP15
#define LOOP16
#define LOOP17
#define LOOP18
#define LOOP19
#define LOOP20


* The PRINTn definitions contain output statements with which changes
* in parameters (typically the ones included in the scan) can be
* written to the output file and screen.  For example,
*       SHOW "MH", MH
* writes out the value of MH in each pass of the loop.

#define PRINT1 SHOW "TB", TB
#define PRINT2 SHOW "MA0", MA0
#define PRINT3 SHOW "M_2", ReM_2
#define PRINT4 SHOW "MUE", ReMUE
#define PRINT5 SHOW "MSusy", MSusy
#define PRINT6 SHOW "At", ReAt
#define PRINT7
#define PRINT8
#define PRINT9
#define PRINT10
#define PRINT11
#define PRINT12
#define PRINT13
#define PRINT14
#define PRINT15
#define PRINT16
#define PRINT17
#define PRINT18
#define PRINT19
#define PRINT20


* Cuts and increments.  May be defined in terms of sqrtS, the CMS energy.
* Note: cuts are numbered by the integration they belong to, which is
* from the last particle downward, e.g. COSTH1CUT refers to cos(th) of
* the *last* particle.

#define SQRTSSTEP 10

#define COSTH1CUT .01D0
#ifdef TESTRUN
#define XCOSTH1STEP 1/3D0
#else
#define XCOSTH1STEP .02D0
#endif

#define XMREM1STEP .1D0
#define E1MIN ESOFTMAX

* The following are veto cuts, i.e. the cross-section is kept if the
* cut conditions apply (set to zero otherwise).  They are formulated
* in terms of the CUT_... macros in cuts.F, e.g.
*   #define CUT1 CUT_kT(3) > 10
* The CUT1..20 together make up a single logical expression in Fortran
* and may include logical operators, e.g.
*   #define CUT1 CUT_kT(3) > 10
*   #define CUT2 .and. CUT_kT(4) > 10

#define CUT1 .TRUE.
#define CUT2
#define CUT3
#define CUT4
#define CUT5
#define CUT6
#define CUT7
#define CUT8
#define CUT9
#define CUT10
#define CUT11
#define CUT12
#define CUT13
#define CUT14
#define CUT15
#define CUT16
#define CUT17
#define CUT18
#define CUT19
#define CUT20


* The integration parameters:
* - the desired relative and absolute accuracy

#define RELACCURACY 1D-3
#define ABSACCURACY 1D-10

* The multidimensional integration routines also need

#define METHOD DIVONNE
#define USERDATA bogus
#define NVEC 1
#define VERBOSE 1
#define SEED 0
#define MINEVAL 1000
#define MAXEVAL 100000
#define STATEFILE ""
#define SPIN -1

* for Vegas:
#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0

* for Suave:
#define NNEW 1000
#define NMIN 2
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

