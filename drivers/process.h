* process.h
* defines all process-dependent parameters
* this file is part of FormCalc
* last modified 29 Apr 15 th


* When using Dirac fermions (FermionChains -> Chiral|VA) and
* the trace technique (HelicityME), the following flag should be 
* defined to compute unpolarized cross-sections efficiently,
* i.e. without actually summing up the different helicities.
* This has no effect on the result, only on the speed of the
* calculation.
* Note: DIRACFERMIONS must NOT be defined when using Weyl fermions,
* i.e. FermionChains -> Weyl in CalcFeynAmp.

#ifndef DIRACFERMIONS
c#define DIRACFERMIONS
#endif


* The combinatorial factor for identical particles in the final state:
* 1/n! for n identical particles, 1 otherwise

#ifndef IDENTICALFACTOR
#define IDENTICALFACTOR 1
#endif


* Possibly a colour factor, e.g.
* - an additional averaging factor if any of the incoming particles
*   carry colour,
* - the overall colour factor resulting from the external particles
*   if that cannot computed by FormCalc (e.g. if the model has no
*   colour indices, as SMew.mod).

#ifndef COLOURFACTOR
#define COLOURFACTOR 1
#endif


* The scale at which the interaction takes place
* (= the factorization scale for an hadronic process).

#ifndef FSCALE
#define FSCALE sqrtS
#endif


* Whether to include soft-photon bremsstrahlung.
* ESOFTMAX is the maximum energy a soft photon may have and may be
* defined in terms of sqrtS, the CMS energy.

#ifndef BREMSSTRAHLUNG
c#define BREMSSTRAHLUNG
#endif
#ifndef ESOFTMAX
#define ESOFTMAX .1D0*sqrtS
#endif


* NCOMP is the number of components of the result vector.  Currently
* the components are 1 = tree-level result, 2 = one-loop result.

#ifndef NCOMP
#define NCOMP 2
#endif


* Choose the appropriate luminosity for the collider:
* - lumi_parton.F for a "parton collider" (e.g. e+ e- -> X),
* - lumi_hadron.F for a hadron collider (e.g. p pbar -> X),
* - lumi_photon.F for a photon collider (gamma gamma -> X)

#ifndef LUMI
#define LUMI "lumi_parton.F"
#endif

* for lumi_parton.F: whether to force the decaying particle to
* be on-shell, independent of the command-line choices for sqrtS;
* the value specifies the maximum value of |sqrtS - sum_masses_in|

#ifndef FORCE_ONSHELL
c#define FORCE_ONSHELL 1D-9
#endif

* for lumi_hadron.F: PARTON1 and PARTON2 identify the
* incoming partons by their PDG code, where
* 0 = gluon
* 1 = down   3 = strange   5 = bottom
* 2 = up     4 = charm     6 = top

#ifndef PARTON1
#define PARTON1 1
#endif
#ifndef PARTON2
#define PARTON2 1
#endif
#ifndef PDFSET
#define PDFSET "cteq5l.LHgrid"
#endif
#ifndef PDFMEM
#define PDFMEM 0
#endif


* Include process specifications

#include "squaredme/specs.h"

