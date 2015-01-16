* process.h
* defines all process-dependent parameters
* this file is part of FormCalc
* last modified 19 Nov 13 th


* When using Dirac fermions (FermionChains -> Chiral|VA) and
* the trace technique (HelicityME), the following flag should be 
* defined to compute unpolarized cross-sections efficiently,
* i.e. without actually summing up the different helicities.
* This has no effect on the result, only on the speed of the
* calculation.
* Note: DIRACFERMIONS must NOT be defined when using Weyl fermions,
* i.e. FermionChains -> Weyl in CalcFeynAmp.

c#define DIRACFERMIONS


* The combinatorial factor for identical particles in the final state:
* 1/n! for n identical particles, 1 otherwise

#define IDENTICALFACTOR 1


* Possibly a colour factor, e.g.
* - an additional averaging factor if any of the incoming particles
*   carry colour,
* - the overall colour factor resulting from the external particles
*   if that cannot computed by FormCalc (e.g. if the model has no
*   colour indices, as SMew.mod).

#define COLOURFACTOR 1


* The scale at which the interaction takes place
* (= the factorization scale for an hadronic process).

#define FSCALE sqrtS


* Whether to include soft-photon bremsstrahlung.
* ESOFTMAX is the maximum energy a soft photon may have and may be
* defined in terms of sqrtS, the CMS energy.

c#define BREMSSTRAHLUNG
#define ESOFTMAX .1D0*sqrtS


* NCOMP is the number of components of the result vector.  Currently
* the components are 1 = tree-level result, 2 = one-loop result.

#define NCOMP 2


* Choose the appropriate luminosity for the collider:
* - lumi_parton.F for a "parton collider" (e.g. e+ e- -> X),
* - lumi_hadron.F for a hadron collider (e.g. p pbar -> X),
* - lumi_photon.F for a photon collider (gamma gamma -> X)

#define LUMI "lumi_parton.F"

* for lumi_parton.F: whether to force the decaying particle to
* be on-shell, independent of the command-line choices for sqrtS;
* the value specifies the maximum value of |sqrtS - sum_masses_in|

c#define FORCE_ONSHELL 1D-9

* for lumi_hadron.F: PARTON1 and PARTON2 identify the
* incoming partons by their PDG code, where
* 0 = gluon
* 1 = down   3 = strange   5 = bottom
* 2 = up     4 = charm     6 = top

#define PARTON1 1
#define PARTON2 1
#define PDFSET "cteq5l.LHgrid"
#define PDFMEM 0


* Include process specifications

#include "squaredme/specs.h"

