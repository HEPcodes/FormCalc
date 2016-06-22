* process.h
* defines all process-dependent parameters
* this file is part of FormCalc
* last modified 8 Jun 16 th


* Possibly a colour factor, e.g.
* - an additional averaging factor if any of the incoming particles
*   carry colour (e.g. 1/8D0 for an incoming gluon),
* - the overall colour factor resulting from the external particles
*   if that cannot computed by FormCalc (e.g. if the model has no
*   colour indices, as SMew.mod).
* You MUST adapt partonic.h if not all subprocesses share this factor.

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

#ifndef PHOTONRADIATION
c#define PHOTONRADIATION SOFT
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

* For lumi_parton.F: whether to force the decaying particle to
* be on-shell, independent of the command-line choices for sqrtS;
* the value specifies the maximum value of |sqrtS - sum_masses_in|

#ifndef FORCE_ONSHELL
c#define FORCE_ONSHELL 1D-9
#endif

* for lumi_hadron.F:

#ifndef PDFSET
#define PDFSET "cteq5l.LHgrid"
#endif
#ifndef PDFMEM
#define PDFMEM 0
#endif


* At which stage to join the partonic XS:
* 0 = add the differential partonic XS, then integrate sum,
* 1 = integrate each differential partonic XS, then add up.

#ifndef JOIN_PARTONIC
#define JOIN_PARTONIC 0
#endif

