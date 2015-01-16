* process.h
* defines all process-dependent parameters for num.F
* this file is part of FormCalc
* last modified 30 Mar 00 th

* Definition of the external particles.
* The TYPEn may be one of SCALAR, FERMION, PHOTON, or VECTOR.
* (PHOTON is equivalent to VECTOR, except that longitudinal
* modes are not allowed)

#define TYPE1 FERMION
#define MASS1 ME
#define CHARGE1 1

#define TYPE2 FERMION
#define MASS2 ME
#define CHARGE2 -1

#define TYPE3 FERMION
#define MASS3 MT
#define CHARGE3 -2/3D0

#define TYPE4 FERMION
#define MASS4 MT
#define CHARGE4 2/3D0

* The combinatorical factor for identical particles in the final state:
* .5D0 for identical particles, 1 otherwise

#define IDENTICALFACTOR 1

* Possibly a colour factor if there are quarks in the final state
* (SM.mod only, SMc.mod gets the factor right by its colour indices)

#define COLOURFACTOR 1

* Whether to include soft-photon bremsstrahlung

*#define BREMSSTRAHLUNG
#define ESOFTMAX .1D0*Ecms

* Possibly some wave-function renormalization

*#define WF_RENORMALIZATION (nW*dWFW1 + nZ*dWFZ1)

* Include the model initialization routine (model_ini)

#include "sm_ini.F"

