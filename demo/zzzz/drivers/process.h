* process.h
* defines all process-dependent parameters for num.F
* this file is part of FormCalc
* last modified 30 Mar 00 th

* Definition of the external particles.
* The TYPEn may be one of SCALAR, FERMION, PHOTON, or VECTOR.
* (PHOTON is equivalent to VECTOR, except that longitudinal
* modes are not allowed)

#define TYPE1 VECTOR
#define MASS1 MZ
#define CHARGE1 0

#define TYPE2 VECTOR
#define MASS2 MZ
#define CHARGE2 0

#define TYPE3 VECTOR
#define MASS3 MZ
#define CHARGE3 0

#define TYPE4 VECTOR
#define MASS4 MZ
#define CHARGE4 0

* The combinatorical factor for identical particles in the final state:
* .5D0 for identical particles, 1 otherwise

#define IDENTICALFACTOR .5D0

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

