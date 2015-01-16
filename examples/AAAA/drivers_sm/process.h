* process.h
* defines all process-dependent parameters for num.F
* this file is part of FormCalc
* last modified 28 Aug 01 th

* Definition of the external particles.
* The TYPEn may be one of SCALAR, FERMION, PHOTON, or VECTOR.
* (PHOTON is equivalent to VECTOR, except that longitudinal
* modes are not allowed)

#define TYPE1 PHOTON
#define MASS1 0
#define CHARGE1 0

#define TYPE2 PHOTON
#define MASS2 0
#define CHARGE2 0

#define TYPE3 PHOTON
#define MASS3 0
#define CHARGE3 0

#define TYPE4 PHOTON
#define MASS4 0
#define CHARGE4 0

* The combinatorical factor for identical particles in the final state:
* .5D0 for identical particles, 1 otherwise

#define IDENTICALFACTOR .5D0

* Possibly a colour factor if the external particles carry colour
* (SM.mod only, SMc.mod gets the factor right by its colour indices)

#define COLOURFACTOR 1

* Whether to include soft-photon bremsstrahlung
* ESOFTMAX is the maximum energy a soft photon may have and may be
* defined in terms of Ecms, the CMS energy.

*#define BREMSSTRAHLUNG
#define ESOFTMAX .1D0*Ecms

* Possibly some wave-function renormalization
* (if calculating in the background-field method)

*#define WF_RENORMALIZATION (nW*dWFW1 + nZ*dWFZ1)

* Include the kinematics-dependent part of the code

#include "2to2.F"

