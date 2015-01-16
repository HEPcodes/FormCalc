* process.h
* defines all process-dependent parameters
* this file is part of FormCalc
* last modified 17 Jan 03 th

* Definition of the external particles.
* The TYPEn may be one of SCALAR, FERMION, PHOTON, or VECTOR.
* (PHOTON is equivalent to VECTOR, except that longitudinal
* modes are not allowed)

#define TYPE1 VECTOR
#define MASS1 MW
#define CHARGE1 -1

#define TYPE2 VECTOR
#define MASS2 MW
#define CHARGE2 1

#define TYPE3 VECTOR
#define MASS3 MW
#define CHARGE3 -1

#define TYPE4 VECTOR
#define MASS4 MW
#define CHARGE4 1

* The combinatorial factor for identical particles in the final state:
* 1/n! for n identical particles, 1 otherwise

#define IDENTICALFACTOR 1

* Possibly a colour factor, e.g.
* - an additional averaging factor if any of the incoming particles
*   carry colour,
* - the overall colour factor resulting from the external particles
*   if that cannot computed by FormCalc (e.g. if the model has no
*   colour indices, as SM.mod).

#define COLOURFACTOR 1

* Whether to include soft-photon bremsstrahlung.
* ESOFTMAX is the maximum energy a soft photon may have and may be
* defined in terms of sqrtS, the CMS energy.

#define BREMSSTRAHLUNG
#define ESOFTMAX .1D0*sqrtS

* Possibly some wave-function renormalization
* (e.g. if calculating in the background-field method)

#define WF_RENORMALIZATION 4*dWFW1

* Include the kinematics-dependent part of the code

#include "2to2.F"

