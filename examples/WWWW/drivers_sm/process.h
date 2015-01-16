* process.h
* defines all process-dependent parameters for num.F
* this file is part of FormCalc
* last modified 19 Jun 01 th

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

* The combinatorical factor for identical particles in the final state:
* .5D0 for identical particles, 1 otherwise

#define IDENTICALFACTOR 1

* Possibly a colour factor if the external particles carry colour
* (SM.mod only, SMc.mod gets the factor right by its colour indices)

#define COLOURFACTOR 1

* Whether to include soft-photon bremsstrahlung

#define BREMSSTRAHLUNG
#define ESOFTMAX .1D0*Ecms

* Possibly some wave-function renormalization
* (if calculating in the background-field method)

*#define WF_RENORMALIZATION (nW*dWFW1 + nZ*dWFZ1)

