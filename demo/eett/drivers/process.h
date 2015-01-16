*	process.h
*	defines all process-dependent parameters for num.F
*	last modified 20 Feb 99 th

* The types of the external particles:
* SCALAR, FERMION, PHOTON, or VECTOR (PHOTON is equivalent to VECTOR,
* except that longitudinal modes are not allowed)

#define TYPE1 FERMION
#define TYPE2 FERMION
#define TYPE3 FERMION
#define TYPE4 FERMION

* and their masses

#define MASS1 ME
#define MASS2 ME
#define MASS3 MT
#define MASS4 MT

* The combinatorical factor for identical particles in the final state:
* 1/2 for identical particles, 1 otherwise

#define IDENTICALFACTOR 1D0

* Possibly a color factor if there are quarks in the final state

#define COLORFACTOR 3D0

* Whether to include soft-photon bremsstrahlung

*#define BREMSSTRAHLUNG
*#include "softphot.m"

* This following line can be used to initialize model parameters
* (or call a subroutine to do so)

#define MODELINI call sm_ini
#include "sm_ini.F"

