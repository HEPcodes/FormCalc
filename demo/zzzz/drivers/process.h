*	process.h
*	defines all process-dependent parameters for num.F
*	last modified 23 Feb 99 th

* whether to run in debugging mode:

#define DEBUG

* the types of the external particles:
* SCALAR, FERMION, PHOTON, or VECTOR (PHOTON is equivalent to VECTOR,
* except that longitudinal modes are not allowed)

#define TYPE1 VECTOR
#define TYPE2 VECTOR
#define TYPE3 VECTOR
#define TYPE4 VECTOR

* and their masses:

#define MASS1 MZ
#define MASS2 MZ
#define MASS3 MZ
#define MASS4 MZ

* the combinatorical factor for identical particles in the final state:
* 1/2 for identical particles, 1 otherwise

#define IDENTICALFACTOR .5D0

* possibly a color factor if there are quarks in the final state:

#define COLORFACTOR 1D0

* whether to include soft-photon bremsstrahlung:

*#define BREMSSTRAHLUNG
*#include "softphot.F"

* This following line can be used to initialize model parameters
* (or call a subroutine to do so)

#define MODELINI call sm_ini
#include "sm_ini.F"

