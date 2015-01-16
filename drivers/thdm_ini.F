* thdm_ini.F
* all parameters for the THDM are initialized here
* this file is part of FormCalc
* Abdesslam Arhrib, 2 Jun 03
* last modified 7 Feb 05 th


* The input parameters for thdm_ini are:

* SA = sin alpha
* TB = tan beta
* Mh0 = Mass(H1)
* MHH = Mass(H2)
* MA0 = Mass(H3)
* MHp = Mass(Hp)
* Lambda5

* The following preprocessor variables can be defined
* (they are usually set in run.F):

* MODEL_TYPE_I, MODEL_TYPE_II
*   Which type of two-higgs-doublet model to select.

* NO_EXCLUSION_LIMITS
*   If defined, the experimental exclusion limits are ignored.


* Exclusion limits (only checked if NO_EXCLUSION_LIMITS is not set):

* PDG
#define DeltaRho_UPPERBOUND .003D0


#define ModelDefaults SMDefaults
#define ModelConstIni SMConstIni
#include "sm_ini.F"
#undef ModelDefaults
#undef ModelConstIni


************************************************************************

	subroutine ModelDefaults
	implicit none

#include "model.h"

	call SMDefaults
	end


************************************************************************

	subroutine ModelConstIni(*)
	implicit none

#include "model.h"

	double precision Lambda1, Lambda2, Lambda3, Lambda4, Lambda6
	double precision DeltaRho

#if !defined(MODEL_TYPE_I) && !defined(MODEL_TYPE_II)
	stop "Neither MODEL_TYPE_I nor MODEL_TYPE_II defined"
#elif defined(MODEL_TYPE_I) && defined(MODEL_TYPE_II)
	stop "Both MODEL_TYPE_I and MODEL_TYPE_II defined"
#endif

	call SMConstIni(*999)

	TB2 = TB**2
	SB = TB/sqrt(1 + TB2)
	CB = SB/TB
	SB2 = SB**2
	CB2 = CB**2
	S2B = 2*SB*CB
	C2B = CB2 - SB2

	SA2 = SA**2
	CA2 = 1 - SA2
	CA = sqrt(CA2)
	S2A = 2*SA*CA
	C2A = CA2 - SA2
	SAB = SA*CB + CA*SB
	CAB = CA*CB - SA*SB
	SBA = SB*CA - CB*SA
	CBA = CB*CA + SB*SA

#ifdef MODEL_TYPE_I
	Yuk1 = CA/SB
	Yuk2 = SA/SB
	Yuk3 = -1/TB
#else
	Yuk1 = -SA/CB
	Yuk2 = CA/CB
	Yuk3 = TB
#endif

	Mh02 = Mh0**2
	MH = Mh0
	MH2 = Mh02
	MHH2 = MHH**2
	MA02 = MA0**2
	MHp2 = MHp**2

#ifndef NO_EXCLUSION_LIMITS
	Lambda1 = .25D0*Lambda5*(1 - TB2) +
     &    Alfa*Pi/(2*MW2*SW2)/CB2*
     &      (MHH2*CA2 + Mh02*SA2 - .5D0*(MHH2 - Mh02)*S2A/TB)
	Lambda2 = .25D0*Lambda5*(1 - 1/TB2) +
     &    Alfa*Pi/(2*MW2*SW2)/SB2*
     &      (Mh02*CA2 + MHH2*SA2 - .5D0*(MHH2 - Mh02)*S2A*TB)
        Lambda3 = -.25D0*Lambda5 +
     &    Alfa*Pi/(2*MW2*SW2)*(MHH2 - Mh02)*S2A/S2B
        Lambda4 = (2*Alfa*Pi)/(MW2*SW2)*MHp2
        Lambda6 = (2*Alfa*Pi)/(MW2*SW2)*MA02

	if( max(abs(Lambda1), abs(Lambda2), abs(Lambda3),
     &          Lambda4, abs(Lambda5), Lambda6) .gt. 8*pi ) then
	  print *, "thdm_ini: Perturbativity violated"
	  return 1
	endif

	if( abs(MA02 - MHH2) .gt. 1D-13 .and.
     &      abs(MHp2 - MHH2) .gt. 1D-13 ) then
	  DeltaRho = GF/(8*sqrt2*pi**2)*(
     &      MHp2*(1 - MA02/(MHp2 - MA02)*log(MHp2/MA02)) +
     &      CBA**2*Mh02*(MA02/(MA02 - Mh02)*log(MA02/Mh02) -
     &                   MHp2/(MHp2 - Mh02)*log(MHp2/Mh02)) +
     &      SBA**2*MHH2*(MA02/(MA02 - MHH2)*log(MA02/MHH2) -
     &                   MHp2/(MHp2 - MHH2)*log(MHp2/MHH2)) )

	  if( DeltaRho .gt. DeltaRho_UPPERBOUND ) then
	    print *, "thdm_ini: DeltaRho = ", DeltaRho,
     &        " > DeltaRho_UPPERBOUND = ", DeltaRho_UPPERBOUND
	    return 1
	  endif
	endif
#endif

	if( thdm_digest ) return
	thdm_digest = .TRUE.

	print *, "THDM parameters"
	print *, "----------------------------"
#ifdef MODEL_TYPE_I
	print *, "using type-I model"
#else
	print *, "using type-II model"
#endif
	print *, "TB =", TB
	print *, "SA =", SA
	print *, "Mh0 =", Mh0
	print *, "MHH =", MHH
	print *, "MA0 =", MA0
	print *, "MHp =", MHp
	print *, "Lambda5 =", Lambda5
	print *, "============================"
	return

999	return 1
	end

