* mssm_ini.F
* all parameters for the MSSM are initialized here
* this file is part of FormCalc
* Christian Schappacher, 17 Feb 99
* last modified 14 Jan 03 th


* The input parameters for mssm_ini are:

* Higgs sector:
*   TB = tan beta
*   MA0 = Mass(A_0)

* Chargino/Neutralino sector:
*   M_2 
*   Mue

* Sfermion sector:
*   MSusy
*   Atau, At, Ab
* If one wants to specify the breaking and mixing in the
* sfermion sector more precisely, the following preprocessor
* variables can be used:
*   MSQ, MSU, MSD, MSL, MSE (= MSusy by default)


* The following preprocessor variables can be defined
* (they are usually set in process.h):

* NO_SQUARK_MIXING
*   If defined, sets Atau, At, and Ab such that the off-diagonal entries
*   of the sfermion mixing matrices vanish, i.e. gauge eigenstates =
*   mass eigenstates.

* COMPLEX_PARAMETERS
*   If defined, uses a simpler (one-loop) approximation of the
*   Higgs masses which is valid for all parameters.  Normally,
*   the more precise (two-loop) approximation is taken, which is
*   however valid only for real parameters.

* SM_ONLY
*   If defined, reverse-engineers the mixing in the Higgs sector
*   (alpha and beta) such that the MSSM Higgs sector looks like
*   a SM Higgs sector.

* NO_EXCLUSION_LIMITS
*   If defined, the experimental exclusion limits are ignored.

* TREE_LEVEL_MH
*   If defined, uses only the tree-level Higgs masses.


#ifndef MSQ
#define MSQ MSusy
#endif

#ifndef MSU
#define MSU MSusy
#endif

#ifndef MSL
#define MSL MSusy
#endif

#ifndef MSD
#define MSD MSusy
#endif

#ifndef MSE
#define MSE MSusy
#endif

#define MSQ2 MSQ**2
#define MSU2 MSU**2
#define MSD2 MSD**2
#define MSE2 MSE**2
#define MSL2 MSL**2


* Exclusion limits (only checked if NO_EXCLUSION_LIMITS is not set):

* hep-ex/9910020 (L3)
#define MStop_LOWERBOUND 80

* hep-ex/9910020 (L3)
#define MSbottom_LOWERBOUND 70

* hep-ex/9910020 (L3)
#define MSquark_LOWERBOUND 150

* hep-ex/9910006 (L3)
#define MSlepton_LOWERBOUND 70

* PDG
#define DeltaRho_UPPERBOUND .003D0

#ifdef COMPLEX_PARAMETERS
* Kane & Wang, Phys Lett B488 (2000) 383
#define Mh0_LOWERBOUND 85
#else
* http://lephiggs.web.cern.ch/LEPHIGGS/papers/osaka_note.ps
#define Mh0_LOWERBOUND 91
#endif

* hep-ex/9910007 (L3) and hep-ex/9909051 (OPAL)
#define MNeutralino_LOWERBOUND 30

* hep-ex/9910007 (L3) and hep-ex/9909051 (OPAL)
#define MChargino_LOWERBOUND 90

* hep-ex/9910020 (L3) and PDG
#define MGl_LOWERBOUND 175


#define ModelConstIni SMConstIni
#define ModelVarIni SMVarIni
#define ModelDigest SMDigest
#include "sm_ini.F"
#undef ModelConstIni
#undef ModelVarIni
#undef ModelDigest

#include "diag.F"


************************************************************************

	subroutine ModelConstIni(*)
	implicit none

#include "model.h"

	double precision alpha, beta

	double precision Alfas_MT, sig1, sig12, sig2
	common /radcor/ Alfas_MT, sig1, sig12, sig2

	double precision ALPHAS2
	external ALPHAS2

	call SMConstIni(*999)

	Alfas_MT = ALPHAS2(MT)

* GUT relation (Higgs Hunter's Guide p. 287)
	M_1 = 5/3D0*SW2/CW2*M_2

#ifdef NO_SQUARK_MIXING
	At = dconjg(MUE)/TB
	Ab = dconjg(MUE)*TB
	Atau = Ab
#endif

	TB2 = TB**2
	SB = TB/sqrt(1 + TB2)
	CB = SB/TB
	SB2 = SB**2
	CB2 = CB**2
	S2B = 2*SB*CB
	C2B = CB2 - SB2

	call SfermionMasses(*999)
	call CharginoMasses(*999)
	call NeutralinoMasses(*999)
	call HiggsMasses(*999)

#ifdef SM_ONLY

* With SM_ONLY defined, the values of the MSSM parameters are set
* such that the SM is reproduced.  See Higgs Hunter's Guide p. 356.
* Note one tricky thing: the Higgs mass Mh0 is first calculated as
* it would be in the MSSM (to get the same value when comparing
* MSSM <-> SM), but then beta is reverse-engineered to be consistent
* with the SM relations.  In other words, we start with the MSSM-
* value of beta, calculate Mh0, and then redefine beta to fulfill
* cos(2 beta) = Mh0/MZ.  In principle, beta could also be fixed such
* that cos(2 beta) = -Mh0/MZ, but this is insignificant here since
* we only need beta to be consistent with the values of Mh0 and MZ.

	C2B = Mh0/MZ
	beta = .5D0*acos(C2B)
	SB = sin(beta)
	CB = cos(beta)
	TB = SB/CB
	S2B = 2*SB*CB

	SA = -CB
	CA = SB
	SA2 = CB2
	CA2 = SB2
	S2A = -S2B
	C2A = -C2B
	SAB = -C2B
	CAB = S2B
	SBA = 1
	CBA = 0
#else
	alpha = (-(MA02 + MZ2)*SB*CB - sig12)/
     &        (MZ2*CB2 + MA02*SB2 - sig1 - Mh02)
	if( alpha .le. 0 ) then
	  alpha = atan(alpha)
	else
	  alpha = atan(alpha) - pi
	endif
	SA = sin(alpha)
	CA = cos(alpha)
	SA2 = SA**2
	CA2 = CA**2
	S2A = 2*SA*CA
	C2A = CA2 - SA2
	SAB = SA*CB + CA*SB
	CAB = CA*CB - SA*SB
	SBA = SB*CA - CB*SA
	CBA = CB*CA + SB*SA
#endif
	return

999	return 1
	end


*-----------------------------------------------------------------------
* Sfermion sector:

	subroutine SfermionMasses(*)
	implicit none

#include "model.h"

	integer j

	double precision Qe, Qu, Qd
	parameter (Qe = -1, Qu = 2/3D0, Qd = -1/3D0)

#ifndef NO_EXCLUSION_LIMITS
	double precision f1f2, abs2, DeltaRho
	double complex c_
	external f1f2
	abs2(c_) = dble(c_)**2 + dimag(c_)**2
#endif

	do j = 1, 3
* Sneutrino masses:
	  MSf2(1, 1, j) = MSL2 + .5D0*MZ2*C2B
	  if( MSf2(1, 1, j) .lt. 0 ) return 1
	  MSf(1, 1, j) = sqrt(MSf2(1, 1, j))
	  MSf(2, 1, j) = -1
	  MSf2(2, 1, j) = -1

* Slepton mass matrix, family j
	  Af(2, j) = Atau
	  call Eigen2x2(MSf(1, 2, j), MSf2(1, 2, j), USf(1, 1, 2, j),
     &      MSL2 + MZ2*(-.5D0 - Qe*SW2)*C2B + MLE2(j),
     &      MSE2 + Qe*MZ2*SW2*C2B + MLE2(j),
     &      MLE(j)*(dconjg(Af(2, j)) - MUE*TB), *999)

* Up-type squark mass matrix, family j
	  Af(3, j) = At
	  call Eigen2x2(MSf(1, 3, j), MSf2(1, 3, j), USf(1, 1, 3, j),
     &      MSQ2 + MZ2*(.5D0 - Qu*SW2)*C2B + MQU2(j),
     &      MSU2 + Qu*MZ2*SW2*C2B + MQU2(j),
     &      MQU(j)*(dconjg(Af(3, j)) - MUE/TB), *999)

* Down-type squark mass matrix, family j
	  Af(4, j) = Ab
	  call Eigen2x2(MSf(1, 4, j), MSf2(1, 4, j), USf(1, 1, 4, j),
     &      MSQ2 + MZ2*(-.5D0 - Qd*SW2)*C2B + MQD2(j),
     &      MSD2 + Qd*MZ2*SW2*C2B + MQD2(j),
     &      MQD(j)*(dconjg(Af(4, j)) - MUE*TB), *999)
	enddo

#ifndef NO_EXCLUSION_LIMITS
	if( MSf(1, 4, 3) .lt. MSbottom_LOWERBOUND ) then
	  print *, "mssm_ini: MSf(1, 4, 3) = ", MSf(1, 4, 3),
     &      " < MSbottom_LOWERBOUND = ", MSbottom_LOWERBOUND
	  return 1
	endif

	if( MSf(1, 3, 3) .lt. MStop_LOWERBOUND ) then
	  print *, "mssm_ini: MSf(1, 3, 3) = ", MSf(1, 3, 3),
     &      " < MStop_LOWERBOUND = ", MStop_LOWERBOUND
	  return 1
	endif

	if( MSf(1, 2, 1) .lt. MSlepton_LOWERBOUND ) then
	  print *, "mssm_ini: MSf(1, 2, 1) = ", MSf(1, 2, 1),
     &      " < MSlepton_LOWERBOUND = ", MSlepton_LOWERBOUND
	  return 1
	endif

	if( MSf(1, 3, 1) .lt. MSquark_LOWERBOUND ) then
	  print *, "mssm_ini: MSf(1, 3, 1) = ", MSf(1, 3, 1),
     &      " < MSquark_LOWERBOUND = ", MSquark_LOWERBOUND
	  return 1
	endif

	DeltaRho =
     &     USf(1,1,3,3)*dconjg(USf(1,2,3,3))
     &     * USf(2,2,3,3)*dconjg(USf(2,1,3,3)) 
     &     * f1f2(MSf2(1,3,3), MSf2(2,3,3))
     &   + USf(1,1,4,3)*dconjg(USf(1,2,4,3))
     &     * USf(2,2,4,3)*dconjg(USf(2,1,4,3)) 
     &     * f1f2(MSf2(1,4,3), MSf2(2,4,3))
     &   + abs2(USf(1,1,3,3)) * abs2(USf(1,1,4,3))
     &     * f1f2(MSf2(1,3,3), MSf2(1,4,3))
     &   + abs2(USf(1,1,3,3)) * abs2(USf(2,1,4,3))
     &     * f1f2(MSf2(1,3,3), MSf2(2,4,3))
     &   + abs2(USf(2,1,3,3)) * abs2(USf(1,1,4,3))
     &     * f1f2(MSf2(2,3,3), MSf2(1,4,3))
     &   + abs2(USf(2,1,3,3)) * abs2(USf(2,1,4,3))
     &     * f1f2(MSf2(2,3,3), MSf2(2,4,3))

	if( DeltaRho .gt. DeltaRho_UPPERBOUND ) then
	  print *, "mssm_ini: DeltaRho = ", DeltaRho,
     &      " > DeltaRho_UPPERBOUND = ", DeltaRho_UPPERBOUND
	  return 1
	endif
#endif
	return

999	return 1
	end


*-----------------------------------------------------------------------
* this function is needed for the calculation of Delta rho

	double precision function f1f2(m1, m2)
	implicit none
	double precision m1, m2

#include "model.h"

	double precision f1, f2, logm1m2, c1

	double precision Alfas_MT, sig1, sig12, sig2
	common /radcor/ Alfas_MT, sig1, sig12, sig2

	double precision Li2
	external Li2

	logm1m2 = log(m1/m2)
	c1 = 2*m1*m2/(m1 - m2)*logm1m2
	f1 = m1 + m2 - c1
	f2 = m1 + m2 - c1*(2 + m1/m2*logm1m2)
     &   + (m1 + m2)*m1**2/(m1 - m2)**2 * logm1m2**2
     &   - 2*(m1 - m2) * dble(Li2(1 - m1/m2))

	f1f2 = 3*GF/(8*sqrt2*pi**2) * (f1 + 2/(3*pi)*Alfas_MT*f2)
	end


*-----------------------------------------------------------------------
* Chargino sector:

	subroutine CharginoMasses(*)
	implicit none

#include "model.h"

	double complex X(2, 2)

	X(1, 1) = M_2
	X(1, 2) = sqrt2*MW*SB
	X(2, 1) = sqrt2*MW*CB
	X(2, 2) = MUE
	call SingularValues(2, 2, X, 2, MCha, UCha, 2, VCha, 2)

#ifndef NO_EXCLUSION_LIMITS
	if( MCha(1) .lt. MChargino_LOWERBOUND ) then
	  print *, "mssm_ini: MCha(2) = ", MCha(2),
     &      " < MChargino_LOWERBOUND = ", MChargino_LOWERBOUND
	  return 1
	endif
#endif

	MCha2(1) = MCha(1)**2
	MCha2(2) = MCha(2)**2
	end


*-----------------------------------------------------------------------
* Neutralino sector:

	subroutine NeutralinoMasses(*)
	implicit none

#include "model.h"

	integer i, j
	double complex Y(4, 4), UNeu(4, 4), VNeu(4, 4), phase

	Y(1, 1) = M_1
	Y(2, 1) = 0
	Y(3, 1) = -MZ*SW*CB
	Y(4, 1) = MZ*SW*SB
	Y(1, 2) = Y(2, 1)
	Y(2, 2) = M_2
	Y(3, 2) = MZ*CW*CB
	Y(4, 2) = -MZ*CW*SB
	Y(1, 3) = Y(3, 1)
	Y(2, 3) = Y(3, 2)
	Y(3, 3) = 0
	Y(4, 3) = -MUE
	Y(1, 4) = Y(4, 1)
	Y(2, 4) = Y(4, 2)
	Y(3, 4) = Y(4, 3)
	Y(4, 4) = 0
	call SingularValues(4, 4, Y, 4, MNeu, UNeu, 4, VNeu, 4)

#ifndef NO_EXCLUSION_LIMITS
	if( MNeu(1) .lt. MNeutralino_LOWERBOUND ) then
	  print *, "mssm_ini: MNeu(1) = ", MNeu(1),
     &      " < MNeutralino_LOWERBOUND = ", MNeutralino_LOWERBOUND
	  return 1
	endif
#endif

* ZNeu_{ij} = (diag(sqrt(V . U^T)^*) . V)_{ij}
*           = sqrt(sum_k V_{ik} U_{ik})^* V_{ij}
	do i = 1, 4
	  MNeu2(i) = MNeu(i)**2
	  phase = 0
	  do j = 1, 4
	    phase = phase + VNeu(i, j)*UNeu(i, j)
	  enddo
	  phase = dconjg(sqrt(phase))
	  do j = 1, 4
	    ZNeu(i, j) = phase*VNeu(i, j)
	  enddo
	enddo
	end


*-----------------------------------------------------------------------
* Higgs masses:
* Mh0, MHH, (MHp) are tree-level masses if sig1 = sig12 = sig2 = 0,
* otherwise they are radiatively corrected

	subroutine HiggsMasses(*)
	implicit none

#include "model.h"

	double precision MTr2, M_S2, MtLR, MbLR, Lambda
	double precision Xt, topstop, toplog, mix
	double precision Pb, Pf, Pg, Pg_, P1h, P2h, P2h_
	double precision pref, pref2, tmp, tmp2, tmplog

	double precision Alfas_MT, sig1, sig12, sig2
	common /radcor/ Alfas_MT, sig1, sig12, sig2

	integer Nc, Ng
	double precision Qt, Qb
	parameter (Nc = 3, Ng = 3, Qt = 2/3D0, Qb = -1/3D0)

	MA02 = MA0**2

#ifdef TREE_LEVEL_MH

	sig1 = 0
	sig12 = 0
	sig2 = 0

	MHp2 = MA02 + MW2

#else

* this is the running top quark mass:
	MTr2 = MT2/(1 + 4/(3*pi)*Alfas_MT)**2

#ifdef COMPLEX_PARAMETERS

* Simple approximation: (for complex and/or real case)
* ~~~~~~~~~~~~~~~~~~~~~
	sig2 = -3*GF/(2*sqrt2*pi**2)*MTr2**2/SB2*
     &    log(1 + (MSQ2*MSU2/MTr2 + MSQ2 + MSU2)/MTr2)
	sig1 = 0
	sig12 = 0

#else

* Heinemeyer/Weiglein approximation: hep-ph/9903404 (only for real case)
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	M_S2 = sqrt(MSQ2*MSU2 + MT2*(MSQ2 + MSU2 + MT2))
	MtLR = Af(3, 3) - dconjg(MUE)/TB
	MbLR = Af(4, 3) - dconjg(MUE)*TB
	mix = MtLR**2/M_S2

* The one-loop top/stop contributions:
	Lambda = 1/8D0 - SW2*(1/3D0 - 4/9D0*SW2)
	topstop = MT2/M_S2
	toplog = log(topstop)
	tmp = MZ2/MT2
	sig1 = GF*sqrt2/pi**2*MZ2**2*Lambda*CB2 * toplog
	sig12 = -GF*sqrt2/pi**2*MZ2/TB *
     &    (-3/8D0*MT2 + MZ2*Lambda*SB2) * toplog
	sig2 = GF*sqrt2/(8*pi**2)*MT2**2/SB2 *
     &    ( tmp*(-2 + 11/10D0*tmp)
     &       + (12 + SB2*tmp*(-6 + 8*Lambda*SB2*tmp)) * toplog
     &       + mix*(-12 + 4*tmp + 6*topstop
     &         + mix*(1 + topstop*(-4 + 3*topstop)
     &           + mix*topstop
     &             *(3/5D0 + topstop*(-12/5D0 + 2*topstop)
     &             + mix*topstop
     &               *(3/7D0 + topstop*(-12/7D0 + 3/2D0*topstop))
     &       ))) )

* The dominant two-loop top/stop contributions:
	toplog = log(MTr2/M_S2)
	tmp = (MSf2(2,3,3) - MSf2(1,3,3))/MTr2
     &    * (USf(1,1,3,3)*USf(2,1,3,3))**2
	tmplog = log(MSf2(2,3,3)/MSf2(1,3,3))
	Xt = tmp*(
     &    tmp*(2 - (MSf2(2,3,3) + MSf2(1,3,3))/
     &             (MSf2(2,3,3) - MSf2(1,3,3))*tmplog)
     &    + 2*tmplog)
	tmplog = log(MSf(2,3,3)*MSf(1,3,3)/MTr2)
	sig2 = sig2 + MTr2**2/SB2 * (
     &    + GF*sqrt2/pi**3*Alfas_MT*(toplog*(3*toplog - 6)
     &        - 6*sqrt(mix) + mix*(-3*toplog + 3/4D0*mix))
     &    - 9*GF**2/(16*pi**4)*MTr2* tmplog*(Xt + tmplog) )

* The one-loop leading-log contributions from all other sectors:
* [NB: MSQ is used here instead of MSusy.]
	Pb = Nc*(1 + 4*Qb*SW2*(1 + 2*Qb*SW2))
	Pf = Nc*(Ng - 1)*(2 - SW2*(4 - 8*(Qt**2 + Qb**2)*SW2))
     $    + Ng*(2 - SW2*(4 - 8*SW2))
	Pg = -44 + SW2*(106 - 62*SW2)
	Pg_ = 10 + SW2*(34 - 26*SW2)
	P1h = (1 - 2*SW2*CW2 - 9*C2B**2)*C2B**2
	P2h = -10 + 2*SW2*CW2
	P2h_ = 8 - SW2*(22 - 10*SW2)
	tmplog = log(MSQ2/MZ2)
	tmp = MB2/(MZ2*CB2)
	pref = GF*MZ2**2/(12*sqrt2*pi**2)
	pref2 = GF*Nc/(4*sqrt2*pi**2)*MB2/MSQ2
	sig1 = sig1
     &    - pref*CB2
     &      *(6*Nc*tmp*(2*tmp - 1) + Pb + Pf + Pg + P2h)*tmplog
     &    - pref2*Af(4, 3)
     &      *(4*MB2*MbLR/CB2*(1 - Af(4, 3)*MbLR/(12*MSQ2))
     &        - MZ2*(MbLR + Af(4, 3)/3D0))
	sig12 = sig12
     &    - pref*SB*CB*(3*Nc*tmp - Pb - Pf - Pg_ - P2h_)*tmplog
     &    + .5D0*pref2
     &      *(4*MB2*MUE*MbLR/CB2*(1 - Af(4, 3)*MbLR/(6*MSQ2))
     &        - MZ2*TB*(MbLR*(Af(4, 3) + MUE/TB)
     &                  + (MUE**2 + Af(4, 3)**2)/3D0))
	sig2 = sig2
     &    - pref*SB2*(Pb + Pf + Pg + P2h)*tmplog
     &    + pref2*MUE
     &      *(MB2*MUE*MbLR**2/(3*CB2*MSQ2)
     &        + MZ2*TB*(MbLR + MUE*TB/3D0))
	if( MA0 .gt. MZ ) then
	  tmplog = log(MA02/MZ2)
	  sig1 = sig1 - pref*CB2*(P1h - P2h)*tmplog
	  sig12 = sig12 - pref*SB*CB*(P1h + P2h_)*tmplog
	  sig2 = sig2 - pref*SB2*(P1h - P2h)*tmplog
	endif

#endif COMPLEX_PARAMETERS

	if( MA0 .le. 250 .and. MSQ .ge. 250 ) then
	  MHp2 = MA02 + MW2
     &      + Alfa*MW2/(12*pi*SW2)
     &        *(Nc*(Ng - 1) + Ng - 9 + 15*SW2/CW2)*log(MSQ2/MW2)
     &      + Nc*Alfa/(8*pi*SW2*MW2)
     &        *(2*MTr2*MB2/(SB2*CB2)
     &          - MW2*(MTr2/SB2 + MB2/CB2 - 2/3D0*MW2))
     &        *log(MSQ2/MTr2)
	else
	  MHp2 = MA02 + MW2
	endif

#endif TREE_LEVEL_MH

	tmp = (MA02 + MZ2 - sig2 - sig1)/2D0
	tmp2 = sqrt(((MA02 + MZ2)**2 + (sig1 - sig2)**2)/4D0
     &    + (.5D0*(MA02 - MZ2)*(sig1 - sig2) - MA02*MZ2*C2B)*C2B
     &    + sig12*((MA02 + MZ2)*S2B + sig12))

	Mh02 = tmp - tmp2
	if( Mh02 .lt. 0 ) then
	  print *, "mssm_ini: negative Mh0^2 detected."
	  return 1
	endif
	Mh0 = sqrt(Mh02)

	MH = Mh0
	MH2 = Mh02

#ifndef NO_EXCLUSION_LIMITS
	if( Mh0 .lt. Mh0_LOWERBOUND ) then
	  print *, "mssm_ini: Mh0 = ", Mh0,
     &      " < Mh0_LOWERBOUND = ", Mh0_LOWERBOUND
	  return 1
	endif
#endif

	MHH2 = tmp + tmp2
	MHH = sqrt(MHH2)

	MHp = sqrt(MHp2)
	end


************************************************************************

	subroutine ModelVarIni(sqrtS, *)
	implicit none
	double precision sqrtS

#include "model.h"

	double precision ALPHAS2
	external ALPHAS2

	call SMVarIni(sqrtS, *999)

* Gluino mass: GUT relation
* we use 1/Alfa(s) approx 1/Alfa(MZ2) - 20/(9*pi)*log(s/MZ2)

	MGl = Alfas*(127.934D0 - 40/(9*pi)*log(sqrtS/MZ))*SW2*M_2
#ifndef NO_EXCLUSION_LIMITS
	if( MGl .lt. MGl_LOWERBOUND ) then
	  print *, "mssm_ini: MGl = ", MGl,
     &      " < MGl_LOWERBOUND = ", MGl_LOWERBOUND
	  return 1
	endif
#endif
	MGl2 = MGl**2
	return

999	return 1
	end


************************************************************************

	subroutine ModelDigest
	implicit none

#include "model.h"

	integer i, j

	call SMDigest

	print *, "SUSY parameters"
	print *, "----------------------------"
	print *, "TB =", TB
	print *, "M_1 =", M_1
	print *, "M_2 =", M_2
	print *, "MUE =", MUE
	print *, "MSusy =", MSusy
	print *, "Atau =", Af(2, 3)
	print *, "At =", Af(3, 3)
	print *, "Ab =", Af(4, 3)
	print *, "----------------------------"
	print *, "Higgs masses:"
#ifdef COMPLEX_PARAMETERS
	print *, "(using simple approximation for Higgs mass)"
#else
	print *, "(using 2-loop approximation for Higgs mass)"
#endif
	print *, "Mh0 =", Mh0
	print *, "MHH =", MHH
	print *, "MHp =", MHp
	print *, "MA0 =", MA0
	print *, "----------------------------"
	print *, "Sfermion masses:"
#ifdef NO_SQUARK_MIXING
	print *, "(diagonal sfermion mass matrix: L-R basis)"
#else
	print *, "(non-diagonal sfermion mass matrix: 1-2 basis)"
#endif
3	format(" MSf(1|2, ", I1, ", ", I1, ") =", 2G24.15)
	do j = 1, 4
	  do i = 1, 3
	    print 3, j, i, MSf(1, j, i), MSf(2, j, i)
	  enddo
	enddo
	print *, "----------------------------"
	print *, "Chargino masses:"
	print *, "MCha(1) =", MCha(1)
	print *, "MCha(2) =", MCha(2)
	print *, "----------------------------"
	print *, "Neutralino masses:"
	print *, "MNeu(1) =", MNeu(1)
	print *, "MNeu(2) =", MNeu(2)
	print *, "MNeu(3) =", MNeu(3)
	print *, "MNeu(4) =", MNeu(4)
	print *, "----------------------------"
	print *, "Gluino mass:"
	print *, "MGl =", MGl
	print *, "============================"
	end
