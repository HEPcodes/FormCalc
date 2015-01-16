* mssm_ini.F
* all parameters for the MSSM are initialized here
* this file is part of FormCalc
* Christian Schappacher, 17 Feb 99
* last modified 25 Jan 05 th


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
* (they are usually set in run.F):

* NO_SQUARK_MIXING
*   If defined, sets Atau, At, and Ab such that the off-diagonal entries
*   of the sfermion mixing matrices vanish, i.e. gauge eigenstates =
*   mass eigenstates.

* HIGGS_MASSES can take the values
* TREE -- uses only tree-level Higgs masses,
* SIMPLE -- uses a simple one-loop approximation of the
*    Higgs masses which is valid for all parameters.
* FEYNHIGGS -- invokes FeynHiggs to compute the Higgs masses,
* (undefined) -- uses a fairly precise two-loop approximation,
*    which is valid though only for real parameters.

* SM_HIGGS
*   If defined, reverse-engineers the mixing in the Higgs sector
*   (alpha and beta) such that the MSSM Higgs sector looks like
*   a SM Higgs sector.

* NO_EXCLUSION_LIMITS
*   If defined, the experimental exclusion limits are ignored.

* NO_GUT_RELATION
*   If defined, does not use the GUT relation which relates
*   the gaugino mass paramters.  In this case, M_1 and MGl are
*   also input parameters.


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


#define TREE 100
#define SIMPLE 200
#define FEYNHIGGS 300

#define ModelDefaults SMDefaults
#define ModelConstIni SMConstIni
#define ModelVarIni SMVarIni
#define ModelDigest SMDigest
#include "sm_ini.F"
#undef ModelDefaults
#undef ModelConstIni
#undef ModelVarIni
#undef ModelDigest


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

	double precision ALPHAS2
	external ALPHAS2

	call SMConstIni(*999)

	AlfasMT = ALPHAS2(MT)

	TB2 = TB**2
	CB2 = 1/(1 + TB2)
	SB2 = 1 - CB2
	CB = sqrt(CB2)
	SB = CB*TB
	S2B = 2*TB*CB2
	C2B = CB2 - SB2

	call SfermionMasses(*999)
	call CharginoMasses(*999)
	call NeutralinoMasses(*999)

* FeynHiggs needs the gluino mass, hence we have to postpone
* the Higgs sector to ModelVarIni.
#if HIGGS_MASSES == FEYNHIGGS
#else
	call HiggsMasses(*999)
#endif
	return

999	return 1
	end


*-----------------------------------------------------------------------
* Sfermion sector:

	subroutine SfermionMasses(*)
	implicit none

#include "model.h"

	integer t, g
	double precision i3part, qpart

	double precision Qe, Qu, Qd
	parameter (Qe = -1, Qu = 2/3D0, Qd = -1/3D0)

#ifndef NO_EXCLUSION_LIMITS
	double precision f1f2, DeltaRho
	external f1f2
#endif

	i3part = .5D0*MZ2*C2B
	qpart = SW2*MZ2*C2B
	DSf(1,1) = +i3part
	DSf(1,2) = 0
	DSf(2,2) = Qe*qpart
	DSf(1,2) = -i3part - DSf(2,2)
	DSf(2,3) = Qu*qpart
	DSf(1,3) = +i3part - DSf(2,3)
	DSf(2,4) = Qd*qpart
	DSf(1,4) = -i3part - DSf(2,4)

	do g = 1, 3
	  MSS(1,2,g) = MSL
	  MSS(2,2,g) = MSE
	  MSS(1,3,g) = MSQ
	  MSS(2,3,g) = MSU
	  MSS(1,4,g) = MSQ
	  MSS(2,4,g) = MSD

	  MSS2(1,2,g) = MSS(1,2,g)**2
	  MSS2(2,2,g) = MSS(2,2,g)**2
	  MSS2(1,3,g) = MSS(1,3,g)**2
	  MSS2(2,3,g) = MSS(2,3,g)**2
	  MSS2(1,4,g) = MSS(1,4,g)**2
	  MSS2(2,4,g) = MSS(2,4,g)**2

#ifdef NO_SQUARK_MIXING
#define Eigen2x2 NoMixing
	  Af(2,g) = dconjg(MUE)*TB
	  Af(3,g) = dconjg(MUE)/TB
	  Af(4,g) = Af(2,g)
	  Xf(2,g) = 0
	  Xf(3,g) = 0
	  Xf(4,g) = 0
#else
	  Af(2,g) = Atau
	  Xf(2,g) = Af(2,g) - dconjg(MUE)*TB
	  Af(3,g) = At
	  Xf(3,g) = Af(3,g) - dconjg(MUE)/TB
	  Af(4,g) = Ab
	  Xf(4,g) = Af(4,g) - dconjg(MUE)*TB
#endif

	  call NoMixing(MSf(1,1,g), MSf2(1,1,g), USf(1,1,1,g),
     &      MSS2(1,2,g) + DSf(1,1), 0D0, dcmplx(0D0), *999)

	  do t = 2, 4
	    call Eigen2x2(MSf(1,t,g), MSf2(1,t,g), USf(1,1,t,g),
     &        MSS2(1,t,g) + DSf(1,t) + Mf2(t,g),
     &        MSS2(2,t,g) + DSf(2,t) + Mf2(t,g),
     &        Mf(t,g)*dconjg(Xf(t,g)), *999)
	  enddo
	enddo

#ifndef NO_EXCLUSION_LIMITS
	if( MSf(1,4,3) .lt. MSbottom_LOWERBOUND ) then
	  print *, "mssm_ini: MSf(1,4,3) = ", MSf(1,4,3),
     &      " < MSbottom_LOWERBOUND = ", MSbottom_LOWERBOUND
	  return 1
	endif

	if( MSf(1,3,3) .lt. MStop_LOWERBOUND ) then
	  print *, "mssm_ini: MSf(1,3,3) = ", MSf(1,3,3),
     &      " < MStop_LOWERBOUND = ", MStop_LOWERBOUND
	  return 1
	endif

	if( MSf(1,2,1) .lt. MSlepton_LOWERBOUND ) then
	  print *, "mssm_ini: MSf(1,2,1) = ", MSf(1,2,1),
     &      " < MSlepton_LOWERBOUND = ", MSlepton_LOWERBOUND
	  return 1
	endif

	if( MSf(1,3,1) .lt. MSquark_LOWERBOUND ) then
	  print *, "mssm_ini: MSf(1,3,1) = ", MSf(1,3,1),
     &      " < MSquark_LOWERBOUND = ", MSquark_LOWERBOUND
	  return 1
	endif

#define abs2(c) (dble(c)**2 + dimag(c)**2)

	DeltaRho =
     &     USf(1,1,3,3)*USfC(1,2,3,3)
     &     * USf(2,2,3,3)*USfC(2,1,3,3)
     &     * f1f2(MSf2(1,3,3), MSf2(2,3,3))
     &   + USf(1,1,4,3)*USfC(1,2,4,3)
     &     * USf(2,2,4,3)*USfC(2,1,4,3)
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

	subroutine NoMixing(m, msq, U, M11, M22, M12, *)
	implicit none
	double precision m(*), msq(*)
	double complex U(2,*)
	double precision M11, M22
	double complex M12

	msq(1) = M11
	m(1) = sqrt(M11)
	msq(2) = M22
	m(2) = sqrt(M22)
	U(1,1) = 1
	U(2,2) = 1
	U(2,1) = 0
	U(1,2) = 0
	end


*-----------------------------------------------------------------------
* this function is needed for the calculation of Delta rho

	double precision function f1f2(m1, m2)
	implicit none
	double precision m1, m2

#include "model.h"

	double precision f1, f2, logm1m2, c1

	double precision Li2
	external Li2

	logm1m2 = log(m1/m2)
	c1 = 2*m1*m2/(m1 - m2)*logm1m2
	f1 = m1 + m2 - c1
	f2 = m1 + m2 - c1*(2 + m1/m2*logm1m2)
     &   + (m1 + m2)*m1**2/(m1 - m2)**2 * logm1m2**2
     &   - 2*(m1 - m2) * Li2(1 - m1/m2)

	f1f2 = 3*GF/(8*sqrt2*pi**2) * (f1 + 2/(3*pi)*AlfasMT*f2)
	end


*-----------------------------------------------------------------------
* Chargino sector:

	subroutine CharginoMasses(*)
	implicit none

#include "model.h"

	double complex X(2,2)

	X(1,1) = M_2
	X(2,1) = sqrt2*MW*CB
	X(1,2) = sqrt2*MW*SB
	X(2,2) = MUE

	call SingularValues(2, 2, X, 2, MCha, UCha, 2, VCha, 2)

#ifndef NO_EXCLUSION_LIMITS
	if( MCha(1) .lt. MChargino_LOWERBOUND ) then
	  print *, "mssm_ini: MCha(1) = ", MCha(1),
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

	double complex Y(4,4)

#ifndef NO_GUT_RELATION
* GUT relation (Higgs Hunter's Guide p. 287)
	M_1 = 5/3D0*SW2/CW2*M_2
#endif

	Y(1,1) = M_1
	Y(2,2) = M_2
	Y(3,1) = -MZ*SW*CB
	Y(1,3) = Y(3,1)
	Y(4,1) = MZ*SW*SB
	Y(1,4) = Y(4,1)
	Y(3,2) = MZ*CW*CB
	Y(2,3) = Y(3,2)
	Y(4,2) = -MZ*CW*SB
	Y(2,4) = Y(4,2)
	Y(4,3) = -MUE
	Y(3,4) = Y(4,3)
	Y(2,1) = 0
	Y(1,2) = 0
	Y(3,3) = 0
	Y(4,4) = 0

	call PseudoEigensystem(4, Y, 4, MNeu, ZNeu, 4)

#ifndef NO_EXCLUSION_LIMITS
	if( MNeu(1) .lt. MNeutralino_LOWERBOUND ) then
	  print *, "mssm_ini: MNeu(1) = ", MNeu(1),
     &      " < MNeutralino_LOWERBOUND = ", MNeutralino_LOWERBOUND
	  return 1
	endif
#endif

	MNeu2(1) = MNeu(1)**2
	MNeu2(2) = MNeu(2)**2
	MNeu2(3) = MNeu(3)**2
	MNeu2(4) = MNeu(4)**2
	end


*-----------------------------------------------------------------------
* Higgs sector:

	subroutine HiggsMasses(*)
	implicit none

#include "model.h"

	double precision avgMH2, deltaMH2

#if HIGGS_MASSES == FEYNHIGGS
	integer mssmpart, fieldren, tanbren, higgsmix, p2approx
	integer looplevel, tl_running_mt, tl_bot_resum
	parameter (mssmpart = 4)
	parameter (fieldren = 0)
	parameter (tanbren = 0)
	parameter (higgsmix = 2)
	parameter (p2approx = 0)
	parameter (looplevel = 2)
	parameter (tl_running_mt = 1)
	parameter (tl_bot_resum = 1)

	double precision scalefactor, lambda_t, lambda_b
	double precision Qtau, Qt, Qb
	parameter (scalefactor = 1)
	parameter (lambda_t = 0)
	parameter (lambda_b = 0)
	parameter (Qtau = -1)
	parameter (Qt = -1)
	parameter (Qb = -1)

	integer error
	double precision MHiggs(4), SAeff
	double complex UHiggs(3,3)
#else
	double precision sig11, sig12, sig22
	double precision MTr2, Qsusy, Lambda
	double precision Xt, topstop, toplog, mix
	double precision Pb, Pf, Pg, Pg_, P1h, P2h, P2h_
	double precision pref, pref2, tmp, tmplog, TA

	double precision MtLR, MbLR
	equivalence (Xf(3,3), MtLR)
	equivalence (Xf(4,3), MbLR)

	integer Nc, Ng
	double precision Qt, Qb
	parameter (Nc = 3, Ng = 3, Qt = 2/3D0, Qb = -1/3D0)
#endif

	MA02 = MA0**2

	avgMH2 = .5D0*(MA02 + MZ2)
	deltaMH2 = sqrt(avgMH2**2 - MA02*MZ2*C2B**2)
	Mh02tree = avgMH2 - deltaMH2
	MHH2tree = avgMH2 + deltaMH2

	MHp2tree = MA02 + MW2

#if HIGGS_MASSES == FEYNHIGGS
	call FHSetFlags(error, mssmpart, fieldren, tanbren,
     &    higgsmix, p2approx, looplevel, tl_running_mt, tl_bot_resum)
	if( error .ne. 0 ) return 1

	call FHSetPara(error, scalefactor,
     &    MT, MB, MW, MZ, TB, MA0, -1D0,
     &    MSS(1,2,3), MSS(2,2,3), MSS(1,3,3), MSS(2,3,3), MSS(2,4,3),
     &    MSS(1,2,2), MSS(2,2,2), MSS(1,3,2), MSS(2,3,2), MSS(2,4,2),
     &    MSS(1,2,1), MSS(2,2,1), MSS(1,3,1), MSS(2,3,1), MSS(2,4,1),
     &    MUE,
     &    Af(2,3), Af(3,3), Af(4,3),
     &    Af(2,2), Af(3,2), Af(4,2),
     &    Af(2,1), Af(3,1), Af(4,1),
     &    dcmplx(M_1), dcmplx(M_2), dcmplx(MGl),
     &    lambda_t, lambda_b,
     &    Qtau, Qt, Qb)
        if( error .ne. 0 ) return 1

	call FHHiggsCorr(error, MHiggs, SA, UHiggs)
	if( error .ne. 0 ) return 1

        Mh0 = MHiggs(1)
	Mh02 = Mh0**2

	MHH = MHiggs(2)
	MHH2 = MHH**2

	MA0 = MHiggs(3)
	MA02 = MA0**2

	MHp = MHiggs(4)
	MHp2 = MHp**2

	SA2 = SA**2
	CA2 = 1 - SA2
	CA = sqrt(CA2)

#else
* no FeynHiggs

	MHp2 = MHp2tree

#if HIGGS_MASSES == TREE

	sig11 = 0
	sig12 = 0
	sig22 = 0

#else

* this is the running top quark mass:
	MTr2 = MT2/(1 + 4/(3*pi)*AlfasMT)**2

#define SusyScale(mt2) sqrt((MSS2(1,3,3) + mt2)*(MSS2(2,3,3) + mt2))

#if HIGGS_MASSES == SIMPLE

* Simple approximation: (for complex and/or real case)
* ~~~~~~~~~~~~~~~~~~~~~
	Qsusy = SusyScale(MTr2)
	sig11 = 0
	sig12 = 0
	sig22 = -3*GF/(sqrt2*pi**2)*MTr2**2/SB2*log(Qsusy/MTr2)

#else

* Heinemeyer/Weiglein approximation: hep-ph/9903404 (only for real case)
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Qsusy = SusyScale(MT2)
	mix = MtLR**2/Qsusy

* The one-loop top/stop contributions:
	Lambda = 1/8D0 - SW2*(1/3D0 - 4/9D0*SW2)
	topstop = MT2/Qsusy
	toplog = log(topstop)
	tmp = MZ2/MT2
	sig11 = GF*sqrt2/pi**2*MZ2**2*Lambda*CB2 * toplog
	sig12 = -GF*sqrt2/pi**2*MZ2/TB *
     &    (-3/8D0*MT2 + MZ2*Lambda*SB2) * toplog
	sig22 = GF*sqrt2/(8*pi**2)*MT2**2/SB2 *
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
	toplog = log(MTr2/Qsusy)
	tmp = (MSf2(2,3,3) - MSf2(1,3,3))/MTr2
     &    * (USf(1,1,3,3)*USf(2,1,3,3))**2
	tmplog = log(MSf2(2,3,3)/MSf2(1,3,3))
	Xt = tmp*(
     &    tmp*(2 - (MSf2(2,3,3) + MSf2(1,3,3))/
     &             (MSf2(2,3,3) - MSf2(1,3,3))*tmplog)
     &    + 2*tmplog)
	tmplog = log(MSf(2,3,3)*MSf(1,3,3)/MTr2)
	sig22 = sig22 + MTr2**2/SB2 * (
     &    + GF*sqrt2/pi**3*AlfasMT*(toplog*(3*toplog - 6)
     &        - 6*sqrt(mix) + mix*(-3*toplog + 3/4D0*mix))
     &    - 9*GF**2/(16*pi**4)*MTr2* tmplog*(Xt + tmplog) )

* The one-loop leading-log contributions from all other sectors:
* [NB: MSQ is used here instead of MSusy.]
	Pb = Nc*(1 + 4*Qb*SW2*(1 + 2*Qb*SW2))
	Pf = Nc*(Ng - 1)*(2 - SW2*(4 - 8*(Qt**2 + Qb**2)*SW2))
     &    + Ng*(2 - SW2*(4 - 8*SW2))
	Pg = -44 + SW2*(106 - 62*SW2)
	Pg_ = 10 + SW2*(34 - 26*SW2)
	P1h = (1 - 2*SW2*CW2 - 9*C2B**2)*C2B**2
	P2h = -10 + 2*SW2*CW2
	P2h_ = 8 - SW2*(22 - 10*SW2)
	tmplog = log(MSS2(1,3,3)/MZ2)
	tmp = MB2/(MZ2*CB2)
	pref = GF*MZ2**2/(12*sqrt2*pi**2)
	pref2 = GF*Nc/(4*sqrt2*pi**2)*MB2/MSS2(1,3,3)
	sig11 = sig11
     &    - pref*CB2
     &      *(6*Nc*tmp*(2*tmp - 1) + Pb + Pf + Pg + P2h)*tmplog
     &    - pref2*Af(4,3)
     &      *(4*MB2*MbLR/CB2*(1 - Af(4,3)*MbLR/(12*MSS2(1,3,3)))
     &        - MZ2*(MbLR + Af(4,3)/3D0))
	sig12 = sig12
     &    - pref*SB*CB*(3*Nc*tmp - Pb - Pf - Pg_ - P2h_)*tmplog
     &    + .5D0*pref2
     &      *(4*MB2*MUE*MbLR/CB2*(1 - Af(4,3)*MbLR/(6*MSS2(1,3,3)))
     &        - MZ2*TB*(MbLR*(Af(4,3) + MUE/TB)
     &                  + (MUE**2 + Af(4,3)**2)/3D0))
	sig22 = sig22
     &    - pref*SB2*(Pb + Pf + Pg + P2h)*tmplog
     &    + pref2*MUE
     &      *(MB2*MUE*MbLR**2/(3*CB2*MSS2(1,3,3))
     &        + MZ2*TB*(MbLR + MUE*TB/3D0))
	if( MA0 .gt. MZ ) then
	  tmplog = log(MA02/MZ2)
	  sig11 = sig11 - pref*CB2*(P1h - P2h)*tmplog
	  sig12 = sig12 - pref*SB*CB*(P1h + P2h_)*tmplog
	  sig22 = sig22 - pref*SB2*(P1h - P2h)*tmplog
	endif

#endif

	if( MA0 .le. 250 .and. MSS(1,3,3) .ge. 250 ) then
	  MHp2 = MHp2 +
     &      Alfa*MW2/(12*pi*SW2)*
     &        (Nc*(Ng - 1) + Ng - 9 + 15*SW2/CW2)*
     &        log(MSS2(1,3,3)/MW2) +
     &      Nc*Alfa/(8*pi*SW2*MW2)*
     &        (2*MTr2*MB2/(SB2*CB2)
     &          - MW2*(MTr2/SB2 + MB2/CB2 - 2/3D0*MW2))*
     &        log(MSS2(1,3,3)/MTr2)
	endif

#endif

	avgMH2 = .5D0*(MA02 + MZ2 - sig22 - sig11)
	deltaMH2 = sqrt(.25D0*((MA02 + MZ2)**2 + (sig11 - sig22)**2) +
     &    (.5D0*(MA02 - MZ2)*(sig11 - sig22) - MA02*MZ2*C2B)*C2B +
     &    sig12*((MA02 + MZ2)*S2B + sig12))

	Mh02 = avgMH2 - deltaMH2
	if( Mh02 .lt. 0 ) then
	  print *, "mssm_ini: negative Mh0^2 detected."
	  return 1
	endif
	Mh0 = sqrt(Mh02)

	MHH2 = avgMH2 + deltaMH2
	MHH = sqrt(MHH2)

	MHp = sqrt(MHp2)

	TA = (-(MA02 + MZ2)*SB*CB - sig12)/
     &    (MZ2*CB2 + MA02*SB2 - sig11 - Mh02)
	CA2 = 1/(1 + TA**2)
	SA2 = 1 - CA2
	SA = -sqrt(SA2)
	CA = SA/TA

#endif

#ifndef NO_EXCLUSION_LIMITS
	if( Mh0 .lt. Mh0_LOWERBOUND ) then
	  print *, "mssm_ini: Mh0 = ", Mh0,
     &      " < Mh0_LOWERBOUND = ", Mh0_LOWERBOUND
	  return 1
	endif
#endif

	MH = Mh0
	MH2 = Mh02

#ifdef SM_HIGGS
* With SM_HIGGS defined, the values of the MSSM parameters are set
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
	if( C2B .gt. 1 ) goto 999
	CB2 = .5D0*(1 + C2B)
	SB2 = 1 - CB2
	TB2 = SB2/CB2
	TB = sqrt(TB2)
	CB = sqrt(CB2)
	SB = CB*TB
	S2B = 2*TB*CB2

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
	S2A = 2*SA*CA
	C2A = CA2 - SA2
	SAB = SA*CB + CA*SB
	CAB = CA*CB - SA*SB
	SBA = SB*CA - CB*SA
	CBA = CB*CA + SB*SA
#endif
	end


************************************************************************

	subroutine ModelVarIni(sqrtS, *)
	implicit none
	double precision sqrtS

#include "model.h"

	call SMVarIni(sqrtS, *999)

#ifndef NO_GUT_RELATION
* Gluino mass: GUT relation
* we use 1/Alfa(s) approx 1/Alfa(MZ2) - 20/(9*pi)*log(s/MZ2)
	MGl = Alfas*(127.934D0 - 40/(9*pi)*log(sqrtS/MZ))*SW2*M_2
#endif

#ifndef NO_EXCLUSION_LIMITS
	if( MGl .lt. MGl_LOWERBOUND ) then
	  print *, "mssm_ini: MGl = ", MGl,
     &      " < MGl_LOWERBOUND = ", MGl_LOWERBOUND
	  return 1
	endif
#endif
	MGl2 = MGl**2

#if HIGGS_MASSES == FEYNHIGGS
	call HiggsMasses(*999)
#endif
	return

999	return 1
	end


************************************************************************

	subroutine ModelDigest
	implicit none

#include "model.h"

	call SMDigest

	print *, "MSSM parameters"
	print *, "----------------------------"
	print *, "TB =", TB
	print *, "M_1 =", M_1
	print *, "M_2 =", M_2
	print *, "MUE =", MUE
	print *, "MSusy =", MSusy
	print *, "Atau =", Af(2,3)
	print *, "At =", Af(3,3)
	print *, "Ab =", Af(4,3)
	print *, "----------------------------"
	print *, "Higgs masses:"
#if HIGGS_MASSES == FEYNHIGGS
	print *, "(using FeynHiggs)"
#elif HIGGS_MASSES == TREE
	print *, "(using tree-level Higgs masses)"
#elif HIGGS_MASSES == SIMPLE
	print *, "(using simple approximation)"
#else
	print *, "(using 2-loop approximation)"
#endif
	print *, "Mh0 =", Mh0
	print *, "MHH =", MHH
	print *, "MA0 =", MA0
	print *, "MHp =", MHp
	print *, "----------------------------"
	print *, "Sfermion masses:"
#ifdef NO_SQUARK_MIXING
	print *, "(diagonal sfermion mass matrix: L-R basis)"
#else
	print *, "(non-diagonal sfermion mass matrix: 1-2 basis)"
#endif
	print *, "Sneutrino  e:  ", MSf(1,1,1)
	print *, "           mu: ", MSf(1,1,2)
	print *, "           tau:", MSf(1,1,3)
	print *, "Slepton    e:  ", MSf(1,2,1), MSf(2,2,1)
	print *, "           mu: ", MSf(1,2,2), MSf(2,2,2)
	print *, "           tau:", MSf(1,2,3), MSf(2,2,3)
	print *, "u-Squark   u:  ", MSf(1,3,1), MSf(2,3,1)
	print *, "           c:  ", MSf(1,3,2), MSf(2,3,2)
	print *, "           t:  ", MSf(1,3,3), MSf(2,3,3)
	print *, "d-Squark   d:  ", MSf(1,4,1), MSf(2,4,1)
	print *, "           s:  ", MSf(1,4,2), MSf(2,4,2)
	print *, "           b:  ", MSf(1,4,3), MSf(2,4,3)
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
