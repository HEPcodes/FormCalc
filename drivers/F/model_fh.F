* model_fh.F
* parameter initialization via FeynHiggs
* this file is part of FormCalc
* last modified 18 Mar 19 th


#ifndef MB_SCHEME
#define MB_SCHEME ONSHELL
#endif

#define ONSHELL -1
#define MSBAR(o) o


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
#define MNeu_LOWERBOUND 30

* hep-ex/9910007 (L3) and hep-ex/9909051 (OPAL)
#define MCha_LOWERBOUND 90

* hep-ex/9910020 (L3) and PDG
#define MGl_LOWERBOUND 175


************************************************************************

	subroutine ModelDefaults(argc, argv)
	implicit none
	integer argc
	character*128 argv(*)

#include "decl.h"
#include "SLHA.h"

	ComplexType slhadata(nslhadata)
	integer error

	character*16 flags
	data flags /"400202113"/

	if( argc < 1 .or. argc > 2 )
     &    stop "model_fh: need arguments <file> [<flags>]"

	if( argc >= 2 ) flags = argv(2)
	call FHSetFlagsString(error, flags)
	if( error /= 0 ) stop

	call FHReadRecord(error, record, slhadata, argv(1))
* 0 = SLHA, 2 = FH file, all others are true errors
	if( ibclr(error, 1) /= 0 ) stop

* turn on verbose mode:   
	record(iAdmin,iStep) = 1
	end


************************************************************************

	subroutine ModelConstIni(fail)
	implicit none
	integer fail

#include "decl.h"

	RealType invAlfaMZ, scalefactor, asMB, Deltab
	integer i, nf
#ifndef NO_EXCLUSION_LIMITS
	RealType gm2, DeltaRho
	RealType MWMSSM, MWSM, SW2MSSM, SW2SM
	RealType edmeTh, edmn, edmHg
#endif

	call FHLoopRecord(fail, record)
	if( fail /= 0 ) return

	call FHSetRecord(fail, record)
	if( fail /= 0 ) return

	call FHRetrieveSMPara(fail,
     &    invAlfaMZ, AlfasMZ, GF,
     &    ME, MU, MD, MM, MC, MS, ML, MBatMB,
     &    MW, MZ,
     &    CKMlambda, CKMA, CKMrhobar, CKMetabar)
	if( fail /= 0 ) return

	call FHGetSMPara(fail, CKM)
	if( fail /= 0 ) return

	call FHRetrieveOSPara(fail, scalefactor,
     &    MT, TB, MA0, MHp,
     &    MSS(1,3), MSS(2,3), MSS(3,3), MSS(4,3), MSS(5,3),
     &    MSS(1,2), MSS(2,2), MSS(3,2), MSS(4,2), MSS(5,2),
     &    MSS(1,1), MSS(2,1), MSS(3,1), MSS(4,1), MSS(5,1),
     &    MUE,
     &    Afd(2,3), Afd(3,3), Afd(4,3),
     &    Afd(2,2), Afd(3,2), Afd(4,2),
     &    Afd(2,1), Afd(3,1), Afd(4,1),
     &    Mino1, Mino2, Mino3)
	if( fail /= 0 ) return

	call FHGetPara(fail, fv, MSf, USf, MASf, UASf,
     &    MCha, UCha, VCha, MNeu, ZNeu, Deltab, MGl,
     &    MHtree, SA)
	if( fail /= 0 ) return

	call FHGetFV(fail, MSS2, Kf)
	if( fail /= 0 ) return

#ifndef NO_EXCLUSION_LIMITS
	if( MSf(1,4,3) < MSbottom_LOWERBOUND ) then
	  Error(fail, "MSf(1,4,3) < MSbottom_LOWERBOUND")
	  return
	endif

	if( MSf(1,3,3) < MStop_LOWERBOUND ) then
	  Error(fail, "MSf(1,3,3) < MStop_LOWERBOUND")
	  return
	endif

	if( MSf(1,2,1) < MSlepton_LOWERBOUND ) then
	  Error(fail, "MSf(1,2,1) < MSlepton_LOWERBOUND")
	  return
	endif

	if( MSf(1,3,1) < MSquark_LOWERBOUND ) then
	  Error(fail, "MSf(1,3,1) < MSquark_LOWERBOUND")
	  return
	endif

	if( MCha(1) < MCha_LOWERBOUND ) then
	  Error(fail, "MCha(1) < MCha_LOWERBOUND")
	  return
	endif

	if( MNeu(1) < MNeu_LOWERBOUND ) then
	  Error(fail, "MNeu(1) < MNeu_LOWERBOUND")
	  return
	endif
#endif

	call FHHiggsCorr(fail, MHiggs, SAeff, UHiggs, ZHiggs)
	if( fail /= 0 ) return

#ifndef NO_EXCLUSION_LIMITS
	if( Mh0 < Mh0_LOWERBOUND ) then
	  Error(fail, "Mh0 < Mh0_LOWERBOUND")
	  return
	endif

	call FHConstraints(fail, gm2, DeltaRho,
     &    MWMSSM, MWSM, SW2MSSM, SW2SM,
     &    edmeTh, edmn, edmHg)

	if( DeltaRho > DeltaRho_UPPERBOUND ) then
	  Error(fail, "DeltaRho > DeltaRho_UPPERBOUND")
	  return
	endif
#endif

	call AlphaS(asMB, nf, MBatMB)
	MB = MBatMB*sqrt(1 + 4*2/(3*pi)*asMB)

	Mf2 = Mf**2

	MW2 = MW**2
	MZ2 = MZ**2
	CW = MW/MZ
	CW2 = CW**2
	SW2 = (1 - CW)*(1 + CW)
	SW = sqrt(SW2)

	AlfaGF = sqrt2/pi*GF*MW2*SW2

	AlfaMZ = 1/invAlfaMZ
	ELMZ = sqrt(4*pi*AlfaMZ)

	Alfa = Alfa0
	Alfa2 = Alfa**2
	EL = sqrt(4*pi*Alfa)

	MH = Mh0
	MH2 = MH**2

	MHiggs2 = MHiggs**2
	MHtree2 = MHtree**2

	MSf2 = MSf**2
	MASf2 = MASf**2
	MCha2 = MCha**2
	MNeu2 = MNeu**2

	TB2 = TB**2
	CB2 = 1/(1 + TB2)
	SB2 = CB2*TB2
	CB = sqrt(CB2)
	SB = CB*TB
	S2B = 2*TB*CB2
	C2B = (CB - SB)*(CB + SB)

	SA2 = SA**2  
	CA2 = (1 - SA)*(1 + SA)
	CA  = sqrt(CA2)
	S2A = 2*SA*CA
	C2A = (CA - SA)*(CA + SA)
	SAB = SA*CB + CA*SB
	CAB = CA*CB - SA*SB
	SBA = SB*CA - CB*SA
	SBA2 = SBA**2 
	CBA = CB*CA + SB*SA
	CBA2 = CBA**2
	end


************************************************************************

	subroutine ModelVarIni(fail, Q)
	implicit none
	integer fail
	RealType Q

#include "decl.h"

	integer nf

	fail = 0

#ifndef ALPHAS_ORDER
#define ALPHAS_ORDER 3
#endif

	call AlphaS(Alfas, nf, Q, ALPHAS_ORDER)
	Alfas2 = Alfas**2
	GS = sqrt(4*pi*Alfas)

#if MB_SCHEME != ONSHELL
	MB = MBatMB
	call MqRun(MB, MBatMB, Q, MB_SCHEME)
	MB2 = MB**2
#endif
	end


************************************************************************

	subroutine ModelDigest
	implicit none

#include "decl.h"

	integer error

	call FHOutput(error, "stdout", 7, -1D0)
	end

