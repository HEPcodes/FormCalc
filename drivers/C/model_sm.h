#if 0
	model_sm.h
	declarations for model_sm.F
	this file is part of FormCalc
	last modified 5 Jan 16 th
#endif


#define CKMlambda .22535
#define CKMA .811
#define CKMrhobar .131
#define CKMetabar .345

#define MZ 91.1876
#define MZ2 (MZ*MZ)
#define MW 80.385
#define MW2 (MW*MW)
#define CW (MW/MZ)
#define CW2 (CW*CW)
#define SW2 (1 - CW2)

#define GF 1.16637e-5
#define Alfa0 (1/137.035999074)
#define AlfaGF (sqrt2/pi*GF*MW2*SW2)
#define AlfaMZ (1/127.944)
#define Alfa Alfa0
#define Alfa2 (Alfa*Alfa)
#define DeltaAlfa5Had .027547
#define AlfasMZ .1184

#define MH_ 125

#define ME .5109989280e-3
#define ME2 (ME*ME)
#define MM 105.6583715e-3
#define MM2 (MM*MM)
#define ML 1776.82e-3
#define ML2 (ML*ML)

#define MU Mf(3,1)
#define MU2 Mf2(3,1)
#define MC 1.275
#define MC2 (MC*MC)
#define MT 173.21
#define MT2 (MT*MT)

#define MD Mf(4,1)
#define MD2 Mf2(4,1)
#define MS 95e-3
#define MS2 (MS*MS)
#define MB1S 4.66
#define MBatMB 4.18
#define MB Mf(4,3)
#define MB2 Mf2(4,3)

#if NOUNDERSCORE
#define smpara_ smpara
#endif

struct smpara_ {
  ComplexType CKM[3][3];
  RealType Mf[3][4], Mf2[3][4], MH, MH2;
  RealType EL, ELMZ, GS, Alfas, Alfas2, SW;
} smpara_;

#define CKM(i,j) smpara_.CKM[j-1][i-1]
#define CKMC(i,j) Conjugate(CKM(i,j))
#define Mf(t,g) smpara_.Mf[g-1][t-1]
#define Mf2(t,g) smpara_.Mf2[g-1][t-1]
#define MH smpara_.MH
#define MH2 smpara_.MH2
#define EL smpara_.EL
#define ELMZ smpara_.ELMZ
#define GS smpara_.GS
#define Alfas smpara_.Alfas
#define Alfas2 smpara_.Alfas2
#define SW smpara_.SW

