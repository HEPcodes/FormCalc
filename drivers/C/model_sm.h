#if 0
	model_sm.h
	declarations for model_sm.F
	this file is part of FormCalc
	last modified 8 Mar 13 th
#endif


#define CKMlambda .2253
#define CKMA .808
#define CKMrhobar .132
#define CKMetabar .341

#define MZ 91.1876
#define MZ2 (MZ*MZ)
#define MW 80.399
#define MW2 (MW*MW)
#define CW (MW/MZ)
#define CW2 (CW*CW)
#define SW2 (1 - CW2)

#define GF 1.16637e-5
#define Alfa (1/137.035999679)
//#define Alfa (sqrt2/pi*GF*MW2*SW2
#define Alfa2 (Alfa*Alfa)
#define AlfaMZ (1/127.934)
#define AlfasMZ .1184

#define ME .510998910e-3
#define ME2 (ME*ME)
#define MM 105.658367e-3
#define MM2 (MM*MM)
#define ML 1776.82e-3
#define ML2 (ML*ML)

#define MU 53.8e-3
#define MU2 (MU*MU)
#define MC 1.27
#define MC2 (MC*MC)
#define MT 172.
#define MT2 (MT*MT)

#define MD 53.8e-3
#define MD2 (MD*MD)
#define MS 101e-3
#define MS2 (MS*MS)
#define MB 4.67
#define MB2 (MB*MB)
#define MBatMB 4.25

#if NOUNDERSCORE
#define smpara_ smpara
#endif

struct smpara_ {
  ComplexType CKM[3][3];
  RealType Mf[3][4], Mf2[3][4];
  RealType MH, MH2;
  RealType EL, GS, Alfas, Alfas2, SW;
} smpara_;

#define CKM(i,j) smpara_.CKM[j-1][i-1]
#define CKMC(i,j) Conjugate(CKM(i,j))
#define Mf(t,g) smpara_.Mf[g-1][t-1]
#define Mf2(t,g) smpara_.Mf2[g-1][t-1]
#define MH smpara_.MH
#define MH2 smpara_.MH2
#define EL smpara_.EL
#define GS smpara_.GS
#define Alfas smpara_.Alfas
#define Alfas2 smpara_.Alfas2
#define SW smpara_.SW

