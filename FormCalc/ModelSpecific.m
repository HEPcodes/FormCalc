(*
	ModelSpecific.m
		global definitions for specific models
		this file is part of FormCalc
		last modified 21 Dec 18 th

Note: This file is read by FormCalc only if $NoModelSpecific is not True.

*)


powQ[1] = powQ[-1] = False

powQ[other_] := IntegerQ[other]

Sq/: (Sq[v_] = v2_) := (
  v^(n_?EvenQ) ^:= v2^(n/2);
  (v^(n_?powQ) ^:= v2^((n - Sign[n])/2) #^Sign[n])&[ v /. Pattern -> (#1 &) ];
  Square[v] = v2
)

Sq/: (Sq[v_] =.) := (
  v/: v^(n_?EvenQ) =.;
  v/: v^(n_?powQ) =.;
  Square[v] =.
)


(* definitions for the Standard Model *)

Sq[EL] = 4 Pi Alfa;
Sq[Alfa] = Alfa2;
Alfa2/: Alfa2/Alfa = Alfa

Sq[GS] = 4 Pi Alfas;
Sq[Alfas] = Alfas2;
Alfas2/: Alfas2/Alfas = Alfas

Sq[SW] = SW2;
Sq[CW] = CW2;
CW2/: CW2 - 1 = -SW2;
SW2/: SW2 - 1 = -CW2;
CW2/: CW2 + SW2 = 1

Sq[MZ] = MZ2;
Sq[MW] = MW2;
Sq[MH] = MH2;

Sq[ME] = ME2;  Sq[MM] = MM2;  Sq[ML] = ML2;
Sq[MU] = MU2;  Sq[MC] = MC2;  Sq[MT] = MT2;
Sq[MD] = MD2;  Sq[MS] = MS2;  Sq[MB] = MB2

MLE[a__] := Mf[2,a];
MQU[a__] := Mf[3,a];
MQD[a__] := Mf[4,a];
Sq[Mf[a__]] = Mf2[a]

Mf[2,1] = ME;  Mf[2,2] = MM;  Mf[2,3] = ML;
Mf[3,1] = MU;  Mf[3,2] = MC;  Mf[3,3] = MT;
Mf[4,1] = MD;  Mf[4,2] = MS;  Mf[4,3] = MB

Mf2[2,1] = ME2;  Mf2[2,2] = MM2;  Mf2[2,3] = ML2;
Mf2[3,1] = MU2;  Mf2[3,2] = MC2;  Mf2[3,3] = MT2;
Mf2[4,1] = MD2;  Mf2[4,2] = MS2;  Mf2[4,3] = MB2

Conjugate[CKM[a__]] ^:= CKMC[a];
Conjugate[CKMC[a__]] ^:= CKM[a]

(* these symbols represent real quantities, i.e. Conjugate[sym] = sym
   for any of these.  Thinking e.g. of complex masses this looks
   dangerous but then again it's easy to remove any such definition.
   The function that really needs this is SquaredME. *)

Scan[ (RealQ[#] = True)&,
  { EL, Alfa, Alfa2, GS, Alfas, Alfas2,
    SW, CW, SW2, CW2,
    MW, MW2, MZ, MZ2,
    MH, MH2, MG0, MG02, MGp, MGp2,
    ME, ME2, MM, MM2, ML, ML2,
    MU, MU2, MC, MC2, MT, MT2,
    MD, MD2, MS, MS2, MB, MB2, _Mf, _Mf2 } ]

(* Model parameters which are defined using the parameter statement in
   Fortran (i.e. as numeric constants; see model.h) are given some
   numeric value here.  Using this information, the OptTimes function can
   significantly optimize the generated Fortran code.  The idea is to put
   everything that is known as constant at compile time in one place,
   i.e. rearrange products such that they are of the form (const)*(vars),
   then the compiler will usually collect all of these constants into one
   number. *)

Scan[ (N[#] = Random[])&,
  { cI, Alfa, Alfa2, SW2, CW, CW2,
    MW, MW2, MZ, MZ2,
    ME, ME2, MM, MM2, ML, ML2,
    MU, MU2, MC, MC2, MT, MT2,
    MD, MD2, MS, MS2 } ]

SMReduce[foo_][expr_, r___] := foo[expr /.
  SW2 -> 1 - CW2 /.
  {CW -> MW/MZ, CW2 -> MW2/MZ2}, r]

SMShorten[foo_][x__] := SMReduce[foo][x] /.
  MW2 - MZ2 -> -SW2 MZ2 /.
  {MZ/MW -> 1/CW, MZ2/MW2 -> 1/CW2,
   MW/MZ -> CW, MW2/MZ2 -> CW2}

SMSimplify = SMShorten[Simplify];
SMFullSimplify = SMShorten[FullSimplify]

(* definitions for the MSSM *)

If[ ValueQ[$FormCalc],
  SetOptions[ CalcFeynAmp,
    NoExpand -> {USf, USfC, UASf, UASfC,
      UCha, UChaC, VCha, VChaC, ZNeu, ZNeuC,
      UHiggs, UHiggsC, ZHiggs, ZHiggsC},
    NoBracket -> {CKM, CKMC, USf, USfC, UASf, UASfC,
      UCha, UChaC, VCha, VChaC, ZNeu, ZNeuC,
      UHiggs, UHiggsC, ZHiggs, ZHiggsC} ]
]

Af[t_, g_] := Af[t, g, g]

USf[t_, g_][a_, b_] := USf[a, b, t, g];
UASf[t_][a_, b_] := UASf[a, b, t]

Conjugate[USf[a__]] ^:= USfC[a];
Conjugate[USfC[a__]] ^:= USf[a]

Conjugate[UASf[a__]] ^:= UASfC[a];
Conjugate[UASfC[a__]] ^:= UASf[a]

Conjugate[UCha[a__]] ^:= UChaC[a];
Conjugate[UChaC[a__]] ^:= UCha[a]

Conjugate[VCha[a__]] ^:= VChaC[a];
Conjugate[VChaC[a__]] ^:= VCha[a]

Conjugate[ZNeu[a__]] ^:= ZNeuC[a];
Conjugate[ZNeuC[a__]] ^:= ZNeu[a]

Conjugate[UHiggs[a__]] ^:= UHiggsC[a];
Conjugate[UHiggsC[a__]] ^:= UHiggs[a]

Conjugate[ZHiggs[a__]] ^:= ZHiggsC[a];
Conjugate[ZHiggsC[a__]] ^:= ZHiggs[a]

Conjugate[Af[a__]] ^:= AfC[a];
Conjugate[AfC[a__]] ^:= Af[a]

Conjugate[MUE] ^:= MUEC;
Conjugate[MUEC] ^:= MUE

Conjugate[Mino3] ^:= Mino3C;
Conjugate[Mino3C] ^:= Mino3

Conjugate[SqrtEGl] ^:= SqrtEGlC;
Conjugate[SqrtEGlC] ^:= SqrtEGl

SqrtEGl/: SqrtEGl SqrtEGlC = 1

Sq[SqrtEGl] = Mino3/MGl;
Sq[SqrtEGlC] = Mino3C/MGl

Sq[SA] = SA2;
Sq[CA] = CA2;
CA2/: CA2 + SA2 = 1

Sq[TB] = TB2;
Sq[SB] = SB2;
Sq[CB] = CB2;
CB2/: CB2 + SB2 = 1;
CB2/: CB2 - 1 = -SB2;
SB2/: SB2 - 1 = -CB2;
CB2/: CB2 TB2 = SB2;
CB2/: CB2 TB = S2B/2;
SB2/: SB2/TB2 = CB2;
CB/: CB TB = SB;
CB/: CB SB = S2B/2;
CB2/: CB2 SB2 = S2B^2/4;
CB2/: CB2 - SB2 = C2B;
SB2/: SB2 - CB2 = -C2B;
SB2/: SB2/SB = SB;
SB2/: SB2/S2B = TB/2;
CB2/: CB2/CB = CB;
CB2/: CB2/S2B = 1/(2 TB);
S2B/: S2B/SB = 2 CB;
S2B/: S2B/CB = 2 SB;
S2B/: S2B TB = 2 SB2

Sq[CBA] = CBA2;
Sq[SBA] = SBA2;
CBA2/: CBA2 + SBA2 = 1

CA/: CA SA = S2A/2;
CA2/: CA2 - SA2 = C2A;
SA2/: SA2 - CA2 = -C2A;

Sq[MGl] = MGl2;
Sq[MSf[a__]] = MSf2[a];
Sq[MASf[a__]] = MASf2[a];
Sq[MCha[a__]] = MCha2[a];
Sq[MNeu[a__]] = MNeu2[a]

Sq[MHiggs[a__]] = MHiggs2[a];
Sq[MHiggstree[a__]] = MHiggstree2[a]

Sq[Mh0] = Mh02;  Sq[Mh0tree] = Mh0tree2;
Sq[MHH] = MHH2;  Sq[MHHtree] = MHHtree2;
Sq[MA0] = MA02;  Sq[MA0tree] = MA0tree2;
Sq[MHp] = MHp2;  Sq[MHptree] = MHptree2

Scan[ (RealQ[#] = True)&,
  { TB, CB, SB, TB2, CB2, SB2, C2B, S2B,
    CA, SA, CA2, SA2, C2A, S2A,
    CAB, SAB, CBA, SBA, CBA2, SBA2,
    Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2, MGl, MGl2,
    _MSf, _MSf2, _MCha, _MCha2, _MNeu, _MNeu2,
    _MHiggs, _MHiggstree } ]

MSSMReduce[foo_, red_:SMReduce][expr_, r___] :=
  red[foo][expr /. SBA2 -> 1 - CBA2, r] //.
  { CB2^2 + SB2^2 + c_ S2B^2 :> C2B^2 + (c + 1/2) S2B^2 /; c < 0,
    CB2^2 + SB2^2 + c_. S2B^2 :> 1 + (c - 1/2) S2B^2 } /.
  { 1 - S2B^2 -> C2B^2, -1 + S2B^2 -> -C2B^2,
    1 - C2B^2 -> S2B^2, -1 + C2B^2 -> -S2B^2 }

MSSMShorten[foo_, red_:SMSimplify][x__] :=
  MSSMReduce[foo, red][x] /. CBA2 - 1 -> -SBA2

MSSMSimplify = MSSMShorten[Simplify];
MSSMFullSimplify = MSSMShorten[FullSimplify]

MSSMTrig[expr_, simp_:Simplify] := SUSYTrigFullReduce @
  MapOnly[Plus, simp, _ | a | b] @ SUSYTrigFullExpand[expr]

MSSMTrig2[expr_, simp_:Simplify] := MapOnly[Plus,
  SUSYTrigFullSimplify[#, simp]&,
  _ | CB | SB | CB2 | SB2 | C2B | S2B | TB |
  CA | SA | CA2 | SA2 | C2A | S2A |
  CAB | CBA | CBA2 | SBA2] @ expr


SUSYTrigFullExpand[expr_] := expr /. {
  SB -> Sin[b], CB -> Cos[b], SB2 -> Sin[b]^2, CB2 -> Cos[b]^2,
  TB -> Tan[b], TB2 -> Tan[b]^2,
  SA -> Sin[a], CA -> Cos[a], SA2 -> Sin[a]^2, CA2 -> Cos[a]^2,
  C2B -> Cos[2 b], S2B -> Sin[2 b],
  C2A -> Cos[2 a], S2A -> Sin[2 a],
  CAB -> Cos[a + b], SAB -> Sin[a + b],
  CBA -> Cos[b - a], SBA -> Sin[b - a],
  CBA2 -> Cos[b - a]^2, SBA2 -> Sin[b - a]^2 }

SUSYTrigFullReduce[expr_] := expr //. {
  Cos[2 b] -> C2B, Sin[2 b] -> S2B,
  Sec[2 b] -> 1/C2B, Csc[2 b] -> 1/S2B,
  Cos[2 a] -> C2A, Sin[2 a] -> S2A,
  Sec[2 a] -> 1/C2A, Csc[2 a] -> 1/S2A,
  Cos[n_?EvenQ x:a | b] :> 1 - 2 Sin[n/2 x]^2,
  Sin[n_?EvenQ x:a | b] :> 2 Cos[n/2 x] Sin[n/2 x] } /. {
  Cos[b] -> CB, Sin[b] -> SB, Tan[b] -> TB,
  Sec[b] -> 1/CB, Csc[b] -> 1/SB,
  Cos[a] -> CA, Sin[a] -> SA, Tan[a] -> SA/CA,
  Sec[a] -> 1/CA, Csc[a] -> 1/SA }

SUSYTrigFullSimplify[expr_, simp_:Simplify] :=
  SUSYTrigFullReduce[simp[SUSYTrigFullExpand[expr]]]

SUSYTrigExpand[expr_] := expr /. {
  SB -> sb, CB -> cb, SB2 -> sb^2, CB2 -> cb^2,
  TB -> sb/cb, TB2 -> sb^2/cb^2,
  SA -> sa, CA -> ca, SA2 -> sa^2, CA2 -> ca^2,
  C2B -> cb^2 - sb^2, S2B -> 2 cb sb,
  C2A -> ca^2 - sa^2, S2A -> 2 ca sa,
  CAB -> ca cb - sa sb, SAB -> cb sa + ca sb,
  CBA -> ca cb + sa sb, SBA -> ca sb - cb sa,
  CBA2 -> (ca cb + sa sb)^2, SBA2 -> (ca sb - cb sa)^2 }

SUSYTrigReduce[expr_] := expr /.
  {ca -> CA, sa -> SA, cb -> CB, sb -> SB}

SUSYTrigSimplify[expr_, simp_:Simplify] :=
  SUSYTrigReduce[simp[SUSYTrigExpand[expr]]]


MassDim0 = { EL, Alfa, Alfa2, GS, Alfas, Alfas2,
  SW, CW, SW2, CW2,
  TB, CB, SB, TB2, CB2, SB2, C2B, S2B,
  CA, SA, CA2, SA2, C2A, S2A,
  CAB, SAB, CBA, SBA, CBA2, SBA2,
  SqrtEGl, SqrtEGlC,
  _CKM, _CKMC,
  _USf, _USfC, _UASf, _UASfC,
  _UCha, _UChaC, _VCha, _VChaC, _ZNeu, _ZNeuC,
  _UHiggs, _UHiggsC, _ZHiggs, _ZHiggsC }

MassDim1 = { MW, MZ, MH, MG0, MGp,
  ME, MM, ML, MU, MC, MT, MD, MS, MB, _Mf,
  Mh0, MHH, MA0, MHp, MGl,
  MUE, MUEC, Mino3, Mino3C,
  _MSf, _MASf, _MCha, _MNeu, _Af, _AfC }

MassDim2 = { MW2, MZ2, MH2, MG02, MGp2,
  ME2, MM2, ML2, MU2, MC2, MT2, MD2, MS2, MB2, _Mf2,
  Mh02, MHH2, MA02, MHp2, MGl2,
  _MSf2, _MASf2, _MCha2, _MNeu2 }

MassDim[expr_] := expr /.
  (# -> Random[] &)/@ MassDim0 /.
  (# -> Mass Random[] &)/@ MassDim1 /.
  (# -> Mass^2 Random[] &)/@ MassDim2

