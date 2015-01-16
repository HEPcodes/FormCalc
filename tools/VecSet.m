(*
	VecSet.m
		explicit representation of external vectors
		and derived quantities
		this file is part of FormCalc
		last modified 23 Jan 14 th
*)


BeginPackage["VecSet`", "FormCalc`"]

VecSet::usage = "VecSet[i, m, p, e] defines the vectors for external
particle i with mass m and three-momentum p moving in direction
e = {ex, ey, ez}."

ToComponents::usage = "ToComponents[expr, pol] substitutes the vectors
appearing in expr by their components as set by VecSet.  The argument
pol specifies the polarizations of the external particles, they are
either given as a list of integers (-1 = left-handed, +1 = right-handed,
0 = longitudinal) or as a string of \"+\", \"-\", \"0\"."

SInvariant::usage = "SInvariant[i, j] computes the invariant (k_i + k_j)^2."

TInvariant::usage = "TInvariant[i, j] computes the invariant (k_i - k_j)^2."

ThreeMom::usage = "ThreeMom[sqrtS, ma, mb] computes the length of
\\vec p_b in the frame in which \\vec p_a + \\vec p_b vanishes, where
ma and mb are the masses corresponding respectively to p_a and p_b."


{ VecK, VecE, MomSpec, Spi,
  SpecM, SpecK, SpecE, SpecDeltaK, SpecEx, SpecEy, SpecEz,
  SpecPhi, SpecKT, SpecET, SpecKT, SpecPRap, SpecRap }


Begin["`Private`"]

cp[{r_, i_}] := r + I i

cm[{r_, i_}] := r - I i


safelog[n_, d_] := Infinity /; n == 0. || d == 0.

safelog[n_, d_] := Log[n/d]


(* i: the index of the momentum
   m, p: mass and three-momentum of the particle
   ex,ey,ez: the unit three-vector of the momentum *)

VecSet[i_, m_, p_, {ex_,ey_,ez_}, simp_:Simplify] :=
Block[ {p0, deltap, eT, sinth, expIphi, phi, onePez, oneMez},
  If[ TrueQ[m < 10^-7],
    p0 = p;
    deltap = 0,
  (* else *)
    p0 = Sqrt[p^2 + m^2];
    deltap = m^2/(p0 + p) ];

  eT = {ex, ey};
  sinth = eT.eT;
  onePez = 1 + ez;
  oneMez = 1 - ez;

  If[ TrueQ[sinth < 10^-14],
	(* phi is irrelevant when theta = 0 *)
    expIphi = {1, 0};
    phi = 0,
  (* else *)
    sinth = Sqrt[sinth];
    expIphi = eT/sinth;
    phi = ArcTan@@ eT ];

  VecK[i] = {{p0 onePez - deltap ez, p cm[eT]},
             {p cp[eT], p0 oneMez + deltap ez}} //simp;

  If[ m =!= 0,
    VecE[i, 0] = {{(p onePez + deltap ez)/m, p0/m cm[eT]},
                  {p0/m cp[eT], (p oneMez - deltap ez)/m}} //simp ];

  VecE[i, -1] = 1/Sqrt[2] {{-sinth, -oneMez cm[expIphi]},
                           {onePez cp[expIphi], sinth}} //simp;
  VecE[i, +1] = 1/Sqrt[2] {{-sinth, onePez cm[expIphi]},
                           {-oneMez cp[expIphi], sinth}} //simp;

  MomSpec[i] = {
    SpecM -> m,
    SpecK -> p,
    SpecE -> p0,
    SpecDeltaK -> deltap,
    SpecEx -> ex,
    SpecEy -> ey,
    SpecEz -> ez,
    SpecKT -> p sinth,
    SpecET -> Sqrt[(p sinth)^2 + m^2],
    SpecPhi -> phi,
    SpecPRap -> safelog[onePez, sinth],
    SpecRap -> safelog[VecK[i][[1,1]], Sqrt[(p sinth)^2 + m^2]]
  };

  (* this is E^[I phi] cos[th/2] = 1/Sqrt[2] Sqrt[1 + ez] expIphi: *)
  expIphi *= 1/Sqrt[2] Sqrt[onePez];

  (* this is sin[th/2]: *)
  sinth = 1/Sqrt[2] Sqrt[oneMez];

  sump = Sqrt[p0 + p];
  deltap = Sqrt[deltap];

  Spi[i, -1,6,1] = deltap {sinth, -cp[expIphi]} //simp; (* undotted *)
  Spi[i, -1,6,2] = deltap {sinth, -cm[expIphi]} //simp; (* dotted *)

  Spi[i, -1,7,1] = sump {sinth, -cp[expIphi]} //simp;
  Spi[i, -1,7,2] = sump {sinth, -cm[expIphi]} //simp;

  Spi[i, +1,6,1] = sump {cm[expIphi], sinth} //simp;
  Spi[i, +1,6,2] = sump {cp[expIphi], sinth} //simp;

  Spi[i, +1,7,1] = deltap {cm[expIphi], sinth} //simp;
  Spi[i, +1,7,2] = deltap {cp[expIphi], sinth} //simp;

  True
]


ToComponents[expr_, pol_String] :=
  ToComponents[expr, Characters[pol] /. {"+" -> 1, "-" -> -1, "0" -> 0}]

ToComponents[expr_, pol_List] :=
Block[ {ChainHead},
  expr /. WeylChain -> SplitChain /.
    Flatten[MapIndexed[vecrule, pol]] /.
    h:_ChainHead | Pair | Eps | Invariant :> rep[h] /.
    Hel[i_] :> pol[[i]]
]


vecrule[pol_, {n_}] := {
  k[n] -> VecK[n],
  e[n] -> VecE[n, pol],
  ec[n] -> VecE[n, -pol] }


inv[i_, j_, s_] :=
Block[ {vij = VecK[i][[##]] + s VecK[j][[##]] &},
  vij[1,1] vij[2,2] - vij[1,2] (vij[1,2] /. I -> -I)
]

SInvariant[i_Integer, j_Integer] := inv[i, j, +1]

TInvariant[i_Integer, j_Integer] := inv[i, j, -1]


ThreeMom[sqrtS_, ma_, mb_] := Sqrt[(# - mb) (# + mb)]&[
  1/2 (sqrtS - (ma - mb)*(ma + mb)/sqrtS) ]


Attributes[rep] = {HoldAll}

rep[Pair][A_List, B_List] :=
  1/2 (A[[1,1]] B[[2,2]] + A[[2,2]] B[[1,1]] -
       A[[1,2]] B[[2,1]] - A[[2,1]] B[[1,2]])

eps[A_, B_, C_, D_] :=
  (A[[1,1]] B[[2,2]] - A[[2,2]] B[[1,1]]) *
  (C[[2,1]] D[[1,2]] - C[[1,2]] D[[2,1]])

rep[Eps][A_List, B_List, C_List, D_List] := 1/4 (
  eps[A, B, C, D] + eps[C, D, A, B] -
  eps[A, C, B, D] - eps[B, D, A, C] +
  eps[A, D, B, C] + eps[B, C, A, D] )


eps[0] = {}

eps[1] = {{{0, 1}, {-1, 0}}}


bar[om_][{{A11_, A12_}, {A21_, A22_}}, {n_}] :=
  {{A22, -A12}, {-A21, A11}} /; EvenQ[om + n]

bar[_][A_, _] := A


rep[ChainHead[om_, n_]][
    Spinor[iL_, sL_, dL_], eL_, V___, eR_, Spinor[iR_, sR_, dR_] ] :=
  Level[
    { {Spi[iL, sL Hel[iL], 6 + BitAnd[om - 1 - eL, 1], dL]},
      eps[eL],
      MapIndexed[bar[om], {V}],
      eps[eR],
      {Spi[iR, sR Hel[iR], 6 + BitAnd[om + n + eR, 1], dR]} },
    {2}, Dot]


End[]

EndPackage[]

