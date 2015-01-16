(*
	softbrems.m
		calculates the soft-photon bremsstrahlung factor
		which multiplies the Born matrix element
		last modified 26 Aug 99 th

Reference:
	Ansgar Denner, "Techniques for the calculation of electroweak
	radiative corrections at one-loop level and results for
	W-physics at LEP200", Fortschr. d. Physik, 41 (1993) 4
*)

(* you have to specify these here: *)

Mass2[1] = MW2; Charge[1] = -1;
Mass2[2] = MW2; Charge[2] = 1;
Mass2[3] = MW2; Charge[3] = -1;
Mass2[4] = MW2; Charge[4] = 1;

(* you shouldn't need to change things below this point *)

k[1][[i_]] ^:= {EE[1], 0, 0, -Pin}[[i]];
k[2][[i_]] ^:= {EE[2], 0, 0, Pin}[[i]];
k[3][[i_]] ^:= {EE[3], -Pout Sin[th], 0, -Pout Cos[th]}[[i]];
k[4][[i_]] ^:= {EE[4], Pout Sin[th], 0, Pout Cos[th]}[[i]];

If[ Mass2[3] === Mass2[4], EE[4] = EE[3] ];
If[ Mass2[1] === Mass2[2], EE[2] = EE[1] ];
If[ Mass2[1] === Mass2[2] === Mass2[3] === Mass2[4],
  EE[4] = EE[3] = EE[2] = EE[1];
  Pout = Pin ];

PP[1 | 2] = Pin;
PP[3 | 4] = Pout;

Attributes[Pair] = {Orderless};
Pair[i_, i_] = Mass2[i];
Pair[1, 2] = S/2 - Mass2[1]/2 - Mass2[2]/2;
Pair[3, 4] = S/2 - Mass2[3]/2 - Mass2[4]/2;
Pair[1, 3] = -T/2 + Mass2[1]/2 + Mass2[3]/2;
Pair[2, 4] = -T/2 + Mass2[2]/2 + Mass2[4]/2;
Pair[2, 3] = -U/2 + Mass2[2]/2 + Mass2[3]/2;
Pair[1, 4] = -U/2 + Mass2[1]/2 + Mass2[4]/2

pair3[a_, b_] := Sum[k[a][[i]] k[b][[i]], {i, 2, 4}];
pair[a_, b_] := k[a][[1]] k[b][[1]] - pair3[a, b]

Block[ {Set},
  (s_Symbol = v_) := (N[s] := N[v]; N[v]);
  << para.m;
];
N[MH2] = (N[MH] = 100)^2;
N[S] = Pair[1, 1] + Pair[2, 2] + 2 pair[1, 2]//Simplify//N;
N[T] = Pair[1, 1] + Pair[3, 3] - 2 pair[1, 3]//Simplify//N;
N[U] = Pair[1, 1] + Pair[4, 4] - 2 pair[1, 4]//Simplify//N;
N[EE[i_]] := xEE[i//Floor];
N[Pin] = N[Pout] = 300;
N[th] = Pi/3;
N[M1] = 101;
N[M2] = 102;
N[M3] = 103;
N[M4] = 104


xEE[i_] := Sqrt[PP[i]^2 + Mass2[i]];
xPin = PowerExpand[
  Sqrt[(S + Mass2[2] - Mass2[1])^2/4/S - Mass2[2]]//Simplify, S];
xPout = PowerExpand[
  Sqrt[(S + Mass2[4] - Mass2[3])^2/4/S - Mass2[4]]//Simplify, S]


(* Choose the physical root of the two cases for alpha. This is done
   by simply plugging in numbers for the parameters -- approximate ones
   suffice -- and taking the positive root. The energy scale should be
   larger than the largest mass in the model so that kinematical variables
   don't pick up an imaginary part. Assuming a light Higgs, 300 GeV for
   the SM is therefore ok. *)

Attributes[ChooseRoot] = {Orderless}

ChooseRoot[i_, j_] := ChooseRoot[i, j] =
Block[ {roots, n},
  roots = Solve[
    alpha^2 Pair[i, i] - 2 alpha Pair[i, j] + Pair[j, j] == 0,
    alpha ]//Simplify;
  n = N[(alpha k[i][[1]] - k[j][[1]])/k[j][[1]] /. roots];
  roots[[ If[n[[1]] < 0, 2, 1], 1, 2 ]]
]

cabs[u_] := PowerExpand[Sqrt[pair3[u, u]]//Simplify, {S, Pin, Pout}]

f[a_, u_] :=
Block[ {xm, xp},
  xm = k[u][[1]] - cabs[u];
  xp = k[u][[1]] + cabs[u];
  1/4 log[xm/xp]^2 + Li2[1 - a xp] + Li2[1 - a xm]
]


(* ESOFT - the maximum energy a soft photon may have without being
           detected,
   delta - the IR regulator (``photon mass'') _squared_ *)

Attributes[L] = {Orderless}

L[i_, i_] := L[i, i] =
  2 Pi (log[4 ESOFT^2/delta] +
    EE[i]/PP[i] log[(EE[i] - PP[i])/(EE[i] + PP[i])])

L[i_, j_] := L[i, j] =
Block[ {alpha, v, vnum},
  alpha = PowerExpand[ChooseRoot[i, j], S];
  vnum = alpha^2 Pair[i, i] - Pair[j, j]//Simplify;
  v = vnum/2/(alpha k[i][[1]] - k[j][[1]])//FullSimplify;
  Simplify[
    4 Pi alpha/vnum Pair[i, j] (
      1/2 log[alpha^2 Pair[i, i]/Pair[j, j]] log[4 ESOFT^2/delta] +
      f[alpha/v, i] - f[1/v, j] ) ]
]

InOut[1 | 2] = 1;
InOut[3 | 4] = -1

SoftPre := SoftPre = Sum[
  Print["L[", i, ", ", j, "] = ", L[i, j]//Shallow];
  InOut[i] InOut[j] Charge[i] Charge[j] L[i, j]/2,
  {i, 4}, {j, 4}]

(* uncomment the /. EE -> xEE ... if you want the soft photon factor
   in terms of S, T, U instead of EE[i] and Pin, Pout *)

SoftBrems := -4 Pi Alpha/(2 Pi)^3 *
  Collect[SoftPre (* /. EE -> xEE /. Pin -> xPin /. Pout -> xPout *),
    {Li2[__], log[__]}, Simplify]


Optimize[expr_] :=
Block[ {plus, abbr, res, qq, qqc = 0},
  Attributes[plus] = {Flat, Orderless};
  abbr[p__] := (abbr[p] = plus[p] = qq[++qqc]) /; FreeQ[{p}, abbr];
  res = expr /. Plus -> abbr;
  { Simplify/@ Sort[
      Cases[DownValues[plus], (_[_[p__]] :> q_) :> q -> Plus[p]] ],
    res } /. qq[i_] :> ToExpression["qq" <> ToString[i]]
]

Unprotect[Rule];
Format[a_ -> b_, FortranForm] := SequenceForm[a, " = ", b];
Protect[Rule]


Print["writing softphot.F"]

hh = OpenWrite[
  "!r8_" <> Environment["HOSTTYPE"] <> " > softphot.F",
  FormatType -> FortranForm, PageWidth -> 67]

WriteString[hh,
  "#define ESOFT YouHaventDefinedESOFTYet\n\n" <>
  "\tdouble precision function softphot()\n" <>
  "\timplicit double precision (q)\n" <>
  "\timplicit logical (a-p,r-z)\n" <>
  "#include \"kin.h\"\n\n" <>
  "\tdouble precision delta\n" <>
  "\tcommon /ffcut/ delta\n\n" <>
  "\tdouble precision Li2\n" <>
  "\texternal Li2\n\n"]

dSB = Optimize[SoftBrems] /.
  n_?NumberQ :> N[n] /.
  x_^r_Real :> x^Rationalize[r] /. -1. x_ -> -x /.
  EE[r_Real] :> EE[Floor[r]];
(WriteString[hh, "\t"]; Write[hh, #])&/@ dSB[[1]];
WriteString[hh, "\t"];
Write[hh, softphot -> dSB[[2]] ]

WriteString[hh, "\tend\n"];
Close[hh];

