(*
	kin2to2.m
		Various kinematical variables for a 2 -> 2 process,
		uses same conventions as num.F
		this file is part of FormCalc
		last modified 14 Jun 01 th
*)


Pin = Sqrt[(S - Mass[1]^2 + Mass[2]^2)^2/(4 S) - Mass[2]^2]

Pout = Sqrt[(S - Mass[3]^2 + Mass[4]^2)^2/(4 S) - Mass[4]^2]


Vec[ Conjugate[x_] ] := ComplexExpand[Conjugate[Vec[x]]] /; ValueQ[Vec[x]]

Four[x0_, {x3__}] := {x0, x3}

VecSet[n_, dir_, p_, er_] :=
Block[ {etheta, efi, sin2th},
  EE[n] = Sqrt[p^2 + Mass[n]^2];
  Vec[ k[n] ] = Four[EE[n], p er];
  Vec[ ep[n][0] ] =
    If[Mass[n] === 0, Undefined, Four[p, EE[n] er]/Mass[n]];
  sin2th = Simplify[1 - er[[3]]^2];
  Vec[ ep[n][1] ] = Four[0, 1/Sqrt[2] If[ sin2th === 0,
    {1, dir er[[3]] I, 0},
  (* else *)
    efi = {-er[[2]], er[[1]], 0}/Sqrt[sin2th];
    etheta = Cross[efi, er];
    etheta + dir I efi ]];
  Vec[ ep[n][2] ] := Conjugate[Vec[ep[n][1]]];
]

VecSet[1, 1, Pin, {0, 0, 1}];
VecSet[2, 1, Pin, -{0, 0, 1}];
VecSet[3, -1, Pout, {Sin[theta], 0, Cos[theta]}];
VecSet[4, -1, Pout, -{Sin[theta], 0, Cos[theta]}]


(*
T = MomSquare[Vec[k[1]] - Vec[k[3]]];
U = MomSquare[Vec[k[1]] - Vec[k[4]]]
*)


KinSimplify[expr_] := Simplify[expr] /; FreeQ[expr, _^_Rational]

KinSimplify[expr_] :=
Block[ {sym, exp = Simplify[expr]},
  sym = Complement[
    Cases[Cases[exp, r_^_Rational -> r, Infinity],
      x_Symbol /; Context[x] =!= "System`", {-1}],
    {S, Ecms, Mass[1], Mass[2], Mass[3], Mass[4], Mass} ];
(* assuming that S and the masses are positive real quantities: *)
  If[ Length[sym] === 0, PowerExpand[Simplify[exp]], exp ]
]


Attributes[Pair] = {Orderless}

Pair[k[i_], k[i_]] = Mass[i]^2

Pair[k[1], k[2]] =  S/2 - Mass[1]^2/2 - Mass[2]^2/2;
Pair[k[3], k[4]] =  S/2 - Mass[3]^2/2 - Mass[4]^2/2;
Pair[k[1], k[3]] = -T/2 + Mass[1]^2/2 + Mass[3]^2/2;
Pair[k[2], k[4]] = -T/2 + Mass[2]^2/2 + Mass[4]^2/2;
Pair[k[2], k[3]] = -U/2 + Mass[2]^2/2 + Mass[3]^2/2;
Pair[k[1], k[4]] = -U/2 + Mass[1]^2/2 + Mass[4]^2/2

Pair[k[i_], ep[i_][_]] = 0;
Pair[k[i_], e[i_]] = 0

Pair[a_, p_Plus] := Pair[a, #]&/@ p
          
Pair[a_, n_?NumberQ b_] := n Pair[a, b]
      
Pair[0, _] = 0

Pair[a_, b_] :=
  (Pair[a, b] = Pair[Vec[a], Vec[b]]) /; ValueQ[Vec[a]] && ValueQ[Vec[b]]

Pair[a_List, b_List] := KinSimplify[ a[[1]] b[[1]] - Pair3[a, b] ]


Pair3[a_, b_] := Pair3[Vec[a], Vec[b]] /; ValueQ[Vec[a]] && ValueQ[Vec[b]]

Pair3[a_List, b_List] := Take[a, -3] . Take[b, -3]


MomSquare[a_] := Pair[a, a]


Eps[a_, b_, c_, d_] :=
  (Eps[a, b, c, d] = Eps[Vec[a], Vec[b], Vec[c], Vec[d]]) /;
  ValueQ[Vec[a]] && ValueQ[Vec[b]] && ValueQ[Vec[c]] && ValueQ[Vec[d]]

Eps[a_List, b_List, c_List, d_List] :=
  KinSimplify[ I Det[{a, b, c, d}] ]

(* Note: although Eps is defined as -I*(Levi-Civita tensor) in FormCalc,
   we have to multiply with I (not -I) here because Det uses the
   Euclidean metric. *)

