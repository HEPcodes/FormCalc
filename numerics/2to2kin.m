(*
	2to2kin.m
		Various kinematical variables for a 2 -> 2 process,
		uses same conventions as num.F
		this file is part of FormCalc
		last modified 23 Jun 00 th
*)


PowSimp[expr_] := PowerExpand[Simplify[expr]]

Pin2 = (S + Mass[2]^2 - Mass[1]^2)^2/4/S - Mass[2]^2;
Pout2 = (S + Mass[4]^2 - Mass[3]^2)^2/4/S - Mass[4]^2;
Pin = Sqrt[Pin2];
Pout = Sqrt[Pout2]

EE[1] = Sqrt[Pin2 + Mass[1]^2]//PowSimp;
EE[2] = Sqrt[Pin2 + Mass[2]^2]//PowSimp;
EE[3] = Sqrt[Pout2 + Mass[3]^2]//PowSimp;
EE[4] = Sqrt[Pout2 + Mass[4]^2]//PowSimp

(*
T = Mass[1]^2 + Mass[3]^2 - 2 (EE[1] EE[3] - Pin Pout Cos[th]);
U = Mass[2]^2 + Mass[4]^2 - 2 (EE[2] EE[4] + Pin Pout Cos[th])
or, if you prefer:
T = MomSquare[Vec[k[1]] - Vec[k[3]]];
U = MomSquare[Vec[k[1]] - Vec[k[4]]]
*)


Vec[ Conjugate[x_] ] := Conjugate[Vec[x]] /; ValueQ[Vec[x]]

Vec[ k[1]     ] = {EE[1], 0, 0, -Pin};
Vec[ ep[1][0] ] = {-Pin, 0, 0, EE[1]}/Mass[1];
Vec[ ep[1][1] ] = {0, -1, I, 0}/Sqrt[2];
Vec[ ep[1][2] ] = {0, -1, -I, 0}/Sqrt[2]

Vec[ k[2]     ] = {EE[2], 0, 0, Pin};
Vec[ ep[2][0] ] = {-Pin, 0, 0, -EE[2]}/Mass[2];
Vec[ ep[2][1] ] = {0, 1, I, 0}/Sqrt[2];
Vec[ ep[2][2] ] = {0, 1, -I, 0}/Sqrt[2]

Vec[ k[3]     ] = {EE[3], -Pout Sin[th], 0, -Pout Cos[th]};
Vec[ ep[3][0] ] = {Pout, -EE[3] Sin[th], 0, -EE[3] Cos[th]}/Mass[3];
Vec[ ep[3][1] ] = {0, -Cos[th], -I, Sin[th]}/Sqrt[2];
Vec[ ep[3][2] ] = {0, -Cos[th], I, Sin[th]}/Sqrt[2]

Vec[ k[4]     ] = {EE[4], Pout Sin[th], 0, Pout Cos[th]};
Vec[ ep[4][0] ] = {Pout, EE[4] Sin[th], 0, EE[4] Cos[th]}/Mass[4];
Vec[ ep[4][1] ] = {0, Cos[th], -I, -Sin[th]}/Sqrt[2];
Vec[ ep[4][2] ] = {0, Cos[th], I, -Sin[th]}/Sqrt[2]


KinSimplify[expr_] := Simplify[expr] /; FreeQ[expr, _^_Rational]

KinSimplify[expr_] :=
Block[ {sym, exp = Simplify[expr]},
  sym = Complement[
    Cases[Cases[exp, r_^_Rational -> r, Infinity],
      x_Symbol /; Context[x] =!= "System`", {-1}],
    {S, Mass[1], Mass[2], Mass[3], Mass[4], Mass} ];
(* assuming that S and the masses are positive real quantities: *)
  If[ Length[sym] === 0, PowSimp[exp], exp ]
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

Pair[a_, x_Plus] := Block[ {Pair}, Distribute[Pair[a, x]] ]
          
Pair[a_, n_?NumberQ x_] := n Pair[a, x]
      
Pair[0, _] = 0

Pair[a_, b_] := Pair[Vec[a], Vec[b]] /; ValueQ[Vec[a]] && ValueQ[Vec[b]]

Pair[a_List, b_List] :=
  Pair[a, b] = KinSimplify[ a[[1]] b[[1]] - Pair3[a, b] ]


Pair3[a_, b_] := Pair3[Vec[a], Vec[b]] /; ValueQ[Vec[a]] && ValueQ[Vec[b]]

Pair3[a_List, b_List] := Pair3[a, b] = Take[a, -3] . Take[b, -3]


MomSquare[a_] := Pair[a, a]


Eps[a_, b_, c_, d_] := Eps[Vec[a], Vec[b], Vec[c], Vec[d]] /;
  ValueQ[Vec[a]] && ValueQ[Vec[b]] && ValueQ[Vec[c]] && ValueQ[Vec[d]]

Eps[a_List, b_List, c_List, d_List] := Eps[a, b, c, d] =
  KinSimplify[ I Det[{a, b, c, d}] ]

(* Note: although Eps is defined as -I*(Levi-Civita tensor) in FormCalc,
   we have to multiply with I (not -I) here because Det uses the
   Euclidean metric. *)

