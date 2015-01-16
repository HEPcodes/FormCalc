(*
	2to2kin.m
		Various kinematical variables for a 2 -> 2 process,
		uses same conventions as num.F
		last modified 8 Sep 99 th
*)


k[1]     = {EE[1], 0, 0, -Pin};
ep[1][0] = 1/Mass[1] {-Pin, 0, 0, EE[1]};
ep[1][1] = {0, -1, 0, 0};
ep[1][2] = {0, 0, 1, 0};

k[2]     = {EE[2], 0, 0, Pin};
ep[2][0] = 1/Mass[2] {-Pin, 0, 0, -EE[2]};
ep[2][1] = {0, 1, 0, 0};
ep[2][2] = {0, 0, 1, 0};

k[3]     = {EE[3], -Pout Sin[th], 0, -Pout Cos[th]};
ep[3][0] = 1/Mass[3] {Pout, -EE[3] Sin[th], 0, -EE[3] Cos[th]};
ep[3][1] = {0, -Cos[th], 0, Sin[th]};
ep[3][2] = {0, 0, -1, 0};

k[4]     = {EE[4], Pout Sin[th], 0, Pout Cos[th]};
ep[4][0] = 1/Mass[4] {Pout, EE[4] Sin[th], 0, EE[4] Cos[th]};
ep[4][1] = {0, Cos[th], 0, -Sin[th]};
ep[4][2] = {0, 0, -1, 0};

Do[
  ep[i]["+"]=(ep[i][1] + I ep[i][2])/Sqrt[2];
  ep[i]["-"]=(ep[i][1] - I ep[i][2])/Sqrt[2],
{i, 4}];

Pin2 = (S + Mass[2]^2 - Mass[1]^2)^2/4/S - Mass[2]^2;
Pout2 = (S + Mass[4]^2 - Mass[3]^2)^2/4/S - Mass[4]^2;
Pin = Sqrt[Pin2];
Pout = Sqrt[Pout2];

EE[1] = Sqrt[Pin2 + Mass[1]^2];
EE[2] = Sqrt[Pin2 + Mass[2]^2];
EE[3] = Sqrt[Pout2 + Mass[3]^2];
EE[4] = Sqrt[Pout2 + Mass[4]^2];

T = Mass[1]^2 + Mass[3]^2 - 2 (EE[1] EE[3] - Pin Pout Cos[th]);
U = Mass[2]^2 + Mass[4]^2 - 2 (EE[2] EE[4] + Pin Pout Cos[th]);
(* or, if you prefer:
T = MomSquare[k[1] - k[3]];
U = MomSquare[k[1] - k[4]];
*)


KinSimplify[expr_] := Simplify[expr] /; FreeQ[expr, _^_Rational]

KinSimplify[expr_] :=
Block[ {sym, exp = Simplify[expr]},
  sym = Complement[
    Cases[Cases[exp, r_^_Rational -> r, Infinity],
      x_Symbol /; Context[x] =!= "System`", {-1}],
    {S, Mass[1], Mass[2], Mass[3], Mass[4], Mass} ];
(* assuming that S and the masses are positive real quantities: *)
  If[ Length[sym] === 0, Simplify[PowerExpand[exp]], exp ]
]

Pair[a_List, b_List] := Pair[a, b] = 
  KinSimplify[
    a[[1]] b[[1]] - a[[2]] b[[2]] - a[[3]] b[[3]] - a[[4]] b[[4]] ]

MomSquare[a_List] := Pair[a, a]


Eps[a_List, b_List, c_List, d_List] := Eps[a, b, c, d] =
  KinSimplify[ I Det[{a, b, c, d}] ]

(* Note: although Eps is defined as -I*(Levi-Civita tensor) in FormCalc,
   we have to multiply with I (not -I) here because Det uses the
   Euclidean metric. *)

