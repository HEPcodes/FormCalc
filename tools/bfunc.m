(*
	bfunc.m
		explicit formulas for the one- and two-point
		functions and their derivatives
		this file is part of FormCalc
		last modified 1 Apr 15 th

The formulas are more or less directly from the Passarino-Veltman
paper. The regularization parameters are
	Mudim  -- the mass scale squared
	Lambda -- the photon mass squared
	Delta  -- the infinity -2/(4 - D) - EulerGamma + Log[4 Pi]
*)


N[eps, ___] = 10^-50

N[Delta, ___] = 0

N[Mudim, ___] = 1

N[Lambda, ___] = 1


A0[0] = 0

A0[m_] := m (1 - Log[m/Mudim] + Delta)

(*
pvroots[p_, m1_, m2_] := pvroots[p, m1, m2] =
Block[ {x},
  (x /. Solve[-p x^2 + (p + m1 - m2) x - m1 == 0, x]) (*+ {I eps, -I eps}*)
]
*)

pvroots[p_, m1_, m2_] :=
  1/(2 p) {p + m1 - m2 + #, p + m1 - m2 - #}& @
    Sqrt[p (p - m1 - m2) - m1 (p - m1 + m2) - m2 (p + m1 - m2)]


(*pvf[n_, x_] := -x^n Log[(x - 1)/x] - Sum[x^(n - m)/m, {m, n}];*)

pvf[n_, 0] := Block[{x}, Limit[pvf[n, x], x -> 0]]

pvf[n_, x_] := LerchPhi[1/x, 1, n + 1]/x


pvB[n_, 0, m1_, m2_] :=
  (-1)^n (Log[m2/Mudim] + pvsimp[pvf[n, m1/(m1 - m2)]] - Delta)/n

pvB[n_, p_, m1_, m2_] :=
  ((-1)^n (Log[m2/Mudim] + pvsimp[pvf[n, #1] + pvf[n, #2]] - Delta)/n)&@@
  pvroots[p, m1, m2]

pvDB[n_, 0, m1_, m2_] :=
  (-1)^n (pvsimp[m2/(m1 - m2) pvf[n, m1/(m1 - m2)]] - 1/(n + 1))/(m1 - m2)

pvDB[n_, p_, m1_, m2_] :=
  ((-1)^n pvsimp[(1 - #1) pvf[n, #1] - (1 - #2) pvf[n, #2]]/p/(#1 - #2))&@@
  pvroots[p, m1, m2]


B0[0, 0, 0] = BAD	(* divergent, must cancel *)

B0[p_, 0, 0] := 2 - Log[-p/Mudim - I eps] + Delta

B0[p_, m1_, 0] := B0[p, 0, m1]

B0[0, m_, m_] := -Log[m/Mudim] + Delta

B0[p_, m1_, m2_] := pvB[1, p, m1, m2]


B1[p_, 0, 0] := -1/2 B0[p, 0, 0]

B1[p_, m1_, 0] := -B1[p, 0, m1] - B0[p, 0, m1]

B1[0, m_, m_] := (Log[m/Mudim] - Delta)/2

B1[p_, m1_, m2_] := pvB[2, p, m1, m2]


B11[p_, m1_, 0] := B11[p, 0, m1] -
  (A0[m1] - (m1 + 2 p) B0[p, 0, m1] - 4 p B1[p, 0, m1])/(3 p)

B11[p_, m1_, m2_] := pvB[3, p, m1, m2]


B111[p_, m1_, m2_] := pvB[4, p, m1, m2]


B00[p_, m1_, m2_] :=
  (3 (m1 + m2) - p)/18 + m1 B0[p, m1, m2]/3 +
  (A0[m2] + (m1 - m2 + p) B1[p, m1, m2])/6


Derivative[1, 0, 0][B0] = DB0

DB0[m_, 0, m_] := -(1 + Log[Lambda/m]/2)/m

DB0[m_, m_, 0] := -(1 + Log[Lambda/m]/2)/m

DB0[0, m_, m_] := 1/6/m

DB0[p_, m1_, m2_] := pvDB[1, p, m1, m2]


Derivative[1, 0, 0][B1] = DB1

DB1[m_, m_, 0] := (3 + Log[Lambda/m])/2/m

DB1[0, m_, m_] := -1/12/m

DB1[p_, m1_, m2_] := pvDB[2, p, m1, m2]


Derivative[1, 0, 0][B11] = DB11

DB11[p_, m1_, m2_] := pvB[3, p, m1, m2]


Derivative[1, 0, 0][B00] = DB00

DB00[p_, m1_, m2_] := -1/18 + m1 DB0[p, m1, m2]/3 + 
 (B1[p, m1, m2] + (p + m1 - m2) DB1[p, m1, m2])/6

(* note: DB00 looks like it could be IR divergent --
   it's not: the lambda-dependence cancels :-) *)


pvsimp[x_] := Simplify[x]

Null

