(*
	btensor.m
		explicit decompositions of the two-point
		tensor-coefficient functions
		this file is part of FormCalc
		last modified 25 Apr 00 th
*)


A0[0] = 0

B0[0, m_, m_] := A0[m]/m - 1 /; m =!= 0;
B0[0, m1_, m2_] := (A0[m2] - A0[m1])/(m2 - m1)

(*
B1[p_, m1_, m2_] =
  (m2 - m1) (B0[p, m1, m2] - B0[0, m1, m2])/(2 p) -
  B0[p, m1, m2]/2//Simplify
*)

B1[p_, m1_, m2_] =
  (A0[m1] - A0[m2] - (p + m2 - m1) B0[p, m1, m2])/(2 p)//Simplify

B00[p_, m1_, m2_] =
  (3 (m1 + m2) - p)/18 + m1 B0[p, m1, m2]/3 +
  (A0[m2] + (m1 - m2 + p) B1[p, m1, m2])/6//Simplify

B11[p_, m1_, m2_] =
  ((p - 3 (m1 + m2))/6 + A0[m2] - m1 B0[p, m1, m2] -
  2 (m1 - m2 + p) B1[p, m1, m2])/(3 p)//Simplify

Derivative[1, 0, 0][B0] = DB0;
Derivative[1, 0, 0][B1] = DB1;
Derivative[1, 0, 0][B00] = DB00;
Derivative[1, 0, 0][B11] = DB11

DB1[p_, m1_, m2_] = D[B1[p, m1, m2], p]//Simplify;
DB00[p_, m1_, m2_] = D[B00[p, m1, m2], p]//Simplify;
DB11[p_, m1_, m2_] = D[B11[p, m1, m2], p]//Simplify

B0[p_, m1_, m2_] := B0[p, m2, m1] /; !OrderedQ[{m1, m2}];
DB0[p_, m1_, m2_] := DB0[p, m2, m1] /; !OrderedQ[{m1, m2}]

