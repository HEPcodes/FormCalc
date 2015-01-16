(*
	chkfinite.m
		last modified 2 Jul 99 th

These definitions replace the one-loop integrals by their
divergent parts. This is useful for checking UV and IR finiteness,
e.g. if you want to check UV finiteness of an expression expr,
load this file and see if Coefficient[expr, 1/EPS] == 0.

*)

(* UV divergent integrals *)

A0[m_] = 2 m/EPS;
B0[__] = 2/EPS;
B1[__] = -1/EPS;
B00[p_, m1_, m2_] = ((m1 + m2)/2 - p/6)/EPS;
B11[__] = 2/3/EPS;
C0i[cc00, __] = 1/2/EPS;
C0i[cc001 | cc002, __] = -1/6/EPS;
D0i[dd0000, __] = 1/12/EPS;
DB00[__] = -2/9/EPS;
DB11[p_, m1_, m2_] = (m1 - m2)/6/p^2/EPS

(* IR divergent integrals *)

(* not tested yet
DB0[m_, 0, m_] = DB0[m_, m_, 0] =
  DB1[m_, m_, 0] =
  C0i[cc0, m_, 0, m_, 0, m_, m_] = -(1 + Log[LAMBDA^2/m]/2)/m;
C0i[cc0, m1_, p_, m2_, 0, m1_, m2_] =
  C0i[cc0, p_, m2_, m1_, m1_, m2_, 0] =
  C0i[cc0, m2_, m1_, p_, m2_, 0, m1_] =
Block[ {pij = (mj - mi - p)/2, beta = Sqrt[1 - mi mj/pij^2]},
  -Log[LAMBDA^2]/4/pij Log[(1 + beta)/(1 - beta)]/beta ];
*)

C0i[__] = D0i[__] = DB0[__] = DB1[__] = 0;

C2 = CW^-2;
C4 = CW^-4;
S2 = SW^-2;
S4 = SW^-4;
SW = Sqrt[1 - CW^2];
a2 = EL^4/16/Pi^2;
MW2 = MZ2 CW^2

