(*
	USfSimplify.m
		simplify the sfermion mixing matrices
		using unitarity as much as possible
		this file is part of FormCalc
		last modified 23 Feb 15 th
*)


Attributes[USfMatch] = {Orderless, Flat}

USfMatch[USf[a__], USfC[a__]] := UCSf[a]

USfMatch[USf[1,j_, r___], USfC[2,j_, r___]] := UCSf[3,j, r];
USfMatch[USfC[1,j_, r___], USf[2,j_, r___]] := UCSfC[3,j, r]

USfMatch[USf[i_,1, r___], USfC[i_,2, r___]] := UCSf[i,3, r];
USfMatch[USfC[i_,1, r___], USf[i_,2, r___]] := UCSfC[i,3, r]

USfMatch[USf[1,1, r___], USfC[2,2, r___]] := UCSf[3,3, r];
USfMatch[USfC[1,1, r___], USf[2,2, r___]] := UCSfC[3,3, r];
USfMatch[USf[1,2, r___], USfC[2,1, r___]] := UCSf[3,4, r];
USfMatch[USfC[1,2, r___], USf[2,1, r___]] := UCSfC[3,4, r]

USfMatch[USf[a__], USf[a__]] := UUSf[a];
USfMatch[USfC[a__], USfC[a__]] := UUSfC[a];

USfMatch[USf[1,j_, r___], USf[2,j_, r___]] := UUSf[3,j, r];
USfMatch[USfC[1,j_, r___], USfC[2,j_, r___]] := UUSfC[3,j, r]

USfMatch[USf[i_,1, r___], USf[i_,2, r___]] := UUSf[i,3, r];
USfMatch[USfC[i_,1, r___], USfC[i_,2, r___]] := UUSfC[i,3, r]

USfMatch[USf[1,1, r___], USf[2,2, r___]] := UUSf[3,3, r];
USfMatch[USfC[1,1, r___], USfC[2,2, r___]] := UUSfC[3,3, r];
USfMatch[USf[1,2, r___], USf[2,1, r___]] := UUSf[3,4, r];
USfMatch[USfC[1,2, r___], USfC[2,1, r___]] := UUSfC[3,4, r]


USfProd/: USfProd[i__] USfProd[j__] := USfProd[i, j]


USfSimplify[expr_] := expr /.
  { Abs[USf[a__]]^2 :> USf[a] USfC[a],
    Abs[USf[a__]]^-2 :> 1/UCSf[a] } /.
  { u_USf^n_.  :> USfProd@@ Table[u, {n}],
    u_USfC^n_. :> USfProd@@ Table[u, {n}] } /.
  p_Plus u:USfProd[_] :> Distribute[p u] /; !FreeQ[p, USfProd[_]] /.
  USfProd -> USfMatch /.
  USfMatch -> Times


i1[1] = 1;	i2[1] = 1;
i1[2] = 2;	i2[2] = 2;
i1[3] = 1;	i2[3] = 2;
i1[4] = 2;	i2[4] = 1

USfRev[UCSf[i_,j_, r___]] := USf[i1[i],i1[j], r] USfC[i2[i],i2[j], r];
USfRev[UCSfC[i_,j_, r___]] := USfC[i1[i],i1[j], r] USf[i2[i],i2[j], r];
USfRev[UUSf[i_,j_, r___]] := USf[i1[i],i1[j], r] USf[i2[i],i2[j], r];
USfRev[UUSfC[i_,j_, r___]] := USfC[i1[i],i1[j], r] USfC[i2[i],i2[j], r]

usimp[a_?AtomQ] := a

usimp[p_Plus] := usimp/@ p

usimp[p_ x_.] := p usimp[x] /;
  FreeQ[p, UUSf | UUSfC (*| UCSf | UCSfC*)]

usimp[x_] := usimp[x] = Simplify @ USfSimplify @
  Distribute[x /. u:_UCSf | _UCSfC | _UUSf | _UUSfC :> USfRev[u]]

USfSimplifyAll[expr_] := usimp[expr]


(* unitarity rules:
   1a   UCSf[1, 1] + UCSf[1, 2] == 1
   1b   UCSf[2, 1] + UCSf[2, 2] == 1
   1c   UCSf[1, 1] + UCSf[2, 1] == 1
   1d   UCSf[1, 2] + UCSf[2, 2] == 1

   2a   UCSf[3, 1] + UCSf[3, 2] == 0
   2a*  UCSfC[3, 1] + UCSfC[3, 2] == 0
   2b   UCSf[1, 3] + UCSf[2, 3] == 0
   2b*  UCSfC[1, 3] + UCSfC[2, 3] == 0
*)

(* 1a - 1c *)
UCSf[2,1, r___] := UCSf[1,2, r]

(* 1a - 1d *)
UCSf[2,2, r___] := UCSf[1,1, r]

(* 1a *)
UCSf/: UCSf[1,j1:1|2, r___] + n_. UCSf[1,j2:1|2, r___] :=
  1 + (n - 1) UCSf[1,j2, r] /; j1 != j2 && n > 0;
UCSf/: UCSf[1,j:1|2, r___] - 1 := -UCSf[1,3-j, r];

(* 2a *)
UCSf[3,2, r___] := -UCSf[3,1, r];
UCSfC[3,2, r___] := -UCSfC[3,1, r];
(* 2b *)
UCSf[2,3, r___] := -UCSf[1,3, r];
UCSfC[2,3, r___] := -UCSfC[1,3, r]

Attributes[USfChk] = {HoldAll};
USfChk[lhs_ := rhs_] := USfChk[lhs] == USfChk[rhs];
USfChk[lhs_] := USfSimplify @
Block[ {UCSf, UCSfC, UUSf, UUSfC},
  Distribute[lhs /. Pattern -> (#1&) /.
    u:_UCSf | _UCSfC | _UUSf | _UUSfC :> USfRev[u]]
];

UCSf/: UCSf[3,j:1|2, r___] UCSfC[3,j_, r___] := UCSf[1,j, r] UCSf[2,j, r];
UCSf/: UCSf[i:1|2,3, r___] UCSfC[i_,3, r___] := UCSf[i,1, r] UCSf[i,2, r];
UUSf/: UUSf[3,j:1|2, r___] UUSfC[3,j_, r___] := UCSf[1,j, r] UCSf[2,j, r];
UUSf/: UUSf[i:1|2,3, r___] UUSfC[i_,3, r___] := UCSf[i,1, r] UCSf[i,2, r]

UCSf/: UCSf[3,3, r___] UCSfC[3,3, r___] := UCSf[1,1, r] UCSf[2,2, r];
UCSf/: UCSf[3,4, r___] UCSfC[3,4, r___] := UCSf[1,2, r] UCSf[2,1, r];

UCSf/: UCSf[3,3, r___] UCSfC[3,4, r___] := UCSf[1,3, r] UCSf[2,3, r];
UCSf/: UCSfC[3,3, r___] UCSf[3,4, r___] := UCSfC[1,3, r] UCSfC[2,3, r];
UCSf/: UCSf[3,3, r___] UCSf[3,4, r___] := UCSf[3,1, r] UCSf[3,2, r];
UCSfC/: UCSfC[3,3, r___] UCSfC[3,4, r___] := UCSfC[3,1, r] UCSfC[3,2, r];

UCSf/: UCSf[3,3, r___] UCSfC[1,3, r___] := UCSf[1,1, r] UCSf[3,2, r];
UCSf/: UCSf[1,3, r___] UCSfC[3,3, r___] := UCSf[1,1, r] UCSfC[3,2, r];
UCSf/: UCSf[3,3, r___] UCSfC[3,1, r___] := UCSf[1,1, r] UCSf[2,3, r];
UCSf/: UCSf[3,1, r___] UCSfC[3,3, r___] := UCSf[1,1, r] UCSfC[2,3, r];

UCSf/: UCSf[1,3, r___] UCSf[3,4, r___] := UCSf[1,2, r] UCSf[3,1, r];
UCSfC/: UCSfC[1,3, r___] UCSfC[3,4, r___] := UCSf[1,2, r] UCSfC[3,1, r];

UCSf/: UCSf[3,4, r___] UCSfC[3,1, r___] := UCSf[2,1, r] UCSfC[1,3, r];
UCSf/: UCSf[3,1, r___] UCSfC[3,4, r___] := UCSf[2,1, r] UCSf[1,3, r];

(*
UUSf/: UUSf[1,3, r___] UUSfC[3,3, r___] := UCSf[1,1, r] UCSf[3,2, r];
UUSf/: UUSfC[1,3, r___] UUSf[3,3, r___] := UCSf[1,1, r] UCSfC[3,2, r];
UUSf/: UUSf[1,3, r___] UUSfC[3,4, r___] := UCSf[1,2, r] UCSf[3,1, r];
UUSf/: UUSfC[1,3, r___] UUSf[3,4, r___] := UCSf[1,2, r] UCSfC[3,1, r];

UUSf/: UUSf[2,3, r___] UUSfC[3,3, r___] := UCSf[2,2, r] UCSfC[3,1, r];
UUSf/: UUSfC[2,3, r___] UUSf[3,3, r___] := UCSf[2,2, r] UCSf[3,1, r];
UUSf/: UUSf[2,3, r___] UUSfC[3,4, r___] := UCSf[2,1, r] UCSfC[3,2, r];
UUSf/: UUSfC[2,3, r___] UUSf[3,4, r___] := UCSf[2,1, r] UCSf[3,2, r];

UCSf/: UCSf[1,3, r___] UCSfC[3,1, r___] := UCSf[1,1, r] UCSfC[3,4, r];
*)

UUSf/: UUSf[i:1|2,j:1|2, r___] UUSfC[i_,j_, r___] := UCSf[i,j, r]^2

UCSf/: UCSf[3,j:1|2, r___] UUSf[3,j_, r___] := UUSf[1,j, r] UCSf[2,j, r];
UCSfC/: UCSfC[3,j:1|2, r___] UUSfC[3,j_, r___] := UUSfC[1,j, r] UCSf[2,j, r]

UCSf/: UCSf[i:1|2,3, r___] UUSf[i_,3, r___] := UUSf[i,1, r] UCSf[i,2, r];
UCSfC/: UCSfC[i:1|2,3, r___] UUSfC[i_,3, r___] := UUSfC[i,1, r] UCSf[i,2, r]


Conjugate[UCSf[a___]] ^:= UCSfC[a];
Conjugate[UCSfC[a___]] ^:= UCSf[a];
Conjugate[UUSf[a___]] ^:= UUSfC[a];
Conjugate[UUSfC[a___]] ^:= UUSf[a];

