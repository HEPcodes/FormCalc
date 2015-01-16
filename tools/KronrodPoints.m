(*
        KronrodPoints.m
                calculate the abscissas and weights for the
		Gauss-Kronrod quadrature in multigauss.F
                this file is part of FormCalc
                last modified 10 Jan 03 th
*)


<< ../FormCalc.m

Prec = 50

S[n_, i_, j_] := Integrate[
  LegendreP[2 i - 1 - Mod[n, 2], x] LegendreP[n, x] LegendreP[2 j - 1, x],
  {x, -1, 1} ]

KronrodP[n_, x_] :=
Block[ {r, a, coeff},
  r = Floor[(n + 3)/2];
  coeff = Solve[
    Flatten[{
      a[1] == 1,
      Table[Sum[a[i] S[n, i, j], {i, r - j, r}] == 0, {j, r - 1}]
    }],
    Array[a, r] ][[1]];
  Expand[
    Sum[a[i] LegendreP[2 i - 1 - Mod[n, 2], x], {i, 1, r}] /. coeff ]
]

KronrodNodes[n_] := x /. NSolve[KronrodP[n, x] == 0, x, Prec]

GaussNodes[n_] := x /. NSolve[LegendreP[n, x] == 0, x, Prec]


InterpPoly[nodes_, i_, x_] :=
Block[ {xi = nodes[[i]]},
  Expand[ Times@@ ((x - #)/(xi - #)&)/@ Drop[nodes, {i}] ]
]

Weights[nodes_] := Array[
  Integrate[InterpPoly[nodes, #, x], {x, -1, 1}]&,
  Length[nodes] ]


DataSets[n_, wrap_:TableForm] :=
Block[ {g, k, w},
  g = GaussNodes[n];
  k = KronrodNodes[n];
  w = Weights[{g, k}//Flatten];
  SetPrecision[
    { xg -> wrap[g],
      xk -> wrap[k],
      wka -> wrap[Take[w, Length[g]]],
      wkb -> wrap[Drop[w, Length[g]]],
      wg -> wrap[Weights[g]] }, 35 ]
]


TheData[_?OddQ] = {}

TheData[n_] :=
Block[ {data, l},
  data = Last/@ DataSets[n, Identity];
  l = Length[ data[[1]] ]/2;
  Flatten[ {data[[-2, -l - 1]], MapThread[List, Take[#, -l]&/@ data]} ]
]


FortranWrite[ns__, file_] :=
Block[ {hh, datax, str},
  hh = OpenFortran[file];
  data = Flatten[TheData/@ {ns}];
  str = {"     &    ", ToFortran[#], ",\n"}&/@ data;
  WriteString[hh,
    "\tdouble precision krdata(" <> ToString[Length[data]] <> ")\n" <>
    "\tdata krdata /\n" <> Drop[Flatten[str], -1] <> " /\n"];
  Close[hh]
]


(* offset into array krdata made from FortranWrite[4, 6, 8, ..., file],
   given the (lower) # of points:

   index[n_] = 1 + Simplify[
     Sum[5/2 i + 1, {i, 4, n - 2, 2}],
     Element[n/2, Integers]]
*)

