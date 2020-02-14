(*
	GaussData.m
		calculate the abscissas and weights for the
		Gaussian quadrature in gauss.F
		this file is part of FormCalc
		last modified 24 Mar 19 th
*)


<< FormCalc`

GaussData[n_] :=
Block[ {nodes, weight, workingprec = 50, finalprec = 35},
  nodes = x /. NSolve[LegendreP[n, x] == 0, x, workingprec];
  weight[x_] = 2/((1 - x^2) D[LegendreP[n, x], x]^2);
  SetPrecision[{#, weight[#]}&/@ Take[-nodes, n/2], finalprec]
]

gaussdata = Table[GaussData[n]/2, {n, 8, 32, 8}]//Flatten

delim := (delim = ","; "\tdata gaussdata /")

hh = OpenFortran["gaussdata.F"];

WriteString[hh,
  "\tRealType gaussdata(" <> ToFortran[Length[gaussdata]] <> ")\n" <>
  ({delim, "\n     &    ", ToFortran[#]}&)/@ gaussdata <>
  " /\n"
]

Close[hh]


(* offset into array gaussdata for a given # of points:
   index[n_] = 1 + Sum[i, {i, 8, n - 8, 8}]
*)

