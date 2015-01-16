(* calculate the abscissas & weigths for gauss int *)

absL[n_] := absL[n] =
  #[[1, 2]]&/@ NSolve[LegendreP[n, x] == 0, x, 50]

weiL[n_] :=
Block[ {dif2},
  dif2[x_] = D[LegendreP[n, x], x]^2;
  2/((1 - #^2) dif2[#])&/@ absL[n]
]

ff = OpenWrite["gauss.F"];
Do[
  WriteString[ff,
    If[i == 8, "#if", "#elif"], " GAUSSPOINTS==", i,
    "\n\tdata gauss_x /\n     +    ", First[absL[i]]];
  WriteString[ff, ",\n     +    ", #]&/@ Take[absL[i], {2, i/2}];
  WriteString[ff, " /\n\tdata gauss_w /\n     +    ", First[weiL[i]]];
  WriteString[ff, ",\n     +    ", #]&/@ Take[weiL[i], {2, i/2}];
  WriteString[ff," /\n"],
{i, 8, 32, 8}];
WriteString[ff, "#endif\n\n"];
Close[ff];

