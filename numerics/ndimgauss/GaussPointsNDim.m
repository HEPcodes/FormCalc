(*
	GaussPointsNDim.m
		calculate the abscissas & weigths for n-dim Gauss int
		last modified 26 Apr 99 th
*)


(* These are the number of sampling points in each dimension, i.e.
   Length[Points] is the number of dimensions.
   The first element in Points is for the innermost integration.
   Points is usually supplied by the makefile. 
   NB: Points is NOT the same thing as in gausspoints.m (the
   one-dimensional version). *)

If[ !ValueQ[Points], Points = {6, 6, 6, 6, 8, 8, 8} ];

(* working precision: *)

Prec = 50

(* the sampling points are the zeros of the Legendre polynomials: *)

SamplingPoints[n_] := SamplingPoints[n] =
  #[[1, 2]]&/@ NSolve[LegendreP[n, x] == 0, x, Prec]

dif2[n_, x_] := D[LegendreP[n, y], y]^2 /. y -> x

Weights[n_] := 2/((1 - #^2) dif2[n, #])&/@ SamplingPoints[n]

Points /= 2;
maxlen = Max[Points]

ZeroFill[func_] := 
  Flatten[Transpose[
    Join[Take[func[2 #], #], Table[0., {maxlen - #}]]&/@ Points ]]

ThePoints = ZeroFill[SamplingPoints];
TheWeights = ZeroFill[Weights]/2	(* note that the weights are 
					   already taken /= 2 *)

r8path = If[ FileType[$Input] === File, $Input,
	(* or, if loaded from a directory in $Path: *)
  Block[ {full},
    Scan[
      If[ FileType[full = # <> "/" <> $Input] === File, Return[full] ]&,
      $Path ] ]
];
pos = StringPosition[r8path, "/"];
If[ Length[pos] === 0, r8path = "",
  r8path = SetDirectory[StringTake[ r8path, pos[[-1, -1]] ]] <> "/";  
  ResetDirectory[] ]

ff = OpenWrite[ "!" <> r8path <>
  "../r8_" <> Environment["HOSTTYPE"] <> " > GaussNDim.f",
  FormatType -> FortranForm, PageWidth -> 67]

WriteString[ff, "\tinteger points(", Length[Points], ")\n"];
WriteString[ff,
  "\tdouble precision gauss_x(", Length[Points], ", ", maxlen,
  "), gauss_w(", Length[Points], ", ", maxlen, ")\n"];
WriteString[ff, "\tsave points, gauss_x, gauss_w\n"];
WriteString[ff, "\tdata points /", First[Points]];
WriteString[ff, ", ", #]&/@ Rest[Points];
WriteString[ff, "/\n\tdata gauss_x /\n     +    ", First[ThePoints]];
WriteString[ff,",\n     +    ",#]&/@ Rest[ThePoints];
WriteString[ff, " /\n\tdata gauss_w /\n     +    ", First[TheWeights]];
WriteString[ff,",\n     +    ",#]&/@ Rest[TheWeights];
WriteString[ff," /\n"]

Close[ff]
