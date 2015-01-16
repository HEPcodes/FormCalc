(*
	gausspoints.m
		calculate the abscissas & weigths for Gauss int
		last modified 30 Aug 99 th
*)


(* write out abscissas and weights for these numbers of sampling points
   (must be even integers): *)

Points = {8, 16, 24, 32}

(* working precision: *)

Prec = 50

(* the sampling points are the zeros of the Legendre polynomials: *)

SamplingPoints[n_] := SamplingPoints[n] =
  Take[ #[[1, 2]]&/@ NSolve[LegendreP[n, x] == 0, x, Prec], n/2 ]

dif2[n_, x_] := D[LegendreP[n, y], y]^2 /. y -> x

Weights[n_] := 2/((1 - #^2) dif2[n, #])&/@ SamplingPoints[n]

ff = OpenWrite["!./r8_" <> Environment["HOSTTYPE"] <>
               " > fortran/gauss.F"];
c = 0

WriteString[ff,
  delim = Prepend[Table[",\n     +    ", {#/2 - 1}], {}];
  If[++c === 1, "#if", "#elif"], " GAUSSPOINTS == ", #,
  "\n\tdata gauss_x /\n     +    ",
  Transpose[{delim, SamplingPoints[#]}] /. List -> Sequence,
  " /\n\tdata gauss_w /\n     +    ",
  Transpose[{delim, Weights[#]}] /. List -> Sequence,
  " /\n" ]&/@ Points;
WriteString[ff, "#endif\n\n"]

Close[ff];

