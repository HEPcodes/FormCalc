(*
	gausspoints.m
		calculate the abscissas & weigths for Gaussian int.
		this file is part of FormCalc
		last modified 29 Feb 00 th
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

hh = OpenWrite["!./r8_" <> Environment["HOSTTYPE"] <>
               " > drivers/gauss.F"];

theif := (theif = "#elif"; "#if")

( delim = Prepend[Table[",\n     +    ", {#/2 - 1}], {}];
  WriteString[hh,
    theif <> " GAUSSPOINTS == " <> ToString[#] <>
    "\n\tdata gauss_x /\n     +    " <>
    Transpose[{delim, ToString/@ SamplingPoints[#]}] <>
    " /\n\tdata gauss_w /\n     +    " <>
    Transpose[{delim, ToString/@ Weights[#]}] <> " /\n" ]
)&/@ Points;

WriteString[hh, "#endif\n\n"]

Close[hh];

