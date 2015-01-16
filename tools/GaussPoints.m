(*
	GaussPoints.m
		calculate the abscissas and weights for the
		Gaussian quadrature in 2to2.F
		this file is part of FormCalc
		last modified 10 Jan 03 th
*)


<< ../FormCalc.m

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

hh = OpenFortran[ToFileName[$DriversDir, "gausspoints.F"]];

theif := (theif = "#elif"; "#if")


WriteString[hh, "\
\tdouble precision gauss_x(GAUSSPOINTS/2)\n\
\tdouble precision gauss_w(GAUSSPOINTS/2)\n\n"]

( delim = Prepend[Table[",\n     &    ", {#/2 - 1}], {}];
  WriteString[hh,
    theif <> " GAUSSPOINTS == " <> ToString[#] <>
    "\n\tdata gauss_x /\n     &    " <>
    Transpose[{delim,
      ToString/@ SetPrecision[SamplingPoints[#], 35]}] <>
    " /\n\tdata gauss_w /\n     &    " <>
    Transpose[{delim,
      ToString/@ SetPrecision[Weights[#], 35]}] <> " /\n" ]
)&/@ Points;

WriteString[hh, "#endif\n"]

Close[hh]

