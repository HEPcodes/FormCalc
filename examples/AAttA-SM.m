(*
	AAttA-SM.m
		generates the Fortran code for
		gamma gamma -> t-bar t gamma in the electroweak SM
		this file is part of FormCalc
		last modified 11 Jun 03 th

Reference: W. Walter, Diploma thesis, Würzburg 1997.

*)


Needs["FeynArts`"]

Needs["FormCalc`"]


time1 = SessionTime[]

CKM = IndexDelta

Small[ME] = Small[ME2] = 0


process = {V[1], V[1]} -> {-F[3, {3}], F[3, {3}], V[1]}

name = "AAttA-SM"

SetOptions[InsertFields, Model -> "SMc", Restrictions -> NoLightFHCoupling]

SetOptions[CalcFeynAmp, Dimension -> 4]


SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]

(* take the comments out if you want the diagrams painted
DoPaint[diags_, file_] := Paint[diags, DisplayFunction ->
  (Display[ToFileName[MkDir[name <> ".diagrams"], file <> ".ps"], #]&)]
*)


Print["Born"]

tops = CreateTopologies[0, 2 -> 3];
ins = InsertFields[tops, process];
DoPaint[ins, "born"];
born = CalcFeynAmp[CreateFeynAmp[ins]]


col = ColourME[All, born]


abbr = OptimizeAbbr[Abbr[]]


dir = SetupCodeDir[name <> ".fortran", Drivers -> name <> ".drivers"]

WriteSquaredME[born, {}, col, abbr, dir]

WriteRenConst[{}, dir]


Print["time used: ", SessionTime[] - time1]

