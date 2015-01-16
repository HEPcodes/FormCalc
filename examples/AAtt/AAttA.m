(*
	AAttA.m
		generates the Fortran code for
		\gamma \gamma -> \bar t t \gamma in the electroweak SM
		this file is part of FormCalc
		last modified 13 Feb 03 th

Reference: W. Walter, Diploma thesis, Wuerzburg 1997.

*)


<< FeynArts`

<< FormCalc`


time1 = SessionTime[]

CKM = IndexDelta

Small[ME] = Small[ME2] = 0


process = {V[1], V[1]} -> {-F[3, {3}], F[3, {3}], V[1]}

SetOptions[InsertFields, Model -> "SMc", Restrictions -> NoLightFHCoupling]

SetOptions[CalcFeynAmp, Dimension -> 4]


SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]

(* take the comments out if you want the diagrams painted
DoPaint[diags_, file_] := (
  If[ FileType["diagrams"] =!= Directory,
    CreateDirectory["diagrams"] ];
  Paint[diags,
    DisplayFunction -> (Display["diagrams/" <> file <> ".ps", #]&)]
)
*)


Print["Born"]

tops = CreateTopologies[0, 2 -> 3];
ins = InsertFields[tops, process];
DoPaint[ins, "born"];
born = CalcFeynAmp[CreateFeynAmp[ins]]


Hel[_] = 0;
hel = HelicityME[All, born]

col = ColourME[All, born]

abbr = Abbr[]

(*
abbr = OptimizeAbbr[Abbr[]]
*)

WriteSquaredME[born, {}, hel, col, abbr, "fortran_brems",
  Drivers -> "drivers_brems"]

WriteRenConst[{}, "fortran_brems"]


Print["time used: ", SessionTime[] - time1]

