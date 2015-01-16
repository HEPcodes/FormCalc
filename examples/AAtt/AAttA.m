(*
	AAttA.m
		generates the Fortran code for
		\gamma \gamma -> \bar t t \gamma in the electroweak SM
		this file is part of FormCalc
		last modified 30 Aug 01 th

Reference: W. Walter, Diploma thesis, Wuerzburg 1997.

*)


<< FeynArts`
<< ../../FormCalc.m


time1 = SessionTime[]

CKM = IndexDelta

Small[ME] = Small[ME2] = 0


AAtt = {V[1], V[1]} -> {-F[3, {3}], F[3, {3}], V[1]}

SetOptions[InsertFields,
  Model -> "SMc", Restrictions -> NoLightFHCoupling]

SetOptions[CalcFeynAmp, Dimension -> 4]

inss := ins = InsertFields[tops, AAtt]


SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]

(* take the comments out if you want the diagrams painted
DoPaint[diags_, file_] := (
  If[ FileType["diagrams"] =!= Directory,
    CreateDirectory["diagrams"] ];
  Paint[diags,
    DisplayFunction -> (Display["diagrams/" <> file <> ".ps", #]&)]
)
*)


tops = CreateTopologies[0, 2 -> 3];
DoPaint[inss, "born"];
born = CalcFeynAmp[CreateFeynAmp[ins]]

num = Simplify

Hel[_] = 0;
hel = HelicityME[All, born];
col = ColourME[All, born]

abbr = Abbr[]

(*
abbr = OptimizeAbbr[Abbr[]]
*)

WriteSquaredME[born, {}, hel, col, abbr, "fortran_brems",
  Drivers -> "drivers_brems"]

WriteRenConst[0, "fortran_brems"]

Print["time used: ", SessionTime[] - time1]

