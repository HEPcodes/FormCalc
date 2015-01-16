(*
	eettA.m
		generates the Fortran code for
		e^+ e^- -> \bar t t gamma in the electroweak SM
		this file is part of FormCalc
		last modified 12 Feb 03 th

Reference: W. Beenakker, S.C. van der Marck, and W. Hollik,
           Nucl. Phys. B365 (1991) 24.

*)


<< FeynArts`

<< FormCalc`


time1 = SessionTime[]

CKM = IndexDelta

Small[ME] = Small[ME2] = 0


process = {-F[2, {1}], F[2, {1}]} -> {-F[3, {3}], F[3, {3}], V[1]}

SetOptions[InsertFields, Model -> "SMc", Restrictions -> NoLightFHCoupling]

inss := ins = InsertFields[tops, process]


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
ins = InsertFields[tops, process];
DoPaint[ins, "born"];
born = CalcFeynAmp[CreateFeynAmp[ins]]


Hel[_] = 0;
hel = HelicityME[All, born]

col = ColourME[All, born]


WriteSquaredME[born, {}, hel, col, Abbr[], "fortran_brems",
  Drivers -> "drivers_brems"]

WriteRenConst[{}, "fortran_brems"]


Print["time used: ", SessionTime[] - time1]

