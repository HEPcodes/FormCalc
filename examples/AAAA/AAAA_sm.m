(*
	AAAA_sm.m  
		generates the Fortran code for
		A A -> A A in the electroweak Standard Model
		this file is part of FormCalc
		last modified 13 Feb 03 th

Reference: M. Boehm, R. Schuster, Z. Phys. C63 (1994) 219.

*)


<< FeynArts`

<< FormCalc`


time1 = SessionTime[]

CKM = IndexDelta


SetOptions[InsertFields, Model -> "SM"]

process = {V[1], V[1]} -> {V[1], V[1]}


SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]

(* take the comments out if you want the diagrams painted
DoPaint[diags_, file_] := (
  If[ FileType["diagrams"] =!= Directory,
    CreateDirectory["diagrams"] ];
  Paint[diags,
    DisplayFunction -> (Display["diagrams/" <> file <> ".ps", #]&)]
)
*)


Print["Boxes"]

tops = CreateTopologies[1, 2 -> 2, BoxesOnly];
ins = InsertFields[tops, process];
DoPaint[ins, "box"];
box = CalcFeynAmp[CreateFeynAmp[ins]]


abbr = OptimizeAbbr[Abbr[]]


WriteSquaredME[{}, box, abbr, "fortran_sm", Drivers -> "drivers_sm"]

WriteRenConst[{}, "fortran_sm"]


Print["time used: ", SessionTime[] - time1]

