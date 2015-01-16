(*
	WWWW_sm.m  
		generates the Fortran code for
		W^+ W^- -> W^+ W^- in the electroweak Standard Model
		this file is part of FormCalc
		last modified 13 Feb 03 th

Reference: A. Denner, T. Hahn,
           Nucl. Phys. B525 (1998) 27 [hep-ph/9711302].
*)


<< FeynArts`

<< FormCalc`


time1 = SessionTime[]

CKM = IndexDelta


SetOptions[InsertFields, Model -> "SM"]

process = {-V[3], V[3]} -> {-V[3], V[3]}


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

tops = CreateTopologies[0, 2 -> 2];
ins = InsertFields[tops, process];
DoPaint[ins, "born"];
born = CalcFeynAmp[CreateFeynAmp[ins]]


Print["Counter terms"]

tops = CreateCTTopologies[1, 2 -> 2,
  ExcludeTopologies -> {TadpoleCTs, WFCorrectionCTs}];
ins = InsertFields[tops, process];
DoPaint[ins, "counter"];
counter = CreateFeynAmp[ins]


Print["Self energies"]

tops = CreateTopologies[1, 2 -> 2, SelfEnergiesOnly];
ins = InsertFields[tops, process];
DoPaint[ins, "self"];
self = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 2 &]]


Print["Vertices"]

tops = CreateTopologies[1, 2 -> 2, TrianglesOnly];
ins = InsertFields[tops, process];
DoPaint[ins, "vert"];
vert = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 1 &]]


Print["Boxes"]

tops = CreateTopologies[1, 2 -> 2, BoxesOnly];
ins = InsertFields[tops, process];
DoPaint[ins, "box"];
box = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 0 &]]


abbr = OptimizeAbbr[Abbr[]]


WriteSquaredME[born, {self, vert, box}, abbr, "fortran_sm",
  Drivers -> "drivers_sm"]

WriteRenConst[{self, vert, box}, "fortran_sm"]


Print["time used: ", SessionTime[] - time1]

