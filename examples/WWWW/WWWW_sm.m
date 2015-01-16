(*
	WWWW_sm.m  
		generates the Fortran code for
		W^+ W^- -> W^+ W^- in the electroweak Standard Model
		this file is part of FormCalc
		last modified 20 Jul 01 th

Reference: A. Denner, T. Hahn,
           Nucl. Phys. B525 (1998) 27 (hep-ph/9711302).
*)


<< FeynArts`
<< ../../FormCalc.m

time1 = SessionTime[]

CKM = IndexDelta


SetOptions[InsertFields, Model -> "SM"]

WWWW = {-V[3], V[3]} -> {-V[3], V[3]}

inss := ins = InsertFields[tops, WWWW]


SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]

(* take the comments out if you want the diagrams painted
DoPaint[diags_, file_] := (
  If[ FileType["diagrams"] =!= Directory,
    CreateDirectory["diagrams"] ];
  Paint[diags,
    DisplayFunction -> (Display["diagrams/" <> file <> ".ps", #]&)]
)
*)


tops = CreateTopologies[0, 2 -> 2];
DoPaint[inss, "born"];
born = CalcFeynAmp[CreateFeynAmp[ins]]

tops = CreateCTTopologies[1, 2 -> 2,
  ExcludeTopologies -> {TadpoleCTs, WFCorrectionCTs}];
DoPaint[inss, "counter"];
counter = CreateFeynAmp[ins]

tops = CreateTopologies[1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, Triangles, AllBoxes}];
DoPaint[inss, "self"];
self = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 2 &]]

tops = CreateTopologies[1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, SelfEnergies, AllBoxes}];
DoPaint[inss, "vert"];
vert = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 1 &]]

tops = CreateTopologies[1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, SelfEnergies, Triangles}];
DoPaint[inss, "box"];
box = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 0 &]]

abbr = OptimizeAbbr[Abbr[]]

WriteSquaredME[born, {self, vert, box}, abbr, "fortran_sm",
  Drivers -> "drivers_sm"]

WriteRenConst[counter, "fortran_sm"]

Print["time used: ", SessionTime[] - time1];

