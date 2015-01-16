(*
	AAAA_smbgf.m  
		generates the Fortran code for
		gamma gamma -> gamma gamma in the electroweak SM
		using the background-field method
		this file is part of FormCalc
		last modified 28 Aug 01 th

Reference: M. Boehm, R. Schuster, Z. Phys. C63 (1994) 219.

*)


<< FeynArts`
<< ../../FormCalc.m

time1 = SessionTime[]

CKM = IndexDelta


SetOptions[InsertFields,
  Model -> "SMbgf", GenericModel -> "Lorentzbgf"]

AAAA = {V[10], V[10]} -> {V[10], V[10]}

inss := ins = InsertFields[tops, AAAA]


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

WriteSquaredME[born, {self, vert, box}, abbr, "fortran_smbgf",
  Drivers -> "drivers_sm"]

WriteRenConst[counter, "fortran_smbgf"]

Print["time used: ", SessionTime[] - time1];

