(*
	eeWW_sm.m
		generates the Fortran code for
		e^+ e^- -> W^+ W^- in the electroweak Standard Model
		this file is part of FormCalc
		last modified 21 Jun 01 th

Reference: W. Beenakker, A. Denner,
           Int. J. Mod. Phys. A9 (1994) 4837.
*)


<< FeynArts`
<< ../../FormCalc.m


time1 = SessionTime[]

CKM = IndexDelta

Small[ME] = Small[ME2] = 0


eeWW = {-F[2, {1}], F[2, {1}]} -> {-V[3], V[3]}

SetOptions[InsertFields,
  Model -> "SM", Restrictions -> NoLightFHCoupling]

inss := ins = InsertFields[tops, eeWW]


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
ins = DiagramSelect[inss, FreeQ[#, Field[5|6] -> S]&];
DoPaint[ins, "self"];
self = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 2 &]]

tops = CreateTopologies[1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, SelfEnergies, AllBoxes}];
ins = DiagramSelect[inss, FreeQ[#, Field[5] -> S]&];
DoPaint[ins, "vert"];
vert = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 1 &]]

tops = CreateTopologies[1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, SelfEnergies, Triangles}];
DoPaint[inss, "box"];
box = CalcFeynAmp[
  CreateFeynAmp[ins],
  Select[counter, DiagramType[#] == 0 &]]

Hel[_] = 0;
hel = HelicityME[All, born];

WriteSquaredME[born, {self, vert, box}, hel, Abbr[], "fortran_sm",
  Drivers -> "drivers_sm"]

WriteRenConst[counter, "fortran_sm"]

Print["time used: ", SessionTime[] - time1]

