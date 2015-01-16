(*
	AAtt_sm.m
		generates the Fortran code for
		\gamma \gamma -> \bar t t in the electroweak SM
		this file is part of FormCalc
		last modified 30 Aug 01 th

Reference: A. Denner, S. Dittmaier, and M. Strobel,
           Phys. Rev. D53 (1996) 44 (hep-ph/9507372).

*)


<< FeynArts`
<< ../../FormCalc.m


time1 = SessionTime[]

CKM = IndexDelta

Small[ME] = Small[ME2] = 0


AAtt = {V[1], V[1]} -> {-F[3, {3}], F[3, {3}]}

SetOptions[InsertFields,
  Model -> "SMc", Restrictions -> NoLightFHCoupling]

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

Hel[_] = 0;
hel = HelicityME[All, born];
col = ColourME[All, born]

abbr = OptimizeAbbr[Abbr[]]

WriteSquaredME[born, {self, vert, box}, hel, col, abbr, "fortran_sm",
  Drivers -> "drivers_sm"]


WriteRenConst[counter, "fortran_sm"]

Print["time used: ", SessionTime[] - time1]

