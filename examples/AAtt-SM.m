(*
	AAtt-SM.m
		generates the Fortran code for
		gamma gamma -> t-bar t in the electroweak SM
		this file is part of FormCalc
		last modified 11 Jun 03 th

Reference: A. Denner, S. Dittmaier, and M. Strobel,
           Phys. Rev. D53 (1996) 44 [hep-ph/9507372].

*)


Needs["FeynArts`"]

Needs["FormCalc`"]


time1 = SessionTime[]

CKM = IndexDelta

Small[ME] = Small[ME2] = 0


process = {V[1], V[1]} -> {-F[3, {3}], F[3, {3}]}

name = "AAtt-SM"

SetOptions[InsertFields, Model -> "SMc", Restrictions -> NoLightFHCoupling]


SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]

(* take the comments out if you want the diagrams painted
DoPaint[diags_, file_] := Paint[diags, DisplayFunction ->
  (Display[ToFileName[MkDir[name <> ".diagrams"], file <> ".ps"], #]&)]
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


col = ColourME[All, born]

abbr = OptimizeAbbr[Abbr[]]


dir = SetupCodeDir[name <> ".fortran", Drivers -> name <> ".drivers"]

WriteSquaredME[born, {self, vert, box}, col, abbr, dir]

WriteRenConst[counter, dir]


Print["time used: ", SessionTime[] - time1]

