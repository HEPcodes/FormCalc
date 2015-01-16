<< ~/FeynArts/FeynArts.m

CKM = IndexDelta

SetOptions[ InsertFields,
  Model -> "SMc",
  InsertionLevel -> {Generic, Particles},
  Restrictions -> {NoQuarkMixing, NoLightFHCoupling} ]

inss := InsertFields[ tops,
  {F[2, {1}], -F[2, {1}]} -> {F[3, {3}], -F[3, {3}]} ]

amps[ins_] := ToFA1Conventions[ CreateFeynAmp[ins] ]

(* tree graphs *)

tops = CreateTopologies[ 0, 2 -> 2 ];
amps[inss] >> fa/born.amp

(* self-energies *)

tops = CreateTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, Triangles, AllBoxes} ];
amps[inss] >> fa/self.amp

(* vertex corrections *)

tops = CreateTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, SelfEnergies, AllBoxes} ];
	(* we're not interested in the QED corrections, and since
	   they're IR divergent anyway, let's just discard them: *)
ins = DiagramSelect[inss, FreeQ[#, Field[6|8] -> V[1]]&];
amps[ins] >> fa/vert.amp

(* boxes *)

tops = CreateTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, SelfEnergies, Triangles} ];
ins = DiagramSelect[inss, FreeQ[#, Field[6|7] -> V[1]]&];
amps[ins] >> fa/box.amp

(* counter terms *)

tops = CreateCTTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {TadpoleCTs, WFCorrectionCTs} ];
amps[inss] >> fa/counter.amp

