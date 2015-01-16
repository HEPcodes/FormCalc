<< ~/develop/HighEnergyPhysics/FeynArts.m

SetOptions[ InsertFields,
  Model -> "SM",
  InsertionLevel -> {Particles},
  Restrictions -> {NoGeneration2, NoGeneration3} ]

amps :=
  ToFA1Conventions[
    CreateFeynAmp[
      InsertFields[ tops, {V[2], V[2]} -> {V[2], V[2]} ] ] ]

(* tree graphs *)

tops = CreateTopologies[ 0, 2 -> 2 ];
amps >> fa_output/zzzz.born.amp

(* self-energies *)

tops = CreateTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, Triangles, AllBoxes} ];
amps >> fa_output/zzzz.self.amp

(* vertex corrections *)

tops = CreateTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, SelfEnergies, AllBoxes} ];
amps >> fa_output/zzzz.vert.amp

(* boxes *)

tops = CreateTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, SelfEnergies, Triangles} ];
amps >> fa_output/zzzz.box.amp

(* counter terms *)

tops = CreateCTTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {TadpoleCTs, WFCorrectionCTs} ];
amps >> fa_output/zzzz.counter.amp

