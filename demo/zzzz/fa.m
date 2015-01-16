<< ~/develop/HighEnergyPhysics/FeynArts.m

amps :=
  ToFA1Conventions[
    CreateFeynAmp[
      InsertFields[ tops, {V[2], V[2]} -> {V[2], V[2]} ] ] ]

(* tree graphs *)

tops = CreateTopologies[ 0, 2 -> 2 ];
amps >> fa/born.amp

(* self-energies *)

tops = CreateTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, Triangles, AllBoxes} ];
amps >> fa/self.amp

(* vertex corrections *)

tops = CreateTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, SelfEnergies, AllBoxes} ];
amps >> fa/vert.amp

(* boxes *)

tops = CreateTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {Tadpoles, WFCorrections, SelfEnergies, Triangles} ];
amps >> fa/box.amp

(* counter terms *)

tops = CreateCTTopologies[ 1, 2 -> 2,
  ExcludeTopologies -> {TadpoleCTs, WFCorrectionCTs} ];
amps >> fa/counter.amp

