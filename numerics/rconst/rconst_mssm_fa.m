(*
	rconst_mssm_fa.m
		generates the amplitudes for the self-energies
		needed for calculating the renormalization constants
		last modified 14 Jun 99 th
*)

<< Utilities`FilterOptions`

<< ~/develop/HighEnergyPhysics/FeynArts.m

SetOptions[Paint, ColumnsXRows -> 4, Destination -> File]

SetOptions[InsertFields, Model -> MSSM]

CKM = IndexDelta

MakeAmp[name_, fromto_, opt___Rule] :=
Block[ {ins, seq = Sequence[]},
  ins = InsertFields[tops, fromto, FilterOptions[InsertFields, opt]];
(*  Paint[ins]; *)
  Put[
    ToFA1Conventions[
      CreateFeynAmp[ins, FilterOptions[CreateFeynAmp, opt]] ],
    "fa/" <> name ]
]


tops = topsSE =
  CreateTopologies[1, 1 -> 1, ExcludeTopologies -> Tadpoles]

MakeAmp["self.aa", V[1] -> V[1]];
MakeAmp["self.az", V[1] -> V[2]];
MakeAmp["self.zz", V[2] -> V[2]];
MakeAmp["self.ww", V[3] -> V[3]];
MakeAmp["self.pp", S[5] -> S[5], ExcludeParticles -> V[1]];

(* this is for calculating with ME = 0 to avoid mass singularities: *)
SetOptions[InsertFields, ExcludeParticles -> V[1]]

MakeAmp["self.nn", F[1, {g1}] -> F[1, {g2}], Truncated -> True];
MakeAmp["self.ee", F[2, {g1}] -> F[2, {g2}], Truncated -> True];
MakeAmp["self.uu", F[3, {g1, c1}] -> F[3, {g2, c1}], Truncated -> True];
MakeAmp["self.dd", F[4, {g1, c1}] -> F[4, {g2, c1}], Truncated -> True]

