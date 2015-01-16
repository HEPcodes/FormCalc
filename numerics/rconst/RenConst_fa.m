(*
	RenConst_fa.m
		generates the amplitudes for the self-energies
		needed for calculating the renormalization constants
		this file is part of FormCalc
		last modified 3 Mar 00 th
*)

<< Utilities`FilterOptions`

<< ~/FeynArts/FeynArts.m

SetOptions[Paint, ColumnsXRows -> 4, Destination -> File]

SetOptions[InsertFields, Restrictions -> NoQuarkMixing]

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
MakeAmp["self.hh", S[1] -> S[1]];
MakeAmp["self.cc", S[2] -> S[2]];
MakeAmp["self.pp", S[3] -> S[3]]

(* this is for calculating with ME = 0 to avoid mass singularities: *)
SetOptions[InsertFields, ExcludeParticles -> V[1]]

MakeAmp["self.nn", F[1, {g1}] -> F[1, {g1}], Truncated -> True];
MakeAmp["self.ee", F[2, {g1}] -> F[2, {g1}], Truncated -> True];
MakeAmp["self.uu", F[3, {g1}] -> F[3, {g1}], Truncated -> True];
MakeAmp["self.dd", F[4, {g1}] -> F[4, {g1}], Truncated -> True]

SetOptions[InsertFields, ExcludeParticles -> {}]

tops = topsTad =
  CreateTopologies[1, 1 -> 0, ExcludeTopologies -> SelfEnergies]

MakeAmp["tad.h", S[1] -> {}]


SetOptions[InsertFields,
  Model -> SMbgf, GenericModel -> Lorentzbgf]

tops = topsSE

MakeAmp["self.aa.bgf", V[10] -> V[10]];
MakeAmp["self.az.bgf", V[10] -> V[20]];
MakeAmp["self.zz.bgf", V[20] -> V[20]];
MakeAmp["self.ww.bgf", V[30] -> V[30]];
MakeAmp["self.hh.bgf", S[10] -> S[10]];
MakeAmp["self.cc.bgf", S[20] -> S[20]];
MakeAmp["self.pp.bgf", S[30] -> S[30]]

tops = topsTad

MakeAmp["tad.h.bgf", S[10] -> {}]

