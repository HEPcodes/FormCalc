<< ~/develop/HighEnergyPhysics/FeynArts.m

CKM[ i_Integer, i_Integer ] = 1;
CKM[ i_Integer, j_Integer ] = 0

SetOptions[ InsertFields,
  Restrictions -> {NoGeneration2, NoGeneration3},
  InsertionLevel -> {Particles} ]

MakeAmp[ fromto_ ] :=
Block[ {ins, amp},
  ins = InsertFields[tops, fromto];
  amp = CreateFeynAmp[ins];
  ToFA1Conventions[amp]
]

tops = CreateTopologies[1, 1 -> 1, ExcludeTopologies -> Tadpoles];
MakeAmp[V[1] -> V[1]] >> rc/self.aa;
MakeAmp[V[1] -> V[2]] >> rc/self.az;
MakeAmp[V[2] -> V[2]] >> rc/self.zz;
MakeAmp[V[3] -> V[3]] >> rc/self.ww;
MakeAmp[S[1] -> S[1]] >> rc/self.hh;
MakeAmp[S[2] -> S[2]] >> rc/self.cc;
MakeAmp[S[3] -> S[3]] >> rc/self.pp

