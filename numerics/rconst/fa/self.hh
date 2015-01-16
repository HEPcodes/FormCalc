FeynAmpList[Model -> "SM", GenericModel -> "Lorentz", 
   InsertionLevel -> Classes, ExcludeParticles -> {}, LastSelections -> {}, 
   ExcludeFieldPoints -> {}, Restrictions -> 
    ExcludeFieldPoints -> 
     {FieldPoint[0][-F[4, {1}], F[3, {2}], S[3]], 
      FieldPoint[0][-F[4, {1}], F[3, {2}], V[3]], 
      FieldPoint[0][-F[4, {1}], F[3, {3}], S[3]], 
      FieldPoint[0][-F[4, {1}], F[3, {3}], V[3]], 
      FieldPoint[0][-F[4, {2}], F[3, {1}], S[3]], 
      FieldPoint[0][-F[4, {2}], F[3, {1}], V[3]], 
      FieldPoint[0][-F[4, {2}], F[3, {3}], S[3]], 
      FieldPoint[0][-F[4, {2}], F[3, {3}], V[3]], 
      FieldPoint[0][-F[4, {3}], F[3, {1}], S[3]], 
      FieldPoint[0][-F[4, {3}], F[3, {1}], V[3]], 
      FieldPoint[0][-F[4, {3}], F[3, {2}], S[3]], 
      FieldPoint[0][-F[4, {3}], F[3, {2}], V[3]]}, VertexFunctions -> False, 
   ProcessName -> "S1S1", Process -> {{S[1], p1, MH}} -> {{S[1], k1, MH}}][
  FeynAmp[GraphName["S1S1", T1, I1], q1, 
   (-I/32*RelativeCF*DiracTrace[DiracSlash[-q1] + 
        Mass[F[Index[Generic, 3]]], 
       ChiralityProjector[-1]*
         G[1][0][F[Index[Generic, 3]], F[Index[Generic, 4]], S[1]][
          ChiralityProjector[-1]] + 
        ChiralityProjector[1]*
         G[1][0][F[Index[Generic, 3]], F[Index[Generic, 4]], S[1]][
          ChiralityProjector[1]], 
       DiracSlash[k1 - q1] + Mass[F[Index[Generic, 4]]], 
       ChiralityProjector[-1]*
         G[1][0][-F[Index[Generic, 4]], -F[Index[Generic, 3]], S[1]][
          ChiralityProjector[-1]] + 
        ChiralityProjector[1]*
         G[1][0][-F[Index[Generic, 4]], -F[Index[Generic, 3]], S[1]][
          ChiralityProjector[1]]]*
      FeynAmpDenominator[PropagatorDenominator[q1, 
        Mass[F[Index[Generic, 3]]]], 
       PropagatorDenominator[-k1 + q1, Mass[F[Index[Generic, 4]]]]])/Pi^4, 
   {Mass[F[Index[Generic, 3]]], Mass[F[Index[Generic, 4]]], 
     G[1][0][F[Index[Generic, 3]], F[Index[Generic, 4]], S[1]][
      ChiralityProjector[-1]], 
     G[1][0][F[Index[Generic, 3]], F[Index[Generic, 4]], S[1]][
      ChiralityProjector[1]], G[1][0][-F[Index[Generic, 4]], 
       -F[Index[Generic, 3]], S[1]][ChiralityProjector[-1]], 
     G[1][0][-F[Index[Generic, 4]], -F[Index[Generic, 3]], S[1]][
      ChiralityProjector[1]], RelativeCF} -> 
    Insertions[Classes][{MLE[Index[Generation, 1]], 
      MLE[Index[Generation, 1]], 
      (-I/2*EL*MLE[Index[Generation, 1]])/(MW*SW), 
      (-I/2*EL*MLE[Index[Generation, 1]])/(MW*SW), 
      (-I/2*EL*MLE[Index[Generation, 1]])/(MW*SW), 
      (-I/2*EL*MLE[Index[Generation, 1]])/(MW*SW), 
      2*SumOver[Index[Generation, 1], 3]}, 
     {MQU[Index[Generation, 1]], MQU[Index[Generation, 1]], 
      (-I/2*EL*MQU[Index[Generation, 1]])/(MW*SW), 
      (-I/2*EL*MQU[Index[Generation, 1]])/(MW*SW), 
      (-I/2*EL*MQU[Index[Generation, 1]])/(MW*SW), 
      (-I/2*EL*MQU[Index[Generation, 1]])/(MW*SW), 
      6*SumOver[Index[Generation, 1], 3]}, 
     {MQD[Index[Generation, 1]], MQD[Index[Generation, 1]], 
      (-I/2*EL*MQD[Index[Generation, 1]])/(MW*SW), 
      (-I/2*EL*MQD[Index[Generation, 1]])/(MW*SW), 
      (-I/2*EL*MQD[Index[Generation, 1]])/(MW*SW), 
      (-I/2*EL*MQD[Index[Generation, 1]])/(MW*SW), 
      6*SumOver[Index[Generation, 1], 3]}]], 
  FeynAmp[GraphName["S1S1", T1, I2], q1, 
   (I/32*RelativeCF*FeynAmpDenominator[PropagatorDenominator[q1, 
        Mass[S[Index[Generic, 3]]]], 
       PropagatorDenominator[-k1 + q1, Mass[S[Index[Generic, 4]]]]]*
      G[1][0][S[1], -S[Index[Generic, 3]], -S[Index[Generic, 4]]][1]*
      G[1][0][S[1], S[Index[Generic, 3]], S[Index[Generic, 4]]][1])/Pi^4, 
   {Mass[S[Index[Generic, 3]]], Mass[S[Index[Generic, 4]]], 
     G[1][0][S[1], -S[Index[Generic, 3]], -S[Index[Generic, 4]]][1], 
     G[1][0][S[1], S[Index[Generic, 3]], S[Index[Generic, 4]]][1], 
     RelativeCF} -> 
    Insertions[Classes][{MH, MH, ((-3*I)/2*EL*MH^2)/(MW*SW), 
      ((-3*I)/2*EL*MH^2)/(MW*SW), 1}, 
     {MZ, MZ, (-I/2*EL*MH^2)/(MW*SW), (-I/2*EL*MH^2)/(MW*SW), 1}, 
     {MW, MW, (-I/2*EL*MH^2)/(MW*SW), (-I/2*EL*MH^2)/(MW*SW), 2}]], 
  FeynAmp[GraphName["S1S1", T1, I3], q1, 
   (-I/32*RelativeCF*FeynAmpDenominator[PropagatorDenominator[q1, 
        Mass[U[Index[Generic, 3]]]], 
       PropagatorDenominator[-k1 + q1, Mass[U[Index[Generic, 4]]]]]*
      G[1][0][S[1], -U[Index[Generic, 3]], -U[Index[Generic, 4]]][1]*
      G[1][0][S[1], U[Index[Generic, 3]], U[Index[Generic, 4]]][1])/Pi^4, 
   {Mass[U[Index[Generic, 3]]], Mass[U[Index[Generic, 4]]], 
     G[1][0][S[1], -U[Index[Generic, 3]], -U[Index[Generic, 4]]][1], 
     G[1][0][S[1], U[Index[Generic, 3]], U[Index[Generic, 4]]][1], 
     RelativeCF} -> 
    Insertions[Classes][{MZ, MZ, (-I/2*EL*MW)/(CW^2*SW), 
      (-I/2*EL*MW)/(CW^2*SW), 2}, 
     {MW, MW, (-I/2*EL*MW)/SW, (-I/2*EL*MW)/SW, 2}, 
     {MW, MW, (-I/2*EL*MW)/SW, (-I/2*EL*MW)/SW, 2}]], 
  FeynAmp[GraphName["S1S1", T1, I4], q1, 
   (I/32*RelativeCF*FeynAmpDenominator[PropagatorDenominator[q1, 
        Mass[V[Index[Generic, 3]]]], 
       PropagatorDenominator[-k1 + q1, Mass[V[Index[Generic, 4]]]]]*
      MetricTensor[li1, li2]*MetricTensor[li1, li3]*MetricTensor[li2, li4]*
      MetricTensor[li3, li4]*G[1][0][S[1], -V[Index[Generic, 3]], 
        -V[Index[Generic, 4]]][MetricTensor[KI1[2], KI1[3]]]*
      G[1][0][S[1], V[Index[Generic, 3]], V[Index[Generic, 4]]][
       MetricTensor[KI1[2], KI1[3]]])/Pi^4, 
   {Mass[V[Index[Generic, 3]]], Mass[V[Index[Generic, 4]]], 
     G[1][0][S[1], -V[Index[Generic, 3]], -V[Index[Generic, 4]]][
      MetricTensor[KI1[2], KI1[3]]], 
     G[1][0][S[1], V[Index[Generic, 3]], V[Index[Generic, 4]]][
      MetricTensor[KI1[2], KI1[3]]], RelativeCF} -> 
    Insertions[Classes][{MZ, MZ, (I*EL*MW)/(CW^2*SW), (I*EL*MW)/(CW^2*SW), 
      1}, {MW, MW, (I*EL*MW)/SW, (I*EL*MW)/SW, 2}]], 
  FeynAmp[GraphName["S1S1", T1, I5], q1, 
   (-I/16*RelativeCF*FeynAmpDenominator[PropagatorDenominator[q1, 
        Mass[S[Index[Generic, 3]]]], 
       PropagatorDenominator[-k1 + q1, Mass[V[Index[Generic, 4]]]]]*
      FourVector[-k1 - q1, li2]*FourVector[p1 + q1, li1]*
      MetricTensor[li1, li2]*G[-1][0][S[1], -S[Index[Generic, 3]], 
        -V[Index[Generic, 4]]][FourVector[Mom[1] - Mom[2], KI1[3]]]*
      G[-1][0][S[1], S[Index[Generic, 3]], V[Index[Generic, 4]]][
       FourVector[Mom[1] - Mom[2], KI1[3]]])/Pi^4, 
   {Mass[S[Index[Generic, 3]]], Mass[V[Index[Generic, 4]]], 
     G[-1][0][S[1], -S[Index[Generic, 3]], -V[Index[Generic, 4]]][
      FourVector[Mom[1] - Mom[2], KI1[3]]], 
     G[-1][0][S[1], S[Index[Generic, 3]], V[Index[Generic, 4]]][
      FourVector[Mom[1] - Mom[2], KI1[3]]], RelativeCF} -> 
    Insertions[Classes][{MZ, MZ, -EL/(2*CW*SW), -EL/(2*CW*SW), 1}, 
     {MW, MW, (I/2*EL)/SW, (-I/2*EL)/SW, 1}, 
     {MW, MW, (-I/2*EL)/SW, (I/2*EL)/SW, 1}]], 
  FeynAmp[GraphName["S1S1", T2, I1], q1, 
   (RelativeCF*FeynAmpDenominator[PropagatorDenominator[q1, 
        Mass[S[Index[Generic, 3]]]]]*
      G[1][0][S[1], S[1], -S[Index[Generic, 3]], S[Index[Generic, 3]]][1])/
    (32*Pi^4), {Mass[S[Index[Generic, 3]]], 
     G[1][0][S[1], S[1], -S[Index[Generic, 3]], S[Index[Generic, 3]]][1], 
     RelativeCF} -> 
    Insertions[Classes][{MH, ((-3*I)/4*EL^2*MH^2)/(MW^2*SW^2), 1}, 
     {MZ, (-I/4*EL^2*MH^2)/(MW^2*SW^2), 1}, 
     {MW, (-I/4*EL^2*MH^2)/(MW^2*SW^2), 2}]], 
  FeynAmp[GraphName["S1S1", T2, I2], q1, 
   -(RelativeCF*FeynAmpDenominator[PropagatorDenominator[q1, 
         Mass[V[Index[Generic, 3]]]]]*MetricTensor[li1, li2]^2*
       G[1][0][S[1], S[1], -V[Index[Generic, 3]], V[Index[Generic, 3]]][
        MetricTensor[KI1[3], KI1[4]]])/(32*Pi^4), 
   {Mass[V[Index[Generic, 3]]], 
     G[1][0][S[1], S[1], -V[Index[Generic, 3]], V[Index[Generic, 3]]][
      MetricTensor[KI1[3], KI1[4]]], RelativeCF} -> 
    Insertions[Classes][{MZ, (I/2*EL^2)/(CW^2*SW^2), 1}, 
     {MW, (I/2*EL^2)/SW^2, 2}]]]
