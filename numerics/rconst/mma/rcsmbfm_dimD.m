(* This file contains the renormalization constants for the SM.
   It was generated automatically by rconst.m. DO NOT EDIT. *)

dMf1[2, g1_] = (MLE[g1]*(-(C2*EL^2)/(16*Pi^2) - (EL^2*S2)/(32*Pi^2) - 
   (C2*EL^2*S2)/(64*Pi^2) + (C2*EL^2*SW^2)/(8*Pi^2) + 
   (EL^2*S2*MLE2[g1]*Re[B0[MLE2[g1], MH2, MLE2[g1]]])/(32*MW2*Pi^2) + 
   (C2*EL^2*Re[B0[MLE2[g1], MZ2, MLE2[g1]]])/(4*Pi^2) - 
   (C2*EL^2*SW^2*Re[B0[MLE2[g1], MZ2, MLE2[g1]]])/(2*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[B0[MLE2[g1], MZ2, MLE2[g1]]])/(32*MW2*Pi^2) - 
   (EL^2*S2*Re[B1[MLE2[g1], 0, MW2]])/(16*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[B1[MLE2[g1], 0, MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[B1[MLE2[g1], MLE2[g1], MH2]])/(32*MW2*Pi^2) + 
   (C2*EL^2*Re[B1[MLE2[g1], MLE2[g1], MZ2]])/(8*Pi^2) - 
   (C2*EL^2*S2*Re[B1[MLE2[g1], MLE2[g1], MZ2]])/(32*Pi^2) - 
   (C2*EL^2*SW^2*Re[B1[MLE2[g1], MLE2[g1], MZ2]])/(4*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[B1[MLE2[g1], MLE2[g1], MZ2]])/(32*MW2*Pi^2)))/2

dZfL1[2, g1_, g1_] = -(C2*EL^2)/(16*Pi^2) + (EL^2*S2)/(32*Pi^2) + (C2*EL^2*S2)/(64*Pi^2) + 
 (C2*EL^2*SW^2)/(16*Pi^2) + (EL^2*S2*Re[B1[MLE2[g1], 0, MW2]])/(16*Pi^2) + 
 (EL^2*S2*MLE2[g1]*Re[B1[MLE2[g1], MLE2[g1], MH2]])/(64*MW2*Pi^2) + 
 (-(C2*EL^2)/(8*Pi^2) + (C2*EL^2*S2)/(32*Pi^2) + (C2*EL^2*SW^2)/(8*Pi^2) + 
   (EL^2*S2*MLE2[g1])/(64*MW2*Pi^2))*Re[B1[MLE2[g1], MLE2[g1], MZ2]] - 
 MLE2[g1]*((EL^2*S2*MLE2[g1]*Re[DB0[MLE2[g1], MH2, MLE2[g1]]])/
    (32*MW2*Pi^2) + (C2*EL^2*Re[DB0[MLE2[g1], MZ2, MLE2[g1]]])/(4*Pi^2) - 
   (C2*EL^2*SW^2*Re[DB0[MLE2[g1], MZ2, MLE2[g1]]])/(2*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[DB0[MLE2[g1], MZ2, MLE2[g1]]])/(32*MW2*Pi^2) - 
   (EL^2*S2*Re[DB1[MLE2[g1], 0, MW2]])/(16*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[DB1[MLE2[g1], 0, MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[DB1[MLE2[g1], MLE2[g1], MH2]])/(32*MW2*Pi^2) + 
   (C2*EL^2*Re[DB1[MLE2[g1], MLE2[g1], MZ2]])/(8*Pi^2) - 
   (C2*EL^2*S2*Re[DB1[MLE2[g1], MLE2[g1], MZ2]])/(32*Pi^2) - 
   (C2*EL^2*SW^2*Re[DB1[MLE2[g1], MLE2[g1], MZ2]])/(4*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[DB1[MLE2[g1], MLE2[g1], MZ2]])/(32*MW2*Pi^2))

dZfR1[2, g1_, g1_] = (C2*EL^2*SW^2)/(16*Pi^2) + (EL^2*S2*MLE2[g1]*Re[B1[MLE2[g1], 0, MW2]])/
  (32*MW2*Pi^2) + (EL^2*S2*MLE2[g1]*Re[B1[MLE2[g1], MLE2[g1], MH2]])/
  (64*MW2*Pi^2) + ((C2*EL^2*SW^2)/(8*Pi^2) + (EL^2*S2*MLE2[g1])/
    (64*MW2*Pi^2))*Re[B1[MLE2[g1], MLE2[g1], MZ2]] - 
 MLE2[g1]*((EL^2*S2*MLE2[g1]*Re[DB0[MLE2[g1], MH2, MLE2[g1]]])/
    (32*MW2*Pi^2) + (C2*EL^2*Re[DB0[MLE2[g1], MZ2, MLE2[g1]]])/(4*Pi^2) - 
   (C2*EL^2*SW^2*Re[DB0[MLE2[g1], MZ2, MLE2[g1]]])/(2*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[DB0[MLE2[g1], MZ2, MLE2[g1]]])/(32*MW2*Pi^2) - 
   (EL^2*S2*Re[DB1[MLE2[g1], 0, MW2]])/(16*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[DB1[MLE2[g1], 0, MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[DB1[MLE2[g1], MLE2[g1], MH2]])/(32*MW2*Pi^2) + 
   (C2*EL^2*Re[DB1[MLE2[g1], MLE2[g1], MZ2]])/(8*Pi^2) - 
   (C2*EL^2*S2*Re[DB1[MLE2[g1], MLE2[g1], MZ2]])/(32*Pi^2) - 
   (C2*EL^2*SW^2*Re[DB1[MLE2[g1], MLE2[g1], MZ2]])/(4*Pi^2) - 
   (EL^2*S2*MLE2[g1]*Re[DB1[MLE2[g1], MLE2[g1], MZ2]])/(32*MW2*Pi^2))

dMf1[3, g1_] = (MQU[g1]*(-(C2*EL^2)/(24*Pi^2) - (EL^2*S2)/(32*Pi^2) - 
   (C2*EL^2*S2)/(64*Pi^2) + (C2*EL^2*SW^2)/(18*Pi^2) + 
   (EL^2*S2*MQU2[g1]*Re[B0[MQU2[g1], MH2, MQU2[g1]]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[B0[MQU2[g1], MW2, MQD2[g1]]])/(16*MW2*Pi^2) + 
   (C2*EL^2*Re[B0[MQU2[g1], MZ2, MQU2[g1]]])/(6*Pi^2) - 
   (2*C2*EL^2*SW^2*Re[B0[MQU2[g1], MZ2, MQU2[g1]]])/(9*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[B0[MQU2[g1], MZ2, MQU2[g1]]])/(32*MW2*Pi^2) - 
   (EL^2*S2*Re[B1[MQU2[g1], MQD2[g1], MW2]])/(16*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[B1[MQU2[g1], MQD2[g1], MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[B1[MQU2[g1], MQD2[g1], MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[B1[MQU2[g1], MQU2[g1], MH2]])/(32*MW2*Pi^2) + 
   (C2*EL^2*Re[B1[MQU2[g1], MQU2[g1], MZ2]])/(12*Pi^2) - 
   (C2*EL^2*S2*Re[B1[MQU2[g1], MQU2[g1], MZ2]])/(32*Pi^2) - 
   (C2*EL^2*SW^2*Re[B1[MQU2[g1], MQU2[g1], MZ2]])/(9*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[B1[MQU2[g1], MQU2[g1], MZ2]])/(32*MW2*Pi^2)))/2

dZfL1[3, g1_, g1_] = -(C2*EL^2)/(24*Pi^2) + (EL^2*S2)/(32*Pi^2) + (C2*EL^2*S2)/(64*Pi^2) + 
 (C2*EL^2*SW^2)/(36*Pi^2) + 
 ((EL^2*S2)/(16*Pi^2) + (EL^2*S2*MQD2[g1])/(32*MW2*Pi^2))*
  Re[B1[MQU2[g1], MQD2[g1], MW2]] + 
 (EL^2*S2*MQU2[g1]*Re[B1[MQU2[g1], MQU2[g1], MH2]])/(64*MW2*Pi^2) + 
 (-(C2*EL^2)/(12*Pi^2) + (C2*EL^2*S2)/(32*Pi^2) + (C2*EL^2*SW^2)/(18*Pi^2) + 
   (EL^2*S2*MQU2[g1])/(64*MW2*Pi^2))*Re[B1[MQU2[g1], MQU2[g1], MZ2]] - 
 MQU2[g1]*((EL^2*S2*MQU2[g1]*Re[DB0[MQU2[g1], MH2, MQU2[g1]]])/
    (32*MW2*Pi^2) - (EL^2*S2*MQD2[g1]*Re[DB0[MQU2[g1], MW2, MQD2[g1]]])/
    (16*MW2*Pi^2) + (C2*EL^2*Re[DB0[MQU2[g1], MZ2, MQU2[g1]]])/(6*Pi^2) - 
   (2*C2*EL^2*SW^2*Re[DB0[MQU2[g1], MZ2, MQU2[g1]]])/(9*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[DB0[MQU2[g1], MZ2, MQU2[g1]]])/(32*MW2*Pi^2) - 
   (EL^2*S2*Re[DB1[MQU2[g1], MQD2[g1], MW2]])/(16*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[DB1[MQU2[g1], MQD2[g1], MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[DB1[MQU2[g1], MQD2[g1], MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[DB1[MQU2[g1], MQU2[g1], MH2]])/(32*MW2*Pi^2) + 
   (C2*EL^2*Re[DB1[MQU2[g1], MQU2[g1], MZ2]])/(12*Pi^2) - 
   (C2*EL^2*S2*Re[DB1[MQU2[g1], MQU2[g1], MZ2]])/(32*Pi^2) - 
   (C2*EL^2*SW^2*Re[DB1[MQU2[g1], MQU2[g1], MZ2]])/(9*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[DB1[MQU2[g1], MQU2[g1], MZ2]])/(32*MW2*Pi^2))

dZfR1[3, g1_, g1_] = (C2*EL^2*SW^2)/(36*Pi^2) + (EL^2*S2*MQU2[g1]*Re[B1[MQU2[g1], MQD2[g1], MW2]])/
  (32*MW2*Pi^2) + (EL^2*S2*MQU2[g1]*Re[B1[MQU2[g1], MQU2[g1], MH2]])/
  (64*MW2*Pi^2) + ((C2*EL^2*SW^2)/(18*Pi^2) + (EL^2*S2*MQU2[g1])/
    (64*MW2*Pi^2))*Re[B1[MQU2[g1], MQU2[g1], MZ2]] - 
 MQU2[g1]*((EL^2*S2*MQU2[g1]*Re[DB0[MQU2[g1], MH2, MQU2[g1]]])/
    (32*MW2*Pi^2) - (EL^2*S2*MQD2[g1]*Re[DB0[MQU2[g1], MW2, MQD2[g1]]])/
    (16*MW2*Pi^2) + (C2*EL^2*Re[DB0[MQU2[g1], MZ2, MQU2[g1]]])/(6*Pi^2) - 
   (2*C2*EL^2*SW^2*Re[DB0[MQU2[g1], MZ2, MQU2[g1]]])/(9*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[DB0[MQU2[g1], MZ2, MQU2[g1]]])/(32*MW2*Pi^2) - 
   (EL^2*S2*Re[DB1[MQU2[g1], MQD2[g1], MW2]])/(16*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[DB1[MQU2[g1], MQD2[g1], MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[DB1[MQU2[g1], MQD2[g1], MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[DB1[MQU2[g1], MQU2[g1], MH2]])/(32*MW2*Pi^2) + 
   (C2*EL^2*Re[DB1[MQU2[g1], MQU2[g1], MZ2]])/(12*Pi^2) - 
   (C2*EL^2*S2*Re[DB1[MQU2[g1], MQU2[g1], MZ2]])/(32*Pi^2) - 
   (C2*EL^2*SW^2*Re[DB1[MQU2[g1], MQU2[g1], MZ2]])/(9*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[DB1[MQU2[g1], MQU2[g1], MZ2]])/(32*MW2*Pi^2))

dMf1[4, g1_] = (MQD[g1]*(-(C2*EL^2)/(48*Pi^2) - (EL^2*S2)/(32*Pi^2) - 
   (C2*EL^2*S2)/(64*Pi^2) + (C2*EL^2*SW^2)/(72*Pi^2) + 
   (EL^2*S2*MQD2[g1]*Re[B0[MQD2[g1], MH2, MQD2[g1]]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[B0[MQD2[g1], MW2, MQU2[g1]]])/(16*MW2*Pi^2) + 
   (C2*EL^2*Re[B0[MQD2[g1], MZ2, MQD2[g1]]])/(12*Pi^2) - 
   (C2*EL^2*SW^2*Re[B0[MQD2[g1], MZ2, MQD2[g1]]])/(18*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[B0[MQD2[g1], MZ2, MQD2[g1]]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[B1[MQD2[g1], MQD2[g1], MH2]])/(32*MW2*Pi^2) + 
   (C2*EL^2*Re[B1[MQD2[g1], MQD2[g1], MZ2]])/(24*Pi^2) - 
   (C2*EL^2*S2*Re[B1[MQD2[g1], MQD2[g1], MZ2]])/(32*Pi^2) - 
   (C2*EL^2*SW^2*Re[B1[MQD2[g1], MQD2[g1], MZ2]])/(36*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[B1[MQD2[g1], MQD2[g1], MZ2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*Re[B1[MQD2[g1], MQU2[g1], MW2]])/(16*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[B1[MQD2[g1], MQU2[g1], MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[B1[MQD2[g1], MQU2[g1], MW2]])/(32*MW2*Pi^2)))/2

dZfL1[4, g1_, g1_] = -(C2*EL^2)/(48*Pi^2) + (EL^2*S2)/(32*Pi^2) + (C2*EL^2*S2)/(64*Pi^2) + 
 (C2*EL^2*SW^2)/(144*Pi^2) + 
 (EL^2*S2*MQD2[g1]*Re[B1[MQD2[g1], MQD2[g1], MH2]])/(64*MW2*Pi^2) + 
 (-(C2*EL^2)/(24*Pi^2) + (C2*EL^2*S2)/(32*Pi^2) + (C2*EL^2*SW^2)/(72*Pi^2) + 
   (EL^2*S2*MQD2[g1])/(64*MW2*Pi^2))*Re[B1[MQD2[g1], MQD2[g1], MZ2]] + 
 ((EL^2*S2)/(16*Pi^2) + (EL^2*S2*MQU2[g1])/(32*MW2*Pi^2))*
  Re[B1[MQD2[g1], MQU2[g1], MW2]] - 
 MQD2[g1]*((EL^2*S2*MQD2[g1]*Re[DB0[MQD2[g1], MH2, MQD2[g1]]])/
    (32*MW2*Pi^2) - (EL^2*S2*MQU2[g1]*Re[DB0[MQD2[g1], MW2, MQU2[g1]]])/
    (16*MW2*Pi^2) + (C2*EL^2*Re[DB0[MQD2[g1], MZ2, MQD2[g1]]])/(12*Pi^2) - 
   (C2*EL^2*SW^2*Re[DB0[MQD2[g1], MZ2, MQD2[g1]]])/(18*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[DB0[MQD2[g1], MZ2, MQD2[g1]]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[DB1[MQD2[g1], MQD2[g1], MH2]])/(32*MW2*Pi^2) + 
   (C2*EL^2*Re[DB1[MQD2[g1], MQD2[g1], MZ2]])/(24*Pi^2) - 
   (C2*EL^2*S2*Re[DB1[MQD2[g1], MQD2[g1], MZ2]])/(32*Pi^2) - 
   (C2*EL^2*SW^2*Re[DB1[MQD2[g1], MQD2[g1], MZ2]])/(36*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[DB1[MQD2[g1], MQD2[g1], MZ2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*Re[DB1[MQD2[g1], MQU2[g1], MW2]])/(16*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[DB1[MQD2[g1], MQU2[g1], MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[DB1[MQD2[g1], MQU2[g1], MW2]])/(32*MW2*Pi^2))

dZfR1[4, g1_, g1_] = (C2*EL^2*SW^2)/(144*Pi^2) + 
 (EL^2*S2*MQD2[g1]*Re[B1[MQD2[g1], MQD2[g1], MH2]])/(64*MW2*Pi^2) + 
 ((C2*EL^2*SW^2)/(72*Pi^2) + (EL^2*S2*MQD2[g1])/(64*MW2*Pi^2))*
  Re[B1[MQD2[g1], MQD2[g1], MZ2]] + 
 (EL^2*S2*MQD2[g1]*Re[B1[MQD2[g1], MQU2[g1], MW2]])/(32*MW2*Pi^2) - 
 MQD2[g1]*((EL^2*S2*MQD2[g1]*Re[DB0[MQD2[g1], MH2, MQD2[g1]]])/
    (32*MW2*Pi^2) - (EL^2*S2*MQU2[g1]*Re[DB0[MQD2[g1], MW2, MQU2[g1]]])/
    (16*MW2*Pi^2) + (C2*EL^2*Re[DB0[MQD2[g1], MZ2, MQD2[g1]]])/(12*Pi^2) - 
   (C2*EL^2*SW^2*Re[DB0[MQD2[g1], MZ2, MQD2[g1]]])/(18*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[DB0[MQD2[g1], MZ2, MQD2[g1]]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[DB1[MQD2[g1], MQD2[g1], MH2]])/(32*MW2*Pi^2) + 
   (C2*EL^2*Re[DB1[MQD2[g1], MQD2[g1], MZ2]])/(24*Pi^2) - 
   (C2*EL^2*S2*Re[DB1[MQD2[g1], MQD2[g1], MZ2]])/(32*Pi^2) - 
   (C2*EL^2*SW^2*Re[DB1[MQD2[g1], MQD2[g1], MZ2]])/(36*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[DB1[MQD2[g1], MQD2[g1], MZ2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*Re[DB1[MQD2[g1], MQU2[g1], MW2]])/(16*Pi^2) - 
   (EL^2*S2*MQD2[g1]*Re[DB1[MQD2[g1], MQU2[g1], MW2]])/(32*MW2*Pi^2) - 
   (EL^2*S2*MQU2[g1]*Re[DB1[MQD2[g1], MQU2[g1], MW2]])/(32*MW2*Pi^2))

dMZsq1 = -(CW^2*EL^2*MZ2*S2)/(24*Pi^2) + (CW^2*EL^2*S2*Re[A0[MW2]])/(4*Pi^2) - 
 (EL^2/(16*Pi^2) - (CW^2*EL^2*S2)/(32*Pi^2) - (C2*EL^2*SW^2)/(32*Pi^2))*
  Re[A0[MW2]] + (C2*EL^2*S2*(Re[A0[MH2]] + Re[A0[MZ2]]))/(64*Pi^2) + 
 (C4*EL^2*MW2*S2*Re[B0[MZ2, MH2, MZ2]])/(16*Pi^2) + 
 (C2*EL^2*MW2*S2*Re[B0[MZ2, MW2, MW2]])/(8*Pi^2) - 
 (CW^2*EL^2*MZ2*S2*Re[B0[MZ2, MW2, MW2]])/(2*Pi^2) + 
 (3*C2*EL^2*S2*Re[B00[MZ2, 0, 0]])/(16*Pi^2) - 
 (C2*EL^2*S2*Re[B00[MZ2, MH2, MZ2]])/(16*Pi^2) - 
 (CW^2*EL^2*S2*Re[B00[MZ2, MW2, MW2]])/(2*Pi^2) + 
 (EL^2/(8*Pi^2) - (CW^2*EL^2*S2)/(16*Pi^2) - (C2*EL^2*SW^2)/(16*Pi^2))*
  Re[B00[MZ2, MW2, MW2]] - (3*C2*EL^2*MZ2*S2*Re[B1[MZ2, 0, 0]])/(32*Pi^2) + 
 Sum[((EL^2*S2)/(8*Pi^2) - (C2*EL^2*S2)/(32*Pi^2) - (CW^2*EL^2*S2)/(8*Pi^2) - 
     (C2*EL^2*SW^2)/(8*Pi^2))*Re[A0[MLE2[Gen1]]] - 
   ((EL^2*S2)/(24*Pi^2) + (C2*EL^2*S2)/(96*Pi^2) + (CW^2*EL^2*S2)/(24*Pi^2) + 
     (C2*EL^2*SW^2)/(24*Pi^2))*Re[A0[MQD2[Gen1]]] + 
   ((EL^2*S2)/(12*Pi^2) - (C2*EL^2*S2)/(96*Pi^2) - (CW^2*EL^2*S2)/(6*Pi^2) - 
     (C2*EL^2*SW^2)/(6*Pi^2))*Re[A0[MQU2[Gen1]]] - 
   ((EL^2*MLE2[Gen1])/(4*Pi^2) - (C2*EL^2*MLE2[Gen1])/(8*Pi^2) - 
     (EL^2*S2*MLE2[Gen1])/(8*Pi^2) + (C2*EL^2*S2*MLE2[Gen1])/(32*Pi^2) + 
     (CW^2*EL^2*S2*MLE2[Gen1])/(8*Pi^2) + (C2*EL^2*SW^2*MLE2[Gen1])/(8*Pi^2))*
    Re[B0[MZ2, MLE2[Gen1], MLE2[Gen1]]] - 
   ((EL^2*MQD2[Gen1])/(12*Pi^2) + (C2*EL^2*MQD2[Gen1])/(24*Pi^2) + 
     (EL^2*S2*MQD2[Gen1])/(24*Pi^2) + (C2*EL^2*S2*MQD2[Gen1])/(96*Pi^2) + 
     (CW^2*EL^2*S2*MQD2[Gen1])/(24*Pi^2) + (C2*EL^2*SW^2*MQD2[Gen1])/
      (24*Pi^2))*Re[B0[MZ2, MQD2[Gen1], MQD2[Gen1]]] - 
   ((EL^2*MQU2[Gen1])/(3*Pi^2) - (C2*EL^2*MQU2[Gen1])/(12*Pi^2) - 
     (EL^2*S2*MQU2[Gen1])/(12*Pi^2) + (C2*EL^2*S2*MQU2[Gen1])/(96*Pi^2) + 
     (CW^2*EL^2*S2*MQU2[Gen1])/(6*Pi^2) + (C2*EL^2*SW^2*MQU2[Gen1])/(6*Pi^2))*
    Re[B0[MZ2, MQU2[Gen1], MQU2[Gen1]]] - 
   ((EL^2*S2)/(4*Pi^2) - (C2*EL^2*S2)/(16*Pi^2) - (CW^2*EL^2*S2)/(4*Pi^2) - 
     (C2*EL^2*SW^2)/(4*Pi^2))*Re[B00[MZ2, MLE2[Gen1], MLE2[Gen1]]] + 
   ((EL^2*S2)/(12*Pi^2) + (C2*EL^2*S2)/(48*Pi^2) + (CW^2*EL^2*S2)/(12*Pi^2) + 
     (C2*EL^2*SW^2)/(12*Pi^2))*Re[B00[MZ2, MQD2[Gen1], MQD2[Gen1]]] - 
   ((EL^2*S2)/(6*Pi^2) - (C2*EL^2*S2)/(48*Pi^2) - (CW^2*EL^2*S2)/(3*Pi^2) - 
     (C2*EL^2*SW^2)/(3*Pi^2))*Re[B00[MZ2, MQU2[Gen1], MQU2[Gen1]]] + 
   MZ2*((EL^2*S2)/(8*Pi^2) - (C2*EL^2*S2)/(32*Pi^2) - 
     (CW^2*EL^2*S2)/(8*Pi^2) - (C2*EL^2*SW^2)/(8*Pi^2))*
    Re[B1[MZ2, MLE2[Gen1], MLE2[Gen1]]] - 
   MZ2*((EL^2*S2)/(24*Pi^2) + (C2*EL^2*S2)/(96*Pi^2) + 
     (CW^2*EL^2*S2)/(24*Pi^2) + (C2*EL^2*SW^2)/(24*Pi^2))*
    Re[B1[MZ2, MQD2[Gen1], MQD2[Gen1]]] + 
   MZ2*((EL^2*S2)/(12*Pi^2) - (C2*EL^2*S2)/(96*Pi^2) - 
     (CW^2*EL^2*S2)/(6*Pi^2) - (C2*EL^2*SW^2)/(6*Pi^2))*
    Re[B1[MZ2, MQU2[Gen1], MQU2[Gen1]]], {Gen1, 3}]

dMWsq1 = (EL^2*MW2)/(8*Pi^2) - (EL^2*MW2*S2)/(8*Pi^2) + (CW^2*EL^2*MW2*S2)/(8*Pi^2) - 
 MW2*(EL^2/(24*Pi^2) + (CW^2*EL^2*S2)/(24*Pi^2)) + 
 (EL^2*S2*Re[A0[MH2]])/(64*Pi^2) + (5*EL^2*S2*Re[A0[MW2]])/(32*Pi^2) + 
 ((EL^2*S2)/(64*Pi^2) + (CW^2*EL^2*S2)/(8*Pi^2))*Re[A0[MZ2]] - 
 (EL^2*MW2*Re[B0[MW2, 0, MW2]])/(4*Pi^2) - 
 (CW^2*EL^2*MW2*S2*Re[B0[MW2, MW2, MZ2]])/(2*Pi^2) - 
 ((EL^2*MW2)/(8*Pi^2) - (CW^2*EL^2*MW2*S2)/(16*Pi^2) - 
   (C2*EL^2*MW2*SW^2)/(16*Pi^2))*Re[B0[MW2, MW2, MZ2]] + 
 (EL^2*MW2*S2*(Re[B0[MW2, MH2, MW2]] + Re[B0[MW2, MW2, MZ2]]))/(16*Pi^2) - 
 (EL^2*Re[B00[MW2, 0, MW2]])/(2*Pi^2) - (EL^2*S2*Re[B00[MW2, MH2, MW2]])/
  (16*Pi^2) - ((EL^2*S2)/(16*Pi^2) + (CW^2*EL^2*S2)/(2*Pi^2))*
  Re[B00[MW2, MZ2, MW2]] + Sum[-(EL^2*S2*Re[A0[MLE2[Gen1]]])/(16*Pi^2) - 
   (3*EL^2*S2*Re[A0[MQD2[Gen1]]])/(16*Pi^2) - 
   (3*EL^2*S2*MQU2[Gen1]*Re[B0[MW2, MQD2[Gen1], MQU2[Gen1]]])/(16*Pi^2) + 
   (EL^2*S2*Re[B00[MW2, 0, MLE2[Gen1]]])/(8*Pi^2) + 
   (3*EL^2*S2*Re[B00[MW2, MQU2[Gen1], MQD2[Gen1]]])/(8*Pi^2) - 
   (EL^2*MW2*S2*Re[B1[MW2, 0, MLE2[Gen1]]])/(16*Pi^2) - 
   (3*EL^2*MW2*S2*Re[B1[MW2, MQU2[Gen1], MQD2[Gen1]]])/(16*Pi^2), {Gen1, 3}]

dMHsq1 = (-3*EL^2*MW2*S2)/(16*Pi^2) - (C4*EL^2*MW2*S2)/(16*Pi^2) - 
 (C2*EL^2*MZ2*S2)/(32*Pi^2) + (3*EL^2*MH2*S2*Re[A0[MH2]])/(128*MW2*Pi^2) + 
 ((3*EL^2*S2)/(32*Pi^2) + (EL^2*MH2*S2)/(64*MW2*Pi^2))*Re[A0[MW2]] + 
 ((3*C2*EL^2*S2)/(64*Pi^2) + (EL^2*MH2*S2)/(128*MW2*Pi^2))*Re[A0[MZ2]] + 
 (9*EL^2*MH2^2*S2*Re[B0[MH2, MH2, MH2]])/(128*MW2*Pi^2) - 
 ((EL^2*MH2*S2)/(16*Pi^2) - (EL^2*MH2^2*S2)/(64*MW2*Pi^2) - 
   (3*EL^2*MW2*S2)/(16*Pi^2))*Re[B0[MH2, MW2, MW2]] - 
 ((C2*EL^2*MH2*S2)/(32*Pi^2) - (EL^2*MH2^2*S2)/(128*MW2*Pi^2) - 
   (3*C4*EL^2*MW2*S2)/(32*Pi^2))*Re[B0[MH2, MZ2, MZ2]] + 
 Sum[-(EL^2*S2*MLE2[Gen1]*Re[A0[MLE2[Gen1]]])/(16*MW2*Pi^2) - 
   (3*EL^2*S2*MQD2[Gen1]*Re[A0[MQD2[Gen1]]])/(16*MW2*Pi^2) - 
   (3*EL^2*S2*MQU2[Gen1]*Re[A0[MQU2[Gen1]]])/(16*MW2*Pi^2) - 
   (EL^2*S2*MLE2[Gen1]^2*Re[B0[MH2, MLE2[Gen1], MLE2[Gen1]]])/(8*MW2*Pi^2) - 
   (3*EL^2*S2*MQD2[Gen1]^2*Re[B0[MH2, MQD2[Gen1], MQD2[Gen1]]])/
    (8*MW2*Pi^2) - (3*EL^2*S2*MQU2[Gen1]^2*
     Re[B0[MH2, MQU2[Gen1], MQU2[Gen1]]])/(8*MW2*Pi^2) - 
   (EL^2*MH2*S2*MLE2[Gen1]*Re[B1[MH2, MLE2[Gen1], MLE2[Gen1]]])/
    (16*MW2*Pi^2) - (3*EL^2*MH2*S2*MQD2[Gen1]*
     Re[B1[MH2, MQD2[Gen1], MQD2[Gen1]]])/(16*MW2*Pi^2) - 
   (3*EL^2*MH2*S2*MQU2[Gen1]*Re[B1[MH2, MQU2[Gen1], MQU2[Gen1]]])/
    (16*MW2*Pi^2), {Gen1, 3}]

dSWsq1 = -((CW^2*dMWsq1)/MW2) + (CW^2*dMZsq1)/MZ2

dCWsq1 = -1 + dSWsq1

dSW1 = 1/2 + dSWsq1 + SW^(-1)

dZAA1 = EL^2/(24*Pi^2) + (EL^2*Re[B0[0, MW2, MW2]])/(2*Pi^2) + 
 (3*EL^2*Re[DB00[0, MW2, MW2]])/(4*Pi^2) + 
 Sum[(EL^2*Re[B1[0, MLE2[Gen1], MLE2[Gen1]]])/(4*Pi^2) + 
   (EL^2*Re[B1[0, MQD2[Gen1], MQD2[Gen1]]])/(12*Pi^2) + 
   (EL^2*Re[B1[0, MQU2[Gen1], MQU2[Gen1]]])/(3*Pi^2) - 
   (EL^2*Re[DB00[0, MLE2[Gen1], MLE2[Gen1]]])/(2*Pi^2) - 
   (EL^2*Re[DB00[0, MQD2[Gen1], MQD2[Gen1]]])/(6*Pi^2) - 
   (2*EL^2*Re[DB00[0, MQU2[Gen1], MQU2[Gen1]]])/(3*Pi^2), {Gen1, 3}]

dZAZ1 = 2 + CW^(-1) + dCWsq1 + SW^(-1)

dZZA1 = 0

dZZZ1 = dZAA1 - (dCWsq1*(CW^2 - SW^2))/(CW^2*SW^2)

dZW1 = dZAA1 - dCWsq1/SW^2

dZH1 = dZW1 + dMWsq1/MW2

dZphi1 = dZH1

dZchi1 = dZH1

dWFZ1 = -dZZZ1 + (CW^2*EL^2*S2)/(24*Pi^2) + (CW^2*EL^2*S2*Re[B0[MZ2, MW2, MW2]])/
  (2*Pi^2) + (3*C2*EL^2*S2*Re[B1[MZ2, 0, 0]])/(32*Pi^2) - 
 (C4*EL^2*MW2*S2*Re[DB0[MZ2, MH2, MZ2]])/(16*Pi^2) - 
 (C2*EL^2*MW2*S2*Re[DB0[MZ2, MW2, MW2]])/(8*Pi^2) + 
 (CW^2*EL^2*MZ2*S2*Re[DB0[MZ2, MW2, MW2]])/(2*Pi^2) - 
 (3*C2*EL^2*S2*Re[DB00[MZ2, 0, 0]])/(16*Pi^2) + 
 (C2*EL^2*S2*Re[DB00[MZ2, MH2, MZ2]])/(16*Pi^2) + 
 (CW^2*EL^2*S2*Re[DB00[MZ2, MW2, MW2]])/(2*Pi^2) - 
 (EL^2/(8*Pi^2) - (CW^2*EL^2*S2)/(16*Pi^2) - (C2*EL^2*SW^2)/(16*Pi^2))*
  Re[DB00[MZ2, MW2, MW2]] + (3*C2*EL^2*MZ2*S2*Re[DB1[MZ2, 0, 0]])/(32*Pi^2) + 
 Sum[-(((EL^2*S2)/(8*Pi^2) - (C2*EL^2*S2)/(32*Pi^2) - 
      (CW^2*EL^2*S2)/(8*Pi^2) - (C2*EL^2*SW^2)/(8*Pi^2))*
     Re[B1[MZ2, MLE2[Gen1], MLE2[Gen1]]]) + 
   ((EL^2*S2)/(24*Pi^2) + (C2*EL^2*S2)/(96*Pi^2) + (CW^2*EL^2*S2)/(24*Pi^2) + 
     (C2*EL^2*SW^2)/(24*Pi^2))*Re[B1[MZ2, MQD2[Gen1], MQD2[Gen1]]] - 
   ((EL^2*S2)/(12*Pi^2) - (C2*EL^2*S2)/(96*Pi^2) - (CW^2*EL^2*S2)/(6*Pi^2) - 
     (C2*EL^2*SW^2)/(6*Pi^2))*Re[B1[MZ2, MQU2[Gen1], MQU2[Gen1]]] + 
   ((EL^2*MLE2[Gen1])/(4*Pi^2) - (C2*EL^2*MLE2[Gen1])/(8*Pi^2) - 
     (EL^2*S2*MLE2[Gen1])/(8*Pi^2) + (C2*EL^2*S2*MLE2[Gen1])/(32*Pi^2) + 
     (CW^2*EL^2*S2*MLE2[Gen1])/(8*Pi^2) + (C2*EL^2*SW^2*MLE2[Gen1])/(8*Pi^2))*
    Re[DB0[MZ2, MLE2[Gen1], MLE2[Gen1]]] + 
   ((EL^2*MQD2[Gen1])/(12*Pi^2) + (C2*EL^2*MQD2[Gen1])/(24*Pi^2) + 
     (EL^2*S2*MQD2[Gen1])/(24*Pi^2) + (C2*EL^2*S2*MQD2[Gen1])/(96*Pi^2) + 
     (CW^2*EL^2*S2*MQD2[Gen1])/(24*Pi^2) + (C2*EL^2*SW^2*MQD2[Gen1])/
      (24*Pi^2))*Re[DB0[MZ2, MQD2[Gen1], MQD2[Gen1]]] + 
   ((EL^2*MQU2[Gen1])/(3*Pi^2) - (C2*EL^2*MQU2[Gen1])/(12*Pi^2) - 
     (EL^2*S2*MQU2[Gen1])/(12*Pi^2) + (C2*EL^2*S2*MQU2[Gen1])/(96*Pi^2) + 
     (CW^2*EL^2*S2*MQU2[Gen1])/(6*Pi^2) + (C2*EL^2*SW^2*MQU2[Gen1])/(6*Pi^2))*
    Re[DB0[MZ2, MQU2[Gen1], MQU2[Gen1]]] + 
   ((EL^2*S2)/(4*Pi^2) - (C2*EL^2*S2)/(16*Pi^2) - (CW^2*EL^2*S2)/(4*Pi^2) - 
     (C2*EL^2*SW^2)/(4*Pi^2))*Re[DB00[MZ2, MLE2[Gen1], MLE2[Gen1]]] - 
   ((EL^2*S2)/(12*Pi^2) + (C2*EL^2*S2)/(48*Pi^2) + (CW^2*EL^2*S2)/(12*Pi^2) + 
     (C2*EL^2*SW^2)/(12*Pi^2))*Re[DB00[MZ2, MQD2[Gen1], MQD2[Gen1]]] + 
   ((EL^2*S2)/(6*Pi^2) - (C2*EL^2*S2)/(48*Pi^2) - (CW^2*EL^2*S2)/(3*Pi^2) - 
     (C2*EL^2*SW^2)/(3*Pi^2))*Re[DB00[MZ2, MQU2[Gen1], MQU2[Gen1]]] - 
   MZ2*((EL^2*S2)/(8*Pi^2) - (C2*EL^2*S2)/(32*Pi^2) - 
     (CW^2*EL^2*S2)/(8*Pi^2) - (C2*EL^2*SW^2)/(8*Pi^2))*
    Re[DB1[MZ2, MLE2[Gen1], MLE2[Gen1]]] + 
   MZ2*((EL^2*S2)/(24*Pi^2) + (C2*EL^2*S2)/(96*Pi^2) + 
     (CW^2*EL^2*S2)/(24*Pi^2) + (C2*EL^2*SW^2)/(24*Pi^2))*
    Re[DB1[MZ2, MQD2[Gen1], MQD2[Gen1]]] - 
   MZ2*((EL^2*S2)/(12*Pi^2) - (C2*EL^2*S2)/(96*Pi^2) - 
     (CW^2*EL^2*S2)/(6*Pi^2) - (C2*EL^2*SW^2)/(6*Pi^2))*
    Re[DB1[MZ2, MQU2[Gen1], MQU2[Gen1]]], {Gen1, 3}]

dWFW1 = -dZW1 + EL^2/(24*Pi^2) + (CW^2*EL^2*S2)/(24*Pi^2) + 
 (EL^2*Re[B0[MW2, 0, MW2]])/(2*Pi^2) + (CW^2*EL^2*S2*Re[B0[MW2, MW2, MZ2]])/
  (2*Pi^2) + (EL^2*MW2*Re[DB0[MW2, 0, MW2]])/(4*Pi^2) + 
 (CW^2*EL^2*MW2*S2*Re[DB0[MW2, MW2, MZ2]])/(2*Pi^2) + 
 ((EL^2*MW2)/(8*Pi^2) - (CW^2*EL^2*MW2*S2)/(16*Pi^2) - 
   (C2*EL^2*MW2*SW^2)/(16*Pi^2))*Re[DB0[MW2, MW2, MZ2]] - 
 (EL^2*MW2*S2*(Re[DB0[MW2, MH2, MW2]] + Re[DB0[MW2, MW2, MZ2]]))/(16*Pi^2) + 
 (EL^2*Re[DB00[MW2, 0, MW2]])/(2*Pi^2) + (EL^2*S2*Re[DB00[MW2, MH2, MW2]])/
  (16*Pi^2) + ((EL^2*S2)/(16*Pi^2) + (CW^2*EL^2*S2)/(2*Pi^2))*
  Re[DB00[MW2, MZ2, MW2]] + 
 Sum[(EL^2*S2*Re[B1[MW2, 0, MLE2[Gen1]]])/(16*Pi^2) + 
   (3*EL^2*S2*Re[B1[MW2, MQU2[Gen1], MQD2[Gen1]]])/(16*Pi^2) + 
   (3*EL^2*S2*MQU2[Gen1]*Re[DB0[MW2, MQD2[Gen1], MQU2[Gen1]]])/(16*Pi^2) - 
   (EL^2*S2*Re[DB00[MW2, 0, MLE2[Gen1]]])/(8*Pi^2) - 
   (3*EL^2*S2*Re[DB00[MW2, MQU2[Gen1], MQD2[Gen1]]])/(8*Pi^2) + 
   (EL^2*MW2*S2*Re[DB1[MW2, 0, MLE2[Gen1]]])/(16*Pi^2) + 
   (3*EL^2*MW2*S2*Re[DB1[MW2, MQU2[Gen1], MQD2[Gen1]]])/(16*Pi^2), {Gen1, 3}]

dZe1 = 1/2 - dZAA1 - (dZZA1*SW)/CW

dTad1 = (EL*MW^3)/(8*Pi^2*SW) + (C2*EL*MW*MZ2)/(16*Pi^2*SW) - 
 (3*EL*MH2*Re[A0[MH2]])/(64*MW*Pi^2*SW) - 
 ((EL*MH2)/(32*MW*Pi^2*SW) + (3*EL*MW)/(16*Pi^2*SW))*Re[A0[MW2]] - 
 ((EL*MH2)/(64*MW*Pi^2*SW) + (3*C2*EL*MW)/(32*Pi^2*SW))*Re[A0[MZ2]] + 
 Sum[(EL*MLE2[Gen1]*Re[A0[MLE2[Gen1]]])/(8*MW*Pi^2*SW) + 
   (3*EL*MQD2[Gen1]*Re[A0[MQD2[Gen1]]])/(8*MW*Pi^2*SW) + 
   (3*EL*MQU2[Gen1]*Re[A0[MQU2[Gen1]]])/(8*MW*Pi^2*SW), {Gen1, 3}]

GammaHMH = (3*EL^2*MH2*S2*Im[A0[MH2]])/(128*MW2*Pi^2) + 
 ((3*EL^2*S2)/(32*Pi^2) + (EL^2*MH2*S2)/(64*MW2*Pi^2))*Im[A0[MW2]] + 
 ((3*C2*EL^2*S2)/(64*Pi^2) + (EL^2*MH2*S2)/(128*MW2*Pi^2))*Im[A0[MZ2]] + 
 (9*EL^2*MH2^2*S2*Im[B0[MH2, MH2, MH2]])/(128*MW2*Pi^2) - 
 ((EL^2*MH2*S2)/(16*Pi^2) - (EL^2*MH2^2*S2)/(64*MW2*Pi^2) - 
   (3*EL^2*MW2*S2)/(16*Pi^2))*Im[B0[MH2, MW2, MW2]] - 
 ((C2*EL^2*MH2*S2)/(32*Pi^2) - (EL^2*MH2^2*S2)/(128*MW2*Pi^2) - 
   (3*C4*EL^2*MW2*S2)/(32*Pi^2))*Im[B0[MH2, MZ2, MZ2]] + 
 Sum[-(EL^2*S2*Im[A0[MLE2[Gen1]]]*MLE2[Gen1])/(16*MW2*Pi^2) - 
   (EL^2*MH2*S2*Im[B1[MH2, MLE2[Gen1], MLE2[Gen1]]]*MLE2[Gen1])/
    (16*MW2*Pi^2) - (EL^2*S2*Im[B0[MH2, MLE2[Gen1], MLE2[Gen1]]]*
     MLE2[Gen1]^2)/(8*MW2*Pi^2) - (3*EL^2*S2*Im[A0[MQD2[Gen1]]]*MQD2[Gen1])/
    (16*MW2*Pi^2) - (3*EL^2*MH2*S2*Im[B1[MH2, MQD2[Gen1], MQD2[Gen1]]]*
     MQD2[Gen1])/(16*MW2*Pi^2) - 
   (3*EL^2*S2*Im[B0[MH2, MQD2[Gen1], MQD2[Gen1]]]*MQD2[Gen1]^2)/
    (8*MW2*Pi^2) - (3*EL^2*S2*Im[A0[MQU2[Gen1]]]*MQU2[Gen1])/(16*MW2*Pi^2) - 
   (3*EL^2*MH2*S2*Im[B1[MH2, MQU2[Gen1], MQU2[Gen1]]]*MQU2[Gen1])/
    (16*MW2*Pi^2) - (3*EL^2*S2*Im[B0[MH2, MQU2[Gen1], MQU2[Gen1]]]*
     MQU2[Gen1]^2)/(8*MW2*Pi^2), {Gen1, 3}]
