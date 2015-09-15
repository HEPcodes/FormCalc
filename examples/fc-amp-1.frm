#-
#:SmallSize 5000000
#:LargeSize 20000000
#:WorkSpace 50000000
#:MaxTermSize 300000
#:TermsInSmall 30000
off stats;
format 75;

#define Dim "D"
#define InsertionPolicy "Default"
#define IntMax "2"
#define ExtMax "2"
#define MomElim "Automatic"
#define DotExpand "0"
#define NoCostly "0"
#define HaveFermions "1"
#define FermionOrder "Automatic"
#define FermionChains "Weyl"
#define Gamma5Test "0"
#define Evanescent "0"
#define HaveSUN "1"
#define SUNN "3"
#define SortDen "1"
#define CombineDen "Automatic"
#define PaVeReduce "False"
#define CancelQ2 "1"
#define OPP "6"
#define OPPMethod "Ninja"

s Alfa, Alfa2, Alfas, Alfas2, C2B, CA, CA2, CAB, CB, CB2, CBA, CBA2,
  cG1, cG2, cG3, cG4, cR1, CW, CW2, D, dMZsq1, dZAA1, dZAZ1, dZG01,
  dZZA1, dZZZ1, EL, GS, MA0, MA02, MA0tree, MA0tree2, MB, MB2, MC, MC2,
  MD, MD2, ME, ME2, MGl, MGl2, MH, Mh0, Mh02, Mh0tree, Mh0tree2, MH2,
  MHH, MHH2, MHHtree, MHHtree2, MHp, MHp2, MHptree, MHptree2, Mino3,
  Mino3C, ML, ML2, MM, MM2, MS, MS2, MT, MT2, MU, MU2, MUE, MUEC, MW,
  MW2, MZ, MZ2, pG1, pG2, pG3, pG4, pG5, pG6, pG7, pG8, Pi, pM1, pM2,
  pM3, pM4, S, S2B, SA, SA2, SAB, SB, SB2, SBA, SBA2, SqrtEGl,
  SqrtEGlC, SW, SW2, T, TB, TB2, U;
cf Af, AfC, Den, dirM, Hel, IndexDelta, intM, List, MASf, MASf2, MCha,
  MCha2, Mf, Mf2, MHiggs, MHiggs2, MHiggstree, MHiggstree2, MNeu,
  MNeu2, MSf, MSf2, powM, SumOver, SUNT, UCha, UChaC, USf, USfC, VCha,
  VChaC, ZNeu, ZNeuC;
f Spinor;
d D;
i Cha5 = 2, Cha6 = 2, Col3 = 3, Col4 = 3, Gen5 = 3, Lor1, Lor2,
  Lor3, Lor4, Lor5, Lor6, Lor7, Lor8, Neu5 = 4, Neu6 = 4, Sfe5 = 2,
  Sfe6 = 2, sM = 0;

#define Vectors "q1, eta1, eta2, eta3, eta4, e1, e2, e3, e4, ec1, ec2, ec3, ec4, z1, z2, z3, z4, zc1, zc2, zc3, zc4, k1, k2, k3, k4"
v `Vectors';

#define Tensors "eT1, eT2, eT3, eT4, eTc1, eTc2, eTc3, eTc4"
t `Tensors';

table HEL(1:4);
fill HEL(1) = Hel(1);
fill HEL(2) = Hel(2);
fill HEL(3) = Hel(3);
fill HEL(4) = Hel(4);


L T1G1 = (i_*cR1*pG5*d_(Lor1, Lor2)*d_(Lor2, Lor4)*
   d_(Lor3, Lor4)*Den(-k3 - k4, pM3^2)*Den(k3 + k4, pM2^2)*
   intM(Den(q1, pM1^2))*Spinor(k1, 0, -1)*g_(100)*
    (pG2*g_(100, Lor1)*(g6_(100)/2) + 
     pG1*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
   Spinor(k4, MT, 1)*g_(101)*(pG4*g_(101, Lor3)*(g6_(101)/2) + 
     pG3*g_(101, Lor3)*(g7_(101)/2))*Spinor(k3, MT, -1))/
  (32*Pi^2);

L T1G2 = -(i_*d_(Lor1, Lor2)*d_(Lor3, Lor4)*d_(Lor5, Lor6)*
    (pG5*d_(Lor2, Lor6)*d_(Lor4, Lor5) + pG5*d_(Lor2, Lor5)*
      d_(Lor4, Lor6) + pG6*d_(Lor2, Lor4)*d_(Lor5, Lor6))*
    Den(-k3 - k4, pM2^2)*Den(k3 + k4, pM1^2)*intM(Den(q1, MW2))*
    Spinor(k1, 0, -1)*g_(100)*(pG2*g_(100, Lor1)*(g6_(100)/2) + 
      pG1*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
    Spinor(k4, MT, 1)*g_(101)*(pG4*g_(101, Lor3)*(g6_(101)/2) + 
      pG3*g_(101, Lor3)*(g7_(101)/2))*Spinor(k3, MT, -1)*
    SUNT(Col3, Col4))/(16*Pi^2);

L T2G1 = -(cR1*d_(Lor1, Lor2)*Den(-k3 - k4, pM3^2)*
    Den(k3 + k4, pM4^2)*intM(Den(q1, pM1^2), Den(-k3 - k4 + q1, pM2^2))*
    Spinor(k1, 0, -1)*g_(100)*(pG4*g_(100, Lor1)*(g6_(100)/2) + 
      pG3*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
    Spinor(k4, MT, 1)*g_(101)*((pG8*g6_(101))/2 + 
      (pG7*g7_(101))/2)*Spinor(k3, MT, -1)*
    g_(1)*(pM1 - g_(1, q1))*((pG6*g6_(1))/2 + (pG5*g7_(1))/2)*
     (pM2 + g_(1, k3) + g_(1, k4) - g_(1, q1))*
     (cG2*g_(1, Lor2)*(g6_(1)/2) + 
      cG1*g_(1, Lor2)*(g7_(1)/2)))/(32*Pi^2);

L T2G2 = (cR1*pG3*pG6*d_(Lor1, Lor2)*Den(-k3 - k4, pM1^2)*
   Den(k3 + k4, pM4^2)*intM(Den(q1, pM2^2), Den(-k3 - k4 + q1, pM3^2))*
   Spinor(k1, 0, -1)*g_(100)*(pG2*g_(100, Lor1)*(g6_(100)/2) + 
     pG1*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
   Spinor(k4, MT, 1)*g_(101)*((pG5*g6_(101))/2 + 
     (pG4*g7_(101))/2)*Spinor(k3, MT, -1)*(k3(Lor2) + k4(Lor2) - 
    2*q1(Lor2)))/(32*Pi^2);

L T2G3 = -(pG5*pG6*d_(Lor1, Lor2)*Den(-k3 - k4, pM1^2)*
    Den(k3 + k4, pM2^2)*intM(Den(q1, MW2), Den(-k3 - k4 + q1, MW2))*
    Spinor(k1, 0, -1)*g_(100)*(pG2*g_(100, Lor1)*(g6_(100)/2) + 
      pG1*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
    Spinor(k4, MT, 1)*g_(101)*((pG4*g6_(101))/2 + 
      (pG3*g7_(101))/2)*Spinor(k3, MT, -1)*(-k3(Lor2) - k4(Lor2) + 
     q1(Lor2))*SUNT(Col3, Col4))/(16*Pi^2);

L T2G4 = (pG3*pG5*d_(Lor1, Lor2)*d_(Lor3, Lor4)*
   d_(Lor4, Lor6)*d_(Lor5, Lor6)*Den(-k3 - k4, pM1^2)*
   Den(k3 + k4, pM2^2)*intM(Den(q1, MW2), Den(-k3 - k4 + q1, MW2))*
   Spinor(k1, 0, -1)*g_(100)*(pG2*g_(100, Lor1)*(g6_(100)/2) + 
     pG1*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
   Spinor(k4, MT, 1)*g_(101)*((pG4*g6_(101))/2 + 
     (pG4*g7_(101))/2)*Spinor(k3, MT, -1)*
   (d_(Lor3, Lor5)*(-k3(Lor2) - k4(Lor2) + 2*q1(Lor2)) + 
    d_(Lor2, Lor5)*(2*k3(Lor3) + 2*k4(Lor3) - q1(Lor3)) + 
    d_(Lor2, Lor3)*(-k3(Lor5) - k4(Lor5) - q1(Lor5)))*
   SUNT(Col3, Col4))/(16*Pi^2);

L T2G5 = (pG3*pG6*d_(Lor1, Lor2)*d_(Lor2, Lor3)*
   d_(Lor3, Lor4)*Den(-k3 - k4, pM1^2)*Den(k3 + k4, pM3^2)*
   intM(Den(q1, pM2^2), Den(-k3 - k4 + q1, pM4^2))*
   Spinor(k1, 0, -1)*g_(100)*(pG2*g_(100, Lor1)*(g6_(100)/2) + 
     pG1*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
   Spinor(k4, MT, 1)*g_(101)*((pG5*g6_(101))/2 + 
     (pG4*g7_(101))/2)*Spinor(k3, MT, -1)*(k3(Lor4) + k4(Lor4) + 
    q1(Lor4))*SUNT(Col3, Col4))/(16*Pi^2);

L T2G6 = (cR1*d_(Lor1, Lor2)*d_(Lor3, Lor4)*
   Den(-k3 - k4, pM4^2)*Den(k3 + k4, pM3^2)*intM(Den(q1, pM1^2), 
    Den(-k3 - k4 + q1, pM2^2))*Spinor(k1, 0, -1)*g_(100)*
    (pG6*g_(100, Lor1)*(g6_(100)/2) + 
     pG5*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
   Spinor(k4, MT, 1)*g_(101)*(pG8*g_(101, Lor3)*(g6_(101)/2) + 
     pG7*g_(101, Lor3)*(g7_(101)/2))*Spinor(k3, MT, -1)*
   g_(1)*(pM1 - g_(1, q1))*(cG2*g_(1, Lor4)*(g6_(1)/2) + 
     cG1*g_(1, Lor4)*(g7_(1)/2))*(pM2 + g_(1, k3) + 
     g_(1, k4) - g_(1, q1))*(cG4*g_(1, Lor2)*(g6_(1)/2) + 
     cG3*g_(1, Lor2)*(g7_(1)/2)))/(32*Pi^2);

L T2G7 = -(cR1*pG5*pG6*d_(Lor1, Lor2)*d_(Lor3, Lor4)*
    Den(-k3 - k4, pM4^2)*Den(k3 + k4, pM3^2)*intM(Den(q1, pM1^2), 
     Den(-k3 - k4 + q1, pM2^2))*Spinor(k1, 0, -1)*g_(100)*
     (pG2*g_(100, Lor1)*(g6_(100)/2) + 
      pG1*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
    Spinor(k4, MT, 1)*g_(101)*(pG4*g_(101, Lor3)*(g6_(101)/2) + 
      pG3*g_(101, Lor3)*(g7_(101)/2))*Spinor(k3, MT, -1)*
    (k3(Lor2) + k4(Lor2) - 2*q1(Lor2))*(-k3(Lor4) - k4(Lor4) + 
     2*q1(Lor4)))/(32*Pi^2);

L T2G8 = (pG5*pG6*d_(Lor1, Lor2)*d_(Lor3, Lor4)*
   Den(-k3 - k4, pM2^2)*Den(k3 + k4, pM1^2)*intM(Den(q1, MW2), 
    Den(-k3 - k4 + q1, MW2))*Spinor(k1, 0, -1)*g_(100)*
    (pG2*g_(100, Lor1)*(g6_(100)/2) + 
     pG1*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
   Spinor(k4, MT, 1)*g_(101)*(pG4*g_(101, Lor3)*(g6_(101)/2) + 
     pG3*g_(101, Lor3)*(g7_(101)/2))*Spinor(k3, MT, -1)*
   (-k3(Lor2) - k4(Lor2) + q1(Lor2))*q1(Lor4)*SUNT(Col3, Col4))/
  (16*Pi^2);

L T2G9 = -(pG5*pG6*d_(Lor1, Lor2)*d_(Lor3, Lor4)*
    d_(Lor5, Lor6)*d_(Lor7, Lor8)*Den(-k3 - k4, pM2^2)*
    Den(k3 + k4, pM1^2)*intM(Den(q1, MW2), Den(-k3 - k4 + q1, MW2))*
    Spinor(k1, 0, -1)*g_(100)*(pG2*g_(100, Lor1)*(g6_(100)/2) + 
      pG1*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
    Spinor(k4, MT, 1)*g_(101)*(pG4*g_(101, Lor3)*(g6_(101)/2) + 
      pG3*g_(101, Lor3)*(g7_(101)/2))*Spinor(k3, MT, -1)*
    (d_(Lor5, Lor7)*(-k3(Lor2) - k4(Lor2) + 2*q1(Lor2)) + 
     d_(Lor2, Lor7)*(2*k3(Lor5) + 2*k4(Lor5) - q1(Lor5)) + 
     d_(Lor2, Lor5)*(-k3(Lor7) - k4(Lor7) - q1(Lor7)))*
    (d_(Lor6, Lor8)*(k3(Lor4) + k4(Lor4) - 2*q1(Lor4)) + 
     d_(Lor4, Lor8)*(-2*k3(Lor6) - 2*k4(Lor6) + q1(Lor6)) + 
     d_(Lor4, Lor6)*(k3(Lor8) + k4(Lor8) + q1(Lor8)))*
    SUNT(Col3, Col4))/(16*Pi^2);

L T2G10 = (pG5*pG6*d_(Lor1, Lor2)*d_(Lor2, Lor5)*
   d_(Lor3, Lor4)*d_(Lor4, Lor6)*d_(Lor5, Lor6)*
   Den(-k3 - k4, pM3^2)*Den(k3 + k4, pM2^2)*intM(Den(q1, pM1^2), 
    Den(-k3 - k4 + q1, pM4^2))*Spinor(k1, 0, -1)*g_(100)*
    (pG2*g_(100, Lor1)*(g6_(100)/2) + 
     pG1*g_(100, Lor1)*(g7_(100)/2))*Spinor(k2, 0, 1)*
   Spinor(k4, MT, 1)*g_(101)*(pG4*g_(101, Lor3)*(g6_(101)/2) + 
     pG3*g_(101, Lor3)*(g7_(101)/2))*Spinor(k3, MT, -1)*
   SUNT(Col3, Col4))/(16*Pi^2);

L T3G1 = i_*d_(Lor1, Lor2)*Den(-k3 - k4, pM1^2)*
  Den(k3 + k4, MZ2)*(pG3*(-k3(Lor1) - k4(Lor1)) + 
   pG4*(k3(Lor1) + k4(Lor1)))*Spinor(k1, 0, -1)*g_(100)*
   (pG2*g_(100, Lor2)*(g6_(100)/2) + 
    pG1*g_(100, Lor2)*(g7_(100)/2))*Spinor(k2, 0, 1)*
  Spinor(k4, MT, 1)*g_(101)*((EL*MT*g6_(101))/(4*MW*SW) - 
    (EL*MT*g7_(101))/(4*MW*SW))*Spinor(k3, MT, -1)*SUNT(Col3, Col4);

L T3G2 = -(i_*d_(Lor1, Lor2)*d_(Lor3, Lor4)*
   Den(-k3 - k4, pM1^2)*Den(k3 + k4, pM2^2)*(pG6*d_(Lor1, Lor3) + 
    pG7*d_(Lor1, Lor3)*(-k3.k3 - 2*k3.k4 - k4.k4) + 
    pG5*(-k3(Lor1) - k4(Lor1))*(k3(Lor3) + k4(Lor3)))*
   Spinor(k1, 0, -1)*g_(100)*(pG2*g_(100, Lor2)*(g6_(100)/2) + 
     pG1*g_(100, Lor2)*(g7_(100)/2))*Spinor(k2, 0, 1)*
   Spinor(k4, MT, 1)*g_(101)*(pG4*g_(101, Lor4)*(g6_(101)/2) + 
     pG3*g_(101, Lor4)*(g7_(101)/2))*Spinor(k3, MT, -1)*
   SUNT(Col3, Col4));


#define Legs "4"
#define Invariants "S, T, U"
#define OnShell "True"
#define Antisymmetrize "1"

#define MomSum "k1 + k2 - k3 - k4"
#define k1 "-k2 + k3 + k4"
#define k2 "-k1 + k3 + k4"
#define k3 "k1 + k2 - k4"
#define k4 "k1 + k2 - k3"

#procedure Neglect
id ME = 0;
id ME2 = 0;
#endprocedure

#procedure Square
id EL^2 = 4*Alfa*Pi;
id EL^-2 = 1/(4*Alfa*Pi);
id GS^2 = 4*Alfas*Pi;
id GS^-2 = 1/(4*Alfas*Pi);
id SqrtEGl^2 = Mino3/MGl;
id SqrtEGl^-2 = 1/(Mino3/MGl);
id SqrtEGlC^2 = Mino3C/MGl;
id SqrtEGlC^-2 = 1/(Mino3C/MGl);
id Alfa^2 = Alfa2;
id Alfa^-2 = Alfa2^-1;
id Alfas^2 = Alfas2;
id Alfas^-2 = Alfas2^-1;
id CA^2 = CA2;
id CA^-2 = CA2^-1;
id CB^2 = CB2;
id CB^-2 = CB2^-1;
id CBA^2 = CBA2;
id CBA^-2 = CBA2^-1;
id CW^2 = CW2;
id CW^-2 = CW2^-1;
id MA0^2 = MA02;
id MA0^-2 = MA02^-1;
id MA0tree^2 = MA0tree2;
id MA0tree^-2 = MA0tree2^-1;
id MB^2 = MB2;
id MB^-2 = MB2^-1;
id MC^2 = MC2;
id MC^-2 = MC2^-1;
id MD^2 = MD2;
id MD^-2 = MD2^-1;
id ME^2 = ME2;
id ME^-2 = ME2^-1;
id MGl^2 = MGl2;
id MGl^-2 = MGl2^-1;
id MH^2 = MH2;
id MH^-2 = MH2^-1;
id Mh0^2 = Mh02;
id Mh0^-2 = Mh02^-1;
id Mh0tree^2 = Mh0tree2;
id Mh0tree^-2 = Mh0tree2^-1;
id MHH^2 = MHH2;
id MHH^-2 = MHH2^-1;
id MHHtree^2 = MHHtree2;
id MHHtree^-2 = MHHtree2^-1;
id MHp^2 = MHp2;
id MHp^-2 = MHp2^-1;
id MHptree^2 = MHptree2;
id MHptree^-2 = MHptree2^-1;
id ML^2 = ML2;
id ML^-2 = ML2^-1;
id MM^2 = MM2;
id MM^-2 = MM2^-1;
id MS^2 = MS2;
id MS^-2 = MS2^-1;
id MT^2 = MT2;
id MT^-2 = MT2^-1;
id MU^2 = MU2;
id MU^-2 = MU2^-1;
id MW^2 = MW2;
id MW^-2 = MW2^-1;
id MZ^2 = MZ2;
id MZ^-2 = MZ2^-1;
id SA^2 = SA2;
id SA^-2 = SA2^-1;
id SB^2 = SB2;
id SB^-2 = SB2^-1;
id SBA^2 = SBA2;
id SBA^-2 = SBA2^-1;
id SW^2 = SW2;
id SW^-2 = SW2^-1;
id TB^2 = TB2;
id TB^-2 = TB2^-1;
id Mf(?a)^2 = Mf2(?a);
id Mf(?a)^-2 = Mf2(?a)^-1;
id MSf(?a)^2 = MSf2(?a);
id MSf(?a)^-2 = MSf2(?a)^-1;
id MASf(?a)^2 = MASf2(?a);
id MASf(?a)^-2 = MASf2(?a)^-1;
id MCha(?a)^2 = MCha2(?a);
id MCha(?a)^-2 = MCha2(?a)^-1;
id MNeu(?a)^2 = MNeu2(?a);
id MNeu(?a)^-2 = MNeu2(?a)^-1;
id MHiggs(?a)^2 = MHiggs2(?a);
id MHiggs(?a)^-2 = MHiggs2(?a)^-1;
id MHiggstree(?a)^2 = MHiggstree2(?a);
id MHiggstree(?a)^-2 = MHiggstree2(?a)^-1;
#endprocedure

#procedure eiei
id e1.ec1 = -1;
id e2.ec2 = -1;
id e3.ec3 = -1;
id e4.ec4 = -1;
#endprocedure

#procedure eiki
id e1.k1 = 0;
id ec1.k1 = 0;
id e2.k2 = 0;
id ec2.k2 = 0;
id e3.k3 = 0;
id ec3.k3 = 0;
id e4.k4 = 0;
id ec4.k4 = 0;
#endprocedure

#procedure kikj
id k1.k1 = ME2;
id k1.k2 = -ME2 + S/2;
id k1.k3 = ME2/2 + MT2/2 - T/2;
id k1.k4 = ME2/2 + MT2/2 - U/2;
id k2.k2 = ME2;
id k2.k3 = ME2/2 + MT2/2 - U/2;
id k2.k4 = ME2/2 + MT2/2 - T/2;
id k3.k3 = MT2;
id k3.k4 = -MT2 + S/2;
id k4.k4 = MT2;
#endprocedure

#procedure InvSimplify(foo)
id `foo'(ARG?) = `foo'(ARG, ARG*replace_(S, 2*ME2 + 2*MT2 - T - U));
#call Fewest(`foo')
id `foo'(ARG?) = `foo'(ARG, ARG*replace_(T, 2*ME2 + 2*MT2 - S - U));
#call Fewest(`foo')
id `foo'(ARG?) = `foo'(ARG, ARG*replace_(U, 2*ME2 + 2*MT2 - S - T));
#call Fewest(`foo')
#endprocedure

#procedure ConstBracket
ab Af, AfC, Alfa, Alfa2, Alfas, Alfas2, C2B, CA, CA2, CAB, CB, CB2,
  CBA, CBA2, CW, CW2, dMZsq1, dZAA1, dZAZ1, dZG01, dZZA1, dZZZ1, EL,
  GS, List, MA0, MA02, MA0tree, MA0tree2, MASf, MASf2, MB, MB2, MC,
  MC2, MCha, MCha2, MD, MD2, ME, ME2, Mf, Mf2, MGl, MGl2, MH, Mh0,
  Mh02, Mh0tree, Mh0tree2, MH2, MHH, MHH2, MHHtree, MHHtree2, MHiggs,
  MHiggs2, MHiggstree, MHiggstree2, MHp, MHp2, MHptree, MHptree2,
  Mino3, Mino3C, ML, ML2, MM, MM2, MNeu, MNeu2, MS, MS2, MSf, MSf2, MT,
  MT2, MU, MU2, MUE, MUEC, MW, MW2, MZ, MZ2, Pi, S2B, SA, SA2, SAB, SB,
  SB2, SBA, SBA2, SqrtEGl, SqrtEGlC, SW, SW2, TB, TB2, UCha, UChaC,
  USf, USfC, VCha, VChaC, ZNeu, ZNeuC;
#endprocedure

#procedure Insertions
L iT1G1 = T1G1*(replace_(pM1, MA0, pM2, MZ, pM3, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (2*i_*Alfa*Pi)/(CW2*SW2), cR1, SUNT(Col3, Col4)) + 
   replace_(pM1, Mh0, pM2, MZ, pM3, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (2*i_*Alfa*Pi)/(CW2*SW2), cR1, SUNT(Col3, Col4)) + 
   replace_(pM1, MHH, pM2, MZ, pM3, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (2*i_*Alfa*Pi)/(CW2*SW2), cR1, SUNT(Col3, Col4)) + 
   replace_(pM1, MHp, pM2, 0, pM3, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 8*i_*Alfa*Pi, cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MHp, pM2, 0, pM3, MZ, pG1, 
    i_*EL, pG2, i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, 
    (2*i_*EL*SW)/(3*CW), pG5, (4*i_*Alfa*Pi*(CW2 - SW2))/(CW*SW), 
    cR1, 2*SUNT(Col3, Col4)) + replace_(pM1, MHp, pM2, MZ, pM3, 0, 
    pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    (4*i_*Alfa*Pi*(CW2 - SW2))/(CW*SW), cR1, 2*SUNT(Col3, Col4)) + 
   replace_(pM1, MHp, pM2, MZ, pM3, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (2*i_*Alfa*Pi*(CW2 - SW2)^2)/(CW2*SW2), cR1, 2*SUNT(Col3, Col4)) + 
   replace_(pM1, MW, pM2, 0, pM3, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 8*i_*Alfa*Pi, cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MW, pM2, 0, pM3, MZ, pG1, 
    i_*EL, pG2, i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, 
    (2*i_*EL*SW)/(3*CW), pG5, (4*i_*Alfa*Pi*(CW2 - SW2))/(CW*SW), 
    cR1, 2*SUNT(Col3, Col4)) + replace_(pM1, MW, pM2, MZ, pM3, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    (4*i_*Alfa*Pi*(CW2 - SW2))/(CW*SW), cR1, 2*SUNT(Col3, Col4)) + 
   replace_(pM1, MW, pM2, MZ, pM3, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (2*i_*Alfa*Pi*(CW2 - SW2)^2)/(CW2*SW2), cR1, 2*SUNT(Col3, Col4)) + 
   replace_(pM1, MZ, pM2, MZ, pM3, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (2*i_*Alfa*Pi)/(CW2*SW2), cR1, SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(1, 1, Gen5), pM2, MZ, pM3, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (2*i_*Alfa*Pi)/(CW2*SW2), cR1, 2*SumOver(Gen5, 3)*
     SUNT(Col3, Col4)) + replace_(pM1, MSf(Sfe5, 2, Gen5), pM2, 0, 
    pM3, 0, pG1, i_*EL, pG2, i_*EL, pG3, (-2*i_*EL)/3, pG4, 
    (-2*i_*EL)/3, pG5, 8*i_*Alfa*Pi, cR1, 2*SumOver(Gen5, 3)*
     SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + replace_(pM1, 
    MSf(Sfe5, 2, Gen5), pM2, 0, pM3, MZ, pG1, i_*EL, pG2, i_*EL, 
    pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), 
    pG5, (-4*i_*Alfa*Pi*((-1 + 2*SW2)*USf(Sfe5, 1, 2, Gen5)*
        USfC(Sfe5, 1, 2, Gen5) + 2*SW2*USf(Sfe5, 2, 2, Gen5)*
        USfC(Sfe5, 2, 2, Gen5)))/(CW*SW), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(Sfe5, 2, Gen5), pM2, MZ, pM3, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    (-4*i_*Alfa*Pi*((-1 + 2*SW2)*USf(Sfe5, 1, 2, Gen5)*
        USfC(Sfe5, 1, 2, Gen5) + 2*SW2*USf(Sfe5, 2, 2, Gen5)*
        USfC(Sfe5, 2, 2, Gen5)))/(CW*SW), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(Sfe5, 2, Gen5), pM2, MZ, pM3, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (2*i_*Alfa*Pi*((1 - 2*SW2)^2*USf(Sfe5, 1, 2, Gen5)*
        USfC(Sfe5, 1, 2, Gen5) + 4*SW2^2*USf(Sfe5, 2, 2, Gen5)*
        USfC(Sfe5, 2, 2, Gen5)))/(CW2*SW2), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(Sfe5, 3, Gen5), pM2, 0, pM3, 0, pG1, i_*EL, 
    pG2, i_*EL, pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    (32*i_*Alfa*Pi)/9, cR1, 6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MSf(Sfe5, 3, Gen5), pM2, 0, 
    pM3, MZ, pG1, i_*EL, pG2, i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/
     (6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (-8*i_*Alfa*Pi*((-3 + 4*SW2)*USf(Sfe5, 1, 3, Gen5)*
        USfC(Sfe5, 1, 3, Gen5) + 4*SW2*USf(Sfe5, 2, 3, Gen5)*
        USfC(Sfe5, 2, 3, Gen5)))/(9*CW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(Sfe5, 3, Gen5), pM2, MZ, pM3, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    (-8*i_*Alfa*Pi*((-3 + 4*SW2)*USf(Sfe5, 1, 3, Gen5)*
        USfC(Sfe5, 1, 3, Gen5) + 4*SW2*USf(Sfe5, 2, 3, Gen5)*
        USfC(Sfe5, 2, 3, Gen5)))/(9*CW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(Sfe5, 3, Gen5), pM2, MZ, pM3, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (2*i_*Alfa*Pi*((3 - 4*SW2)^2*USf(Sfe5, 1, 3, Gen5)*
        USfC(Sfe5, 1, 3, Gen5) + 16*SW2^2*USf(Sfe5, 2, 3, Gen5)*
        USfC(Sfe5, 2, 3, Gen5)))/(9*CW2*SW2), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(Sfe5, 4, Gen5), pM2, 0, pM3, 0, pG1, i_*EL, 
    pG2, i_*EL, pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    (8*i_*Alfa*Pi)/9, cR1, 6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MSf(Sfe5, 4, Gen5), pM2, 0, 
    pM3, MZ, pG1, i_*EL, pG2, i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/
     (6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (-4*i_*Alfa*Pi*((-3 + 2*SW2)*USf(Sfe5, 1, 4, Gen5)*
        USfC(Sfe5, 1, 4, Gen5) + 2*SW2*USf(Sfe5, 2, 4, Gen5)*
        USfC(Sfe5, 2, 4, Gen5)))/(9*CW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(Sfe5, 4, Gen5), pM2, MZ, pM3, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    (-4*i_*Alfa*Pi*((-3 + 2*SW2)*USf(Sfe5, 1, 4, Gen5)*
        USfC(Sfe5, 1, 4, Gen5) + 2*SW2*USf(Sfe5, 2, 4, Gen5)*
        USfC(Sfe5, 2, 4, Gen5)))/(9*CW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(Sfe5, 4, Gen5), pM2, MZ, pM3, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (2*i_*Alfa*Pi*((3 - 2*SW2)^2*USf(Sfe5, 1, 4, Gen5)*
        USfC(Sfe5, 1, 4, Gen5) + 4*SW2^2*USf(Sfe5, 2, 4, Gen5)*
        USfC(Sfe5, 2, 4, Gen5)))/(9*CW2*SW2), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)));

L iT1G2 = T1G2*(replace_(pM1, 0, pM2, 0, pG1, i_*EL, pG2, 
    i_*EL, pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    4*i_*Alfa*Pi, pG6, -8*i_*Alfa*Pi) + replace_(pM1, 0, pM2, MZ, 
    pG1, i_*EL, pG2, i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), 
    pG4, (2*i_*EL*SW)/(3*CW), pG5, (4*i_*Alfa*CW*Pi)/SW, pG6, 
    (-8*i_*Alfa*CW*Pi)/SW) + replace_(pM1, MZ, pM2, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, (4*i_*Alfa*CW*Pi)/SW, 
    pG6, (-8*i_*Alfa*CW*Pi)/SW) + replace_(pM1, MZ, pM2, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (4*i_*Alfa*CW2*Pi)/SW2, pG6, (-8*i_*Alfa*CW2*Pi)/SW2));

L iT2G1 = T2G1*(replace_(pM1, MCha(Cha5), pM2, MCha(Cha5), pM3, 
    MA0, pM4, 0, cG1, i_*EL, cG2, i_*EL, pG3, i_*EL, pG4, i_*EL, 
    pG5, -((EL*powM(2, -1/2)*(SB*UChaC(Cha5, 2)*VChaC(Cha5, 1) + 
        CB*UChaC(Cha5, 1)*VChaC(Cha5, 2)))/SW), pG6, 
    (EL*powM(2, -1/2)*(SB*UCha(Cha5, 2)*VCha(Cha5, 1) + 
       CB*UCha(Cha5, 1)*VCha(Cha5, 2)))/SW, pG7, (EL*MT)/(2*MW*SW*TB), 
    pG8, -(EL*MT)/(2*MW*SW*TB), cR1, 2*SumOver(Cha5, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MCha(Cha5), pM2, MCha(Cha5), 
    pM3, Mh0, pM4, 0, cG1, i_*EL, cG2, i_*EL, pG3, i_*EL, pG4, 
    i_*EL, pG5, (i_*EL*powM(2, -1/2)*
      (SA*UChaC(Cha5, 2)*VChaC(Cha5, 1) - CA*UChaC(Cha5, 1)*
        VChaC(Cha5, 2)))/SW, pG6, (i_*EL*powM(2, -1/2)*
      (SA*UCha(Cha5, 2)*VCha(Cha5, 1) - CA*UCha(Cha5, 1)*VCha(Cha5, 2)))/
     SW, pG7, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG8, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), cR1, 2*SumOver(Cha5, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MCha(Cha5), pM2, MCha(Cha5), 
    pM3, MHH, pM4, 0, cG1, i_*EL, cG2, i_*EL, pG3, i_*EL, pG4, 
    i_*EL, pG5, -((i_*EL*powM(2, -1/2)*
       (CA*UChaC(Cha5, 2)*VChaC(Cha5, 1) + SA*UChaC(Cha5, 1)*
         VChaC(Cha5, 2)))/SW), pG6, 
    -((i_*EL*powM(2, -1/2)*(CA*UCha(Cha5, 2)*VCha(Cha5, 1) + 
        SA*UCha(Cha5, 1)*VCha(Cha5, 2)))/SW), pG7, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG8, -(i_*EL*MT*SA)/(2*MW*SB*SW), 
    cR1, 2*SumOver(Cha5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MCha(Cha5), pM2, MCha(Cha5), pM3, MZ, pM4, 0, cG1, 
    i_*EL, cG2, i_*EL, pG3, i_*EL, pG4, i_*EL, pG5, 
    (EL*powM(2, -1/2)*(CB*UChaC(Cha5, 2)*VChaC(Cha5, 1) - 
       SB*UChaC(Cha5, 1)*VChaC(Cha5, 2)))/SW, pG6, 
    (EL*powM(2, -1/2)*(-(CB*UCha(Cha5, 2)*VCha(Cha5, 1)) + 
       SB*UCha(Cha5, 1)*VCha(Cha5, 2)))/SW, pG7, (EL*MT)/(2*MW*SW), pG8, 
    -(EL*MT)/(2*MW*SW), cR1, 2*SumOver(Cha5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MCha(Cha5), pM2, MCha(Cha6), pM3, MA0, pM4, MZ, cG1, 
    -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - UCha(Cha6, 1)*
         UChaC(Cha5, 1) - (UCha(Cha6, 2)*UChaC(Cha5, 2))/2))/(CW*SW)), 
    cG2, -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - 
        VCha(Cha5, 1)*VChaC(Cha6, 1) - (VCha(Cha5, 2)*VChaC(Cha6, 2))/
         2))/(CW*SW)), pG3, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, 
    -((i_*EL*SW)/CW), pG5, 
    -((EL*powM(2, -1/2)*(SB*UChaC(Cha6, 2)*VChaC(Cha5, 1) + 
        CB*UChaC(Cha6, 1)*VChaC(Cha5, 2)))/SW), pG6, 
    (EL*powM(2, -1/2)*(SB*UCha(Cha5, 2)*VCha(Cha6, 1) + 
       CB*UCha(Cha5, 1)*VCha(Cha6, 2)))/SW, pG7, (EL*MT)/(2*MW*SW*TB), 
    pG8, -(EL*MT)/(2*MW*SW*TB), cR1, 2*SumOver(Cha5, 2)*SumOver(Cha6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MCha(Cha5), pM2, MCha(Cha6), 
    pM3, Mh0, pM4, MZ, cG1, 
    -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - UCha(Cha6, 1)*
         UChaC(Cha5, 1) - (UCha(Cha6, 2)*UChaC(Cha5, 2))/2))/(CW*SW)), 
    cG2, -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - 
        VCha(Cha5, 1)*VChaC(Cha6, 1) - (VCha(Cha5, 2)*VChaC(Cha6, 2))/
         2))/(CW*SW)), pG3, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, 
    -((i_*EL*SW)/CW), pG5, (i_*EL*powM(2, -1/2)*
      (SA*UChaC(Cha6, 2)*VChaC(Cha5, 1) - CA*UChaC(Cha6, 1)*
        VChaC(Cha5, 2)))/SW, pG6, (i_*EL*powM(2, -1/2)*
      (SA*UCha(Cha5, 2)*VCha(Cha6, 1) - CA*UCha(Cha5, 1)*VCha(Cha6, 2)))/
     SW, pG7, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG8, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), cR1, 2*SumOver(Cha5, 2)*
     SumOver(Cha6, 2)*SUNT(Col3, Col4)) + replace_(pM1, MCha(Cha5), 
    pM2, MCha(Cha6), pM3, MHH, pM4, MZ, cG1, 
    -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - UCha(Cha6, 1)*
         UChaC(Cha5, 1) - (UCha(Cha6, 2)*UChaC(Cha5, 2))/2))/(CW*SW)), 
    cG2, -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - 
        VCha(Cha5, 1)*VChaC(Cha6, 1) - (VCha(Cha5, 2)*VChaC(Cha6, 2))/
         2))/(CW*SW)), pG3, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, 
    -((i_*EL*SW)/CW), pG5, 
    -((i_*EL*powM(2, -1/2)*(CA*UChaC(Cha6, 2)*VChaC(Cha5, 1) + 
        SA*UChaC(Cha6, 1)*VChaC(Cha5, 2)))/SW), pG6, 
    -((i_*EL*powM(2, -1/2)*(CA*UCha(Cha5, 2)*VCha(Cha6, 1) + 
        SA*UCha(Cha5, 1)*VCha(Cha6, 2)))/SW), pG7, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG8, -(i_*EL*MT*SA)/(2*MW*SB*SW), 
    cR1, 2*SumOver(Cha5, 2)*SumOver(Cha6, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MCha(Cha5), pM2, MCha(Cha6), pM3, MZ, pM4, MZ, cG1, 
    -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - UCha(Cha6, 1)*
         UChaC(Cha5, 1) - (UCha(Cha6, 2)*UChaC(Cha5, 2))/2))/(CW*SW)), 
    cG2, -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - 
        VCha(Cha5, 1)*VChaC(Cha6, 1) - (VCha(Cha5, 2)*VChaC(Cha6, 2))/
         2))/(CW*SW)), pG3, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, 
    -((i_*EL*SW)/CW), pG5, 
    (EL*powM(2, -1/2)*(CB*UChaC(Cha6, 2)*VChaC(Cha5, 1) - 
       SB*UChaC(Cha6, 1)*VChaC(Cha5, 2)))/SW, pG6, 
    (EL*powM(2, -1/2)*(-(CB*UCha(Cha5, 2)*VCha(Cha6, 1)) + 
       SB*UCha(Cha5, 1)*VCha(Cha6, 2)))/SW, pG7, (EL*MT)/(2*MW*SW), pG8, 
    -(EL*MT)/(2*MW*SW), cR1, 2*SumOver(Cha5, 2)*SumOver(Cha6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, Mf(3, Gen5), pM2, Mf(3, Gen5), 
    pM3, MA0, pM4, 0, cG1, (-2*i_*EL)/3, cG2, (-2*i_*EL)/3, pG3, 
    i_*EL, pG4, i_*EL, pG5, (EL*Mf(3, Gen5))/(2*MW*SW*TB), pG6, 
    -(EL*Mf(3, Gen5))/(2*MW*SW*TB), pG7, (EL*MT)/(2*MW*SW*TB), pG8, 
    -(EL*MT)/(2*MW*SW*TB), cR1, 6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(3, Gen5), pM2, Mf(3, Gen5), pM3, MA0, pM4, MZ, 
    cG1, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), cG2, (2*i_*EL*SW)/(3*CW), 
    pG3, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, -((i_*EL*SW)/CW), 
    pG5, (EL*Mf(3, Gen5))/(2*MW*SW*TB), pG6, -(EL*Mf(3, Gen5))/
     (2*MW*SW*TB), pG7, (EL*MT)/(2*MW*SW*TB), pG8, 
    -(EL*MT)/(2*MW*SW*TB), cR1, 6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(3, Gen5), pM2, Mf(3, Gen5), pM3, Mh0, pM4, 0, cG1, 
    (-2*i_*EL)/3, cG2, (-2*i_*EL)/3, pG3, i_*EL, pG4, i_*EL, 
    pG5, -(i_*CA*EL*Mf(3, Gen5))/(2*MW*SB*SW), pG6, 
    -(i_*CA*EL*Mf(3, Gen5))/(2*MW*SB*SW), pG7, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG8, -(i_*CA*EL*MT)/(2*MW*SB*SW), 
    cR1, 6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(3, Gen5), pM2, Mf(3, Gen5), pM3, Mh0, pM4, MZ, 
    cG1, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), cG2, (2*i_*EL*SW)/(3*CW), 
    pG3, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, -((i_*EL*SW)/CW), 
    pG5, -(i_*CA*EL*Mf(3, Gen5))/(2*MW*SB*SW), pG6, 
    -(i_*CA*EL*Mf(3, Gen5))/(2*MW*SB*SW), pG7, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG8, -(i_*CA*EL*MT)/(2*MW*SB*SW), 
    cR1, 6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(3, Gen5), pM2, Mf(3, Gen5), pM3, MHH, pM4, 0, cG1, 
    (-2*i_*EL)/3, cG2, (-2*i_*EL)/3, pG3, i_*EL, pG4, i_*EL, 
    pG5, -(i_*EL*SA*Mf(3, Gen5))/(2*MW*SB*SW), pG6, 
    -(i_*EL*SA*Mf(3, Gen5))/(2*MW*SB*SW), pG7, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG8, -(i_*EL*MT*SA)/(2*MW*SB*SW), 
    cR1, 6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(3, Gen5), pM2, Mf(3, Gen5), pM3, MHH, pM4, MZ, 
    cG1, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), cG2, (2*i_*EL*SW)/(3*CW), 
    pG3, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, -((i_*EL*SW)/CW), 
    pG5, -(i_*EL*SA*Mf(3, Gen5))/(2*MW*SB*SW), pG6, 
    -(i_*EL*SA*Mf(3, Gen5))/(2*MW*SB*SW), pG7, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG8, -(i_*EL*MT*SA)/(2*MW*SB*SW), 
    cR1, 6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(3, Gen5), pM2, Mf(3, Gen5), pM3, MZ, pM4, 0, cG1, 
    (-2*i_*EL)/3, cG2, (-2*i_*EL)/3, pG3, i_*EL, pG4, i_*EL, 
    pG5, (EL*Mf(3, Gen5))/(2*MW*SW), pG6, -(EL*Mf(3, Gen5))/(2*MW*SW), 
    pG7, (EL*MT)/(2*MW*SW), pG8, -(EL*MT)/(2*MW*SW), cR1, 
    6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + replace_(pM1, Mf(3, Gen5), 
    pM2, Mf(3, Gen5), pM3, MZ, pM4, MZ, cG1, (i_*EL*(-3 + 4*SW2))/
     (6*CW*SW), cG2, (2*i_*EL*SW)/(3*CW), pG3, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, -((i_*EL*SW)/CW), pG5, 
    (EL*Mf(3, Gen5))/(2*MW*SW), pG6, -(EL*Mf(3, Gen5))/(2*MW*SW), pG7, 
    (EL*MT)/(2*MW*SW), pG8, -(EL*MT)/(2*MW*SW), cR1, 
    6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + replace_(pM1, MNeu(Neu5), 
    pM2, MNeu(Neu6), pM3, MA0, pM4, MZ, cG1, 
    (i_*EL*(-(ZNeu(Neu6, 3)*ZNeuC(Neu5, 3)) + ZNeu(Neu6, 4)*
        ZNeuC(Neu5, 4)))/(2*CW*SW), cG2, 
    (i_*EL*(ZNeu(Neu5, 3)*ZNeuC(Neu6, 3) - ZNeu(Neu5, 4)*
        ZNeuC(Neu6, 4)))/(2*CW*SW), pG3, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, -((i_*EL*SW)/CW), pG5, 
    (EL*((SB*ZNeuC(Neu5, 3) - CB*ZNeuC(Neu5, 4))*(SW*ZNeuC(Neu6, 1) - 
         CW*ZNeuC(Neu6, 2)) + SB*(SW*ZNeuC(Neu5, 1) - CW*ZNeuC(Neu5, 2))*
        ZNeuC(Neu6, 3) + (-(CB*SW*ZNeuC(Neu5, 1)) + 
         CB*CW*ZNeuC(Neu5, 2))*ZNeuC(Neu6, 4)))/(2*CW*SW), pG6, 
    -(EL*((SB*ZNeu(Neu5, 3) - CB*ZNeu(Neu5, 4))*(SW*ZNeu(Neu6, 1) - 
          CW*ZNeu(Neu6, 2)) + SB*(SW*ZNeu(Neu5, 1) - CW*ZNeu(Neu5, 2))*
         ZNeu(Neu6, 3) + (-(CB*SW*ZNeu(Neu5, 1)) + CB*CW*ZNeu(Neu5, 2))*
         ZNeu(Neu6, 4)))/(2*CW*SW), pG7, (EL*MT)/(2*MW*SW*TB), pG8, 
    -(EL*MT)/(2*MW*SW*TB), cR1, SumOver(Neu5, 4)*SumOver(Neu6, 4)*
     SUNT(Col3, Col4)) + replace_(pM1, MNeu(Neu5), pM2, MNeu(Neu6), 
    pM3, Mh0, pM4, MZ, cG1, (i_*EL*(-(ZNeu(Neu6, 3)*ZNeuC(Neu5, 3)) + 
       ZNeu(Neu6, 4)*ZNeuC(Neu5, 4)))/(2*CW*SW), cG2, 
    (i_*EL*(ZNeu(Neu5, 3)*ZNeuC(Neu6, 3) - ZNeu(Neu5, 4)*
        ZNeuC(Neu6, 4)))/(2*CW*SW), pG3, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, -((i_*EL*SW)/CW), pG5, 
    -(i_*EL*((SA*ZNeuC(Neu5, 3) + CA*ZNeuC(Neu5, 4))*
         (SW*ZNeuC(Neu6, 1) - CW*ZNeuC(Neu6, 2)) + 
        SA*(SW*ZNeuC(Neu5, 1) - CW*ZNeuC(Neu5, 2))*ZNeuC(Neu6, 3) + 
        CA*(SW*ZNeuC(Neu5, 1) - CW*ZNeuC(Neu5, 2))*ZNeuC(Neu6, 4)))/
     (2*CW*SW), pG6, 
    (i_*EL*(-((SA*ZNeu(Neu5, 3) + CA*ZNeu(Neu5, 4))*
         (SW*ZNeu(Neu6, 1) - CW*ZNeu(Neu6, 2))) + 
       (-(SA*SW*ZNeu(Neu5, 1)) + CW*SA*ZNeu(Neu5, 2))*ZNeu(Neu6, 3) + 
       (-(CA*SW*ZNeu(Neu5, 1)) + CA*CW*ZNeu(Neu5, 2))*ZNeu(Neu6, 4)))/
     (2*CW*SW), pG7, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG8, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), cR1, SumOver(Neu5, 4)*
     SumOver(Neu6, 4)*SUNT(Col3, Col4)) + replace_(pM1, MNeu(Neu5), 
    pM2, MNeu(Neu6), pM3, MHH, pM4, MZ, cG1, 
    (i_*EL*(-(ZNeu(Neu6, 3)*ZNeuC(Neu5, 3)) + ZNeu(Neu6, 4)*
        ZNeuC(Neu5, 4)))/(2*CW*SW), cG2, 
    (i_*EL*(ZNeu(Neu5, 3)*ZNeuC(Neu6, 3) - ZNeu(Neu5, 4)*
        ZNeuC(Neu6, 4)))/(2*CW*SW), pG3, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, -((i_*EL*SW)/CW), pG5, 
    (i_*EL*((CA*ZNeuC(Neu5, 3) - SA*ZNeuC(Neu5, 4))*
        (SW*ZNeuC(Neu6, 1) - CW*ZNeuC(Neu6, 2)) + 
       CA*(SW*ZNeuC(Neu5, 1) - CW*ZNeuC(Neu5, 2))*ZNeuC(Neu6, 3) + 
       (-(SA*SW*ZNeuC(Neu5, 1)) + CW*SA*ZNeuC(Neu5, 2))*ZNeuC(Neu6, 4)))/
     (2*CW*SW), pG6, (i_*EL*((CA*ZNeu(Neu5, 3) - SA*ZNeu(Neu5, 4))*
        (SW*ZNeu(Neu6, 1) - CW*ZNeu(Neu6, 2)) + 
       CA*(SW*ZNeu(Neu5, 1) - CW*ZNeu(Neu5, 2))*ZNeu(Neu6, 3) + 
       (-(SA*SW*ZNeu(Neu5, 1)) + CW*SA*ZNeu(Neu5, 2))*ZNeu(Neu6, 4)))/
     (2*CW*SW), pG7, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG8, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), cR1, SumOver(Neu5, 4)*
     SumOver(Neu6, 4)*SUNT(Col3, Col4)) + replace_(pM1, MNeu(Neu5), 
    pM2, MNeu(Neu6), pM3, MZ, pM4, MZ, cG1, 
    (i_*EL*(-(ZNeu(Neu6, 3)*ZNeuC(Neu5, 3)) + ZNeu(Neu6, 4)*
        ZNeuC(Neu5, 4)))/(2*CW*SW), cG2, 
    (i_*EL*(ZNeu(Neu5, 3)*ZNeuC(Neu6, 3) - ZNeu(Neu5, 4)*
        ZNeuC(Neu6, 4)))/(2*CW*SW), pG3, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG4, -((i_*EL*SW)/CW), pG5, 
    -(EL*((CB*ZNeuC(Neu5, 3) + SB*ZNeuC(Neu5, 4))*(SW*ZNeuC(Neu6, 1) - 
          CW*ZNeuC(Neu6, 2)) + CB*(SW*ZNeuC(Neu5, 1) - 
          CW*ZNeuC(Neu5, 2))*ZNeuC(Neu6, 3) + 
        SB*(SW*ZNeuC(Neu5, 1) - CW*ZNeuC(Neu5, 2))*ZNeuC(Neu6, 4)))/
     (2*CW*SW), pG6, (EL*((CB*ZNeu(Neu5, 3) + SB*ZNeu(Neu5, 4))*
        (SW*ZNeu(Neu6, 1) - CW*ZNeu(Neu6, 2)) + 
       CB*(SW*ZNeu(Neu5, 1) - CW*ZNeu(Neu5, 2))*ZNeu(Neu6, 3) + 
       SB*(SW*ZNeu(Neu5, 1) - CW*ZNeu(Neu5, 2))*ZNeu(Neu6, 4)))/
     (2*CW*SW), pG7, (EL*MT)/(2*MW*SW), pG8, -(EL*MT)/(2*MW*SW), cR1, 
    SumOver(Neu5, 4)*SumOver(Neu6, 4)*SUNT(Col3, Col4)));

L iT2G2 = T2G2*(replace_(pM1, MA0, pM2, Mh0, pM3, MA0, pM4, MZ, 
    pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (CBA*EL)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW*TB), pG5, 
    -(EL*MT)/(2*MW*SW*TB), pG6, -(i_*C2B*EL*MW*SAB)/(2*CW2*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MA0, pM2, Mh0, pM3, MZ, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (EL*SBA)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW*TB), pG5, 
    -(EL*MT)/(2*MW*SW*TB), pG6, -(i_*EL*MW*S2B*SAB)/(2*CW2*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MA0, pM2, MHH, pM3, MA0, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, -(EL*SBA)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW*TB), pG5, 
    -(EL*MT)/(2*MW*SW*TB), pG6, (i_*C2B*CAB*EL*MW)/(2*CW2*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MA0, pM2, MHH, pM3, MZ, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (CBA*EL)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW*TB), pG5, 
    -(EL*MT)/(2*MW*SW*TB), pG6, (i_*CAB*EL*MW*S2B)/(2*CW2*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MA0, pM2, MSf(Sfe5, 2, Gen5), 
    pM3, MSf(Sfe5, 2, Gen5), pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    i_*EL, pG4, (EL*MT)/(2*MW*SW*TB), pG5, -(EL*MT)/(2*MW*SW*TB), pG6, 
    -(EL*Mf(2, Gen5)*(-((MUEC + TB*Af(2, Gen5, Gen5))*
          USf(Sfe5, 2, 2, Gen5)*USfC(Sfe5, 1, 2, Gen5)) + 
        (MUE + TB*AfC(2, Gen5, Gen5))*USf(Sfe5, 1, 2, Gen5)*
         USfC(Sfe5, 2, 2, Gen5)))/(2*MW*SW), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MA0, pM2, MSf(Sfe5, 2, Gen5), pM3, 
    MSf(Sfe6, 2, Gen5), pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, 
    -(i_*EL*((-1 + 2*SW2)*USf(Sfe6, 1, 2, Gen5)*USfC(Sfe5, 1, 2, 
          Gen5) + 2*SW2*USf(Sfe6, 2, 2, Gen5)*USfC(Sfe5, 2, 2, Gen5)))/
     (2*CW*SW), pG4, (EL*MT)/(2*MW*SW*TB), pG5, -(EL*MT)/(2*MW*SW*TB), 
    pG6, -(EL*Mf(2, Gen5)*(-((MUEC + TB*Af(2, Gen5, Gen5))*
          USf(Sfe5, 2, 2, Gen5)*USfC(Sfe6, 1, 2, Gen5)) + 
        (MUE + TB*AfC(2, Gen5, Gen5))*USf(Sfe5, 1, 2, Gen5)*
         USfC(Sfe6, 2, 2, Gen5)))/(2*MW*SW), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MA0, pM2, MSf(Sfe5, 3, Gen5), 
    pM3, MSf(Sfe5, 3, Gen5), pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    (-2*i_*EL)/3, pG4, (EL*MT)/(2*MW*SW*TB), pG5, 
    -(EL*MT)/(2*MW*SW*TB), pG6, 
    -(EL*Mf(3, Gen5)*(-((MUEC*TB + Af(3, Gen5, Gen5))*
          USf(Sfe5, 2, 3, Gen5)*USfC(Sfe5, 1, 3, Gen5)) + 
        (MUE*TB + AfC(3, Gen5, Gen5))*USf(Sfe5, 1, 3, Gen5)*
         USfC(Sfe5, 2, 3, Gen5)))/(2*MW*SW*TB), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MA0, pM2, MSf(Sfe5, 3, Gen5), pM3, 
    MSf(Sfe6, 3, Gen5), pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*((-3 + 4*SW2)*USf(Sfe6, 1, 3, Gen5)*USfC(Sfe5, 1, 3, 
         Gen5) + 4*SW2*USf(Sfe6, 2, 3, Gen5)*USfC(Sfe5, 2, 3, Gen5)))/
     (6*CW*SW), pG4, (EL*MT)/(2*MW*SW*TB), pG5, -(EL*MT)/(2*MW*SW*TB), 
    pG6, -(EL*Mf(3, Gen5)*(-((MUEC*TB + Af(3, Gen5, Gen5))*
          USf(Sfe5, 2, 3, Gen5)*USfC(Sfe6, 1, 3, Gen5)) + 
        (MUE*TB + AfC(3, Gen5, Gen5))*USf(Sfe5, 1, 3, Gen5)*
         USfC(Sfe6, 2, 3, Gen5)))/(2*MW*SW*TB), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MA0, pM2, MSf(Sfe5, 4, Gen5), 
    pM3, MSf(Sfe5, 4, Gen5), pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    (i_*EL)/3, pG4, (EL*MT)/(2*MW*SW*TB), pG5, -(EL*MT)/(2*MW*SW*TB), 
    pG6, -(EL*Mf(4, Gen5)*(-((MUEC + TB*Af(4, Gen5, Gen5))*
          USf(Sfe5, 2, 4, Gen5)*USfC(Sfe5, 1, 4, Gen5)) + 
        (MUE + TB*AfC(4, Gen5, Gen5))*USf(Sfe5, 1, 4, Gen5)*
         USfC(Sfe5, 2, 4, Gen5)))/(2*MW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MA0, pM2, MSf(Sfe5, 4, Gen5), pM3, 
    MSf(Sfe6, 4, Gen5), pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, 
    -(i_*EL*((-3 + 2*SW2)*USf(Sfe6, 1, 4, Gen5)*USfC(Sfe5, 1, 4, 
          Gen5) + 2*SW2*USf(Sfe6, 2, 4, Gen5)*USfC(Sfe5, 2, 4, Gen5)))/
     (6*CW*SW), pG4, (EL*MT)/(2*MW*SW*TB), pG5, -(EL*MT)/(2*MW*SW*TB), 
    pG6, -(EL*Mf(4, Gen5)*(-((MUEC + TB*Af(4, Gen5, Gen5))*
          USf(Sfe5, 2, 4, Gen5)*USfC(Sfe6, 1, 4, Gen5)) + 
        (MUE + TB*AfC(4, Gen5, Gen5))*USf(Sfe5, 1, 4, Gen5)*
         USfC(Sfe6, 2, 4, Gen5)))/(2*MW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, Mh0, pM2, MHp, pM3, MHp, pM4, 
    0, pG1, i_*EL, pG2, i_*EL, pG3, i_*EL, pG4, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, -(i_*CA*EL*MT)/(2*MW*SB*SW), 
    pG6, -((i_*EL*MW*((C2B*SAB)/(2*CW2) + SBA))/SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, Mh0, pM2, MHp, pM3, MHp, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (i_*EL*(CW2 - SW2))/(2*CW*SW), pG4, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, -(i_*CA*EL*MT)/(2*MW*SB*SW), 
    pG6, -((i_*EL*MW*((C2B*SAB)/(2*CW2) + SBA))/SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, Mh0, pM2, MW, pM3, MW, pM4, 0, 
    pG1, i_*EL, pG2, i_*EL, pG3, i_*EL, pG4, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, -(i_*CA*EL*MT)/(2*MW*SB*SW), 
    pG6, (i_*C2B*EL*MW*SAB)/(2*CW2*SW), cR1, 2*SUNT(Col3, Col4)) + 
   replace_(pM1, Mh0, pM2, MW, pM3, MW, pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(CW2 - SW2))/(2*CW*SW), pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), 
    pG5, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG6, (i_*C2B*EL*MW*SAB)/
     (2*CW2*SW), cR1, 2*SUNT(Col3, Col4)) + replace_(pM1, Mh0, pM2, 
    MSf(1, 1, Gen5), pM3, MSf(1, 1, Gen5), pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    -(i_*EL)/(2*CW*SW), pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG6, (i_*EL*MZ*SAB)/(2*CW*SW), cR1, 
    2*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + replace_(pM1, Mh0, pM2, 
    MSf(Sfe5, 2, Gen5), pM3, MSf(Sfe5, 2, Gen5), pM4, 0, pG1, i_*EL, 
    pG2, i_*EL, pG3, i_*EL, pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG6, 
    (i_*EL*(((CB*MW*MZ*SAB*(-1 + 2*SW2) + 2*CW*SA*Mf2(2, Gen5))*
          USf(Sfe5, 1, 2, Gen5) + CW*(CA*MUEC + SA*Af(2, Gen5, Gen5))*
          Mf(2, Gen5)*USf(Sfe5, 2, 2, Gen5))*USfC(Sfe5, 1, 2, Gen5) + 
       (CW*(CA*MUE + SA*AfC(2, Gen5, Gen5))*Mf(2, Gen5)*
          USf(Sfe5, 1, 2, Gen5) - 2*CB*MW*MZ*SAB*SW2*USf(Sfe5, 2, 2, 
           Gen5) + 2*CW*SA*Mf2(2, Gen5)*USf(Sfe5, 2, 2, Gen5))*
        USfC(Sfe5, 2, 2, Gen5)))/(2*CB*CW*MW*SW), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mh0, pM2, MSf(Sfe5, 2, Gen5), pM3, 
    MSf(Sfe6, 2, Gen5), pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, 
    -(i_*EL*((-1 + 2*SW2)*USf(Sfe6, 1, 2, Gen5)*USfC(Sfe5, 1, 2, 
          Gen5) + 2*SW2*USf(Sfe6, 2, 2, Gen5)*USfC(Sfe5, 2, 2, Gen5)))/
     (2*CW*SW), pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG6, 
    (i_*EL*(((CB*MW*MZ*SAB*(-1 + 2*SW2) + 2*CW*SA*Mf2(2, Gen5))*
          USf(Sfe5, 1, 2, Gen5) + CW*(CA*MUEC + SA*Af(2, Gen5, Gen5))*
          Mf(2, Gen5)*USf(Sfe5, 2, 2, Gen5))*USfC(Sfe6, 1, 2, Gen5) + 
       (CW*(CA*MUE + SA*AfC(2, Gen5, Gen5))*Mf(2, Gen5)*
          USf(Sfe5, 1, 2, Gen5) - 2*CB*MW*MZ*SAB*SW2*USf(Sfe5, 2, 2, 
           Gen5) + 2*CW*SA*Mf2(2, Gen5)*USf(Sfe5, 2, 2, Gen5))*
        USfC(Sfe6, 2, 2, Gen5)))/(2*CB*CW*MW*SW), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, Mh0, pM2, MSf(Sfe5, 3, Gen5), 
    pM3, MSf(Sfe5, 3, Gen5), pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    (-2*i_*EL)/3, pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG6, 
    -(i_*EL*(((MW*MZ*SAB*SB*(-3 + 4*SW2) + 6*CA*CW*Mf2(3, Gen5))*
           USf(Sfe5, 1, 3, Gen5) + 3*CW*(MUEC*SA + CA*Af(3, Gen5, Gen5))*
           Mf(3, Gen5)*USf(Sfe5, 2, 3, Gen5))*USfC(Sfe5, 1, 3, Gen5) + 
        (3*CW*(MUE*SA + CA*AfC(3, Gen5, Gen5))*Mf(3, Gen5)*
           USf(Sfe5, 1, 3, Gen5) - 4*MW*MZ*SAB*SB*SW2*USf(Sfe5, 2, 3, 
            Gen5) + 6*CA*CW*Mf2(3, Gen5)*USf(Sfe5, 2, 3, Gen5))*
         USfC(Sfe5, 2, 3, Gen5)))/(6*CW*MW*SB*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mh0, pM2, MSf(Sfe5, 3, Gen5), pM3, 
    MSf(Sfe6, 3, Gen5), pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*((-3 + 4*SW2)*USf(Sfe6, 1, 3, Gen5)*USfC(Sfe5, 1, 3, 
         Gen5) + 4*SW2*USf(Sfe6, 2, 3, Gen5)*USfC(Sfe5, 2, 3, Gen5)))/
     (6*CW*SW), pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG6, 
    -(i_*EL*(((MW*MZ*SAB*SB*(-3 + 4*SW2) + 6*CA*CW*Mf2(3, Gen5))*
           USf(Sfe5, 1, 3, Gen5) + 3*CW*(MUEC*SA + CA*Af(3, Gen5, Gen5))*
           Mf(3, Gen5)*USf(Sfe5, 2, 3, Gen5))*USfC(Sfe6, 1, 3, Gen5) + 
        (3*CW*(MUE*SA + CA*AfC(3, Gen5, Gen5))*Mf(3, Gen5)*
           USf(Sfe5, 1, 3, Gen5) - 4*MW*MZ*SAB*SB*SW2*USf(Sfe5, 2, 3, 
            Gen5) + 6*CA*CW*Mf2(3, Gen5)*USf(Sfe5, 2, 3, Gen5))*
         USfC(Sfe6, 2, 3, Gen5)))/(6*CW*MW*SB*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, Mh0, pM2, MSf(Sfe5, 4, Gen5), 
    pM3, MSf(Sfe5, 4, Gen5), pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    (i_*EL)/3, pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG6, 
    (i_*EL*(((CB*MW*MZ*SAB*(-3 + 2*SW2) + 6*CW*SA*Mf2(4, Gen5))*
          USf(Sfe5, 1, 4, Gen5) + 3*CW*(CA*MUEC + SA*Af(4, Gen5, Gen5))*
          Mf(4, Gen5)*USf(Sfe5, 2, 4, Gen5))*USfC(Sfe5, 1, 4, Gen5) + 
       (3*CW*(CA*MUE + SA*AfC(4, Gen5, Gen5))*Mf(4, Gen5)*
          USf(Sfe5, 1, 4, Gen5) - 2*CB*MW*MZ*SAB*SW2*USf(Sfe5, 2, 4, 
           Gen5) + 6*CW*SA*Mf2(4, Gen5)*USf(Sfe5, 2, 4, Gen5))*
        USfC(Sfe5, 2, 4, Gen5)))/(6*CB*CW*MW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mh0, pM2, MSf(Sfe5, 4, Gen5), pM3, 
    MSf(Sfe6, 4, Gen5), pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, 
    -(i_*EL*((-3 + 2*SW2)*USf(Sfe6, 1, 4, Gen5)*USfC(Sfe5, 1, 4, 
          Gen5) + 2*SW2*USf(Sfe6, 2, 4, Gen5)*USfC(Sfe5, 2, 4, Gen5)))/
     (6*CW*SW), pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG6, 
    (i_*EL*(((CB*MW*MZ*SAB*(-3 + 2*SW2) + 6*CW*SA*Mf2(4, Gen5))*
          USf(Sfe5, 1, 4, Gen5) + 3*CW*(CA*MUEC + SA*Af(4, Gen5, Gen5))*
          Mf(4, Gen5)*USf(Sfe5, 2, 4, Gen5))*USfC(Sfe6, 1, 4, Gen5) + 
       (3*CW*(CA*MUE + SA*AfC(4, Gen5, Gen5))*Mf(4, Gen5)*
          USf(Sfe5, 1, 4, Gen5) - 2*CB*MW*MZ*SAB*SW2*USf(Sfe5, 2, 4, 
           Gen5) + 6*CW*SA*Mf2(4, Gen5)*USf(Sfe5, 2, 4, Gen5))*
        USfC(Sfe6, 2, 4, Gen5)))/(6*CB*CW*MW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MHH, pM2, MHp, pM3, MHp, pM4, 
    0, pG1, i_*EL, pG2, i_*EL, pG3, i_*EL, pG4, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, -(i_*EL*MT*SA)/(2*MW*SB*SW), 
    pG6, -((i_*(CBA - (C2B*CAB)/(2*CW2))*EL*MW)/SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MHH, pM2, MHp, pM3, MHp, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (i_*EL*(CW2 - SW2))/(2*CW*SW), pG4, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, -(i_*EL*MT*SA)/(2*MW*SB*SW), 
    pG6, -((i_*(CBA - (C2B*CAB)/(2*CW2))*EL*MW)/SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MHH, pM2, MW, pM3, MW, pM4, 0, 
    pG1, i_*EL, pG2, i_*EL, pG3, i_*EL, pG4, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, -(i_*EL*MT*SA)/(2*MW*SB*SW), 
    pG6, -(i_*C2B*CAB*EL*MW)/(2*CW2*SW), cR1, 2*SUNT(Col3, Col4)) + 
   replace_(pM1, MHH, pM2, MW, pM3, MW, pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(CW2 - SW2))/(2*CW*SW), pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), 
    pG5, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG6, -(i_*C2B*CAB*EL*MW)/
     (2*CW2*SW), cR1, 2*SUNT(Col3, Col4)) + replace_(pM1, MHH, pM2, 
    MSf(1, 1, Gen5), pM3, MSf(1, 1, Gen5), pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    -(i_*EL)/(2*CW*SW), pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG6, -(i_*CAB*EL*MZ)/(2*CW*SW), 
    cR1, 2*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, MHH, pM2, MSf(Sfe5, 2, Gen5), pM3, 
    MSf(Sfe5, 2, Gen5), pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    i_*EL, pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG6, 
    (i_*EL*(((CAB*CB*MW*MZ*(1 - 2*SW2) - 2*CA*CW*Mf2(2, Gen5))*
          USf(Sfe5, 1, 2, Gen5) + CW*(MUEC*SA - CA*Af(2, Gen5, Gen5))*
          Mf(2, Gen5)*USf(Sfe5, 2, 2, Gen5))*USfC(Sfe5, 1, 2, Gen5) + 
       (CW*(MUE*SA - CA*AfC(2, Gen5, Gen5))*Mf(2, Gen5)*
          USf(Sfe5, 1, 2, Gen5) + 2*CAB*CB*MW*MZ*SW2*USf(Sfe5, 2, 2, 
           Gen5) - 2*CA*CW*Mf2(2, Gen5)*USf(Sfe5, 2, 2, Gen5))*
        USfC(Sfe5, 2, 2, Gen5)))/(2*CB*CW*MW*SW), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MHH, pM2, MSf(Sfe5, 2, Gen5), pM3, 
    MSf(Sfe6, 2, Gen5), pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, 
    -(i_*EL*((-1 + 2*SW2)*USf(Sfe6, 1, 2, Gen5)*USfC(Sfe5, 1, 2, 
          Gen5) + 2*SW2*USf(Sfe6, 2, 2, Gen5)*USfC(Sfe5, 2, 2, Gen5)))/
     (2*CW*SW), pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG6, 
    (i_*EL*(((CAB*CB*MW*MZ*(1 - 2*SW2) - 2*CA*CW*Mf2(2, Gen5))*
          USf(Sfe5, 1, 2, Gen5) + CW*(MUEC*SA - CA*Af(2, Gen5, Gen5))*
          Mf(2, Gen5)*USf(Sfe5, 2, 2, Gen5))*USfC(Sfe6, 1, 2, Gen5) + 
       (CW*(MUE*SA - CA*AfC(2, Gen5, Gen5))*Mf(2, Gen5)*
          USf(Sfe5, 1, 2, Gen5) + 2*CAB*CB*MW*MZ*SW2*USf(Sfe5, 2, 2, 
           Gen5) - 2*CA*CW*Mf2(2, Gen5)*USf(Sfe5, 2, 2, Gen5))*
        USfC(Sfe6, 2, 2, Gen5)))/(2*CB*CW*MW*SW), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MHH, pM2, MSf(Sfe5, 3, Gen5), 
    pM3, MSf(Sfe5, 3, Gen5), pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    (-2*i_*EL)/3, pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG6, 
    -(i_*EL*(((CAB*MW*MZ*SB*(3 - 4*SW2) + 6*CW*SA*Mf2(3, Gen5))*
           USf(Sfe5, 1, 3, Gen5) + 3*CW*(-(CA*MUEC) + 
            SA*Af(3, Gen5, Gen5))*Mf(3, Gen5)*USf(Sfe5, 2, 3, Gen5))*
         USfC(Sfe5, 1, 3, Gen5) + 
        (3*CW*(-(CA*MUE) + SA*AfC(3, Gen5, Gen5))*Mf(3, Gen5)*
           USf(Sfe5, 1, 3, Gen5) + 4*CAB*MW*MZ*SB*SW2*USf(Sfe5, 2, 3, 
            Gen5) + 6*CW*SA*Mf2(3, Gen5)*USf(Sfe5, 2, 3, Gen5))*
         USfC(Sfe5, 2, 3, Gen5)))/(6*CW*MW*SB*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MHH, pM2, MSf(Sfe5, 3, Gen5), pM3, 
    MSf(Sfe6, 3, Gen5), pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*((-3 + 4*SW2)*USf(Sfe6, 1, 3, Gen5)*USfC(Sfe5, 1, 3, 
         Gen5) + 4*SW2*USf(Sfe6, 2, 3, Gen5)*USfC(Sfe5, 2, 3, Gen5)))/
     (6*CW*SW), pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG6, 
    -(i_*EL*(((CAB*MW*MZ*SB*(3 - 4*SW2) + 6*CW*SA*Mf2(3, Gen5))*
           USf(Sfe5, 1, 3, Gen5) + 3*CW*(-(CA*MUEC) + 
            SA*Af(3, Gen5, Gen5))*Mf(3, Gen5)*USf(Sfe5, 2, 3, Gen5))*
         USfC(Sfe6, 1, 3, Gen5) + 
        (3*CW*(-(CA*MUE) + SA*AfC(3, Gen5, Gen5))*Mf(3, Gen5)*
           USf(Sfe5, 1, 3, Gen5) + 4*CAB*MW*MZ*SB*SW2*USf(Sfe5, 2, 3, 
            Gen5) + 6*CW*SA*Mf2(3, Gen5)*USf(Sfe5, 2, 3, Gen5))*
         USfC(Sfe6, 2, 3, Gen5)))/(6*CW*MW*SB*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MHH, pM2, MSf(Sfe5, 4, Gen5), 
    pM3, MSf(Sfe5, 4, Gen5), pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    (i_*EL)/3, pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG6, 
    -(i_*EL*(((CAB*CB*MW*MZ*(-3 + 2*SW2) + 6*CA*CW*Mf2(4, Gen5))*
           USf(Sfe5, 1, 4, Gen5) + 3*CW*(-(MUEC*SA) + 
            CA*Af(4, Gen5, Gen5))*Mf(4, Gen5)*USf(Sfe5, 2, 4, Gen5))*
         USfC(Sfe5, 1, 4, Gen5) + 
        (3*CW*(-(MUE*SA) + CA*AfC(4, Gen5, Gen5))*Mf(4, Gen5)*
           USf(Sfe5, 1, 4, Gen5) - 2*CAB*CB*MW*MZ*SW2*USf(Sfe5, 2, 4, 
            Gen5) + 6*CA*CW*Mf2(4, Gen5)*USf(Sfe5, 2, 4, Gen5))*
         USfC(Sfe5, 2, 4, Gen5)))/(6*CB*CW*MW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MHH, pM2, MSf(Sfe5, 4, Gen5), pM3, 
    MSf(Sfe6, 4, Gen5), pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, 
    -(i_*EL*((-3 + 2*SW2)*USf(Sfe6, 1, 4, Gen5)*USfC(Sfe5, 1, 4, 
          Gen5) + 2*SW2*USf(Sfe6, 2, 4, Gen5)*USfC(Sfe5, 2, 4, Gen5)))/
     (6*CW*SW), pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG6, 
    -(i_*EL*(((CAB*CB*MW*MZ*(-3 + 2*SW2) + 6*CA*CW*Mf2(4, Gen5))*
           USf(Sfe5, 1, 4, Gen5) + 3*CW*(-(MUEC*SA) + 
            CA*Af(4, Gen5, Gen5))*Mf(4, Gen5)*USf(Sfe5, 2, 4, Gen5))*
         USfC(Sfe6, 1, 4, Gen5) + 
        (3*CW*(-(MUE*SA) + CA*AfC(4, Gen5, Gen5))*Mf(4, Gen5)*
           USf(Sfe5, 1, 4, Gen5) - 2*CAB*CB*MW*MZ*SW2*USf(Sfe5, 2, 4, 
            Gen5) + 6*CA*CW*Mf2(4, Gen5)*USf(Sfe5, 2, 4, Gen5))*
         USfC(Sfe6, 2, 4, Gen5)))/(6*CB*CW*MW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MZ, pM2, Mh0, pM3, MA0, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (CBA*EL)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW), pG5, 
    -(EL*MT)/(2*MW*SW), pG6, -(i_*EL*MW*S2B*SAB)/(2*CW2*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MZ, pM2, Mh0, pM3, MZ, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (EL*SBA)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW), pG5, 
    -(EL*MT)/(2*MW*SW), pG6, (i_*C2B*EL*MW*SAB)/(2*CW2*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MZ, pM2, MHH, pM3, MA0, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, -(EL*SBA)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW), pG5, 
    -(EL*MT)/(2*MW*SW), pG6, (i_*CAB*EL*MW*S2B)/(2*CW2*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MZ, pM2, MHH, pM3, MZ, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (CBA*EL)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW), pG5, 
    -(EL*MT)/(2*MW*SW), pG6, -(i_*C2B*CAB*EL*MW)/(2*CW2*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MZ, pM2, MSf(Sfe5, 2, Gen5), 
    pM3, MSf(Sfe5, 2, Gen5), pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    i_*EL, pG4, (EL*MT)/(2*MW*SW), pG5, -(EL*MT)/(2*MW*SW), pG6, 
    -(EL*Mf(2, Gen5)*((-(MUEC*TB) + Af(2, Gen5, Gen5))*
         USf(Sfe5, 2, 2, Gen5)*USfC(Sfe5, 1, 2, Gen5) + 
        (MUE*TB - AfC(2, Gen5, Gen5))*USf(Sfe5, 1, 2, Gen5)*
         USfC(Sfe5, 2, 2, Gen5)))/(2*MW*SW), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MZ, pM2, MSf(Sfe5, 2, Gen5), pM3, MSf(Sfe6, 2, Gen5), 
    pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, 
    -((i_*EL*SW)/CW), pG3, 
    -(i_*EL*((-1 + 2*SW2)*USf(Sfe6, 1, 2, Gen5)*USfC(Sfe5, 1, 2, 
          Gen5) + 2*SW2*USf(Sfe6, 2, 2, Gen5)*USfC(Sfe5, 2, 2, Gen5)))/
     (2*CW*SW), pG4, (EL*MT)/(2*MW*SW), pG5, -(EL*MT)/(2*MW*SW), pG6, 
    -(EL*Mf(2, Gen5)*((-(MUEC*TB) + Af(2, Gen5, Gen5))*
         USf(Sfe5, 2, 2, Gen5)*USfC(Sfe6, 1, 2, Gen5) + 
        (MUE*TB - AfC(2, Gen5, Gen5))*USf(Sfe5, 1, 2, Gen5)*
         USfC(Sfe6, 2, 2, Gen5)))/(2*MW*SW), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MZ, pM2, MSf(Sfe5, 3, Gen5), 
    pM3, MSf(Sfe5, 3, Gen5), pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    (-2*i_*EL)/3, pG4, (EL*MT)/(2*MW*SW), pG5, -(EL*MT)/(2*MW*SW), 
    pG6, (EL*Mf(3, Gen5)*((-MUEC + TB*Af(3, Gen5, Gen5))*
        USf(Sfe5, 2, 3, Gen5)*USfC(Sfe5, 1, 3, Gen5) + 
       (MUE - TB*AfC(3, Gen5, Gen5))*USf(Sfe5, 1, 3, Gen5)*
        USfC(Sfe5, 2, 3, Gen5)))/(2*MW*SW*TB), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MZ, pM2, MSf(Sfe5, 3, Gen5), pM3, MSf(Sfe6, 3, Gen5), 
    pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, 
    -((i_*EL*SW)/CW), pG3, 
    (i_*EL*((-3 + 4*SW2)*USf(Sfe6, 1, 3, Gen5)*USfC(Sfe5, 1, 3, 
         Gen5) + 4*SW2*USf(Sfe6, 2, 3, Gen5)*USfC(Sfe5, 2, 3, Gen5)))/
     (6*CW*SW), pG4, (EL*MT)/(2*MW*SW), pG5, -(EL*MT)/(2*MW*SW), pG6, 
    (EL*Mf(3, Gen5)*((-MUEC + TB*Af(3, Gen5, Gen5))*
        USf(Sfe5, 2, 3, Gen5)*USfC(Sfe6, 1, 3, Gen5) + 
       (MUE - TB*AfC(3, Gen5, Gen5))*USf(Sfe5, 1, 3, Gen5)*
        USfC(Sfe6, 2, 3, Gen5)))/(2*MW*SW*TB), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MZ, pM2, MSf(Sfe5, 4, Gen5), 
    pM3, MSf(Sfe5, 4, Gen5), pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    (i_*EL)/3, pG4, (EL*MT)/(2*MW*SW), pG5, -(EL*MT)/(2*MW*SW), pG6, 
    -(EL*Mf(4, Gen5)*((-(MUEC*TB) + Af(4, Gen5, Gen5))*
         USf(Sfe5, 2, 4, Gen5)*USfC(Sfe5, 1, 4, Gen5) + 
        (MUE*TB - AfC(4, Gen5, Gen5))*USf(Sfe5, 1, 4, Gen5)*
         USfC(Sfe5, 2, 4, Gen5)))/(2*MW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MZ, pM2, MSf(Sfe5, 4, Gen5), pM3, MSf(Sfe6, 4, Gen5), 
    pM4, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, 
    -((i_*EL*SW)/CW), pG3, 
    -(i_*EL*((-3 + 2*SW2)*USf(Sfe6, 1, 4, Gen5)*USfC(Sfe5, 1, 4, 
          Gen5) + 2*SW2*USf(Sfe6, 2, 4, Gen5)*USfC(Sfe5, 2, 4, Gen5)))/
     (6*CW*SW), pG4, (EL*MT)/(2*MW*SW), pG5, -(EL*MT)/(2*MW*SW), pG6, 
    -(EL*Mf(4, Gen5)*((-(MUEC*TB) + Af(4, Gen5, Gen5))*
         USf(Sfe5, 2, 4, Gen5)*USfC(Sfe6, 1, 4, Gen5) + 
        (MUE*TB - AfC(4, Gen5, Gen5))*USf(Sfe5, 1, 4, Gen5)*
         USfC(Sfe6, 2, 4, Gen5)))/(2*MW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)));

L iT2G3 = T2G3*(replace_(pM1, Mh0, pM2, 0, pG1, i_*EL, pG2, 
    i_*EL, pG3, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG4, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, -(i_*EL*MW*SBA)/(2*SW), pG6, 
    -(i_*EL)) + replace_(pM1, Mh0, pM2, 0, pG1, i_*EL, pG2, 
    i_*EL, pG3, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG4, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, -(i_*EL*MW*SBA)/(2*SW), pG6, 
    i_*EL) + replace_(pM1, Mh0, pM2, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), 
    pG5, -(i_*EL*MW*SBA)/(2*SW), pG6, -((i_*CW*EL)/SW)) + 
   replace_(pM1, Mh0, pM2, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG4, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, -(i_*EL*MW*SBA)/(2*SW), pG6, 
    (i_*CW*EL)/SW) + replace_(pM1, MHH, pM2, 0, pG1, i_*EL, pG2, 
    i_*EL, pG3, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG4, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, -(i_*CBA*EL*MW)/(2*SW), pG6, 
    -(i_*EL)) + replace_(pM1, MHH, pM2, 0, pG1, i_*EL, pG2, 
    i_*EL, pG3, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG4, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, -(i_*CBA*EL*MW)/(2*SW), pG6, 
    i_*EL) + replace_(pM1, MHH, pM2, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), 
    pG5, -(i_*CBA*EL*MW)/(2*SW), pG6, -((i_*CW*EL)/SW)) + 
   replace_(pM1, MHH, pM2, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG4, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, -(i_*CBA*EL*MW)/(2*SW), pG6, 
    (i_*CW*EL)/SW) + replace_(pM1, MZ, pM2, 0, pG1, i_*EL, pG2, 
    i_*EL, pG3, (EL*MT)/(2*MW*SW), pG4, -(EL*MT)/(2*MW*SW), pG5, 
    -(EL*MW)/(2*SW), pG6, -(i_*EL)) + replace_(pM1, MZ, pM2, 0, pG1, 
    i_*EL, pG2, i_*EL, pG3, (EL*MT)/(2*MW*SW), pG4, 
    -(EL*MT)/(2*MW*SW), pG5, (EL*MW)/(2*SW), pG6, i_*EL) + 
   replace_(pM1, MZ, pM2, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, (EL*MT)/(2*MW*SW), pG4, 
    -(EL*MT)/(2*MW*SW), pG5, -(EL*MW)/(2*SW), pG6, -((i_*CW*EL)/SW)) + 
   replace_(pM1, MZ, pM2, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, (EL*MT)/(2*MW*SW), pG4, 
    -(EL*MT)/(2*MW*SW), pG5, (EL*MW)/(2*SW), pG6, (i_*CW*EL)/SW));

L iT2G4 = T2G4*(replace_(pM1, Mh0, pM2, 0, pG1, i_*EL, pG2, 
    i_*EL, pG3, i_*EL, pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, 
    (i_*EL*MW*SBA)/SW) + replace_(pM1, Mh0, pM2, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*CW*EL)/SW, pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, 
    (i_*EL*MW*SBA)/SW) + replace_(pM1, MHH, pM2, 0, pG1, i_*EL, 
    pG2, i_*EL, pG3, i_*EL, pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, 
    (i_*CBA*EL*MW)/SW) + replace_(pM1, MHH, pM2, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*CW*EL)/SW, pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, 
    (i_*CBA*EL*MW)/SW));

L iT2G5 = T2G5*(replace_(pM1, MA0, pM2, Mh0, pM3, MZ, pM4, MZ, 
    pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, -(CBA*EL)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW*TB), pG5, 
    -(EL*MT)/(2*MW*SW*TB), pG6, (i_*EL*MW*SBA)/(CW2*SW)) + 
   replace_(pM1, MA0, pM2, MHH, pM3, MZ, pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (EL*SBA)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW*TB), pG5, 
    -(EL*MT)/(2*MW*SW*TB), pG6, (i_*CBA*EL*MW)/(CW2*SW)) + 
   replace_(pM1, Mh0, pM2, MW, pM3, 0, pM4, MW, pG1, i_*EL, pG2, 
    i_*EL, pG3, -(i_*EL*SBA)/(2*SW), pG4, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, -(i_*CA*EL*MT)/(2*MW*SB*SW), 
    pG6, i_*EL*MW) + replace_(pM1, Mh0, pM2, MW, pM3, 0, pM4, MW, 
    pG1, i_*EL, pG2, i_*EL, pG3, (i_*EL*SBA)/(2*SW), pG4, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, -(i_*CA*EL*MT)/(2*MW*SB*SW), 
    pG6, i_*EL*MW) + replace_(pM1, Mh0, pM2, MW, pM3, MZ, pM4, MW, 
    pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, -(i_*EL*SBA)/(2*SW), pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG6, -((i_*EL*MW*SW)/CW)) + 
   replace_(pM1, Mh0, pM2, MW, pM3, MZ, pM4, MW, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*SBA)/(2*SW), pG4, -(i_*CA*EL*MT)/(2*MW*SB*SW), pG5, 
    -(i_*CA*EL*MT)/(2*MW*SB*SW), pG6, -((i_*EL*MW*SW)/CW)) + 
   replace_(pM1, MHH, pM2, MW, pM3, 0, pM4, MW, pG1, i_*EL, pG2, 
    i_*EL, pG3, -(i_*CBA*EL)/(2*SW), pG4, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, -(i_*EL*MT*SA)/(2*MW*SB*SW), 
    pG6, i_*EL*MW) + replace_(pM1, MHH, pM2, MW, pM3, 0, pM4, MW, 
    pG1, i_*EL, pG2, i_*EL, pG3, (i_*CBA*EL)/(2*SW), pG4, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, -(i_*EL*MT*SA)/(2*MW*SB*SW), 
    pG6, i_*EL*MW) + replace_(pM1, MHH, pM2, MW, pM3, MZ, pM4, MW, 
    pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, -(i_*CBA*EL)/(2*SW), pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG6, -((i_*EL*MW*SW)/CW)) + 
   replace_(pM1, MHH, pM2, MW, pM3, MZ, pM4, MW, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*CBA*EL)/(2*SW), pG4, -(i_*EL*MT*SA)/(2*MW*SB*SW), pG5, 
    -(i_*EL*MT*SA)/(2*MW*SB*SW), pG6, -((i_*EL*MW*SW)/CW)) + 
   replace_(pM1, MZ, pM2, Mh0, pM3, MZ, pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    -(EL*SBA)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW), pG5, 
    -(EL*MT)/(2*MW*SW), pG6, (i_*EL*MW*SBA)/(CW2*SW)) + 
   replace_(pM1, MZ, pM2, MHH, pM3, MZ, pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    -(CBA*EL)/(2*CW*SW), pG4, (EL*MT)/(2*MW*SW), pG5, 
    -(EL*MT)/(2*MW*SW), pG6, (i_*CBA*EL*MW)/(CW2*SW)) + 
   2*replace_(pM1, MZ, pM2, MW, pM3, 0, pM4, MW, pG1, i_*EL, pG2, 
     i_*EL, pG3, EL/(2*SW), pG4, (EL*MT)/(2*MW*SW), pG5, 
     -(EL*MT)/(2*MW*SW), pG6, i_*EL*MW) + 
   2*replace_(pM1, MZ, pM2, MW, pM3, MZ, pM4, MW, pG1, 
     -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
     EL/(2*SW), pG4, (EL*MT)/(2*MW*SW), pG5, -(EL*MT)/(2*MW*SW), pG6, 
     -((i_*EL*MW*SW)/CW)));

L iT2G6 = T2G6*(replace_(pM1, 0, pM2, 0, pM3, MZ, pM4, MZ, cG1, 
    -(i_*EL)/(2*CW*SW), cG2, 0, cG3, -(i_*EL)/(2*CW*SW), cG4, 0, 
    pG5, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG6, -((i_*EL*SW)/CW), 
    pG7, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG8, (2*i_*EL*SW)/(3*CW), 
    cR1, 6*SUNT(Col3, Col4)) + replace_(pM1, MCha(Cha5), pM2, 
    MCha(Cha5), pM3, 0, pM4, 0, cG1, i_*EL, cG2, i_*EL, cG3, 
    i_*EL, cG4, i_*EL, pG5, i_*EL, pG6, i_*EL, pG7, 
    (-2*i_*EL)/3, pG8, (-2*i_*EL)/3, cR1, 2*SumOver(Cha5, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MCha(Cha5), pM2, MCha(Cha5), 
    pM3, 0, pM4, MZ, cG1, 
    -((i_*EL*(SW2 - UCha(Cha5, 1)*UChaC(Cha5, 1) - 
        (UCha(Cha5, 2)*UChaC(Cha5, 2))/2))/(CW*SW)), cG2, 
    -((i_*EL*(SW2 - VCha(Cha5, 1)*VChaC(Cha5, 1) - 
        (VCha(Cha5, 2)*VChaC(Cha5, 2))/2))/(CW*SW)), cG3, i_*EL, cG4, 
    i_*EL, pG5, i_*EL, pG6, i_*EL, pG7, (i_*EL*(-3 + 4*SW2))/
     (6*CW*SW), pG8, (2*i_*EL*SW)/(3*CW), cR1, 
    2*SumOver(Cha5, 2)*SUNT(Col3, Col4)) + replace_(pM1, MCha(Cha5), 
    pM2, MCha(Cha5), pM3, MZ, pM4, 0, cG1, i_*EL, cG2, i_*EL, cG3, 
    -((i_*EL*(SW2 - UCha(Cha5, 1)*UChaC(Cha5, 1) - 
        (UCha(Cha5, 2)*UChaC(Cha5, 2))/2))/(CW*SW)), cG4, 
    -((i_*EL*(SW2 - VCha(Cha5, 1)*VChaC(Cha5, 1) - 
        (VCha(Cha5, 2)*VChaC(Cha5, 2))/2))/(CW*SW)), pG5, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG6, -((i_*EL*SW)/CW), pG7, 
    (-2*i_*EL)/3, pG8, (-2*i_*EL)/3, cR1, 2*SumOver(Cha5, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MCha(Cha5), pM2, MCha(Cha6), 
    pM3, MZ, pM4, MZ, cG1, 
    -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - UCha(Cha5, 1)*
         UChaC(Cha6, 1) - (UCha(Cha5, 2)*UChaC(Cha6, 2))/2))/(CW*SW)), 
    cG2, -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - 
        VCha(Cha6, 1)*VChaC(Cha5, 1) - (VCha(Cha6, 2)*VChaC(Cha5, 2))/
         2))/(CW*SW)), cG3, 
    -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - UCha(Cha6, 1)*
         UChaC(Cha5, 1) - (UCha(Cha6, 2)*UChaC(Cha5, 2))/2))/(CW*SW)), 
    cG4, -((i_*EL*(SW2*IndexDelta(Cha5, Cha6) - 
        VCha(Cha5, 1)*VChaC(Cha6, 1) - (VCha(Cha5, 2)*VChaC(Cha6, 2))/
         2))/(CW*SW)), pG5, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG6, 
    -((i_*EL*SW)/CW), pG7, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG8, 
    (2*i_*EL*SW)/(3*CW), cR1, 2*SumOver(Cha5, 2)*SumOver(Cha6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, Mf(2, Gen5), pM2, Mf(2, Gen5), 
    pM3, 0, pM4, 0, cG1, i_*EL, cG2, i_*EL, cG3, i_*EL, cG4, 
    i_*EL, pG5, i_*EL, pG6, i_*EL, pG7, (-2*i_*EL)/3, pG8, 
    (-2*i_*EL)/3, cR1, 2*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(2, Gen5), pM2, Mf(2, Gen5), pM3, 0, pM4, MZ, cG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), cG2, -((i_*EL*SW)/CW), cG3, 
    i_*EL, cG4, i_*EL, pG5, i_*EL, pG6, i_*EL, pG7, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG8, (2*i_*EL*SW)/(3*CW), cR1, 
    2*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + replace_(pM1, Mf(2, Gen5), 
    pM2, Mf(2, Gen5), pM3, MZ, pM4, 0, cG1, i_*EL, cG2, i_*EL, cG3, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), cG4, -((i_*EL*SW)/CW), pG5, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG6, -((i_*EL*SW)/CW), pG7, 
    (-2*i_*EL)/3, pG8, (-2*i_*EL)/3, cR1, 2*SumOver(Gen5, 3)*
     SUNT(Col3, Col4)) + replace_(pM1, Mf(2, Gen5), pM2, Mf(2, Gen5), 
    pM3, MZ, pM4, MZ, cG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), cG2, 
    -((i_*EL*SW)/CW), cG3, -((i_*EL*(-1/2 + SW2))/(CW*SW)), cG4, 
    -((i_*EL*SW)/CW), pG5, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG6, 
    -((i_*EL*SW)/CW), pG7, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG8, 
    (2*i_*EL*SW)/(3*CW), cR1, 2*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(3, Gen5), pM2, Mf(3, Gen5), pM3, 0, pM4, 0, cG1, 
    (-2*i_*EL)/3, cG2, (-2*i_*EL)/3, cG3, (-2*i_*EL)/3, cG4, 
    (-2*i_*EL)/3, pG5, i_*EL, pG6, i_*EL, pG7, (-2*i_*EL)/3, 
    pG8, (-2*i_*EL)/3, cR1, 6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(3, Gen5), pM2, Mf(3, Gen5), pM3, 0, pM4, MZ, cG1, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), cG2, (2*i_*EL*SW)/(3*CW), cG3, 
    (-2*i_*EL)/3, cG4, (-2*i_*EL)/3, pG5, i_*EL, pG6, i_*EL, 
    pG7, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG8, (2*i_*EL*SW)/(3*CW), 
    cR1, 6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(3, Gen5), pM2, Mf(3, Gen5), pM3, MZ, pM4, 0, cG1, 
    (-2*i_*EL)/3, cG2, (-2*i_*EL)/3, cG3, (i_*EL*(-3 + 4*SW2))/
     (6*CW*SW), cG4, (2*i_*EL*SW)/(3*CW), pG5, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG6, -((i_*EL*SW)/CW), pG7, 
    (-2*i_*EL)/3, pG8, (-2*i_*EL)/3, cR1, 6*SumOver(Gen5, 3)*
     SUNT(Col3, Col4)) + replace_(pM1, Mf(3, Gen5), pM2, Mf(3, Gen5), 
    pM3, MZ, pM4, MZ, cG1, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), cG2, 
    (2*i_*EL*SW)/(3*CW), cG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), cG4, 
    (2*i_*EL*SW)/(3*CW), pG5, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG6, 
    -((i_*EL*SW)/CW), pG7, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG8, 
    (2*i_*EL*SW)/(3*CW), cR1, 6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(4, Gen5), pM2, Mf(4, Gen5), pM3, 0, pM4, 0, cG1, 
    (i_*EL)/3, cG2, (i_*EL)/3, cG3, (i_*EL)/3, cG4, (i_*EL)/3, 
    pG5, i_*EL, pG6, i_*EL, pG7, (-2*i_*EL)/3, pG8, 
    (-2*i_*EL)/3, cR1, 6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, Mf(4, Gen5), pM2, Mf(4, Gen5), pM3, 0, pM4, MZ, cG1, 
    -(i_*EL*(-3 + 2*SW2))/(6*CW*SW), cG2, -(i_*EL*SW)/(3*CW), cG3, 
    (i_*EL)/3, cG4, (i_*EL)/3, pG5, i_*EL, pG6, i_*EL, pG7, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG8, (2*i_*EL*SW)/(3*CW), cR1, 
    6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + replace_(pM1, Mf(4, Gen5), 
    pM2, Mf(4, Gen5), pM3, MZ, pM4, 0, cG1, (i_*EL)/3, cG2, 
    (i_*EL)/3, cG3, -(i_*EL*(-3 + 2*SW2))/(6*CW*SW), cG4, 
    -(i_*EL*SW)/(3*CW), pG5, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG6, 
    -((i_*EL*SW)/CW), pG7, (-2*i_*EL)/3, pG8, (-2*i_*EL)/3, cR1, 
    6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + replace_(pM1, Mf(4, Gen5), 
    pM2, Mf(4, Gen5), pM3, MZ, pM4, MZ, cG1, -(i_*EL*(-3 + 2*SW2))/
     (6*CW*SW), cG2, -(i_*EL*SW)/(3*CW), cG3, 
    -(i_*EL*(-3 + 2*SW2))/(6*CW*SW), cG4, -(i_*EL*SW)/(3*CW), pG5, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG6, -((i_*EL*SW)/CW), pG7, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG8, (2*i_*EL*SW)/(3*CW), cR1, 
    6*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + replace_(pM1, MNeu(Neu5), 
    pM2, MNeu(Neu6), pM3, MZ, pM4, MZ, cG1, 
    (i_*EL*(-(ZNeu(Neu5, 3)*ZNeuC(Neu6, 3)) + ZNeu(Neu5, 4)*
        ZNeuC(Neu6, 4)))/(2*CW*SW), cG2, 
    (i_*EL*(ZNeu(Neu6, 3)*ZNeuC(Neu5, 3) - ZNeu(Neu6, 4)*
        ZNeuC(Neu5, 4)))/(2*CW*SW), cG3, 
    (i_*EL*(-(ZNeu(Neu6, 3)*ZNeuC(Neu5, 3)) + ZNeu(Neu6, 4)*
        ZNeuC(Neu5, 4)))/(2*CW*SW), cG4, 
    (i_*EL*(ZNeu(Neu5, 3)*ZNeuC(Neu6, 3) - ZNeu(Neu5, 4)*
        ZNeuC(Neu6, 4)))/(2*CW*SW), pG5, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG6, -((i_*EL*SW)/CW), pG7, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG8, (2*i_*EL*SW)/(3*CW), cR1, 
    SumOver(Neu5, 4)*SumOver(Neu6, 4)*SUNT(Col3, Col4)));

L iT2G7 = T2G7*(replace_(pM1, Mh0, pM2, MA0, pM3, MZ, pM4, MZ, 
    pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), 
    pG5, (CBA*EL)/(2*CW*SW), pG6, (CBA*EL)/(2*CW*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, Mh0, pM2, MZ, pM3, MZ, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), 
    pG5, (EL*SBA)/(2*CW*SW), pG6, (EL*SBA)/(2*CW*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MHH, pM2, MA0, pM3, MZ, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), 
    pG5, -(EL*SBA)/(2*CW*SW), pG6, -(EL*SBA)/(2*CW*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MHH, pM2, MZ, pM3, MZ, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), 
    pG5, (CBA*EL)/(2*CW*SW), pG6, (CBA*EL)/(2*CW*SW), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MHp, pM2, MHp, pM3, 0, pM4, 0, 
    pG1, i_*EL, pG2, i_*EL, pG3, (-2*i_*EL)/3, pG4, 
    (-2*i_*EL)/3, pG5, i_*EL, pG6, -(i_*EL), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MHp, pM2, MHp, pM3, 0, pM4, 
    MZ, pG1, i_*EL, pG2, i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/
     (6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, i_*EL, pG6, 
    -(i_*EL*(CW2 - SW2))/(2*CW*SW), cR1, 2*SUNT(Col3, Col4)) + 
   replace_(pM1, MHp, pM2, MHp, pM3, MZ, pM4, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, (i_*EL*(CW2 - SW2))/
     (2*CW*SW), pG6, -(i_*EL), cR1, 2*SUNT(Col3, Col4)) + 
   replace_(pM1, MHp, pM2, MHp, pM3, MZ, pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (i_*EL*(CW2 - SW2))/(2*CW*SW), pG6, -(i_*EL*(CW2 - SW2))/
     (2*CW*SW), cR1, 2*SUNT(Col3, Col4)) + replace_(pM1, MW, pM2, MW, 
    pM3, 0, pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, (-2*i_*EL)/3, 
    pG4, (-2*i_*EL)/3, pG5, i_*EL, pG6, -(i_*EL), cR1, 
    2*SUNT(Col3, Col4)) + replace_(pM1, MW, pM2, MW, pM3, 0, pM4, MZ, 
    pG1, i_*EL, pG2, i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), 
    pG4, (2*i_*EL*SW)/(3*CW), pG5, i_*EL, pG6, 
    -(i_*EL*(CW2 - SW2))/(2*CW*SW), cR1, 2*SUNT(Col3, Col4)) + 
   replace_(pM1, MW, pM2, MW, pM3, MZ, pM4, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, (i_*EL*(CW2 - SW2))/
     (2*CW*SW), pG6, -(i_*EL), cR1, 2*SUNT(Col3, Col4)) + 
   replace_(pM1, MW, pM2, MW, pM3, MZ, pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (i_*EL*(CW2 - SW2))/(2*CW*SW), pG6, -(i_*EL*(CW2 - SW2))/
     (2*CW*SW), cR1, 2*SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(1, 1, Gen5), pM2, MSf(1, 1, Gen5), pM3, MZ, pM4, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), 
    pG5, -(i_*EL)/(2*CW*SW), pG6, (i_*EL)/(2*CW*SW), cR1, 
    2*SumOver(Gen5, 3)*SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(Sfe5, 2, Gen5), pM2, MSf(Sfe5, 2, Gen5), pM3, 0, 
    pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, (-2*i_*EL)/3, pG4, 
    (-2*i_*EL)/3, pG5, i_*EL, pG6, -(i_*EL), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + 
   replace_(pM1, MSf(Sfe5, 2, Gen5), pM2, MSf(Sfe5, 2, Gen5), pM3, 0, 
    pM4, MZ, pG1, i_*EL, pG2, i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/
     (6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, i_*EL, pG6, 
    (i_*EL*((-1 + 2*SW2)*USf(Sfe5, 1, 2, Gen5)*USfC(Sfe5, 1, 2, 
         Gen5) + 2*SW2*USf(Sfe5, 2, 2, Gen5)*USfC(Sfe5, 2, 2, Gen5)))/
     (2*CW*SW), cR1, 2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MSf(Sfe5, 2, Gen5), pM2, 
    MSf(Sfe5, 2, Gen5), pM3, MZ, pM4, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    -(i_*EL*((-1 + 2*SW2)*USf(Sfe5, 1, 2, Gen5)*USfC(Sfe5, 1, 2, 
          Gen5) + 2*SW2*USf(Sfe5, 2, 2, Gen5)*USfC(Sfe5, 2, 2, Gen5)))/
     (2*CW*SW), pG6, -(i_*EL), cR1, 2*SumOver(Gen5, 3)*
     SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + replace_(pM1, 
    MSf(Sfe5, 2, Gen5), pM2, MSf(Sfe6, 2, Gen5), pM3, MZ, pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    -(i_*EL*((-1 + 2*SW2)*USf(Sfe6, 1, 2, Gen5)*USfC(Sfe5, 1, 2, 
          Gen5) + 2*SW2*USf(Sfe6, 2, 2, Gen5)*USfC(Sfe5, 2, 2, Gen5)))/
     (2*CW*SW), pG6, (i_*EL*((-1 + 2*SW2)*USf(Sfe5, 1, 2, Gen5)*
        USfC(Sfe6, 1, 2, Gen5) + 2*SW2*USf(Sfe5, 2, 2, Gen5)*
        USfC(Sfe6, 2, 2, Gen5)))/(2*CW*SW), cR1, 
    2*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MSf(Sfe5, 3, Gen5), pM2, 
    MSf(Sfe5, 3, Gen5), pM3, 0, pM4, 0, pG1, i_*EL, pG2, i_*EL, pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, (-2*i_*EL)/3, pG6, 
    (2*i_*EL)/3, cR1, 6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MSf(Sfe5, 3, Gen5), pM2, 
    MSf(Sfe5, 3, Gen5), pM3, 0, pM4, MZ, pG1, i_*EL, pG2, i_*EL, 
    pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), 
    pG5, (-2*i_*EL)/3, pG6, 
    -(i_*EL*((-3 + 4*SW2)*USf(Sfe5, 1, 3, Gen5)*USfC(Sfe5, 1, 3, 
          Gen5) + 4*SW2*USf(Sfe5, 2, 3, Gen5)*USfC(Sfe5, 2, 3, Gen5)))/
     (6*CW*SW), cR1, 6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MSf(Sfe5, 3, Gen5), pM2, 
    MSf(Sfe5, 3, Gen5), pM3, MZ, pM4, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    (i_*EL*((-3 + 4*SW2)*USf(Sfe5, 1, 3, Gen5)*USfC(Sfe5, 1, 3, 
         Gen5) + 4*SW2*USf(Sfe5, 2, 3, Gen5)*USfC(Sfe5, 2, 3, Gen5)))/
     (6*CW*SW), pG6, (2*i_*EL)/3, cR1, 6*SumOver(Gen5, 3)*
     SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + replace_(pM1, 
    MSf(Sfe5, 3, Gen5), pM2, MSf(Sfe6, 3, Gen5), pM3, MZ, pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (i_*EL*((-3 + 4*SW2)*USf(Sfe6, 1, 3, Gen5)*USfC(Sfe5, 1, 3, 
         Gen5) + 4*SW2*USf(Sfe6, 2, 3, Gen5)*USfC(Sfe5, 2, 3, Gen5)))/
     (6*CW*SW), pG6, 
    -(i_*EL*((-3 + 4*SW2)*USf(Sfe5, 1, 3, Gen5)*USfC(Sfe6, 1, 3, 
          Gen5) + 4*SW2*USf(Sfe5, 2, 3, Gen5)*USfC(Sfe6, 2, 3, Gen5)))/
     (6*CW*SW), cR1, 6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*
     SumOver(Sfe6, 2)*SUNT(Col3, Col4)) + replace_(pM1, 
    MSf(Sfe5, 4, Gen5), pM2, MSf(Sfe5, 4, Gen5), pM3, 0, pM4, 0, pG1, 
    i_*EL, pG2, i_*EL, pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, 
    pG5, (i_*EL)/3, pG6, -(i_*EL)/3, cR1, 6*SumOver(Gen5, 3)*
     SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + replace_(pM1, 
    MSf(Sfe5, 4, Gen5), pM2, MSf(Sfe5, 4, Gen5), pM3, 0, pM4, MZ, pG1, 
    i_*EL, pG2, i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, 
    (2*i_*EL*SW)/(3*CW), pG5, (i_*EL)/3, pG6, 
    (i_*EL*((-3 + 2*SW2)*USf(Sfe5, 1, 4, Gen5)*USfC(Sfe5, 1, 4, 
         Gen5) + 2*SW2*USf(Sfe5, 2, 4, Gen5)*USfC(Sfe5, 2, 4, Gen5)))/
     (6*CW*SW), cR1, 6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*
     SUNT(Col3, Col4)) + replace_(pM1, MSf(Sfe5, 4, Gen5), pM2, 
    MSf(Sfe5, 4, Gen5), pM3, MZ, pM4, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    -(i_*EL*((-3 + 2*SW2)*USf(Sfe5, 1, 4, Gen5)*USfC(Sfe5, 1, 4, 
          Gen5) + 2*SW2*USf(Sfe5, 2, 4, Gen5)*USfC(Sfe5, 2, 4, Gen5)))/
     (6*CW*SW), pG6, -(i_*EL)/3, cR1, 6*SumOver(Gen5, 3)*
     SumOver(Sfe5, 2)*SUNT(Col3, Col4)) + replace_(pM1, 
    MSf(Sfe5, 4, Gen5), pM2, MSf(Sfe6, 4, Gen5), pM3, MZ, pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    -(i_*EL*((-3 + 2*SW2)*USf(Sfe6, 1, 4, Gen5)*USfC(Sfe5, 1, 4, 
          Gen5) + 2*SW2*USf(Sfe6, 2, 4, Gen5)*USfC(Sfe5, 2, 4, Gen5)))/
     (6*CW*SW), pG6, (i_*EL*((-3 + 2*SW2)*USf(Sfe5, 1, 4, Gen5)*
        USfC(Sfe6, 1, 4, Gen5) + 2*SW2*USf(Sfe5, 2, 4, Gen5)*
        USfC(Sfe6, 2, 4, Gen5)))/(6*CW*SW), cR1, 
    6*SumOver(Gen5, 3)*SumOver(Sfe5, 2)*SumOver(Sfe6, 2)*
     SUNT(Col3, Col4)));

L iT2G8 = T2G8*(replace_(pM1, 0, pM2, 0, pG1, i_*EL, pG2, 
    i_*EL, pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, -(i_*EL), 
    pG6, -(i_*EL)) + replace_(pM1, 0, pM2, 0, pG1, i_*EL, pG2, 
    i_*EL, pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, i_*EL, 
    pG6, i_*EL) + replace_(pM1, 0, pM2, MZ, pG1, i_*EL, pG2, 
    i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, 
    (2*i_*EL*SW)/(3*CW), pG5, -((i_*CW*EL)/SW), pG6, -(i_*EL)) + 
   replace_(pM1, 0, pM2, MZ, pG1, i_*EL, pG2, i_*EL, pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (i_*CW*EL)/SW, pG6, i_*EL) + replace_(pM1, MZ, pM2, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, -(i_*EL), pG6, 
    -((i_*CW*EL)/SW)) + replace_(pM1, MZ, pM2, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, i_*EL, pG6, 
    (i_*CW*EL)/SW) + replace_(pM1, MZ, pM2, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    -((i_*CW*EL)/SW), pG6, -((i_*CW*EL)/SW)) + 
   replace_(pM1, MZ, pM2, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, 
    (2*i_*EL*SW)/(3*CW), pG5, (i_*CW*EL)/SW, pG6, 
    (i_*CW*EL)/SW));

L iT2G9 = T2G9*(replace_(pM1, 0, pM2, 0, pG1, i_*EL, pG2, 
    i_*EL, pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, i_*EL, 
    pG6, -(i_*EL)) + replace_(pM1, 0, pM2, MZ, pG1, i_*EL, pG2, 
    i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, 
    (2*i_*EL*SW)/(3*CW), pG5, i_*EL, pG6, -((i_*CW*EL)/SW)) + 
   replace_(pM1, MZ, pM2, 0, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, 
    pG5, (i_*CW*EL)/SW, pG6, -(i_*EL)) + replace_(pM1, MZ, pM2, 
    MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), 
    pG5, (i_*CW*EL)/SW, pG6, -((i_*CW*EL)/SW)));

L iT2G10 = T2G10*(replace_(pM1, Mh0, pM2, MZ, pM3, MZ, pM4, MZ, 
    pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), 
    pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), 
    pG5, (i_*EL*MW*SBA)/(CW2*SW), pG6, (i_*EL*MW*SBA)/(CW2*SW)) + 
   replace_(pM1, MHH, pM2, MZ, pM3, MZ, pM4, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
    (i_*CBA*EL*MW)/(CW2*SW), pG6, (i_*CBA*EL*MW)/(CW2*SW)) + 
   2*replace_(pM1, MW, pM2, 0, pM3, 0, pM4, MW, pG1, i_*EL, pG2, 
     i_*EL, pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, i_*EL*MW, 
     pG6, i_*EL*MW) + 2*replace_(pM1, MW, pM2, 0, pM3, MZ, pM4, MW, 
     pG1, i_*EL, pG2, i_*EL, pG3, (i_*EL*(-3 + 4*SW2))/(6*CW*SW), 
     pG4, (2*i_*EL*SW)/(3*CW), pG5, i_*EL*MW, pG6, 
     -((i_*EL*MW*SW)/CW)) + 2*replace_(pM1, MW, pM2, MZ, pM3, 0, 
     pM4, MW, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, 
     -((i_*EL*SW)/CW), pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
     -((i_*EL*MW*SW)/CW), pG6, i_*EL*MW) + 
   2*replace_(pM1, MW, pM2, MZ, pM3, MZ, pM4, MW, pG1, 
     -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
     (i_*EL*(-3 + 4*SW2))/(6*CW*SW), pG4, (2*i_*EL*SW)/(3*CW), pG5, 
     -((i_*EL*MW*SW)/CW), pG6, -((i_*EL*MW*SW)/CW)));

L iT3G1 = T3G1*(replace_(pM1, 0, pG1, i_*EL, pG2, i_*EL, 
    pG3, (dZZA1*MZ)/4, pG4, -(dZZA1*MZ)/4) + replace_(pM1, MZ, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (MZ*(dZG01 + dZZZ1 + dMZsq1/MZ2))/4, pG4, 
    (MZ*(-dZG01 - dZZZ1 - dMZsq1/MZ2))/4));

L iT3G2 = T3G2*(replace_(pM1, 0, pM2, 0, pG1, i_*EL, pG2, 
    i_*EL, pG3, (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, 
    -(i_*dZAA1), pG6, 0, pG7, i_*dZAA1) + 
   replace_(pM1, 0, pM2, MZ, pG1, i_*EL, pG2, i_*EL, pG3, 
    -((i_*EL*(1/2 - (2*SW2)/3))/(CW*SW)), pG4, (2*i_*EL*SW)/(3*CW), 
    pG5, i_*(-dZAZ1/2 - dZZA1/2), pG6, (i_*dZZA1*MZ2)/2, pG7, 
    i_*(dZAZ1/2 + dZZA1/2)) + replace_(pM1, MZ, pM2, 0, pG1, 
    -((i_*EL*(-1/2 + SW2))/(CW*SW)), pG2, -((i_*EL*SW)/CW), pG3, 
    (-2*i_*EL)/3, pG4, (-2*i_*EL)/3, pG5, i_*(-dZAZ1/2 - dZZA1/2), 
    pG6, (i_*dZZA1*MZ2)/2, pG7, i_*(dZAZ1/2 + dZZA1/2)) + 
   replace_(pM1, MZ, pM2, MZ, pG1, -((i_*EL*(-1/2 + SW2))/(CW*SW)), 
    pG2, -((i_*EL*SW)/CW), pG3, -((i_*EL*(1/2 - (2*SW2)/3))/
      (CW*SW)), pG4, (2*i_*EL*SW)/(3*CW), pG5, -(i_*dZZZ1), pG6, 
    i_*(dMZsq1 + dZZZ1*MZ2), pG7, i_*dZZZ1));

.sort
drop;

L Result = iT1G1 + iT1G2 + iT2G1 + iT2G10 + iT2G2 + 
  iT2G3 + iT2G4 + iT2G5 + iT2G6 + iT2G7 + iT2G8 + iT2G9 + 
  iT3G1 + iT3G2;

#endprocedure


* CalcFeynAmp.frm
* the FORM part of the CalcFeynAmp function
* this file is part of FormCalc
* last modified 24 Sep 14 th


#procedure Contract
repeat once e_([i]?, [j]?, [k]?, [LA]?)*e_([I]?, [J]?, [K]?, [LA]?) =
#if "`Dim'" != "4"
  (1 + Dminus4) *
#endif
  ( d_([i], [I]) * (d_([j], [J])*d_([k], [K]) - d_([j], [K])*d_([k], [J])) +
    d_([i], [J]) * (d_([j], [K])*d_([k], [I]) - d_([j], [I])*d_([k], [K])) +
    d_([i], [K]) * (d_([j], [I])*d_([k], [J]) - d_([j], [J])*d_([k], [I])) );
#endprocedure

***********************************************************************

#if "`MomElim'" == "Automatic"
#define MomRange "1, `Legs'"
#elseif `MomElim'
#define MomRange "`MomElim', `MomElim'"
#endif

#procedure MomConserv(foo, r1, r2)
#do rep1 = 1, `r1'
#ifdef `MomRange'
id `foo'([x]?) = `foo'(nterms_([x]), [x]);

#do rep2 = 1, `r2'
#do i = `MomRange'
#ifdef `k`i''
id `foo'([n]?, [x]?) = `foo'([n], [x]) * NOW([x]);
argument NOW;
id k`i' = `k`i'';
#call eiki
endargument;

id NOW(0) = 0;
id NOW([x]?) = `foo'(nterms_([x]), [x]);
once `foo'(?a) = `foo'(?a);
also `foo'(?a) = 1;
#endif
#enddo
#enddo

id `foo'([n]?, [x]?) = `foo'([x]);
#endif

argument `foo';
#call kikj
#call Square
endargument;

#call InvSimplify(`foo')
id `foo'(0) = 0;
#enddo
#endprocedure

***********************************************************************

#procedure DotSimplify
#call eiki

id GA(?g) = GF(?g);

id q1 = q1 * QTAG;
id e_([mu]?, [nu]?, [ro]?, [si]?) =
  e_([mu], [nu], [ro], [si]) * TMP([mu], [nu], [ro], [si]) * ETAG;
id [t]?(?i) = [t](?i) * TMP(?i);
chainout TMP;
id TMP([p1]?) = 1;
id TMP([mu]?)^2 = 1;
id TMP([mu]?) = TAG;
id ETAG^[n]?{>1} = ETAG;

ab k1,...,k`Legs';
.sort
on oldFactArg;

collect dotM, dotM, 50;
makeinteger dotM;

id ETAG = 1;
id QTAG = TAG;

b dotM;
.sort
keep brackets;

#call MomConserv(dotM, 1, 2)

#if `DotExpand' == 1

id dotM([x]?) = [x];

.sort
off oldFactArg;

id TAG = 1;

#else

factarg dotM;
chainout dotM;
makeinteger dotM;
id dotM(1) = 1;

ab `Vectors', `Invariants', dotM;
.sort
off oldFactArg;

collect dotM, dotM;

#call InvSimplify(dotM)
id dotM(0) = 0;

repeat id TAG * dotM([x]?) = TAG * [x];
id TAG = 1;

*makeinteger dotM;
*id dotM(dotM(?x)) = dotM(?x);

argument dotM;
id dotM([x]?) = dotM(nterms_([x]), [x]);
id dotM(1, [x]?) = [x];
id dotM([n]?, [x]?) = dotM([x]);
argument dotM;
toPolynomial;
endargument;
toPolynomial;
endargument;

makeinteger dotM;
id dotM(1) = 1;
id dotM([x]?^[n]?) = dotM([x])^[n];
id dotM([x]?INVS) = [x];

toPolynomial onlyfunctions dotM;

.sort

#endif

id GF(?g) = GA(?g);
#endprocedure

***********************************************************************

#procedure ChainOrder(contr)
label 1;

#if "`OnShell'" == "True"
* Apply Dirac equation to right spinor
repeat;
  id GA([om]?, ?a, [p1]?, ?b) * Spinor([p1]?, [m1]?, [s1]?, ?s1) =
    ( 2*GD([om], ?a, [p1]) * distrib_(-1, 1, GD, GD, ?b) +
      sign_(nargs_(?b)) * [s1]*[m1] * GA([om], ?a, ?b) ) *
    Spinor([p1], [m1], [s1], ?s1);
  id GD(?a, [p1]?) * GD([mu]?) * GD(?b) =
    d_([p1], [mu]) * GA(?a, ?b) * TAG;
endrepeat;

* Apply Dirac equation to left spinor
repeat;
  id Spinor([p1]?, [m1]?, [s1]?, ?s1) * GA([om]?{6,7}[[n]], ?a, [p1]?, ?b) =
    Spinor([p1], [m1], [s1], ?s1) * sign_(nargs_(?a)) *
    ( [s1]*[m1] * GA({7,6}[[n]], ?a, ?b) -
      2*distrib_(-1, 1, GD, GD, ?a) * GD([p1], [om], ?b) );
  id GD([mu]?) * GD(?a) * GD([p1]?, [om]?, ?b) =
    d_([p1], [mu]) * GA([om], ?a, ?b) * TAG;
endrepeat;
#endif

* Eliminate contractions within each Dirac chain using the
* formulas from M. Veltman's Gammatrica [Nucl Phys B319 (1989) 253]
id GA([om]?, [mu]?, ?a) = GA([om]) * GB([mu], ?a);
while( count(GB,1) );
#if `contr'
  repeat;
    id GB([LA]?, [LA]?, ?a) = d_([LA], [LA]) * GB(?a);
    also GB([LA]?, [mu]?, [LA]?, ?a) = (2 - d_([LA], [LA])) * GB([mu], ?a);
    also GB([LA]?, [mu]?, [nu]?, [LA]?, ?a) =
#if "`Dim'" != "4"
      Dminus4 * GB([mu], [nu], ?a) +
#endif
      4*d_([mu], [nu]) * GB(?a) * TAG;
    also GB([LA]?, [mu]?, [nu]?, [ro]?, ?b, [LA]?, ?a) =
      -sign_(nargs_(?b)) * (
#if "`Dim'" != "4"
        Dminus4 * GB([mu], [nu], [ro], ?b, ?a) +
#endif
        2*GB([ro], [nu], [mu], ?b, ?a) +
        2*GD([mu], [nu], [ro]) * distrib_(-1, 1, GD, GD, ?b) * GD(?a) );
    id GD(?a) * GD([mu]?) * GD(?b) * GD(?c) = GB([mu], ?a, ?b, ?c);
  endrepeat;
#endif
  id GB([mu]?, ?a) = GC([mu]) * GB(?a);
  id GB() = 1;
endwhile;

* Order the gamma matrices canonically
repeat;
  id GC([p1]?) * GC([p1]?) = [p1].[p1];
  disorder GC([mu]?) * GC([nu]?) = 2*d_([mu], [nu]) * TAG - GC([nu]) * GC([mu]);
endrepeat;
chainin GC;
id GA(?a) * GC(?b) = GA(?a, ?b);

id ifmatch->1 TAG = 1;
#endprocedure

***********************************************************************

#procedure ChainSimplify(contr)
#call ChainOrder(`contr')

#if "`OnShell'" == "True"

#do i = 1, `Legs'
#ifdef `k`i''
b `Fermionic';
.sort
keep brackets;

id k`i' = `k`i'';
#call ChainOrder(`contr')
#endif
#enddo

#endif
#endprocedure

***********************************************************************

#procedure FierzBefore
#if `Evanescent' == 1
b Spinor;
.sort
keep brackets;

id Spinor(?a) * GA(?g) * Spinor(?b) =
  Spinor(?a) * GA(?g) * Spinor(?b) *
  Evanescent(DiracChain(Spinor(?a), ?g, Spinor(?b)));
repeat id Evanescent([x]?) * Evanescent([y]?) = Evanescent([x] * [y]);
#endif
#endprocedure

***********************************************************************

#procedure FierzPre(ord)
b `Fermionic', ORD;
.sort
keep brackets;

id Spinor([p1]?, ?a) * GA(?g) * Spinor([p2]?, ?b) =
  CH(Spinor([p1], ?a), ?g, Spinor([p2], ?b))
#if `ord'
  * ORD([p1]) * ORD([p2]);

repeat id ORD([p1]?) * ORD(?a, [p1]?, ?b) = ORD(?a, [p1], ?b);

chainin ORD
#endif
  ;

* these relations are obtained by Fierzing twice

id CH(Spinor(?a), 6, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), 7, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), 6, Spinor(?b)) *
    CH(Spinor(?c), 7, Spinor(?d));

id CH(Spinor(?a), [om]?, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), [om]?, [ro]?, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), [om], Spinor(?b)) *
    CH(Spinor(?c), [om], [ro], Spinor(?d));
also CH(Spinor(?a), [omA]?, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), [omB]?, [ro]?, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), [omA], [ro], [mu], Spinor(?b)) *
    CH(Spinor(?c), [omB], [mu], Spinor(?d));

id CH(Spinor(?a), 6, [ro]?, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), 7, [si]?, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), 6, [ro], Spinor(?b)) *
    CH(Spinor(?c), 7, [si], Spinor(?d));

id CH(Spinor(?a), [om]?, [mu]?, [nu]?, [ro]?, Spinor(?b)) *
   CH(Spinor(?c), [om]?, [mu]?, [nu]?, [ro]?, Spinor(?d)) =
  16*CH(Spinor(?a), [om], [mu], Spinor(?b)) *
     CH(Spinor(?c), [om], [mu], Spinor(?d));
also CH(Spinor(?a), [omA]?, [mu]?, [nu]?, [ro]?, Spinor(?b)) *
   CH(Spinor(?c), [omB]?, [mu]?, [nu]?, [ro]?, Spinor(?d)) =
  4*CH(Spinor(?a), [omA], [mu], Spinor(?b)) *
    CH(Spinor(?c), [omB], [mu], Spinor(?d));
#endprocedure

***********************************************************************

#procedure FierzPost
id CH([x]?, ?g, [y]?) = [x] * GA(?g) * [y];

id D = 4;
#call ChainSimplify(1)

.sort

id D = 4;
#endprocedure

***********************************************************************
* The following general Fierz identity is from hep-ph/0412245.

#procedure FierzIdentity(lhs, rhs)
once ifnomatch->2 `lhs'
  CH(Spinor([p1]?, ?a), [omA]?, ?A, [s2]?) *
  CH([s1]?, [omB]?, ?B, Spinor([p2]?, ?b)) = `rhs' sum_(JJ, 1, 5,
    CHI([omA])*g_(1, ?A) *
    DUAL(JJ, N21_?, N22_?) *
    CHI([omB])*g_(1, ?B) * sum_(KK, 1, 5,
      DUAL(KK, N31_?, N32_?) *
      Spinor([p1], ?a) * BASIS(KK, 3, N31_?, N32_?) * Spinor([p2], ?b)) *
    [s1] * BASIS(JJ, 2, N21_?, N22_?) * [s2]);

trace4, 1;

id g6_([i]?) = 2*GA(6);
id g7_([i]?) = 2*GA(7);
repeat id GA(?g) * g_([i]?, [mu]?) = GA(?g, [mu]);

* Chisholm's identity backwards to get rid of all e_
repeat;
  once GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [ro]?, [LA]?) =
    sign_([om] + nargs_(?a)) * (
      GA([om], ?a, [mu], [nu], [ro], ?b) -
      d_([mu], [nu]) * GA([om], ?a, [ro], ?b) +
      d_([mu], [ro]) * GA([om], ?a, [nu], ?b) -
      d_([nu], [ro]) * GA([om], ?a, [mu], ?b) );
endrepeat;

renumber;

#call ChainOrder(1)
#endprocedure

***********************************************************************

#procedure FierzUnordered
#call FierzPre(0)

#do rep = 1, 1
#call FierzIdentity(,)

label 2;
if( count(CH,1) ) redefine rep "0";

b CH;
.sort
keep brackets;
#enddo

#call FierzPost
#endprocedure

***********************************************************************

#procedure FierzOrdered
#call FierzPre(1)

#do rep = 1, 1
id ifmatch->2 ORD([p1]?, [p2]?, ?p) *
     CH(Spinor([p1]?, ?m1), ?g, Spinor([p2]?, ?m2)) =
  ORD(?p) * CH(Spinor([p1], ?m1), ?g, Spinor([p2], ?m2));

* charge conjugation to get first spinor in front
* the rules for this are exactly as in HelicityME
id ORD([p1]?, [p3]?, ?p) *
   CH(Spinor([p2]?, [m2]?, [s2]?), [x]?, ?g,
      Spinor([p1]?, [m1]?, [s1]?)) =
  -ORD([p1], [p3], ?p) * sign_(nargs_(?g)) *
  CH(Spinor([p1], [m1], -[s1]),
     (sign_(nargs_(?g))*(2*[x] - 13) + 13)/2, reverse_(?g),
     Spinor([p2], [m2], -[s2]));

* charge conjugation to get second spinor in back
id ORD([p3]?, [p2]?, ?p) *
   CH(Spinor([p2]?, [m2]?, [s2]?), [x]?, ?g,
      Spinor([p1]?, [m1]?, [s1]?)) =
  -ORD([p3], [p2], ?p) * sign_(nargs_(?g)) *
  CH(Spinor([p1], [m1], -[s1]),
     (sign_(nargs_(?g))*(2*[x] - 13) + 13)/2, reverse_(?g),
     Spinor([p2], [m2], -[s2]));

* Fierz to get second spinor together with first
#call FierzIdentity(ORD([p1]?\, [p2]?\, ?p) *, ORD(?p) *)

id Spinor(?a) * GA(?g) * Spinor(?b) = CH(Spinor(?a), ?g, Spinor(?b));

label 2;

id ORD() = 1;
if( count(ORD,1) ) redefine rep "0";

b CH, ORD;
.sort
keep brackets;
#enddo

#call FierzPost
#endprocedure

***********************************************************************

#procedure DiracFinal

#if `Antisymmetrize' == 0

id GA([om]?, ?g) = CC([om]) * GD(?g);

#else

* introduce antisymmetrized Dirac chains
id GA([om]?, ?g) = CC([om]) *
  sum_(KK, 0, nargs_(?g), 2, distrib_(-1, KK, DD, GD, ?g));

id DD() = 1;
id DD([mu]?, [nu]?) = d_([mu], [nu]);
repeat;
  once DD(?a) = g_(1, ?a)/4;
  trace4, 1;
endrepeat;
id D = 4;

antisymm GD;

#endif

b CC, GD, Evanescent, `Fermionic';
.sort
keep brackets;

#if "`FermionChains'" == "VA"
id CC([om]?) = CC(1)/2 + sign_([om]) * CC(5)/2;
#endif

#if `Antisymmetrize' == 1
id CC([w]?) * GD([mu]?, [nu]?, ?r) = CC(-[w]) * GD([mu], [nu], ?r);
#endif

id Spinor(?a) * CC([w]?) * GD(?g) * Spinor(?b) =
  ABB(1, DiracChain(Spinor(?a), [w], ?g, Spinor(?b)), ?g);

id CC([w]?) * GD(?g) = ABB(1, DiracChain([w], ?g), ?g);

repeat id ABB(1, [x]?, ?g) * ABB(1, [y]?, ?h) =
  ABB(1, [x]*[y], ?g, ?h);

#if 1
id ABB(1, DiracChain(?a, -6, [mu]?, [nu]?, ?b) *
          DiracChain(?c, -7, [mu]?, [nu]?, ?d), ?g) = 0;
#endif

#if "`FermionOrder'" != "None" && `Evanescent' == 1
ab ABB;
.sort
collect TMP;
normalize TMP;
id Evanescent([x]?) * TMP([y]?) = Evanescent([x], [y]) + [y];
id TMP([y]?) = [y];
argument Evanescent;
  id ABB([x]?, [y]?, ?r) = [y];
  id DiracChain(?a) = DiracChain(?a) * TMP(?a);
  chainout TMP;
  id TMP([mu]?index_) = ORD([mu]);
  id ORD([mu]?fixed_) = 1;
  id TMP(?a) = 1;
  id ORD([mu])^[n]? = ORD([mu]);
  repeat;
    once ORD([mu]?) = replace_([mu], N100_?);
    renumber;
  endrepeat;
endargument;
id Evanescent([x]?, [x]?) = 0;
#endif

#call Abbreviate
#endprocedure

***********************************************************************

#procedure Abbreviate
.sort

#call DotSimplify

*----------------------------------------------------------------------

if( count(cutM,1) );

id q1.[p1]? = qfM(q1.[p1]);
id e_(q1, [p1]?, [p2]?, [p3]?) = qfM(e_(q1, [p1], [p2], [p3]));
id abbM(fermM(WeylChain(?a, q1, ?b))) = qfM(WeylChain(?a, q1, ?b));

endif;

.sort

*----------------------------------------------------------------------

id [p1]?.[p2]? = ABB(0, [p1].[p2], [p1], [p2]);

id e_([mu]?, [nu]?, [ro]?, [si]?) =
  ABB(0, Eps([mu], [nu], [ro], [si]), [mu], [nu], [ro], [si]);

id d_([mu]?, [nu]?) = ABB(0, MetricTensor([mu], [nu]), [mu], [nu]);

id [t]?(?a) = ABB(0, [t](?a), ?a);

id [p1]?([mu]?) = ABB(0, [p1]([mu]), [p1]);

repeat;
  once ABB([s1]?, [x]?, ?a, [mu]?!fixed_, ?b) *
       ABB([s2]?, [y]?, ?c, [mu]?, ?d) =
    ABB([s1] + [s2], [x]*[y], ?a, ?b, ?c, ?d) * replace_([mu], N100_?);
  also once ABB([s1]?, [x]?, ?a, [mu]?!fixed_, ?b, [mu]?, ?c) =
    ABB([s1], [x], ?a, ?b, ?c) * replace_([mu], N100_?);
  renumber;
endrepeat;

id ABB(0, [x]?, ?a) = abbM([x]);
id ABB([i]?, [x]?, ?a) = fermM([x]);

#if "`FermionChains'" != "Weyl"
repeat id fermM([x]?) * fermM([y]?) = fermM([x] * [y]);
argument fermM;
toPolynomial;
endargument;
id fermM([x]?) = Mat(fermM([x]));
argument Mat;
toPolynomial;
endargument;
#endif

argument abbM;
toPolynomial;
endargument;

b addM, mulM;
moduleoption polyfun=abbM;
.sort

b abbM;
.sort
on oldFactArg;
keep brackets;

id abbM([x]?) = abbM(nterms_([x]), [x]);
id abbM(1, [x]?) = TMP([x]);
also abbM([n]?, [x]?) = abbM([x]);

factarg abbM;
chainout abbM;
id TMP([x]?) = abbM([x]);

makeinteger abbM;
id abbM(1) = 1;

b abbM;
.sort
off oldFactArg;
keep brackets;

toPolynomial onlyfunctions abbM;
#endprocedure

***********************************************************************

#procedure CollectTerms
collect dotM;

moduleoption polyfun=dotM;
.sort

makeinteger dotM;
id dotM([x]?) = dotM(nterms_([x]), [x]);
id dotM(1, [x]?) = mulM([x]);
id dotM([n]?, [x]?) = mulM(dotM([x]));

argument mulM;
toPolynomial;
endargument;

moduleoption polyfun=mulM;
.sort
on oldFactArg;

#call Factor(mulM)

b mulM;
.sort
off oldFactArg;
keep brackets;

argument mulM;
toPolynomial;
endargument;

id mulM([x]?symbol_) = [x];

toPolynomial onlyfunctions mulM;
#endprocedure

***********************************************************************

#procedure MomReduce
id MOM(0) = 0;

#ifdef `MomSum'
* Apply momentum conservation to generate as few terms as possible

id MOM([p1]?) = MOM([p1], nterms_([p1]),
  [p1] + (`MomSum'), nterms_([p1] + (`MomSum')),
  [p1] - (`MomSum'), nterms_([p1] - (`MomSum')));

symm MOM (2,1) (4,3) (6,5);
#endif
#endprocedure

***********************************************************************

#procedure MomSquare
#call MomReduce
id MOM([p1]?, ?p) = MOM([p1].[p1]);
argument MOM;
#call kikj
endargument;
id MOM([x]?) = [x];
#endprocedure

***********************************************************************

#procedure Fewest(foo)
argument `foo';
#call Neglect
endargument;
id `foo'([x]?, [y]?) = `foo'([x], nterms_([x])*2 - 1, [y], nterms_([y])*2);
symm `foo' (2,1), (4,3);
id `foo'([x]?, ?a) = `foo'([x]);
#endprocedure

***********************************************************************

#procedure Factor(foo)
id `foo'(?x) = mulM(`foo'(?x));
argument mulM;
factarg `foo';
chainout `foo';
makeinteger `foo';
id `foo'([x]?) = `foo'(nterms_([x]), [x]);
id `foo'(1, [x]?) = [x];
id `foo'([n]?, [x]?) = `foo'([x]);
endargument;
makeinteger mulM;
#endprocedure

***********************************************************************

#procedure IndexSimplify
collect mulM;

repeat id powM([x]?, [y]?)^2 = powM([x], 2*[y]);
repeat id powM([x]?, [y]?) * powM([x]?, [z]?) = powM([x], [y] + [z]);
id powM([x]?, [y]?int_) = [x]^[y];

moduleoption polyfun=mulM;
.sort
#ifdef `Inserted'
#if `NoCostly' == 1
on oldFactArg;
#endif
#endif

argument mulM;
#call Neglect
endargument;

argument;
argument;
#call Square
endargument;
#call Square
endargument;

id mulM(0) = 0;
#call Factor(mulM)

.sort
off oldFactArg;

repeat;
  once IndexSum([x]?, [i]?, ?n) =
    TMP([x] * replace_([i], N100_?) * SumOver(N100_?, ?n, Renumber));
  renumber;
  id TMP([x]?) = [x];
endrepeat;

#ifdef `Inserted'
repeat;
  once SumOver([i]?, ?a, Renumber) =
    TMP(N100_?) * SumOver(N100_?, ?a) * replace_([i], N100_?);
  renumber;
endrepeat;

id IndexEps([i]?, [j]?, [k]?) = EPS([i], [j], [k]);

repeat;
  id EPS([I]?, [J]?, [K]?) * EPS([I]?, [J]?, [K]?) *
    SumOver([I]?, 3) * SumOver([J]?, 3) * SumOver([K]?, 3) = 6;
  id EPS([I]?, [J]?, [k]?) * EPS([I]?, [J]?, [c]?) *
    SumOver([I]?, 3) * SumOver([J]?, 3) = 2*IndexDelta([k], [c]);
  id EPS([I]?, [j]?, [k]?) * EPS([I]?, [b]?, [c]?) *
    SumOver([I]?, 3) =
    IndexDelta([j], [b])*IndexDelta([k], [c]) -
    IndexDelta([j], [c])*IndexDelta([k], [b]);
  repeat;
    id IndexDelta([I]?, [I]?) = 1;
    symm IndexDelta;
    once ifmatch->1 IndexDelta([i]?, [J]?) * SumOver([J]?, [x]?) =
      replace_([J], [i]);
    once IndexDelta([I]?, [j]?) * SumOver([I]?, [x]?) =
      replace_([I], [j]);
    label 1;
  endrepeat;
endrepeat;

id IndexDelta([x]?int_, [y]?int_) = delta_([x], [y]);

id TMP([x]?int_) = 1;
repeat id TMP([I]?)^2 = TMP([I]);

renumber 1;

#do i = 1, 9
once TMP([I]?) = replace_([I], Ind`i');
#enddo

id EPS([i]?, [j]?, [k]?) = IndexEps([j], [k], [i]);

argument mulM;
toPolynomial;
endargument;

toPolynomial onlyfunctions addM, powM, IndexDelta, IndexEps;
#endif
#endprocedure

***********************************************************************

#procedure DoInsertions
.sort
hide;

#call Insertions
#call Neglect

.sort
#endprocedure


***********************************************************************
*** main program starts here

#if "`InsertionPolicy'" == "Begin"
#call DoInsertions
#define Inserted
#else
#call Neglect
#endif

#call eiki

#call ConstBracket
.sort

*----------------------------------------------------------------------

#define SUNObjs "SUNSum, SUNT, SUNTSum, SUNF, SUNEps"

#define Fermionic "Spinor, GA, e_, `Tensors'"

* variables appearing in the CalcFeynAmp input and output
s I, Pi, D, Dminus4, `Invariants';
s Gamma5Test, Finite, MuTilde, MuTildeSq, Renumber;
s tnj, xnj, b0nj, b1nj, b2nj;
cf SumOver, PowerOf, Mat, Den, A0, IGram, List;
cf MetricTensor, Eps, DiracChain, WeylChain, Evanescent;
cf IndexDelta, IndexEps, IndexSum, `SUNObjs', SUNTr(c);
f Spinor, g5M, g6M, g7M;
i Col1,...,Col`Legs', Ind1,...,Ind9;
v nul, vTnj, v0nj, v1nj, v2nj, v3nj, v4nj;

* variables that make it into Mma but don't appear in the output
extrasymbols array subM;
cf addM, mulM, powM, dotM, abbM, fermM, sunM;
cf intM, extM, paveM, cutM, numM, qfM, qcM;
s dm4M, njM;

* patterns
s [x], [y], [z], [w], [n], [h];
s [k1], [k2], [k1k2];
s <[m0]>,...,<[m20]>, [mk];
s <[s0]>,...,<[s20]>;
v <[p0]>,...,<[p20]>, [pk], <[q1]>,...,<[q20]>;
i <[i0]>,...,<[i20]>;
i [mu], [nu], [ro], [si], [LA];
i [om], [omA], [omB];
i [i], [j], [k], [l], [I], [J], [K];
i [a], [b], [c], [d];
cf [f];
t [t];

* variables internal to FORM
s TAG, ETAG, QTAG, CUTRAT, JJ, KK;
cf TMP, MOM, ABB, SUNX, ORD, NOW, CH, SIGN(antisymm);
cf NEQ, NN, FF, DROP, JGRAM, D1, D2, E1, E2, HDEL;
t NUM, EQ, DD, EPS(antisymm);
nt GA, GB, GC, GD;
f GF, CC, WC;
auto s ARG;
set MOMS: k1,...,k`Legs';
set COLS: Col1,...,Col`Legs';
set INVS: `Invariants';

ntable BASIS(1:5, [i]?, [mu]?, [nu]?);
ntable DUAL(1:5, [mu]?, [nu]?);
ntable CHI(6:7);

*#define sig(i,mu,nu) "i_/2*(g_(`~i',`~mu',`~nu') - g_(`~i',`~nu',`~mu'))"
#define sig(i,mu,nu) "i_*(g_(`~i',`~mu',`~nu') - g_(`~i')*d_(`~mu',`~nu'))"

fill BASIS(1) = g6_([i])/2;
fill BASIS(2) = g7_([i])/2;
fill BASIS(3) = g_([i], 6_, [mu])/2;
fill BASIS(4) = g_([i], 7_, [mu])/2;
fill BASIS(5) = 1/2*(g6_([i]) + g7_([i]))/2*`sig([i], [mu], [nu])';

* DUAL includes the 1/2 to cancel the 2 = Tr BASIS(i) DUAL(i)
fill DUAL(1) = 1/2*g6_(1)/2;
fill DUAL(2) = 1/2*g7_(1)/2;
fill DUAL(3) = 1/2*g_(1, 7_, [mu])/2;
fill DUAL(4) = 1/2*g_(1, 6_, [mu])/2;
fill DUAL(5) = 1/2*1/2*`sig(1, [mu], [nu])';

fill CHI(6) = g6_(1)/2;
fill CHI(7) = g7_(1)/2;

*----------------------------------------------------------------------

#call IndexSimplify

ab `SUNObjs';
.sort

collect SUNX;
makeinteger SUNX;
id SUNX(1) = 1;

id Eps([mu]?, [nu]?, [ro]?, [si]?) = e_([mu], [nu], [ro], [si]);


if( count(dirM,1) ) ;

repeat;
  repeat id dirM([x]?, [i]?, [j]?) * dirM([y]?, [j]?, [k]?) =
    dirM([x]*[y], [i], [k]);
  id dirM([x]?, [i]?, [j]?) * dirM([y]?, [k]?, [j]?) =
    dirM([x]*CC([y]), [i], [k]);
  id dirM([x]?, [j]?, [i]?) * dirM([y]?, [j]?, [k]?) =
    dirM(CC([x])*[y], [i], [k]);
  argument dirM, 1;
    argument CC;
      id g_(sM, [mu]?) = -CC([mu]);
      chainin CC;
      id CC(?a) = CC(reverse_(?a));
      id g5_(sM) * CC(?a) = g_(sM, ?a, 5_);
      id g6_(sM) * CC(?a) = g_(sM, ?a, 6_);
      id g7_(sM) * CC(?a) = g_(sM, ?a, 7_);
      id CC(?a) = g_(sM, ?a);
    endargument;
    id CC([x]?) = [x];
  endargument;
endrepeat;

id dirM([x]?, [i]?) * dirM([y]?, [i]?, [j]?) *
     dirM(Spinor(?k, [s2]?)*gi_(sM), [j]?) =
  dirM([x]*[y]*Spinor(?k, -[s2]), [i], [j]);

$fline = 1;
id dirM([x]?, [i]?, [i]?) = -CH([x]);
while( count(CH,1) );
  once CH([x]?) = TMP([x]*replace_(sM, $fline));
  $fline = $fline + 1;
endwhile;

$fline = 100;
while( count(dirM,1) );
  once dirM([x]?, ?i) = TMP([x]*replace_(sM, $fline)) * ORD(?i);
  $fline = $fline + 1;
endwhile;

id TMP([x]?) = [x];

if( count(ORD,1) );
  redefine HaveFermions "1";
  chainin ORD;
  id ORD(?a) = SIGN(?a)*sign_(nargs_(?a)/2);
  mul replace_(SIGN, ORD);
  id ORD(?a) = 1;
endif;

endif;

id g6M([i]?{<100}) = gi_([i]) + g5M([i]);
id g7M([i]?{<100}) = gi_([i]) - g5M([i]);
repeat;
  once g5M([i]?{<100}) = e_(N100_?, N101_?, N102_?, N103_?) *
    g_([i], N100_?, N101_?, N102_?, N103_?)/24;
  renumber;
endrepeat;

mul replace_(g5M, g5_, g6M, g6_, g7M, g7_);

#do i = 1, 10
trace4, `i';
#enddo

#call eiki

moduleoption local $fline;
.sort

#call Contract

id intM(?a) = intM(?a) * NN(nargs_(?a));

*----------------------------------------------------------------------

#if `HaveFermions' == 1

id MuTilde * NN([n]?{<`OPP'}) = 0;
id MuTilde^2 = qfM(MuTildeSq);
id MuTilde = 0;

b g_, `Fermionic';
.sort
keep brackets;

id gi_([i]?) = g6_([i])/2 + g7_([i])/2;
id g5_([i]?) = g6_([i])/2 - g7_([i])/2;
id g_([i]?, [mu]?) = TMP([i]) * g_([i], [mu]);
repeat id TMP([i]?) * TMP([i]?) = TMP([i]);
id g_([i]?, [mu]?) * TMP([i]?) = g_([i], 6_, [mu])/2 + g_([i], 7_, [mu])/2;

id g6_([i]?) = 2*GA(6)
#if `Gamma5Test' == 1
  + 2*(GA(6) - GA(7)) * Gamma5Test * Dminus4
#endif
  ;
id g7_([i]?) = 2*GA(7)
#if `Gamma5Test' == 1
  - 2*(GA(6) - GA(7)) * Gamma5Test * Dminus4
#endif
  ;

repeat id GA(?g) * g_([i]?, [mu]?) = GA(?g, [mu]);

#call ChainSimplify(1)
#endif

*----------------------------------------------------------------------

#if `CancelQ2' == 1

b q1, intM, NN;
.sort
keep brackets;

id ifmatch->1 NN([n]?{>=`OPP'}) = 1;
id NN(?i) = 1;

* cancel q^2's in the numerator

repeat;
  once q1.q1 * intM(?a, Den(q1, [m1]?), ?b) =
    TAG * intM(?a, ?b) + [m1] * intM(?a, Den(q1, [m1]), ?b);
  once TAG * q1.q1 * intM(?a, Den([p1]?!{q1}, 0), ?b) =
    replace_(q1, 2*q1 - [p1]) * q1.q1 * intM(?a, Den([p1], 0), ?b);
  also once TAG * q1.q1 * intM(Den([p1]?!{q1}, [m1]?), ?a) =
    replace_(q1, 2*q1 - [p1]) * q1.q1 * intM(Den([p1], [m1]), ?a);
  id TAG = 1;
endrepeat;

id intM() = 0;

label 1;

#endif

*----------------------------------------------------------------------

.sort

id intM(Den([p1]?, 0)) = 0;

id Den([p1]?, [m1]?) * [p1]?.[p1]? = 1 + [m1]*Den([p1], [m1]);

*----------------------------------------------------------------------

b q1, intM;
.sort
keep brackets;

#if `OPP' <= `IntMax'
id intM(?d) = intM(nargs_(?d), ?d);
id intM([n]?{<`OPP'}, ?d) = intM(?d);
also intM([n]?, ?d) = cutM(?d)
#if "`OPPMethod'" == "AnaRat"
  + CUTRAT * intM(?d)
#endif
  ;
#endif

argument intM;
id Den([p1]?, ?m) = Den(?m)*MOM([p1]);
endargument;

#if `SortDen' == 1
symm intM;
#endif

once intM(Den(?m1)*MOM([p1]?)) = intM(Den(0,[p1],?m1)) *
  ORD(0) * replace_(q1, 2*q1 - [p1]);
#do n = 2, `IntMax'
also once intM(<Den(?m1)*MOM([p1]?)>,...,<Den(?m`n')*MOM([p`n']?)>) =
  intM(<Den(0,[p1],?m1)>*...*<Den({`n'-1},[p`n'],?m`n')>) *
    NN(`n') *
    ORD(<paveM(1)*([p2]-[p1])>+...+<paveM({`n'-1})*([p`n']-[p1])>) *
    replace_(q1, 2*q1 - [p1]);
#enddo

*----------------------------------------------------------------------

b q1, NUM, ORD, NN, intM, D, Dminus4, CUTRAT;
.sort
keep brackets;

if( count(ORD,1) ) totensor q1, NUM;

#if "`Dim'" == "4"
* add local terms for dimred/CDR as given in Appendix B of
* hep-ph/9806451 (note: 1/(16 Pi^2) already included in intM)

#if 0
id NUM([mu]?, [nu]?, [ro]?, [si]?) * NN(4) * intM([x]?) =
  NUM([mu], [nu], [ro], [si]) * NN(4) * intM([x]) -
  5/144 * NEQ([mu], [nu], [ro], [si]) * Finite +
  1/8 * distrib_(1, 2, EQ, NEQ, [mu], [nu], [ro], [si]) * Finite;
#endif

#if 0
also NUM([mu]?, [nu]?, [ro]?) * NN(3) *
       intM(<Den(0,[p0]?,[m0]?)>*...*<Den(2,[p2]?,[m2]?)>) =
  NUM([mu], [nu], [ro]) * NN(3) *
    intM(<Den(0,[p0],[m0])>*...*<Den(2,[p2],[m2])>) +
  1/36 * NEQ([mu], [nu], [ro], [p2] - [p0]) * Finite;
#endif

#if 0
also NUM([mu]?, [nu]?) * NN(3) * intM([x]?) =
  NUM([mu], [nu]) * NN(3) * intM([x]) -
  1/8 * NEQ([mu], [nu]) * Finite;
#endif

#if 0
id EQ([mu]?, [mu]?) = 1;
id EQ(?a) = 0;

symm NEQ;
id NEQ(?a, [mu]?, [mu]?, ?b) = 0;
id NEQ(?a) = dd_(?a);
#endif

id CUTRAT * intM(?a) = 0;
id CUTRAT = 1;
#endif


* decompose into Lorentz-covariant tensors

* The following statement introduces the g_{\mu\nu}'s in a smart way.
* Lifted from: S.A. Larin, T. van Ritbergen, and J.A.M. Vermaseren,
* The optimization of a huge FORM program,
* in: Proceedings Oberammergau 1993, ISBN 9-810-21699-8.

id NUM(?i) = sum_(KK, 0, nargs_(?i), 2,
  paveM(0)^KK * distrib_(1, KK, dd_, NUM, ?i));

id ORD(0) * NUM([mu]?, ?i) = 0;
repeat id ORD([p1]?) * NUM([mu]?, ?i) = ORD([p1]) * d_([p1], [mu]) * NUM(?i);

#call eiki

id ORD(?p) = 1;
id NUM() = 1;

chainin paveM;

*----------------------------------------------------------------------

id D = Dminus4 + 4;

#if "`PaVeReduce'" != "False"

id NN([i]?) * Dminus4 = Dminus4;

#do rep = 1, 1

b NN, paveM, intM;
.sort
keep brackets;

id ifnomatch->1 NN([i]?) * intM([x]?) = NN([i]) * intM([x]) * [x];

* symmetrize the coefficients for N > 4
* hep-ph/0509141 Eq. (6.14+15)
if( match(NN([n]?{>4})) );
  id paveM(0,0,[i1]?,?i) = paveM([i1],0,0,?i);
  id paveM([i]?,[j]?,?i) = TMP([i],[j],TAG,?i) + paveM([j],[i],?i);
  repeat id TMP([i]?,?i,TAG,[j]?,?j) =
    TMP([i],?i,[j],TAG,?j) + paveM([j],?i,[i],?j);
  id TMP(?i,TAG) = paveM(?i);
endif;

* hep-ph/0509141 Eq. (7.13)
id NN(6) * paveM([i1]?,?i) = NN(6) *
  deltap_([i1], 0) * sum_(KK, 1, 5, IGram([i1],KK) * DROP(KK, ?i));

* hep-ph/0509141 Eq. (6.13)
also NN(5) * paveM(0,0,?i) = NN(5) * (
  Dminus4 * paveM(0,0,?i) +
  sum_(KK, 1, 4, JGRAM(KK,0) * DROP(KK, 0,0,?i)) );

* hep-ph/0509141 Eq. (6.12)
also NN(5) * paveM([i1]?,?i) = NN(5) * (
  Dminus4 * paveM([i1],?i) +
  sum_(KK, 0, 4, JGRAM([i1],KK) * DROP(KK, ?i)) -
  2*sum_(KK, 1, 4, E1([i1], KK) * distrib_(1, 1, E2, paveM, ?i)) );

* hep-ph/0509141 Eq. (5.10)
also NN([n]?) * paveM(0,0,?i) = NN([n])/(3 + nargs_(0,0,?i) - [n]) * (
  -Dminus4 * paveM(0,0,?i) +
  FF(0) * paveM(?i) +
  sum_(KK, 1, [n] - 1, FF(KK) * paveM(KK,?i))/2 -
  DROP(0, ?i)/2 );

* hep-ph/0509141 Eq. (5.11+8)
also NN([n]?) * paveM([i1]?,?i) = NN([n]) *
  sum_(KK, 1, [n] - 1, IGram([i1],KK) * (
    DROP(KK, ?i) -
    FF(KK) * paveM(?i) -
    2*D1(KK) * distrib_(1, 1, D2, paveM, ?i) ));

id D1([k]?) * D2([i2]?) * paveM(?i) =
  delta_([k], [i2]) * paveM(0,0,?i);
id E1([i1]?, [k]?) * E2([i2]?) * paveM(?i) =
  JGRAM([i1],[k], 0,[i2]) * DROP([k], 0,0,?i);

* hep-ph/0509141 Eq. (2.28)
id NN([n]?) * JGRAM([s1]?,0) = -NN([n]) *
  sum_(KK, 1, [n], IGram([s1],KK) * FF(KK));
id NN([n]?) * JGRAM([s1]?,[s2]?) = NN([n]) * (
  2*FF(0) * IGram([s1],[s2]) +
  sum_(JJ, 1, [n], sum_(KK, 1, [n],
    IGram([s1],KK, [s2],JJ) * FF(KK) * FF(JJ))) );
* hep-ph/0509141 Eq. (2.29)
id NN([n]?) * JGRAM([s1]?,[s2]?, 0,[s4]?) =
  -NN([n]) * sum_(KK, 1, [n], IGram([s1],[s2], KK,[s4]));

id FF(0) * Den(0, [p0]?, [m0]?) = [m0] * Den(0, [p0], [m0]);
id FF([k]?) * Den(0, [p0]?, [m0]?) * Den([k]?, [pk]?, [mk]?) =
  (MOM([pk] - [p0]) - [mk] + [m0]) *
  Den(0, [p0], [m0]) * Den([k], [pk], [mk]);

id IGram([s1]?,[s2]?) * intM([x]?) =
  sign_([s1] + [s2]) *
  IGram(1, DROP([s1]) * [x]) *
  IGram(1, DROP([s2]) * [x]) *
  IGram(2, [x]);
also IGram([s1]?,[s2]?, [s3]?,[s4]?) * intM([x]?) =
  sign_([s1] + [s2] + [s3] + [s4]) *
  sig_([s1] - [s3]) * IGram(1, DROP([s1]) * DROP([s3]) * [x]) *
  sig_([s4] - [s2]) * IGram(1, DROP([s2]) * DROP([s4]) * [x]) *
  IGram(2, [x]);
id intM(?x) = 1;

argument IGram;
id DROP([k]?) * Den([k]?, ?p) = 1;
endargument;

id IGram(?i, Den(?p)) = IGram(?i);
#do n = 2, `IntMax'
id IGram(?i, <Den([i1]?,[p1]?,[m1]?)>*...*<Den([i`n']?,[p`n']?,[m`n']?)>) =
  IGram(?i, <[p2]-[p1]>,...,<[p`n']-[p1]>);
#enddo
id IGram(1, ?n1) * IGram(1, ?n2) * IGram(2, ?d) =
  IGram(MOM(?n1) * MOM(?n2), MOM(?d)^2)/2;

#call MomSquare

argument IGram;
id MOM() = 1;
repeat id MOM(?a, 0, ?b) = MOM(?a, nul, ?b);
#do n = 2, `IntMax'
id MOM(<[p2]?>,...,<[p`n']?>) = e_(<[p2]>,...,<[p`n']>);
#enddo
contract;
endargument;

id IGram([x]?, [p1]?.[p1]?) * [p1]?.[p1]? = [x];

#if "`PaVeReduce'" == "True"
id IGram([x]?, [y]?) = IGram(TMP([x]), TMP([y]));
argument IGram;
argument TMP;
#call kikj
endargument;
#call InvSimplify(TMP)
id TMP([x]?) = [x];
endargument;

id IGram([x]?, [y]?) = [x] * IGram([y]);
factarg IGram;
chainout IGram;
id IGram(0) = IGram(0);
also IGram([x]?number_) = 1/[x];
also IGram([x]?symbol_) = 1/[x];
#endif

id NN([n]?) * DROP([k]?, ?i) = NN([n] - 1) * paveM() *
  (deltap_([k], 0) * DROP([k], ?i) - DROP(0, ?i));

* hep-ph/0509141 Eq. (2.8)
repeat;
  id paveM(?n) * DROP(0, 1,?i) * NN([n]?) =
    (-paveM(?n) - sum_(KK, 1, [n] - 1, paveM(?n, KK))) *
    DROP(0, ?i) * NN([n]);
  also paveM(?n) * DROP(0, [h]?,?i) =
    paveM(?n, [h] - theta_([h] - 1)) * DROP(0, ?i);
endrepeat;

repeat id paveM(?n) * DROP([n]?, [h]?,?i) =
  deltap_([n], [h]) * paveM(?n, [h] - theta_([h] - [n])) *
  DROP([n], ?i);

id DROP([k]?) * Den([k]?, ?q) = 1;

#do n = 1, `IntMax'
id NN(`n') * <Den([i1]?,?p1)>*...*<Den([i`n']?,?p`n')> =
  NN(`n') * intM(<Den(0,?p1)>*...*<Den({`n'-1},?p`n')>);
#enddo

id NN(1) = 1;
id NN([i]?) * Dminus4 = Dminus4;
id paveM() = 1;

if( count(paveM,1, NN,1) == 2 ) redefine rep "0";

label 1;

.sort

#enddo

#endif

*----------------------------------------------------------------------

id NN(?i) = 1;

b intM, Den;
.sort
keep brackets;

id intM(Den(0, [p1]?, [m1]?)) = A0([m1]);
#do n = 2, `IntMax'
also intM(<Den(0,[p1]?,[m1]?)>*...*<Den({`n'-1},[p`n']?,[m`n']?)>) =
  intM(`n',
#do i = 1, {`n'/2}
    <MOM([p{`i'+1}]-[p1])>,...,<MOM([p`n']-[p{`n'-`i'}])>,
#if {2*`i'} < `n'
    <MOM([p1]-[p{`n'-`i'+1}])>,...,<MOM([p`i']-[p`n'])>,
#endif
#enddo
    <[m1]>,...,<[m`n']>);
#enddo

id Den([p1]?, ?m) = Den(MOM([p1]), ?m);

argument intM, Den;
#call MomSquare
endargument;

b A0, paveM, intM;
.sort
keep brackets;

id paveM(?i) * A0([m1]?) = theta_(sign_(nargs_(?i))) *
  ([m1]/2)^(nargs_(?i)/2)/fac_(nargs_(?i)/2 + 1) *
  (A0([m1]) + [m1]*sum_(KK, 2, nargs_(?i)/2 + 1, 1/KK));

id intM([n]?, ?r) * paveM(?i) = paveM([n], List(?i), ?r);
symm intM:4 3, 4;
id intM([n]?, ?r) = paveM([n], List(0), ?r);

*----------------------------------------------------------------------

.sort

#call Contract

#call DotSimplify

*----------------------------------------------------------------------

#if `HaveFermions' == 1
* Dirac algebra on open fermion chains again

#if "`FermionChains'" == "Weyl"

* Chisholm's identity backwards to get rid of all Eps
*repeat id GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [ro]?, [LA]?) =
*  1/4 * sign_([om]) * sign_(nargs_(?a)) * (
*    GA([om], ?a, [mu], [nu], [ro], ?b) -
*    GA([om], ?a, [ro], [nu], [mu], ?b) );
repeat;
  once GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [ro]?, [LA]?) =
    sign_([om] + nargs_(?a)) * (
      GA([om], ?a, [mu], [nu], [ro], ?b) -
      d_([mu], [nu]) * GA([om], ?a, [ro], ?b) +
      d_([mu], [ro]) * GA([om], ?a, [nu], ?b) -
      d_([nu], [ro]) * GA([om], ?a, [mu], ?b) );
endrepeat;

*#elseif "`Dim'" == "4"
#elseif 0

* this is Chisholm's identity:
repeat;
  once GA([om]?, [mu]?, [nu]?, [ro]?, ?a) =
    sign_([om]) * GA([om], N100_?, ?a) * e_([mu], [nu], [ro], N100_?) +
    d_([mu], [nu]) * GA([om], [ro], ?a) -
    d_([mu], [ro]) * GA([om], [nu], ?a) +
    d_([nu], [ro]) * GA([om], [mu], ?a);
  renumber;
endrepeat;

#call Contract

#endif

b `Fermionic';
.sort
keep brackets;

#call ChainSimplify(1)

#endif

*----------------------------------------------------------------------

#if "`Dim'" != 4
b D, Dminus4, A0, paveM, cutM, CUTRAT;
.sort
keep brackets;

id D = Dminus4 + 4;

#if "`Dim'" == "D"

#if "`OPPMethod'" == "AnaRat"
id Dminus4 * cutM(?d) = 0;
#else
id Dminus4 * cutM(?d) = dm4M * cutM(?d);
#endif

* add local terms for dimreg
also Dminus4 * A0([m1]?) = -2*[m1]*Finite;
also Dminus4 * paveM(2, List(0), ?a) = -2*Finite;
also Dminus4 * paveM(2, List(1), ?a) = Finite;
also Dminus4 * paveM(2, List(0,0), [k1]?, [m1]?, [m2]?) =
  1/6*([k1] - 3*[m1] - 3*[m2])*Finite;
also Dminus4 * paveM(2, List(1,1), ?a) = -2/3*Finite;
also Dminus4 * paveM(2, List(0,0,1), [k1]?, [m1]?, [m2]?) =
  -1/12*([k1] - 2*[m1] - 4*[m2])*Finite;
also Dminus4 * paveM(2, List(1,1,1), ?a) = 1/2*Finite;
also Dminus4 * paveM(3, List(0,0), ?a) = -1/2*Finite;
also Dminus4 * paveM(3, List(0,0,[i]?), ?a) = 1/6*Finite;
also Dminus4 * paveM(3, List(0,0,0,0), [k1]?, [k2]?, [k1k2]?, [m1]?, [m2]?, [m3]?) =
  1/48*([k1] + [k2] + [k1k2] - 4*([m1] + [m2] + [m3]))*Finite;
also Dminus4 * paveM(3, List(0,0,[i]?,[i]?), ?a) = -1/12*Finite;
also Dminus4 * paveM(3, List(0,0,[i]?,[j]?), ?a) = -1/24*Finite;
also Dminus4 * paveM(4, List(0,0,0,0), ?a) = -1/12*Finite;
also Dminus4 * paveM(4, List(0,0,0,0,[i]?), ?a) = 1/48*Finite;

also Dminus4 = 0;

#redefine Dim "4"

#endif

id CUTRAT * paveM(?a) = 0;
id CUTRAT = 1;

#endif

*----------------------------------------------------------------------

#if `HaveFermions' == 0

#call Abbreviate

#else

#if "`FermionChains'" == "Weyl"

if( count(cutM,1) );
* HDEL arguments shall vanish if the external helicities
* select the 'small' components of a massless spinor:
* Spinor[+1] = {Spinor[+1, 6] propto Sqrt[p0 + p],
*               Spinor[+1, 7] propto Sqrt[p0 - p] -> 0},
* Spinor[-1] = {Spinor[-1, 6] propto Sqrt[p0 - p] -> 0,
*               Spinor[-1, 7] propto Sqrt[p0 + p]},
* i.e. zero if mismatch in 'oddity' ([+1, 7], [-1, 6])

* [s1]*HEL([n]) != sign_([om] - 1)
* [s1]*HEL([n]) == sign_([om])
* [s1]*HEL([n]) - sign_([om]) == 0
* HEL([n]) - [s1]*sign_([om]) == 0
id Spinor([p1]?MOMS[[n]], 0, [s1]?) * GA([om]?, ?g) =
  HDEL(mulM(HEL([n]) - [s1]*sign_([om]))) *
  Spinor([p1], 0, [s1]) * GA([om], ?g);

* [s1]*HEL([n]) != sign_([om] + nargs_(?g))
* [s1]*HEL([n]) == -sign_([om] + nargs_(?g))
* [s1]*HEL([n]) + sign_([om] + nargs_(?g)) == 0
* HEL([n]) + [s1]*sign_([om] + nargs_(?g)) == 0
id GA([om]?, ?g) * Spinor([p1]?MOMS[[n]], 0, [s1]?) =
  HDEL(mulM(HEL([n]) + [s1]*sign_([om] + nargs_(?g)))) *
  GA([om], ?g) * Spinor([p1], 0, [s1]);

id HDEL(mulM([x]?number_)) = delta_([x]);
endif;

b Spinor;
.sort
keep brackets;

id Spinor(?a) * GA(?g) * Spinor(?b) =
  CH(Spinor(?a, 2, 0), ?g, Spinor(?b, 1, 0));

repeat;
  once CH([s1]?, [x]?, ?a, [LA]?, ?b) *
       CH([s2]?, [y]?, ?c, [LA]?, ?d) =
    WC(sign_(nargs_(?a, ?c) + [x] + [y]), [s1], [x], ?a) *
    WC(?b) * WC([s2], [y], ?c) * WC(?d);

* Fierz 1: <A|sig_mu|B> <C|sigbar^mu|D> = 2 <A|D> <C|B>
  id WC(-1, ?a) * WC(?b) * WC(?c) * WC(?d) =
    2 * CH(?a, ?d) * CH(?c, ?b);

* Fierz 2: <A|sig(bar)_mu|B> <C|sig(bar)^mu|D> = 2 <A|eps|C> <B|eps|D>
  also WC(1, ?a) * WC(?b, Spinor(?s1, [s1]?)) *
         WC(Spinor(?s2, [s2]?), [x]?, ?c) * WC(?d) =
    2 * CH(?a, reverse_(?c), Spinor(?s2, 1 - [s2])) *
      CH(Spinor(?s1, 1 - [s1]),
        7 - mod_([x] + nargs_(?b, ?c), 2),
        reverse_(?b), ?d);

* due to the canonical ordering of the Dirac chains this
* is the only(?) case we need of Fierz on the same chain:
  repeat id CH(?a, [LA]?, [LA]?, ?b) = 4*CH(?a, ?b);
endrepeat;

id CH([s1]?, ?g, [s2]?) = [s1] * GA(?g) * [s2];

#call ChainSimplify(0)
*id CH(?a, [p1]?, [p1]?, ?b) = [p1].[p1] * CH(?a, ?b);

#call eiki
#call kikj

id Spinor(?s1) * GA(?g) * Spinor(?s2) =
  abbM(fermM(WeylChain(Spinor(?s1), ?g, Spinor(?s2))));

#call Abbreviate

#else

#switch "`FermionOrder'"

#case "None"
.sort
#call DiracFinal
#break

#case "Fierz"
* Fierz twice for simplification
#call FierzBefore
#call FierzUnordered
#call FierzUnordered
#call DiracFinal
#break

#case "Automatic"
* lexicographical ordering
#call FierzBefore
#call FierzOrdered
#call DiracFinal
#break

#case "Colour"
* postponed until after SUNT simplification
#break

#default
#do i = {`FermionOrder'}
#ifndef `order'
#define order "k{`i'}"
#else
#redefine order "`order',k{`i'}"
#endif
#enddo
mul ORD(`order');
#call FierzBefore
#call FierzOrdered
#call DiracFinal

#endswitch

#endif

#endif

*----------------------------------------------------------------------

#ifdef `Inserted'
.sort
#else
#call DoInsertions
#endif

id A0(0) = 0;

id Den(0, [x]?) = -Den([x], 0);
id Den([x]?, 0) * [x]? = 1;

#ifndef `Inserted'
#define Inserted

#call ConstBracket
.sort

#call IndexSimplify
.sort
#endif

*----------------------------------------------------------------------

#if `HaveSUN' == 1
* simplification of SU(N) structures

* The algorithm implemented here is an extension of the one given in
* J.A.M. Vermaseren, The use of computer algebra in QCD,
* in: Proceedings Schladming 1996, ISBN 3-540-62478-3.

* The idea is to transform all SU(N) objects to generators, SUNT.
* In the output, only two types of objects can appear:
* - chains of SUNTs (with external colour indices), or
* - traces of SUNTs.
* A chain of SUNTs is denoted by SUNT(a, b, ..., i, j), where
* a, b, ... are gluon indices and i and j are colour indices.
* SUNT(i, j) is the special case of the identity in colour space.
* A trace over SUNTs is marked by both colour indices being zero,
* i.e. SUNT(a, b, ..., 0, 0).

id SUNX([x]?) = [x];

b `SUNObjs';
.sort
keep brackets;

if( count(SUNF,1) );

  repeat;
    once SUNF(?a, [a]?, [b]?, [c]?, [d]?) =
      SUNF(?a, [a], [b], N100_?) * SUNF(N100_?, [c], [d]) * SUNSum(N100_?);
    renumber;
  endrepeat;

* f^{abc} = 2 i Tr(T^c T^b T^a - T^a T^b T^c)

  id SUNF([a]?, [b]?, [c]?) =
    2*i_*(SUNT([c], [b], [a], 0, 0) - SUNT([a], [b], [c], 0, 0));

endif;


repeat;
  once SUNT(?a, 0, 0) = SUNT(?a, N100_?, N100_?) * SUNSum(N100_?);
  renumber;
endrepeat;

repeat;
  once SUNT(?a, [a]?, [b]?, [i]?, [j]?) =
    SUNT(?a, [a], [i], N100_?) * SUNT([b], N100_?, [j]) * SUNSum(N100_?);
  renumber;
endrepeat;


* T^a_{ij} T^a_{kl} =
*   1/2 (delta_{il} delta_{jk} - 1/N delta_{ij} delta_{kl})

id SUNT([a]?, [i]?, [j]?) * SUNT([a]?, [k]?, [l]?) * SUNSum([a]?, ?a) =
  1/2 * SUNT([i], [l]) * SUNT([j], [k]) -
  1/2/(`SUNN') * SUNT([i], [j]) * SUNT([k], [l]);

id SUNTSum([i]?, [j]?, [k]?, [l]?) =
  1/2 * SUNT([i], [l]) * SUNT([j], [k]) -
  1/2/(`SUNN') * SUNT([i], [j]) * SUNT([k], [l]);


id SUNEps([i]?, [j]?, [k]?) = EPS([i], [j], [k]);


* cleaning up, step 1: get rid of the deltas

repeat;
  id EPS([I]?, [j]?, [k]?) * EPS([I]?, [b]?, [c]?) *
    SUNSum([I]?, [x]?) = 
    SUNT([j], [b])*SUNT([k], [c]) -
    SUNT([j], [c])*SUNT([k], [b]);
  repeat;
    id SUNT([I]?, [I]?) * SUNSum([I]?, ?a) = `SUNN';
    symm SUNT:2 1, 2;
    once ifmatch->1 SUNT([I]?, [j]?) * SUNSum([I]?, ?a) = replace_([I], [j]);
    once SUNT([i]?, [J]?) * SUNSum([J]?, ?a) = replace_([J], [i]);
    label 1;
  endrepeat;
endrepeat;

id SUNT([x]?int_, [y]?int_) = delta_([x], [y]);
id SUNT([a]?, [i]?, [i]?) * SUNSum([i]?, ?a) = 0;

id EPS([i]?, [j]?, [k]?) = sunM(SUNEps([j], [k], [i]));

* cleaning up, step 2: bead up the SUNTs into chains

repeat;
  once SUNT(?a, [a]?, [i]?, [j]?) = TMP(?a, [a], [i], [j]);
  repeat;
    id TMP(?a, [i]?, [j]?) * SUNT(?b, [j]?, [k]?) * SUNSum([j]?, ?c) =
      TMP(?a, ?b, [i], [k]);
    id SUNT(?a, [i]?, [j]?) * TMP(?b, [j]?, [k]?) * SUNSum([j]?, ?c) =
      TMP(?a, ?b, [i], [k]);
  endrepeat;

  id TMP(?a, [i]?, [i]?) * SUNSum([i]?, ?b) = SUNTr(?a);
  id SUNTr(?a) = TMP(?a, 0, 0);

* special case of Tr(T^a T^b) = 1/2 delta_{ab}
*  id TMP([a]?, [a]?, 0, 0) = 1/2;
  id TMP([x]?int_, [y]?int_, 0, 0) = 1/2*delta_([x], [y]);

  id TMP(?a) = sunM(SUNT(?a));
endrepeat;

symm SUNT;

id SUNT(?a) = sunM(SUNT(?a));

#if "`FermionOrder'" == "Colour"
id sunM(SUNT(?a, [i]?COLS[[x]], [j]?COLS[[y]])) =
  sunM(SUNT(?a, [i], [j])) * ORD(MOMS[[x]], MOMS[[y]]);
#endif

repeat id sunM([x]?) * sunM([y]?) = sunM([x] * [y]);

id SUNSum([i]?, [x]?) = [x];

* the Mat(...) are kept at the almost outermost level (only SumOver
* comes before), i.e. the amplitude is of the form Sum[c[i] Mat[i], i];
* this is required for the calculation of the squared amplitude

id sunM([x]?) = Mat(sunM([x]));

argument Mat;
toPolynomial;
endargument;

#if "`FermionOrder'" == "Colour"
#call FierzBefore
#call FierzOrdered
#call DiracFinal
#endif

.sort

#endif

*----------------------------------------------------------------------

b SumOver, PowerOf, Mat, Den, IGram, HDEL,
  mulM, A0, paveM, cutM, extM, qfM, Dminus4, dm4M;
.sort

#call CollectTerms

*----------------------------------------------------------------------

#if `OPP' > `IntMax'

.sort

repeat id Mat([x]?) * Mat([y]?) = Mat([x] * [y]);

#else

b cutM, extM, A0, paveM, qfM, Dminus4, dm4M,
  SumOver, PowerOf, Mat, HDEL, Den;
.sort

collect mulM;
makeinteger mulM;

b cutM, extM, A0, paveM, qfM, Dminus4, dm4M,
  SumOver, PowerOf, Mat, HDEL;
.sort

collect qcM;

*if( count(cutM,1) );

makeinteger qcM;
argument qcM;
toPolynomial;
endargument;

*endif;

argument qfM;
id WeylChain([s1]?, ?g, [s2]?) = [s1] * GA(?g) * [s2];
endargument;

moduleoption polyfun=qfM;
.sort

makeinteger qfM;

moduleoption polyfun=Mat;
.sort

makeinteger Mat;
id Mat(1) = 1;

*----------------------------------------------------------------------
* extend the denominators of the cutM

#if "`CombineDen'" != "False"

b cutM, extM, Mat;
.sort

collect TMP, TMP;

id cutM(?o) * extM(intM(?o), intM(?n), ?e) = cutM(?n) * extM(?e)
#if "`CombineDen'" == "Automatic"
  * intM(?o) * TAG^nargs_(?o) * ETAG^(nargs_(?e) - 1)
#endif
  ;
also extM(?a) = 1;

#if "`CombineDen'" == "Automatic"

b cutM, TAG, Mat;
.sort

collect TMP, TMP;

#define cost(n, r) "1 + \
  2*min_(`~n', `~r')^3 + \
  max_(`~n' - max_(`~r', 2), 0)^4/4"

$min = 1;
$max = 0;
$sum = 0;

id TAG^[n]?$n = 1;

if( $n > 0 );

argument TMP;
$q = 0;
argument TMP;
argument qfM;
id q1 = q1 * QTAG;
id QTAG^[n]?$e = 1;
$q = max_($q, $e);
endargument;
endargument;

id ETAG^[n]?$e = 1;

$min = $min*$e;
$max = max_($max, `cost($n + $e, $q + 2*$e)');
$sum = $sum + `cost($n, $q)';
endargument;

endif;

id TMP([x]?) = [x] * ORD($min, $max - $sum);

id ORD(0, [x]?neg0_) = 1;
also intM(?o) * cutM(?n) * extM(?e) = cutM(?o);
id ORD(?x) = 1;

#endif

id intM(?o) = 1;
id extM(?e, [p1]?) * TMP([x]?) = extM(?e) * TMP([x] * replace_(q1, [p1]));
id TMP([x]?) = [x];

b qfM, cutM, extM;
#if "`CombineDen'" == "Automatic"
moduleoption local $n, $q, $e, $max, $sum;
#endif
.sort

#endif

argument cutM;   
id Den([p1]?, ?m1) = Den(?m1)*MOM([p1]);
endargument;

#if `SortDen' == 1
symm cutM;
#endif

#do n = 2, `ExtMax'
also once cutM(<Den(?m1)*MOM([p1]?)>,...,<Den(?m`n')*MOM([p`n']?)>) =
  replace_(q1, 2*q1 - [p1]) *
  cutM(`n', List(
#do i = 1, {`n'/2}
      <MOM([p{`i'+1}]-[p1])>,...,<MOM([p`n']-[p{`n'-`i'}])>,
#if {2*`i'} < `n'
      <MOM([p1]-[p{`n'-`i'+1}])>,...,<MOM([p`i']-[p`n'])>,
#endif
#enddo
      <?m1>,...,<?m`n'>),
    <MOM([p2]-[p1])>,...,<MOM([p`n']-[p1])>, <?m1>,...,<?m`n'>);
#enddo

argument cutM;
argument List;
#call MomSquare
endargument;
#call MomReduce
id MOM([p1]?, ?p) = [p1];
endargument;

argument extM;
id Den([p1]?, [m1]?) = MOM([p1]) - [m1];
#call MomSquare
endargument;

b qfM;
.sort

collect ABB;
makeinteger ABB;
id qfM([x]?) = [x];

#call eiki
#call ChainOrder(0)
#call eiki

#call DotSimplify

id MuTildeSq = qfM(MuTildeSq);
id q1.[p1]? = qfM(q1.[p1]);
id e_(q1, [p1]?, [p2]?, [p3]?) = qfM(e_(q1, [p1], [p2], [p3]));
id Spinor(?s1) * GA(?g) * Spinor(?s2) =
  abbM(fermM(WeylChain(Spinor(?s1), ?g, Spinor(?s2))));
id abbM(fermM(WeylChain(?a, q1, ?b))) = qfM(WeylChain(?a, q1, ?b));

argument abbM;
toPolynomial;
endargument;

b ABB, qfM;
.sort

collect qcM;

moduleoption polyfun=qfM;
.sort

factarg qfM;
chainout qfM;
makeinteger qfM;
id qfM(1) = 1;

id ABB([x]?) = [x];

moduleoption polyfun=qcM;
.sort

on oldfactarg;

factarg qcM;
chainout qcM;
makeinteger qcM;
id qcM(1) = 1;
chainin qcM;

id extM(?a) = qfM(?a);
chainout qfM;

#call InvSimplify(qfM)

moduleoption polyfun=HDEL;
.sort

off oldfactarg;

makeinteger HDEL;

id cutM([n]?, ?a) * HDEL([x]?) = cutM([n], [x], ?a);
also cutM([n]?, ?a) = cutM([n], 1, ?a);

id HDEL(?h) = 1;

*----------------------------------------------------------------------

if( count(cutM,1) );

id SumOver([i]?, ?a) = SumOver([i], ?a) * ORD([i]);
chainin ORD;
id qcM(?a) * ORD(?i) = qcM(?a, List(?i)) * ORD(?i);

else;

id qcM(?a) = mulM(?a);

endif;

splitarg Mat;
repeat id Mat(?a, [x]?, [y]?) = Mat(?a, [x]) + Mat([y]);
makeinteger Mat;

b cutM, A0, paveM, SumOver, PowerOf, Mat, ORD;
.sort

collect numM, numM;

if( count(cutM,1) );

makeinteger numM;

#if "`OPPMethod'" == "Ninja"
argument numM;
repeat id qfM([x]?) * qfM([y]?) = qfM([x] * [y]);
argument qfM;
id WeylChain([s1]?, ?g, [s2]?) = [s1] * GA(?g) * [s2];

id MuTildeSq = TMP(0) * MuTildeSq +
  TMP(1) * tnj^2 * vTnj.vTnj;
id q1 = TMP(0) * q1 +
  TMP(1) * tnj * vTnj +
  TMP(2) * (v0nj + tnj * v3nj + (b0nj + MuTildeSq)/tnj * v4nj) +
  TMP(3) * (v1nj + xnj * v2nj + tnj * v3nj +
              (b0nj + b1nj*xnj + b2nj*xnj^2 + MuTildeSq)/tnj * v4nj);

id v2nj.v3nj = 0;
id v2nj.v4nj = 0;
id v3nj.v3nj = 0;
id v4nj.v4nj = 0;
id v3nj.v4nj = 1/2;

repeat id TMP([i]?) * TMP([i]?) = TMP([i]);
id TMP([i]?) * TMP([j]?) = 0;
id TMP([i]?) = njM^[i];

id Spinor(?s1) * GA(?g) * Spinor(?s2) =
  WeylChain(Spinor(?s1), ?g, Spinor(?s2));

toPolynomial;
endargument;
endargument;
#endif

id numM(qcM(?a)) = TMP(qcM(?a)) * numM(1);
also numM([x]?) * ORD(?i) = numM([x], List(?i));

id ORD(?i) = 1;

argument numM,1, TMP;
toPolynomial;
endargument;

id cutM([n]?, ?a) * numM(?x) = cutM([n], numM([n], ?x), ?a);

id TMP([x]?) = [x];
id numM([x]?, ?i) = [x];

argument cutM,2,3,4;
toPolynomial;
endargument;

else;

id numM([x]?) = [x];

endif;

#endif

*----------------------------------------------------------------------

.sort

#write "%X"

b SumOver, PowerOf, Mat, Den, IGram, A0, paveM, cutM;
print;

.end


