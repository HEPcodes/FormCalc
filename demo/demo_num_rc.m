source = "rc/";
dest = "fortran/"

<< para.m

A0[0] = 0;
Derivative[1, 0, 0][B0] = DB0;
Derivative[1, 0, 0][B1] = DB1;
Derivative[1, 0, 0][B11] = DB11;
Derivative[1, 0, 0][B00] = DB00

abfunc = A0 | B0 | B1 | B00 | B11 | DB0 | DB1 | DB00 | DB11

InsertFermions[x_] := (
  (x /. {ME2 -> fME^2, MD2 -> fMD^2, MU2 -> fMU^2}) +
  (x /. {ME2 -> fMM^2, MD2 -> fMS^2, MU2 -> fMC^2}) +
  (x /. {ME2 -> fML^2, MD2 -> fMB^2, MU2 -> fMT^2})
)

SelfEnergy[x_] :=
  Expand[
    Get[source <> x <> ".m"] +
    InsertFermions[Get[source <> x <> "F.m"]] ]

exco[rc_] := Collect[rc//Expand, abfunc[__]]

re[rc_] := exco[rc] /. x:abfunc[__] :> dble[x] /.
  Complex[a_, b_] -> a

im[rc_] := Select[exco[rc], !FreeQ[#, abfunc]&] /.
  x:abfunc[__] :> dimag[x] /. Complex[a_,b_] -> b

saa = -SelfEnergy["saa"];
saz = -SelfEnergy["saz"];
sww = -SelfEnergy["sww"];
szz = -SelfEnergy["szz"];
shh = SelfEnergy["shh"];
spp = SelfEnergy["spp"];
scc = SelfEnergy["scc"];

v[GammaHMH] = im[shh /. K2 -> MH2]

v[dMZsq1] = re[szz /. K2 -> MZ2];
v[dMWsq1] = re[sww /. K2 -> MW2];
v[dMHsq1] = re[shh /. K2 -> MH2]

v[dSWsq1] = -CW^2*dMWsq1/MW2 + CW^2*dMZsq1/MZ2;
v[dCWsq1] = -dSWsq1;
v[dSW1]   = 1/2*dSWsq1/SW

v[dZAA1]  = -re[D[saa,K2] /. K2 -> 0];
v[dZAZ1]  = re[(-2*saz/MZ2) /. K2 -> MZ2];
v[dZZA1]  = re[(2*saz/MZ2) /. K2 -> 0];
v[dZZZ1]  = re[-D[szz,K2] /. K2 -> MZ2];
v[dZchi1] = re[-D[scc,K2] /. K2 -> MZ2];
v[dZW1]   = re[-D[sww,K2] /. K2 -> MW2];
v[dZphi1] = re[-D[spp,K2] /. K2 -> MW2];
v[dZH1]   = re[-D[shh,K2] /. K2 -> MH2]

v[dWFZ]   = 0;
v[dWFW]   = 0

v[dZe1]     = -1/2*(dZAA1 + SW/CW*dZZA1);

v[shhren]   = exco[shh - dMHsq1 + dZH1*(K2 - MH2)]


fi = OpenWrite[dest <> "rcsm.h"];
WriteString[fi,
  "\tdouble precision GammaHMH\n",
  "\tdouble precision dMZsq1, dMWsq1, dMHsq1, dSWsq1, dCWsq1, dSW1\n",
  "\tdouble precision dZAA1, dZAZ1, dZZA1, dZZZ1, dZW1, dZH1\n",
  "\tdouble precision dZphi1, dZchi1, dWFZ, dWFW, dZe1\n",
  "\tcommon /renconst/ GammaHMH,\n",
  "     +    dMZsq1, dMWsq1, dMHsq1, dSWsq1, dCWsq1, dSW1,\n",
  "     +    dZAA1, dZAZ1, dZZA1, dZZZ1, dZW1, dZH1,\n",
  "     +    dZphi1, dZchi1, dWFZ, dWFW, dZe1\n"];
Close[fi];

Nv[x_] := NN[v[x]] /. 1.*a_ -> a /. -1.*a_ -> -a /. 0 -> 0.

fi = OpenWrite["!sed 's/\\\\//g' > " <> dest <> "rcsm.F",
  FormatType -> FortranForm, PageWidth -> 65]

externals =
  "#include \"kin.h\"\n" <>
  "#include \"rcsm.h\"\n\n" <>
  "\tdouble complex A0, B0, B1, B00, B11, DB0, DB1, DB00, DB11\n" <>
  "\texternal A0, B0, B1, B00, B11, DB0, DB1, DB00, DB11\n"

WriteString[ fi,
  "* This file contains the renormalization constants and\n",
  "* the H-H-Selfenergy for the SM\n",
  "* It was generated automatically by ", $Input, "\n",
  "* DO NOT EDIT.\n\n",
  "\tsubroutine calc_renconst\n",
  "\timplicit none\n",
  externals
]

(WriteString[fi, "\n\n\t", #, " = "];
 Write[fi, #//Nv];
 WriteString[fi, "\tprint *, \"", #, " =\",", #])&/@
  { GammaHMH,
    dMZsq1, dMWsq1, dMHsq1, dSWsq1, dCWsq1, dSW1,
    dZAA1, dZAZ1, dZZA1, dZZZ1, dZW1, dZH1,
    dZphi1, dZchi1, dWFZ, dWFW, dZe1 }

WriteString[fi, "\n\tend\n\n",
  "#ifdef DYSON\n",
  "\tdouble complex function shh(K2)\n",
  "\timplicit none\n",
  "\tdouble precision K2\n",
  externals,
  "\n\tshh="];
Write[fi, shhren//Nv];
WriteString[fi, "\tend\n#endif\n"]

Close[fi];

