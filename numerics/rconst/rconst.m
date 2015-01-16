(*
	rconst.m
		calculates the renormalization constants
		last modified 2 Jul 99 th
*)


$HKSign = -1	(* 1 = Haber-Kane, -1 = Denner conventions 
		   $HKSign is the sign of the SU(2) covariant
		   derivative *)

<< ~/FormCalc/FormCalc.m

$OnShell = False

Pair[k[i_], k[j_]] = K2

(* take only the transverse part of vector-boson SEs *)

Pair[e[i_], k[j_]] = 0;
Pair[e[i_], e[j_]] = 1


NN[expr_] := Chop[N[expr, 35], 10^-30] //.
  {1. a_ -> a, -1. a_ -> -a} /. 0 -> 0. /.
  rc:$IntArgFuncs[__] :> (rc /. r_Real :> Floor[r])

	(* these functions must always have integer arguments *)
$IntArgFuncs = SumOver | dZfL1 | dZfR1 | dMf1 | Af |
  MSf | MSf2 | MCha | MCha2 | MNeu | MNeu2 |
  USf | USfC | UCha | UChaC | VCha | VChaC | ZNeu | ZNeuC


(* The next Block is admittedly a dirty trick: the idea is to get the
   numerical definitions in -- not for the symbols directly, but only
   for N[sym], so as not to interfere with FormCalc. *)

Block[ {Set},
  (sym_ = val_) := ((N[#1, _] := #2)&@@ {sym, NN[val]}; NN[val]);
  << ~/FormCalc/numerics/para.m
]

mla = MLA = 0


(* The function Chunks is from NumPrep -- it cuts a Fortran
   expression into blocks which the compiler can handle *)

$BlockSize = 700

Chunks[expr_, var_, addtovar_] :=
Block[ {theexpr = expr},
  If[ LeafCount[expr] > $BlockSize,
    If[ Head[expr] === Plus && !MemberQ[var, _Symbol],
      theexpr = List@@ expr /. Plus -> PlusChop;
      ChopUp[Plus@@ theexpr, var, addtovar];
      Return[] ];
    theexpr = expr /. Plus -> PlusChop
  ];
  WriteString[fi, "\n\t", var//FortranForm, " = "];
  If[addtovar,
    WriteString[fi, var//FortranForm, " +\n     -  "] ];
  Write[fi, theexpr];
]


PlusChop[expr__] := ChopUp[Plus[expr], Unique["tmp"], False] /;
  LeafCount[{expr}] > $BlockSize

PlusChop[expr__] := Plus[expr]

ChopUp[expr_, var_, addtovar_] :=
Block[ {i, l = ToString[var//FortranForm]},
  i = "\n\t" <> l <> " = ";
  l = i <> l <> " +\n     -  ";
  If[addtovar, i = l];
  Scan[(WriteString[fi, i]; Write[fi, #]; i = l)&,
    Flatten[CodeBlocks@@ expr]];
  var
]


CodeBlocks[a_, b_, r___] :=
  CodeBlocks[a + b, r] /; LeafCount[{a, b}] < $BlockSize

CodeBlocks[a_, r___] := {a, CodeBlocks[r]}

CodeBlocks[] = Sequence[]


WriteOut[expr_, var_] :=
Block[ {addtovar, loops, Conjugate = dconjg},
  addtovar = False;
  Scan[
    ( loops = Cases[Collect[#, _SumOver, simp], _SumOver];
      Apply[
        WriteString[fi, "\n\tdo ", #1, " = 1, ", #2]&,
        loops, 1];
      Chunks[# /. _SumOver -> 1, var, addtovar];
      WriteString[fi,
        StringJoin[Table["\tenddo\n", {Length[loops]}]]];
      addtovar = True )&,
    Flatten[{expr}] ];
]


simp[expr_] := Collect[expr, _dble | _dimag]

abfunc = A0 | B0 | B1 | B00 | B11 | DB0 | DB1 | DB00 | DB11

externals =
  "\timplicit logical (a-s,u-z)\n" <>
  "\timplicit double complex (t)\n" <>
  "#include \"kin.h\"\n" <>
  "#include \"rcsm.h\"\n\n" <>
  "\tinteger IndexDelta\n" <>
  "\texternal IndexDelta\n" <>
  "\tdouble complex A0, B0, B1, B00, B11, DB0, DB1, DB00, DB11\n" <>
  "\texternal A0, B0, B1, B00, B11, DB0, DB1, DB00, DB11\n"

SetAttributes[{re, im}, Listable]

re[rc_] := rc /. x:abfunc[__] :> dble[x]

im[rc_] := (rc /. x:abfunc[__] :> dimag[x]) - (rc /. abfunc[__] -> 0)

AddScalar[{l1_, lr___}, x_] := {l1 + x, lr}

dble[0] = dimag[0] = 0

Mass[1][_] = 0;
Mass[2] = MLE;
Mass[3] = MQU;
Mass[4] = MQD


SelfEnergy[srcfile_] := 
Block[ {amps, res},
  Print[""];
  Print["processing ", srcfile];
  res = ProcessFile["fa/" <> srcfile,
    Classification -> IndexSumsOnly];
  res //. Abbreviations[]
]


CalcRC[] := (
  sff[2] = SelfEnergy["self.ee"];
  sff[3] = SelfEnergy["self.uu"];
  sff[4] = SelfEnergy["self.dd"];
  saa = -SelfEnergy["self.aa"];
  saz = -SelfEnergy["self.az"];
  sww = -SelfEnergy["self.ww"];
  szz = -SelfEnergy["self.zz"];
  shh = SelfEnergy["self.hh"];
  spp = SelfEnergy["self.pp"];
  scc = SelfEnergy["self.cc"];
  tad = SelfEnergy["tad.h"];

  Do[
    sffL = Coefficient[sff[i], omp ** ga[k[1]]];
    sffR = Coefficient[sff[i], omm ** ga[k[1]]];
    sffS = Coefficient[sff[i], omp]/Mass[i][g1];
    com = Expand[sffL + sffR + 2 sffS];
    v[dMf1[i, g1]] = Mass[i][g1]/2 re[com /. K2 -> Mass[i][g1]^2];
	(* assuming we're only interested in the
	   flavour-diagonal RCs for the moment *)
    com = Mass[i][g1]^2 re[D[com, K2] /. K2 -> Mass[i][g1]^2];
    v[dZfL1[i, g1, g1]] = -re[sffL /. K2 -> Mass[i][g1]^2] - com;
    v[dZfR1[i, g1, g1]] = -re[sffR /. K2 -> Mass[i][g1]^2] - com,
  {i, 2, 4}];

  v[dMZsq1] = re[szz /. K2 -> MZ2];
  v[dMWsq1] = re[sww /. K2 -> MW2];
  v[dMHsq1] = re[shh /. K2 -> MH2];
  v[dZAA1] = -re[D[saa, K2] /. K2 -> 0];
  v[dZAZ1] = re[-2 saz/MZ2 /. K2 -> MZ2];
  v[dZZA1] = re[2 saz/MZ2 /. K2 -> 0];
  v[dZZZ1] = re[-D[szz, K2] /. K2 -> MZ2];
  v[dZchi1] = re[-D[scc, K2] /. K2 -> MZ2];
  v[dZW1] = re[-D[sww, K2] /. K2 -> MW2];
  v[dZphi1] = re[-D[spp, K2] /. K2 -> MW2];
  v[dZH1] = re[-D[shh, K2] /. K2 -> MH2];
  v[dWFZ1] = 0;
  v[dWFW1] = 0;
  v[dTad1] = re[-tad];
  v[GammaHMH] = im[shh /. K2 -> MH2];

	(* derived ones *)
  v[dSWsq1] = -CW^2 dMWsq1/MW2 + CW^2 dMZsq1/MZ2;
  v[dCWsq1] = -dSWsq1;
  v[dSW1] = -$HKSign 1/2 dSWsq1/SW;
  v[dZe1] = -1/2 (dZAA1 - $HKSign SW/CW dZZA1);

	(* and the renormalized Higgs SE, which is needed
	   for Dyson summation *)
  v[shhren] = AddScalar[shh, -dMHsq1 + dZH1 (K2 - MH2)];

  WriteFortran[""];

  If[ FileType["fa/self.aa.bfm"] === File,
    saa = -SelfEnergy["self.aa.bfm"];
    saz = -SelfEnergy["self.az.bfm"];
    sww = -SelfEnergy["self.ww.bfm"];
    szz = -SelfEnergy["self.zz.bfm"];
    shh =  SelfEnergy["self.hh.bfm"];
    spp =  SelfEnergy["self.pp.bfm"];
    scc =  SelfEnergy["self.cc.bfm"];
    tad =  SelfEnergy["tad.h.bfm"];

    v[dMZsq1] = re[szz /. K2 -> MZ2];
    v[dMWsq1] = re[sww /. K2 -> MW2];
    v[dMHsq1] = re[shh /. K2 -> MH2];
    v[dZAA1] = -re[D[saa, K2] /. K2 -> 0];
    v[dZAZ1] = -$HKSign 2 dCWsq1/(SW CW);
    v[dZZA1] = 0;
    v[dZZZ1] = dZAA1 - (CW^2 - SW^2)/(SW^2 CW^2) dCWsq1;
    v[dZW1] = dZAA1 - dCWsq1/SW^2;
    v[dZH1] = dZW1 + dMWsq1/MW2;
    v[dZphi1] = dZH1;
    v[dZchi1] = dZH1;
    v[dWFZ1] = -(re[D[szz, K2] /. K2 -> MZ2] + dZZZ1)/2;
    v[dWFW1] = -(re[D[sww, K2] /. K2 -> MW2] + dZW1)/2;
    v[dTad1] = re[-tad];
    v[GammaHMH] = im[shh /. K2 -> MH2];

    v[shhren] = AddScalar[shh, -dMHsq1 + dZH1 (K2 - MH2)];

    WriteFortran["bfm"];
  ];
)


WriteFortran[bfm_] :=
Block[ {fi},
  fi = OpenWrite["fortran/rcsm.h"];
  WriteString[fi,
    "\tdouble complex dZfL1(4, 3, 3), dZfR1(4, 3, 3), dMf1(4, 3)\n",
    "\tdouble complex dMZsq1, dMWsq1, dMHsq1, dSWsq1, dCWsq1, dSW1\n",
    "\tdouble complex dZAA1, dZAZ1, dZZA1, dZZZ1, dZW1, dZH1\n",
    "\tdouble complex dZphi1, dZchi1, dWFZ1, dWFW1, dZe1, dTad1\n",
    "\tdouble complex GammaHMH\n",
    "\tcommon /renconst/ dZfL1, dZfR1, dMf1,\n",
    "     +    dMZsq1, dMWsq1, dMHsq1, dSWsq1, dCWsq1, dSW1,\n",
    "     +    dZAA1, dZAZ1, dZZA1, dZZZ1, dZW1, dZH1,\n",
    "     +    dZphi1, dZchi1, dWFZ1, dWFW1,\n",
    "     +    dZe1, dTad1, GammaHMH\n"];
  Close[fi];

  fi = OpenWrite[
    StringForm[
      "!`1`numerics/r8_`2` > fortran/rcsm`3`_dim`4`.F",
      $FormCalcDir, $Platform,
      bfm, $Dimension,
    ]//ToString,
    FormatType -> FortranForm,
    PageWidth -> 65];

  WriteString[ fi,
    "* This file contains the renormalization constants and\n",
    "* the Higgs self-energy for the SM\n",
    "* This file was generated automatically by rconst.m.\n",
    "* DO NOT EDIT.\n\n",
    "\tsubroutine calc_renconst\n",
    externals,
    "\n\tinteger g1",
    StringJoin[{", ", ToString[#]}&/@
      Cases[DownValues[v], SumOver[i_, _] -> i, Infinity]//Union],
    "\n" ];

  WriteString[fi, "\n\tdo g1 = 1, 3"];
  ( WriteString[fi, "\n"];
    WriteOut[v[#]//NN, #];
    WriteString[fi, "\n\tprint *, '",
      StringReplace[ToString[FortranForm[#]], "g1" -> "',g1,'"],
      " =', ", FortranForm[#]] )&/@
      Flatten[Table[{dMf1[i, g1], dZfL1[i, g1, g1], dZfR1[i, g1, g1]},
        {i, 2, 4}]];
  WriteString[fi, "\n\tenddo"];

  ( WriteString[fi, "\n"];
    WriteOut[v[#]//NN, #];
    WriteString[fi, "\n\tprint *, '", #, " =', ", #] )&/@
      { dMZsq1, dMWsq1, dMHsq1, dSWsq1, dCWsq1, dSW1,
        dZAA1, dZAZ1, dZZA1, dZZZ1, dZW1, dZH1,
        dZphi1, dZchi1, dWFZ1, dWFW1,
        dZe1, dTad1, GammaHMH };

  WriteString[fi, "\n\tend\n\n",
    "#ifdef DYSON\n",
    "\tdouble complex function shh(K2)\n",
    "\tdouble precision K2\n",
    externals];
  WriteOut[v[shhren]//NN, "shh"];
  WriteString[fi, "\tend\n#endif\n"];
  Close[fi];
]



$Dimension = D;
CalcRC[]

$Dimension = 4;
CalcRC[]

