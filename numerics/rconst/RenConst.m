(*
	RenConst.m
		calculates the renormalization constants
		this file is part of FormCalc
		last modified 3 Mar 00 th
*)

<< ../../FormCalc.m;
<< ../NumPrep.m

$OnShell = False

Pair[k[i_], k[j_]] = K2

(* take only the transverse part of vector-boson SEs: *)

Pair[e[i_], k[j_]] = 0;
Pair[e[i_], e[j_]] = 1


AddScalar[{l1_, lr___}, x_] := {l1 + x, lr}


IntCollect[p__] := Plus[p] /; FreeQ[{p}, Re]

IntCollect[p__] := Collect[Plus[p], _Re, Simplify]

Simp[expr_] := Simplify[expr] /. Plus -> IntCollect

FastSimp[expr_] := expr /. Plus -> IntCollect


SelfEnergy[srcfile_] := 
Block[ {amps, res},
  Print[""];
  Print["processing ", srcfile];
  res = ProcessFile["fa/" <> srcfile, Classification -> IndexSumsOnly];
  res //. Abbreviations[]
]


SetOptions[WriteSummedExpr, PrintResult -> IfDebug]

WriteFortran[rcs_, file_] :=
Block[ {fi, s},
  fi = OpenWrite[
    StringForm[
      "!`1`r8_`2` > `3`_`4`dim.F",
      $NumPrepDir, $Platform, file, $Dimension ]//ToString,
    FormatType -> FortranForm,
    PageWidth -> 65];

  WriteString[ fi, "\
* Note: this file was generated automatically by RenConst.m.\n\
* Changes by hand will be overwritten the next time RenConst.m is run.\n\n\
\tsubroutine calc_renconst\n\
\timplicit character (a-s,u-z)\n\
\timplicit double complex (t)\n\n\
#include \"model.h\"\n\
#include \"looptools.h\"\n\
#include \"renconst.h\"\n\n\
\tinteger IndexDelta\n\
\texternal IndexDelta\n\n\
\tinteger g1" <>
    ({", ", ToString[#]}&)/@
      Union[Cases[rcs, SumOver[i_, _] -> i, Infinity]] <> "\n\n" ];

  WriteDoLoops[fi, ToDoLoops[rcs], WriteSummedExpr];

  WriteString[fi, "\n\tend\n"];
  Close[fi];
]


WriteMma[rcs_, file_] :=
Block[ {fi, SumOver, Sum},
  fi = OpenWrite[file <> "_" <> ToString[$Dimension] <> "dim.m"];

  WriteString[ fi, "\
(* This file contains the renormalization constants for the SM.\n\
   It was generated automatically by RenConst.m. DO NOT EDIT. *)\n"];

  SumOver/: SumOver[i__] x_ := Sum[x, {i}];

  ( WriteString[fi, "\n", #[[1]] /. g1 -> g1_, " = "];
    Write[ fi, Plus@@ #[[2]] ] )&/@ rcs;

  Close[fi];
]


IndexRange[g1] = 3

Mass[2] = MLE;
Mass[3] = MQU;
Mass[4] = MQD

mla = MLA = 0

	(* $HKSign is the sign of the SU(2) covariant derivative
	   1 = Haber-Kane, -1 = Denner conventions *)

RenConstSM :=
Block[ {sff, saa, saz, sww, szz, shh, spp, scc, tad, com, Zcom, $HKSign = -1},
  sff[1] = SelfEnergy["self.nn"];
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

  { Table[
	(* Note: it seems weird that the left-handed component sffL
	   is taken as the coefficient from (omp ** ga[...]): this
	   is because omp ** ga[...] = ga[...] ** omm *)
      sffL = Coefficient[ sff[i], Mat[omp ** ga[k[1]]] ];
      sffR = Coefficient[ sff[i], Mat[omm ** ga[k[1]]] ];
      sffS = Coefficient[ sff[i], Mat[omp] ]/Mass[i][g1];
      com = Expand[sffL + sffR + 2 sffS];
      Zcom = Mass[i][g1]^2 ReTilde[D[com, K2] /. K2 -> Mass[i][g1]^2];
      { dMf1[i, g1] ->
          Mass[i][g1]/2 ReTilde[com /. K2 -> Mass[i][g1]^2],
		(* assuming we're only interested in the
		   flavour-diagonal RCs for the moment: *)
        dZfL1[i, g1, g1] ->
          -ReTilde[sffL /. K2 -> Mass[i][g1]^2] - Zcom,
        dZfR1[i, g1, g1] ->
          -ReTilde[sffR /. K2 -> Mass[i][g1]^2] - Zcom } /.
        Mass[1][_] -> 0,
    {i, 1, 4}],

    dMZsq1 -> ReTilde[szz /. K2 -> MZ2],
    dMWsq1 -> ReTilde[sww /. K2 -> MW2],
    dMHsq1 -> ReTilde[shh /. K2 -> MH2],
    dZAA1  -> ReTilde[-D[saa, K2] /. K2 -> 0],
    dZAZ1  -> ReTilde[-2 saz/MZ2 /. K2 -> MZ2],
    dZZA1  -> ReTilde[2 saz/MZ2 /. K2 -> 0],
    dZZZ1  -> ReTilde[-D[szz, K2] /. K2 -> MZ2],
    dZchi1 -> ReTilde[-D[scc, K2] /. K2 -> MZ2],
    dZW1   -> ReTilde[-D[sww, K2] /. K2 -> MW2],
    dZphi1 -> ReTilde[-D[spp, K2] /. K2 -> MW2],
    dZH1   -> ReTilde[-D[shh, K2] /. K2 -> MH2],
    dWFZ1  -> 0,
    dWFW1  -> 0,
    dTad1  -> ReTilde[-tad],

	(* derived ones *)
    dSWsq1 -> -CW2 dMWsq1/MW2 + CW2 dMZsq1/MZ2,
    dCWsq1 -> -dSWsq1,
    dSW1   -> -$HKSign 1/2 dSWsq1/SW,
    dZe1   -> -1/2 (dZAA1 - $HKSign SW/CW dZZA1)
  }//Flatten
]


RenConstBgF :=
Block[ {saa, saz, sww, szz, shh, spp, scc, tad, $HKSign = -1},
  saa = -SelfEnergy["self.aa.bgf"];
  saz = -SelfEnergy["self.az.bgf"];
  sww = -SelfEnergy["self.ww.bgf"];
  szz = -SelfEnergy["self.zz.bgf"];
  shh = SelfEnergy["self.hh.bgf"];
  spp = SelfEnergy["self.pp.bgf"];
  scc = SelfEnergy["self.cc.bgf"];
  tad = SelfEnergy["tad.h.bgf"];

  { dMZsq1 -> ReTilde[szz /. K2 -> MZ2],
    dMWsq1 -> ReTilde[sww /. K2 -> MW2],
    dMHsq1 -> ReTilde[shh /. K2 -> MH2],
    dZAA1  -> ReTilde[-D[saa, K2] /. K2 -> 0],
    dZAZ1  -> -$HKSign 2 dCWsq1/(SW CW),
    dZZA1  -> 0,
    dZZZ1  -> dZAA1 - (CW2 - SW2)/(SW2 CW2) dCWsq1,
    dZW1   -> dZAA1 - dCWsq1/SW2,
    dZH1   -> dZW1 + dMWsq1/MW2,
    dZphi1 -> dZH1,
    dZchi1 -> dZH1,
    dWFZ1  -> -AddScalar[ReTilde[D[szz, K2] /. K2 -> MZ2], dZZZ1],
    dWFW1  -> -AddScalar[ReTilde[D[sww, K2] /. K2 -> MW2], dZW1],
    dTad1  -> ReTilde[-tad]
  }
]


RenConstMSSM :=
Block[ {sff, saa, saz, sww, szz, shh, spp, scc, tad, com, Zcom, $HKSign = 1},
  sff[1] = SelfEnergy["self.nn.mssm"];
  sff[2] = SelfEnergy["self.ee.mssm"];
  sff[3] = SelfEnergy["self.uu.mssm"];
  sff[4] = SelfEnergy["self.dd.mssm"];
  saa = -SelfEnergy["self.aa.mssm"];
  saz = -SelfEnergy["self.az.mssm"];
  sww = -SelfEnergy["self.ww.mssm"];
  szz = -SelfEnergy["self.zz.mssm"];
  spp = SelfEnergy["self.pp.mssm"];

  { Table[
	(* Note: it seems weird that the left-handed component sffL
	   is taken as the coefficient from (omp ** ga[...]): this
	   is because omp ** ga[...] = ga[...] ** omm *)
      sffL = Coefficient[ sff[i], Mat[omp ** ga[k[1]]] ];
      sffR = Coefficient[ sff[i], Mat[omm ** ga[k[1]]] ];
      sffS = Coefficient[ sff[i], Mat[omp] ]/Mass[i][g1];
      com = Expand[sffL + sffR + 2 sffS];
      Zcom = Mass[i][g1]^2 ReTilde[D[com, K2] /. K2 -> Mass[i][g1]^2];
      { dMf1[i, g1] ->
          Mass[i][g1]/2 ReTilde[com /. K2 -> Mass[i][g1]^2],
		(* assuming we're only interested in the
		   flavour-diagonal RCs for the moment: *)
        dZfL1[i, g1, g1] ->
          -ReTilde[sffL /. K2 -> Mass[i][g1]^2] - Zcom,
        dZfR1[i, g1, g1] ->
          -ReTilde[sffR /. K2 -> Mass[i][g1]^2] - Zcom } /.
        Mass[1][_] -> 0,
    {i, 1, 4}],

    dMZsq1 -> ReTilde[szz /. K2 -> MZ2],
    dMWsq1 -> ReTilde[sww /. K2 -> MW2],
    dZAA1  -> ReTilde[-D[saa, K2] /. K2 -> 0],
    dZAZ1  -> ReTilde[-2 saz/MZ2 /. K2 -> MZ2],
    dZZA1  -> ReTilde[2 saz/MZ2 /. K2 -> 0],
    dZZZ1  -> ReTilde[-D[szz, K2] /. K2 -> MZ2],
    dZW1   -> ReTilde[-D[sww, K2] /. K2 -> MW2],
    dZphi1 -> ReTilde[-D[spp, K2] /. K2 -> MHp2],

	(* derived ones *)
    dSWsq1 -> -CW2 dMWsq1/MW2 + CW2 dMZsq1/MZ2,
    dCWsq1 -> -dSWsq1,
    dSW1   -> -$HKSign 1/2 dSWsq1/SW,
    dZe1   -> -1/2 (dZAA1 - $HKSign SW/CW dZZA1)
  }//Flatten
]



CalcRC :=
Block[ {rc},
  rc = Simp[RenConstSM];
  WriteFortran[rc, "fortran/rc_SM"];
  WriteMma[rc, "mma/rc_SM"];

  If[ FileType["fa/self.aa.bgf"] === File,
    rc = OnePassOrder[
      Fold[#1 /. (#2[[1]] -> _) -> #2 &, rc, Simp[RenConstBgF]] ];
    WriteFortran[rc, "fortran/rc_SMbgf"];
    WriteMma[rc, "mma/rc_SMbgf"];
  ];

  If[ FileType["fa/self.aa.mssm"] === File,
    rc = FastSimp[RenConstMSSM];
    WriteFortran[rc, "fortran/rc_MSSM"];
    WriteMma[rc, "mma/rc_MSSM"];
  ];
]



$Dimension = D

CalcRC


$Dimension = 4

CalcRC

