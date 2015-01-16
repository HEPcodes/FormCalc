$Resonance = False;
  (* whether to handle the Higgs resonance *)
$ModifyingFunction = 1;
  (* additional multiplicative factor - e.g. MH2/S - 
     in applying the pole scheme *)
$Dyson = False;
  (* whether to apply Dyson summation *)

$BlockSize = 700;
  (* max LeafCount allowed per block *)

$AbbrFile = "abbr";
  (* where to look for the abbreviations *)
$BornFile = "born.m";
  (* appends this file automatically to the oneloop list *)
$RenConst = "rcsm.F";
  (* which file to use for the renormalization constants *)
$F77 = "f77 -r8 -extend_source";
  (* which compiler to use for compile Mma generated f77 code *)
$LoopToolsDir = "$(HOME)/LoopTools/";
  (* where LoopTools is installed *)

$NumPrepDir = StringReplace[
  If[ FileType[$Input] === File, $Input,
        (* if FormCalc was loaded from a directory in $Path: *)
    Block[ {full},
      Scan[
        If[ FileType[full = # <> "/" <> $Input] === File, Return[full] ]&,
        $Path ] ]
  ],
  {"NumPrep.m" -> "", "~" -> HomeDirectory[]} ]

$NumberMarks = False;

Get[$NumPrepDir <> "para.m"];

Off[FileDate::nffil];

DEN[S, MH2] = reso;
DEN[T, MH2] = resoT;
DEN[U, MH2] = resoU;
DEN[p_, m_] := (p - m)^-1;

out = (#//NN//FortranForm)&;

fmass[g_, a_, b_, c_] :=
  "\n\tgen = " <> ToString[g] <>
  "\n\tME2 = " <> ToString[a^2//out] <>
  "\n\tMD2 = " <> ToString[b^2//out] <>
  "\n\tMU2 = " <> ToString[c^2//out] <> "\n";
fma1 = fmass[1, fME, fMD, fMU];
fma2 = fmass[2, fMM, fMS, fMC];
fma3 = fmass[3, fML, fMB, fMT];

k[n_Integer] := n;
e[n_Integer] := ToExpression["i" <> ToString[n]];

A0[0] = 0;

Derivative[1, 0, 0][B0] = DB0;
Derivative[1, 0, 0][B1] = DB1;
Derivative[1, 0, 0][B11] = DB11;
Derivative[1, 0, 0][B00] = DB00;

ResoReplace[expr_] := expr /; FreeQ[expr, reso];
ResoReplace[expr_Times] :=
Block[ {stu, amp, yy, dsigma},
  amp = expr/reso;
  yy = amp /. U -> stu - S - T /. S -> MH2 /. stu -> S + T + U;
  yy *= $ModifyingFunction;
  If[ FreeQ[amp, reso],
    yy*resoc + (amp - yy)*reso,
  (* else *)
    dsigma = D[amp/reso, S] /. S -> MH2;
    dsigma *= $ModifyingFunction;
    dsigma*resoc + (amp - yy - dsigma)*reso
  ]
];
ResoReplace[plu_Plus] := ResoReplace/@ plu;

op[Plus] = "+";
op[Times] = "*";

Chunks[expr_] := expr /; LeafCount[expr] <= $BlockSize;
Chunks[expr_, pr_] := (
  WriteString[hh, "\n\t", pr, " = "];
  Write[hh,expr]; ) /; LeafCount[expr] <= $BlockSize;
Chunks[hd_[args__], r___] :=
  Chunks[#, r]&/@ hd[args] /; hd =!= Plus && hd =!= Times;
Chunks[expr_,r___] :=
Block[ {ex = Chunks/@ expr},
  If[LeafCount[ex] <= $BlockSize, Chunks, ChopChunks][ex,r]
];

SetAttributes[ChopChunks, HoldRest];
ChopChunks[expr_, var_:Unique["tmp"]] :=
Block[ {na = {}, re = {}, i = 0, l, hd = Head[expr], tmp = ToString[var]},
  Scan[(l = LeafCount[#];
    If[ (i += l) <= $BlockSize, re = {re, #},
      na = {na, hd@@ Flatten[re]}; re = {#}; i = l ])&, expr];
  i = "\n\t" <> tmp <> " = ";
  l = i <> tmp <> op[hd] <> "\n     -  ";
  Scan[(WriteString[hh, i]; Write[hh, #]; i = l)&,
    Flatten[{na, hd@@ Flatten[re]}]];
  tmp
];

Pieces[li_List, n_] :=
Block[ {k},
  If[ Length[li] < n, Return[{li}] ];
  k = Mod[Length[li], n];
  If[k === 0,
    Partition[li, n],
    Append[Partition[li, n], Take[li, -k]]
  ]
];

AddList[li_] :=
Block[ {r, j = 0},
  If[ li === {}, Return["0D0"] ];
  r = First[li] <> "()";
  Scan[(r = r <> " + ";
    If[ StringLength[r] - j > 45,
      j = StringLength[r = r <> "\n     -    "] ];
    r = r <> If[StringMatchQ[#, "*FT*"], "3D0*", ""] <> # <> "()")&,
    Rest[li]];
  r
];

WriteVars[vars_, type_, common_] :=
Block[ {v, s},
  If[ vars =!= {},
    s = StringReplace[ToString[ Pieces[First/@
      Select[vars, StringTake[ToString[#], 1] =!= "t" &], 6] ],
      {"{{" -> "}, {", "}}" -> ""}];
    v = StringReplace[s,
      {"}, {" -> "\n\t" <> type, "[gen]" -> "(3)"}];
    s = StringReplace[s,
      {"}, {" -> "\n\tcommon " <> common, "[gen]" -> ""}];
    WriteString[hh, v, s, "\n"]
  ]
];

WriteAbbr[abbr_, subr_] := (
  WriteString[hh, "\tsubroutine ", subr, "\n"];
  If[ abbr =!= {},
    WriteString[ hh,
      "\timplicit logical (a-s,u-z)\n",
      "\timplicit double complex (t)\n",
      "#include <vars.h>\n",
      "\tdouble complex Pair, Eps\n",
      "\texternal Pair, Eps\n",
      "\tdouble complex A0, B0, DB0, B1, DB1, B00, DB00, B11, DB11\n",
      "\texternal A0, B0, DB0, B1, DB1, B00, DB00, B11, DB11\n",
      "\tinteger Cget, Dget\n",
      "\texternal Cget, Dget\n\n"];
    (WriteString[hh, "\t", #[[1]]//FortranForm, " = "];
      Write[hh, #[[2]] ])&/@ abbr;
  ];
  WriteString[hh, "\tend\n\n"];
);

abfunc = A0 | B0 | B1 | B00 | B11 | DB0 | DB1 | DB00 | DB11;

C0[args__] = C0i[cc0, args];
D0[args__] = D0i[dd0, args];

abr[bca:abfunc[args__]] :=
Block[ {uu = Unique["ab"]},
  If[!FreeQ[{args}, ME2 | MD2 | MU2], uu = uu[gen]];
  cresolve[uu] = bca;
  abr[bca] = uu
];
abr[C0i[i_, args__]] :=
Block[ {uu = Unique["cd"]},
  If[!FreeQ[{args}, ME2 | MD2 | MU2], uu = uu[gen]];
  iresolve[uu] = Cget[args]//NN;
  abr[C0i[id_, args]] = Cval[id,uu];
  Cval[i, uu]
];
abr[D0i[i_, args__]] :=
Block[ {uu = Unique["cd"]},
  If[!FreeQ[{args}, ME2 | MD2 | MU2], uu = uu[gen]];
  iresolve[uu] = Dget[args]//NN;
  abr[D0i[id_, args]] = Dval[id, uu];
  Dval[i, uu]
];

MakeCode[indir_, outdir_, mask_:"*.m"] :=
Block[{i, j, hh, mf, amp, procs, borns, d2,
cbca, ibca, oname, fi, na, abbr},
  fi = Select[
    FileNames[Flatten[{mask, $BornFile}], indir],
    StringMatchQ[#, "*.m"]& ];
  If[$Dyson && StringMatchQ[StringJoin[fi], "*" <> indir <> "/self*"],
    Print["Warning: Applying Dyson summation to self energies!"]];

  mf = OpenWrite[outdir <> "Makefile"];
  WriteString[mf, "DIR = ", $NumPrepDir,
    "../\nLTDIR = ", $LoopToolsDir,
    "\nLIBS = -L$(LTDIR)$(HOSTTYPE) -lff\n",
    "FC = ", $F77, " -I. -I.. -I$(LTDIR)include",
    If[$Dyson, " -DDYSON", ""],
    " -g\nOBJS = renconst.o abbr.o"];

  procs = {};
  filelist = {};
  Do[
    oname = StringReplace[fi[[i]], {indir -> "", "/" -> ""}];
    If[ oname === $BornFile, oname = "born",
      oname = StringDrop[oname, -2] ];
    Print["processing class ", oname];

(********
 CAUTION: Using preprocessed files saves a lot of time,
          but can be fatal:
          *MUST* delete .mi-file when changing $Resonance
          or $ModifyingFunction or anything in para.m!!!
 ********)

    d2 = FileDate[fi[[i]] <> "i"];
    If[ d2 =!= $Failed && FromDate[fi[[i]]//FileDate] <= FromDate[d2],
      Print["using preprocessed file ", fi[[i]], "i"];
      amp = Get[fi[[i]] <> "i"],
    (* else *)
      amp = Get[fi[[i]] ];
      If[ $Resonance, amp = ResoReplace[amp] ];
      amp = Chop[NN[amp] /. {1. -> 1, -1. -> -1}, 10^-100] /. 0 -> 0.;
      Put[amp, fi[[i]] <> "i"]
    ];
    amp = amp /. bca:(abfunc | C0i | D0i)[__] -> abr[bca];
    AppendTo[procs, oname];
    WriteString[mf, " \\\n ", oname, ".o"];
    Print["  writing ", na = outdir <> oname <> ".F"];
    hh = OpenWrite["!sed 's/[\"\\]//g' > "<>na,
      FormatType -> FortranForm, PageWidth -> 67];
    WriteString[hh,
      "#include <defs.h>\n\n\tdouble complex function ", oname, "()\n",
      "\timplicit logical (a-s,u-z)\n",
      "\timplicit double complex (t)\n",
      "#include <vars.h>"];
    Chunks[amp, oname];
    WriteString[hh, "\n#ifdef DEBUG\n\tprint *, '", oname, " =',",
      oname, "\n#endif\n\tend\n"];
    Close[hh],
  {i, Length[fi]}];

  WriteString[mf,
    "\n\n.SUFFIXES: .f .F .o",
    "\n\nnum$(MH): $(DIR)num.F $(LTDIR)include/tools.F $(DIR)kin.h ",
    "feyn.F $(OBJS)\n",
    "\t$(FC) -DMH=$(MH) -o num$(MH) $(DIR)num.F $(OBJS) $(LIBS)\n\n",
    "renconst.o: ", $RenConst,
    "\n\t$(FC) -c -o renconst.o ", $RenConst,
    "\n\n.F.o: $(DIR)kin.h\n\t$(FC) -c $<\n\n",
    "clean:\n\trm -f *.o num$(MH)\n\n"];
  Close[mf];

  Print["writing ", na = outdir <> "vars.h"];
  hh = OpenWrite[na];
  WriteString[hh,
    "#include <kin.h>\n#include <rcsm.h>\n\n",
    "\tdouble complex Cval(13,1)\n\tcommon /cpave/ Cval\n\n",
    "\tdouble complex Dval(46,1)\n\tcommon /dpave/ Dval\n"];
  abbr = Get[indir <> $AbbrFile];
  WriteVars[
    If[ FreeQ[abbr, scale], abbr,
      WriteString[hh,
        "\n\tdouble precision scale\n\tcommon /sme/ scale\n"];
      DeleteCases[abbr, scale -> _] ],
    "double complex ", "/sme/ "];
  cbca = DownValues[cresolve] /. 
    cresolve -> Identity /.
    Literal -> Identity /. HoldPattern -> Identity;
  WriteVars[cbca, "double complex ", "/bca/ "];
  ibca = DownValues[iresolve] /.
    iresolve -> Identity /.
    Literal -> Identity /. HoldPattern -> Identity;
  WriteVars[ibca, "integer ", "/bca/ "];
  Close[hh];

  Print["writing ", na = outdir <> "abbr.F"];
  hh = OpenWrite["!sed 's/\\\\//g' > " <> na,
    FormatType -> FortranForm, PageWidth -> 67];
  WriteAbbr[Select[abbr, FreeQ[#, Plus]&]//ReleaseHold, "calc_sme"];
  WriteAbbr[Select[abbr, !FreeQ[#, Plus]&], "calc_abbr"];
  cbca = Join[cbca, ibca];
  ibca = Select[cbca, !FreeQ[#, gen]&];
  cbca = Join[Select[cbca, FreeQ[#, gen]&],
    ibca /. {gen -> 1,
      ME2 -> NN[fME^2], MD2 -> NN[fMD^2], MU2 -> NN[fMU^2]},
    ibca /. {gen -> 2,
      ME2 -> NN[fMM^2], MD2 -> NN[fMS^2], MU2 -> NN[fMC^2]},
    ibca /. {gen -> 3,
      ME2 -> NN[fML^2], MD2 -> NN[fMB^2], MU2 -> NN[fMT^2]}
  ];
  WriteAbbr[Select[cbca, FreeQ[#, T | U]&], "calc_bca"];
  WriteAbbr[Select[cbca, !FreeQ[#, T | U]&], "calc_bca2"];
  Close[hh];

  Print["writing ", na = outdir <> "feyn.F"];
  hh = OpenWrite[na];
  j = Count[cbca, Dget[__], {2}];
  i = Count[cbca, Cget[__], {2}] + 4*j;
  WriteString[hh,
    "*#define LAMBDA 1D10\n*#define MUDIM 1D10\n\n",
    "#define CSTORE ", 2*i + 1,
    "\n#define DSTORE ", 2*j + 1,
    "\n#define KINFAC (", out[kinfactor /. EE^2 -> EE2], ")\n"];
  If[ValueQ@@ #,
    WriteString[hh, "#define ", #//FortranForm, " ",
      ReleaseHold[#]//out, "\n"]]&/@
    Thread[HoldForm[{Alpha,
      MZ, GammaZ, MW, GammaW, MH,
      fME, fMD, fMU,
      fMM, fMS, fMC,
      fML, fMB, fMT}]];
  borns = Select[procs, StringMatchQ[#, "born*"]&];
  Which[
    Length[borns] === 0,
      WriteString[hh, "\n#define born() 0D0\n"],
    Length[borns] > 1,
      WriteString[hh, "\n#define BORN ",
        DeleteCases[borns, "born"][[1]], "\n"] ];
  procs = Complement[procs, borns];
  If[ Length[procs] === 0,
    WriteString[hh, "\n#define oneloop() 0D0\n"],
  (* else *)
    WriteString[hh,
      "\n\tdouble complex function oneloop()\n#include <kin.h>\n"];
    WriteString[hh, "\n\tdouble complex ", #,
      "\n\texternal ", #]&/@ procs;
    ferm = Select[procs,
      (StringMatchQ[#, "*F*"] && !StringMatchQ[#, "*FT*"])&];
    If[ ferm === {},
      WriteString[hh, "\n\n\toneloop = ", AddList[procs]],
      ferm = AddList[ferm];
      WriteString[hh, "\n",
        fma1, "\toneloop = ", AddList[procs],
        fma2, "\toneloop = oneloop + ", ferm,
        fma3, "\toneloop = oneloop + ", ferm]
    ];
    WriteString[hh, "\n\tend\n"]
  ];
  Close[hh];
];

Print["NumPrep loaded"];
