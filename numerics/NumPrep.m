(*
	NumPrep.m
		prepares Fortran code from FormCalc output
		this file is part of FormCalc
		last modified 14 May 00 th

NumPrep is a program for generating Fortran code out of FormCalc output.
It has been tested only for 2 -> 2 processes.

The philosophy behind NumPrep is that the user should NOT HAVE TO MODIFY
the code after NumPrep did its job since that is typically the point
where most bugs creep in. NumPrep thus produces all modules of code, the
necessary include files, and a Makefile.

*)


Print[""];
Print["NumPrep 2"];
Print["by Thomas Hahn"];
Print["last revision: 14 May 00"];


NumPrep`$NumPrepDir =
  If[ FileType[$Input] === File, $Input,
        (* or, if NumPrep was loaded from a directory in $Path: *)
    Block[ {full},
      Scan[
        If[ FileType[full = # <> "/" <> $Input] === File, Return[full] ]&,
        $Path ] ]
  ];
Block[ {pos = StringPosition[NumPrep`$NumPrepDir, "/"]},
  NumPrep`$NumPrepDir = If[ Length[pos] === 0, ".",
    StringTake[ NumPrep`$NumPrepDir, pos[[-1, -1]] ] ];
];
	(* expand relative paths *)
NumPrep`$NumPrepDir = SetDirectory[NumPrep`$NumPrepDir] <> "/";
ResetDirectory[]


BeginPackage["FormCalc`"]

(* declare some FormCalc symbols: *)

{ A0, B0, B1, B00, B11, DB0, DB1, DB00, DB11,
  C0i, Cget, D0i, Dget,
  IndexDelta, SumOver, Pair, Eps,
  DEN, S, T, U, Sf, Tf, Uf,
  e, k, s, Hel, Scale, Mat, Spinor }

EndPackage[]


BeginPackage["NumPrep`"]

(* the main function MakeCode and its options *)

MakeCode::usage = "MakeCode[indir, outdir, mask] is the main function of
NumPrep. It takes two mandatory arguments: fromdir points to a directory
containing the input files, and outdir to the directory where the Fortran
files are put. mask is an optional parameter (default value *.m) which
specifies a mask to select the input files."

Legs::usage = "Legs is an option of MakeCode. It specifies the number of
external legs to generate code for."

AbbrFile::usage = "AbbrFile is an option of MakeCode and specifies the
file in which the abbreviations were saved."

MatFile::usage = "MatFile is an option of MakeCode and specifies the file
in which the matrix elements (the Mat[i, j]) were saved."

BornMask::usage = "BornMask is an option of MakeCode. It singles out the
code modules which contain Born amplitudes."

RenConst::usage = "RenConst is an option of MakeCode. It tells MakeCode
how the file containing the renormalization constants is called."

Drivers::usage = "Drivers points to a directory containing customised
versions of the driver programs for running the generated Fortran code.
This directory need not contain all driver programs: files not contained
therein are taken from the default directory, " <> $NumPrepDir <>
"drivers."


(* utility routines, mainly for generating Fortran code *)

ReTilde::usage = "ReTilde[expr] is a special version of Re used for
computing renormalization constants. It takes the real part only of the
loop integrals occurring in expr."

ImTilde::usage = "ImTilde[expr] is a special version of Im used for
computing renormalization constants. It takes the imaginary part only of
the loop integrals occurring in expr."

ToFortran::usage = "ToFortran[expr] returns the Fortran form of expr as
a string."

SizeSplit::usage = "SizeSplit[expr, n] splits expr into subexpressions
whose LeafCount is at most n. The list of subexpressions is arranged such
that Join@@ SizeSplit[expr, n] === expr."

WriteExpr::usage = "WriteExpr[file, {var1 -> expr1, var2 -> expr2, ...},
addto] writes a list of variable assignments in Fortran format to
file. If an expression's leaf count is larger than $BlockSize, the
expression is cut into blocks, viz.\n\
\tvar = part1\n\
\tvar = var + part2\n\
\t...\n
addto is an optional argument: if True, the variables var1, var2, ... are
always added to, i.e. \"var = var + expr\"."

NewLine::usage = "NewLine is an option of WriteExpr. If set to True,
WriteExpr writes out an empty line (\\n) after each statement."

Optimize::usage = "Optimize is an option of WriteExpr. With Optimize ->
True, WriteExpr introduces variables for subexpressions which are used
more than once."

WriteSummedExpr::usage = "WriteSummedExpr[file, var -> expr] writes
\"var = expr\" in Fortran format to file. Unlike WriteExpr, it writes out
code to perform any summations over indices occurring in expr."

PrintResult::usage = "PrintResult is an option of WriteSummedExpr. With
PrintResult -> True, WriteSummedExpr puts a print statement in the Fortran
code after the variable assignment to print out the result. With
PrintResult -> IfDebug, the print statement is wrapped in \"#ifdef DEBUG\"
and \"#endif\"."

IfDebug::usage = "IfDebug is a possible selection for the option
PrintResult of WriteSummedExpr. It specifies that the print statement
written out by WriteSummedExpr is wrapped in \"#ifdef DEBUG\" and
\"#endif\"."

ToDoLoops::usage = "ToDoLoops[list, ifunc] splits list into patches which
must be summed over the same set of indices. ifunc is an optional
argument: ifunc[expr] must return the indices occurring in expr."

DoLoop::usage = "DoLoop[ind, expr] is a symbol introduced by ToDoLoops
indicating that expr is to be summed over the set of indices ind."

WriteDoLoops::usage = "WriteDoLoops[file, list, func] writes the do-loops
occurring in list to file. For writing the actual expressions, func[file,
item] is called."

IndexRange::usage = "IndexRange[i] gives the range of index i."

OnePassOrder::usage = "OnePassOrder[defs] orders a list of definitions
such that each definition comes before its first use. In Mathematica
lingo, the effect of OnePassOrder can be shown as\n
\texpr //. defs === expr /. defs"


(* symbols used in the Fortran code *)

Ctree::usage = "Ctree[Fi] is the coefficient of Fi of the tree-level
amplitude."

Cloop::usage = "Cloop[Fi] is the coefficient of Fi of the one-loop
amplitude."

MomSquare::usage = "MomSquare[a, b] calculates the difference of the
momenta a and b squared, i.e. Pair[k[a] - k[b], k[a] - k[b]]."

dconjg::usage = "dconjg[z] takes the complex conjugate of z in Fortran."

dble::usage = "dble[z] takes the real part of z in Fortran."

dimag::usage = "dimag[z] takes the imaginary part of z in Fortran."

Cval::usage = "Cval is the array containing the cached three-point
integrals in LoopTools."

Dval::usage = "Dval is the array containing the cached four-point
integrals in LoopTools."


(* system variables *)

$NumPrepDir::usage = "$NumPrepDir is the path where NumPrep and its
companion files are located."

$GenerateCodeFor::usage = "$GenerateCodeFor specifies the platform for
which Fortran code is produced. This affects only the calling sequence for
the Fortran compiler, and is closely connected with the $F77 variable. If
you want the code to be produced always for the same architecture on which
NumPrep is running on, put $GenerateCodeFor = $SystemID."

$ExtraLibs::usage = "$ExtraLibs specifies additional libraries needed for
linking the final executable. By default, this is the CERNlib."

$LoopToolsDir::usage = "$LoopToolsDir is the path where the include files
for LoopTools are located."

$RenConstDir::usage = "$RenConstDir is the path where the Fortran files 
for calculating the renormalization constants is located."

$BlockSize::usage = "$BlockSize is the maximum LeafCount a single Fortran
statement written out by WriteExpr may have. Any expression with
LeafCount > $BlockSize will be chopped up before being written to the
Fortran file."

$FileSize::usage = "$FileSize gives the maximum LeafCount the expressions
in a single Fortran file may have. If the expressions grow larger than
$FileSize, the file is split into several pieces."

$F77::usage = "$F77 specifies the command line for the Fortran compiler."


Begin["`Private`"]

PrependTo[$ContextPath, "FormCalc`"]

$NumberMarks = False

Off[CopyFile::filex]


(* using the convention k[n_] = 4 (n - 1) + 1, we have: *)

MandelstamDef[T] := MomSquare[1, 9];		(* (k[1] - k[3])^2 *)
MandelstamDef[U] := MomSquare[1, 13];		(* (k[1] - k[4])^2 *)
MandelstamDef[Sf] := MomSquare[9, -13];		(* (k[3] + k[4])^2 *)
MandelstamDef[Tf] := MomSquare[5, 13];		(* (k[2] - k[4])^2 *)
MandelstamDef[Uf] := MomSquare[5, 9];		(* (k[2] - k[3])^2 *)

MandelstamVars = {T, U, Sf, Tf, Uf}


ToFortran[s_String] = s

ToFortran[s_] := ToString[FortranForm[s]]


ToSymbol[x__] := ToExpression[ StringJoin[ToString/@ Flatten[{x}]] ]


ToList[a_Max] := List@@ a

ToList[a_List] = a

ToList[a_] = {a}


(* convert FormCalc's abbreviations (symbols) to arrays which are easier
   to handle in Fortran: *)

ToArray[s_Symbol] :=
Block[ {c = Characters[ToString[s]], h},
  {h, c} = ToExpression[StringJoin[#]]&/@
    {Select[c, LetterQ], Select[c, DigitQ]};
  s = h[c]
]

ToArray[s_] = s


Unprotect[Max];
Max[s_[i_Integer], s_[j_Integer]] := s[Max[i, j]];
Protect[Max]


MatType[_Mat, _] = 1

MatType[h_[i_], h_[j_]] := MatType[h][i, j]

MatType[h_] := MatType[h] = ToSymbol["Mat", h]

Mat[0] = 0


A0[0] = 0

loopint = A0 | B0 | B1 | B00 | B11 | DB0 | DB1 | DB00 | DB11 | C0i | D0i

loopintForArgType = Flatten[loopint | Cget | Dget]

loopintForResType = Flatten[loopint | Cval | Dval]


(* These are special versions of Re and Im used for computing
   renormalization constants (see A. Denner, Forts. Phys. 41 (1993) 307).
   The real and imaginary part is taken only of the loop integrals. *)

ReTilde[expr_] := expr /. int:loopintForResType[__] :> Re[int]

ImTilde[expr_] :=
  (expr /. int:loopintForResType[__] :> Im[int]) -
  (expr /. loopintForResType[__] -> 0)


r8cmd = "!" <> $NumPrepDir <> "r8_" <> Environment["HOSTTYPE"] <> " > "

OpenWriteFortran[file_] :=
  OpenWrite[r8cmd <> outdir <> file <> ".F",
    FormatType -> FortranForm, PageWidth -> 67]


(* The following routines are concerned with breaking a large
   expression into pieces the Fortran compiler will compile.
   This is controlled by two variables:

   - $BlockSize is the maximum LeafCount a single Fortran statement
     may have. The cutting-up of expressions into such blocks is
     performed by the function WriteExpr.

   - $FileSize is the maximum LeafCount a whole file may have.
     The function which separates large expressions into file-size
     fragments is SizeSplit. *)

SizeSplit[expr_, n_] :=
Block[ {maxsize = n, cb},
  List@@ Flatten[ Operate[Coalesce, expr] ]
]

Coalesce[h_][a_, b_, r___] :=
  Coalesce[h][{a, b}, r] /; LeafCount[{a, b}] < maxsize

Coalesce[h_][a_, r___] := cb[ h@@ Flatten[{a}], Coalesce[h][r] ]

Coalesce[_][] = Sequence[]


Options[WriteExpr] = {NewLine -> False, Optimize -> False}

WriteExpr[fi_, expr_, addtovar_:False, opt___Rule] :=
Block[ {hh = fi, newline, optim, addto,
Conjugate = dconjg, Re = dble, Im = dimag},
  {newline, optim} = {NewLine, Optimize} /. {opt} /. Options[WriteExpr];
  If[ newline, newline := WriteString[hh, "\n"] ];
  If[ addtovar,
    Cases[{expr}, (var_ -> _) :> (addto[var] = True), Infinity] ];
  WriteAssign[ If[optim, RemoveRedundancy[expr], expr] ];
]

RemoveRedundancy[expr_] :=
Block[ {dups, objs, vars, tmpdefs = {}, optexpr = expr},
  Block[ {Plus},
    Do[
      objs = Cases[optexpr, p_Plus /; LeafCount[N[p]] > 10, {-i}];
      dups = Union[objs];
      dups = Select[dups,
        Length[Position[objs, #, {1}, 2, Heads -> False]] > 1 &,
        Length[objs] - Length[dups]];
      vars = Table[Unique["tmp"], {Length[dups]}];
      tmpdefs = {tmpdefs, Hold@@ Thread[vars -> dups]};
      MapThread[Set, {dups, vars}];
      optexpr = optexpr,
    {i, 3, Depth[expr] - 1}];
  ];
  Flatten[{ReleaseHold[tmpdefs], optexpr}]
]


Attributes[WriteAssign] = {Listable}

WriteAssign[var_ -> expr_, addtovar_] :=
  WriteBlock[var, expr] /; LeafCount[expr] < $BlockSize

WriteAssign[var_ -> expr_Plus] :=
  ChopUp[var, Plus@@ (List@@ expr /. Plus -> PlusChop)]

WriteAssign[var_ -> expr_] :=
  WriteBlock[var, expr /. Plus -> PlusChop]


PlusChop[expr__] := ChopUp[Unique["tmp"], Plus[expr]] /;
  LeafCount[{expr}] > $BlockSize

PlusChop[expr__] := Plus[expr]

ChopUp[var_, expr_] :=
  (Scan[WriteBlock[var, #]&, SizeSplit[expr, $BlockSize]]; var)

WriteBlock[var_, expr_] :=
Block[ {for},
  for = expr /. int:loopintForArgType[__] :> Nargs/@ int /.
    p:_Integer^_Rational :> N[p, 35] /.
    Times -> OptTimes;
  If[ TrueQ[addto[var]], for += var, addto[var] = True ];
  Write[hh, var -> for];
  newline
]

Nargs[i_Integer] := N[i]

Nargs[x_] = x

Unprotect[Rule];
Format[a_ -> b_, FortranForm] := SequenceForm[a, " = ", b];
Protect[Rule]


OptTimes[a__] :=
Block[ {p, const, var},
  p = Position[N[{a}], _Real, 1, Heads -> False];
  const = Times@@ {a}[[ Flatten[p] ]];
  If[ NumberQ[const], Return[Times[a]] ];
  var = Times@@ Delete[{a}, p];
  If[ NumberQ[var], Return[Times[a]] ];
  If[ MatchQ[const, _?Negative _.],
    -HoldForm[HoldForm[#1] #2]&[-const, var],
    HoldForm[HoldForm[#1] #2]&[const, var] ]
]


Options[WriteSummedExpr] = {PrintResult -> False}

WriteSummedExpr[hh_, var_ -> expr_, opt___Rule] :=
Block[ {SumOver, IndexRange, addto, pr, svar, si},
  loops = ToDoLoops[ {SeparateSums[expr]}//Flatten ];
  If[ loops[[1, 0]] === DoLoop,
    WriteString[hh, "\t", ToFortran[var], " = 0\n"];
    addto = True,
  (* else *)
    addto := (addto = True; False) ];
  SumOver[i_, r_] := (IndexRange[i] = r; 1);
  WriteDoLoops[hh, loops,
    WriteExpr[#1, var -> #2, addto, Optimize -> True]&];
  pr = PrintResult /. {opt} /. Options[WriteSummedExpr];
  If[ pr =!= False,
    svar = ToFortran[var];
    WriteString[hh,
      If[pr === IfDebug, "\n#ifdef DEBUG", ""] <>
      "\n\tprint *, '" <>
      StringReplace[svar, Cases[var,
        i_Symbol :> (si = ToFortran[i]; si -> "'," <> si <> ",'")]] <>
      " =', " <> svar <> 
      If[pr === IfDebug, "\n#endif\n\n", "\n\n"] ]
  ];
]


Attributes[SeparateSums] = {Listable}

SeparateSums[x_] := x /; FreeQ[x, SumOver]

SeparateSums[x_] := If[Head[#] === Plus, List@@ #, #]& @
  Collect[x, Cases[x, _SumOver, Infinity]//Union]


FindIndices[var_ -> _] := Union[Cases[var, _Symbol]]

FindIndices[t_Times] := Cases[t, SumOver[i_, _] -> i]

FindIndices[_] = {}

ToDoLoops[h_[li__], indices_:FindIndices] :=
Block[ {do},
  do[_] = {};
  Apply[(do[#1] = {do[#1], #2})&, {indices[#], #}&/@ {li}, 1];
  Cases[ DownValues[do],
    _[_[_[ind_List]], a_] :> DoLoop[ind, h@@ Flatten[a]] ]
]

ToDoLoops[x_, ___] := {x}


DoLoop[_[], {a__}] = a

DoLoop[_[], a_] = a


Attributes[WriteDoLoops] = {Listable}

WriteDoLoops[hh_, DoLoop[ind_, expr_], func_] := (
  WriteString[hh, {"\n\tdo ", ToFortran[#], " = 1, ",
    ToFortran[IndexRange[#]]}&/@ ind <> "\n"];
  WriteDoLoops[hh, expr, func];
  WriteString[hh, StringJoin[Table["\tenddo\n", {Length[ind]}]] ];
)

WriteDoLoops[hh_, expr_, func_] := func[hh, expr]



(* NumPrep pulls the calculation of several kinds of objects out of the
   actual computation of the form factor. That is, these objects are
   calculated first and stored in variables. After this, the form factors
   are computed. Because the replacement of more complex objects by
   scalar variables shortens the expressions, these objects are referred
   to as abbreviations.

   NumPrep treats these objects as abbreviations:
   1. variables replacing some function call: Pair, Eps, cint, iint,
   2. (true) abbreviations, i.e. products or sums of variables
      of type 1: Abb, AbbSum,
   3. matrix elements: Mat.

   Calculating these abbreviations in a clever way is key to a decent
   performance of the generated code. Therefore, NumPrep splits the
   abbreviations into these categories:
   1. objects that depend only on model constants
      -> subroutine abbr_const,
   2. objects that in addition depend on S
      -> subroutine abbr_s,
   3. objects that in addition depend on other phase-space variables
      (angles etc.)
      -> subroutine abbr_angle,
   4. objects that depend on the helicities
      -> subroutine abbr_hel.
   The final user-callable subroutine squared_me takes care to
   invoke these abbr_nnn subroutines only when necessary. *)

Category[ rul:_[_, def_] ] := {{}, {}, {}, rul} /;
  !FreeQ[def, helvars]

Category[ rul:_[_, def_] ] := {{}, {}, rul, {}} /;
  !FreeQ[def, T | U | Sf | Tf | Uf | Pair | Eps]

Category[ rul:_[_, def_] ] := {{}, rul, {}, {}} /;
  !FreeQ[def, S]

Category[ rul_ ] := {rul, {}, {}, {}}


Dependencies[l_] := {l}

Dependencies[f__, l_] :=
Block[ {deps = First/@ l, f2 = {f}, f3, pl = {}, p},
  f3 = Map[Last, f2, {2}];
  While[ True,
    p = Union[ Take[#, 2]&/@
      Position[f3, Alternatives@@ deps, {2, Infinity}, Heads -> False] ];
    If[ Length[p] === 0, Break[] ];
    pl = Union[pl, p];
    deps = Apply[f2[[##, 1]]&, p, 1];
  ];
  Append[
    Dependencies@@ Delete[f2, pl],
    Join[l, Apply[f2[[##]]&, pl, 1]] ]
]


DependenceCategories[ab_] :=
  OnePassOrder/@ Dependencies@@ Flatten/@ Transpose[Category/@ ab]


(* OnePassOrder orders the abbreviations such that the definition of
   each abbreviation comes before its first use. In Mathematica lingo,
   the effect of OnePassOrder can be shown as
	expr //. Abbreviations[] ===
	  expr /. OnePassOrder[Abbreviations[]]
   In principle, FormCalc should produce abbreviations which are already
   in that order. However, during the calculation of the fermionic ME new
   abbreviations may have been introduced. *)

OnePassOrder[li_] :=
Block[ {c = 0, l = Length[li], fi, Dep, Ticket, posmap},
  fi = Alternatives@@ First/@ li;
  Ticket[a_, b_] := (Dep[a] = b; Ticket[a, b] = ++c) /; FreeQ[b, Dep];
  posmap = Apply[Ticket[#1, #2 /. x:fi -> Dep[x]]&, li, 1];
  While[c < l, posmap = posmap];
  li[[ Level[Sort[MapIndexed[List, posmap]], {3}] ]]
]



(* NumPrep-generated code is written out in five parts:

				     |	handled by NumPrep-function
   ----------------------------------|-----------------------------------
   1. the code modules		     |	WriteCode (-> WriteCodeModule)
   2. the abbreviations		     |	WriteAbbr (-> WriteAbbrModule)
   3. the variable declarations      |	MakeCode (-> CommonDecl, VarDecl)
   4. the makefile		     |	MakeCode
   5. the squared_me subroutine	     |	MakeCode
   6. copy the driver files          |  MakeCode 
*)


VarDecl[_, _[]] = ""

VarDecl[decl_, vars_] :=
Block[ {llen = Infinity, dl = StringLength[decl], l, s},
  Fold[
    ( l = StringLength[s = ToFortran[#2]];
      {#1, If[(llen += l + 2) > 64, llen = dl + l; decl, ", "], s} )&,
		(* 64 = 70 - 8 for the tab + 1 for \t + 1 for \n *)
    "",
    vars ]
]


CommonDecl[_[], __] = {}

CommonDecl[vars_, type_, common_] :=
Block[ {v, pindex, phead},
  v = Select[DeleteCases[TakeVar/@ vars, _[0]],
    StringTake[ToString[#], 1] =!= "t" &];
  pindex = Position[v, _[_Symbol..], 1, Heads -> False];
  phead = Position[v, _[__], 1, Heads -> False];
  VarDecl["\n\t" <> type, MapAt[IndexRange/@ # &, v, pindex]] <>
    VarDecl["\n\tcommon " <> common, MapAt[Head, v, phead]] <> "\n"
]

TakeVar[v_ -> _] = v

TakeVar[v_] = v


SubroutineDecl[name_, cext_, iext_, expr_] :=
Block[ {ce, ie},
  ce = Select[cext, !FreeQ[expr, #]&];
  ie = Select[iext, !FreeQ[expr, #]&];
  "\
\tsubroutine " <> name <> "\n\
\timplicit character (a-s,u-z)\n\
\timplicit double complex (t)\n\n\
#include \"vars.h\"" <>
    VarDecl["\n\tdouble complex ", ce] <>
    VarDecl["\n\tinteger ", ie] <>
    VarDecl["\n\texternal ", Join[ce, ie]] <> "\n\n"
]


WriteAbbr[_[], _] = {}

WriteAbbr[abbr_, name_] :=
Block[ {ab, namex},
  ab = ToDoLoops[abbr];
  If[ LeafCount[ab] > $FileSize,
    MapIndexed[
      ( namex = name <> FromCharacterCode[#2[[1]] + 96];
        WriteAbbrModule[#1, namex] )&,
      SizeSplit[ab, $FileSize] ],
  (* else *)
    {WriteAbbrModule[ab, name]} ]
]

WriteAbbrModule[abbr_, name_] :=
Block[ {hh},
  Print["writing ", name, ".F"];
  hh = OpenWriteFortran[name];
  WriteString[hh, SubroutineDecl[name, {Pair, Eps}, {}, abbr]];
  WriteDoLoops[hh, abbr, WriteExpr];
  WriteString[hh, "\tend\n"];
  Close[hh];
  name
]


WriteCode[file_] :=
Block[ {name, namex, amp, coeffh},
  name = StringReplace[file, ".m" -> ""];
  amp = Get[indir <> file];
  mandel = Select[mandel, FreeQ[amp, #]&];
  amp = amp /. int:loopint[__] :> abbint[int] /.
    DEN[p_, m_] -> 1/(p - m);
  coeffh = If[ StringMatchQ[name, bornmask],
    amp = amp /. unusedtree; Ctree,
    amp = amp /. unusedloop; Cloop ];

  If[LeafCount[amp] > $FileSize && Head[amp] === Plus,
    MapIndexed[
      ( namex = name <> FromCharacterCode[#2[[1]] + 96];
        Indices[namex] = Indices[name];
        WriteCodeModule[#1, namex] )&,
      SizeSplit[amp, $FileSize] ],
  (* else *)
    {WriteCodeModule[amp, name]} ]
]

WriteCodeModule[amp_, name_] :=
Block[ {hh},
  Print["writing ", name, ".F"];
  hh = OpenWriteFortran[name];
  WriteString[hh, SubroutineDecl[name, {Pair, Eps}, {IndexDelta}, amp]];
  WriteExpr[hh, ToCoeff[amp], True, NewLine -> True, Optimize -> True];
  WriteString[hh, "\tend\n"];
  Close[hh];
  name
]


ToCoeff[x_] := (maxmat[coeffh] = {Mat[1]}; coeffh[1] -> x) /;
  FreeQ[x, Mat]

ToCoeff[p_Plus] := ToCoeff/@ List@@ p

ToCoeff[Mat[m_] x_] :=
  ( maxmat[coeffh] = Max[maxmat[coeffh], Level[m, {-2}]];
    Level[m, {-1}, coeffh] -> x )


	(* LoopComponents gives back e.g.
		1. {F[4], SUN[3]}
		2. {F[j1], SUN[j2]}
		3. Ctree[4, 3]
		4. Ctree[j1, j2]
		5. "\tdo j1 = 1, 4\n\tdo j2 = 1, 3\n"
		6. "\tenddo\n\tenddo\n"			*)

LoopComponents[arr_] := 0 /; maxmat[arr] === {}

LoopComponents[arr_] :=
  ReplacePart[
    Transpose[LoopVar/@ ToList[maxmat[arr]]],
    arr, {{3, 0}, {4, 0}} ]

LoopVar[mat:_[1]] := { mat, mat, 1, 1, "", "" }

LoopVar[mat_[max_]] :=
Block[ {var = Unique["j"]},
  { mat[max], mat[var], max, var,
    "\tdo " <> ToString[var] <> " = 1, " <> ToString[max] <> "\n",
    "\tenddo\n" }
]


GetFile[file_] := Flatten[Get/@ file]

DefNeeded[m_[i_, j_]] := (Needed[m[x_, y_] -> _] := x <= i && y <= j)


Options[MakeCode] = {
  Legs -> 4,
  AbbrFile -> "abbr",
  MatFile -> "mat*",
  BornMask -> "born*",
  RenConst -> "rc_SM_Ddim.F",
  Drivers -> "drivers"
}

MakeCode[idir_, odir_, mask_String:"*.m", opt___Rule] :=
Block[ {legs, abbr, mat, bornmask, rconst, drivers,
SumOver, Indices, IndexRange, Hel, Global`F, Global`SUN,
abbint, cints, iints, helvars, mandel,
indir, outdir, files, hh, n,
unusedtree, unusedloop, maxmat, Needed, ntree, ntree2, nloop,
codemod, abbrmod, Conjugate = dconjg},

  {legs, abbr, mat, bornmask, rconst, drivers} =
    {Legs, AbbrFile, MatFile, BornMask, RenConst, Drivers} /.
    {opt} /. Options[MakeCode];

  SumOver[index_, range_] := (
    IndexRange[index] = range;
    AppendTo[Indices[name], index];
    1
  );

	(* Hel[n] = 0, 1 (= +), 2 (= -) *)
  Hel[n_Integer] := Hel[n] = ToSymbol["Hel", n];

(* abbint introduces abbreviations for the loop integrals.
   They fall into two categories:
   1. A0, B0, B1, B00, B11, DB0, DB1, DB00, DB11 (cint..),
   2. C0i and D0i (iint..).
   The reason is that for the latter the LoopTools functions Cget and Dget
   can be used to compute all tensor coefficients at once (which is much
   more efficient). Unlike the other integrals, whose results are double
   complex numbers, Cget and Dget return an integer pointing into a cache
   array. We follow here the conventions of LoopTools, where the actual
   tensor coefficients are retrieved from the arrays Cval and Dval. *)

  abbint[ C0i[i_, args__] ] :=
  Block[ {uu = Unique["iint"], ind},
    ind = Select[Indices[name], !FreeQ[{args}, #]&];
    If[ Length[ind] =!= 0, uu = uu@@ ind ];
    iints = {iints, uu -> Cget[args]};
    abbint[C0i[id_, args]] = Cval[id, uu];
    Cval[i, uu]
  ];

  abbint[ D0i[i_, args__] ] :=
  Block[ {uu = Unique["iint"], ind},
    ind = Select[Indices[name], !FreeQ[{args}, #]&];
    If[ Length[ind] =!= 0, uu = uu@@ ind ];
    iints = {iints, uu -> Dget[args]};
    abbint[D0i[id_, args]] = Dval[id, uu];
    Dval[i, uu]
  ];

  abbint[ func_ ] :=
  Block[ {uu = Unique["cint"], ind},
    ind = Select[Indices[name], !FreeQ[func, #]&];
    If[ Length[ind] =!= 0, uu = uu@@ ind ];
    cints = {cints, uu -> func};
    abbint[func] = uu
  ];


  outdir = Check[SetDirectory[odir] <> "/", Abort[]];
  ResetDirectory[];

  indir = Check[SetDirectory[idir] <> "/", Abort[]];
  mat = FileNames[mat];
  abbr = FileNames[abbr];
  files = Complement[
    FileNames[Flatten[{mask, bornmask <> ".m"}]],
    abbr, mat ];
  mat = GetFile[mat] /. m_Mat :> ToArray/@ m;
  abbr = Select[GetFile[abbr], AtomQ[ #[[1]] ]&];
  ResetDirectory[];

  n = First/@ DeleteCases[mat, _ -> 0];
  unusedtree = Thread[
    Select[ Union[#[[1, 2]]&/@ mat], FreeQ[n, Mat[_, #]]& ] -> 0 ];
  unusedloop = Thread[
    Select[ Union[#[[1, 1]]&/@ mat], FreeQ[n, Mat[#, _]]& ] -> 0 ];
  maxmat[_] = {};

(* The most economical way to handle the vectors is to put all of
   them in one array and use indices into that array in calls to
   Pair and Eps. NumPrep imposes solely that this array be ordered
   such that
	vec(1) = k1	vec(5) = k2
	vec(2) = e1(0)	vec(6) = e2(0)
	vec(3) = e1(+)	...
	vec(4) = e1(-)
   How the array is called is completely up to the driver program
   which also has to supply the Pair and Eps routines. Complex
   conjugated vectors are referenced by negative indices, e.g.
   vec(-3) = e1(+)^*. *)

  abbr = abbr /. Conjugate[ep_e] -> -ep /.
    { k[j_Integer] -> 4 (j - 1) + 1,
	(* the spin reference vector of fermions is equivalent to the
	   longitudinal polarization vector of bosons *)
      s[j_Integer] -> 4 (j - 1) + 2,
      e[j_Integer] -> 4 (j - 1) + 2 + Hel[j] };

(* Part 1: the code modules *)

  cints = iints = Indices[_] = {};
  mandel = Select[MandelstamVars, FreeQ[mat, #]&];
  codemod = Flatten[WriteCode/@ files];
  mandel = Complement[MandelstamVars, mandel];

(* Part 2: the variable declarations *)

  Print["writing vars.h"];

  iints = Flatten[iints];
  cints = Flatten[cints];
  helvars = Alternatives@@ Hel/@ Range[legs];

  ntree = LoopComponents[Ctree];
  nloop = LoopComponents[Cloop];
  ntree2 = LoopComponents[ If[ntree === 0, Cloop, Ctree] ];
  If[ nloop === 0, nloop = {{}, {}, {}, {}, "", ""} ];
  n = MapThread[MatType,
    {ToList[Max[ ntree2[[1]], nloop[[1]] ]], ntree2[[1]]}];
  Scan[DefNeeded, n];

  hh = OpenWrite[outdir <> "vars.h"];

  WriteString[hh, "\
#include \"model.h\"\n\
#include \"looptools.h\"\n\
#include \"renconst.h\"\n" <>
    CommonDecl[Flatten[{S, mandel}], "double precision ", "/kinvars/ "] <>
    If[ FreeQ[abbr, Scale], "",
      CommonDecl[{Scale}, "double precision ", "/kinvars/ "] ] <>
    CommonDecl[helvars, "integer ", "/kinvars/ "] <>
    CommonDecl[DeleteCases[First/@ abbr, Scale],
      "double complex ", "/abbrev/ "] <>
    CommonDecl[cints, "double complex ", "/loopint/ "] <>
    CommonDecl[iints, "integer ", "/loopint/ "] <>
    CommonDecl[
      #[[1, 1, 1]]&/@ DownValues[IndexRange],
      "integer ", "/indices/ "] <>
    CommonDecl[
      Flatten[{ntree2[[3]], nloop[[3]], DeleteCases[n, 1]}]//Union,
      "double complex ", "/coeff/ " ]
  ];
  Close[hh];

(* Part 3: the abbreviations *)

  mat = Select[mat /. Mat -> MatType, Needed];
  abbr = DependenceCategories[ Join[abbr, mat, cints, iints] ];
  abbrmod = MapThread[WriteAbbr,
    {abbr, {"abbr_const", "abbr_s", "abbr_angle", "abbr_hel"}}];

(* Part 4: the makefile *)

  Print["writing GNUmakefile"];

  hh = OpenWrite[outdir <> "GNUmakefile"];

  WriteString[hh, "\
LTDIR = " <> $LoopToolsDir <> "\n\n\
LOOPTOOLS = -L$(LTDIR)lib -looptools\n\n\
EXTRALIBS = " <> $ExtraLibs <> "\n\n\
FC = " <> $F77 <> " -I. -I$(LTDIR)include -g\n\n\
OBJS =" <> ({" \\\n ", #, ".o"}&)/@ Flatten[{abbrmod, codemod}] <> "\n\n\
default: run\n\n\
renconst.o: " <> rconst <> " model.h\n\
\t$(FC) -c -o $@ $< # -DDEBUG\n\n\
$(OBJS): %.o: %.F model.h vars.h\n\
\t$(FC) -c $<\n\n\
clean:\n\
\trm -f *.o\n\n\
%:: %.F num.F model.h process.h squared_me.F $(OBJS) renconst.o\n\
\t$(FC) -o $@ $< renconst.o $(OBJS) $(LOOPTOOLS) $(EXTRALIBS)\n\n"];
  Close[hh];

(* Part 5: the control file squared_me.F *)

  Print["writing squared_me.F"];

  helvars = {#, # <> "from", # <> "to"}&/@ ToString/@ List@@ helvars;
  n = Flatten[Rest/@ helvars];

  hh = OpenWriteFortran["squared_me"];

  WriteString[hh, "\
\tsubroutine squared_me(treeres, loopres, Ecms2" <>
    VarDecl[",\n     +    ", n] <> ", reset)\n\
\timplicit integer (j)\n\
\timplicit character (a-i,k-z)\n\n\
\tdouble precision treeres, loopres\n\
\tdouble precision Ecms2" <> VarDecl["\n\tinteger ", n] <> "\n\
\tlogical reset\n\n\
#include \"vars.h\"\n\n" <>
    If[ Length[mandel] === 0, "", "\
\tdouble precision MomSquare\n\
\texternal MomSquare\n\n"
    ] <> "\
\tdouble complex c, m\n\
\tinteger Cptr0, Dptr0, Cptr1, Dptr1\n\
\tdouble precision prevS\n\
\tsave Cptr0, Dptr0, prevS\n\
\tdata prevS /-1/\n\n\
\tS = Ecms2\n"];

  WriteExpr[hh, # -> MandelstamDef[#]&/@ mandel];

  WriteString[hh, "\n\
\tif(reset .or. S .ne. prevS) then\n\
\t  if(reset .or. prevS .lt. 0) then\n\
\t    call setcachelast(Ccache, 0)\n\
\t    call setcachelast(Dcache, 0)\
" <> ({"\n\t    call ", #}&)/@ abbrmod[[1]] <> "\n\
\t    Cptr0 = getcachelast(Ccache)\n\
\t    Dptr0 = getcachelast(Dcache)\n\
\t    reset = .FALSE.\n\
\t  else\n\
\t    call setcachelast(Ccache, Cptr0)\n\
\t    call setcachelast(Dcache, Dptr0)\n\
\t  endif\
" <> ({"\n\t  call ", #}&)/@ abbrmod[[2]] <> "\n\
\t  Cptr1 = getcachelast(Ccache)\n\
\t  Dptr1 = getcachelast(Dcache)\n\
\t  prevS = S\n\
\tendif\n\
" <> ({"\n\tcall ", #}&)/@ abbrmod[[3]] <> "\n\n\
\ttreeres = 0\n\
\tloopres = 0\n\
" <> Apply[{"\n\tdo 1 ", #1, " = ", #2, ", ", #3}&, helvars, 1] <> "\n\
" <> ({"\n\tcall ", #}&)/@ abbrmod[[4]] <> "\n\n" <>
  If[ ntree === 0, "", ntree[[5]] <> "\
\t" <> ToFortran[ ntree[[4]] ] <> " = 0\n\
" <> ntree[[6]] <> "\n"] <> 
  If[ nloop[[1]] === {}, "", nloop[[5]] <> "\
\t" <> ToFortran[ nloop[[4]] ] <> " = 0\n\
" <> nloop[[6]] <> "\n" ]
  ];

  WriteDoLoops[hh, ToDoLoops[codemod, Indices],
    WriteString[#1, "\tcall " <> #2 <> "\n"]& ];

  WriteString[hh, "\n\
" <> ntree2[[5]] <> "\
\tc = dconjg(" <> ToFortran[ ntree2[[4]] ] <> ")\n\n" <>
  If[ ntree === 0, "", "\
\tm = 0\n\
" <> ntree[[5]] <> "\
\tm = m + " <>
  ToFortran[ ntree[[4]] *
    Inner[MatType, ntree[[2]], ntree2[[2]], Times] ] <> "\n\
" <> ntree[[6]] <> "\
\ttreeres = treeres + dble(c*m)\n\n" ] <>
  If[ nloop[[1]] === {}, "", "\
\tm = 0\n\
" <> nloop[[5]] <> "\
\tm = m + " <>
  ToFortran[ nloop[[4]] *
    Inner[MatType, nloop[[2]], ntree2[[2]], Times] ] <> "\n\
" <> nloop[[6]] <> "\
\tloopres = loopres + " <> If[ntree === 0, "", "2*"] <> "dble(c*m)\n\n" ] <>
  ntree2[[6]] <> "\n\
1\tcontinue\n\n\
\tcall setcachelast(Ccache, Cptr1)\n\
\tcall setcachelast(Dcache, Dptr1)\n\
\tend\n"];

  Close[hh];

(* Part 6: copy the driver files *)

  files = {};

  If[ FileType[drivers] === Directory,
    SetDirectory[drivers];
    CopyFile[#, outdir <> #]&/@ (files = FileNames[]);
    ResetDirectory[]
  ];

  SetDirectory[$NumPrepDir <> "drivers"];
  CopyFile[#, outdir <> #]&/@ Complement[FileNames[], files];
  ResetDirectory[];

  SetDirectory[$RenConstDir];
  CopyFile[#, outdir <> #]&/@
    Complement[FileNames[{rconst, "*.h"}], files];
  ResetDirectory[];
]

End[]


$F77::undef = "Warning: I don't know the correct f77 compiler switches for
your system. I'll use g77 this time, but if you have a native Fortran
compiler, you should update NumPrep.m."

(* IMPORTANT: you must choose here the correct compiler flags for your
   Fortran compiler, *in particular* (= at least) one to make the
   compiler disregard the stubborn 72-characters-per-line limit of
   Fortran. Practically all compilers have such a flag, only
   unfortunately, the naming of this flag varies widely across different
   computer systems.

   Choices for some systems are listed below in the $F77 = ... line. If
   you don't find your system among them, please look up the appropriate
   options in your f77 man page and insert them there. *)

$F77 := $F77 = Switch[ $GenerateCodeFor,
  "DEC-AXP", "f77 -O -extend_source -warn truncated_source",
  "HP-RISC", "fort77 -O2 +es +U77",
  "Solaris", "f77 -e",
  "PGF77",   "pgf77 -O2 -Mextend -g77libs",
  "Linux",   "g77 -O0",
  _,         Message[f77::undef]; "g77 -O0" ]

$GenerateCodeFor = "DEC-AXP"

$ExtraLibs = "-L$(CERN)/$(CERN_LEVEL)/lib -lpdflib -lmathlib -lpacklib"

$BlockSize = 700

$FileSize = 30 $BlockSize

(* this assumes that LoopTools has been installed in /usr/local/ *)
$LoopToolsDir = "/usr/local/"

(* for the standard LoopTools installation in a user's home directory use
$LoopToolsDir = "$(HOME)/LoopTools/$(HOSTTYPE)/" *)

$RenConstDir = $NumPrepDir <> "rconst/fortran/"

EndPackage[]


(* Model parameters which are defined using the parameter statement in
   Fortran (i.e. as numeric constants; see model.h) are given some
   numeric value here. Using this information, the OptTimes function can
   significantly optimize the generated Fortran code. The idea is to put
   everything that is known as constant at compile time in one place,
   i.e. rearrange products such that they are of the form (const)*(vars),
   then the compiler will usually collect all of these constants into one
   number. *)

Scan[ (N[#] = 1)&,
  { Alfa, Alfa2,
    MW, MW2, MZ, MZ2, SW2, CW, CW2,
    ME, ME2, MM, MM2, ML, ML2,
    MU, MU2, MC, MC2, (* MT, MT2, *)
    MD, MD2, MS, MS2, MB, MB2 } ]


Null

