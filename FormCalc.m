(*

This is FormCalc, Version 3.1
Copyright by Thomas Hahn 1996-2002
last modified 4 Feb 02 by Thomas Hahn

Release notes:

FormCalc is free software, but is not in the public domain.
Instead it is covered by the GNU library general public license.
In plain English this means:

1. We don't promise that this software works.
   (But if you find any bugs, please let us know!)

2. You can use this software for whatever you want.
   You don't have to pay us.

3. You may not pretend that you wrote this software.
   If you use it in a program, you must acknowledge
   somewhere in your publication that you've used
   our code.

If you're a lawyer, you can find the legal stuff at
http://www.fsf.org/copyleft/lgpl.html.

The user guide for this program can be found at
http://www.feynarts.de/formcalc.

If you find any bugs, or want to make suggestions, or
just write fan mail, address it to:
	Thomas Hahn
	Institute for Theoretical Physics
	University of Karlsruhe
	D-76128 Karlsruhe, Germany
	e-mail: hahn@feynarts.de

To join the FormCalc mailing list, send a mail (any text) to
	hahn-formcalc-subscribe@particle.uni-karlsruhe.de

Have fun!

*)

Print[""];
Print["FormCalc 3.1"];
Print["by Thomas Hahn"];
Print["last revised 4 Feb 02"];


(* symbols from the model files live in Global` *)

{ DiracMatrix, DiracSlash, ChiralityProjector,
  DiracSpinor, MajoranaSpinor, PolarizationVector,
  MetricTensor, FourVector, ScalarProduct, Lorentz,
  SUNT, SUNF, SUNTSum, Colour, Gluon }

(* symbols from FeynArts *)

Begin["FeynArts`"]

{ FeynAmp, FeynAmpList, Process,
  GraphID, Number, Insertions, Classes, Particles,
  PropagatorDenominator, FeynAmpDenominator,
  FourMomentum, Internal,
  Index, IndexDelta, IndexSum, SumOver, External,
  MatrixTrace, FermionChain, NonCommutative,
  CreateTopologies, ExcludeTopologies, Tadpoles,
  Paint, InsertFields, CreateFeynAmp }

End[]


(* these symbols are used in the model file, and hence live in Global` *)

RenConst::usage = "RenConst[rc] := ... defines the renormalization
constant rc."

SelfEnergy::usage = "SelfEnergy[from -> to, mass] calculates the
self-energy with incoming particle from and outgoing particle to, taken at
k^2 = mass^2."

DSelfEnergy::usage = "DSelfEnergy[from -> to, mass] calculates the
derivative with respect to k^2 of the self-energy with incoming particle
from and outgoing particle to, taken at k^2 = mass^2."

LVectorCoeff::usage = "LVectorCoeff[expr] returns the coefficient of
DiracChain[6, k[1]] (= k1slash omega_-) in expr."

RVectorCoeff::usage = "RVectorCoeff[expr] returns the coefficient of
DiracChain[7, k[1]] (= k1slash omega_+) in expr."

LScalarCoeff::usage = "LScalarCoeff[expr] returns the coefficient of
DiracChain[7] (= omega_-) in expr."

RScalarCoeff::usage = "RScalarCoeff[expr] returns the coefficient of
DiracChain[6] (= omega_+) in expr."

ReTilde::usage = "ReTilde[expr] takes the real part of loop integrals
occurring in expr."

ImTilde::usage = "ImTilde[expr] takes the imaginary part of loop integrals
occurring in expr."


If[ !MemberQ[$ContextPath, "FeynArts`"],
  PrependTo[$ContextPath, "FeynArts`"] ]

BeginPackage["FormCalc`"]

(* some internal symbols must be visible for ReadForm *)

{ ReadForm, ClearCache, FormSetup, FormSubst,
  abb, fme, sun, pave, sqrt2, d$, e$, q1 }


(* symbols appearing in the output *)

Amp::usage = "Amp[proc][expr1, expr2, ...] is the result of the
calculation of diagrams of the process proc. The result is divided into
parts expr1, expr2, ..., such that index sums (marked by SumOver) apply
to the whole of each part."

Den::usage = "Den[p2, m2] stands for 1/(p2 - m2). Note that in contrast
to PropagatorDenominator, p2 and m2 are the momentum and mass *squared*."

Delta::usage = "Delta[a, b] is the Kronecker delta with indices a and b."

DiracChain::usage = "DiracChain[objs] is a chain of Dirac matrices
contracted with the given objects. The integers 5, 6, and 7 appearing as
arguments denote gamma_5, (1 + gamma_5)/2, and (1 - gamma_5)/2,
respectively."

FeynArts`Spinor::usage = "Spinor[p, m, s] is a spinor with momentum p and
mass m, i.e. a solution of the Dirac equation (pslash + s m) Spinor[p, m,
s] = 0. On screen, particle spinors (s = 1) are printed as u[p, m],
antiparticle spinors (s = -1) as v[p, m]."

e::usage = "e[n] is the nth polarization vector."

k::usage = "k[n] is the nth momentum."

SUNN::usage = "SUNN specifies the N in SU(N), i.e. the number of colours."

FeynArts`S::usage = "S is the Mandelstam variable s. If k1 and k2 are the
incoming momenta, S = (k1 + k2)^2."

T::usage = "T is the Mandelstam variable t. If k1 denotes the first
incoming and k3 the first outgoing momentum, T = (k1 - k3)^2."

FeynArts`U::usage = "U is the Mandelstam variable u. If k1 denotes the
first incoming and k4 the second outgoing momentum, U = (k1 - k4)^2."

Sf::usage = "Sf is an extended Mandelstam variable for a 2 -> 3 reaction.
If k3 and k4 denote the first two outgoing momenta, Sf = (k3 + k4)^2."

Tf::usage = "Tf is an extended Mandelstam variable for a 2 -> 3 reaction.
If k2 denotes the second incoming and k4 the second outgoing momentum,
Tf = (k2 - k4)^2."

Uf::usage = "Uf is an extended Mandelstam variable for a 2 -> 3 reaction.
If k2 denotes the second incoming and k3 the first outgoing momentum,
Uf = (k2 - k3)^2."


(* CalcFeynAmp and its options *)

CalcFeynAmp::usage = "CalcFeynAmp[amps] calculates the Feynman amplitudes
given in amps. The resulting expression is broken up into categories which
are returned in a list."

CalcLevel::usage = "CalcLevel is an option of CalcFeynAmp. It specifies
the level (Classes or Particles) at which to calculate the amplitudes.
Automatic takes Classes level, if available, otherwise Particles."

Dimension::usage = "Dimension is an option of CalcFeynAmp. It specifies
the space-time dimension in which to perform the calculation and can take
the values D, where dimensional regularization is used, and 4, where
constrained differential renormalization is used. The latter method is
equivalent to dimensional reduction at the one-loop level."

OnShell::usage = "OnShell is an option of CalcFeynAmp. It specifies
whether FORM should be instructed to put the external particles on their
mass shell, i.e. to apply ki^2 = mi^2."

Mandelstam::usage = "Mandelstam is an option of CalcFeynAmp. It specifies
whether FORM should be instructed to introduce kinematical invariants,
like the Mandelstam variables for a 2 -> 2 process."

Transverse::usage = "Transverse is an option of CalcFeynAmp. It specifies
whether FORM should be instructed to apply the transversality relations
for polarization vectors (ei.ki = 0)."

Normalized::usage = "Normalized is an option of CalcFeynAmp. It specifies
whether FORM should be instructed to exploit the normalization of the
polarization vectors (ei.ei^* = -1)."

MomSimplify::usage = "MomSimplify is an option of CalcFeynAmp. It
specifies whether FORM should be instructed to try all possible
permutations of the momentum conservation equation on dot products and
fermion chains in order to find the shortest possible combination. This
might be slow, however."

VADecompose::usage = "VADecompose is an option of CalcFeynAmp. It
specifies whether the usual chiral fermion chains should be decomposed
into vector and axial-vector parts."

NoExpand::usage = "NoExpand is an option of CalcFeynAmp. NoExpand ->
{sym1, sym2, ...} specifies that sums containing any of sym1, sym2, ...
are not expanded during the FORM calculation."

AbbrScale::usage = "AbbrScale is an option of CalcFeynAmp. It specifies a
factor with which the abbreviations introduced by CalcFeynAmp are scaled
for every momentum they include. This way, the abbreviations can be made
dimensionless which is of advantage in some applications, e.g. in treating
resonances."

EditCode::usage = "EditCode is a debugging option of CalcFeynAmp and
HelicityME. It edits the temporary file passed to FORM using the editor
command in $Editor."

RetainFile::usage = "RetainFile is a debugging option of CalcFeynAmp and
HelicityME. When set to True, the temporary file containing the FORM
code is not removed after running FORM."


(* abbreviationing-related functions *)

Abbr::usage = "Abbr[] returns a list of all abbreviations introduced
so far."

OptimizeAbbr::usage = "OptimizeAbbr[abbr] optimizes the abbreviations
in abbr by eliminating common subexpressions."

Pair::usage = "Pair[a, b] represents the contraction of the two
four-vectors or Lorentz indices a and b."

Eps::usage = "Eps[a, b, c, d] represents -I times the antisymmetric
Levi-Civita tensor epsilon_{abcd}. The sign convention is
epsilon^{0123} = +1."

num::usage = "num is the function in which FORM collects the numerical
prefactor of a term."

kin::usage = "kin is the function in which FORM collects kinematical
factors."


(* miscellaneous functions *)

MandelstamSimplify::usage = "MandelstamSimplify attempts to simplify
expressions containing S, T, and U using the Mandelstam relation
S + T + U = sum of external masses squared, which holds for on-shell
2 -> 2 processes. For all other reactions, MandelstamSimplify is
equivalent to Simplify."

ClearProcess::usage = "ClearProcess[] is necessary to clear internal
definitions before calculating a process with a different kinematical
set-up."

DiagramType::usage = "DiagramType[diag] returns the number of
denominators not containing the integration momentum."

FermionicQ::usage = "FermionicQ[diag] gives True for a diagram containing
fermions and False otherwise."

Small::usage = "Small[sym] = 0 tells FormCalc to put sym = 0 except
when it appears in negative powers or in loop integrals."


(* FeynCalc compatibility functions *)

FeynCalcGet::usage = "FeynCalcGet[mask] reads files produced with
FeynCalc. mask is taken as input to the Mathematica function FileNames, so
it might be FeynCalcGet[\"file.m\"] or FeynCalcGet[\"*.m\"] or
FeynCalcGet[\"*.m\", \"~/feyncalcfiles\"]."

FeynCalcPut::usage = "FeynCalcPut[expr, file] writes expr to file in
FeynCalc format."

C0::usage = "C0[p1, p2, p1p2, m1, m2, m3] is the scalar three-point
Passarino-Veltman function as used by FeynCalc. It is converted to C0i in
FormCalc."

D0::usage = "D0[p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the scalar
four-point Passarino-Veltman function as used by FeynCalc. It is converted
to D0i in FormCalc."

PaVe::usage = "PaVe[ind, {pi}, {mi}] is the generalized Passarino-Veltman
function as used by FeynCalc. It is converted to C0i or D0i in FormCalc."


(* one-loop integrals *)

A0::usage = "A0[m] is the one-point one-loop scalar integral. m is the
mass squared."

B0::usage = "B0[p, m1, m2] is the two-point one-loop scalar integral. p is
the external momentum squared and m1 and m2 are the masses squared."

B1::usage = "B1[p, m1, m2] is the coefficient of k_\\mu in the two-point
one-loop tensor integral B_\\mu. p is the external momentum k squared and
m1 and m2 are the masses squared."

B00::usage = "B00[p, m1, m2] is the coefficient of g_{\\mu\\nu} in the
two-point one-loop tensor integral B_{\\mu\\nu}. p is the external
momentum squared and m1 and m2 are the masses squared."

B11::usage = "B11[p, m1, m2] is the coefficient of k_\\mu k_\\nu in the
two-point one-loop tensor integral B_{\\mu\\nu}. p is the external
momentum k squared and m1 and m2 are the masses squared."

DB0::usage = "DB0[p, m1, m2] is the derivative of B0[p, m1, m2] with
respect to p."

DB1::usage = "DB1[p, m1, m2] is the derivative of B1[p, m1, m2] with
respect to p."

DB00::usage = "DB00[p, m1, m2] is the derivative of B00[p, m1, m2] with
respect to p."

DB11::usage = "DB11[p, m1, m2] is the derivative of B11[p, m1, m2] with
respect to p."

C0i::usage = "C0i[id, p1, p2, p1p2, m1, m2, m3] is the generic three-point
loop integral which includes both scalar and tensor coefficients,
specified by id. For example, C0i[cc0, ...] is the scalar function C_0,
C0i[cc112, ...] the tensor coefficient function C_{112} etc. The arguments
of C0i are p1 = k1^2, p2 = k2^2, p1p2 = (k1 + k2)^2, where k1...k3 are the
external momenta, and m1...m3 are the masses squared."

D0i::usage = "D0i[id, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the
generic four-point loop integral which includes both scalar and tensor
coefficients, specified by id. For example, D0i[dd0, ...] is the scalar
function D_0, D0i[dd1233, ...] the tensor function D_{1233} etc. The
arguments of D0i are p1 = k1^2, p2 = k2^2, p3 = k3^2, p4 = k4^2, p1p2 =
(k1 + k2)^2, p2p3 = (k2 + k3)^2, where k1...k4 are the external momenta,
and m1...m4 are the masses squared."

{ A0i, B0i, E0i }

{ cc0, cc1, cc2, cc00, cc11, cc12, cc22, cc001, cc002, cc111, cc112,
  cc122, cc222,
  dd0, dd1, dd2, dd3, dd00, dd11, dd12, dd13, dd22, dd23, dd33, dd001,
  dd002, dd003, dd111, dd112, dd113, dd122, dd123, dd133, dd222, dd223,
  dd233, dd333, dd0000, dd0011, dd0012, dd0013, dd0022, dd0023, dd0033,
  dd1111, dd1112, dd1113, dd1122, dd1123, dd1133, dd1222, dd1223,
  dd1233, dd1333, dd2222, dd2223, dd2233, dd2333, dd3333 }


(* finiteness checks *)

UVDivergentPart::usage = "UVDivergentPart[expr] returns expr with all loop
integrals replaced by their UV-divergent part. The divergence itself is
denoted by Divergence, so to assert that expr is UV-finite, one can check
if FreeQ[Expand[UVDivergentPart[expr], Divergence], Divergence] is True."

Divergence::usage = "Divergence stands for the dimensionally regularized
divergence 2/(4 - D) of loop integrals. It is used by the function
UVDivergentPart."


(* matrix elements *)

HelicityME::usage = "HelicityME[plain, conj] calculates the helicity
matrix elements for all combinations of spinor chains that appear in the
expression (plain conj^*). Terms of this kind arise in the calculation of
the squared matrix element, where typically plain is the one-loop result
and conj the Born expression. The arguments do not necessarily have to
be amplitudes since they are only used to determine which spinor chains
to select from the abbreviations. The symbol All can be used to select
all spinor chains currently defined in the abbreviations."

ColourME::usage = "ColourME[plain, conj] calculates the colour matrix
elements. ColourME is very similar to HelicityME, except that it computes
the matrix elements for SU(N) objects, not for spinor chains."

All::usage = "All as an argument of HelicityME and ColourME indicates
that all spinor chains or SUNT objects currently defined in the
abbreviations should be used instead of just those appearing in the
argument."

AbbrToUse::usage = "AbbrToUse is an option of HelicityME. It specifies
which abbreviations are used to calculate the helicity matrix elements."

Hel::usage = "Hel[i] is the helicity of the ith external particle. It can
take the values +1, 0, -1, where 0 stands for an unpolarized particle."

s::usage = "s[n] is the nth helicity reference vector."

Mat::usage = "Mat[Fi SUNi] is a matrix element in an amplitude, i.e. an
amplitude is a linear combination of Mat objects.\n
Mat[Fi, Fj] appears in the squared matrix element and stands for the
product of the two arguments, Fi Fj^*. Such expressions are calculated
by HelicityME and ColourME."

Lor::usage = "Lor[n] is a contracted Lorentz index in a product of Dirac
chains."

SquaredME::usage = "SquaredME[plain, conj] returns the matrix element
(plain Conjugate[conj)]. This performs a nontrivial task only for
fermionic amplitudes: the product of two fermionic amplitudes\n
    M1 = a1 F1 + a2 F2 + ... and\n
    M2 = b1 F1 + b2 F2 + ... is returned as\n
    M1 M2^* = a1 b1^* Mat[F1, F1] + a2 b1^* Mat[F2, F1] + ...\n
The special case of plain === conj can be written as SquaredME[plain]
which is of course equivalent to SquaredME[plain, plain]."

RealQ::usage = "RealQ[sym] is True if sym represents a real quantity
which means in particular that Conjugate[sym] = sym."


(* writing out Fortran code *)

WriteSquaredME::usage = "WriteSquaredME[tree, loop, me, abbr, ..., dir]
writes out Fortran code to compute the squared matrix element for a
process whose tree-level and one-loop contributions are given in the first
and second argument, respectively. All further arguments except the last
specify the necessary matrix elements and abbreviations. The last argument
dir finally gives the path to write the generated code to."

SymbolPrefix::usage = "SymbolPrefix is an option of WriteSquaredME. It
specifies a string which is prepended to global symbols in the generated
Fortran code to prevent collision of names when several processes are
linked together."

RenConstFile::usage = "RenConstFile is an option of WriteSquaredME and
WriteRenConst. It specifies the name of the Fortran program used to
calculate the renormalization constants."

Drivers::usage = "Drivers is an option of WriteSquaredME. Drivers points
to a directory containing customized versions of the driver programs for
running the generated Fortran code. This directory need not contain all
driver programs: files not contained therein are taken from the default
directory $DriversDir."


(* renormalization constants *)

CalcRenConst::usage = "CalcRenConst[expr] calculates the renormalization
constants appearing in expr."

WriteRenConst::usage = "WriteRenConst[expr, dir] calculates the
renormalization constants appearing in expr and generates a Fortran
program from the results. The resulting files (the Fortran program itself
and the corresponding declarations) are written to the directory dir.  
The names of the files are determined by the RenConstFile option."

InsertFieldsHook::usage = "InsertFieldsHook[tops, proc] is the function
called by SelfEnergy and DSelfEnergy to insert fields into the topologies
tops for the process proc. It is normally equivalent to InsertFields, but
may be redefined to change the diagram content of certain self-energies."

ClearSE::usage = "ClearSE[] clears the internal definitions of
already-calculated self-energies."

$PaintSE = "$PaintSE determines whether SelfEnergy paints the diagrams it
generates to compute the self-energies. $PaintSE can be True, False, or a
string which indicates that the output should be saved in a PostScript
file instead of being displayed on screen, and is prepended to the
filename."


(* low-level Fortran output functions *)

ToFortran::usage = "ToFortran[expr] returns the Fortran form of expr as
a string."

OpenFortran::usage = "OpenFortran[file] opens file for writing in
Fortran format."

TimeStamp::usage = "TimeStamp[] returns a string with the current date
and time."

WriteExpr::usage = "WriteExpr[file, {var1 -> expr1, var2 -> expr2, ...},
addto] writes a list of variable assignments in Fortran format to
file. If an expression's leaf count is larger than $BlockSize, the
expression is cut into blocks, viz.\n\
\tvar = part1\n\
\tvar = var + part2\n\
\t...\n
addto is an optional argument: if True, the variables var1, var2, ... are
always added to, i.e. \"var = var + expr\"."

Newline::usage = "Newline is an option of WriteExpr. If set to True,
WriteExpr writes out an empty line (\\n) after each statement."

Optimize::usage = "Optimize is an option of WriteExpr. With Optimize ->
True, WriteExpr introduces variables for subexpressions which are used
more than once."

WriteSummedExpr::usage = "WriteSummedExpr[file, var -> exprlist] writes
Fortran code to file which computes each expression in exprlist, performs
the necessary index summations, and stores the sum in var. The members of
exprlist must be arranged such that their index sums (marked by SumOver)
always apply to the whole expression."

SplitSums::usage = "SplitSums[expr] splits expr into a list of expressions
such that index sums (marked by SumOver) always apply to the whole of each
part."

ToDoLoops::usage = "ToDoLoops[list, ifunc] splits list into patches which
must be summed over the same set of indices. ifunc is an optional
argument: ifunc[expr] must return the indices occurring in expr."

DoLoop::usage = "DoLoop[ind, expr] is a symbol introduced by ToDoLoops
indicating that expr is to be summed over the set of indices ind."

WriteDoLoops::usage = "WriteDoLoops[file, list, func] writes the do-loops
occurring in list to file. For writing the actual expressions, func[file,
item] is called."

Dim::usage = "Dim[i] returns the highest value the index i takes on."

OnePassOrder::usage = "OnePassOrder[r] orders a list of interdependent
rules such that the definition of each item (item -> ...) comes before its
use in the right-hand sides of other rules."

SubroutineDecl::usage = "SubroutineDecl[name, cfuncs, ifuncs, expr]
returns a string with the declaration of the Fortran subroutine name. The
members of cfuncs and ifuncs, if they appear in expr, are declared
respectively as double complex and integer external functions."

VarDecl::usage = "VarDecl[\"\\n\\tTYPE \", vars] returns a string with the
declaration of vars as variables of type TYPE in Fortran."

CommonDecl::usage = "CommonDecl[\"\\n\\tTYPE \", vars, com] returns a
string with the declaration of vars as variables of type TYPE and members
of the common block com in Fortran."


(* symbols used in the Fortran code *)

Ctree::usage = "Ctree[Fi] is the ith form factor (the coefficient of Fi)
of the tree-level amplitude."

Cloop::usage = "Cloop[Fi] is the ith form factor (the coefficient of Fi)
of the one-loop amplitude."

MomSquare::usage = "MomSquare[a, b] is a Fortran function which computes
the square of the difference of the momenta indexed by a and b, i.e.
Pair[k[a] - k[b], k[a] - k[b]]."

dconjg::usage = "dconjg[z] takes the complex conjugate of z in Fortran."

dble::usage = "dble[z] takes the real part of z in Fortran."

dimag::usage = "dimag[z] takes the imaginary part of z in Fortran."

Cval::usage = "Cval is the array containing the cached three-point
integrals in LoopTools."

Cget::usage = "Cget computes all three-point coefficients in LoopTools."

Dval::usage = "Dval is the array containing the cached four-point
integrals in LoopTools."

Dget::usage = "Dget computes all four-point coefficients in LoopTools."


(* system variables *)

$Editor::usage = "$Editor specifies the editor command line used in
debugging FORM code."

$FormCalc::usage = "$FormCalc contains the version number of FormCalc."

$FormCalcDir::usage = "$FormCalcDir points to the directory where
FormCalc lives."

$Platform::usage = "$Platform is a string that identifies the platform
FormCalc is running on. It is the value of the environment variable
HOSTTYPE in the Bourne shell (sh) and is used to distinguish the
executables of different platforms."

$FormCmd::usage = "$FormCmd gives the name of the actual FORM executable.
It may contain a path."

$DriversDir::usage = "$DriversDir is the path where the driver programs
for the generated Fortran code are located."

$BlockSize::usage = "$BlockSize is the maximum LeafCount a single Fortran
statement written out by WriteExpr may have. Any expression with
LeafCount > $BlockSize will be chopped up before being written to the
Fortran file."

$FileSize::usage = "$FileSize gives the maximum LeafCount the expressions
in a single Fortran file may have. If the expressions grow larger than
$FileSize, the file is split into several pieces."


Begin["`Private`"]

$ContextPath = {"FormCalc`", "FeynArts`", "Global`", "System`"}


FullFileName[file_, path_] :=
Block[ {full},
  Catch[
    Scan[
      If[ FileType[full = ToFileName[#, file]] === File, Throw[full] ] &,
      path ];
    file ]
]

$FormCalc = 3.1

$FormCalcDir =
  SetDirectory[DirectoryName[ FullFileName[$Input, $Path] ]]

ResetDirectory[]


$Platform = Environment["HOSTTYPE"]

If[ Head[$Platform] =!= String, $Platform = "" ]

Check[
  Install[ToFileName[$FormCalcDir, "ReadForm_" <> $Platform]],
  ReadForm::notcompiled = "\
The ReadForm executable could not be installed.
Did you run the compile script first?";
  Message[ReadForm::notcompiled];
  Abort[] ]


$NumberMarks = False

Off[General::spell1, General::spell, CopyFile::filex]


(* generic functions *)

General::noopt = "Warning: `2` is not a valid option of `1`."

ParseOpt[opt___, func_] :=
Block[ {names = First/@ Options[func]},
  Message[func::noopt, func, #]&/@
    Complement[First/@ {opt}, names];
  names /. {opt} /. Options[func]
]


ToSymbol[x__] := ToExpression[ StringJoin[ToString/@ Flatten[{x}]] ]


DiagramType[a_FeynAmp] := Count[a[[3]] /. _FeynAmpDenominator -> 1,
  _PropagatorDenominator, Infinity]


FermionicQ[a_FeynAmp] := !FreeQ[a[[3]], FermionChain | MatrixTrace]


(* preparations for FORM *)

Attributes[delta] = Attributes[Delta] = {Orderless}

delta[c1:Index[Colour, _], c2_] := SUNT[c1, c2]

delta[g1:Index[Gluon, _], g2_] := 2 SUNT[g1, g2, 0, 0]

delta[x__] := Delta[x]


sum[expr_, {i_, r_}] := r expr /; FreeQ[expr, i]

sum[expr_, {i_, r_}] :=
Block[ {dummy = Unique[ToString[i]]},
  FormIndices = {FormIndices, dummy};
  (expr /. i -> dummy) SumOver[dummy, r]
]


momlist = Array[ToSymbol["Global`k", #]&, 8];
epslist = Array[ToSymbol["Global`e", #]&, 8];
epsclist = Array[ToSymbol["Global`e", #, "c"]&, 8]

MapThread[
  (pol[_, #1, li_] := #2[li]; Conjugate[#2[li_]] ^= #3[li])&,
  {momlist, epslist, epsclist} ]

Conjugate[pol] ^= pol

vector[p_Plus, li_] := vector[#, li]&/@ p

vector[n_. p_Symbol, li_] := n p[li]


Attributes[scalar] = {Orderless}

scalar[0, _] = 0

scalar[p_Plus, x_] := scalar[#, x]&/@ p

scalar[n1_. p1_Symbol, n2_. p2_Symbol] := n1 n2 p1.p2


prop[0, m_Symbol] = -1/m^2

prop[p_, m_] := prop[-p, m] /; !FreeQ[p, -q1]

prop[p_, m_] := Den[p, m^2]

Attributes[loop] = {Orderless}

loop[a___, d_[p_, m1_], _[p_, m2_], b___] :=
  (loop[a, d[p, m1], b] - loop[a, d[p, m2], b])/Factor[m1 - m2]

loop[_[q1, m_], d___] := I Pi^2 *
  {A0i, B0i, C0i, D0i, E0i}[[ Length[{d}] + 1 ]]@@
    Flatten[{First/@ {d} - q1, m, Last/@ {d}}]


slash[p_Plus] := slash/@ p

slash[n_. p_Symbol] := n ga[p]

noncomm[p_Plus] := noncomm/@ p

noncomm[g_] := g ga[] /; FreeQ[g, ga | Spinor]

noncomm[g_] = g

noncomm[g__] := NonCommutativeMultiply[g]


DeclareProcess[proc_,
  flags:{onshell_, mandel_, transv_, norm_, momsimp_}] :=
Block[ {moms, signs, masses, legs, eps, epsc,
s, ki, kj, kiki = {}, kikj = {}, eiki = {}, eiei = {}},

  CurrentProcess = proc;
  CurrentFlags = flags;
  FormSymbols = FormProcs = MomRules = {};
  Clear[STUSum];

  signs = Flatten[{1&/@ proc[[1]], -1&/@ proc[[2]]}];
  masses = Level[proc, {2}]^2;
  legs = Length[masses];
  moms = Take[momlist, legs];
  eps = Take[epslist, legs];
  epsc = Take[epsclist, legs];

  FormVectors = Flatten[{eps, epsc, moms}];
  FromFormRules = MapThread[Rule,
    { FormVectors,
      Flatten[{{#, Conjugate[#]}&[Array[e, legs]], Array[k, legs]}] }];
  PrependTo[FormVectors, q1];
  Block[ {Set},
    (FormExec[cmd_] := Block[#, ReadForm[cmd]])& @
      Flatten[{Apply[Set, FromFormRules, 1], Dot = p$, sqrt2 = Sqrt[2]}]
  ];

  If[ onshell,
    kiki = MapThread[#1.#1 -> #2 &, {moms, FormSymbols = masses}] ];

  If[ mandel,
    Switch[ Length[moms],
      3,
        kikj = Array[
          Dot@@ Drop[moms, {#}] -> Plus@@
            (-Times@@ Drop[signs, {#}]/2 MapAt[-#&, #.#&/@ moms, #])&,
          3 ] /. kiki,
      4 | 5,
        kikj = If[ Length[moms] === 5, {S, T, U, Sf, Tf, Uf},
          If[onshell, STUSum = Plus@@ masses /.
            s_Symbol :> Small[s] /. Small -> Identity];
          {S, T, U, S, T, U} ];
        FormSymbols = {FormSymbols, kikj};
        kikj = MapThread[
          ( ki = moms[[#2]]; kj = moms[[#3]];
            s = signs[[#2]] signs[[#3]]/2;
            ki.kj -> s #1 - s ki.ki - s kj.kj )&,
          {kikj, {1, 1, 1, 3, 2, 2}, {2, 3, 4, 4, 4, 3}} ] /. kiki;
        If[ Length[moms] === 5,
          s = -signs[[-1]] Drop[signs, -1];
          ki = Drop[moms, -1];
          kikj = {kikj, Block[ {Dot},
            Attributes[Dot] = {Orderless, Listable};
            #1.k5 -> Plus@@ (s #1.ki /. kiki /. kikj)&/@ ki ]} ]
    ]
  ];

  ki = {""};
  Switch[ Length[moms],
    0 | 1,
      Null,
    2,
      MomRules = moms[[2]] -> moms[[1]];
      moms[[2]] = moms[[1]],
    _,
      kj = signs moms;
      FormProcs = {"#define MomSum \"", ToString[Plus@@ kj], "\"\n"};
      kj = Array[{-Plus@@ Drop[kj, {#}], kj[[#]]}&, Length[kj]];
      kj = Join[kj, -kj];
      MomRules = Apply[Rule, kj, 1];
      If[ momsimp,
        ki = Flatten[{ki, Cases[kj, {p2_, p1_Symbol} :>
          "id " <> ToString[p1] <> " = " <> ToString[p2] <> ";"]}];
        If[!(transv || onshell || mandel), ki = Drop[ki, {2, -2}]]
      ]
  ];

  If[ transv, eiki =
    MapThread[Dot[##] -> 0 &, Flatten/@ {{eps, epsc}, {moms, moms}}] ];

  If[ norm, eiei = MapThread[Dot[##] -> -1 &, {eps, epsc}] ];

  s = If[onshell, 1, Length[ki]];
  FormProcs = FormProcs <>
    "#define Legs \"" <> ToString[legs] <>
    "\"\n\n#procedure eiki\n" <>
    Id[eiki] <>
    "#endprocedure\n\n#procedure eiei\n" <>
    Id[eiei] <>
    "#endprocedure\n\n#procedure kikj\n" <>
    Id[kiki, kikj] <>
    "#endprocedure\n\n#procedure MomConserv(proc)\n" <>
    ({"#call `proc'(", #, callkikj[--s]}&)/@ ki <>
    "#endprocedure\n\n";
] /; proc =!= CurrentProcess || flags =!= CurrentFlags


callkikj[0] = ",#call kikj);\n"

callkikj[_] = ",);\n"

Id[eqs__] := Apply[{"id ", ToString[#1//InputForm], " = ",
  ToString[#2//InputForm], ";\n"}&, Flatten[{eqs}], 1]

CurrentProcess = Sequence[]


GList = Array[ToSymbol["Global`G", #]&, 30]

lhs[gr_ -> _] = gr

lhs[gr_] = gr

LevelSelect[lev_][id_, _, amp_, coup_ -> ins_] :=
Block[ {gen = amp, ins2, t, u, p, r, c = 0,
coupvar = Take[GList, Length[coup]]},
  ins2 = Apply[TrivialSums,
    lhs/@ Cases[{ins}, Insertions[lev][rulz__] :> rulz, Infinity], 1];
  t = Transpose[ins2];
  u = Union/@ t;
  p = Position[u, {_}, 1, Heads -> False];
  Apply[(coupvar[[#]] = u[[#, 1]])&, p, 1];
  t = ReplacePart[t, 0, p];
  While[ Length[t] > 1,
    u = t[[1]];
    t = Rest[t];
    ++c;
    If[ u =!= 0,
      If[Length[r = Position[t, -u, 1, Heads -> False]] =!= 0,
        t = ReplacePart[t, 0, r];
        coupvar = ReplacePart[coupvar, -coupvar[[c]], r += c];
        p = Join[p, r] ];
      If[Length[r = Position[t, u, 1, Heads -> False]] =!= 0,
        t = ReplacePart[t, 0, r];
        coupvar = ReplacePart[coupvar, coupvar[[c]], r += c];
        p = Join[p, r] ]
    ]
  ];
  r = amp /. Thread[coup -> coupvar];
  If[ Length[coupvar = Delete[coupvar, p]] === 0,
    {AmpID[id], Length[ins2] r},
    {AmpID[id], r, coupvar -> (Delete[#, p]&)/@ ins2} ]
]

LevelSelect[_][id_, _, amp_] :=
  Flatten[{AmpID[Select[id, FreeQ[#, Number]&]], TrivialSums[amp]}]


AmpID[g_] :=
Block[ {s, c},
  s = StringJoin@@
    Apply[{StringTake[ToString[#1], 1], ToString[#2]}&, g, 1];
  If[ (c = ++uniq[s]) === 96, s, s <> FromCharacterCode[c] ]
]


TrivialSums[ins___, r_] := {ins, r} /; FreeQ[r, SumOver]

TrivialSums[ins___, r_] :=
Block[ {test = {ins, r /. _SumOver -> 1, gen}},
  {ins, r /. SumOver -> CarryOut}
]

CarryOut[i_, r__] := SumOver[i, r] /; !FreeQ[test, i]

CarryOut[_, v_] = v

CarryOut[_, v_, _] := Sqrt[v]


InsList = Array[ToSymbol["Global`Ins", #]&, 500]

SUNobjs = SUNT | SUNTSum | SUNF


(* the main function CalcFeynAmp *)

Options[CalcFeynAmp] = {
  CalcLevel -> Automatic,
  Dimension -> D,
  OnShell -> True,
  Mandelstam -> True,
  Transverse -> True,
  Normalized -> True,
  MomSimplify -> True,
  VADecompose -> False,
  NoExpand -> {},
  AbbrScale -> 1,
  EditCode -> False,
  RetainFile -> False }

CalcFeynAmp::syntax = "Wrong syntax: CalcFeynAmp expects FeynAmpList
objects as arguments."

CalcFeynAmp::incomp = "Calculation of incompatible process(es) attempted.
If you want to calculate a new process, run ClearProcess[] first."

CalcFeynAmp[fal:FeynAmpList[__][___].., opt___Rule] :=
Block[ {lev, dim, onshell, mandel, transv, norm, momsimp, va, edit,
retain, procs, proc, uniq, FourMomentum, c = 0,
Index, FormIndices = {}, Dim, deny, hide, smalls, temp, hh,
amps, ampnoins, ampins, allins, ins, iabbr, insgen = {},
traces = False, ampsum = 0, Global`NoIns, Global`Result, res},

  procs = Cases[Head/@ {fal}, _[Process, p_] -> p, {2}];
  proc = Union[Apply[#3 &, procs, {3}], {CurrentProcess}];
  If[ Length[proc] =!= 1, Message[CalcFeynAmp::incomp]; Abort[] ];

  {lev, dim, onshell, mandel, transv, norm, momsimp,
    va, noexp, scale, edit, retain} = ParseOpt[opt, CalcFeynAmp];
  lev = Flatten[{lev}][[-1]] /.
    Automatic :> If[FreeQ[{fal}, Classes], Particles, Classes];

  DeclareProcess[proc[[1]],
    TrueQ/@ {onshell, mandel, transv, norm, momsimp}];
  Apply[(#2 = momlist[[++c]])&, procs[[1]], {2}];
  FourMomentum[Internal, 1] = q1;

  uniq[_] = 95;
  amps = Apply[LevelSelect[lev], Level[Hold[fal], {2}], 1] /.
    PolarizationVector -> pol /.
    MomRules /. {
      _MetricTensor^2 -> dim,
      MetricTensor -> "d_",
      ScalarProduct -> scalar,
      IndexDelta -> delta,
      IndexSum -> sum,
      PropagatorDenominator -> prop,
      FeynAmpDenominator -> loop,
      ChiralityProjector[c_] :> ga[(13 - c)/2],
      DiracMatrix -> ga,
      DiracSlash -> slash,
      (DiracSpinor | MajoranaSpinor)[(s_Integer:1) p_, m_] ->
        Spinor[p, m, s],
      NonCommutative -> noncomm,
      FourVector -> vector };

  amps = DeleteCases[amps, {_, 0, ___}];
  If[ Length[amps] === 0, Return[Amp[CurrentProcess][0]] ];

  Index[type_, n_] :=
  Block[ {i = ToSymbol[StringTake[ToString[type], 3], n]},
    FormIndices = {FormIndices, i};
    If[type === Lorentz, Dim[i] = i];
    Index[type, n] = i
  ];

  (amps = amps /. {2^(1/2) -> sqrt2, 2^(-1/2) -> 1/sqrt2,
            Complex[a_, b_] -> a + "i_" b}) /.
    SumOver[x1_, x2_, ___] :> (Dim[x1] = x1 == x2; 1);
  Dim[x_] = x == 0;

  FormIndices = Union[Flatten[ {FormIndices,
    Cases[amps, SUNobjs[s__] :> Cases[{s}, _Symbol], Infinity]} ]];
	(*      ^^^^^^^ possibly some colour or gluon
	   indices have been fixed by the user in FeynArts;
	   these would by default be declared as symbols
	   and not be seen by FORM *)

  noexp = Alternatives@@ Flatten[{noexp}];
  If[ Length[noexp] =!= 0,
    deny = Alternatives@@ Flatten[{FormVectors, FormFunc}];
    c = 0;
    hide[p_] := (hide[p] = -hide[-p]) /; MatchQ[p[[1]], -_];
    hide[p_] := hide[p] = ToSymbol["nx", ++c];
    amps = amps /. p_Plus :> hide[p] /; !FreeQ[p, noexp] && FreeQ[p, deny]
  ];

  ampnoins = Select[amps, Length[#] === 2 &];
  ampins = Select[amps, Length[#] === 3 &];
  allins = Union[Flatten[#[[-1, 2]]&/@ ampins]];
  ins = DeleteCases[allins, _Symbol];
  iabbr = Take[InsList, Length[ins]];
  ins = Thread[ins -> iabbr];

  OpenFormTemp;
  If[ Head[dim] === Symbol, WriteString[hh, "s ", dim, ";\n"] ];
  WriteString[hh, "d ", dim, ";\n"];
  DeclareVars[amps, SUNN];
  WriteString[hh, "f Spinor;\n\n.global\n\n"];

  Print["> ", Length[ampins], " amplitudes with insertions"];
  Print["> ", Length[ampnoins], " amplitudes without insertions"];

  If[ Length[ampnoins] =!= 0,
    Apply[FormWrite, ampnoins, 1];
    Write[hh, "g ", Global`NoIns == ampsum, ";"];
    WriteString[hh, "\n"];
    ampsum = Global`NoIns ];
  Apply[FormWrite, ampins, 1];

  WriteString[hh, ".sort\n\n" <>
    If[traces, "trace4,1;\n\n", ""] <>
    FormProcs <> "\
#procedure TrivialSubst\n\
id sqrt2^2 = 2;\n\
id sqrt2^-2 = 1/2;\n" <>
    FormSubst <>
    smalls <> "\
#endprocedure\n\n\
#procedure Insertions\n"];

  If[ Length[ampins] =!= 0,
    WriteString[hh, FormDecl["s ", iabbr] <> ".global\n\n"];
    Apply[
      ( Write[hh, "g ", #1 == #2, ";"];
        WriteString[hh, ".store\n\n"] )&, Flatten[insgen], 1 ] ];
  Write[hh, "l ", Global`Result == ampsum, ";"];
  If[ Length[ampins] =!= 0,
    WriteString[hh, "\
.sort\n\n\
#call FillIns\n\
argument;\n\
#call FillIns\n\
endargument;\n\
#endprocedure\n\n\
#procedure FillIns\n"];
    Apply[Write[hh, "id ", #2 == #1, ";"]&, ins, 1]
  ];

  WriteString[hh, "\
#endprocedure\n\n\
#define Dim \"" <> ToString[dim /. D -> 0] <> "\"\n\
#define OnShell \"" <> ToBin[onshell] <> "\"\n\
#define VADecompose \"" <> ToBin[va] <> "\"\n\
#define Scaled \"" <> ToBin[scale =!= 1] <> "\"\n\
#define FermionChains \"" <> ToBin[!FreeQ[amps, FermionChain]] <> "\"\n\
#define SUNobjs \"" <> ToBin[!FreeQ[amps, SUNobjs]] <> "\"\n\
#define SUNN \"" <> ToString[SUNN] <> "\"\n\n\
#include " <> ToFileName[$FormCalcDir, "OneLoop.h"] <> "\n"];

  If[scale === 1, scale = Sequence[]];
  Block[ {Set},
    (res := Block[#, RunForm])& @
      Cases[DownValues[hide] /. "i_" -> I,
        _[p_, x_Symbol] :> (x = p[[1, 1]])]
  ];
  Amp[CurrentProcess]@@ res[[1]]
]

CalcFeynAmp[___] := (Message[CalcFeynAmp::syntax]; Abort[])


ToBin[True] = "1"

ToBin[___] = "0"


(* FORM interfacing *)

chain[expr__] := (++fline; NonCommutativeMultiply[expr] /. ga -> om)

trace[expr__] := (traces = True;
  NonCommutativeMultiply[expr] /. {
    a_. ga[li_] ** ga[1] + a_. ga[li_] ** ga[-1] :> a "g_"[fline, li],
    a_. ga[1] + a_. ga[-1] :> a "g_"[fline],
    ga -> om }
)

FormWrite[na_, amp_] :=
Block[ {fline = 1},
  Write[hh, "g " <> na == amp /.
    MatrixTrace -> trace /. FermionChain -> chain /.
    NonCommutativeMultiply[a_] -> a, ";"];
  WriteString[hh, ".store\n\n"];
  ampsum += na;
]

FormWrite[na_, amp_, gm_ -> rulz_] :=
Block[ {insna = "Ins" <> na, fline = 1},
  Write[hh, "g ", na@@ gm == amp /.
    MatrixTrace -> trace /. FermionChain -> chain /.
    NonCommutativeMultiply[a_] -> a, ";"];
  WriteString[hh, "\n"];
  insgen = {insgen,
    insna -> Plus@@ Apply[na, Map[Replace[#, ins]&, rulz, {2}], 1]};
  ampsum += insna;
]


FormDecl[_, _[]] = {}

FormDecl[type_, _[f_, v___]] :=
Block[ {l, ll, s},
  ll = StringLength[s = type <> ToString[f]];
  { Fold[
      ( l = StringLength[s = ToString[#2]] + 2;
        {#1, If[(ll += l) > 70, ll = l; ",\n  ", ", "], s} )&,
      s, {v} ],
    ";\n" }
]


FormPattern[ _[_[ sym_Symbol ]], val_ ] :=
  (vars = {vars, sym}; Id[sym -> val])

FormPattern[ _[_[ (h_Symbol)[args__] ]], val_ ] :=
Block[ {patt, lhs, c1 = 0, c2 = 96},
  lhs = h[args] /. Pattern -> (#2 &) /.
    _Blank :> "patt" <> ToString[++c1] <> "?" /.
    _BlankSequence | _BlankNullSequence :> "?" <> FromCharacterCode[++c2];
  vars = {vars, ToSymbol["patt", #]&/@ Range[c1], Cases[lhs, _Symbol]};
  func = {func, h};
  Id[lhs -> val]
]


FormFunc = {ga, Spinor, Den, A0i, B0i, C0i, D0i, E0i,
  SumOver, SUNT, SUNF, SUNTSum}

DeclareVars[expr__] :=
Block[ {theexpr, vars, func = {}, const},
  theexpr = {expr, FormSymbols};
  vars = Cases[theexpr, x_Symbol /; Context[x] =!= "System`",
    Infinity, Heads -> True];
  smalls = StringJoin[ Apply[FormPattern, DownValues[Small], 1] ];
  vars = Complement[Flatten[{vars, Pi}], FormVectors,
    {ga, Spinor, MatrixTrace, FermionChain}];
  func = Union[ Flatten[{func, Sqrt}],
    Select[Flatten[{vars, Conjugate}],
      !FreeQ[theexpr, HoldPattern[Blank[#]]]&] ];
  vars = Complement[vars, FormIndices];
  const = Complement[vars, FormFunc,
    {S, T, U, Sf, Tf, Uf, External, Delta}, GList];
  WriteString[hh,
    FormDecl["i ", Dim/@ FormIndices] <>
    FormDecl["v ", FormVectors] <>
    FormDecl["s ", Complement[vars, func]] <>
    FormDecl["cf ", func] <>
    "#define Const \"" <> StringTake[ToString[const], {2, -2}] <> "\"\n"
  ]
]


toform = "!" <> ToFileName[$FormCalcDir, "ToForm_" <> $Platform] <> " > "

tempnum = 1

OpenFormTemp := (
  While[
    temp = $TemporaryPrefix <> ToString[tempnum] <> ".frm";
    FileType[temp] =!= None, ++tempnum];
  Print[""];
  Print["preparing FORM code in ", temp];
  hh = OpenWrite[toform <> temp,
    FormatType -> InputForm, PageWidth -> 73];
  WriteString[hh, FormSetup];
)

RunForm := (
  Close[hh];
  If[edit, Pause[1]; Run[StringForm[$Editor, temp]]; Pause[3]];
  WriteString["stdout", "running FORM... "];
  res = FormExec["!" <> $FormCmd <> " " <> temp];
  Print["ok"];
  If[!retain, DeleteFile[temp]];
  res
)


(* things to do when the amplitude comes back from FORM *)

num = Identity

kin = MandelstamSimplify

d$ = MetricTensor

pave[A0i[0], args__] := A0[args];
pave[B0i[0], args__] := B0[args];
pave[B0i[1], args__] := B1[args];
pave[B0i[0, 0], args__] := B00[args];
pave[B0i[1, 1], args__] := B11[args];
pave[n_[i__], args__] := n[paveid[n, i], args]

paveid[n_, i__] := paveid[n, i] =
Block[ {t = ToLowerCase[StringTake[ToString[n], 1]]},
  ToSymbol[t, t, i]
]

A0[0] = 0

B0[p_, m1_, m2_] := B0[p, m2, m1] /; !OrderedQ[{m1, m2}]

DB0[p_, m1_, m2_] := DB0[p, m2, m1] /; !OrderedQ[{m1, m2}]

Derivative[1, 0, 0][B0] = DB0;
Derivative[1, 0, 0][B1] = DB1;
Derivative[1, 0, 0][B11] = DB11;
Derivative[1, 0, 0][B00] = DB00


(* abbreviationing business *)

Mat[0] = 0

DiracChain[s___Spinor, 1, r___] := DiracChain[s, r]

sun[x_] := sun[x] = Unique["SUN"]

fme[x_] := ferm[x, scale]

ferm[x__] := ferm[x] = Unique["F"]


p$[a_, x_^n_] := p$[a, x]^n
	(* different operator priority in FORM *)

p$[x__] := pair[Pair[x], scale]

pair[x__] := pair[x] = Unique["Pair"]


e$[x__] := eps[Eps[x], scale]

eps[x__] := eps[x] = Unique["Eps"]


abb[p_Plus] := abbsum[abb/@ p]

abb[n_?NumberQ s_] := n abb[s]

abb[s_?AtomQ] = s

abb[s_] := abb[s] = Unique["Abb"]

abbsum[s_] := abbsum[s] = Unique["AbbSum"]


Abbr[] := Flatten[
  Apply[dv, DownValues/@ {ferm, sun, pair, eps, abb, abbsum}, {2}] ]

Attributes[dv] = {HoldAll}

dv[x_, _] := {} /; !FreeQ[x, Pattern]

dv[_[_[x_]], s_Symbol] := s -> x

dv[_[_[x_, sc_]], s_Symbol] := s -> sc^Count[5 x, _k, {2}] x


ClearProcess[] := (
  CurrentProcess = Sequence[];
  ClearCache[];
  Apply[zap, DownValues/@ {ferm, sun, pair, eps, abb, abbsum}, {2}]; )

Attributes[zap] = {HoldAll}

zap[p_, s_Symbol] := (Unset@@ p; Remove[s]) /; FreeQ[p, Pattern]


Attributes[set] = {Flat, Orderless}

Overlap[] = set[]

Overlap[x__] := Intersection[x]


PlusCSE[{}] = {}

PlusCSE[rul_] :=
Block[ {pl, var, def, i, com, tmp, new = {}},
  Attributes[pl] = {Flat, Orderless};
  Apply[
    (Set@@ {pl@@ -#2, -#1}; Set@@ {pl@@ #2, #1})&,
    Sort[rul, Length[ #1[[2]] ] < Length[ #2[[2]] ] &], 1 ];
  {var, def} = Transpose[
    Cases[DownValues[pl], _[_[_[p__]], s_Symbol] :> {s, set[p]}] ];
  Do[
    While[
      com = Intersection[ def[[i]], def[[i + 1]] ];
      If[ Length[com] < Length[ def[[i]] ]/2,
        tmp = Intersection[ def[[i]], Thread[-def[[i + 1]], set] ];
        If[Length[tmp] > Length[com], com = tmp] ];
      Length[com] > 3,
    (* while body: *)
      tmp = Ceiling[Length[com]/2];
      tmp = Overlap@@ Select[
        set[Intersection[#, com], Intersection[Thread[-#, set], com]]&/@
          Drop[def, {i, i + 1}],
        Length[#] > tmp & ];
      If[Length[tmp] > 3, com = tmp];
      tmp = ToSymbol["help", ++c];
      new = {new, tmp -> Plus@@ com};
      def = def /. {com -> tmp, Thread[-com, set] -> -tmp}
    ],
  {i, Length[def] - 1}];
  Flatten[{new, Thread[var -> Apply[Plus, def, 1]]}]
]


TimesCSE[{}] = {}

TimesCSE[rul_] :=
Block[ {tm, var, def, i, com, tmp, new = {}},
  Attributes[tm] = {Flat, Orderless};
  Apply[ Set@@ {tm@@ #2, #1}&,
    Sort[rul, Length[ #1[[2]] ] < Length[ #2[[2]] ] &], 1 ];
  {var, def} = Transpose[
    Cases[DownValues[tm], _[_[_[t__]], s_Symbol] :> {s, set[t]}] ];
  Do[
    While[
      com = Intersection[ def[[i]], def[[i + 1]] ];
      Length[com] > 3,
    (* while body: *)
      tmp = Ceiling[Length[com]/2];
      tmp = Overlap@@ Select[
        Intersection[#, com]&/@ Drop[def, {i, i + 1}],
        Length[#] > tmp & ];
      If[Length[tmp] > 3, com = tmp];
      tmp = ToSymbol["help", ++c];
      new = {new, tmp -> Times@@ com};
      def = def /. com -> tmp
    ],
  {i, Length[def] - 1}];
  Flatten[{new, Thread[var -> Apply[Times, def, 1]]}]
]


AbbrCat[rul:_[_, _Plus]] := {{}, {}, rul}

AbbrCat[rul:_[_, t_Times]] := {{}, rul, {}} /; FreeQ[t, DiracChain]

AbbrCat[rul_] := {rul, {}, {}}

OptimizeAbbr[rul:{__Rule}] :=
Block[ {c = 0},
  MapThread[ #1[#2]&,
    { {Identity, TimesCSE, PlusCSE},
      Flatten/@ Transpose[AbbrCat/@ rul] } ]//Flatten
]


Cost[x_] :=
Block[ {sx = Simplify[x]},
  {Plus@@ ((Length[#] - 1) Head[#]&)/@
    Level[sx /. Plus -> 1 /. -1 -> 1 /. Times | Power -> 2, {0, -2}], sx}
]

MandelstamSimplify[x_] := x /; FreeQ[x, Plus]

MandelstamSimplify[x_] :=
  Sort[Cost/@ {x, x /. T -> STUSum - S - U,
    x /. U -> STUSum - S - T, x /. S -> STUSum - T - U}][[1, 2]] /;
  ValueQ[STUSum]

MandelstamSimplify[x_] := Simplify[x]


(* UV and IR finiteness checks *)

loopint = A0 | B0 | B1 | B00 | B11 | DB0 | DB1 | DB00 | DB11 | C0i | D0i

loopintArg = Flatten[loopint | Cget | Dget]


UVDivergentPart[expr_] := expr /. int:loopint[__] :> UVDiv[int]

UVDiv[A0[m_]] = m Divergence

UVDiv[_B0] = Divergence

UVDiv[_B1] = -1/2 Divergence

UVDiv[B00[p_, m1_, m2_]] = ((m1 + m2)/2 - p/6)/2 Divergence

UVDiv[_DB00] = -1/12 Divergence

UVDiv[_B11] = 1/3 Divergence

UVDiv[C0i[cc00, __]] = 1/4 Divergence

UVDiv[C0i[cc001 | cc002, __]] = -1/12 Divergence

UVDiv[D0i[dd0000, __]] = 1/24 Divergence

UVDiv[_] = 0


(* FeynCalc-compatibility functions *)

C0[p1_, p2_, p1p2_, m1_, m2_, m3_] :=
  C0i[cc0, p1, p2, p1p2, m1, m2, m3]

PaVe[i__Integer, {p__}, {m1_, m2_, m3_}] :=
  C0i[ToSymbol["cc", Sort[{i}]], p, m1, m2, m3]

D0[p1_, p2_, p3_, p4_, p1p2_, p2p3_, m1_, m2_, m3_, m4_] :=
  D0i[dd0, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4]

PaVe[i__Integer, {p__}, {m1_, m2_, m3_, m4_}] :=
  D0i[ToSymbol["dd", Sort[{i}]], p, m1, m2, m3, m4]


FeynCalcGet[mask___] :=
Block[ {Global`OneLoopResult, Global`GraphName},
  _Global`GraphName = 0;
  Plus@@ ((Get[#]; Global`OneLoopResult[0])&)/@
    FileNames[mask] /. ep_Eps :> I ep
]


FeynCalcPut[expr_, file_] :=
Block[ {PaVe, C0i, D0i, C0, D0},
  C0i[cc0, args__] := C0[args];
  D0i[dd0, args__] := D0[args];
  C0i[i_, p__, m1_, m2_, m3_] := PaVe[
    Sequence@@ (ToExpression/@ Drop[Characters[ToString[i]], 2]),
    {p}, {m1, m2, m3} ];
  D0i[i_, p__, m1_, m2_, m3_, m4_] := PaVe[
    Sequence@@ (ToExpression/@ Drop[Characters[ToString[i]], 2]),
    {p}, {m1, m2, m3, m4} ];
  Put[expr /. _Amp -> List, file]
]


(* helicity matrix elements *)

om[5] := "g5_"[fline]

om[6] := "g6_"[fline]/2

om[7] := "g7_"[fline]/2

om[rho[k[n_], 0, sign_]] :=
  ("g_"[fline] + sign Hel[n] "g5_"[fline]) ** "g_"[fline, k[n]]/2

om[rho[k[n_], m_, sign_]] :=
  ("g_"[fline] + Hel[n] "g5_"[fline] ** "g_"[fline, e[n]]) **
    ("g_"[fline, k[n]] + sign m "g_"[fline])/2

om[rhoc[k[n_], 0, sign_]] :=
  -"g_"[fline, k[n]] ** ("g_"[fline] + sign Hel[n] "g5_"[fline])/2

om[rhoc[k[n_], m_, sign_]] :=
  (-"g_"[fline, k[n]] + sign m "g_"[fline]) **
    ("g_"[fline] - Hel[n] "g_"[fline, e[n]] ** "g5_"[fline])/2

om[li___] := "g_"[fline, li]


Spinor[p_, m_, s_] := Spinor[p, Small[m], s] /; Head[Small[m]] =!= Small


ToTrace[fi_ -> plain_, fj_ -> conj_] :=
Block[ {me, fline = 0},
  me = plain conj //. {
    DiracChain[a___, Spinor[km__]] DiracChain[Spinor[km__], b___] :>
      DiracChain[a, rho[km], b],
	(* If the spinors at the ends don't match directly, we
	   have to reverse one chain. This is a charge conjugation,
	   not a hermitian conjugation like the one HelicityME
	   does for the "conj" part. The rules:
	     a) reverse the chain and exchange u <-> v
	     b) gamma_mu -> -gamma_mu
	     c) add a global minus sign to compensate for the change
	        in the permutation of the external fermions.
	   For details see the Denner/Eck/Hahn/Kueblbeck paper. *)
    DiracChain[Spinor[k1_, m1_, _], a___, Spinor[k2_, m2_, s2_]] *
      DiracChain[Spinor[k1_, m1_, s1_], b___] :>
      -gasign[a] DiracChain[
        Spinor[k2, m2, -s2],
        Sequence@@ Reverse[{a} /. {rho -> rhoc, rhoc -> rho}],
        rho[k1, m1, s1], b ]
  } /.
    DiracChain[Spinor[km__], a___, Spinor[km__]] :>
      (++fline; om/@ (rho[km] ** a));
  Mat[fi, fj] -> {fline, me}
]

gasign[a___] := -1 /;
  OddQ[Length[DeleteCases[{a}, _Integer | _rho | _rhoc]]]

gasign[___] = 1


ToHel[k[n_], __] := {{}, {}} /; Head[Hel[n]] =!= Hel

ToHel[k[n_], 0, s_] :=
Block[ {h = Heli[n]},
  {Hel[n] -> h - s, h -> Hel[n] + s}
]

ToHel[k_, __] := ToHel[k, 0, 0]


Heli[n_] := Heli[n] = ToSymbol["Hel", n]


SelectArg[All] := fabbr

SelectArg[expr_] := Union[Select[fabbr, !FreeQ[ expr, #[[1]] ]&]]


uv[1] = "u"

uv[-1] = "v"

Format[DiracChain[Spinor[p1_, m1_, s1_], c___, Spinor[p2_, m2_, s2_]]] :=
  DiracChain[Overscript[uv[s1], "_"][p1, m1], c, uv[s2][p2, m2]]


ConjChain[s___Spinor, 5, g___] := -Reverse[DiracChain[s, 5, g]]

ConjChain[s___Spinor, 6, g___] := Reverse[DiracChain[s, 7, g]]

ConjChain[s___Spinor, 7, g___] := Reverse[DiracChain[s, 6, g]]

ConjChain[g__] := Reverse[DiracChain[g]]


Options[HelicityME] = {
  AbbrToUse :> Abbr[],
  EditCode -> False,
  RetainFile -> False }

HelicityME::noprocess =
"No process defined so far. HelicityME works only after CalcFeynAmp."

HelicityME::nomat = "Warning: No matrix elements to compute."

HelicityME[plain_, conj_, opt___?OptionQ] :=
Block[ {abbr, edit, retain,
fabbr, part, tohel, fromhel, hels,
FormIndices = {}, ind, c = 0,
temp, hh, res, e, traces, smalls},

  If[{CurrentProcess} === {}, Message[HelicityME::noprocess]; Abort[]];
  {abbr, edit, retain} = ParseOpt[opt, HelicityME];

  Update[Spinor];
  fabbr = Select[abbr, !FreeQ[#, Spinor]&];
  abbr = SelectArg/@ {plain, conj};
  If[ Times@@ Length/@ abbr === 0,
    Message[HelicityME::nomat];
    Return[{}] ];

  part = Cases[abbr, _Spinor, Infinity]//Union;
  {tohel, fromhel} = Flatten/@ Transpose[Apply[ToHel, part, 1]];
  part = #[[1, 1]]&/@ part;
  hels = First/@ fromhel;

  ind = Map[# -> "N" <> ToString[++c] <> "_?" &,
    Union[Cases[#, _Lor, Infinity]]&/@ abbr, {2}];

  traces = Flatten[Outer[ ToTrace, abbr[[1]] /. ind[[1]],
    abbr[[2]] /. ind[[2]] /. DiracChain -> ConjChain /.
      ep_Eps -> -ep /. ep_e :> Conjugate[ep] ]];
  traces = traces /. tohel /. Reverse/@ FromFormRules /. Eps -> "e_";

  OpenFormTemp;
  DeclareVars[Last/@ traces, hels];

  Print["> ", Length[traces], " helicity matrix elements"];

  WriteString[hh, "\
i J1, J2, J3, J4;\n\
v P1, P2;\n\
s scale;\n\
cf abb, kin;\n\n\
.global\n\n" <>
    FormProcs <> "\
#procedure DotSimplify(momsubst, moresimp)\n\
`momsubst'\n\
#call eiki\n\
`moresimp'\n\
.sort\n\
#endprocedure\n\n\
#procedure Simplify()\n\
contract 0;\n\
#call eiei\n\
#call MomConserv(DotSimplify)\n\
polyfun abb;\n" <>
    smalls <>
    If[scale === 1, "", "\
$pow = count_(<k1,1>,...,<k`Legs',1>);\n\
multiply scale^$pow;\n"
    ] <> "\
id P1?.P2? = abb(P1.P2);\n\
id e_(J1?, J2?, J3?, J4?) = abb(e_(J1, J2, J3, J4));\n\
.sort\n\
polyfun;\n\
normalize abb;\n\
id abb(1) = 1;\n" <>
    FormDecl["b ", Append[hels, abb]] <> "\
.sort\n\
collect kin;\n\
normalize kin;\n\
id kin(1) = 1;\n" <>
    FormDecl["b ", hels] <> "\
print;\n\
.store\n\
#endprocedure\n\n"];

  Apply[
    ( Write[hh, "l " <> "Mat" <> ToString/@ List@@ #1 <> " = ",
        #2[[2]], ";"];
      Array[ WriteString[hh, "trace4,", #, ";\n"]&, #2[[1]] ];
      WriteString[hh, "#call Simplify()\n\n"] )&,
    traces, 1 ];

  WriteString[hh, "\n.end\n"];

  (e[#] = s[#])&/@ part;

  Thread[First/@ traces -> Apply[Plus, RunForm, 1]] /. fromhel
]


(* colour matrix elements *)

sunT[a___, 0, 0] := Block[ {c = Unique["col"]}, sunT[a, c, c] ]

sunT[t1___, a_, t2___, a_, t3___, i_, j_] :=
  1/2 (sunT[t1, t3, i, j] sunT[t2, 0, 0] - 1/SUNN sunT[t1, t2, t3, i, j])

sunT/: sunT[t1___, a_, t2___, i_, j_] sunT[t3___, a_, t4___, k_, l_] :=
  1/2 (sunT[t1, t4, i, l] sunT[t3, t2, k, j] -
        1/SUNN sunT[t1, t2, i, j] sunT[t3, t4, k, l])

sunT[i_, i_] := SUNN

sunT/: sunT[a___, i_, j_]^2 := ((1 - 1/SUNN)/2)^Length[{a}] SUNN

sunT/: sunT[a___, i_, j_] sunT[b___, j_, k_] := sunT[a, b, i, k]


ColourFactor[fi_ -> plain_, fj_ -> conj_] :=
  Mat[fi, fj] -> Simplify[plain conj /. SUNT -> sunT]


Options[ColourME] = {AbbrToUse :> Abbr[]}

ColourME::nomat = HelicityME::nomat

ColourME[plain_, conj_, opt___?OptionQ] :=
Block[ {abbr, fabbr},
  {abbr} = ParseOpt[opt, ColourME];
  fabbr = Select[abbr, !FreeQ[#, SUNT]&];
  abbr = SelectArg/@ {plain, conj};
  If[ Times@@ Length/@ abbr === 0,
    Message[ColourME::nomat];
    Return[{}] ];

  Outer[ ColourFactor, abbr[[1]],
    abbr[[2]] /. t_SUNT :> RotateLeft[Reverse[t], 2] ]//Flatten
]


UniquefyIndices[conj_, plain__] :=
Block[ {ind},
  ind = Intersection@@
    (Union[Cases[#, SumOver[x_] -> x, Infinity]]&)/@ {conj, plain};
  conj /. Thread[ind -> (ToSymbol[#, "c"]&)/@ ind]
]


SquaredME[amp_] := SquaredME[amp, amp]

SquaredME[Amp[_][plain__], Amp[_][conj__]] :=
  Plus[plain] Conjugate[UniquefyIndices[Plus[conj], plain]] /;
  FreeQ[{plain, conj}, Mat]

SquaredME[Amp[_][plain__], Amp[_][conj__]] :=
  Flatten[Outer[ToSquared,
    Collect[Plus[plain], _Mat, Hold],
    Collect[UniquefyIndices[Plus[conj], plain], _Mat, Hold]
  ]]


ToSquared[Mat[m1_] x1_., Mat[m2_] x2_.] :=
  ReleaseHold[x1 Conjugate[x2]] ToMat[m1, m2]

ToMat[m1_Symbol, m2_Symbol] := Mat[m1, m2]

ToMat[m1_, m2_] := Inner[Mat, m1, m2, Times]


Unprotect[Conjugate]

Format[ Conjugate[x_] ] := SequenceForm[x, Superscript["*"]]

Format[ Conjugate[t_Times] ] :=
  SequenceForm["(", t, ")", Superscript["*"]]

Conjugate[p:_Plus | _Den] := Conjugate/@ p

Conjugate[sym_?RealQ] := sym

Protect[Conjugate]


(RealQ[#] = True)&/@ {S, T, U, Sf, Tf, Uf}


(* Fortran code generation *)

Attributes[WriteFF] = {HoldFirst, Listable}

WriteFF[s_Symbol, array_] :=
  FFPut[s, array, Block[{s}, ToString[s]]]

WriteFF[amp_, array_] :=
  FFPut[amp, array, ToString[array] <> ToString[++modnum]]

ProcessCheck[p_] := (
  proc = p;
  header = "\n\
* this file is part of the process " <> ToString[p] <> "\n\
* generated by WriteSquaredME " <> TimeStamp[] <> "\n\n\
#include \"prefix.h\"\n\n";
)

ProcessCheck[p_, p_] = 0

ProcessCheck[__] := Message[WriteSquaredME::incomp]

FFPut[Amp[p_][amp__], array_, file_] := (
  ProcessCheck[p, proc];
  mandel = Select[mandel, FreeQ[{amp}, #]&];
  FFPut[#, array, file]&/@ {amp}
)

FFPut[0, _, _] = {}

FFPut[amp_, array_, file_] :=
Block[ {ind},
  ind = Cases[amp, SumOver[i_, r_] :> (Dim[i] = r; i)];
  mods = FFMod[
    amp /. unused[array] -> 0 /. fcs /. {
      _SumOver -> 1,
      Den[p_, m_] -> 1/(p - m),
      int:loopint[__] :> abbint[int] },
    array,
    file <> ({"_", ToString[#]}&)/@ ind ];
  Indices[mods /. List -> Alternatives] = ind;
  mods
]

FFMod[amp_Plus, array_, mod_] :=
Block[ {c = 96},
  FFMod[#, array, mod <> FromCharacterCode[++c]]&/@
    SizeSplit[amp, $FileSize]
] /; LeafCount[amp] > $FileSize

FFMod[amp_, array_, mod_] :=
Block[ {file = mod <> ".F", hh},
  hh = OpenFortran[ModName[file]];
  WriteString[hh,
    "* " <> file <> header <>
    SubroutineDecl[mod, {Pair, Eps}, {Delta}, amp] <>
    "#include \"vars.h\"\n\n"];
  WriteExpr[hh, FFList[amp, array], True,
    Newline -> True, Optimize -> True];
  WriteString[hh, "\tend\n"];
  Close[hh];
  mod
]

FFList[0, _] := Sequence[]

FFList[amp_, array_] :=
  (maxmat[array] = {Mat[1]}; array[1] -> amp) /; FreeQ[amp, Mat]

FFList[amp_Plus, array_] := FFList[#, array]&/@ List@@ amp

FFList[x_. Mat[m_], array_] :=
  ( maxmat[array] = MaxDims[maxmat[array], Level[m, {-2}]];
    Level[m, {-1}, array] -> x )


(* Calculating the abbreviations in a clever way is key to a decent
   performance of the generated code. Therefore, the abbreviations are
   split into three categories:
   1. objects that depend only on model constants and S
      -> subroutine abbr_s,
   2. objects that depend on other phase-space variables (angles etc.)
      -> subroutine abbr_angle,
   3. objects that depend on the helicities
      -> subroutine abbr_hel.
   The final user-callable subroutine squared_me takes care to
   invoke these abbr_nnn subroutines only when necessary. *)

Category[rul_] := {{}, {}, rul} /; !FreeQ[rul[[2]], hels]

Category[rul_] := {{}, rul, {}} /;
  !FreeQ[rul[[2]], T | U | Sf | Tf | Uf | Pair | Eps]

Category[rul_] := {rul, {}, {}}


Dependencies[li_] := {li}

Dependencies[f__, li_] :=
Block[ {pos, c = 0, cc = -1},
  pos = {f};
  Block[ #,
    Apply[(#1 = Indeterminate)&, li, 1];
    While[ c != cc,
      cc = c;
      pos = Apply[dep, pos, {2}] ]
  ]&[ Union@@ Apply[sym, {f, li}, {2}] ];
  pos = Position[pos, {}, {2}, Heads -> False];
  Append[
    Dependencies@@ Delete[{f}, pos],
    Flatten[{li, Extract[{f}, pos]}] ]
]

sym[s_[__], _] = s

sym[s_, _] = s

dep[s_, Indeterminate] := (++c; s = Indeterminate; {})

dep[s___] := {s}


OnePassOrder::recurs =
"Recursive definition in list. Returning list unordered."

OnePassOrder[li_] :=
Block[ {c = 0, cc = -1, l = Length[li], Dep, Ticket, posmap},
  Attributes[Dep] = Attributes[Ticket] = {HoldFirst};
  Ticket[a_, b_] := (a = Random[]; ++c) /; FreeQ[b, Dep];
  Block[ #,
    Apply[(#1 = Dep[#1])&, li, 1];
    posmap = Apply[Ticket, li, 1]
  ]&[ Union[Apply[sym, li, 1]] ];
  While[ c < l,
    If[cc === c, Message[OnePassOrder::recurs]; Return[li]];
    cc = c;
    posmap = posmap ];
  li[[ Level[Sort[MapIndexed[List, posmap]], {3}] ]]
]


DepCats[abbr__] := OnePassOrder/@ Dependencies@@
  Flatten/@ Transpose[Category/@ Flatten[{abbr}]]


AbbrMod[_[], _] = {}

AbbrMod[abbr_, mod_] :=
Block[ {c = 96},
  Flatten[AbbrMod[#, mod <> FromCharacterCode[++c]]&/@
    SizeSplit[abbr, $FileSize]]
] /; LeafCount[abbr] > $FileSize

AbbrMod[abbr_, mod_] :=
Block[ {file = mod <> ".F", hh},
  hh = OpenFortran[ModName[file]];
  WriteString[hh,
    "* " <> file <> header <>
    SubroutineDecl[mod, {Pair, Eps}, {}, abbr] <>
    "#include \"vars.h\"\n\n"];
  WriteDoLoops[hh, abbr, WriteExpr];
  WriteString[hh, "\tend\n"];
  Close[hh];
  {mod}
]


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
  v = Select[DeleteCases[lhs/@ vars, _[0]],
    StringTake[ToString[#], 1] =!= "t" &];
  pindex = Position[v, _[_Symbol..], 1, Heads -> False];
  phead = Position[v, _[__], 1, Heads -> False];
  VarDecl["\n\t" <> type,
    MapAt[Dim/@ # &, v, pindex]] <>
  VarDecl["\n\tcommon /" <> common <> "/ ",
    MapAt[Head, v, phead]] <>
  "\n"
]


SubroutineDecl[name_, cext_, iext_, expr_] :=
Block[ {ce, ie},
  ce = Select[cext, !FreeQ[expr, #]&];
  ie = Select[iext, !FreeQ[expr, #]&];
  { "\tsubroutine ", name, "\n\
\timplicit character (a-s,u-z)\n\
\timplicit double complex (t)\n",
    VarDecl["\n\tdouble complex ", ce],
    VarDecl["\n\tinteger ", ie],
    VarDecl["\n\texternal ", Join[ce, ie]], "\n" }
]


MaxDims[args__] := List@@ MaxIndex@@ Union[Flatten[{args}]]

Attributes[MaxIndex] = {Flat}

MaxIndex[s_[i__], s_[j__]] :=
  MaxIndex[ s@@ MapThread[Max, {{i}, {j}}] ]


	(* LoopComponents gives back e.g.
		1. {F[4], SUN[3]}
		2. {F[j1], SUN[j2]}
		3. Tree[4, 3]
		4. Tree[j1, j2]
		5. "\tdo j1 = 1, 4\n\tdo j2 = 1, 3\n"
		6. "\tenddo\n\tenddo\n"			*)

LoopComponents[arr_] := 0 /; maxmat[arr] === {}

LoopComponents[arr_] :=
  ReplacePart[
    Transpose[LoopVar/@ maxmat[arr]],
    arr, {{3, 0}, {4, 0}} ]

LoopVar[mat:_[1]] := {mat, mat, 1, 1, "", ""}

LoopVar[mat_[max_]] :=
Block[ {var = Unique["j"]},
  { mat[max], mat[var], max, var,
    "\tdo " <> ToString[var] <> " = 1, " <> ToString[max] <> "\n",
    "\tenddo\n" }
]


MatType[_Mat, _] = 1

MatType[h_[i_], h_[j_]] := MatType[h][i, j]

MatType[h_] := MatType[h] = ToSymbol["Mat", h]


Assort[m_Mat -> x_] := {{}, m -> (x /. Hel -> Heli), {}}

	(* convert FormCalc's abbreviations (symbols) to arrays which
	   are easier to handle in Fortran: *)
Assort[f_ -> x_] :=
Block[ {c = Characters[ToString[f]]},
  {{}, {}, (f -> #1[#2])&@@ (ToExpression[StringJoin[#]]&)/@
    {Select[c, LetterQ], Select[c, DigitQ]}}
] /; !FreeQ[x, DiracChain | SUNT]

Assort[x_] := {x /. Conjugate[ep_e] -> -ep /.
  { k[j_] :> 4 (j - 1) + 1,
	(* the spin reference vector of fermions is equivalent
	   to the longitudinal polarization vector of bosons *)
    s[j_] :> 4 (j - 1) + 2,
    e[j_] :> 4 (j - 1) + 2 + Heli[j] }, {}, {}}


DefNeeded[m_[i_, j_]] := (Needed[m[x_, y_] -> _] := x <= i && y <= j)


MandelstamVars = {T, U, Sf, Tf, Uf}

MapThread[
  (MandelstamDef[#1] := MomSquare[ mom[[#2]], -mom[[#3]] ])&,
  {MandelstamVars, {1, 1, 3, 2, 2}, {3, 4, 4, 4, 3}}
]

MomSquare[a_?Negative, b_] := MomSquare[-a, -b]


CheckDir[dir_] :=
Block[ {full},
  Check[
    If[FileType[dir] === None, CreateDirectory[dir]];
    full = SetDirectory[dir],
    Abort[] ];
  ResetDirectory[];
  full
]


Attributes[WriteSquaredME] = {HoldAll}

Options[WriteSquaredME] = {
  SymbolPrefix -> "",
  RenConstFile -> "renconst.F",
  Drivers -> "drivers" }

WriteSquaredME::baddir = "`1` is not a directory."

WriteSquaredME::incomp = "Warning: writing out Fortran code for
incompatible processes."

WriteSquaredME::empty = "Warning: no amplitudes were specified."

WriteSquaredME[tree_, loop_, outdir_String, opt___Rule] :=
  WriteSquaredME[tree, loop, Abbr[], outdir, opt]

WriteSquaredME[tree_, loop_, abbr__, outdir_String, opt___Rule] :=
Block[ {prefix, rconst, drivers,
dir, abr, mat, fcs, proc = Sequence[],
abbint, cints = {}, iints = {}, cc = 0, ic = 0,
ModName, Indices, hels, mandel, mom,
files, hh, n, unused, maxmat, ntree, ntree2, nloop,
header, ffmods, abbrmods, Conjugate = dconjg},

  {prefix, rconst, drivers} = ParseOpt[opt, WriteSquaredME];

  dir = CheckDir[outdir];
  ModName[mod_] := ModName[mod] = ToFileName[dir, mod];

(* abbint introduces abbreviations for the loop integrals.
   They fall into two categories:
   1. A0, B0, B1, B00, B11, DB0, DB1, DB00, DB11 (cint..),
   2. C0i and D0i (iint..).
   For the latter the LoopTools functions Cget and Dget can be used to
   compute all tensor coefficients at once (which is much more efficient).
   Unlike the other integrals, whose results are double complex numbers,
   Cget and Dget return an integer pointing into a cache array. In the
   conventions of LoopTools 2, the actual tensor coefficients are
   retrieved from the arrays Cval and Dval. *)

  abbint[C0i[i_, args__]] :=
  Block[ {uu = ToSymbol["iint", ++ic]},
    If[Length[#] =!= 0, uu = uu@@ #]& @
      Union[Cases[{args} /. Plus | Times -> Sequence, _Symbol, {2}]];
    iints = {iints, uu -> Cget[args]};
    abbint[C0i[id_, args]] = Cval[id, uu];
    Cval[i, uu]
  ];

  abbint[D0i[i_, args__]] :=
  Block[ {uu = ToSymbol["iint", ++ic]},
    If[Length[#] =!= 0, uu = uu@@ #]& @
      Union[Cases[{args} /. Plus | Times -> Sequence, _Symbol, {2}]];
    iints = {iints, uu -> Dget[args]};
    abbint[D0i[id_, args]] = Dval[id, uu];
    Dval[i, uu]
  ];

  abbint[func_] :=
  Block[ {uu = ToSymbol["cint", ++cc]},
    If[Length[#] =!= 0, uu = uu@@ #]& @
      Union[Cases[func /. Plus | Times -> Sequence, _Symbol, {2}]];
    cints = {cints, uu -> func};
    abbint[func] = uu
  ];

  {abr, mat, fcs} = Flatten/@ Transpose[Assort/@ Flatten[{abbr}]];
  n = First/@ DeleteCases[mat, _ -> 0];
  unused[Ctree] = Alternatives@@
    Select[Union[#[[1, 2]]&/@ mat], FreeQ[n, Mat[_, #]]&];
  unused[Cloop] = Alternatives@@
    Select[Union[#[[1, 1]]&/@ mat], FreeQ[n, Mat[#, _]]&];

(* Part 1: the form factors *)

  mandel = Select[MandelstamVars, FreeQ[mat, #]&];
  maxmat[_] = {};
  ffmods = Flatten[{
    Block[{modnum = 0}, WriteFF[tree, Ctree]],
    Block[{modnum = 0}, WriteFF[loop, Cloop]] }];
  If[ Length[ffmods] === 0,
    Message[WriteSquaredME::empty]; Return[{}] ];
  mandel = Complement[MandelstamVars, mandel];

  n = 0;
  mom = Flatten[{#[[1]], -#[[2]]}]&[Map[4 n++ + 1 &, proc, {2}]];

(* Part 2: the variable declarations *)

  iints = Flatten[iints];
  cints = Flatten[cints];
  hels = Alternatives@@ Array[Heli, Plus@@ Length/@ proc];

  ntree = LoopComponents[Ctree];
  nloop = LoopComponents[Cloop];
  ntree2 = LoopComponents[If[ntree === 0, Cloop, Ctree]];
  If[ nloop === 0, nloop = {{}, {}, {}, {}, "", ""} ];
  n = MapThread[MatType,
    {MaxDims[ ntree2[[1]], nloop[[1]] ], ntree2[[1]]}];
  Scan[DefNeeded, n];

  hh = OpenWrite[ModName["vars.h"]];

  WriteString[hh, "\
#include \"model.h\"\n\
#include \"looptools.h\"\n\
#include \"renconst.h\"\n" <>
    CommonDecl[Flatten[{S, mandel}], "double precision ", "kinvars"] <>
    CommonDecl[hels, "integer ", "kinvars"] <>
    CommonDecl[First/@ abr, "double complex ", "abbrev"] <>
    CommonDecl[cints, "double complex ", "loopint"] <>
    CommonDecl[iints, "integer ", "loopint"] <>
    CommonDecl[
      #[[1, 1, 1]]&/@ DownValues[Dim],
      "integer ", "indices"] <>
    CommonDecl[
      Flatten[{ntree2[[3]], nloop[[3]], DeleteCases[n, 1]}]//Union,
      "double complex ", "coeff" ]
  ];
  Close[hh];

(* Part 3: the abbreviations *)

  mat = Select[mat /. fcs /. Mat -> MatType, Needed];
  abbrmods = MapThread[AbbrMod,
    { ToDoLoops/@ DepCats[abr, mat, cints, iints],
      {"abbr_s", "abbr_angle", "abbr_hel"} }];

(* Part 4: global defs in prefix.h *)

  hh = OpenWrite[ModName["prefix.h"]];

  If[ prefix =!= "",
    WriteString[hh, "#define " <> # <> " " <> prefix <> # <> "\n"]&/@
      Flatten[{ffmods, abbrmods, "squared_me",
        "kinvars", "abbrev", "loopint", "indices", "coeff",
        "calc_renconst", "renconst"}] ];

  Close[hh];

(* Part 5: the makefile *)

  hh = OpenWrite[ModName["GNUmakefile.in"]];

  WriteString[hh, "\
OBJS =" <> ({" \\\n  ", #, ".o"}&)/@ Flatten[{abbrmods, ffmods}] <> "\n\n\
RENCONST = " <> StringReplace[rconst, ".F" -> ".o"] <> "\n\n\
ALLOBJS = $(OBJS) $(RENCONST) squared_me.o\n\n\
STDDEPS = prefix.h model.h renconst.h\n\n\
POSSIBLEDEPS = \\\n\
  2to2.F 2to2.h gauss.F \\\n\
  2to3.F 2to3.h multigauss.F vegas.F \\\n\
  sm_ini.F mssm_ini.F\n\n\n\
default: run\n\n\
clean:\n\
\t$(RM) $(ALLOBJS) squared_me.a\n\n\
squared_me.a: squared_me.a($(ALLOBJS))\n\
squared_me.a($(RENCONST)): $(STDDEPS)\n\
squared_me.a($(OBJS) squared_me.o): vars.h $(STDDEPS)\n\n\
renconst.h:\n\
\ttouch renconst.h\n\n\
(%.o): %.F\n\
\t$(FC) -c $<\n\
\t$(AR) cr $@ $%\n\
\t$(RM) $%\n\n\
%:: %.F num.F process.h $(STDDEPS) $(POSSIBLEDEPS) squared_me.a\n\
\t$(FC) -o $@ $< squared_me.a $(LIBS)\n\n"];

  Close[hh];

(* Part 6: the control file squared_me.F *)

  hels = {#, # <> "from", # <> "to"}&/@ ToString/@ List@@ hels;
  n = Flatten[Rest/@ hels];

  hh = OpenFortran[ModName["squared_me.F"]];

  WriteString[hh, "* squared_me.F" <> header <> "\
\tsubroutine squared_me(tree, loop, sqrtS" <>
    VarDecl[",\n     +    ", n] <> ", reset)\n\
\timplicit integer (j)\n\
\timplicit character (a-i,k-z)\n\n\
\tdouble precision tree, loop\n\
\tdouble precision sqrtS" <> VarDecl["\n\tinteger ", n] <> "\n\
\tlogical reset\n\n\
#include \"vars.h\"\n\n" <>
    If[ Length[mandel] === 0, "", "\
\tdouble precision MomSquare\n\
\texternal MomSquare\n\n"
    ] <> "\
\tdouble complex c, m\n\
\tinteger Cptr, Dptr\n\
\tdouble precision prevS\n\
\tsave prevS\n\
\tdata prevS /-1/\n\n\
\tS = sqrtS**2\n"];

  WriteExpr[hh, # -> MandelstamDef[#]&/@ mandel];

  WriteString[hh, "\n\
\tif(reset .or. S .ne. prevS) then\n\
\t  call setcachelast(Ccache, 0)\n\
\t  call setcachelast(Dcache, 0)\
" <> ({"\n\t  call ", #}&)/@ abbrmods[[1]] <> "\n\
\t  reset = .FALSE.\n\
\t  prevS = S\n\
\tendif\n\
\tCptr = getcachelast(Ccache)\n\
\tDptr = getcachelast(Dcache)\n\
" <> ({"\n\tcall ", #}&)/@ abbrmods[[2]] <> "\n\n\
\ttree = 0\n\
\tloop = 0\n\
" <> Apply[{"\n\tdo 1 ", #1, " = ", #2, ", ", #3}&, hels, 1] <> "\n\
" <> ({"\n\tcall ", #}&)/@ abbrmods[[3]] <> "\n\n" <>
    If[ ntree === 0, "", ntree[[5]] <> "\
\t" <> ToFortran[ ntree[[4]] ] <> " = 0\n\
" <> ntree[[6]] <> "\n"
    ] <>
    If[ nloop[[1]] === {}, "", nloop[[5]] <> "\
\t" <> ToFortran[ nloop[[4]] ] <> " = 0\n\
" <> nloop[[6]] <> "\n"
    ]
  ];

  WriteDoLoops[hh, ToDoLoops[ffmods, Indices],
    WriteString[#1, "\tcall " <> #2 <> "\n"]&];

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
\ttree = tree + dble(c*m)\n\n"
    ] <>
    If[ nloop[[1]] === {}, "", "\
\tm = 0\n\
" <> nloop[[5]] <> "\
\tm = m + " <>
      ToFortran[ nloop[[4]] *
        Inner[MatType, nloop[[2]], ntree2[[2]], Times] ] <> "\n\
" <> nloop[[6]] <> "\
\tloop = loop + " <> If[ntree === 0, "", "2*"] <> "dble(c*m)\n\n"
    ] <>
    ntree2[[6]] <> "\n\
1\tcontinue\n\n\
\tcall setcachelast(Ccache, Cptr)\n\
\tcall setcachelast(Dcache, Dptr)\n\
\tend\n"];

  Close[hh];

(* Part 7: copy the driver files *)

  files = {};

  If[ FileType[drivers] === Directory,
    SetDirectory[drivers];
    CopyFile[#, ModName[#]]&/@ (files = FileNames[]);
    ResetDirectory[]
  ];

  If[ FileType[$DriversDir] === Directory,
    SetDirectory[$DriversDir];
    CopyFile[#, ModName[#]]&/@ Complement[FileNames[], files];
    ResetDirectory[]
  ];

  Cases[DownValues[ModName], _[_, s_String] -> s]
]


(* renormalization constants *)

SelfEnergy[proc_, m_] :=
  CalcSelfEnergy[proc, Options[InsertFields]] /. K2 -> m^2

DSelfEnergy[proc_, m_] :=
  D[CalcSelfEnergy[proc, Options[InsertFields]], K2] /. K2 -> m^2

CalcSelfEnergy[proc_, opt_] := CalcSelfEnergy[proc, opt] =
Block[ {se, num = Simplify, Small},
  Needs["FeynArts`"];
  ClearProcess[];
  se = InsertFieldsHook[
    CreateTopologies[1, Length[Flatten[{#}]]&/@ proc,
      ExcludeTopologies -> Internal],
    proc ];
  OptionalPaint[se, $PaintSE];
  se = CreateFeynAmp[se, Truncated -> !FreeQ[proc, F]];
  Plus@@ CalcFeynAmp[se, OnShell -> False, Transverse -> False] //.
    Abbr[] /. {
    Mat -> Identity,
    Pair[_k, _k] -> K2,
	(* take only the transverse part of vector-boson SEs: *)
    Pair[_e, _k] -> If[MatchQ[proc, _V -> _V], 0, 1],
    Pair[_e, _e] -> 1,
    SUNT[_, _] -> 1,
    SUNT[_, _, 0, 0] -> 1/2 }
]

InsertFieldsHook[tops_, proc_] := InsertFields[tops, proc]


ClearSE[] := (DownValues[CalcSelfEnergy] =
  Select[DownValues[CalcSelfEnergy], #[[1, 1, 1, 0]] =!= Pattern &];)


OptionalPaint[ins_, True] := Paint[ins]

OptionalPaint[ins_, path_String] :=
Block[ {file},
  file = path <>
    ToString/@ (Cases[Process /. List@@ Head[ins],
      _Integer | (s_Symbol /; Context[s] === "FeynArts`"),
      {-1}, Heads -> True] /. -1 -> "-") <>
    "_" <> ToString[Model /. List@@ Head[ins]] <> ".ps";
  Paint[ins, DisplayFunction -> (Display[file, #]&)];
]


(* These are special versions of Re and Im where the real and
   imaginary part is taken only of the loop integrals (see A. Denner, 
   Forts. Phys. 41 (1993) 307). *)

ReTilde[expr_] := expr /. int:loopint[__] :> Re[int]

ImTilde[expr_] :=
  (expr /. int:loopint[__] :> Im[int]) - (expr /. loopint[__] -> 0)


	(* Note: it seems weird that the left-handed vector component
	   is taken as the coefficient of DiracChain[6, k]: this is
	   because DiracChain[6, k] = DiracChain[k, 7]. *)

LVectorCoeff[se_] := Coefficient[se, DiracChain[6, k[1]]]

RVectorCoeff[se_] := Coefficient[se, DiracChain[7, k[1]]]

LScalarCoeff[se_] := Coefficient[se, DiracChain[7]]

RScalarCoeff[se_] := Coefficient[se, DiracChain[6]]


IntCollect[p__] := Plus[p] /; FreeQ[{p}, Re]

IntCollect[p__] := Collect[Plus[p], _Re, Simplify]


ExecRenConst[rc_[args___]] := ExecRenConst[rc[args], Options[rc]]

ExecRenConst[rc_] := ExecRenConst[rc, Options[rc]]

ExecRenConst[rc_, {}] := RenConst[rc]

ExecRenConst[rc_, opts_] :=
Block[ {saveopts = Options[InsertFields], res},
  SetOptions[InsertFields, Sequence@@ opts];
  res = RenConst[rc];
  Options[InsertFields] = saveopts;
  res
]


RenConst::nodef =
"Warning: `` might be renormalization constants, but have no definition."

CalcRenConst[expr__] :=
Block[ {patt, rcs = {}, rcdefs = {}, test = {expr}, new},
  patt = Alternatives@@ (#[[1, 1, 1]]&)/@ DownValues[RenConst];
  While[ Length[new = Complement[Cases[test, patt, Infinity], rcs]] =!= 0,
    rcs = Flatten[{rcs, new}];
    test = ExecRenConst/@ new;
    rcdefs = {rcdefs, Thread[new -> test]} ];

  new = Select[ToExpression/@ Names["Global`d*"],
    FreeQ[rcs, #] && !FreeQ[{expr}, #]&];
  If[ Length[new] =!= 0, Message[RenConst::nodef, new] ];

  Expand[Sort[ TakeApart/@ Flatten[rcdefs] ]] /. Plus -> IntCollect
]

TakeApart = Identity


NElements[rc_[i__]] := rc[Times[i]]

NElements[rc_] := rc[1]


Options[WriteRenConst] = {RenConstFile -> "renconst.F"}

WriteRenConst::norcs = "Warning: no renormalization constants found."

WriteRenConst[rcdefs:{___Rule}..., outdir_String, opt___Rule] :=
Block[ {rconst, include, dir, file, ifile, header, hh, arr,
rcs = Flatten[{rcdefs}]},

  {rconst} = ParseOpt[opt, WriteRenConst];

  dir = CheckDir[outdir];
  include = StringReplace[rconst, ".F" -> ""] <> ".h";
  header = "\n* generated by WriteRenConst " <> TimeStamp[] <> "\n";

  If[ Length[rcs] === 0, Message[WriteRenConst::norcs],
    rcs = OnePassOrder[rcs] ];

  hh = OpenFortran[file = ToFileName[dir, rconst]];
  WriteString[hh, "* " <> rconst <> "\n\
* this file contains renormalization constants" <> header <> "\n\
#include \"prefix.h\"\n\n" <>
    SubroutineDecl["calc_renconst", {}, {Delta}, rcs] <> "\
#include \"model.h\"\n\
#include \"looptools.h\"\n\
#include \"" <> include <> "\"\n" <>
    VarDecl[ "\n\tinteger ",
      Union[Cases[rcs, SumOver[i_, _] -> i, Infinity]] ] <>
    "\n\n"];
  WriteDoLoops[hh, ToDoLoops[rcs], WriteSummedExpr];
  WriteString[hh, "\tend\n"];
  Close[hh];

  rcs = MaxDims[First/@ rcs] /. m_Max :> Dim/@ m;
  arr = NElements/@ rcs;
  hh = OpenFortran[ifile = ToFileName[dir, include]];
  WriteString[hh, "* " <> include <> "\n\
* the declarations for " <> rconst <> header <>
    CommonDecl[rcs, "double complex ", "renconst"] <> "\n\
\tinteger sizeof_rc\n\
\tparameter (sizeof_rc = " <> ToFortran[Level[arr, {2}, Plus]] <> ")\n" <>
    If[ Length[arr] === 0, "", "\
\tdouble complex rc(sizeof_rc)\n\
\tequivalence (" <> ToFortran[ arr[[1, 0]] ] <> ", rc)\n" ]
  ];
  Close[hh];

  {file, ifile}
]

WriteRenConst[expr__, outdir_String, opt___Rule] :=
Block[ {Dim, TakeApart = SplitSums},
  Dim[i_Integer] = i;
  Cases[{expr},
    SumOver[i_, r_, ___] | IndexSum[_, {i_, r_}] :> (Dim[i] = r),
    Infinity];
  WriteRenConst[CalcRenConst[expr], outdir, opt]
]


(* low-level Fortran output functions *)

ToFortran[s_String] = s

ToFortran[s_] := ToString[FortranForm[s]]


tofortran =
  "!" <> ToFileName[$FormCalcDir, "ToFortran_" <> $Platform] <> " > "

OpenFortran[file_] := OpenWrite[tofortran <> file,
  FormatType -> FortranForm, PageWidth -> 67]


TimeStamp[] :=
  ToString[#3] <> " " <>
  {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
   "Sep", "Oct", "Nov", "Dec"}[[#2]] <> " " <>
  ToString[#1] <> " " <>
  ToString[#4] <> ":" <>
  StringTake["0" <> ToString[#5], -2]&@@ Date[]


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
  List@@ Flatten[Operate[Coalesce, expr]]
]

Coalesce[h_][a_, b_, r___] :=
  Coalesce[h][{a, b}, r] /; LeafCount[{a, b}] < maxsize

Coalesce[h_][a_, r___] := cb[ h@@ Flatten[{a}], Coalesce[h][r] ]

Coalesce[_][] = Sequence[]


Options[WriteExpr] = {Newline -> False, Optimize -> False}

WriteExpr[fi_, expr_, addtovar_:False, opt___Rule] :=
Block[ {hh = fi, newline, optim, addto,
Conjugate = dconjg, Re = dble, Im = dimag},
  {newline, optim} = ParseOpt[opt, WriteExpr];
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

WriteAssign[var_ -> 0] := 0 /; addto[var]

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
  for = expr /. int:loopintArg[__] :> Nargs/@ int /.
    p:_Integer^_Rational :> N[p] /.
    Times -> OptTimes;
  If[ TrueQ[addto[var]], for += var, addto[var] = True ];
  Write[hh, var -> for];
  newline
]

Nargs[0] = 0.

Nargs[i_Integer] := N[i]

Nargs[x_] = x


Unprotect[Rule]

Format[a_ -> b_, FortranForm] := SequenceForm[a, " = ", b]

Protect[Rule]


OptTimes[a__] :=
Block[ {p, const, var},
  p = Position[N[{a}], _Real, 1, Heads -> False];
  const = Times@@ {a}[[ Flatten[p] ]];
  If[ IntegerQ[const], Return[Times[a]] ];
  var = Times@@ Delete[{a}, p];
  If[ var === 1, Return[Times[a]] ];
  If[ MatchQ[const, _?Negative _.],
    -HoldForm[HoldForm[#1] #2]&[-const, var],
    HoldForm[HoldForm[#1] #2]&[const, var] ]
]


WriteSummedExpr[hh_, var_ -> li_List, opt___Rule] :=
Block[ {SumOver, Dim, loops, addto, svar = ToFortran[var], si},
  loops = ToDoLoops[li];
  If[ loops[[1, 0]] === DoLoop,
    WriteString[hh, "\t" <> ToFortran[var] <> " = 0\n"];
    addto = True,
  (* else *)
    addto := (addto = True; False) ];
  SumOver[i_, r_] := (Dim[i] = r; 1);
  WriteDoLoops[hh, loops,
    WriteExpr[#1, var -> #2, addto, Optimize -> True]&];
  WriteString[hh,
    "\n#ifdef DEBUG\n\tprint *, '" <>
    StringReplace[svar, Cases[var,
      i_Symbol :> (si = ToFortran[i]; si -> "'," <> si <> ",'")]] <>
    " =', " <> svar <> "\n#endif\n\n" ];
]


SplitSums[var_ -> expr_] := var -> SplitSums[expr]

SplitSums[li_List] := SplitSums[Plus@@ li]

SplitSums[x_] := {x} /; FreeQ[x, SumOver]

SplitSums[x_] :=
Block[ {term},
  term[_] = 0;
  assign[Expand[x, SumOver]];
  term[_] =.;
  #[[1, 1, 1]] Plus@@ Flatten[ #[[2]] ]&/@ DownValues[term]
]

assign[p_Plus] := assign/@ p

assign[t_Times] := (term[#1] = {term[#1], #2})&@@ cull/@ t

assign[other_] := term[1] = {term[1], other}

cull[s_SumOver] := {s, 1}

cull[other_] := {1, other}


FindIndices[var_ -> _] := Union[Cases[var, _Symbol]]

FindIndices[t_Times] := Cases[t, SumOver[i_, _] -> i]

FindIndices[_] = {}

ToDoLoops[h_[li__], func_:FindIndices] :=
Block[ {do},
  do[_] = {};
  Apply[(do[#1] = {do[#1], #2})&, {func[#], #}&/@ {li}, 1];
  Cases[ DownValues[do],
    _[_[_[ind_List]], a_] :> DoLoop[ind, h@@ Flatten[a]] ]
]

ToDoLoops[x_, ___] := Flatten[{x}]


DoLoop[_[], {a___}] = a

DoLoop[_[], a_] = a


Attributes[WriteDoLoops] = {Listable}

WriteDoLoops[hh_, DoLoop[ind_, expr_], func_] := (
  WriteString[hh,
    {"\n\tdo ", ToFortran[#], " = 1, ", ToFortran[Dim[#]]}&/@ ind <>
    "\n"];
  WriteDoLoops[hh, expr, func];
  WriteString[hh, StringJoin[Table["\tenddo\n", {Length[ind]}]] ];
)

WriteDoLoops[hh_, expr_, func_] := func[hh, expr]

End[]


Format[ Continuation[_] ] = "    "
  (* eliminate those `>' in front of continuation lines so one can cut
     and paste more easily *)

$FormCmd = "form3"
  (* the filename of the actual FORM executable; may contain a path *)

FormSetup = "\
#-\n\
#:SmallSize 5000000\n\
#:LargeSize 10000000\n\
#:WorkSpace 500000\n\
#:MaxTermSize 30000\n\
#:TermsInSmall 30000\n\
#:TempDir " <> DirectoryName[$TemporaryPrefix] <> "\n\
off stats;\n\
format 255;\n\n"

$Editor = "xterm -geometry 80x30 -e pico `1` &"
  (* editor to use when debugging FORM code *)

$BlockSize = 700

$FileSize = 30 $BlockSize

$DriversDir = ToFileName[{$FormCalcDir, "drivers"}]

$PaintSE = False

EndPackage[]


(* global definitions for specific models *)

(* definitions for the Standard Model *)

FormSubst = "\
id CW^2 = CW2;\n\
id CW^-2 = CW2^-1;\n\
id SW^2 = SW2;\n\
id SW^-2 = SW2^-1;\n\
id EL^2 = 4*Pi*Alfa;\n\
id Alfa^2 = Alfa2;\n\
id GS^2 = 4*Pi*Alfas;\n\
id Alfas^2 = Alfas2;\n"

Conjugate[CKM[a__]] ^:= CKMC[a];
Conjugate[CKMC[a__]] ^:= CKM[a]

EL/: EL^2 = 4 Pi Alfa;
EL/: EL^4 = 16 Pi^2 Alfa2;
Alfa^(n_?EvenQ) ^= Alfa2^(n/2)

GS/: GS^2 = 4 Pi Alfas;
GS/: GS^4 = 16 Pi^2 Alfas2;
Alfas^(n_?EvenQ) ^= Alfas2^(n/2)

SW^(n_?EvenQ) ^= SW2^(n/2);
CW^(n_?EvenQ) ^= CW2^(n/2)

CW2/: CW2 + SW2 = 1

MZ^(n_?EvenQ) ^= MZ2^(n/2);
MW^(n_?EvenQ) ^= MW2^(n/2);
MH^(n_?EvenQ) ^= MH2^(n/2)

ME^(n_?EvenQ) ^= ME2^(n/2);
MM^(n_?EvenQ) ^= MM2^(n/2);
ML^(n_?EvenQ) ^= ML2^(n/2);
MLE[a__]^(n_?EvenQ) ^= MLE2[a]^(n/2)

MU^(n_?EvenQ) ^= MU2^(n/2);
MC^(n_?EvenQ) ^= MC2^(n/2);
MT^(n_?EvenQ) ^= MT2^(n/2);
MQU[a__]^(n_?EvenQ) ^= MQU2[a]^(n/2)

MD^(n_?EvenQ) ^= MD2^(n/2);
MS^(n_?EvenQ) ^= MS2^(n/2);
MB^(n_?EvenQ) ^= MB2^(n/2);
MQD[a__]^(n_?EvenQ) ^= MQD2[a]^(n/2)

SUNN = 3

(* these symbols represent real quantities, i.e. Conjugate[sym] = sym
   for any of these. Thinking e.g. of complex masses this looks
   dangerous but then again it's easy to remove any such definition.
   The function that really needs this is SquaredME. *)

Scan[ (RealQ[#] = True)&,
  { EL, Alfa, Alfa2, GS, Alfas, Alfas2,
    MW, MW2, MZ, MZ2,
    SW, CW, SW2, CW2,
    MH, MH2, MG0, MG02, MGp, MGp2,
    ME, ME2, MM, MM2, ML, ML2, _MLE, _MLE2,
    MU, MU2, MC, MC2, MT, MT2, _MQU, _MQU2,
    MD, MD2, MS, MS2, MB, MB2, _MQD, _MQD2 } ]

(* Model parameters which are defined using the parameter statement in
   Fortran (i.e. as numeric constants; see model.h) are given some
   numeric value here. Using this information, the OptTimes function can
   significantly optimize the generated Fortran code. The idea is to put
   everything that is known as constant at compile time in one place,
   i.e. rearrange products such that they are of the form (const)*(vars),
   then the compiler will usually collect all of these constants into one
   number. *)

Scan[ (N[#] = Random[])&,
  { Alfa, Alfa2, SW2, CW, CW2,
    MW, MW2, MZ, MZ2,
    ME, ME2, MM, MM2, ML, ML2,
    MU, MU2, MC, MC2, MT, MT2,
    MD, MD2, MS, MS2, MB, MB2 } ]


(* definitions for the MSSM *)

SetOptions[CalcFeynAmp,
  NoExpand -> {USf, USfC, UCha, UChaC, VCha, VChaC, ZNeu, ZNeuC}]

USf[t_, g_][a_, b_] := USf[a, b, t, g]

Conjugate[USf[a__]] ^:= USfC[a];
Conjugate[UCha[a__]] ^:= UChaC[a];
Conjugate[VCha[a__]] ^:= VChaC[a];
Conjugate[ZNeu[a__]] ^:= ZNeuC[a]

Conjugate[USfC[a__]] ^:= USf[a];
Conjugate[UChaC[a__]] ^:= UCha[a];
Conjugate[VChaC[a__]] ^:= VCha[a];
Conjugate[ZNeuC[a__]] ^:= ZNeu[a]

SA^(n_?EvenQ) ^= SA2^(n/2);
CA^(n_?EvenQ) ^= CA2^(n/2);
SB^(n_?EvenQ) ^= SB2^(n/2);
CB^(n_?EvenQ) ^= CB2^(n/2);
TB^(n_?EvenQ) ^= TB2^(n/2)

CA2/: CA2 + SA2 = 1;
CB2/: CB2 + SB2 = 1

MGl^(n_?EvenQ) ^= MGl2^(n/2);
MSf[a__]^(n_?EvenQ) ^= MSf2[a]^(n/2);
MCha[a__]^(n_?EvenQ) ^= MCha2[a]^(n/2);
MNeu[a__]^(n_?EvenQ) ^= MNeu2[a]^(n/2)

Mh0^(n_?EvenQ) ^= Mh02^(n/2);
MHH^(n_?EvenQ) ^= MHH2^(n/2);
MA0^(n_?EvenQ) ^= MA02^(n/2);
MHp^(n_?EvenQ) ^= MHp2^(n/2)

Scan[ (RealQ[#] = True)&,
  { TB, CB, SB, CA, SA, CB2, SB2, C2A, S2A, CAB, SAB, CBA, SBA,
    Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2,
    _MSf, _MSf2, _MCha, _MCha2, _MNeu, _MNeu2 } ]

Null

