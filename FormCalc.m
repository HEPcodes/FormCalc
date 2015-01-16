(*

This is FormCalc, Version 1.2
Copyright by Thomas Hahn 1998
last modified 11 Oct 98 by Thomas Hahn

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
   somewhere in your documentation that you've used
   our code.

If you're a lawyer, you can find the legal stuff at
http://www.fsf.org/copyleft/lgpl.html.

The user guide for this program can be found at
http://www-itp.physik.uni-karlsruhe.de/formcalc

If you find any bugs, or want to make suggestions, or
just write fan mail, address it to:
	Thomas Hahn
	Institut fuer Theoretische Physik
	Universitaet Karlsruhe
	e-mail: hahn@particle.physik.uni-karlsruhe.de

To join the FormCalc mailing list, send a mail (any text) to
	hahn-formcalc-subscribe@particle.physik.uni-karlsruhe.de

Have fun!

*)

Print[""];
Print["FormCalc 1.2"];
Print["by Thomas Hahn"];
Print["last revision: 11 Oct 98"];


BeginPackage["FormCalc`"]

FeynAmp::usage = "FeynAmp[gname, mom, amp] is the FeynArts way of writing
a Feynman amplitude amp with name gname and integration momentum mom."

FeynAmpList::usage = "FeynAmpList[info][amps] is the head of a list of
FeynAmp objects. info contains additional information handed down by
FeynArts about the process, the momenta etc."

Process::usage = "Process -> {proc} contains the process specification
handed down by FeynArts in the information field of a FeynAmpList."

PropagatorDenominator::usage = "PropagatorDenominator[p, m] is the
FeynArts expression for 1/(p^2 - m^2)."

FeynAmpDenominator::usage = "FeynAmpDenominator[prden..] is the head
wrapped around the PropagatorDenominators inside a loop."

Spinor::usage = "Spinor[p, m] is the spinor with momentum p and mass m.
Spinor corresponds to the more conventional way of writing spinors by\n
   spinor left  in SpinorChain + momentum incoming -> \\bar v\n
   spinor left  in SpinorChain + momentum outgoing -> \\bar u\n
   spinor right in SpinorChain + momentum incoming -> u\n
   spinor right in SpinorChain + momentum outgoing -> v."

LeptonSpinor = Spinor

QuarkSpinor = Spinor

DiracSlash::usage = "DiracSlash[p] represents p_mu gamma_mu."

DiracMatrix::usage = "DiracMatrix[mu] represents the Dirac matrix with
Lorentz index mu."

ChiralityProjector::usage = "ChiralityProjector[+-1] represents the
chirality projectors omega_{+-} = (1 +- gamma_5)/2."

MetricTensor::usage = "MetricTensor[mu, nu] represents the metric tensor
with Lorentz indices mu and nu."

IndexDelta::usage = "IndexDelta[i1, i2] is used to force indices to be the
same. For integers i1 and i2, IndexDelta is 1 if i1 = i2 or 0 otherwise.
If i1 is a symbol, i2 in the expression multiplied by IndexDelta is
replaced by i1 and vice versa."

Index::usage = "Index[t, n] represents an index of type t with number n."

SpinorChain::usage = "SpinorChain[Spinor[...], ..., Spinor[...]]
represents an open fermion chain."

DiracTrace::usage = "DiracTrace represents the trace over Dirac matrices
in a fermion loop."

tr::usage = "tr[i, t] represents the Dirac trace over t with running
number i."

ga::usage = "ga[li] represents the gamma matrix with Lorentz index li."

ga5::usage = "ga5 represents gamma_5."

LI::usage = "LI[n] is the nth Lorentz index in spinor chains. Summation
over all LI is implied."

FourVector::usage = "FourVector[p, mu] represents the four-vector p with
Lorentz index mu."

PolarizationVector::usage = "PolarizationVector[p, mu] represents the
polarization vector belonging to the momentum p with Lorentz index mu."

Momentum::usage = "Momentum[p] is an obsolete form to represent the
four-momentum p."

q1::usage = "q1 is the integration momentum in a one-loop amplitude."

ScalarProduct::usage = "ScalarProduct[p1, p2] is the scalar product of two
four-vectors p1 and p2. It is converted to Pair[p1, p2] in FormCalc."

S::usage = "S is the Mandelstam variable s. If p1 and p2 denote the
incoming momenta, S = (p1 + p2)^2."

T::usage = "T is the Mandelstam variable t. If p1 denotes the first
incoming and k1 the first outgoing momentum, T = (p1 - k1)^2."

U::usage = "U is the Mandelstam variable u. If p1 denotes the first
incoming and k2 the second outgoing momentum, U = (p1 - k2)^2."

Sf::usage = "Sf is an extended Mandelstam variable for a 2 -> 3 reaction.
If k1 and k2 denote the first two outgoing momenta, Sf = (k1 + k2)^2."

Tf::usage = "Tf is an extended Mandelstam variable for a 2 -> 3 reaction.
If p2 denotes the second incoming and k2 the second outgoing momentum,
Tf = (p2 - k2)^2."

Uf::usage = "Uf is an extended Mandelstam variable for a 2 -> 3 reaction.
If p2 denotes the second incoming and k1 the first outgoing momentum,
Uf = (p2 - k1)^2."

DEN::usage = "DEN[p2, m2] stands for 1/(p2 - m2). Note that in contrast
to PropagatorDenominator p2 and m2 are the momentum and mass _squared_."

xi::usage = "xi[v] is the gauge parameter of gauge boson v."

OneLoop::usage = "OneLoop[amps] is the basic function for evaluating
Feynman diagrams. It takes as an argument either a single diagram (with
head FeynAmp) or a list of diagrams and returns the results in a collected
and abbreviated form."

ProcessFile::usage = "ProcessFile[amps, outtag] evaluates the diagrams in
amps by running them through OneLoop and splits the result into a bosonic
part which is written to outtag.m and a fermionic part which is written
to outtagF.m. outtag may contain a pathname together with an identifier
for the diagrams, e.g. destpath/box for box diagrams."

Comment::usage = "Comment -> \"string\" is an option of OneLoop and places
an extra comment at the beginning of the temporary file handed to Form.
This may be useful when debugging Form code."

FormCommands::usage = "FormCommands -> \"cmd1\\ncmd2...\" is an option
of OneLoop and ProcessFile and specifies extra Form commands to be
executed for each diagram."

CancelEps::usage = "CancelEps is an option of OneLoop and ProcessFile. It
is used for CP invariant processes which must not contain terms
proportional to Eps. Setting CancelEps -> True will cancel these terms at
an early stage thus speeding up the calculation (they should drop out in
the end in any case). The default is CancelEps -> False."

EditCode::usage = "EditCode is an option of OneLoop and ProcessFile. It
edits the temporary file passed to Form using $Editor and is of course
used only for debugging."

RetainFile::usage = "RetainFile is an option of OneLoop. When set to
True, it prevents OneLoop from deleting the temporary file containing the
Form input after running Form."

ReduceME::usage = "ReduceME is an option of OneLoop and ProcessFile. When
set to True, FormCalc tries to reduce the number of independent matrix
elements by eliminating one momentum via momentum conservation.
Empirically, bosonic amplitudes are shortened but fermionic ones are not.
If getting the shortest possible result is your aim, you might want to set
ReduceME -> Automatic (the default) which is equivalent to True for
bosonic and False for fermionic amplitudes."

ToughStuff::usage = "ToughStuff is an option of ProcessFile. When set to
True, ProcessFile will split the output not only into fermionic and
bosonic parts, but also into different topologies and (if applicable)
different fermion families. It is mainly used for large amplitudes."

FermionFamily::usage = "FermionFamily[patt] = \"string\" returns the
identifier appended to the filename by ProcessFile if ToughStuff -> True
is set, the amplitude is fermionic and it contains patt. Typically used
like FermionFamily[ME | MM | ML] = \"e\" which gives files containing
leptonic amplitudes an extra e in the name."

ReadForm::usage = "ReadForm[\"file\"] reads file which it assumes to
contain Form output into Mathematica.\n
ReadForm[\"!cmd\"] executes cmd and pipes its output back into
Mathematica."

r2::usage = "r2 is a symbol used to substitute Sqrt[2] before handing an
expression to Form. It is used internally only."

Abbreviations::usage = "Abbreviations[] returns a list of all
abbreviations introduced so far. It is typically used at the end of a
calculation to save the abbreviations with a command like
Abbreviations[] >> abbrfile."

UseAbbreviations::usage = "UseAbbreviations[abbr] registers the
abbreviations abbr with FormCalc. abbr may be either the abbreviations
themselves (a list of rules) or a filename containing them. This is useful
when continuing a former FormCalc session. UseAbbreviations requires that
the value of $O1ME that was used in the former session is also set in the
current session.";

OptimizeAbbreviations::usage = "OptimizeAbbreviations[abbr] optimizes the
set of abbreviations returned by Abbreviations[] by eliminating common
subexpressions."

Scale::usage = "Scale is a scale introduced via $O1ME to scale matrix
elements, e.g. to make them dimensionless. The variable Scale itself
appears in the Abbreviations[]."

Pair::usage = "Pair[a, b] represents the Minkovskian scalar product of
the four-vectors a and b."

Eps::usage = "Eps[a, b, c, d] represents -I times the total antisymmetric
Levi-Civita tensor: Eps[a, b, c, d] = -I a[mu] b[nu] c[rho] d[sigma]
epsilon[mu, nu, rho, sigma], where a, b, c, d are four-vectors and the
sign convention is epsilon[0, 1, 2, 3] = +1."

e::usage = "e[i] is the ith polarization vector."

k::usage = "k[i] is the ith momentum."

DiagramType::usage = "DiagramType[diag] returns 2 for a self-energy
diagram, 1 for a vertex diagram, and 0 for a box. In other words, it
returns the number of denominators not containing the integration
momentum."

FermionicQ::usage = "FermionicQ[diag] gives True for a diagram containing
fermions and False otherwise."

Pick::usage = "Pick[amp, no] picks diagrams from amp where no is of the
form {5} which picks diagram # 5, {{5, 9}} which picks diagrams # 5
through 9, or All which picks all diagrams. Combinations are also allowed,
for example, Pick[amp, {6, {10, 12}, 27}] picks diagrams # 6, 10, 11, 12,
and 27 from amp."

FeynCalcGet::usage = "FeynCalcGet[mask] reads files produced with
FeynCalc. mask is taken as input to the Mathematica function FileNames, so
it might be FeynCalcGet[\"file.m\"] or FeynCalcGet[\"*.m\"] or
FeynCalcGet[\"*.m\",\"~/feyncalcfiles\"]."

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

a0::usage = "a0[m] is used internally to represent a pre-form of the
one-point function."

A0::usage = "A0[m] is the one-point scalar Passarino-Veltman function
where m is the mass squared."

b0::usage = "b0[m, p] is used internally to represent a pre-form of the
two-point function."

B0::usage = "B0[p, m1, m2] is the scalar two-point Passarino-Veltman
function where p is the external momentum squared and m1 and m2 are the
masses squared."

B1::usage = "B1[p, m1, m2] is the tensor two-point Passarino-Veltman
function B_1 where p is the external momentum squared and m1 and m2 are
the masses squared."

B00::usage = "B00[p, m1, m2] is the tensor two-point Passarino-Veltman
function B_00 where p is the external momentum squared and m1 and m2 are
the masses squared."

B11::usage = "B11[p, m1, m2] is the tensor two-point Passarino-Veltman
function B_11 where p is the external momentum squared and m1 and m2 are
the masses squared."

c0::usage = "c0[p1, p2, m1, m2, m3] is used internally to represent a
pre-form of the three-point function."

C0i::usage = "C0i[id, p1, p2, p1p2, m1, m2, m3] is the generic three-point
Passarino-Veltman function which includes both scalar and tensor integrals
specified by id. For example, C0i[cc0, ...] is the scalar function C0,
C0i[cc112, ...] the tensor function C_112 etc. Call the external momenta
k1...k3, then the arguments are given as p1 = k1^2, p2 = k2^2,
p1p2 = (k1 + k2)^2, and m1...m3 are the masses squared."

d0::usage = "d0[p1, p2, p3, m1, m2, m3, m4] is used internally to
represent a pre-form of the four-point function."

D0i::usage = "D0i[id, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the
generic Passarino-Veltman four-point function which includes both scalar
and tensor integrals specified by id. For example, D0i[dd0, ...] is the
scalar function D0, D0i[dd1233, ...] the tensor function D_1233 etc.
Call the external momenta k1...k4, then the arguments are given as
p1 = k1^2, p2 = k2^2, p3 = k3^2, p4 = k4^2, p1p2 = (k1 + k2)^2,
p2p3 = (k2 + k3)^2, and m1...m4 are the masses squared."

o1::usage = "o1 is a head wrapped around the prefactor of a standard
matrix element. Its default value is o1 = Identity (doing nothing). To
get the shortest possible amplitude, you can set o1 = Simplify or
similar."

o2::usage = "o2 is a head wrapped around a linear combination of standard
matrix elements. Its default value is o2 = Identity (doing nothing). To
get the shortest possible amplitude, you can set o2 = Simplify or
similar."

$Editor::usage = "$Editor specifies the editor used in debugging Form
code."

$TempFile::usage = "$TempFile is the name for the temporary file used to
write out Form code in OneLoop. Unless RetainFile is set to True, it is
deleted after usage."

$OnShell::usage = "$OnShell specifies whether the external particles in
your calculation are on shell. Caution: This variable must be set *before*
loading amplitudes."

$FormCalcDir::usage = "$FormCalcDir is the path where FormCalc
(specifically, form.h, form.set and pave.frm) resides."

$Platform::usage = "$Platform is a string that identifies the platform
FormCalc is running on. It is the value of the environment variable
HOSTTYPE in the Bourne shell (sh) and is used to distinguish the ReadForm
binaries of different platforms."

$FormCmd::usage = "$FormCmd gives the name of the actual Form executable.
It may contain a path."

$Transversality::usage = "$Transversality specifies whether FormCalc may
assume transversality for polarization vectors, i.e. e[i].k[i] = 0."

$O1ME::usage = "Matrix elements are scaled with $O1ME for every momentum
they include, e.g. e1.k2 is abbreviated by k12/$O1ME, e1.k2*e2.k3 by
k12*k23/$O1ME^2, etc. This way, all matrix elements can be made
dimensionless. For instance, in a 2 -> 2 process, $O1ME = Sqrt[S] makes
all matrix elements of O(1)."

$Dimension::usage = "$Dimension specifies the dimension FormCalc works in.
It can be D for dimensional regularization and 4 for dimensional reduction
and constrained differential renormalization."


Begin["`Private`"]

(* generic functions *)

$FormCalcDir = StringReplace[
  If[ FileType[$Input] === File, $Input,
        (* if FormCalc was loaded from a directory in $Path: *)
    Block[ {full},
      Scan[
        If[ FileType[full = # <> "/" <> $Input] === File, Return[full] ]&,
        $Path ] ]
  ],
  {"FormCalc.m" -> "", "~" -> HomeDirectory[]} ]

$Platform = Environment["HOSTTYPE"]

If[ Head[$Platform] =!= String, $Platform = "" ]


If[ $LinkSupported,
  Install[$FormCalcDir <> "ReadForm_" <> $Platform],
(* else *)
  Print["WARNING: Your Mathematica Kernel does not support MathLink."];
  Print["This means that the main functions OneLoop and ProcessFile"];
  Print["will not work."];
  ReadForm[_] = $Failed ]


Off[General::spell1, General::spell]

ToSymbol[x__] := ToExpression[ StringJoin[ToString/@ Flatten[{x}]] ]


Pick[amp_, graphs_] :=
  amp[[ If[graphs === All, Range[Length[amp]],
    Flatten[ If[Head[#] === List, Range@@ #, #]&/@ graphs ]//Union] ]]

Attributes[DiagramType] = {Listable};
DiagramType[ FeynAmp[___, a_] ] := Exponent[a /. DEN[__] -> DEN, DEN]

FermionicQ[a_] :=
  !FreeQ[a, Spinor | DiracMatrix | DiracSlash | ChiralityProjector |
    DiracTrace | SpinorChain | "g_" | "g5_"]


PropagatorDenominator[p_, m_] :=
  PropagatorDenominator[-p, m] /; !FreeQ[p, -q1]

PropagatorDenominator[0, m_] = -1/m^2

FeynAmpDenominator[ p:PropagatorDenominator[q1, _], b___ ] :=
  den[p, b];
FeynAmpDenominator[ a__, p:PropagatorDenominator[q1, _], b___ ] :=
  den[p, b, a]

den[ a___, PropagatorDenominator[q_, m1_], b___,
  PropagatorDenominator[q_, m2_], c___ ] :=
  1/Factor[m1^2 - m2^2] (den[a, PropagatorDenominator[q, m1], b, c] -
    den[a, PropagatorDenominator[q, m2], b, c]);
den[ PropagatorDenominator[q1, m1_] ] :=
       I Pi^2 a0[m1^2];
den[ PropagatorDenominator[p1_, m1_],
     PropagatorDenominator[p2_, m2_] ] :=
       I Pi^2 b0[p2 - p1, m1^2, m2^2];
den[ PropagatorDenominator[p1_, m1_],
     PropagatorDenominator[p2_, m2_],
     PropagatorDenominator[p3_, m3_] ] :=
       I Pi^2 c0[p2 - p1, p3 - p1, m1^2, m2^2, m3^2];
den[ PropagatorDenominator[p1_, m1_],
     PropagatorDenominator[p2_, m2_],
     PropagatorDenominator[p3_, m3_],
     PropagatorDenominator[p4_, m4_] ] :=
       I Pi^2 d0[p2 - p1, p3 - p1, p4 - p1, m1^2, m2^2, m3^2, m4^2]


Attributes[IndexDelta] = {Orderless}
        
IndexDelta[n_Integer, n_Integer] = 1
        
IndexDelta[_Integer, n_Integer] = 0
        
FeynAmp[ args__, amp_ ] :=
Block[ {res, Index},
  IndexDelta[n_, sym_] := sym = n;
  res = FeynAmp[args, amp];
  IndexDelta[n_, sym_] = .;
  res
] /; !FreeQ[amp, IndexDelta]


(* Fermion stuff: Dirac algebra etc *)

ds[ind_, n_?NumberQ k_] := n ds[ind, k]

ds[ind_, k_Symbol] := "g_"[ind, k]

DiracSlash[p_] := -DiracSlash[-p] /; !FreeQ[p, -q1]

Spinor[-p_, m_] := Spinor[p, -m]


FermionResolve[n__, a_] :=
Block[ {fline = 0, DiracTrace, SpinorChain},
  DiracTrace[x__] := (++fline; tr[fline, HandleFermions[fline, x]]);
  SpinorChain[x__] := (++fline; HandleFermions[fline, x]);
  FeynAmp[n, a] /. d_Dot :> SpinorChain@@ d /; FermionicQ[d]
] /; FermionicQ[a]

FermionResolve[a__] := FeynAmp[a]


HandleFermions[index_, expr__] :=
Block[ {DiracMatrix, Dot = NonCommutativeMultiply},
  DiracMatrix[5] = "g5_"[index];
  DiracMatrix[li_] = "g_"[index, li];
  NonCommutativeMultiply[expr] /.
    f1_. DiracSlash[k_] + m_. :>
      f1 Distribute[ds[index, k]] + m "gi_"[index] /.
    f1_. ChiralityProjector[1] + f2_. ChiralityProjector[-1] + rr_. :>
      (f1/2 + f2/2 + rr) "gi_"[index] + (f1 - f2)/2 "g5_"[index] /.
    f1_. ChiralityProjector[s_] + m_. :>
      (f1/2 + m) "gi_"[index] + s f1/2 "g5_"[index]
]


(* Boson stuff: 4-vectors etc *)

FourVector[ksum_Plus, li_] := FourVector[#, li]&/@ ksum

FourVector[n_?NumberQ k_, li_] := n FourVector[k, li]

FourVector[k_ q_xi, li_] := q FourVector[k, li]

FourVector[k_, li_] := k[li]


Attributes[MetricTensor] = {Orderless}

MetricTensor/: MetricTensor[li__]^2 :=
  MetricTensor[li] ** MetricTensor[li]


Momentum[p_] = p


PolarizationVector[n_?NumberQ p_, li_] := n PolarizationVector[p, li]

PolarizationVector[p_, li_] :=
  PolarizationVector[p,li] = ToSymbol["e", p][li]


ScalarProduct = Pair

Attributes[Pair] = {Orderless}

(* due to different operator priority in Form: *)
Pair[a_, x_^n_] := Pair[a, x]^n

Pair[a_, x_Plus] := Pair[a, #]&/@ x

Pair[a_, n_?NumberQ x_] := n Pair[a, x]

(* must say a_ (not _) here, or def will appear in Abbreviations[] *)
Pair[a_, 0] = 0

Pair[e[i_], e[j_]] := Pair[e[i], e[j]] = ToSymbol["e", i, j];
Pair[e[i_], k[j_]] := Pair[e[i], k[j]] = ToSymbol["k", i, j]

MomSquare[0] = 0

MomSquare[p_] := Pair[p, p]


Eps[___, a_, a_, ___] = 0

(* must say a___ (not ___) here, or def will appear in Abbreviations[] *)
Eps[a___, 0, ___] = 0

Eps[a___, x_Plus, b___] := Eps[a, #, b]&/@ x

Eps[a___, n_?NumberQ x_, b___] := n Eps[a, x, b]

Eps[a__] := Signature[{a}] Eps@@ Sort[{a}] /; !OrderedQ[{a}]

Eps[e[i1_], e[i2_], e[i3_], e[i4_]] := Eps[e[i1], e[i2], e[i3], e[i4]] =
  ToSymbol["e", i1, i2, i3, i4];
Eps[e[i1_], e[i2_], e[i3_], k[i4_]] := Eps[e[i1], e[i2], e[i3], k[i4]] =
  ToSymbol["e", i1, i2, i3, "k", i4];
Eps[e[i1_], e[i2_], k[i3_], k[i4_]] := Eps[e[i1], e[i2], k[i3], k[i4]] =
  ToSymbol["e", i1, i2, "k", i3, i4];
Eps[e[i1_], k[i2_], k[i3_], k[i4_]] := Eps[e[i1], k[i2], k[i3], k[i4]] =
  ToSymbol["e", i1, "k", i2, i3, i4]


ga[0] = 0

Unprotect[NonCommutativeMultiply];
___ ** 0 ** ___ = 0;
Protect[NonCommutativeMultiply]


(* Abbreviationing business *)

FME[x_] := x /; AtomQ[x] || !FreeQ[x, ga[_Symbol]];
FME[x_] := FME[x] = Unique["F"]

SME[x_] :=
Block[ {c = Cases[x, ga[li_Symbol] :> li, Infinity], once},
  SME[ x /. MapIndexed[#1 -> LI@@ #2 &,
    once[l_] := (once[l] = Sequence[]; l); once/@ c] ] /; Length[c] =!= 0
];
SME[x_?AtomQ] := x;
SME[-x_] := -SME[x];
SME[x_] := SME[x] = Unique["O"]

SMEPLUS[x_?AtomQ] := x;
SMEPLUS[x_] := (SMEPLUS[x] =
Block[ {ft = FactorTerms[x]},
  If[ Head[ft] === Times, ft[[1]] SMEPLUS[ ft[[2]] ], ft ]
]) /; Union[Head/@ List@@ x] === {Times};
SMEPLUS[x_] := SMEPLUS[x] = Unique["P"]


dvhold[sym_] :=
  Cases[ DownValues[sym], (ab_ :> p_Symbol) :>
    (p -> HoldForm@@ ab Scale^Count[ab, k[_], {-2}]) /;
    FreeQ[ab, Pattern] ];

dv[sym_] :=
  Cases[ DownValues[sym], (_[_[ab_]] :> p_Symbol) :>
    (p -> ab Scale^Count[ab /. Spinor[__] -> 1, k[_], {-2}]) /;
    FreeQ[ab, Pattern] ];

Abbreviations[] :=
Block[ {ab = Join[dvhold[Pair], dvhold[Eps],
dv[FME], dv[SME], dv[SMEPLUS]]},
  If[ $O1ME === 1, ab /. Scale -> 1,
    Prepend[ab, Scale -> $O1ME^-1]]
]


UseAbbreviations::warning =
"Warning: there are ambiguous definitions."

UseAbbreviations::scale =
"You must set $O1ME to `1` to use these abbreviations."

UseAbbreviations[s_String] := UseAbbreviations[Get[s]]

UseAbbreviations[abbr_] :=
Block[ {sc, ab},
  sc = 1/ReleaseHold[Scale /. abbr /. Scale -> 1];
  If[ sc =!= $O1ME,
    Message[UseAbbreviations::scale, sc];
    Return[$Failed] ];
  ab = DeleteCases[abbr, Scale -> _] /. Scale -> 1;
  ab /. HoldForm -> Identity;
  Apply[
    If[ AtomQ[#1],
      If[ FreeQ[#2, Plus], SME[#2] = #1, SMEPLUS[#2] = #1 ],
      If[#1 =!= #2, Message[UseAbbreviations::warning] ]
    ]&,
    Select[ab, FreeQ[#, HoldForm]&], 1 ]
]


IntSec[] = 0;
IntSec[a_, b_] := 0 /; Head[a] =!= Plus || Head[b] =!= Plus;
IntSec[a__] := Block[ {i = Intersection[a]}, If[Head[i] === Plus, i, 0] ]

TempInsert[li_, t_] := Insert[li, t,
  Max[0, Position[First/@ li, #, {1}, 1]&/@
    DeleteCases[Level[t[[2]], {-1}], -1]] + 1]

OptLevel12[rul_] :=
Block[ {l = 0, pl, rl = Length[rul], i, is, iss, isl, pr,
repl, nurul = {}},
  Print["level 1: redundancy removal"];
  Attributes[pl] = {Flat, Orderless};
  Apply[
    ( WriteString["stdout", "\r", --rl, " "];
      Set@@ {pl@@ -#2, -#1};
      Set@@ {pl@@ #2, #1})&,
    Sort[rul, Length[ #1[[2]] ] < Length[ #2[[2]] ] &], 1 ];
  repl = Cases[DownValues[pl],
    (_[_[p__]] :> o_Symbol) :> (o -> Plus[p])];

  Print["\rlevel 2: common subexpression elimination"];
  rl = Length[repl] - 1;
  Do[
    pr = False;
    While[
      is = IntSec[ repl[[i, 2]], repl[[i + 1, 2]] ];
      If[ Length[is] < Length[ repl[[i, 2]] ]/2,
        iss = IntSec[repl[[i, 2]], -repl[[i + 1, 2]] ];
        If[Length[is] < Length[iss], is = iss]
      ];
      Length[is] > 3,
    (* while body: *)
      isl = Ceiling[Length[is]/2];
      iss = IntSec@@ Select[
        (IntSec[#[[2]], is] + IntSec[-#[[2]], is])&/@
          Drop[repl, {i, i + 1}],
        Length[#] > isl &];
      If[Length[iss] < 4, iss = is];
      is = ToSymbol["t", ++l];
      AppendTo[nurul, is -> iss];
      repl = repl /. iss -> is /. -iss -> -is;
      pr = True;
    ];
    If[pr, WriteString["stdout", "\r", rl - i, " "]],
  {i, rl}];
  WriteString["stdout", "\r     \r"];
  Fold[TempInsert, repl, nurul]
]

OptimizeAbbreviations[file_String] :=
Block[ {abbr = Get[file], new, orig = file <> ".orig"},
  new = OptimizeAbbreviations[abbr];
  If[ FileType[orig] === File, DeleteFile[orig] ];
  RenameFile[file, orig];
  Put[new, file];
]

OptimizeAbbreviations[rul_] :=
  Join[ Select[rul, FreeQ[#, Plus]&],
    OptLevel12[Select[rul, !FreeQ[#, Plus]&]] ]


(* Things to do when amplitude comes back from Form *)

(* for maximum simplification set o1 = o2 = Simplify or similar *)
o1 = o2 = Identity

Global`fme[x_] := FME[x /. FromFormRules];
Global`fme[x__] := FME[NonCommutativeMultiply[x] /. FromFormRules]

Global`d$ = MetricTensor

Global`i$ = I

Global`gi$[_] = Sequence[]

Global`pave4[i__Integer, p1_, p2_, p3_, m1_, m2_, m3_, m4_] := D0i[
  ToSymbol["dd", Sort[{i}]],
  MomSquare[p1], MomSquare[p1 - p2], MomSquare[p2 - p3],
  MomSquare[p3], MomSquare[p2], MomSquare[p1 - p3],
  m1, m2, m3, m4 ]

Global`pave3[i__Integer, p1_, p2_, m1_, m2_, m3_] := C0i[
  ToSymbol["cc", Sort[{i}]],
  MomSquare[p1], MomSquare[p1 - p2], MomSquare[p2],
  m1, m2, m3 ]

Global`B0m[p_, m1_, m2_] := B0[MomSquare[p], m1, m2];
Global`B1m[p_, m1_, m2_] := B1[MomSquare[p], m1, m2];
Global`B00m[p_, m1_, m2_] := B00[MomSquare[p], m1, m2];
Global`B11m[p_, m1_, m2_] := B11[MomSquare[p], m1, m2]

B0[p_, m1_, m2_] := B0[p, m2, m1] /; !OrderedQ[{m1, m2}]

a0[0] = 0


(* Reading and analyzing FeynArts files *)

FormMandelstam = FormMomcons = ""

FormSymbols = FromFormRules = FormVectors = {}


SetPair22[p1_, p2_, s_, kin_] :=
  p1/: Pair[p1, p2] = s kin - s MomSquare[p1] - s MomSquare[p2]

SetPair12[mom_, i_] :=
Block[ {lhs = Pair@@ Drop[mom, {i}], rhs},
  rhs = Plus@@
    Array[ If[# == i, 1, -1]/2 MomSquare[ mom[[#]] ] &, Length[mom]];
  If[MatchQ[lhs, -_], lhs = -lhs; rhs = -rhs];
  TagSet@@ {lhs[[1]], lhs, rhs}
]


FeynAmpList::noFAfile = "Hold it! This is no FeynArts amplitude."

FeynAmpList[h__][a__] :=
Block[ {proc, mom, ps, eps, dp, mdp, dot, x, x1, x2, s, t, u,
Conjugate = Identity, PropagatorDenominator, MetricTensor = "d_"},
  proc = Process /. {h};
  If[ Head[proc] =!= Rule,
    Message[FeynAmpList::noFAfile];
    Return[$Failed] ];

  mom = #[[2]]&/@ Join[proc[[1]], -proc[[2]]];
  proc = Select[Join@@ proc, Head[ #[[2]] ] === Symbol &];
  Scan[Clear, ps = #[[2]]&/@ proc];

  dot[x1___, -p_, x2___] := -dot[x1, p, x2];
  eps = (dot[ x = ToSymbol["e", #], # ] = 0; x)&/@ ps;
  FormVectors = Join[{q1}, eps, ps];
  FromFormRules = Join[
    MapIndexed[#1 -> k@@ #2 &, ps], MapIndexed[#1 -> e@@ #2 &, eps] ];

  FormSymbols = {};
  If[ $OnShell,
    Apply[UpSet, {MomSquare[ #[[2]] ], #[[3]]^2}&/@ proc, 1];
    FormSymbols = Union[Cases[Last/@ proc, x_Symbol :> x^2]];
  ];

	(* introduce kinematical invariants: *)
  Switch[ Length[ps],
    2,
      Set@@ Reverse[ps],
    3,
      SetPair12[mom, 1];
      SetPair12[mom, 2];
      SetPair12[mom, 3],
    4 | 5,
      {s, t, u} =
        If[ Length[ps] === 4, {S, T, U}, {Sf, Tf, Uf} ];
      FormSymbols = Join[FormSymbols, Union[{S, T, U, s, t, u}]];
      SetPair22[ ps[[1]], ps[[2]], 1/2, S ];
      SetPair22[ ps[[3]], ps[[1]], -1/2, T ];
      SetPair22[ ps[[3]], ps[[2]], -1/2, U ];
      SetPair22[ ps[[3]], ps[[4]], 1/2, s ];
      SetPair22[ ps[[4]], ps[[2]], -1/2, t ];
      SetPair22[ ps[[4]], ps[[1]], -1/2, u ]
  ];

  dp = If[ $Transversality,
    Apply["id " <> ToString[#1] <> " . " <> ToString[#2] <> " = 0;\n" &,
      Transpose[{eps, ps}], 1],
    {} ];
  mdp = Cases[ UpValues[#] /. Pair -> Dot, (_[p_] :> sp_) :>
    "id " <> ToString[p//InputForm] <> " = " <>
      ToString[sp//InputForm] <> ";\n" /; FreeQ[p, Pattern] ]&/@ ps;
  If[ Length[ps] > 4,
	(* eliminate last mom *)
    x = Solve[Plus@@ mom == 0, ps[[-1]] ][[1, 1]];
    UpSet@@ {Pair[x[[1]], x1_], Pair[x[[2]], x1]};
    dp = Insert[Reverse[dp],
      {mdp[[-1]], "id " <> ToString[#1] <> " = " <>
        ToString[#2//InputForm] <> ";\n" &@@ x}, 2];
    mdp = Drop[mdp, -1] ];
  FormMandelstam = StringJoin[dp, mdp];

  FormMomcons = "";
  If[ Length[mom] > 2,
    mom = Array[ {Plus@@ Delete[mom, #], -mom[[#]]}&, Length[mom] ];
    ( If[ $Transversality && AtomQ[ #[[2]] ],
        x = #;
        Scan[ If[ (dp = dot[ #, x[[2]] ]) =!= 0,
          FormMomcons = FormMomcons <> "id " <>
            ToString[dp /. dot -> Dot] <> " = " <>
            ToString[Distribute[dot[ #, x[[1]] ]] /. dot -> Dot] <>
            ";\n"] &,
          eps ];
        FormMomcons = FormMomcons <> ".sort\n" ];
      x = Select[#[[1]], AtomQ, 1];
      If[ x =!= 0, TagSet@@ Prepend[#, x]] )&/@ Join[mom, -mom]
  ];

  PropagatorDenominator[pp_, mm_] :=
    PropagatorDenominator[pp, mm] = DEN[MomSquare[pp], mm^2];
  {a} /. FeynAmp -> FermionResolve /. (1/xi[f_])^r_Rational :> xi[f]^-r
]

ToFormRules = {
  2^(1/2) -> r2,
  2^(-1/2) -> 1/r2,
  Complex[a_, b_] -> a + "i_" b
}


(* Form interface and main function *)

FormWrite[ FeynAmp[gname_, ___, ampl_] ] :=
Block[ {Pair = Dot, na = Last[gname], trc = {}},
  WriteString["stdout", na, " "];
  If[ amplist === na || MemberQ[amplist, na],
    na = Unique[ToString[na] <> "n"] ];
  amplist += na;
  WriteString[hh, "* ", gname, "\ng ", na, " = "];
  Write[hh, ampl /. tr[i_, ferm_] :> (AppendTo[trc, i]; ferm), ";"];
  If[ Length[trc] =!= 0,
    WriteString[hh, ".sort\n"];
    Scan[WriteString[hh, "trace4,", #, ";\n"]&, Union[trc]] ];
  If[ canceps && !FreeQ[ampl, "g5_"],
    WriteString[hh, "if( count(e_, 1) > 0 ) discard;\n"] ];
  WriteString[hh, formcmd, "\n.store\n\n"]
]

FormDecl[type_, _[v___, l_]] :=
Block[ {llen = Length[type]},
  WriteString[hh, type];
  (WriteString[hh, #, ", "];
    If[ (llen += StringLength[ToString[#]] + 2) > 70,
      llen = 0; WriteString[hh, "\n  "] ])&/@ {v};
  WriteString[hh, l, ";\n"];
]

Options[OneLoop] = {
  Comment -> "",
  FormCommands -> "",
  CancelEps -> False,
  EditCode -> False,
  RetainFile -> False,
  ReduceME -> Automatic }

OneLoop[ amps_, opt___Rule ] :=
Block[ {hh, amplist, com, edcode, retain, formcmd, canceps, cancqp,
redme, res, vars, func, indx},
  {com, edcode, retain, formcmd, canceps, redme} =
    {Comment, EditCode, RetainFile, FormCommands,
      CancelEps, ReduceME} /. {opt} /. Options[OneLoop];
  If[ redme === Automatic, redme = !FermionicQ[amps] ];
  Print["preparing output for Form in ", $TempFile];
  hh = OpenWrite[
    "!sed -e '/^[^#]/s/\"//g' -e 's/\\[/(/g' -e 's/]/)/g' " <>
    "-e 's/\\*\\*/\\*/g' -e 's/\\*\\\\//g' > " <> $TempFile,
    FormatType -> InputForm, PageWidth -> 73];

  res = Flatten[{amps}] /. ToFormRules /.
    NonCommutativeMultiply[a_] -> a;
	(* NCM _is_ OneIdentity, only Mma doesn't seem to notice *)
  amplist = {FormSymbols, Last/@ res /. tr -> List};
  vars = Cases[amplist,
    x_Symbol /; Context[x] =!= "System`",
    Infinity, Heads -> True]//Union;
  vars = Append[Complement[vars, FormVectors], Pi];
  func = Select[vars, !FreeQ[amplist, #[__]]&];
  indx = Intersection[vars, ToExpression/@ Names["li*"]];
  vars = Complement[vars, func, indx];

  WriteString[hh,
    "#-\n* ", com,
    "\n#define DIM \"", $Dimension /. D -> 0, "\"\n",
    If[$Dimension === D, "s D;\nd D", "d " <> ToString[$Dimension] ],
    ";\n"];
  FormDecl["i ", indx];
  FormDecl["v ", FormVectors];
  FormDecl["s ", vars];
  FormDecl["cf ", DeleteCases[func, Spinor]];
  WriteString[hh,
    "f Spinor;\n.global\n\n#procedure mandelstam()\n",
    FormMandelstam,
    "#endprocedure\n\n#procedure momcons()\n",
    If[redme, FormMomcons, ""],
    "#endprocedure\n\n"];
  amplist = 0;
  Scan[FormWrite, res];
  WriteString[hh, "l amp = "];
  Write[hh, amplist, ";"];
  WriteString[hh, "#include ", $FormCalcDir, "pave.frm\n\n.end\n"];
  Close[hh];
  If[edcode, Pause[1]; Run[StringForm[$Editor, $TempFile]]; Pause[3]];

  WriteString["stdout", "\nrunning Form... "];
  res = Block[ {Dot = Pair, r2 = Sqrt[2]},
    ReadForm[ToString[
      StringForm["!`1` -s `2`form.set `3`",
        $FormCmd, $FormCalcDir, $TempFile] ]] 
  ] /. FromFormRules /. Global`e$ -> Eps /.
    Global`sme -> SME /. Global`smeplus -> SMEPLUS;
  If[ res === $Failed,
    Print["warning: Form reports an error!!!"],
  (* else *)
    Print["ok"];
    If[!retain, DeleteFile[$TempFile]]
  ];
  res
]


(* This is for backward compatibility with FeynCalc *)

C0[p1_, p2_, p1p2_, m1_, m2_, m3_] :=
  C0i[cc0, p1, p2, p1p2, m1, m2, m3]

PaVe[i__Integer, {p__}, {m1_, m2_, m3_}] :=
  C0i[ToSymbol["cc", Sort[{i}]], p, m1, m2, m3]

D0[p1_, p2_, p3_, p4_, p1p2_, p2p3_, m1_, m2_, m3_, m4_] :=
  D0i[dd0, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4]

PaVe[i__Integer, {p__}, {m1_, m2_, m3_, m4_}] :=
  D0i[ToSymbol["dd", Sort[{i}]], p, m1, m2, m3, m4]

StandardMatrixElement[a_] StandardMatrixElement[b_] ^:=
  StandardMatrixElement[a b];
StandardMatrixElement[a_] := StandardMatrixElement[a /.
  StandardMatrixElement -> Identity] /; !FreeQ[a, StandardMatrixElement]

fcpair[e[i_], e[j_]] := fcsme[FCPair[e[i], e[j]]];
fcpair[e[i_], k[j_]] := $O1ME fcsme[FCPair[e[i], k[j]]];
fcpair[a__] := Pair[a]

fceps[a__] := Signature[{a}] (fceps@@ Sort[{a}]) /; !OrderedQ[{a}];
fceps[e[i1_], e[i2_], e[i3_], e[i4_]] :=
  I fcsme[FCEps[e[i1], e[i2], e[i3], e[i4]]];
fceps[e[i1_], e[i2_], e[i3_], k[i4_]] :=
  I $O1ME fcsme[FCEps[e[i1], e[i2], e[i3], k[i4]]];
fceps[e[i1_], e[i2_], k[i3_], k[i4_]] :=
  I $O1ME^2 fcsme[FCEps[e[i1], e[i2], k[i3], k[i4]]];
fceps[e[i1_], k[i2_], k[i3_], k[i4_]] :=
  I $O1ME^3 fcsme[FCEps[e[i1], k[i2], k[i3], k[i4]]]

fcsme/: fcsme[a_] fcsme[b_] := fcsme[a b]

FeynCalcGet[mask___] :=
Block[ {Global`OneLoopResult, Global`GraphName, Global`SW, Global`CW,
Pair = fcpair, Eps = fceps, StandardMatrixElement = Identity},
  Global`SW/: Global`SW^-2 = Global`S2;
  Global`SW/: Global`SW^-4 = Global`S4;
  Global`CW/: Global`CW^-2 = Global`C2;
  Global`CW/: Global`CW^-4 = Global`C4;
  Global`GraphName[__] = 0;
  Plus@@ ((Print[#]; Get[#]; Global`OneLoopResult[0])&)/@
    FileNames[mask] /.
    (p:S | T | U - m_)^(n_?Negative) :> DEN[p, m]^(-n) /.
    (p:S | T | U)^(n_?Negative) :> DEN[p, 0]^(-n) /.
    (m_ - p:S | T | U)^(n_?Negative) :> -DEN[p, m]^(-n)
] /. FCPair -> Pair /. FCEps -> Eps /. fcsme -> SME


FeynCalcPut[expr_, file_] :=
Block[ {PaVe, C0i, D0i, C0, D0},
  C0i[Global`cc0, args__] := C0[args];
  D0i[Global`dd0, args__] := D0[args];
  C0i[id_, p1_, p2_, p1p2_, m1_, m2_, m3_] := PaVe[
    Sequence@@ (ToExpression/@ Characters[StringDrop[ToString[id], 2]]),
    {p1, p2, p1p2}, {m1, m2, m3}];
  D0i[id_, p1_, p2_, p3_, p4_, p1p2_, p2p3_, m1_, m2_, m3_, m4_] := PaVe[
    Sequence@@ (ToExpression/@ Characters[StringDrop[ToString[id], 2]]),
    {p1, p2, p3, p4, p1p2, p2p3}, {m1, m2, m3, m4}];
  Put[expr, file]
]


(* Classifying amplitudes *)

ClassIndex[ FeynAmp[_[_, top_, ___], ___, a_] ] :=
Block[ {s = If[FermionicQ[a], "F", ""]},
  If[tough && !FreeQ[a, d0],
    If[s === "F",
      s = StringJoin[s, Cases[a,
            m_Symbol "gi_"[_] :> FermionFamily[m], Infinity]//Union] ];
    s = s <> ToString[top];
  ];
  s
];

Options[ProcessFile] = {
  FormCommands -> "",
  CancelEps -> False,
  EditCode -> False,
  ToughStuff -> False,
  ReduceME -> Automatic }

ProcessFile[filename_String, outtag_, opt___Rule] :=
  ProcessFile[Get[filename], outtag, opt]

ProcessFile[amps_List, outtag_, opt___Rule] :=
Block[ {class, na, tough, oo},
  tough = ToughStuff /. {opt} /. Options[ProcessFile];
  oo = Sequence@@ Select[{opt}, FreeQ[#, ToughStuff]&];
  class[_] = {};
  Scan[AppendTo[class[ClassIndex[#]], #]&, amps];
  Scan[
    (Print[""];
     Print["processing class ", na = outtag <> #[[1, 1, 1]] <> ".m"];
     Put[OneLoop[#[[2]], Comment -> na, oo], na];
     Share[])&,
    Select[DownValues[class], FreeQ[#, {}]&]];
]

End[]


Format[ Continuation[_] ] = "    "
  (* eliminate those `>' in front of continuation lines so one can cut
     and paste more easily *)

$OnShell = True
  (* whether the external particles are on-shell *)

$Transversality = True
  (* whether to enforce transversality on the external polarization
     vectors, i.e. e[i].k[i] = 0 *)

$Dimension = D

$FormCmd = "form"
  (* the filename of the actual Form executable; may contain a path *)

$Editor = "xterm -geometry 80x30 -e pico `1` &"
  (* which editor to use when debugging Form code *)

$O1ME = 1
  (* factor which scales SMEs, see usage *)

$TempFile = "t" <> ToString[$SessionID] <> ".frm"

EndPackage[]


(* global definitions for specific models *)

(* definitions specially for the Standard Model *)

EL/: EL^4 = 16 Pi^2 a2

GaugeXi = xi

(* MLA = mla = 0; *)
MW^(n_?EvenQ) ^= MW2^(n/2);
MP^(n_?EvenQ) ^= MP2^(n/2);
MZ^(n_?EvenQ) ^= MZ2^(n/2);
MH^(n_?EvenQ) ^= MH2^(n/2);
ME^(n_?EvenQ) ^= ME2^(n/2);
MD^(n_?EvenQ) ^= MD2^(n/2);
MU^(n_?EvenQ) ^= MU2^(n/2);
MM^(n_?EvenQ) ^= MM2^(n/2);
MS^(n_?EvenQ) ^= MS2^(n/2);
MC^(n_?EvenQ) ^= MC2^(n/2);
ML^(n_?EvenQ) ^= ML2^(n/2);
MB^(n_?EvenQ) ^= MB2^(n/2);
MT^(n_?EvenQ) ^= MT2^(n/2)

(* these are the identifiers for fermion classes used by ProcessFile *)

FermionFamily[ME | MM | ML] = "E";    (* leptons *)
FermionFamily[MU | MC | MT] = "U";    (* up-type quarks *)
FermionFamily[MD | MS | MB] = "D";    (* down-type quarks *)
FermionFamily[_] = ""

Null

